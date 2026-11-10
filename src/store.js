/* ============================================================================
   PlasmidBench — local persistence.

   Everything this page keeps lives in IndexedDB in the browser it is opened
   in: the construct library, the primer and plasmid stock, the enzyme
   preferences, the lab profile. Nothing leaves the machine, which is what
   makes the page deployable on a lab server and what keeps a tool that is
   not a GLP record system from pretending to be one.

   WHY NOT localStorage
   It held the library until now, and it is a synchronous string store of
   five to ten megabytes for the whole origin — settings, stock lists and
   sequences together. A single 100 kb plasmid is 100 kb of JSON; a few
   dozen of them with their feature tables, plus an oligo inventory of a
   few thousand rows, reaches the ceiling, and a localStorage write that
   does not fit throws where it stands. IndexedDB is asynchronous, orders of
   magnitude larger, and indexable.

   WHY THE SETTINGS STILL READ SYNCHRONOUSLY
   Around a hundred call sites read a setting in the middle of drawing
   something — which tab was open, which enzyme panel, the reaction
   conditions. Making each of those wait for a transaction would mean
   rewriting the interface around promises for no benefit the person can
   see. So settings are held in memory, read synchronously from there, and
   written through to IndexedDB in the background. The map is seeded from
   localStorage at start-up, which is already synchronous and already holds
   the old values, and is then reconciled with IndexedDB.

   WHAT HAPPENS WHEN IndexedDB IS NOT THERE
   Private windows, blocked site data and some embedded viewers refuse it.
   Then this falls back to localStorage alone and the status line says so,
   rather than appearing to save and losing the lot on close.
   ============================================================================ */
(function (root) {
  'use strict';
  /** @type {any} the one namespace the page hangs everything on */
  const PB = root.PB = root.PB || {};

  const DB_NAME = 'plasmidbench';
  const DB_VERSION = 1;
  const RECORDS = 'records';         /* constructs, keyed by their own id */
  const SETTINGS = 'settings';       /* everything else, as key/value     */

  /* A value small enough to keep mirrored in localStorage as well, so that
     the next load has it before IndexedDB has opened. Sequences and stock
     lists are far bigger than this and live only in IndexedDB. */
  const MIRROR_LIMIT = 32 * 1024;

  let db = null;                     /* the open connection, or null */
  let available = null;              /* null = not yet decided */
  const mem = new Map();             /* the synchronous settings cache */
  /** @type {Promise<any>} */
  let queue = Promise.resolve();     /* writes, one at a time, in order */

  /* ---------------------------------------------------------------- open -- */

  function openDB() {
    return new Promise((resolve, reject) => {
      if (typeof indexedDB === 'undefined' || !indexedDB) { reject(new Error('no IndexedDB')); return; }
      let req;
      try { req = indexedDB.open(DB_NAME, DB_VERSION); }
      catch (e) { reject(e); return; }
      req.onupgradeneeded = ev => {
        const d = req.result;
        if (!d.objectStoreNames.contains(RECORDS)) d.createObjectStore(RECORDS, { keyPath: 'id' });
        if (!d.objectStoreNames.contains(SETTINGS)) d.createObjectStore(SETTINGS, { keyPath: 'key' });
        void ev;
      };
      req.onsuccess = () => resolve(req.result);
      req.onerror = () => reject(req.error || new Error('IndexedDB refused to open'));
      req.onblocked = () => reject(new Error('IndexedDB is blocked by another tab'));
    });
  }

  function tx(storeName, mode) {
    if (!db) throw new Error('store is not open');
    return db.transaction(storeName, mode).objectStore(storeName);
  }

  function wrap(request) {
    return new Promise((resolve, reject) => {
      request.onsuccess = () => resolve(request.result);
      request.onerror = () => reject(request.error);
    });
  }

  /* ------------------------------------------------------------ settings -- */

  /* The keys localStorage used, with their `pb.` prefix. They are read once
     at start-up so that a library saved by an earlier version is not lost,
     and left where they are afterwards: deleting them would strand anyone
     who opens an older copy of the page. */
  function seedFromLocalStorage() {
    const found = [];
    try {
      for (let i = 0; i < localStorage.length; i++) {
        const k = localStorage.key(i);
        if (!k || k.indexOf('pb.') !== 0) continue;
        const key = k.slice(3);
        try { mem.set(key, JSON.parse(localStorage.getItem(k))); found.push(key); }
        catch (e) { /* a value that is not JSON is not ours */ }
      }
    } catch (e) { /* localStorage itself refused; nothing to seed from */ }
    return found;
  }

  function get(key, dflt) {
    if (!mem.has(key)) return dflt;
    const v = mem.get(key);
    return v === undefined || v === null ? dflt : v;
  }

  function set(key, value) {
    mem.set(key, value);
    /* small values stay in localStorage too, so that the next load has them
       before any transaction has completed */
    try {
      const s = JSON.stringify(value);
      if (s.length <= MIRROR_LIMIT) localStorage.setItem('pb.' + key, s);
      else localStorage.removeItem('pb.' + key);
    } catch (e) { /* full or refused: IndexedDB is the real store */ }
    if (!db) return Promise.resolve(false);
    queue = queue.then(() => wrap(tx(SETTINGS, 'readwrite').put({ key: key, value: value })))
                 .then(() => true, () => false);
    return queue;
  }

  /* ------------------------------------------------------------- records -- */

  function allRecords() {
    if (!db) return Promise.resolve([]);
    return wrap(tx(RECORDS, 'readonly').getAll()).catch(() => []);
  }

  function putRecord(rec) {
    if (!db || !rec || !rec.id) return Promise.resolve(false);
    /* stored as a plain object: an IndexedDB value has to be structured
       cloneable, and a record carries nothing else, but a stray function or
       DOM node on it would otherwise fail the whole transaction */
    const flat = JSON.parse(JSON.stringify(rec));
    queue = queue.then(() => wrap(tx(RECORDS, 'readwrite').put(flat)))
                 .then(() => true, () => false);
    return queue;
  }

  function deleteRecord(id) {
    if (!db || !id) return Promise.resolve(false);
    queue = queue.then(() => wrap(tx(RECORDS, 'readwrite').delete(id)))
                 .then(() => true, () => false);
    return queue;
  }

  /* ---------------------------------------------------------------- init -- */

  /* Opens the database, seeds the settings cache, and brings across anything
     an earlier version left in localStorage. Resolves to a short report so
     the caller can say what happened rather than guess. */
  async function init() {
    const seeded = seedFromLocalStorage();
    let migratedRecords = 0, migratedSettings = 0;
    try {
      db = await openDB();
      available = true;
    } catch (e) {
      available = false;
      return { available: false, why: String(e && e.message || e), seeded: seeded.length,
               migratedRecords: 0, migratedSettings: 0 };
    }

    /* settings already in IndexedDB win over the localStorage copy: they are
       the ones that were written last, and the mirror may have been skipped
       for size */
    try {
      const rows = await wrap(tx(SETTINGS, 'readonly').getAll());
      rows.forEach(r => { if (r && r.key) mem.set(r.key, r.value); });
      /* and anything only localStorage had is written across, once */
      const have = new Set(rows.map(r => r.key));
      for (const key of seeded) {
        if (have.has(key)) continue;
        await wrap(tx(SETTINGS, 'readwrite').put({ key: key, value: mem.get(key) }));
        migratedSettings++;
      }
    } catch (e) { /* settings unreadable; the in-memory seed still serves */ }

    /* the construct library used to live under one localStorage key as an
       array; each construct becomes a row of its own here */
    try {
      const existing = await wrap(tx(RECORDS, 'readonly').count());
      const old = mem.get('library.v1');
      if (!existing && Array.isArray(old) && old.length) {
        for (const rec of old) {
          if (!rec || !rec.id) continue;
          await wrap(tx(RECORDS, 'readwrite').put(JSON.parse(JSON.stringify(rec))));
          migratedRecords++;
        }
      }
    } catch (e) { /* leave the old key alone; it is still readable */ }

    return { available: true, seeded: seeded.length,
             migratedRecords: migratedRecords, migratedSettings: migratedSettings };
  }

  /* How much room is left, where the browser will say. Shown in the storage
     panel so that running out is something a person sees coming rather than
     meets as a failed save. */
  async function quota() {
    try {
      if (!navigator.storage || !navigator.storage.estimate) return null;
      const e = await navigator.storage.estimate();
      return { used: e.usage || 0, available: e.quota || 0 };
    } catch (err) { return null; }
  }

  PB.store = {
    init: init,
    get: get,
    set: set,
    allRecords: allRecords,
    putRecord: putRecord,
    deleteRecord: deleteRecord,
    quota: quota,
    isAvailable: () => available === true,
    /* for tests: forget everything and start again */
    _reset: async function () {
      mem.clear();
      if (db) { db.close(); db = null; }
      available = null;
      await new Promise(res => {
        const r = indexedDB.deleteDatabase(DB_NAME);
        r.onsuccess = r.onerror = r.onblocked = () => res();
      });
    }
  };

  if (typeof module !== 'undefined' && module.exports) module.exports = PB.store;
})(typeof window !== 'undefined' ? window : globalThis);
