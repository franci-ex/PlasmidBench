# PlasmidBench and eLabFTW

eLabFTW is the lab's record system. PlasmidBench is not one, and should not
become one. This document says how the two work together, and answers the
question that decides the architecture.

## The division of labour

| system | responsibility |
|---|---|
| **eLabFTW** | the records: plasmids, primers, samples, inventory, experiments, and the audit trail that makes them GLP records at all |
| **OpenCloning** | reproducible cloning design and construct provenance |
| **PlasmidBench** | portable analysis, design and visualisation, offline, aware of what the freezer holds |

## The rule this follows

> **Read freely. Write only by adding, and only through the API.**

Good Laboratory Practice expects records that are attributable, legible,
contemporaneous, original and accurate. eLabFTW provides that. A design tool
that held authoritative records in browser storage and synchronised later
would break both attribution and contemporaneity, so this one does not try.

In practice, and enforced in `src/elab.js`:

- Nothing ever replaces an attachment, patches a body, or deletes anything.
- A result goes back as a **new attachment** plus a **comment** saying what
  made it, when, and with the sequence checksum.
- The original stays exactly as it was. The addition is attributable to the
  person whose key was used and timestamped by eLabFTW's own audit trail.
- A construct is proposed, not asserted: `propose()` creates a resource
  titled `[proposed] …`, whose body says what built it and that it should be
  sequenced before being treated as stock.

`tests/uitest-elab.js` checks this against the wire rather than against
intent: it asserts that no `PATCH`, no `DELETE` and no POST to an existing
upload ever leaves the page.

---

## The question that decides everything: CORS

Section 7.3 of the evaluation is right that this comes first. Two routes are
possible and the server decides which:

- **Option A — file exchange.** Download from eLabFTW, work here, attach the
  result there. Works today, needs no code and no permission.
- **Option B — the browser calls the API directly**, each user with their own
  personal API token held in their own browser.

Option B turns entirely on whether eLabFTW's CORS configuration permits
browser-origin requests. And there is a subtlety that catches people: the API
key goes in an `Authorization` header, which is **not** a CORS-safelisted
header, so every request is preceded by an `OPTIONS` preflight. A server that
allows a plain `GET` but does not answer the preflight naming `Authorization`
will still fail — and will fail in a way the browser reports to the page as
an ordinary network error, with no detail, by design.

### How to get the answer in two seconds

**Tools → eLabFTW → Test the connection.** It runs five steps and reports
each separately, because the useful answer is which one failed:

| step | what it establishes |
|---|---|
| 1. Something answers | sent deliberately **without** the key, and not preflighted, so "the host is unreachable" is separated from "the preflight is refused". A 401 here is a success: something answered. |
| 2. A request carrying the API key | the preflight. **This is the step a default install fails.** |
| 3. Whose key it is | attribution, read from eLabFTW. Nothing here invents an author. |
| 4. Reading resources | whether records can be listed |
| 5. Permission to add | whether a read-write key is present |

The verdicts it gives are `direct browser access works`, `blocked by CORS`,
`the key was rejected`, `nothing answered`, and `the server refused` — and
they are distinguished, because "blocked by CORS" and "the key was rejected"
have completely different remedies and look identical from inside a browser
if you only try one request.

### If it says "blocked by CORS"

The server needs to send, for this page's origin:

```
Access-Control-Allow-Origin: https://your-page-origin
Access-Control-Allow-Headers: Authorization, Content-Type
Access-Control-Allow-Methods: GET, POST, PATCH, OPTIONS
```

and answer `OPTIONS` with a 204. In an nginx front end that is a few lines.
Note that `Access-Control-Allow-Origin: *` is enough here because no cookies
are sent — the key is the whole of the authentication and requests are made
with `credentials: 'omit'`, which deliberately keeps the server's CORS
requirements as loose as possible.

Until that is changed, **use the file route**. The panel describes it, it
needs nothing installed, and the audit trail is identical: the file arrives
through eLabFTW's own upload, attributed to you and timestamped by it. It is
a few more clicks and it is the honest first proof of concept.

Only if both routes fail is a server-side component justified, and then as
the smallest possible proxy rather than a synchronisation service.

---

## The API key

- It is **the person's own**, generated in eLabFTW → Settings → API keys.
- It is **never in the distributed file.** There is nowhere in the source for
  one to live.
- It is typed in and kept **in that browser**, or **for the session only**,
  as they choose. "Forget the key" clears both.
- A **read-only** key is enough for everything except sending results back,
  and the connection test says which kind has been given.
- The lab profile **cannot** export it: that is enforced by pattern in
  `src/profile.js` and checked by a test that plants a key and searches the
  exported file for it.

---

## What travels, and how it is identified

Per section 7.4 of the evaluation:

- **GenBank** for annotated sequences, FASTA for plain sequence.
- **A sequence checksum** on everything: `sha256:` over the topology and the
  cleaned sequence, so two systems can agree they hold the same molecule
  without passing it back and forth. SHA-256 is implemented in the engine
  rather than taken from `crypto.subtle`, because `crypto.subtle` does not
  exist on a `file://` page — which is exactly how a lab machine opens this.
  It is checked against Python's `hashlib`.
- **A download is verified** against the checksum eLabFTW stored for it. A
  file whose bytes do not match is **refused, not opened**.
- **The eLabFTW identifier is recorded on the record** when a file is opened
  from it, so the link back is not only in somebody's memory.

---

## Trying it without a server

`src/../mock-elab.js` is a small stand-in for the API, used by the test
suite. It reproduces the behaviours that matter:

```sh
node mock-elab.js 8099 full      # correctly configured
node mock-elab.js 8099 simple    # a default install: preflight refused
node mock-elab.js 8099 none      # no CORS headers at all
```

It prints two keys, one read-write and one read-only. Point the panel at
`http://127.0.0.1:8099` and the connection test will tell you, correctly,
which case you are in. That is how the integration was developed without
touching a real instance, and it is how a change to it can be checked before
it goes anywhere near one.
