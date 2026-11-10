/* ============================================================================
   The shapes that cross a module boundary.

   This is the only TypeScript file in the project and it holds no code — just
   the types, so that `npx tsc -p tsconfig.json` can hold the engine to them:
   a function called with the wrong arguments, a property read that does not
   exist on what a function returns, a value that may be null used as though
   it cannot. Nothing is compiled and nothing in src/ changes shape; the
   engine stays plain JavaScript that Node can require as it is. Section 4.5
   of the evaluation argues for this over a TypeScript migration.

   It is a .d.ts rather than JSDoc in one of the engine files because each of
   those is a CommonJS module, and a typedef written inside a module is
   visible only there. Declared here, these names are available to every
   file the project checks.

   Two conventions to know before reading:

   * A position is 0-based and a range is end-exclusive, everywhere, in every
     module. GenBank's 1-based inclusive numbers exist only in the text a
     reader consumes and a writer produces.
   * Several functions return an array that also carries named properties —
     E.pcr returns the products with `found` and `sites` on the side. That is
     written out here as an intersection rather than left to be discovered.
   ============================================================================ */

/**
 * One annotated feature. `parts` is present when the location is a join;
 * `mixed` when the parts are not all on the same strand. A feature whose
 * `end` is at or before its `start`, or past the sequence length, crosses
 * the origin of a circular molecule — both forms occur in real files, and
 * Bio.featurePieces resolves either.
 */
declare interface Feature {
  start: number;
  end: number;
  strand: number;
  type: string;
  /** what to call it; kept whole, shortened only where it is drawn */
  label?: string;
  quals?: { [key: string]: any };
  parts?: Array<{ start: number; end: number; strand?: number }>;
  mixed?: boolean;
  /** a site between two bases, written a^b */
  between?: boolean;
  color?: string;
  pinned?: boolean;
  [extra: string]: any;
}

/**
 * A construct. `molecule` is 'protein' for a record read from a protein
 * file, in which case `seq` holds residues and no DNA operation applies.
 * `topologyFromFeature` marks a record whose LOCUS line claimed linear
 * while a feature crossed the origin, so the topology was corrected.
 */
declare interface SeqRecord {
  name: string;
  seq: string;
  topology: string;
  features: Feature[];
  description?: string;
  molecule?: string;
  dropped?: number;
  droppedChars?: string;
  skippedFeatures?: string[];
  topologyFromFeature?: boolean;
  contigOnly?: boolean;
  updatedAt?: number;
  createdAt?: number;
  id?: string;
  builtin?: boolean;
  [extra: string]: any;
}

/**
 * A restriction enzyme. One that cuts inside its own site carries top/bot,
 * each the distance from the start of the site to the break in that strand.
 * A type IIS enzyme cuts beyond its site and carries o1/o2 instead,
 * measured from the end of it.
 */
declare interface Enzyme {
  name: string;
  /** IUPAC codes */
  site: string;
  top?: number;
  bot?: number;
  o1?: number;
  o2?: number;
  typeIIS?: boolean;
  /** REBASE supplier letters */
  suppliers?: string;
  alias?: string;
}

declare interface Site {
  pos: number;
  cutTop: number;
  cutBottom: number;
  strand: number;
  [extra: string]: any;
}

declare interface Product {
  seq: string;
  length?: number;
  [extra: string]: any;
}

/** The products E.pcr found, with the search's own findings on the array. */
declare type PcrResult = Product[] & {
  found?: number;
  sites?: { forward: number; reverse: number };
  truncated?: boolean;
};

/**
 * The features a table reader found, with the locations it refused attached
 * to the array. Refusing a location that cannot be placed on this sequence
 * — a part in another accession, a negative coordinate — and saying which,
 * is not the same as losing it.
 */
declare type FeatureList = Feature[] & { skipped?: string[] };

declare interface SeqTextReport {
  /** DNA-cleaned, U rewritten as T */
  seq: string;
  /** every letter, kept as written: what a protein record needs */
  protein: string;
  dropped: number;
  droppedChars: string;
  hadU: boolean;
  looksProtein: boolean;
}

/**
 * One entry of a feature library. The provenance fields are not optional
 * decoration: a library assembled from several sources has different licence
 * terms and different reliability per entry, and a feature called in a
 * plasmid map has to be traceable to the thing that called it.
 */
declare interface LibraryEntry {
  name: string;
  seq: string;
  type?: string;
  /** where this entry came from — ApE, FPbase, a file someone wrote */
  source?: string | null;
  licence?: string | null;
  accession?: string | null;
  note?: string | null;
  /** other names for the same sequence, folded together on import */
  aliases?: string[];
}

/**
 * An Error with a short machine-readable reason on it. The readers throw
 * these so a caller can tell "that is not JSON" from "that JSON is not a
 * sequence" without matching on the message text.
 */
declare interface CodedError extends Error {
  code?: string;
}
