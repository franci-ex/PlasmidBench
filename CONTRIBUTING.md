# Contributing

The project is a lab tool first. The thing that makes it usable by anyone
other than its author is that every claim it makes can be checked, so that
is what contributions are held to.

Read `AGENTS.md` before changing anything. It has the conventions and the one
rule about tests. `docs/developer-guide.md` has the practical steps.

## Before opening a pull request

```sh
./validate-all.sh
```

All four layers green, or say in the pull request which is not and why.
`npx tsc -p tsconfig.json` must report zero errors.

## What a pull request should say

1. **What changed, and why.**
2. **What it was validated against** — name the oracle, the library, the
   paper, or the arithmetic. "Tests pass" is not an answer; the tests are
   part of the change.
3. **Anything measured, with the numbers.** If you made something faster or
   more sensitive, say how much. If you did not measure it, say that.
4. **Whether a result changed.** Someone designed a plasmid with the old
   number. That belongs at the top of the description, not in a footnote.

## One module per pull request

Small enough that someone can actually read it. This matters more than it
used to: it is easy to produce a thousand plausible lines, and a reviewer
cannot validate what they cannot read.

## The rule about tests

**A test and the code it tests must not come from the same place.** Expected
values come from an outside source — Biopython, REBASE, a published formula,
a vendor's spec sheet, or arithmetic written out in the test.

**Never change a test expectation to make it pass.** Find out which side is
wrong. If it is the test, say in a comment what settled it. There are a dozen
such comments in `tests/` and each one is evidence; please add to them rather
than tidying them away.

## AI assistance

It is welcome and it is part of this project's history — `NOTICE.md` records
that. Two conditions:

- **Say so in the pull request.** Not as a disclaimer; as information a
  reviewer should have.
- **Be able to explain every line.** An unexplainable nearest-neighbour
  calculation cannot be validated by anyone, which makes it worthless in a
  lab whether or not it happens to be right. If you cannot explain it, it is
  not ready.

## Reporting a problem

An issue is most useful with:

- the **file** that caused it, if a file did — a malformed GenBank is a gift,
  not a nuisance, and three real bugs came from exactly that;
- what you expected, and what happened;
- which build (the version in the page footer, or the commit).

If something gave a **wrong scientific answer**, say so plainly and label the
issue `correctness`. Those come before features, always.

## Things deliberately not done

Before proposing one of these, see the reasoning in `AGENTS.md`:

no crawling Addgene · no hand-built feature library · no TypeScript
migration · no new cloning simulators · no authoritative records in this tool
· no credentials in the built file.
