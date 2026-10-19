# TVP-AIDS convergence and speed study protocol

**Status: DRAFT v0 -- not approved or frozen; do not run the study grid.**

This post-0.4.0 study addresses the next TVP boundary without changing
shipped source or public API. It first tests whether the existing optimizer
reaches the best available endpoint on known-truth data. Only after those
reference endpoints are frozen may a separate speed phase compare an
implementation change.

## Claims and non-claims

Phase C asks whether the shipped optimizer, unchanged, reaches the best
log-likelihood found across prespecified starts within a fixed iteration
budget. Results apply only to the prespecified DGPs and P6 anchor. They are
not a convergence guarantee or a general scale claim.

Phase S asks whether a candidate optimization change reduces runtime while
reaching an equivalent endpoint on every Phase C reference cell. No speed
claim is allowed before endpoint equivalence passes. `mleRetcode == 0` alone
is never treated as proof of a maximum.

## Phase C questions

1. What share of known-truth fits return `mleRetcode == 0`, by `n1`, sample
   length, starting value and fixed-H rule?
2. Among those fits, which agree with the best log-likelihood found across
   the prespecified starts for the same dataset?
3. Does each returned endpoint satisfy the necessary check
   `LL(Qfit) >= LL(Qtrue) - tolerance`?
4. How often are fitted Q elements near the positivity boundary? This is
   descriptive, not a per-element accuracy gate.
5. Does the P6 `n1 = 3` example-data anchor converge under a larger,
   iteration-bound budget, and do its starts agree?
6. Are iteration-capped results deterministic?

## Proposed fixed design -- owner approval required before freeze

| Factor | Proposed levels |
|---|---|
| `n1` | 2, 3, 4 (`k_states` 7, 12, 18); `n1 >= 5` deferred |
| `tobs` | 200 for all; 500 additionally for `n1 = 2` |
| study seeds | `n1=2`: 1001-1008; `n1=3`: 1001-1006; `n1=4`: 1001-1004 |
| reserve seeds | 1101-1110, used in order only when a study DGP has a share outside (0,1) |
| true Q | alpha block `1e-5`; gamma and beta blocks `1e-4` |
| true H | `2.5e-5 * eye(n1)` |
| H supplied to fit | oracle true H (primary); static LA-AIDS diagonal H (secondary, `n1=2/3`) |
| Q starts | shipped rule; 10 times shipped; 0.1 times shipped |
| caps | 500 iterations, 600 minutes so the iteration cap binds |
| P6 anchor | `n1=3`, three starts, 2000 iterations; `n1=4`, three starts, 500 iterations |
| total | approximately 129 sequential fits |
| hard budget | 6 hours; staged A (`n1=2`), B (`n1=3`), C (`n1=4`), anchor |

The exact `state0` vectors for `n1=3/4` must be written here and validated
only on reserve seeds before this document is frozen. The `n1=2` vector is
the existing Stage 6 fixture:

```text
.40 | .30 | .05 | -.02 | .04 | -.04 | .03
```

If Stage A projects the remaining work above six hours, reduce `n1=4` to
two seeds first, then `n1=3` to four seeds. No other factor changes.

## Per-fit record and categories

Each row records the cell identifiers; CMLMT retcode, iterations, final
log-likelihood and maximum transformed-gradient magnitude; log-likelihood
at true Q and at the best start; full Q in the job log plus boundary count;
wall time, CMLMT time and seconds per iteration; and warnings/errors.

Proposed endpoint tolerance: `1e-6 * nobs` in total log-likelihood units.
The Stage A pilot must first establish whether CMLMT's stored `Fct` is a
mean or sum; the frozen protocol will state the conversion explicitly.

- `CONV-GLOBAL`: rc 0, within tolerance of best found LL, and the true-Q
  necessary condition passes.
- `CONV-LOCAL`: rc 0 but worse than best found LL beyond tolerance.
- `CONV-FAILC3`: rc 0 but worse than true-Q LL beyond tolerance.
- `CAPPED`: rc 2.
- `OTHER-RC`, `ERROR`, or `TIMEOUT`: reported literally.

Q recovery summaries are descriptive. A single realization does not
identify each Q element accurately, especially near a variance boundary.

## Stop rules

Before every stage, require a clean recorded HEAD, strict SSLib pin, TVP
Stage 6 at 30/30, and the P5 likelihood cross-check at 12/12. The first
harness job must exactly reproduce `quaidsTVPFit()`'s `Qfit` and `bFinal`
on the Stage 6 fixture through the private component chain.

Stop the current stage on any GAUSS error, outer timeout or
`CONV-FAILC3`. Do not change settings and rerun. Any protocol or harness
change after freeze requires a documented amendment and a new hash. Phase C
contains no `src/` edits.

## Smallest implementation slice

After owner approval and protocol freeze:

1. Implement a Stage A-only harness under this directory, reusing the P6
   sequential-runner pattern and tagged output.
2. Run the harness-equivalence gate and the `n1=2` cells only.
3. Record partial results and re-project Stages B/C before continuing.

## Phase S entry gate

Profile the Phase C harness before altering code: count likelihood
evaluations per optimizer iteration and measure filter versus Q-update time.
Numerical gradients are a plausible cause of per-iteration cost growth, but
that is currently an inference, not evidence. An analytic-gradient candidate
must first match numerical gradients at prespecified Q points, then preserve
Phase C categories and log-likelihood endpoints before timing is compared.

## Boundary order after this study

1. Q identification and behavior at variance boundaries.
2. Boundary-aware Q inference and standard errors.
3. Fixed non-diagonal H estimation scope, then separately estimated H.
4. IV/endogeneity treatment of total expenditure after an explicit
   methodology decision and validated inference machinery.

