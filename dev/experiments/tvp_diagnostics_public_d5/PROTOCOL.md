# D5 public TVP diagnostic evidence protocol

This bounded release-candidate check exercises the shipped
`quaidsTVPMultiStartDiagnostics()` API on the frozen D1/P6 anchor. It does
not expand the convergence study, tune the estimator, or claim global
optimality or identification.

## Frozen design

- Source candidate: v0.5.0 integration commit `600ec16` or its D5 evidence
  descendant, with no estimator-source changes.
- Dependency: the repository's strict SSLib 2.1.2 pin must pass before and
  after execution.
- Data: `quaidsExampleData(500, 204)`, first four goods, with the retained
  shares renormalized to sum to one.
- Static seed: homogeneous linear QUAIDS fit; `H` is diagonal with entries
  `homogSse[1:3]/nobs`.
- TVP seed: `q0 = 0.1 * mean(diag(H)) * ones(12, 1)`.
- Optimizer caps: `maxIters=2000`, `maxTime=600` minutes per arm (the
  `quaidsTVPControl.maxTime`/CMLMT unit). The runner imposes the effective
  bounded wall-clock limit: 2,400 seconds for the complete GAUSS process.
- Starts: the public API's fixed `0.1|1|10` multipliers.
- One public-API execution. The frozen private gate already ran this anchor
  twice with byte-identical payloads; D5 is a bounded parity check, not a
  new determinism study.

## Prespecified acceptance

All three fits must return `rc=0`, valid Q, and valid likelihoods. The best
start must be index 1 (multiplier 0.1), all three collapse counts must be
12, and the maximum likelihood gap must not exceed the public tolerance
`5e-4`. Likelihoods must agree with the frozen private values within
`1e-9`, and Q sums within `1e-12`:

```
LL:    -6841.8654840494210  -6841.8655741471166  -6841.8654911245685
Q sum:  0.0013269734775710979 0.0013406784859562100 0.0013243386530964532
```

The public interpretation must be status 5, label `SCALE-CAUTION`, with
`startStable=1` and `anyScaleCaution=1`. Start agreement and scale caution
are intentionally simultaneous.

The runner verifies that the frozen private D1/P6 payload remains at SHA-256
`b33836572df51a196449208059026b6f26eaa295cdd74bc2742af9a8761090d4`
and that the tracked private evidence has no diff before and after the run.

