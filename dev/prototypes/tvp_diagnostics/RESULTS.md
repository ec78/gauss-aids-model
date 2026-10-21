# TVP-AIDS Phase D1 private tiered diagnostic prototype: initial results

The prototype is private and dev-only. It changes no shipped source, public
API, documentation, package metadata, version or release behavior.

## Result

The initial gate passes.

- The cheap tier reports solver status, Q validity and caller-start-relative
  scale collapse from one `quaidsTVPOut` without calling the optimizer again.
- The opt-in tier runs fixed 0.1x, 1x and 10x starts and compares total filter
  likelihoods at tolerance `1e-6*nobs`. Q-vector equality is descriptive only.
- Nine pure-logic GAUSS checks pass. They cover status precedence, invalid Q,
  the strict scale marker (`< 1e-3`, not `<=`), the inclusive LL tolerance,
  a gap above tolerance, a nonzero retcode and deterministic tie selection.
- On the Stage-6 n1=2/tobs=120/seed=7 fixture, the 1x arm is bit-exact with
  `quaidsTVPFit()`: Qfit difference 0, retcode equal and independently
  re-evaluated total-LL difference 0. The three-start maximum LL gap is
  `1.0573e-11`, below tolerance `1.2e-4`.
- On the frozen representative n1=2/tobs=200/seed=1001/oracle-H probe, all
  three fits return rc 0. Maximum LL gap is `5.0022e-12`, below tolerance
  `2e-4`, so the multi-start result is `START-STABLE-OBSERVED`. At the same
  time the three arms report 1, 1 and 2 start-relative scale-collapse
  elements. This demonstrates the intended separation: optimizer agreement
  does not erase a scale warning.
- The tagged probe output reproduced exactly across two consecutive runs.

## Interpretation

The interface separation is feasible and the first live case behaves as D0
predicted. The single-fit tier is cheap but intentionally weak; its q0-relative
marker cannot be called a boundary or identification diagnostic. The
multi-start tier is materially more expensive and supports only observed
agreement among the starts tried, not global optimality.

Hessian, profile, score-threshold, fixed-point-residual and inference features
remain excluded. The optional P6 anchor probe has not been run in this initial
gate; it is the next bounded validation if the owner wants D1 extended.

Two implementation boundaries remain before any public API proposal:
`scaleRatio` needs explicit input validation, and the multi-start result needs
an explicit invalid-Q/invalid-LL state rather than relying only on fit status.
Neither affects the frozen initial measurements above.

## Reproduction

From this directory, under the repository's pinned SSLib installation:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_d1.ps1 -Equivalence -ProbeN12
```

The accepted measurement boundary is `4a7a939`; it retains the original
prototype freeze at `419319c` and includes the exact probe job in
`FREEZE.sha256`. The tagged
equivalence output is preserved in `results/equiv.txt`. Both exact
probe payloads and their matching SHA-256 hashes are preserved in
`results/probe_n12_run1.txt`, `results/probe_n12_run2.txt` and
`results/probe_n12.sha256`. The representative job is
`results/gen/probe_n12.e`. The frozen configuration and amendments are in
`PROTOCOL.md`.
