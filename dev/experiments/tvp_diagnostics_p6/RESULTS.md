# TVP-AIDS D1-P6 private anchor results

The frozen bounded gate passes. This experiment is private and dev-only; it
changes no shipped source, public API, documentation, examples, package
metadata, version or release behavior.

## Result

- The accepted D1 unit suite passed 9/9 before the anchor runs.
- The strict installed SSLib 2.1.2 pin passed before and after execution.
- Both complete anchor runs exited successfully and produced byte-identical
  tagged payloads with SHA-256
  `b33836572df51a196449208059026b6f26eaa295cdd74bc2742af9a8761090d4`.
- The static homogeneous linear fit converged. All three TVP fits returned
  `rc=0`.
- In fixed D1 start order 0.1, 1 and 10, total filter likelihoods were
  `-6841.8654840494210`, `-6841.8655741471166` and
  `-6841.8654911245685`. The maximum gap was `9.00977e-5`, below the
  `5e-4` tolerance, so the result is `START-STABLE-OBSERVED`; the 0.1 arm
  has the best observed likelihood.
- Q sums were `0.0013269734775710979`, `0.0013406784859562100` and
  `0.0013243386530964532`. Both likelihoods and Q sums reproduce the prior
  frozen convergence study within their prespecified tolerances.
- All 12 Q elements in every arm are below `1e-3` times that arm's own
  starting value: collapse counts are `12|12|12`.

## Interpretation

The expensive multi-start tier is reproducible on the P6 anchor and gives
observed likelihood agreement across the three tested starts. At the same
time, the cheap start-relative marker raises the strongest possible scale
caution. These statements coexist: start agreement does not establish
identification or remove the scale warning.

This result does not establish global optimality, a statistical boundary,
Q identification, Q standard errors or inferential validity. Timing and
iteration counts are not exposed by the unchanged D1 result schema and were
not retained as evidence.

## Evidence

- Frozen protocol and harness: `PROTOCOL.md`, `FREEZE.sha256`
- Accepted D1 byte baseline: `D1_BASELINE.sha256`
- Exact payloads: `results/probe_p6_run1.txt`,
  `results/probe_p6_run2.txt`
- Portable payload hashes: `results/probe_p6.sha256`

Reproduce from the repository root with:

```powershell
powershell -ExecutionPolicy Bypass -File .\dev\experiments\tvp_diagnostics_p6\run_d1_p6.ps1 -TimeoutSeconds 3600
```
