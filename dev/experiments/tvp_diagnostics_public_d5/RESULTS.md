# D5 public TVP diagnostic evidence

The frozen bounded check passed on GAUSS 26.1.4 with the strict SSLib 2.1.2
pin. It changes no estimator source, public API, documentation, example, or
package metadata.

## Result

- The static homogeneous linear fit converged.
- All three public multi-start fits returned `rc=0`, with valid Q and valid
  common-scale likelihoods.
- For starts `0.1|1|10`, likelihoods were
  `-6841.8654840494210|-6841.8655741471166|-6841.8654911245685` and Q sums
  were `0.0013269734775710979|0.0013406784859562100|0.0013243386530964532`.
  They match the frozen private D1/P6 evidence within the prespecified
  `1e-9` and `1e-12` tolerances.
- The 0.1 start was best. The maximum likelihood gap was `9.00977e-5`, below
  the `5e-4` tolerance, so `startStable=1`.
- Every arm had 12 start-relative scale collapses. The public result therefore
  reports status 5, `SCALE-CAUTION`, with `anyScaleCaution=1` while retaining
  the separate observed start-stability result.
- Strict SSLib 2.1.2 pin and private-evidence integrity checks passed before
  and after execution. The frozen private payload remains unchanged at its
  portable SHA-256 `b3383657...0d4`.

## Interpretation

This confirms that the shipped public result reproduces the frozen private
anchor and communicates the intended precedence: scale caution is the overall
status even when the three tested starts agree. It is descriptive evidence,
not a claim of global optimality, Q identification, boundary-aware inference,
or general convergence across untested data and starts.

## Evidence

- Frozen design and acceptance: `PROTOCOL.md`, `FREEZE.sha256`
- Exact tagged payload: `results/public_p6.txt`
- Payload hash: `results/public_p6.sha256`

Reproduce from the repository root with:

```powershell
powershell -ExecutionPolicy Bypass -File .\dev\experiments\tvp_diagnostics_public_d5\run_public_p6.ps1
```
