# Public TVP diagnostics D3 implementation protocol

Status: **D3 accepted; public integration remains D4.**

## Frozen decisions

The owner approved the D2 contract on 2026-10-22, including the policy that
valid converged likelihood agreement with any start-relative Q collapse
returns overall status `5` (`SCALE-CAUTION`) while `startStable` remains `1`.

D3 implements only the computational contract. It does not update package
version, release metadata, public API inventory, command documentation,
examples or release claims. Those changes remain one synchronized D4 gate.

## Implementation surface

- `src/quaids.sdf`: plain-matrix/string result structs
  `quaidsTVPDiagOut` and `quaidsTVPMultiStartOut`.
- `src/quaidstvpfit.src`: the cheap `quaidsTVPDiagnostics()` summary and the
  opt-in `quaidsTVPMultiStartDiagnostics()` fixed 0.1x/1x/10x fit, plus private
  validation, normalization, likelihood-comparison and status helpers.
- `tests/quaidstvpdiag_test.e`: pure contract checks and exact Stage-6 1x-arm
  parity.
- `tests/guard_error_cases/tvpdiag_*.e`: fail-loud ordering checks for each
  public tier, base q0 and derived-start scaling.
- `tests/run_source_tests.ps1`: SSLib-gated registration only.

The implementation validates all three derived starts before the first
optimizer call. Returned Q is normalized to a conformable all-missing column
before result assignment, and invalid likelihood is normalized to a missing
scalar before storage. Missing/NaN/Inf values therefore never reach division,
filtering, `minc`, `maxc` or `maxindc`.

## Compatibility boundary

- Existing `quaidsTVPFit()` behavior is unchanged.
- The existing optional TVP include chain is unchanged.
- `package.json` source/dependency lists are unchanged; SSLib is not made a
  core dependency.
- `dev/prototypes/tvp_diagnostics` and
  `dev/experiments/tvp_diagnostics_p6` are read-only frozen evidence.
- No inference, identification, boundary or global-optimum claim is added.

## D3 verification record

- SSLib strict pin: 2.1.2 at commit
  `a32885be37544c41e6379ffda1505f19b66855cb`.
- Focused diagnostic test: 27/27 passed.
- Four public validation guards: expected diagnostics observed.
- Source runner with bootstrap/pubtable/curvature/direct-evidence exclusions:
  passed, including Stage 6 (30/30), fixed-Q cross-check (12/12), all core
  tests and registered guards.
- D1/P6 files: no diff from implementation base `ae72c46`.
- Frozen P6 deterministic reproduction: a reduced 600-second trial timed out
  during the first probe without changing accepted evidence. The subsequent
  full-budget run passed both probes byte-for-byte with the original hash
  `b33836572df51a196449208059026b6f26eaa295cdd74bc2742af9a8761090d4`.
  All three arms returned rc 0; the maximum LL gap was `9.0098e-5` versus
  tolerance `5e-4`; collapse counts were 12/12/12; and prior LL/Q totals
  matched within the frozen tolerances.

## Remaining D3 gate

Accept D3 after `git diff --check` passes, the independent read-only review's
focused coverage findings are resolved, and the implementation/test changes
are committed together on the isolated branch.
