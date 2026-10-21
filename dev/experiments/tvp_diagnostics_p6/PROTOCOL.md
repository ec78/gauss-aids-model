# TVP-AIDS D1-P6 private anchor continuation

Status: **FROZEN before measurement.**

This private, dev-only continuation applies the accepted D1 tiered diagnostic
unchanged to the one remaining prespecified anchor. It changes no shipped
source, public API, documentation, example, package metadata or version.

## Fixed fixture

- `quaidsExampleData(500, 204)`, first four goods (`n1=3`), with shares
  renormalized after subsetting;
- homogeneous linear `quaidsFit()` static fit and diagonal `H` from
  `homogSse[1:3]/nobs`, matching example 14 and the prior P6 studies;
- `q0 = 0.1*mean(diag(H))*ones(12,1)`;
- unchanged D1 starts `0.1 | 1 | 10`, `maxIters=2000`, `maxTime=600`;
- common total filter likelihood and D1 tolerance `1e-6*nobs = 5e-4`.

The runner executes the job twice sequentially. Each job has a 3600-second
outer timeout with process-tree cleanup. Tagged payloads exclude timing and
optimizer prose, trim line-ending whitespace, use UTF-8 without a BOM, and
must be byte-identical. Every fitted Q element is emitted on its own tagged
line so GAUSS line wrapping cannot hide data.

## Acceptance

Both runs must exit zero and report the complete payload, static-fit
convergence, retcodes `0|0|0`, `START-STABLE-OBSERVED`, best start 0.1,
collapse counts `12|12|12`, and a maximum LL gap no larger than `5e-4`.
As a reproduction check, likelihoods must match the prior frozen convergence
study within `1e-9` and Q sums within `1e-12`, in D1 start order:

| multiplier | prior total LL | prior Q sum |
|---:|---:|---:|
| 0.1 | -6841.8654840494210 | 0.0013269734775710979 |
| 1 | -6841.8655741471166 | 0.0013406784859562100 |
| 10 | -6841.8654911245685 | 0.0013243386530964532 |

`D1_BASELINE.sha256` protects every accepted D1 file byte-for-byte.
`FREEZE.sha256` protects this protocol, the baseline manifest, runner and
probe using LF-normalized hashes. The runner verifies both manifests, D1's
own normalized freeze and the strict SSLib 2.1.2 pin before execution, then
rechecks the manifests and pin afterward.

This is a bounded reproducibility and observed-start-agreement check. It is
not evidence of global optimality, parameter identification, boundary
inference or valid Q uncertainty estimates. Timings and iteration counts are
not exposed by the unchanged D1 result schema and are deliberately excluded.
