# Stage D execution record (2026-10-13)

- `attempt1_default_outroot_failure.log`: `run_sym_stage.ps1 -Stage dev`
  and `-Stage dev_det` with the default `-OutRoot` (tail of each run's
  output). Both failed during parameter binding, before any guard, solve or
  directory creation:

  ```
  Join-Path : Cannot bind argument to parameter 'Path' because it is an empty string.
  At ...\run_sym_stage.ps1:27 char:35
  +     [string]$OutRoot = (Join-Path $PSScriptRoot "symmetry"),
  ```

  Cause: under Windows PowerShell 5.1, `$PSScriptRoot` is empty in the
  parameter defaults of an advanced script (this one has
  `[Parameter(Mandatory)]`). Afterwards `symmetry/` held only the
  implementation stamp: no `dev/` or `dev_det/` directory and no seed output.
- `attempt2_explicit_outroot.log`: the same two stages with the same
  directory passed explicitly (`-OutRoot <repo>\dev\prototypes\aids_direct_lm\symmetry`).
  These produced `dev/` and `dev_det/`; `dev_out/` is the analysis.

The default was repaired after Stage D (the default now resolves after
parameter binding). The repair is in the execution wrapper only and is
numerically irrelevant: no Stage D output was rerun or changed.
