# printQuaidsPreflight

## Purpose

Prints a compact console report from a
[quaidsPreflight](quaidsPreflight.md) result.

## Format

```gauss
call printQuaidsPreflight(pOut);
```

## Parameters

- `pOut` - a `quaidsPreflightOut` structure returned by
  [quaidsPreflight](quaidsPreflight.md).

## Returns

Nothing. The procedure prints data, design, IV, cluster, and convergence-risk
diagnostics to the console.

## Remarks

Separated from [quaidsPreflight](quaidsPreflight.md) so preflight checks can
be used silently in scripts, tests, and higher-level workflows.

## Examples

```gauss
pOut = quaidsPreflight(w, intcpt, prices, totexp, instr, aCtl, 0);

if not pOut.ok;
    call printQuaidsPreflight(pOut);
endif;
```

## Source

`quaidsdiagnostics.src`

## See Also

[quaidsPreflight](quaidsPreflight.md), [quaidsFit](quaidsFit.md)
