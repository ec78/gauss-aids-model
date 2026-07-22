# getDefaultQuaidsControl

## Purpose

Alias for [quaidsControlCreate](quaidsControlCreate.md), named to match the
`getDefault...Control` naming convention used in the `gauss-qardl` sibling
library.

## Format

```gauss
struct quaidsControl aCtl;
aCtl = getDefaultQuaidsControl();
```

## Parameters

None.

## Returns

Identical to [quaidsControlCreate](quaidsControlCreate.md) -- see that page
for the full field list and defaults.

## Remarks

Purely a naming-convention alias; both procs return the same defaults.
Prefer whichever name reads better in your own code.

## Examples

```gauss
struct quaidsControl aCtl;
aCtl = getDefaultQuaidsControl();
```

## Source

`quaidsutil.src`

## See Also

[quaidsControlCreate](quaidsControlCreate.md)
