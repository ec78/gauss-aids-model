# quaidsTVPControlCreate

## Purpose

Returns a default `quaidsTVPControl` structure for
[quaidsTVPFit](quaidsTVPFit.md) -- mirrors
[quaidsControlCreate](quaidsControlCreate.md)'s own minimalism.

## Format

```gauss
tvpCtl = quaidsTVPControlCreate();
```

## Parameters

None.

## Returns

`tvpCtl` is a `quaidsTVPControl` structure:

- `smooth` - `1` by default: also run the RTS smoother (see
  [quaidsTVPFit](quaidsTVPFit.md)). Set to `0` to return only the
  filtered (real-time-causal) path, skipping the full-sample smoother.
- `othnam` - `0` by default: auto-generated `W1..Wn` good names. Set to
  a legacy character matrix of `n` names to override. Unlike
  `quaidsControl.othnam` (`string`-typed, defaults to `""`), this field is
  `matrix`-typed and defaults to `0` -- a `string`-typed field was found
  to reject a real character-matrix assignment (`error G0071 : Type
  mismatch`), so this field uses the same `matrix` type every output
  struct's own name field (`wnam`/`xnam`/etc.) already uses successfully.
  **A second gotcha**: a value built with `$|` (vertical string
  concatenation -- e.g. `"Food" $| "Housing" $| ...`) is still rejected
  even by this `matrix`-typed field; only the legacy character-matrix form
  `$+`/`ftocv()` produce is accepted. Coerce with a leading `0$+` first:
  `tvpCtl.othnam = 0$+someDollarPipeBuiltNames;`.

## Remarks

Deliberately minimal: `H` and `q0` (the observation covariance and
starting state-innovation-covariance diagonal) are NOT control-struct
fields -- they are required, data-dependent arguments to
[quaidsTVPFit](quaidsTVPFit.md) itself, not tunable defaults.

## Examples

```gauss
tvpCtl = quaidsTVPControlCreate();
tvpCtl.smooth = 0;   // filtered path only, skip the RTS smoother

tvOut = quaidsTVPFit(w, prices, totexp, H, q0, tvpCtl);
```

## Source

`quaidstvpfit.src`

## See Also

[quaidsTVPFit](quaidsTVPFit.md), [quaidsControlCreate](quaidsControlCreate.md)
