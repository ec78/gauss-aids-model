# quaidsFull

## Purpose

Estimates a demand system from a GAUSS dataframe by column name, instead of
assembling the `w`/`intcpt`/`prices`/`totexp`/`instr` matrices by hand.

## Format

```gauss
struct quaidsOut qOut;
qOut = quaidsFull(data, shareVars, priceVars, totexpVar, instrVars, extraVars, aCtl);
```

## Parameters

- `data` (*dataframe*) - as returned by `loadd()` or built with
  `asDF()`/`dfaddcol()`. `quaidsFull()` does not load files itself.
- `shareVars` (*Nx1 string array*) - budget-share column names, in the
  order the goods should appear.
- `priceVars` (*Nx1 string array*) - log-price column names, in the
  **same** order as `shareVars` -- `shareVars[i]` and `priceVars[i]` must
  refer to the same good. There is no name-matching magic; matching is by
  position only.
- `totexpVar` (*string*) - the log-total-expenditure column name.
- `instrVars` (*string array*) - instrument column name(s) for log total
  expenditure.
- `extraVars` (*string array, or scalar `0`*) - extra intercept-shifter
  column names (demographics etc.), or `0` for none -- matches
  `quaidsFit()`'s `intcpt == 0` convention. Checked via `type(extraVars)`,
  not value equality (comparing a string array to `0` with `==` is itself
  a type-mismatch trap in GAUSS).
- `aCtl` (*`quaidsControl` structure*).

## Returns

`qOut` is a `quaidsOut` structure, identical to calling
[quaidsFit](quaidsFit.md) directly on the matrices selected from `data`.

## Remarks

Not a `"y ~ x1 + x2"` formula string, deliberately: AIDS/QUAIDS is a
multi-equation system (N budget shares against N parallel log prices),
which doesn't fit GAUSS's single-equation formula grammar. Column-name
lists, matched positionally, are the natural fit instead.

Verified numerically identical to the equivalent `quaidsFit(matrices...)`
call, including the `extraVars == 0` path, by
`tests/quaids_formula_parity_test.e` (17 checks).

## Examples

```gauss
data = loadd("mydata.csv");

shareVars = "W1"$|"W2"$|"W3"$|"W4"$|"W5";
priceVars = "P1"$|"P2"$|"P3"$|"P4"$|"P5";

struct quaidsControl aCtl;
aCtl = quaidsControlCreate();

struct quaidsOut qOut;
qOut = quaidsFull(data, shareVars, priceVars, "TOTEXP", "Z1", "X1", aCtl);

// No extra intercept shifters:
qOut = quaidsFull(data, shareVars, priceVars, "TOTEXP", "Z1", 0, aCtl);
```

## Source

`quaidsformula.src`

## See Also

[quaidsFit](quaidsFit.md), [quaidsControlCreate](quaidsControlCreate.md)
