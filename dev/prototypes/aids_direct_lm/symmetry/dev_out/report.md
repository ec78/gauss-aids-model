# Symmetry protocol v1 -- dev analysis

- block: dev; seeds: 30
- verdict: None (informational: REJECT)

## Integrity

- no problem found

## Runtime (never a decision input)

```
{
 "exceeds_cap": false,
 "max_seed_s": 43.332999999998144,
 "median_seed_s": 13.66350000000093,
 "p90_seed_s": 24.298799999998423,
 "projected_stage_c_s": 2322.795000000158,
 "slow_arms": []
}
```

## Gates

| Gate | Statistic | Threshold | Pass |
|---|---|---|---|
| C0 | 0 | no integrity problem | True |
| C1 | 30 | 161 | False |
| C2 | 0 | 1 | True |
| C3 | 0 | 1 | True |
| C4 | 0 | 3 | True |
| C5 | 0 | 3 | True |
| G6 | None | all G6 conditions | None |

## Descriptive (never pass/fail)

```
{
 "W_H": 0,
 "W_P": 0,
 "W_S1": 0,
 "borderline": [
  [
   1001,
   "H",
   0.2604739079892493
  ]
 ],
 "k_H": 30,
 "k_P": 30,
 "k_S1": 30,
 "lrExceedance": {
  "chi2_df": 6,
  "crit95": 12.591587243743977,
  "note": "preview only; validates no test",
  "rate": 0.0
 },
 "n1Flags": [],
 "none": {
  "H": 0,
  "P": 0,
  "S": 0,
  "S1": 0
 },
 "pairedPfErrSminusH": {
  "median": -0.026220345201090955,
  "n": 30,
  "note": "context only",
  "wilcoxonP": 3.2391399145126343e-06
 },
 "sensitivity": {
  "H": {
   "0.25": {
    "correct": 29,
    "none": 0,
    "wrong": 1
   },
   "0.5": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   },
   "1": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   },
   "2": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   },
   "5": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   }
  },
  "P": {
   "0.25": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   },
   "0.5": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   },
   "1": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   },
   "2": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   },
   "5": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   }
  },
  "S": {
   "0.25": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   },
   "0.5": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   },
   "1": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   },
   "2": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   },
   "5": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   }
  },
  "S1": {
   "0.25": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   },
   "0.5": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   },
   "1": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   },
   "2": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   },
   "5": {
    "correct": 30,
    "none": 0,
    "wrong": 0
   }
  }
 },
 "support": {
  "1001": 16,
  "1002": 16,
  "1003": 16,
  "1004": 16,
  "1005": 16,
  "1006": 16,
  "1007": 16,
  "1008": 16,
  "1009": 16,
  "1010": 16,
  "1011": 16,
  "1012": 16,
  "1013": 16,
  "1014": 16,
  "1015": 16,
  "1016": 16,
  "1017": 16,
  "1018": 16,
  "1019": 16,
  "1020": 16,
  "1021": 16,
  "1022": 16,
  "1023": 16,
  "1024": 16,
  "1025": 16,
  "1026": 16,
  "1027": 14,
  "1028": 16,
  "1029": 16,
  "1030": 16
 }
}
```

## Stage D stop rule

- stop (revise): False
