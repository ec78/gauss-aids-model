"""pilot_posthoc_pooled.py -- POST-HOC EXPLORATORY analysis (not part of
the frozen PILOT_PROTOCOL.md comparison). Pools the genuinely converged
candidates of several methods/starts per seed and applies the SAME
operational rule (largest critExact, converged only). Run from this
directory:  python pilot_posthoc_pooled.py > pilot_out/pilot_posthoc_pooled_output.txt
"""
import pilot_analyze as pa

df = pa.load_cells("pilot_out")
configs = {
    "cold LM, Stone only": [(2, [1])],
    "continuation, Stone only": [(3, [1])],
    "cold LM + continuation, Stone only (2 solves)": [(2, [1]), (3, [1])],
    "cold LM 16 starts": [(2, range(1, 17))],
    "continuation 16 starts": [(3, range(1, 17))],
    "cold LM 16 + continuation Stone (17 solves)": [(2, range(1, 17)), (3, [1])],
    "cold LM 16 + continuation 16 (32 solves)": [(2, range(1, 17)), (3, range(1, 17))],
}
print("POST-HOC EXPLORATORY: pooled selection by largest critExact among genuinely converged candidates")
print("columns: config | bS<=.5/1/2/5 | raw<=.5/1/2/5 | never | wrong(bS>1) seeds | wrong(raw>1) seeds | median / max wall per seed (s)")
for model in (0, 1):
    print(f"\n== {pa.MODELS[model]} ==")
    for name, parts in configs.items():
        bs, raw, never, wb, wr, walls = [], [], 0, [], [], []
        for seed in range(1, 31):
            g = df[(df.model == model) & (df.seed == seed)]
            mask = False
            for m, starts in parts:
                mask = mask | ((g.method == m) & g.start.isin(list(starts)))
            c = g[mask]
            walls.append(c.wallTime.sum())
            best = None
            for _, r in c.iterrows():
                if r.conv and r.critExact == r.critExact and (best is None or r.critExact > best.critExact):
                    best = r
            if best is None:
                never += 1
                continue
            bs.append(best.bsErr); raw.append(best.rawStructErr)
            if best.bsErr > 1: wb.append(seed)
            if best.rawStructErr > 1: wr.append(seed)
        thr = lambda v: "/".join(str(sum(x <= c for x in v)) for c in pa.THRESHOLDS)
        walls.sort()
        print(f"{name:48s} | {thr(bs):12s} | {thr(raw):12s} | {never} | {wb} | {wr} | {walls[15]:.1f} / {walls[-1]:.1f}")
