# 1-filter.py  GSE249874 RAW mtx 流式过滤（gastritis 6 samples）
# 结构（已实测）: 36,601 genes x 122,307,840 droplets, nnz=324,107,674
# barcodes 带样本后缀 -1..-18；胃炎层 = 后缀 -7,-8,-9 (GS-HP-N) / -10,-11,-12 (GS-HP-P)
# 内存策略: 不保存全量条码字符串；后缀用 uint8 数组；QC 后只取胃炎层真实细胞的字符串
import gzip, numpy as np, time
from collections import Counter

RAW = "raw"
KEEP_SUFFIX = (7, 8, 9, 10, 11, 12)
SAMPLE_NAME = {7: "GS-HP-N_1", 8: "GS-HP-N_2", 9: "GS-HP-N_3",
               10: "GS-HP-P_1", 11: "GS-HP-P_2", 12: "GS-HP-P_3"}

def suf_of(bc: str):
    suf = bc.rsplit("-", 1)[1] if "-" in bc else ""
    return int(suf) if suf.isdigit() else 0

t0 = time.time()
# ---- Pass0: 流式读条码 -> 后缀 uint8 数组 ----
print("[pass0] streaming barcodes ...", flush=True)
sufs = []
with gzip.open(f"{RAW}/barcodes_full", "rt") as f:
    for line in f:
        sufs.append(suf_of(line.strip()))
sufs = np.array(sufs, dtype=np.uint8)
nb = sufs.shape[0]
print(f"[pass0] barcodes={nb} elapsed={time.time()-t0:.0f}s", flush=True)
print("suffix hist:", dict(sorted(Counter(sufs.tolist()).items())), flush=True)

# ---- features ----
genes = []
with gzip.open(f"{RAW}/features_full", "rt") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        genes.append(p[1] if len(p) > 1 else p[0])
ng = len(genes)
print("features:", ng, flush=True)
TARGET = ["C1QA", "C1QB", "CD68", "LYZ", "NFKBIA", "CTNNBIP1", "DTL",
          "PTPRC", "HBB", "EPCAM", "MUC5AC", "GZMB", "ACTA2", "PECAM1"]
trow = {}
for k, g in enumerate(TARGET):
    if g in genes:
        trow[genes.index(g)] = k
print("targets found:", [TARGET[v] for v in trow.values()], flush=True)

# ---- Pass1: 全矩阵流式 -> 每液滴 UMI/基因数 ----
umi = np.zeros(nb, dtype=np.float32)
nfeat = np.zeros(nb, dtype=np.int32)
tg_sum = np.zeros(len(TARGET), dtype=np.float64)
path = f"{RAW}/matrix_full"
with gzip.open(path, "rt") as f:
    banner = f.readline()
    dims = None
    for _ in range(3):
        line = f.readline()
        if not line.startswith("%"):
            dims = line.split(); break
    nr, nc, nnz = int(dims[0]), int(dims[1]), int(dims[2])
    print(f"mtx dims: {nr} x {nc}, nnz={nnz}", flush=True)
    assert nc == nb, f"barcode mismatch: mtx {nc} vs {nb}"
    cnt = 0
    for line in f:
        cnt += 1
        if cnt % 40000000 == 0:
            print(f"  pass1 {cnt/nnz*100:.0f}%  {time.time()-t0:.0f}s", flush=True)
        p = line.split()
        if len(p) < 3:
            continue
        j = int(p[1]) - 1
        v = float(p[2])
        umi[j] += v
        nfeat[j] += 1
        k = trow.get(int(p[0]) - 1)
        if k is not None:
            tg_sum[k] += v
print(f"[pass1] done {time.time()-t0:.0f}s", flush=True)
np.save("qc_umi.npy", umi); np.save("qc_nfeat.npy", nfeat)
for k, g in enumerate(TARGET):
    if g in genes:
        print(f"  total UMI {g}: {tg_sum[k]:.0f}", flush=True)

# ---- QC：只对胃炎层判定 ----
keep_mask = (nfeat >= 200) & (umi >= 1000)
for s in KEEP_SUFFIX:
    m = (sufs == s)
    print(f"sample suffix -{s}: droplets={int(m.sum())}, QC-pass={int((m & keep_mask).sum())}", flush=True)

sel = np.nonzero(keep_mask & np.isin(sufs, KEEP_SUFFIX))[0]
print("kept gastritis cells:", len(sel), f"{time.time()-t0:.0f}s", flush=True)
sel_set = set((sel + 1).tolist())            # 1-based 列号

# ---- 取回保留条码字符串 ----
kept_bc = {}
with gzip.open(f"{RAW}/barcodes_full", "rt") as f:
    for i, line in enumerate(f):
        if (i + 1) in sel_set:
            kept_bc[i + 1] = line.strip()
print("barcodes fetched:", len(kept_bc), flush=True)

# ---- Pass2: 流式写出过滤后 mtx ----
sel_sorted = np.sort(sel)
newidx = np.full(nc, -1, dtype=np.int64)
newidx[sel_sorted] = np.arange(1, len(sel_sorted) + 1)
nnz2 = int(nfeat[sel_sorted].sum())
with open("sample_map.csv", "w") as f:
    f.write("barcode,sample,group\n")
    for c1 in sel_sorted.tolist():
        b = kept_bc[c1]; s = suf_of(b)
        f.write(f"{b},{SAMPLE_NAME[s]},{'HpP' if s >= 10 else 'HpN'}\n")
print("filtered nnz:", nnz2, flush=True)
with gzip.open(path, "rt") as f, open("filtered_matrix.mtx", "w") as out:
    out.write("%%MatrixMarket matrix coordinate integer general\n%\n")
    out.write(f"{nr} {len(sel_sorted)} {nnz2}\n")
    buf = []
    for line in f:
        p = line.split()
        if len(p) < 3:
            continue
        c = int(p[1])
        ni = newidx[c - 1]
        if ni > 0:
            buf.append(f"{p[0]} {ni} {p[2]}\n")
            if len(buf) >= 1000000:
                out.write("".join(buf)); buf.clear()
    out.write("".join(buf))
print(f"[pass2] done {time.time()-t0:.0f}s", flush=True)

import subprocess
subprocess.run(["gzip", "-f", "filtered_matrix.mtx"], check=True)
with open("filtered_barcodes.tsv", "w") as f:
    for c1 in sel_sorted.tolist():
        f.write(kept_bc[c1] + "\n")
with open("filtered_features.tsv", "w") as f:
    for g in genes:
        f.write(f"{g}\t{g}\tGene Expression\n")
print("ALL DONE", f"{time.time()-t0:.0f}s", flush=True)
