# 54_fill_tableS2.py -- fill Table S2 Egger-intercept P cells and add a
# weighted-mode P column from 53_egger_mode_backfill.csv (mock-review backfill).
# Idempotent: safe to re-run; only fills '-' cells and inserts the mode column once.
import csv

BACKFILL = "D:/MRout/70-MR/53_egger_mode_backfill.csv"
TABLES = "D:/GEO-limma,肌少症1/70-MR/paper/Tables_draft.md"

val = {}
for r in csv.DictReader(open(BACKFILL, encoding="utf-8")):
    val[r["pair"].strip().lower()] = (float(r["intercept_p"]), float(r["mode_p"]))

ALIAS = {
    "seropositivity -> right grip": "seropositivity -> grip_right",
    "seropositivity -> left grip": "seropositivity -> grip_left",
    "seropositivity -> walking pace": "seropositivity -> walking",
    "seropositivity -> low grip (>=60y)": "seropositivity -> lowgrip60",
    "igg level -> alm": "igg level (titer) -> alm",
    "nfkbia -> right grip": "nfkbia -> grip_right",
    "nfkbia -> left grip": "nfkbia -> grip_left",
    "nfkbia -> walking pace": "nfkbia -> walking",
    "nfkbia -> low grip (>=60y)": "nfkbia -> lowgrip60",
    "nfkbia -> alm": "nfkbia -> alm",
    "c1qa -> right grip": "c1qa -> grip_right",
    "c1qa -> walking pace": "c1qa -> walking",
    "c1qa -> low grip (>=60y)": "c1qa -> lowgrip",
    "ctnnbip1 -> right grip": "ctnnbip1 -> grip_right",
    "ctnnbip1 -> walking pace": "ctnnbip1 -> walking",
    "ctnnbip1 -> low grip (>=60y)": "ctnnbip1 -> lowgrip",
}
def lookup(pair):
    p = pair.strip().lower()
    p = ALIAS.get(p, p)
    return val.get(p)

fmt3 = lambda x: ("%.3g" % x)

lines = open(TABLES, encoding="utf-8").read().split("\n")
out = []
in_s2 = False
fi = fm = 0
for ln in lines:
    if ln.startswith("## Table S2"):
        in_s2 = True
    elif in_s2 and ln.startswith("## ") and not ln.startswith("## Table S2"):
        in_s2 = False
    if in_s2 and ln.startswith("|"):
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        if cells[0] == "Line":
            if "Weighted mode P" not in cells:
                cells.insert(3, "Weighted mode P")
            out.append("| " + " | ".join(cells) + " |")
            continue
        if cells[0] and set(cells[0]) <= set("-: "):
            if len(cells) == 6:
                cells.insert(3, "---")
            out.append("| " + " | ".join(cells) + " |")
            continue
        # data row; after first pass there are 7 cells, first pass 6
        pair = cells[1].strip().lower()
        v = lookup(pair)
        if len(cells) == 6:  # first pass: insert mode column
            if v and cells[2] in ("-", ""):
                cells[2] = fmt3(v[0]); fi += 1
            cells.insert(3, fmt3(v[1]) if v else "-")
            if v: fm += 1
        else:                # idempotent pass: fill remaining '-', refresh mode cells
            if v:
                if cells[2] in ("-", ""):
                    cells[2] = fmt3(v[0]); fi += 1
                if v[1] == v[1]:  # always refresh mode cell to match latest backfill
                    if cells[3] != fmt3(v[1]):
                        cells[3] = fmt3(v[1]); fm += 1
        out.append("| " + " | ".join(cells) + " |")
        continue
    out.append(ln)

open(TABLES, "w", encoding="utf-8").write("\n".join(out))
print("intercept filled this pass:", fi, "| mode filled this pass:", fm)
print("remaining '-' intercept cells:", [l.split("|")[2].strip() for l in out
      if l.startswith("|") and "->" in l and len(l.split("|")) >= 4 and l.split("|")[2].strip() == "-"])
