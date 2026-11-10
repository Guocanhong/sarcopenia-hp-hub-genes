# -*- coding: utf-8 -*-
"""Strict Vancouver renumbering — dynamic version (v2).

Computes the old->new mapping from FIRST-APPEARANCE ORDER in the main text
(Introduction -> Methods -> Results -> Discussion), resolving BOTH numbered
{{CITE N}} tokens, literal ^N^ superscripts, and un-numbered {{CITE: PMID ...}}
tags (via PMID/DOI parsed from the reference list itself). Then applies to:
  1. {{CITE N(+M,): ...}} tokens in all four drafts
  2. literal ^N^ superscripts in all drafts
  3. the numbered list in References_draft.md (reordered + renumbered)
  4. the {{CITE}} -> number mapping table in References_draft.md (regenerated)
Un-numbered {{CITE: ...}} tags are left as-is (19_build_submission_md.py resolves
them dynamically from the reordered list).
"""
import io, os, re

PAPER = r"D:/GEO-limma,肌少症1/70-MR/paper"
ORDER = ["Introduction_draft.md", "Methods_draft.md", "Results_draft.md", "Discussion_draft.md"]

def rd(fn):
    return io.open(os.path.join(PAPER, fn), encoding="utf-8").read()

def wr(fn, t):
    io.open(os.path.join(PAPER, fn), "w", encoding="utf-8", newline="\n").write(t)

# ---------- 0) sanity: literal superscripts present in drafts ----------
for fn in ORDER:
    lits = sorted(set(int(x) for x in re.findall(r"\^(\d+)\^", rd(fn))))
    if lits:
        print(fn, "literal ^N^ refs:", lits)

# ---------- 0b) PMID/DOI maps from the CURRENT reference list ----------
refs_raw = rd("References_draft.md")
list_part = refs_raw.split("## Numbered reference list")[1].split("## {{CITE}}")[0]
entries = re.findall(r"^(\d+)\. (.*)$", list_part, re.M)
old_nums = [int(n) for n, _ in entries]
N = len(old_nums)
assert old_nums == list(range(1, N + 1)), "refs list not sequential: %s" % old_nums

PMID2N, DOI2N = {}, {}
for n_s, body in entries:
    n = int(n_s)
    for pm in re.findall(r"PMID (\d{7,9})", body):
        PMID2N[pm] = n
    for dm in re.findall(r"DOI (10\.\d{4,5}/[^\s;]+)", body):
        DOI2N[dm.rstrip(".")] = n
print("PMID map size:", len(PMID2N), "| DOI map size:", len(DOI2N))

# ---------- 1) first-appearance order over the 4 main-text drafts ----------
TOKEN_NUM = re.compile(r"\{\{CITE ([\d,+]+):")
TOKEN_UNNUM = re.compile(r"\{\{CITE: ([^}]*)\}\}")
LITERAL = re.compile(r"\^(\d+)\^")

def nums_from_body(body):
    out = []
    for pm in re.findall(r"PMID (\d{7,9})", body):
        if pm in PMID2N and PMID2N[pm] not in out:
            out.append(PMID2N[pm])
    for dm in re.findall(r"10\.\d{4,5}/[^\s;,)]+", body):
        dm = dm.rstrip(".")
        if dm in DOI2N and DOI2N[dm] not in out:
            out.append(DOI2N[dm])
    return out

seq = []
for fn in ORDER:
    t = rd(fn)
    pos = []
    for m in TOKEN_NUM.finditer(t):
        pos.append((m.start(), [int(x) for x in m.group(1).replace("+", ",").split(",")]))
    for m in TOKEN_UNNUM.finditer(t):
        pos.append((m.start(), nums_from_body(m.group(1))))
    for m in LITERAL.finditer(t):
        pos.append((m.start(), [int(m.group(1))]))
    pos.sort()
    for _, nums in pos:
        for n in nums:
            if n not in seq:
                seq.append(n)
print("\nfirst-appearance sequence:", seq)

missing = [n for n in old_nums if n not in seq]
extra = [n for n in seq if n not in old_nums]
assert not missing, "refs never cited in text: %s" % missing
assert not extra, "cited numbers not in refs list: %s" % extra
assert len(seq) == N, "cited %d unique but list has %d" % (len(seq), N)

MAP = {old: i + 1 for i, old in enumerate(seq)}
INVERT = {v: k for k, v in MAP.items()}
print("map built: %d refs | unchanged positions: %d" % (N, sum(1 for k in MAP if MAP[k] == k)))

def map_group(g):
    return ",".join(str(MAP[int(x)]) for x in g.replace("+", ",").split(","))

# ---------- 2) apply to drafts (numbered tokens + literal superscripts) ----------
for fn in ORDER:
    t = rd(fn)
    t = TOKEN_NUM.sub(lambda m: "{{CITE %s:" % map_group(m.group(1)), t)
    t = LITERAL.sub(lambda m: "^%d^" % MAP[int(m.group(1))], t)
    wr(fn, t)

# ---------- 3) reorder + renumber References_draft.md ----------
lines = list_part.split("\n")
entry_re = re.compile(r"^(\d+)\. ")
out, buf, cur = [], [], None
def flush():
    global cur, buf
    if cur is not None:
        out.append("%d. %s" % (MAP[cur], "".join(buf).strip()))
for ln in lines:
    m = entry_re.match(ln)
    if m:
        flush()
        cur = int(m.group(1))
        buf = [ln[len(m.group(0)):]]
    elif cur is not None and ln.strip():
        buf.append("\n" + ln)
flush()
out.sort(key=lambda e: int(e.split(".", 1)[0]))
got = [int(e.split(".", 1)[0]) for e in out]
assert got == list(range(1, N + 1)), "renumber failed: %s" % got

head = refs_raw.split("## Numbered reference list")[0]
tail = "## {{CITE}}" + refs_raw.split("## {{CITE}}", 1)[1]
refs_new = head + "## Numbered reference list\n\n" + "\n".join(out) + "\n\n" + tail

# ---------- 4) regenerate the mapping table ----------
bodies = {MAP[int(n)]: b for n, b in entries}
sec_of = {}
for fn in ORDER:
    t = rd(fn)
    pos = []
    for m in TOKEN_NUM.finditer(t):
        pos.append((m.start(), [int(x) for x in m.group(1).replace("+", ",").split(",")]))
    for m in TOKEN_UNNUM.finditer(t):
        pos.append((m.start(), nums_from_body(m.group(1))))
    for m in LITERAL.finditer(t):
        pos.append((m.start(), [int(m.group(1))]))
    pos.sort()
    for _, nums in pos:
        for n in nums:
            sec_of.setdefault(n, fn.replace("_draft.md", "").replace("Draft", ""))

tbl = ["## {{CITE}} → reference-number mapping (for typesetting)", "",
       "**Regenerated by `code/23_strict_renumber.py` after strict first-appearance renumbering (%d refs).**" % N, "",
       "| Ref # | First author, year | First cited in |", "|---|---|---|"]
for n in range(1, N + 1):
    body = bodies[n]
    au = re.search(r"^([A-ZÀ-Ž][\w'’\-]+)", body)
    first_author = au.group(1) if au else body[:20]
    yr = (re.findall(r"\.\s?(19[89]\d|20[0-2]\d);", body) or ["?"])[0]
    tbl.append("| %d | %s (%s) | %s |" % (n, first_author, yr, sec_of.get(n, "?")))
tbl.append("")

t2 = refs_new.split("## {{CITE}}", 1)[1]
map_hdr = "## {{CITE}} → reference-number mapping"
start = refs_new.find(map_hdr)
end = refs_new.find("## Remaining pending items")
assert 0 < start < end, "table markers not found (start=%d end=%d)" % (start, end)
refs_final = refs_new[:start] + "\n".join(tbl) + "\n" + refs_new[end:]
wr("References_draft.md", refs_final)

print("\ndone: %d refs renumbered strictly by first appearance" % N)
for n in range(1, N + 1):
    body = bodies[n]
    au = re.search(r"^([A-ZÀ-Ž][\w'’\-]+)", body)
    print("  %2d <- old %-2d  %s" % (n, INVERT[n], (au.group(1) if au else body[:20]) + " " + ((re.findall(r"\.\s?(19[89]\d|20[0-2]\d);", body) or ["?"])[0])))
