# -*- coding: utf-8 -*-
"""22_strict_renumber.py — renumber all references strictly by first appearance (Vancouver).
First run (2026-11-06): steps 1-2 applied; rerun executes steps 3-4 only (MAP is NOT idempotent)."""
import re, io, os

PAPER = r"D:/GEO-limma,肌少症1/70-MR/paper"
BUILD = r"D:/GEO-limma,肌少症1/70-MR/code/19_build_submission_md.py"

# old -> new  (derived from first-appearance scan of manuscript_combined.md, 2026-11-06)
MAP = {n: n for n in range(1, 11)}
MAP.update({
    18: 11, 13: 12, 11: 13, 12: 14, 16: 15, 14: 16, 15: 17,
    19: 18, 17: 19, 20: 20,
    24: 21, 25: 22, 21: 23, 23: 24, 22: 25, 28: 26,
    30: 27, 29: 28, 31: 29,
    27: 30, 26: 31,
    32: 32, 33: 33, 34: 34,
})
assert sorted(MAP.values()) == list(range(1, 35)), "MAP not a permutation of 1..34"

STEPS_1_2_DONE = True  # set False only to re-run from a pristine git/backup state

if not STEPS_1_2_DONE:
    def remap_cite_tokens(text):
        return re.sub(r"\{\{CITE ([0-9+]+):",
                      lambda m: "{{CITE " + "+".join(str(MAP[int(x)]) for x in m.group(1).split("+")) + ":",
                      text)
    for fn in ["Introduction_draft.md", "Methods_draft.md", "Results_draft.md", "Discussion_draft.md"]:
        p = os.path.join(PAPER, fn)
        t = io.open(p, encoding="utf-8").read()
        io.open(p, "w", encoding="utf-8", newline="\n").write(remap_cite_tokens(t))
        print(fn, "tokens remapped")

    p = os.path.join(PAPER, "References_draft.md")
    t = io.open(p, encoding="utf-8").read()
    m = re.search(r"(## Numbered reference list\n)(.*?)(\n## )", t, re.S)
    head, body, tail = m.group(1), m.group(2), m.group(3)
    entries = {}
    for line in body.split("\n"):
        mm = re.match(r"^(\d+)\. (.*)$", line)
        if mm:
            entries[int(mm.group(1))] = mm.group(2)
    new_lines = []
    for new_n in range(1, len(entries) + 1):
        old_n = [k for k, v in MAP.items() if v == new_n][0]
        new_lines.append(f"{new_n}. {entries[old_n]}")
    t = t[:m.start()] + head + "\n".join(new_lines) + tail + t[m.end():]
    io.open(p, "w", encoding="utf-8", newline="\n").write(t)
    print("References_draft.md reordered:", len(entries), "entries")

# ---------- 3) Discussion drafting-notes prose numbers ----------
p = os.path.join(PAPER, "Discussion_draft.md")
t = io.open(p, encoding="utf-8").read()
for old, new in [
    ("**refs 6, 21** (Mayr 2000", "**refs 6, 23** (Mayr 2000"),
    ("**refs 25, 26** (Jackson 2009", "**refs 22, 31** (Jackson 2009"),
    ("**refs 21, 22** (Singh 2024", "**refs 23, 25** (Singh 2024"),
    ("**ref 24** (Cai 2004 Cell", "**ref 21** (Cai 2004 Cell"),
    ("**ref 23** (Hoffmann 2002", "**ref 25** (Hoffmann 2002"),
    ("**ref 24** (same paper: MISR", "**ref 21** (same paper: MISR"),
    ("**refs 3, 28**.", "**refs 3, 26**."),
    ("**ref 27** (Burgess 2016", "**ref 31** (Burgess 2016"),
]:
    assert old in t, "note pattern missing: " + old
    t = t.replace(old, new)
io.open(p, "w", encoding="utf-8", newline="\n").write(t)
print("Discussion notes updated")

# ---------- 4) build script: PMID2N / DOI2N values, SUP list, Hemani pattern ----------
p = BUILD
t = io.open(p, encoding="utf-8").read()

t = re.sub(r'"(\d{8})":"(\d+)"', lambda mm: '"%s":"%d"' % (mm.group(1), MAP[int(mm.group(2))]), t)
t = re.sub(r'"(10\.[^"]+)":"(\d+)"', lambda mm: '"%s":"%d"' % (mm.group(1), MAP[int(mm.group(2))]), t)
t = re.sub(r'SUP\(\["(\d+)"\]\)', lambda mm: 'SUP(["%d"])' % MAP[int(mm.group(1))], t)
t = t.replace("CITE 17:", "CITE 19:")  # inside the Hemani author-year regex pattern (unique)
io.open(p, "w", encoding="utf-8", newline="\n").write(t)
print("build script maps remapped")
print("DONE")
