# -*- coding: utf-8 -*-
"""Plan-B renumbering: remove ref 11 (companion manuscript) and shift
refs 12-28 -> 11-27 across paper drafts.

- References_draft.md: leading list numbers (>=12 -> -1) and usage-map table cells
- Discussion/Methods drafts: numbers >=12 inside {{CITE n}} / {{CITE a+b}} tokens
- Introduction/Abstract/S1: untouched by this script (resolver-map based / manual)
"""
import re, io

PAPER = r"D:/GEO-limma,肌少症1/70-MR/paper"

def shift_num(n):
    n = int(n)
    return str(n - 1) if n >= 12 else str(n)

def shift_group(g):
    return "+".join(shift_num(x) for x in g.split("+")) if "+" in g else ",".join(shift_num(x) for x in g.split(","))

# --- {{CITE n:}} tokens in Discussion / Methods ---
for fn in ["Discussion_draft.md", "Methods_draft.md"]:
    p = PAPER + "/" + fn
    t = io.open(p, encoding="utf-8").read()
    t2 = re.sub(r"\{\{CITE ([0-9+,]+):", lambda m: "{{CITE " + shift_group(m.group(1)) + ":", t)
    io.open(p, "w", encoding="utf-8", newline="\n").write(t2)
    print(fn, "cite tokens shifted")

# --- References_draft.md: list entries "N. " and usage-map "| N |" / "| 8, 9 |" ---
p = PAPER + "/References_draft.md"
t = io.open(p, encoding="utf-8").read()
# remove companion entry (line 11) but preserve a note for future re-insertion
t = t.replace(
    "11. [Companion manuscript — **DOI placeholder pending**] Guo C, et al. Integrative transcriptomic analysis of *Helicobacter pylori* infection and sarcopenia-related phenotypes [title TBC]. Under review; Research Square \"In Review\" preprint (manuscript MD-D-26-16765). DOI to insert when the preprint posting goes live.\n",
    "<!-- Plan-B (2026-11-06): former ref 11 (companion manuscript, MD-D-26-16765) removed for the submit-first version; refs 12-28 shifted to 11-27. To revert to the companion-citation version, re-insert as 11 and shift 11-27 back up by one: the companion entry read: Guo C, et al. Integrative transcriptomic analysis of Helicobacter pylori infection and sarcopenia-related phenotypes [title TBC]. Under review; Research Square In Review preprint (MD-D-26-16765). -->\n",
)
t = re.sub(r"^(\d+)\. ", lambda m: shift_num(m.group(1)) + ". ", t, flags=re.M)
# usage-map table cells: last cell of each row
def fix_row(m):
    return m.group(1) + shift_group(m.group(2)) + m.group(3)
t = re.sub(r"(\|[^\n]*\| )([0-9,]+)( \|)", fix_row, t)
io.open(p, "w", encoding="utf-8", newline="\n").write(t)
print("References_draft.md renumbered")

# --- build script maps ---
p = PAPER + "/../code/19_build_submission_md.py"
t = io.open(p, encoding="utf-8").read()
t = re.sub(r'":(\d+)"', lambda m: '":' + shift_num(m.group(1)) + '"', t)  # PMID2N
# DOI2N: only the gepi.21998 entry is >=12
t = t.replace('"10.1002/gepi.21998":"27"', '"10.1002/gepi.21998":"26"')
# Hemani special-case regex: CITE 18 -> CITE 17
t = t.replace(r"\{\{CITE 18:", r"\{\{CITE 17:")
# title-page reference count
t = t.replace("28 references.", "27 references.")
io.open(p, "w", encoding="utf-8", newline="\n").write(t)
print("build script maps shifted")
