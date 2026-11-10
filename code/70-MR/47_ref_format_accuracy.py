# -*- coding: utf-8 -*-
# 47_ref_format_accuracy.py -- full reference audit fixes (2026-11-07):
#  A. ref 4 article number 56 -> 67 (PubMed: Gut Pathog 2021;13(1):67)
#  B. AMA author truncation (>6 authors -> first 3 + et al): refs 10,14,16,17,19,23,26,27
#  C. add ref 39 (Burgess & Thompson 2011, weak instruments) -- supports the
#     conservative-lower-bound claim misattributed to Burgess 2016
#  D. body fixes: EWGSOP2 -> EWGSOP (2010) x2; ^8,9^ sentence accuracy x2;
#     haemorrhage -> hemorrhage; LIMIT3 cite -> new Burgess 2011 token
# Then run 23_strict_renumber.py (separately). Idempotent.
import io, re, os

P = r"D:/GEO-limma,肌少症1/70-MR/paper"

def rd(f):
    return io.open(os.path.join(P, f), encoding="utf-8").read()

def wr(f, s):
    io.open(os.path.join(P, f), "w", encoding="utf-8", newline="").write(s)

def sub(s, old, new, tag):
    if new in s:
        print("SKIP:", tag); return s
    assert old in s, "anchor missing: " + tag
    assert s.count(old) == 1, "anchor not unique: " + tag
    print("OK:", tag)
    return s.replace(old, new, 1)

# ---------- A/B: References_draft.md ----------
s = rd("References_draft.md")

s = sub(s,
    "Gut Pathog. 2021;13(1):56. PMID 34782007",
    "Gut Pathog. 2021;13(1):67. PMID 34782007 *(article number corrected to 67, verified via PubMed/eutils 2026-11-07)*",
    "ref 4 pages")

s = sub(s,
    "10. Mei Z, Wan Y, Mao X, Qi W, Zhou G, Chen M, Li X, Pan H, Wang D. A Mendelian",
    "10. Mei Z, Wan Y, Mao X, et al. A Mendelian",
    "ref 10 authors")

s = sub(s,
    "14. Chong AHW, Mitchell RE, Hemani G, Davey Smith G, Yolken RH, Richmond RC, Paternoster L. Genetic analyses",
    "14. Chong AHW, Mitchell RE, Hemani G, et al. Genetic analyses",
    "ref 14 authors")

s = sub(s,
    "16. Pei YF, Liu YZ, Yang XL, Zhang H, Feng GJ, Wei XT, Zhang L. The genetic architecture",
    "16. Pei YF, Liu YZ, Yang XL, et al. The genetic architecture",
    "ref 16 authors")

s = sub(s,
    "17. Jones G, Trajanoska K, Santanasto AJ, Stringa N, Kuo CL, Atkins JL, et al. Genome-wide",
    "17. Jones G, Trajanoska K, Santanasto AJ, et al. Genome-wide",
    "ref 17 authors")

s = sub(s,
    "19. Hemani G, Zheng J, Elsworth B, Wade KH, Haberland V, Baird D, et al. The MR-Base platform",
    "19. Hemani G, Zheng J, Elsworth B, et al. The MR-Base platform",
    "ref 19 authors")

s = sub(s,
    "23. Jackson L, Britton J, Lewis SA, McKeever TM, Atherton J, Fullerton D, Fogarty AW. A population-based",
    "23. Jackson L, Britton J, Lewis SA, et al. A population-based",
    "ref 23 authors")

s = sub(s,
    "26. Singh MK, Shin Y, Han S, Ha J, Tiwari PK, Kim SS, Kang I. Molecular chaperonin",
    "26. Singh MK, Shin Y, Han S, et al. Molecular chaperonin",
    "ref 26 authors")

s = sub(s,
    "27. Cai D, Frantz JD, Tawa NE Jr, Melendez PA, Oh BC, Lidov HGW, et al. IKK\u03b2/NF-\u03baB activation",
    "27. Cai D, Frantz JD, Tawa NE Jr, et al. IKK\u03b2/NF-\u03baB activation",
    "ref 27 authors")

# C: append ref 39 inside the "## Numbered reference list" section (before the
# "## {{CITE}}" heading), keeping the list sequential 1..39 for renumbering.
NEWREF = ("39. Burgess S, Thompson SG, CRP CHD Genetics Collaboration. Avoiding bias from "
          "weak instruments in Mendelian randomization studies. Int J Epidemiol. "
          "2011;40(3):755\u2013764. PMID 21414999. DOI 10.1093/ije/dyr036 "
          "*(added 2026-11-07: supports the weak-instrument attenuation claim in "
          "Limitations Third, previously misattributed to the overlap-bias paper; "
          "metadata verified via PubMed/eutils)*\n\n")
if "PMID 21414999" not in s:
    assert "\n## {{CITE}}" in s
    s = s.replace("\n## {{CITE}}", "\n" + NEWREF + "## {{CITE}}", 1)
    print("OK: ref 39 appended")
else:
    print("SKIP: ref 39 present")
wr("References_draft.md", s)

# ---------- D: body fixes ----------
d = rd("Discussion_draft.md")

d = sub(d,
    "(Burgess et al., 2016) {{CITE 21: Burgess S, Davies NM, Thompson SG. Genet Epidemiol "
    "2016;40(7):597\u2013608, DOI 10.1002/gepi.21998 \u2014 weak-instrument attenuation and "
    "overlap-bias bounds}}",
    "(Burgess and Thompson, 2011) {{CITE: Burgess S, Thompson SG, CRP CHD Genetics "
    "Collaboration. Int J Epidemiol 2011;40(3):755\u2013764, PMID 21414999, DOI "
    "10.1093/ije/dyr036 \u2014 weak-instrument attenuation toward the null}}",
    "LIMIT3 cite swap")

d = sub(d,
    "which showed that broad-panel MR screens generate many signals that fail rigorous "
    "validation (Zhang et al., 2026; Machnik et al., 2026)",
    "which documented the fragility of conventional MR pipelines \u2014 pleiotropy-driven "
    "inconsistency when many instruments are used and previously reported causal paths "
    "that fail rigorous re-analysis (Zhang et al., 2026; Machnik et al., 2026)",
    "STRENGTHS v reword")

d = sub(d, "subarachnoid haemorrhage", "subarachnoid hemorrhage", "haemorrhage")
wr("Discussion_draft.md", d)

s = rd("Introduction_draft.md")
s = sub(s,
    "showed this directly and argued for targeted, hypothesis-driven designs with "
    "systematic validation across complementary data sources (Zhang et al., 2026; "
    "Machnik et al., 2026)",
    "showed that conventional pipelines mislead in exactly this regime \u2014 through "
    "pleiotropy when many instruments are used and through previously reported paths "
    "that dissolve on systematic re-analysis \u2014 and argued for robust instrument "
    "selection and validation across complementary data sources (Zhang et al., 2026; "
    "Machnik et al., 2026)",
    "INTRO 8,9 reword")
wr("Introduction_draft.md", s)

s = rd("Methods_draft.md")
s = sub(s, "EWGSOP2 criterion", "EWGSOP (2010) criterion", "Methods EWGSOP")
wr("Methods_draft.md", s)

s = rd("Tables_draft.md")
s = sub(s, "Low grip strength (>=60y, EWGSOP2)", "Low grip strength (>=60y, EWGSOP)", "Table S1 EWGSOP")
wr("Tables_draft.md", s)

# ---------- readback ----------
s = rd("References_draft.md")
for must in ["13(1):67", "Mei Z, Wan Y, Mao X, et al.", "PMID 21414999"]:
    assert must in s, must
for gone in ["13(1):56", "Stringa N", "Wade KH, Haberland", "McKeever TM", "Tiwari PK", "Lidov HGW"]:
    assert gone not in s, gone
d = rd("Discussion_draft.md")
assert "Burgess and Thompson, 2011" in d and "hemorrhage" in d and "myofibers" in d
assert "EWGSOP2" not in rd("Methods_draft.md") and "EWGSOP2" not in rd("Tables_draft.md")
print("readback OK")
