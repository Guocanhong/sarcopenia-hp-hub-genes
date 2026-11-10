# -*- coding: utf-8 -*-
"""Patch 50_build_cover_letter.py: align 6 spots with the manuscript's exact wording.
Line-anchored, idempotent, asserted."""
import io

p = r"D:/GEO-limma,肌少症1/70-MR/code/50_build_cover_letter.py"
s = io.open(p, encoding="utf-8").read()

pairs = [
    # 5) reverse causation: remove invented mechanism, use manuscript framing
    ("'susceptible to reverse causation \\u2014 systemic illness may reduce eating, '",
     "'vulnerable to reverse causation and residual confounding \\u2014 precisely the '"),
    ("'activity, and antibody titers alike. Whether the infection-muscle association '",
     "'biases that a genetic design is built to mitigate. Whether the '\n    'infection-muscle association '"),
    # 1) antigen list: real names
    ("(CagA, GroEL, UreA, OMP, HcpC, NabA)",
     "(CagA, catalase, GroEL, OMP, UreA, VacA)"),
    # 2) decision rule: manuscript rule = FE+WM corroboration; Egger is pleiotropy check
    ("'every signal was classified by a decision rule fixed in advance: only effects '",
     "'every signal was classified by a decision rule fixed in advance: only signals '"),
    ("'supported by the primary random-effects estimator, the weighted median, and '",
     "'corroborated by the fixed-effect and weighted-median estimators are reported '"),
    ("'MR-Egger concordance are reported as supported, whereas random-effects-only '",
     "'as supported, whereas random-effects-only '"),
    # 3) beta units
    ("0.061 SD per log-odds of seropositivity",
     "0.061 SD units per log-odds of seropositivity"),
    # 4) NFKBIA: participation, not direction/activity readout
    ("'higher predicted NFKBIA expression \\u2014 a readout of NF-\\u03baB pathway activity '",
     "'higher predicted NFKBIA expression \\u2014 an NF-\\u03baB-related transcript \\u2014 '"),
    ("'\\u2014 to lower lean mass. We report these signals with their decision-rule '",
     "'\\u2014 to lower lean mass: evidence that the pathway participates in the '\n"
     "    'phenotype rather than of one specific direction of NF-\\u03baB signaling. We '\n"
     "    'report these signals with their decision-rule '"),
    # 6) AI disclosure per Wolters Kluwer policy
    ("'`70-MR`) and will be archived with a Zenodo DOI upon acceptance.',",
     "'`70-MR`) and will be archived with a Zenodo DOI upon acceptance. In line with '\n"
     "    'journal policy, we disclose that Tencent WorkBuddy (a GLM large language '\n"
     "    'model, Tencent) was used to assist with language editing and analysis-code '\n"
     "    'development; the manuscript contains a corresponding disclosure, and the '\n"
     "    'authors reviewed and take full responsibility for all content.',"),
]

n = 0
for old, new in pairs:
    if new in s and old not in s:
        print("SKIP (applied):", old[:50])
        continue
    assert old in s, "NOT FOUND: " + old[:60]
    assert s.count(old) == 1, "NOT UNIQUE: " + old[:60]
    s = s.replace(old, new, 1)
    n += 1

io.open(p, "w", encoding="utf-8", newline="").write(s)
print("patched:", n)

# readback: none of the wrong anchors remain
for bad in ["HcpC", "NabA", "readout of", "systemic illness", "MR-Egger concordance are reported"]:
    assert bad not in s, "still present: " + bad
print("readback clean")
