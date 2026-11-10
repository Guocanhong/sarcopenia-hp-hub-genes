# -*- coding: utf-8 -*-
"""40b_second_pass.py -- second compression round for the abstract (target <=300 words)
+ fix the F6 false-skip (placeholder "(eqtl-a-ENSG...)" still present because the
idempotency check matched an unrelated "through OpenGWAS" occurrence)."""
import io, os, re, sys

PAPER = r"D:/GEO-limma,肌少症1/70-MR/paper"

def read(p):  return io.open(p, encoding="utf-8").read()
def write(p, s): io.open(p, "w", encoding="utf-8", newline="").write(s)

report = []
def patch(fname, old, new, label, count=1):
    path = fname if os.path.isabs(fname) else os.path.join(PAPER, fname)
    s = read(path)
    if old not in s:
        report.append((label, "SKIP (old gone)" if new in s else "FAIL: old not found", fname)); return
    write(path, s.replace(old, new, count))
    report.append((label, "OK (%d/%d)" % (count, s.count(old)), fname))

# --- F6b placeholder removal (unique anchor) ---
patch("Methods_draft.md",
      "expression (N = 31,684) ^15^ through OpenGWAS (eqtl-a-ENSG\u2026).",
      "expression (N = 31,684) ^15^ through OpenGWAS.", "F6b drop placeholder")

# --- F2b abstract, round 2 ---
patch("Abstract_draft.md",
      "but causality remains untested; inflammatory pathways are the presumed mediators.",
      "but causality is untested; inflammation is the presumed mediator.", "F2b Background")
patch("Abstract_draft.md",
      "(iii) blood expression of four hub genes (C1QA, NFKBIA, DTL, CTNNBIP1) using cis-eQTLs (eQTLGen), each tested against appendicular lean mass (ALM), grip strength, usual walking pace and low grip strength.",
      "(iii) four blood hub genes (C1QA, NFKBIA, DTL, CTNNBIP1; cis-eQTLs from eQTLGen), each tested against appendicular lean mass (ALM), grip strength, walking pace and low grip strength.",
      "F2b Methods instruments")
patch("Abstract_draft.md",
      "Inverse-variance weighting was the primary estimator; MR-Egger, weighted median and weighted mode, Cochran's Q, leave-one-out and Steiger analyses served as sensitivity checks, with pre-specified Bonferroni corrections (\u03b1 = 0.05/24 antigen line; \u03b1 = 0.01 for NFKBIA) and a uniform decision rule classifying signals supported only by the random-effects model as suggestive.",
      "Inverse-variance weighting was primary; MR-Egger, weighted median and mode, Cochran's Q, leave-one-out and Steiger analyses were sensitivity checks; Bonferroni corrections (\u03b1 = 0.05/24 antigens; \u03b1 = 0.01 for NFKBIA) and a uniform decision rule classifying random-effects-only signals as suggestive were pre-specified.",
      "F2b Methods estimator sentence")
patch("Abstract_draft.md",
      "A reverse-MR line tested the sarcopenia traits as exposures for seropositivity.",
      "A reverse-MR line tested sarcopenia traits as exposures for seropositivity.", "F2b Methods reverse line")
patch("Abstract_draft.md",
      "Of the six antigens, anti-GroEL antibody levels showed",
      "Anti-GroEL antibody levels showed", "F2b Results GroEL lead-in")
patch("Abstract_draft.md",
      "The fixed-effect and weighted-median estimators did not support this signal, which is therefore reported as suggestive.",
      "The fixed-effect and weighted-median estimators did not support this signal, which is therefore suggestive.",
      "F2b Results decision rule")
patch("Abstract_draft.md",
      "the nominal CTNNBIP1\u2013low grip strength estimate (P = 6.7 \u00d7 10\u207b\u2078) was unsupported by all other estimators and is likewise suggestive.",
      "the nominal CTNNBIP1\u2013low grip strength estimate (P = 6.7 \u00d7 10\u207b\u2078), unsupported by all other estimators, is likewise suggestive.",
      "F2b Results CTNNBIP1")
patch("Abstract_draft.md",
      "and provides suggestive, directionally concordant evidence for the anti-GroEL immune response.",
      "and provides suggestive concordant evidence for the anti-GroEL immune response.", "F2b Conclusions")

# word-count check
ab = read(os.path.join(PAPER, "Abstract_draft.md"))
body = ab.split("## Abstract")[1].split("## Keywords")[0]
body = re.sub(r"\^[\d,]+\^", "", body)
words = len(body.split())
print("final abstract word count (incl. section labels):", words)
for lbl, status, fn in report:
    print("%-38s %-24s" % (lbl, status))
sys.exit(0 if all(not st.startswith("FAIL") for _, st, _ in report) else 1)
