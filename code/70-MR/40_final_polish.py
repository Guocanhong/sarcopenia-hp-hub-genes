# -*- coding: utf-8 -*-
"""40_final_polish.py -- final copy-edit before submission (2026-11-07).
Keeps all viewpoints, technical route, data and parameters unchanged.
F1  build: strip draft-version headings ("... -- Draft v0.x") from sections
F2  abstract: compress to <=300 words (Medicine limit); drop self-report N
    (never carried forward) + tighten decision-rule wording
F3  terminology: general quantitative exposure = "IgG titre" everywhere
    (Methods/Results/Table 1/Table 2/Fig 1 legend); antigen-specific
    exposures stay "antibody levels"; abstract "antibody titres"->"levels"
F4  Results 3.5: split the 90-word closing sentence
F5  remove empty self-praise / over-modifiers: "transparently",
    "appropriately conservative"; fix "ambiguous ... ambiguities"
F6  drop placeholder "(eqtl-a-ENSG...)" from Methods 2.2.2
F7  build: strip italic editorial notes + the **Note:** block from the
    numbered reference list (they stay in References_draft.md)
F8  build: abstract word count excludes the keywords line
"""
import io, os, re, sys

PAPER = r"D:/GEO-limma,肌少症1/70-MR/paper"
CODE = r"D:/GEO-limma,肌少症1/70-MR/code"

def read(p):  return io.open(p, encoding="utf-8").read()
def write(p, s): io.open(p, "w", encoding="utf-8", newline="").write(s)

report = []
def patch(fname, old, new, label, count=1):
    path = fname if os.path.isabs(fname) else os.path.join(PAPER, fname)
    s = read(path)
    if new in s:
        report.append((label, "SKIP (applied)", fname)); return
    n = s.count(old)
    if n == 0:
        report.append((label, "FAIL: old not found", fname)); return
    write(path, s.replace(old, new, count))
    report.append((label, "OK (%d/%d)" % (count, n), fname))

# ================= F2 abstract compression =================
AB_OLD_BG = ("**Background:** Observational studies link chronic *Helicobacter pylori* infection to lower "
 "muscle mass and sarcopenia, but causality has never been tested; inflammatory pathways are the presumed "
 "mediators.")
AB_NEW_BG = ("**Background:** Observational studies link chronic *Helicobacter pylori* infection to lower "
 "muscle mass and sarcopenia, but causality remains untested; inflammatory pathways are the presumed "
 "mediators.")

AB_OLD_ME = ("**Methods:** In a two-sample Mendelian randomization (MR) framework with four complementary "
 "evidence lines, we instrumented (i) *H. pylori* infection status (seropositivity, N = 8,735; self-report, "
 "N = 462,933; IgG titre, N = 4,683), (ii) six antigen-specific antibody titres (CagA, catalase, GroEL, OMP, "
 "UreA, VacA) and (iii) blood expression of four candidate hub genes (C1QA, NFKBIA, DTL, CTNNBIP1) using "
 "cis-eQTLs (eQTLGen, N = 31,684), testing each against appendicular lean mass (ALM), grip strength, usual "
 "walking pace and low grip strength. Inverse-variance weighting was the primary estimator; MR-Egger, "
 "weighted median and weighted mode, Cochran's Q, leave-one-out and Steiger directionality analyses served "
 "as sensitivity checks, with Bonferroni correction (\u03b1 = 0.05/24 within the antigen line; \u03b1 = 0.01 for "
 "NFKBIA) and a uniform decision rule that signals supported only by the random-effects model are reported "
 "as suggestive. A reverse-MR line tested sarcopenia-related traits as exposures for seropositivity.")
AB_NEW_ME = ("**Methods:** In a two-sample Mendelian randomization (MR) framework, we instrumented (i) "
 "*H. pylori* seropositivity (N = 8,735) and circulating anti-*H. pylori* IgG titre (N = 4,683), (ii) six "
 "antigen-specific antibody levels (CagA, catalase, GroEL, OMP, UreA, VacA) and (iii) blood expression of "
 "four hub genes (C1QA, NFKBIA, DTL, CTNNBIP1) using cis-eQTLs (eQTLGen), each tested against appendicular "
 "lean mass (ALM), grip strength, usual walking pace and low grip strength. Inverse-variance weighting was "
 "the primary estimator; MR-Egger, weighted median and weighted mode, Cochran's Q, leave-one-out and Steiger "
 "analyses served as sensitivity checks, with pre-specified Bonferroni corrections (\u03b1 = 0.05/24 antigen "
 "line; \u03b1 = 0.01 for NFKBIA) and a uniform decision rule classifying signals supported only by the "
 "random-effects model as suggestive. A reverse-MR line tested the sarcopenia traits as exposures for "
 "seropositivity.")

AB_OLD_RE = ("**Results:** Genetically predicted seropositivity lowered ALM (\u03b2 = \u22120.061 SD, P = 3.7 \u00d7 "
 "10\u207b\u2074; 11 instruments), with concordant alternative estimators and supportive sensitivity analyses. Of "
 "the six antigens, anti-GroEL antibody levels showed a nominally significant random-effects association with "
 "lower ALM (\u03b2 = \u22120.021 SD, P = 3.1 \u00d7 10\u207b\u2078). The fixed-effect and weighted-median estimators did "
 "not support this signal, which is therefore reported as suggestive under the uniform decision rule. Higher "
 "genetically predicted NFKBIA expression likewise lowered ALM (\u03b2 = \u22120.044 SD, P = 4.9 \u00d7 10\u207b\u2074; "
 "weighted median P = 2.3 \u00d7 10\u207b\u2078). No robust effects were observed for grip strength, walking pace, C1QA "
 "or DTL, or in the reverse direction; one nominal random-effects estimate for CTNNBIP1 and low grip strength "
 "(P = 6.7 \u00d7 10\u207b\u2078) was unsupported by all other estimators and is likewise reported as suggestive.")
AB_NEW_RE = ("**Results:** Genetically predicted seropositivity lowered ALM (\u03b2 = \u22120.061 SD, P = 3.7 \u00d7 "
 "10\u207b\u2074; 11 instruments), with concordant sensitivity analyses. Of the six antigens, anti-GroEL antibody "
 "levels showed a nominally significant random-effects association with lower ALM (\u03b2 = \u22120.021 SD, "
 "P = 3.1 \u00d7 10\u207b\u2078). The fixed-effect and weighted-median estimators did not support this signal, which "
 "is therefore reported as suggestive. Higher genetically predicted NFKBIA expression likewise lowered ALM "
 "(\u03b2 = \u22120.044 SD, P = 4.9 \u00d7 10\u207b\u2074; weighted median P = 2.3 \u00d7 10\u207b\u2078). No robust effects were "
 "observed for grip strength, walking pace, C1QA or DTL, or in the reverse direction; the nominal "
 "CTNNBIP1\u2013low grip strength estimate (P = 6.7 \u00d7 10\u207b\u2078) was unsupported by all other estimators and is "
 "likewise suggestive.")

AB_OLD_CO = ("**Conclusions:** Genetic evidence supports a causal effect of *H. pylori* infection on lower "
 "appendicular lean mass, identifies NFKBIA/NF-\u03baB signalling as a candidate pathway, and provides "
 "suggestive, directionally concordant evidence implicating the anti-GroEL immune response. If confirmed, "
 "*H. pylori* screening and eradication may merit evaluation in sarcopenia prevention.")
AB_NEW_CO = ("**Conclusions:** Genetic evidence supports a causal effect of *H. pylori* infection on lower "
 "appendicular lean mass, identifies NFKBIA/NF-\u03baB signalling as a candidate pathway, and provides "
 "suggestive, directionally concordant evidence for the anti-GroEL immune response. If confirmed, "
 "*H. pylori* screening and eradication may merit evaluation in sarcopenia prevention.")

for lbl, o, n in [("F2 abstract Background", AB_OLD_BG, AB_NEW_BG),
                  ("F2 abstract Methods", AB_OLD_ME, AB_NEW_ME),
                  ("F2 abstract Results", AB_OLD_RE, AB_NEW_RE),
                  ("F2 abstract Conclusions", AB_OLD_CO, AB_NEW_CO)]:
    patch("Abstract_draft.md", o, n, lbl)

# word-count check of the new abstract
_ab_words = len(re.sub(r"\^[\d,]+\^", "", " ".join([AB_NEW_BG, AB_NEW_ME, AB_NEW_RE, AB_NEW_CO])).split())
print("new abstract word count:", _ab_words)
if _ab_words > 300:
    print("WARNING: abstract still above 300 words")

# ================= F3 terminology: IgG titre =================
patch("Methods_draft.md",
      "and circulating anti-*H. pylori* IgG levels were used to estimate",
      "and circulating anti-*H. pylori* IgG titre was used to estimate", "F3 2.1#1 titre")
patch("Methods_draft.md",
      "**Circulating anti-*H. pylori* IgG levels** (quantitative; ieu-b-4905)",
      "**Circulating anti-*H. pylori* IgG titre** (quantitative; ieu-b-4905)", "F3 2.2.1 bullet titre")
patch("Methods_draft.md",
      "the *H. pylori* seropositivity, IgG-level and antigen-level GWAS",
      "the *H. pylori* seropositivity, IgG-titre and antigen-level GWAS", "F3 2.3 titre")
patch("Results_draft.md",
      "circulating anti-*H. pylori* IgG levels yielded 12 instruments",
      "circulating anti-*H. pylori* IgG titre yielded 12 instruments", "F3 3.1 titre")
patch("Results_draft.md",
      "When circulating anti-*H. pylori* IgG levels (a quantitative exposure) were used instead",
      "When circulating anti-*H. pylori* IgG titre (a quantitative exposure) was used instead", "F3 3.2 titre")
patch("Results_draft.md",
      "antibody titre per se (IgG level), the other five antigens",
      "antibody titre per se, the other five antigens", "F3 3.6 drop parenthetical")
patch("Tables_draft.md", "Anti-H. pylori IgG level |", "Anti-H. pylori IgG titre |", "F3 Table 1 titre")
patch("Tables_draft.md", "| IgG level |", "| IgG titre |", "F3 Table 2 titre", count=99)
patch("Tables_draft.md",
      "Seropositivity, IgG-level and antigen-level instruments",
      "Seropositivity, IgG-titre and antigen-level instruments", "F3 Table 1 note titre")
patch(os.path.join(CODE, "19_build_submission_md.py"),
      "circulating anti-*H. pylori* IgG level) on five",
      "circulating anti-*H. pylori* IgG titre) on five", "F3 Fig1 legend titre")

# ================= F5/F6 Methods cleanups =================
patch("Methods_draft.md", "and report it transparently as a limitation",
      "and report it as a limitation", "F5 drop 'transparently'")
patch("Methods_draft.md", "and ambiguous forward-strand ambiguities were resolved",
      "and forward-strand ambiguities were resolved", "F5 fix ambiguous/ambiguities")
patch("Methods_draft.md", "through OpenGWAS (eqtl-a-ENSG\u2026)", "through OpenGWAS", "F6 drop placeholder")

# ================= F4/F7 Results cleanups =================
patch("Results_draft.md",
      "\u22480.29 (walking pace), so the nulls exclude only appreciable reverse effects and should be read as "
      "absence of evidence rather than evidence of absence, with the walking-pace",
      "\u22480.29 (walking pace). The nulls therefore exclude only appreciable reverse effects and should be "
      "read as absence of evidence rather than evidence of absence, with the walking-pace",
      "F4 split 3.5 closing sentence")
patch("Results_draft.md", "is therefore the appropriately conservative primary estimate",
      "is therefore the conservative primary estimate", "F5 drop 'appropriately'")
patch("Results_draft.md", "supporting that the signal is not driven by a single variant",
      "indicating that the signal is not driven by a single variant", "F5 'supporting that'->'indicating that'")

# ================= Discussion =================
patch("Discussion_draft.md",
      "Recent MR studies of sarcopenia-related traits cover lifestyle factors and vitamin D",
      "Recent MR studies of sarcopenia-related traits have covered lifestyle factors, vitamin D",
      "F-list grammar 4.2")

# ================= F1/F7/F8 build script =================
b19 = os.path.join(CODE, "19_build_submission_md.py")
s = read(b19)
if "draft-version heading guard" not in s:
    guard = ('        if re.match(r"^#.*\\u2014 Draft v", ln):  # draft-version heading guard\n'
             "            continue\n")
    anchor = '        if ln.startswith(("*Manuscript:", "*Drafting conventions", "*Generated by")):\n            continue\n'
    assert anchor in s, "strip_meta anchor not found"
    s = s.replace(anchor, anchor + guard, 1)
    report.append(("F1 strip draft-version headings", "OK", "19"))
else:
    report.append(("F1 strip draft-version headings", "SKIP (applied)", "19"))

if "ref note scrub" not in s:
    scrub = ('# ref note scrub: italic editorial annotations and the **Note:** block must not ship\n'
             '_clean = []\n'
             'for _ln in ref_lines:\n'
             '    _ln = re.sub(r"\\s*\\*\\*Note:\\*\\*.*$", "", _ln)\n'
             '    _ln = re.sub(r"\\s*\\*\\([^\\*]*\\)\\*", "", _ln)\n'
             '    if _ln.strip():\n'
             '        _clean.append(_ln.rstrip())\n'
             'ref_lines = _clean\n')
    anchor2 = 'refs_md = "\\n".join(ref_lines).strip()'
    assert anchor2 in s, "refs anchor not found"
    s = s.replace(anchor2, scrub + anchor2, 1)
    report.append(("F7 strip ref editorial notes", "OK", "19"))
else:
    report.append(("F7 strip ref editorial notes", "SKIP (applied)", "19"))

old_wc = 'abstract_plain = re.sub(r"\\^[\\d,]+\\^", "", re.sub(r"^#.*$", "", ab, flags=re.M))'
new_wc = ('abstract_plain = re.sub(r"^Mendelian randomization;.*$", "",\n'
          '              re.sub(r"\\^[\\d,]+\\^", "", re.sub(r"^#.*$", "", ab, flags=re.M)), flags=re.M)  # exclude keywords line')
if old_wc in s:
    s = s.replace(old_wc, new_wc, 1)
    report.append(("F8 abstract word count excludes keywords", "OK", "19"))
elif "exclude keywords line" in s:
    report.append(("F8 abstract word count excludes keywords", "SKIP (applied)", "19"))
else:
    report.append(("F8 abstract word count", "FAIL: anchor not found", "19"))
write(b19, s)

# ================= report =================
print("%-42s %-24s %s" % ("PATCH", "STATUS", "FILE"))
fails = 0
for label, status, fn in report:
    print("%-42s %-24s %s" % (label, status, fn))
    if status.startswith("FAIL"):
        fails += 1
print("\n%d patch(es) failed" % fails)
sys.exit(1 if fails else 0)
