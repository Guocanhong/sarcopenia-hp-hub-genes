# -*- coding: utf-8 -*-
"""38_supervisor_edits.py -- Supervisor-round edits (2026-11-07), 6 items:
P1  title: "NF-kB-mediated inflammation" -> "NF-kB signalling" (all headers + build script)
P2  4.3: two competing readings of NFKBIA + mediation-test infeasibility + drop "independently"
P3a 4.5 Limitations Seventh: serology lines share one sub-study (triangulation weight)
P3b 4.5 Limitations Eighth: ALM is fat-mass-adjusted (AFM residualisation boundary)
P4  3.5 ending: soften reverse-MR conclusion + power numbers (MDE = 2.8 x SE)
P5  3.6: "argued against inverse causation" -> "provided no evidence of inverse causation"
P6  4.1: same softening with power qualifier
Idempotent: skips patches whose NEW text is already present; reports any OLD text not found.
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
        report.append((label, "SKIP (already applied)", fname)); return
    n = s.count(old)
    if n == 0:
        report.append((label, "FAIL: old string not found", fname)); return
    s = s.replace(old, new, count)
    write(path, s)
    report.append((label, "OK (%d of %d occurrence(s))" % (count, n), fname))

# ---------- assertions on reference numbering ----------
refs = read(os.path.join(PAPER, "References_draft.md"))
assert re.search(r"^15\. .*V\u00f5sa.*34475573", refs, re.M), "ref 15 is not Vosa 34475573"
assert re.search(r"^13\. .*Butler-Laporte.*33204752", refs, re.M), "ref 13 is not Butler-Laporte 33204752"

# ---------- P1 title (draft headers, framework, references) ----------
OLD_TITLE = "implicating NF-\u03baB\u2013mediated inflammation in appendicular lean mass"
NEW_TITLE = "implicating NF-\u03baB signalling in appendicular lean mass"
for fn in ["Abstract_draft.md", "Introduction_draft.md", "Methods_draft.md",
           "Results_draft.md", "Discussion_draft.md", "References_draft.md",
           "Framework_outline.md"]:
    patch(fn, OLD_TITLE, NEW_TITLE, "P1 title", count=99)

# build script: TITLE tuple (literal \u escapes in source) + title-page heading line
b19 = os.path.join(CODE, "19_build_submission_md.py")
s = read(b19)
changed = False
tup_old = 'NF-\\u03baB\\u2013mediated "' + "\n" + '         "inflammation in appendicular lean mass'
tup_new = 'NF-\\u03baB signalling "' + "\n" + '         "in appendicular lean mass'
if tup_old in s:
    s = s.replace(tup_old, tup_new); changed = True
hdr_old = r'implicating NF-\u03baB\u2013mediated inflammation in appendicular lean mass'
hdr_new = r'implicating NF-\u03baB signalling in appendicular lean mass'
if hdr_old in s:
    s = s.replace(hdr_old, hdr_new); changed = True
if changed:
    write(b19, s); report.append(("P1 title (build script)", "OK", "19"))
elif hdr_new in s:
    report.append(("P1 title (build script)", "SKIP (already applied)", "19"))
else:
    report.append(("P1 title (build script)", "FAIL: neither old nor new found", "19"))

# ---------- P2 4.3 rewrite ----------
d_old = ("NFKBIA encodes I\u03baB\u03b1, the principal cytoplasmic inhibitor of NF-\u03baB and an NF-\u03baB target gene itself: "
 "transcriptional upregulation of NFKBIA is part of the pathway's autoregulatory feedback loop, so genetically "
 "predicted elevation of NFKBIA expression plausibly marks a state of heightened NF-\u03baB pathway activation rather "
 "than simple pathway suppression (Hoffmann et al., 2002) {{CITE 28: Hoffmann A, et al. Science 2002;298(5596):1241\u20131245, "
 "PMID 12424381 \u2014 I\u03baB\u03b1 strong negative-feedback module}}. The direction of our estimate (greater predicted NFKBIA "
 "expression, lower ALM) is therefore consistent with a catabolic, NF-\u03baB\u2013dominant muscle phenotype. We note, "
 "alternatively, that in myogenic contexts I\u03baB\u03b1 overexpression has been shown to block NF-\u03baB activation and "
 "ameliorate inflammatory atrophy (Cai et al., 2004) {{CITE 27: PMID 15479644 \u2014 I\u03baB\u03b1 super-repressor (MISR) mice; "
 "NF-\u03baB inhibition reduced denervation- and tumor-induced muscle loss and improved survival}}, and our eQTL "
 "instruments reflect whole-blood expression rather than muscle tissue; distinguishing \"NFKBIA as a marker of "
 "pathway activation\" from \"NFKBIA as a causal mediator within muscle\" will require tissue-resolved expression data.")
d_new = ("NFKBIA encodes I\u03baB\u03b1, the principal cytoplasmic inhibitor of NF-\u03baB and an NF-\u03baB target gene itself, and the "
 "signal admits two competing readings that the present design cannot discriminate. Under the activation-marker "
 "reading, transcriptional upregulation of NFKBIA is part of the pathway's autoregulatory feedback loop, so "
 "genetically predicted elevation of blood expression plausibly marks a state of heightened NF-\u03baB pathway "
 "activation rather than simple pathway suppression (Hoffmann et al., 2002) {{CITE 28: Hoffmann A, et al. Science "
 "2002;298(5596):1241\u20131245, PMID 12424381 \u2014 I\u03baB\u03b1 strong negative-feedback module}}; on this reading, the estimate "
 "(greater predicted expression, lower ALM) indexes a catabolic, NF-\u03baB\u2013dominant muscle phenotype. Under the "
 "inhibitor reading, I\u03baB\u03b1 acts on the pathway directly: in myogenic contexts, I\u03baB\u03b1 overexpression blocks "
 "NF-\u03baB activation and ameliorates inflammatory atrophy (Cai et al., 2004) {{CITE 27: PMID 15479644 \u2014 I\u03baB\u03b1 "
 "super-repressor (MISR) mice; NF-\u03baB inhibition reduced denervation- and tumor-induced muscle loss and improved "
 "survival}}, so a deleterious effect of higher predicted expression is not expected on simple inhibitory grounds. "
 "The gene-expression line therefore establishes that genetically predicted NFKBIA expression is causally "
 "associated with ALM \u2014 evidence that the pathway participates in the phenotype \u2014 but does not establish which "
 "direction of NF-\u03baB signalling the estimate reflects, nor does it test the link from infection to pathway "
 "activation itself: a proximate mediation test (seropositivity instruments against NFKBIA expression as the "
 "outcome) is not currently feasible, because the phase-I eQTLGen resource tested trans effects only for a "
 "targeted set of about 10,000 pre-selected trait-associated variants (V\u00f5sa et al., 2021) {{CITE 15: V\u00f5sa U, et al. "
 "Nat Genet 2021;53(9):1300\u20131310, PMID 34475573}} and the cis-restricted outcome files used here cannot accommodate "
 "genome-wide seropositivity instruments; a genome-wide trans-eQTL resource for whole blood is not yet available. "
 "The instruments also reflect whole-blood expression rather than muscle tissue, so distinguishing \"NFKBIA as a "
 "marker of pathway activation\" from \"NFKBIA as a causal mediator within muscle\" will require tissue-resolved "
 "expression data.")
patch("Discussion_draft.md", d_old, d_new, "P2 4.3 two readings")

d2_old = ("Either interpretation, however, places the causal signal squarely on the inflammatory signalling axis "
          "that the infection and antigen lines independently implicate.")
d2_new = ("Either interpretation, however, places the causal signal on the inflammatory signalling axis that the "
          "infection and antigen lines implicate.")
patch("Discussion_draft.md", d2_old, d2_new, "P2 4.3 closing sentence")

# ---------- P3a Limitations: shared serology sub-study ----------
l7_old = ("and the bias from overlap is bounded by the instruments' small variance explained (Burgess et al., 2016) "
 "{{CITE 21: Burgess S, Davies NM, Thompson SG. Genet Epidemiol 2016;40(7):597\u2013608, DOI 10.1002/gepi.21998 \u2014 bias "
 "proportional to overlap fraction \u00d7 1/F}}.")
l7_new = l7_old + (" Relatedly, the seropositivity and antigen-spectrum exposure GWAS derive from the same UK Biobank "
 "serology sub-study (Butler-Laporte et al., 2020) {{CITE 13: Butler-Laporte G, et al. Nat Commun 2020;11:5371, "
 "PMID 33204752}}, so the directional agreement between the two serology lines reflects phenotypic decomposition "
 "within one sample rather than fully independent replication, and the triangulation weight of the antigen line "
 "is correspondingly reduced.")
patch("Discussion_draft.md", l7_old, l7_new, "P3a shared serology sub-study")

# ---------- P3b Limitations: AFM residualisation boundary ----------
l8_old = "replication against DXA-derived lean-mass GWAS as these resources scale is a priority."
l8_new = l8_old + (" In addition, the source GWAS residualised ALM for appendicular fat mass (Methods 2.2.3), so the "
 "analysed phenotype is fat-mass-adjusted lean mass rather than total ALM; because fat mass is itself influenced "
 "by inflammatory and metabolic exposures, the estimates should not be transferred unqualified to unadjusted ALM "
 "or to the syndromic sarcopenia definitions used in observational studies.")
patch("Discussion_draft.md", l8_old, l8_new, "P3b AFM residualisation boundary")

# ---------- P4 3.5 ending ----------
r_old = ("Weighing all sensitivity evidence, we conclude that reverse causation is unlikely to explain the "
 "forward-line findings, while acknowledging the walking-pace \u2192 seropositivity signal as the residual uncertainty "
 "in this assessment.")
r_new = ("Weighing all sensitivity evidence, the reverse line provides no evidence that reverse causation explains "
 "the forward-line findings, but its power is limited: at the standard errors of the primary random-effects "
 "estimates, 80% power at \u03b1 = 0.05 is reached only for log-OR effects of \u22480.13 (grip strength) and \u22480.29 (walking "
 "pace), so the nulls exclude only appreciable reverse effects and should be read as absence of evidence rather "
 "than evidence of absence, with the walking-pace \u2192 seropositivity signal remaining the residual uncertainty in "
 "this assessment.")
patch("Results_draft.md", r_old, r_new, "P4 3.5 reverse conclusion + power")

# ---------- P5 3.6 ----------
r2_old = "and reverse MR argued against inverse causation."
r2_new = "and reverse MR provided no evidence of inverse causation."
patch("Results_draft.md", r2_old, r2_new, "P5 3.6 wording")

# ---------- P6 4.1 ----------
d3_old = ("Reverse MR provided no evidence that muscle phenotypes causally influence infection status, arguing "
          "against inverse causation.")
d3_new = ("Reverse MR provided no evidence that muscle phenotypes causally influence infection status; given the "
 "limited power of that line (80% power at \u03b1 = 0.05 only for log-OR effects of \u22480.13 or larger), it excludes "
 "appreciable reverse effects and should be read as absence of evidence rather than evidence of absence.")
patch("Discussion_draft.md", d3_old, d3_new, "P6 4.1 wording + power")

# ---------- report ----------
print("%-38s %-34s %s" % ("PATCH", "STATUS", "FILE"))
fails = 0
for label, status, fn in report:
    print("%-38s %-34s %s" % (label, status, fn))
    if status.startswith("FAIL"):
        fails += 1
print("\n%d patch(es) failed" % fails)
sys.exit(1 if fails else 0)
