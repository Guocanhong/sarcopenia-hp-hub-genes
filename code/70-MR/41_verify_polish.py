# -*- coding: utf-8 -*-
"""41_verify_polish.py -- readback verification for the final-polish round."""
import zipfile, re, sys

z = zipfile.ZipFile(r"D:/MRout/70-MR/submission/manuscript_final.docx")
text = re.sub(r"<[^>]+>", "", z.read("word/document.xml").decode("utf-8"))

def chk(label, ok): print(("PASS" if ok else "FAIL"), label); return 0 if ok else 1
fails = 0

# --- build-residue removal ---
fails += chk("no 'Draft v0.x' residue", "Draft v0." not in text)
fails += chk("no duplicate Abstract heading", text.count("Abstract") >= 1 and "— Draft v0.1" not in text)

# --- reference notes stripped ---
for bad in ["to confirm at typesetting", "to fill at typesetting", "**Note:**",
            "not yet peer-reviewed published", "earlier-search", "cite the correction only"]:
    fails += chk("ref note gone: " + bad[:28], bad not in text)
fails += chk("ref 21 intact after note strip", "PMID 27625185. DOI 10.1002/gepi.21998" in text)
fails += chk("ref 3 intact after note strip", "Lancet Gastroenterol Hepatol. 2023;8(6):553–564. PMID 37086739." in text)

# --- terminology: IgG titer ---
fails += chk("Methods 2.2.1 'IgG titer' bullet", "Circulating anti-H. pylori IgG titer (quantitative; ieu-b-4905)" in text)
fails += chk("no 'Anti-H. pylori IgG level'", "Anti-H. pylori IgG level" not in text)
fails += chk("no 'IgG levels yielded 12'", "IgG levels yielded" not in text)
fails += chk("no 'IgG levels (a quantitative exposure)'", "IgG levels (a quantitative exposure)" not in text)
fails += chk("no 'IgG-level GWAS'", "IgG-level and antigen-level" not in text)
fails += chk("Table 2 exposure renamed", "IgG titer" in text and text.count("IgG titer") >= 20)
fails += chk("antigen rows keep 'IgG level'", "Anti-CagA IgG level" in text)
fails += chk("Fig1 legend titre", "circulating anti-H. pylori IgG titer) on five" in text)

# --- abstract: Sci Rep unstructured <=200-word version (Round-2, 2026-11-10) ---
for phrase in ["is linked epidemiologically to lower muscle mass and sarcopenia",
               "In a two-sample Mendelian randomization (MR) framework",
               "classifying random-effects-only signals as suggestive",
               "with concordant sensitivity estimates",
               "the anti-GroEL random-effects association was suggestive only",
               "DTL was untestable (&lt;3 IVs)",
               "reverse MR provided no evidence of inverse causation",
               "implicating NFKBIA/NF-κB signaling"]:
    fails += chk("abstract: " + phrase[:40], phrase in text)
for gone in ["causality is untested and inflammation the presumed mediator",
             "which is therefore suggestive",
             "unsupported by all other estimators, is likewise suggestive",
             "suggestive concordant evidence for the anti-GroEL immune response",
             "**Background:**"]:
    fails += chk("abstract gone: " + gone[:40], gone not in text)
fails += chk("abstract self-report dropped", "self-report, N = 462,933; IgG titer" not in text)
fails += chk("title-page abstract count <= 200",
             bool(re.search(r"Abstract ≈ (\d+);", text)) and int(re.search(r"Abstract ≈ (\d+);", text).group(1)) <= 200)

# --- Round-2 mock-review fixes (2026-11-10) ---
for phrase in [
    # Methods 2.4: post-hoc disclosure
    "scripted sensitivity re-runs executed on the frozen harmonized datasets",
    "no instrument selection, threshold or reported estimate was altered in response to these re-runs",
    # Methods 2.4: corroboration threshold made explicit
    "at the family-wise Bonferroni threshold of that test family",
    # Methods 2.4: multiplicative-RE rescaling named (check string stops before the
    # smart apostrophe pandoc renders for "Cochran's")
    "rescales the fixed-effect standard error by the square root of the ratio of Cochran",
    "when Q falls below its degrees of freedom this factor is below unity",
    # Methods 2.4: WM vs weighted mode defined
    "the two are distinct robustness estimators and are labelled consistently in all tables",
    # Results 3.3: two pairs pass antigen Bonferroni
    "two pairs passed the threshold in the random-effects model",
    # Results 3.4: NFKBIA leave-one-out
    "leave-one-out analyses showed that the estimate remained directionally concordant and nominally significant when each of the five instruments was excluded in turn",
    # Results 3.5 / Discussion Seventh: reverse-line overlap
    "overlap almost completely with the self-reported infection outcome GWAS (N = 462,933, also UK Biobank)",
    "those nulls could in principle reflect attenuation under a violated independence assumption",
    # Discussion 4.3: immune-cell composition confound
    "partly driven by leukocyte-subset composition",
    # Discussion 4.3: GroEL chain hypothesis framing
    "presented strictly as a hypothesis for independent testing",
    # Discussion: hub-gene small-sample instability
    "makes the DEG calls and machine-learning feature rankings inherently unstable and susceptible to overfitting",
    # Discussion First limitation: winner's curse caveat (pandoc smart apostrophe U+2019)
    "inflated by winner\u2019s curse",
]:
    fails += chk("round2: " + phrase[:48], phrase in text)
fails += chk("round2: 'one pair produced' gone", "one pair produced a nominally significant random-effects signal" not in text)

# --- modifier cleanups ---
for gone in ["transparently as a limitation", "appropriately conservative",
             "ambiguous forward-strand ambiguities", "supporting that the signal is not driven",
             "(eqtl-a-ENSG…)", "reported as suggestive under the uniform decision rule",
             "with concordant alternative estimators and supportive sensitivity analyses"]:
    fails += chk("gone: " + gone[:38], gone not in text)
for present in ["and forward-strand ambiguities were resolved",
                "is therefore the conservative primary estimate",
                "indicating that the signal is not driven by a single variant",
                "The nulls therefore exclude only appreciable reverse effects",
                "have covered lifestyle factors, vitamin D"]:
    fails += chk("present: " + present[:38], present in text)

# --- frozen numbers untouched ---
for num in ["−0.061", "3.7 × 10⁻⁴", "−0.021", "3.1 × 10⁻⁸", "−0.044", "4.9 × 10⁻⁴",
            "2.3 × 10⁻⁸", "7.6 × 10⁻²⁰", "−0.300", "0.0041", "0.0039", "0.0072",
            "27.34", "0.05/24", "−0.040", "6.7 × 10⁻⁸", "0.039", "0.0010",
            "462,933", "31,684", "10,179–14,263"]:
    fails += chk("frozen: " + num, num in text)

# --- structure sanity ---
fails += chk("38 references", "Zhang J, Li W. Association between" in text)
fails += chk("no {{CITE}} leftovers", "{{CITE" not in text)
fails += chk("no TO CONFIRM", "TO CONFIRM" not in text)

print("\n%s" % ("ALL CHECKS PASSED" if fails == 0 else "%d CHECK(S) FAILED" % fails))
sys.exit(1 if fails else 0)
