# -*- coding: utf-8 -*-
"""39_verify_supervisor_edits.py -- readback verification of the final docx."""
import zipfile, re, sys

z = zipfile.ZipFile(r"D:/MRout/70-MR/submission/manuscript_final.docx")
xml = z.read("word/document.xml").decode("utf-8")
text = re.sub(r"<[^>]+>", "", xml)

def chk(label, ok): print(("PASS" if ok else "FAIL"), label); return 0 if ok else 1

fails = 0

# 1. title
fails += chk("new title present", "implicating NF-\u03baB signaling in appendicular lean mass" in text)
fails += chk("old 'mediated inflammation' gone", "mediated inflammation" not in text)

# 2. P2 two readings + mediation infeasibility
for phrase in [
    "two competing readings that the present design cannot discriminate",
    "Under the activation-marker reading",
    "Under the inhibitor reading",
    "evidence that the pathway participates in the phenotype",
    "not currently feasible",
    "targeted set of about 10,000 pre-selected trait-associated variants",
    "genome-wide trans-eQTL resource for whole blood is not yet available",
]:
    fails += chk("4.3: " + phrase[:48], phrase in text)
fails += chk("4.3 'independently implicate' gone", "independently implicate" not in text)
fails += chk("4.3 'squarely' gone", "squarely" not in text)

# 3. P3 limitations
fails += chk("4.5 shared serology sub-study", "same UK Biobank serology sub-study" in text)
fails += chk("4.5 triangulation weight reduced", "triangulation weight of the antigen line is correspondingly reduced" in text)
fails += chk("4.5 fat-mass-adjusted ALM", "fat-mass-adjusted lean mass rather than total ALM" in text)
fails += chk("4.5 not transferred unqualified", "should not be transferred unqualified" in text)

# 4. P4-P6 reverse wording + power
fails += chk("3.5 power 0.13/0.29", "\u22480.13 (grip strength) and \u22480.29 (walking" in text)
fails += chk("3.5 absence of evidence", "absence of evidence rather than evidence of absence" in text)
fails += chk("old 3.5 'we conclude ... unlikely' gone", "we conclude that reverse causation is unlikely" not in text)
fails += chk("3.6 'provided no evidence of inverse causation'", "provided no evidence of inverse causation" in text)
fails += chk("old 'argued against inverse causation' gone", "argued against inverse causation" not in text)
fails += chk("4.1 power qualifier", "80% power at \u03b1 = 0.05 only for log-OR effects of \u22480.13 or larger" in text)

# 5. frozen numbers untouched
for num in ["\u22120.061", "3.7 \u00d7 10\u207b\u2074", "\u22120.021", "3.1 \u00d7 10\u207b\u2078",
            "\u22120.044", "4.9 \u00d7 10\u207b\u2074", "2.3 \u00d7 10\u207b\u2078", "7.6 \u00d7 10\u207b\u00b2\u2070",
            "\u22120.300", "0.0041", "0.0039", "0.0072", "27.34", "1.7 \u00d7 10\u207b\u2075",
            "0.05/24", "\u22120.040", "6.7 \u00d7 10\u207b\u2078", "0.039", "0.0010"]:
    fails += chk("frozen: " + num, num in text)

# 6. no unresolved markers / no literal superscript stars
fails += chk("no {{CITE}} leftovers", "{{CITE" not in text)
fails += chk("no TO CONFIRM", "TO CONFIRM" not in text)

# 7. citation resolution of the two new author-year cites
# (docx superscripts appear as bare numbers in tag-stripped text)
fails += chk("Vosa cite resolved to 15",
             "pre-selected trait-associated variants 15 and the cis-restricted" in text)
fails += chk("Butler-Laporte cite resolved to 13 (Limitations)",
             "same UK Biobank serology sub-study 13, so the directional agreement" in text)

# 8. reference count
nrefs = len(re.findall(r"\^", text))  # rough; skip
print("\n%s" % ("ALL CHECKS PASSED" if fails == 0 else "%d CHECK(S) FAILED" % fails))
sys.exit(1 if fails else 0)
