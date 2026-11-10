# -*- coding: utf-8 -*-
"""Re-apply edits lost to parallel-write races (GroEL downgrade, 2026-11-07)."""
import io, os
BASE = 'D:/GEO-limma,肌少症1/70-MR/paper'

def patch(fn, pairs):
    p = os.path.join(BASE, fn)
    t = io.open(p, encoding='utf-8').read()
    for old, new in pairs:
        n = t.count(old)
        if n == 1:
            t = t.replace(old, new)
        print('  [%s] %s -> count=%d' % (fn, old[:60].replace('\n', ' '), n))
    io.open(p, 'w', encoding='utf-8').write(t)

# ---------- Abstract ----------
patch('Abstract_draft.md', [
 ("with Bonferroni correction (\u03b1 = 0.05/24 within the antigen line; \u03b1 = 0.01 for NFKBIA). A reverse-MR line tested",
  "with Bonferroni correction (\u03b1 = 0.05/24 within the antigen line; \u03b1 = 0.01 for NFKBIA) and a uniform decision rule that signals supported only by the random-effects model are reported as suggestive. A reverse-MR line tested"),
 ("Of the six antigens, only anti-GroEL antibodies lowered ALM (\u03b2 = \u22120.021 SD, P = 3.1 \u00d7 10\u207b\u2078) \u2014 the sole signal surviving Bonferroni correction.",
  "Of the six antigens, anti-GroEL antibody levels showed a nominally significant random-effects association with lower ALM (\u03b2 = \u22120.021 SD, P = 3.1 \u00d7 10\u207b\u2078), which was not supported by the fixed-effect or weighted-median estimators and is therefore reported as suggestive under the uniform decision rule."),
])

# ---------- Results ----------
patch('Results_draft.md', [
 ("One additional signal approached but did not robustly survive correction: genetic liability to higher anti-OMP antibody levels was associated with higher odds of low grip strength in the random-effects model (log-OR = 0.039, SE = 0.012, P = 0.0010, 9 instruments), which nominally passes the Bonferroni threshold; however, the fixed-effect estimate was weaker (P = 0.028) and the weighted-median and MR-Egger estimates were null (P = 0.12 and P = 0.34 respectively). We therefore report the OMP signal as suggestive only.",
  "The same decision rule applies to the OMP antigen: genetic liability to higher anti-OMP antibody levels was associated with higher odds of low grip strength in the random-effects model (log-OR = 0.039, SE = 0.012, P = 0.0010, 9 instruments), which nominally passes the Bonferroni threshold; however, the fixed-effect estimate was weaker (P = 0.028) and the weighted-median and MR-Egger estimates were null (P = 0.12 and P = 0.34 respectively). The OMP signal is likewise reported as suggestive only."),
 ("Taken together, the antigen-spectrum analysis localises the forward-line infection signal to the anti-GroEL (chaperonin-60/HSP60) immune response: the general seropositivity signal and the GroEL-specific signal converge on ALM, while antibody responses to the other five antigens do not reproduce it.",
  "Taken together, the antigen-spectrum analysis provides suggestive, directionally concordant evidence that the forward-line infection signal may be carried in part by the anti-GroEL (chaperonin-60/HSP60) immune response: the general seropositivity signal and the GroEL-specific estimate converge on ALM in sign and outcome, while antibody responses to the other five antigens do not reproduce it. With only four harmonised GroEL instruments and no estimator corroboration beyond the random-effects model, this antigen-level localisation remains provisional (Discussion 4.3)."),
 ("  4. GroEL\u2192walking FE P = 0.0013 is a positive-direction point estimate inconsistent with ALM \u2014 interpreted as null (wording above).",
  "  4. GroEL\u2192walking FE P = 0.0013 is a positive-direction point estimate inconsistent with ALM \u2014 interpreted as null (wording above).\n  5. **GroEL downgrade (2026-11-07, post-review):** GroEL\u2192ALM is RE-only supported (RE SE 0.0039 < FE SE 0.0072 at Q P = 0.83; FE P = 0.003 > 0.0021; WM P = 0.017) and is now reported as suggestive under the uniform Methods 2.4 rule \u2014 same treatment as OMP and CTNNBIP1. Title and Abstract updated accordingly (\"NF-\u03baB\u2013mediated inflammation\"); 3.6 and Discussion 4.1\u20134.6 synced.\n  6. **Reverse Bonferroni corrected (2026-11-07):** walking\u2192seropositivity P = 0.0041 DOES pass \u03b1 = 0.0083; 3.5 now states this explicitly with WM P = 0.054 and no-Steiger caveat (Table 5 walking rows have no Steiger entry; recomputation pending)."),
])

# ---------- Methods ----------
patch('Methods_draft.md', [
 ("(range 21.1\u2013857.2 across all instrument sets; per-set ranges in Table 1)",
  "(range 20.9\u2013857.2 across all instrument sets; per-set ranges in Table 1)"),
])

# ---------- Discussion ----------
patch('Discussion_draft.md', [
 ("indicating that the signal is outcome-specific (ALM), antigen-specific (GroEL) and pathway-specific (NFKBIA).",
  "indicating that the robust signals are outcome-specific (ALM) and pathway-specific (NFKBIA), with the antigen-level GroEL estimate providing suggestive directional concordance."),
 ("Among six *H. pylori* antigens, only the GroEL-specific immune response produced a robust causal signal.",
  "Among six *H. pylori* antigens, only the GroEL-specific immune response produced a nominally significant random-effects signal, directionally concordant with the general serology."),
 ("## 4.3 Mechanistic interpretation: a GroEL\u2013HSP60\u2013NF-\u03baB axis for muscle loss\n\nA coherent mechanistic chain integrates the three convergent signals.",
  "## 4.3 Mechanistic interpretation: a proposed GroEL\u2013HSP60\u2013NF-\u03baB pathway for muscle loss\n\nA coherent, though partly hypothesised, mechanistic chain integrates the convergent signals."),
 ("Within the limits of a two-sample MR design, the convergent genetic evidence supports a causal pathway from *H. pylori* infection \u2014 specifically the anti-GroEL immune response \u2014 through an NF-\u03baB\u2013centric inflammatory programme to reduced appendicular lean mass.",
  "Within the limits of a two-sample MR design, the genetic evidence supports a causal effect of *H. pylori* infection on reduced appendicular lean mass, with causal support for an NF-\u03baB\u2013centric inflammatory programme (NFKBIA) and suggestive, directionally concordant evidence implicating the anti-GroEL immune response."),
])
print('ALL DONE')
