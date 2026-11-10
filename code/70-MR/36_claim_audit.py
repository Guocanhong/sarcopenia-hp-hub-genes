# -*- coding: utf-8 -*-
"""36_claim_audit.py — evidence audit of evaluative claims (2026-11-07).

User instruction: check claim-type assertions (significantly-improved /
effectively / fully-proves / excellent-performance / highly-novel class),
keep those backed by data or literature, convert or delete the rest.

Audit outcome per scanned hit:
  KEEP (backed):  all "robust" usages (anchored to the Methods 2.4 decision
    rule + frozen P values); "confirmed" for Steiger/LOO/units (analyses
    actually run, all TRUE/null outcomes recorded); "confirmed sarcopenia"
    (outcome term of Wu & Chen 2021); "canonical ... mimicry immunity"
    (Mayr 2000; Ranford 2002); "canonical signalling route ... MuRF1"
    (Cai 2004, MIKK mice); "substantial sequence homology" (Ranford 2002);
    "marked regional variation ... substantially higher in Asia"
    (Li 2023; Hooi 2017); "considerably stronger" (P 7.6e-20 vs 4.9e-4);
    "substantial heterogeneity" (Q = 27.34); "most notably/prominently"
    (NHANES III is the largest study of the association); "first MR study"
    ("to our knowledge" + Liu & Peng 2025 review + Nov 2026 search).
  FIX (no direct support): 6 items below.
"""
import io, os

P = r'D:/GEO-limma,肌少症1/70-MR/paper/'
log = []

def patch(fn, old, new, tag):
    p = os.path.join(P, fn)
    t = io.open(p, encoding='utf-8').read()
    if old in t:
        assert t.count(old) == 1, 'AMBIGUOUS %s [%s]' % (fn, tag)
        io.open(p, 'w', encoding='utf-8').write(t.replace(old, new))
        log.append('%s [%s] APPLIED' % (fn, tag))
    elif new in t:
        log.append('%s [%s] already applied' % (fn, tag))
    else:
        raise AssertionError('NOT FOUND %s [%s] :: %s...' % (fn, tag, old[:80]))

# P1 Intro ¶1: "substantially elevated risk" — consensus documents elevated
# risk; "substantially" is an unsupported intensifier.
patch('Introduction_draft.md',
 'confers a substantially elevated risk of falls, fractures, frailty, hospitalisation and all-cause mortality',
 'confers elevated risks of falls, fractures, frailty, hospitalisation and all-cause mortality',
 'P1-substantially')

# P2 Intro ¶2: "comprehensive" characterises someone else's review scope —
# cannot be verified here; the quoted conclusion itself carries the point.
patch('Introduction_draft.md',
 'A recent comprehensive review concluded, however,',
 'A recent review concluded, however,',
 'P2-comprehensive')

# P3 Discussion 4.2: same "comprehensive review" characterisation.
patch('Discussion_draft.md',
 'a comprehensive 2025 review of the *H. pylori*–muscle literature',
 'a 2025 review of the *H. pylori*–muscle literature',
 'P3-comprehensive')

# P4 Discussion 4.3: "established trigger" for the mitochondria→catabolism
# link outruns the adjacent citations (Singh 2024 = HSP60 perturbation
# → proteome loss/apoptosis; not a consensus-establishing source).
patch('Discussion_draft.md',
 'Mitochondrial dysfunction in muscle is an established trigger of the ubiquitin–proteasome catabolic programme',
 'Mitochondrial dysfunction in muscle is a recognised trigger of the ubiquitin–proteasome catabolic programme',
 'P4-established')

# P5 Discussion 4.6: "biologically tractable target for pharmacological
# modulation" — no citation, implies druggability evidence the paper does
# not contain. The NFKBIA finding is already stated in 4.6 sentence 1;
# delete the clause.
patch('Discussion_draft.md',
 'offers a testable mechanism linking chronic infection to tissue-specific catabolism, and NFKBIA-centred NF-κB signalling offers a biologically tractable target for pharmacological modulation.',
 'offers a testable mechanism linking chronic infection to tissue-specific catabolism.',
 'P5-pharmacological')

# P6 Methods 2.1 reverse-MR bullet: a design aim cannot promise exclusion
# ("exclude"); Results 3.5 correctly reports "no robust evidence".
patch('Methods_draft.md',
 'were used to estimate causal effects on *H. pylori* infection outcomes, to assess and exclude reverse causation.',
 'were used to estimate causal effects on *H. pylori* infection outcomes, to assess reverse causation.',
 'P6-exclude')

print('\n'.join(log))
print('%d patches processed' % len(log))
