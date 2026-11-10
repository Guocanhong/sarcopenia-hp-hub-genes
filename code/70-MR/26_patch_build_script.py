# -*- coding: utf-8 -*-
"""Update 19_build_submission_md.py (title, H1, fig legends, Table 5 note).
Source file stores literal \\uXXXX escapes, so match with raw strings."""
import io
p = '19_build_submission_md.py'
t = io.open(p, encoding='utf-8').read()

pairs = [
 (r'implicating the GroEL\u2013NFKBIA "' + '\n         ' + r'"inflammatory axis in appendicular lean mass")',
  r'implicating NF-\u03baB\u2013mediated "' + '\n         ' + r'"inflammation in appendicular lean mass")'),
 (r'# Causal effects of *Helicobacter pylori* infection on sarcopenia-related traits: a two-sample Mendelian randomization study implicating the GroEL\u2013NFKBIA inflammatory axis in appendicular lean mass',
  r'# Causal effects of *Helicobacter pylori* infection on sarcopenia-related traits: a two-sample Mendelian randomization study implicating NF-\u03baB\u2013mediated inflammation in appendicular lean mass'),
 (r'Both nominal signals (grip \u2192 seropositivity; walking pace \u2192 seropositivity) failed sensitivity checks and neither survives Bonferroni correction across the six reverse tests (\u03b1 = 0.0083). Steiger direction could not be confirmed for grip \u2192 seropositivity (FALSE).',
  r'The grip \u2192 seropositivity signal failed the Steiger directionality check (FALSE); the walking pace \u2192 seropositivity signal survives Bonferroni correction across the six reverse tests (\u03b1 = 0.0083) but is unsupported by the weighted-median estimator (P = 0.054), absent for self-reported infection, and lacks a computable directionality test for the walking-pace pairs (Results 3.5).'),
 (r'The red square marks the GroEL \u2192 ALM signal, the only estimate surviving correction (P = 3.1 \u00d7 10\u207b\u2078).',
  r'The red square marks the GroEL \u2192 ALM signal, the only antigen estimate nominally passing the Bonferroni threshold in the random-effects model (P = 3.1 \u00d7 10\u207b\u2078); under the uniform decision rule (Methods 2.4) it is reported as suggestive (Results 3.3).'),
 (r'Nominal signals (grip \u2192 seropositivity, P = 0.010; walking pace \u2192 seropositivity, P = 0.0041) failed Steiger/sensitivity checks and do not survive Bonferroni correction (Results 3.5).',
  r'Nominal signals (grip \u2192 seropositivity, P = 0.010; walking pace \u2192 seropositivity, P = 0.0041) failed Steiger or weighted-median corroboration (Results 3.5).'),
]
ok = True
for old, new in pairs:
    n = t.count(old)
    print('count=%d :: %s' % (n, old[:60]))
    if n != 1:
        ok = False
    t = t.replace(old, new)
io.open(p, 'w', encoding='utf-8').write(t)
print('BUILD SCRIPT UPDATED' if ok else 'SOME PATTERNS NOT FOUND')
