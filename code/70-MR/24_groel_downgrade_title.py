# -*- coding: utf-8 -*-
import io, os
old = 'implicating the GroEL\u2013NFKBIA inflammatory axis in appendicular lean mass'
new = 'implicating NF-\u03baB\u2013mediated inflammation in appendicular lean mass'
base = 'D:/GEO-limma,\u808c\u5c11\u75c71/70-MR/paper'
files = ['Abstract_draft.md', 'Introduction_draft.md', 'Methods_draft.md',
         'Results_draft.md', 'Discussion_draft.md', 'References_draft.md',
         'Framework_outline.md']
for fn in files:
    p = os.path.join(base, fn)
    t = io.open(p, encoding='utf-8').read()
    n = t.count(old)
    if n:
        io.open(p, 'w', encoding='utf-8').write(t.replace(old, new))
    print(fn, n)
