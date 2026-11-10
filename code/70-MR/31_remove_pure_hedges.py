# -*- coding: utf-8 -*-
"""31_remove_pure_hedges.py — 删除只规避争议、不贡献论证的语句/短语。
2026-11-07，用户指令："如果某一语句只用来规避争议，而无法对论点的推导论证起到贡献，予以删除"。

三处删除（均不含引用，无重编号连锁）：
1. Results 3.5："; the nominal signal should therefore be interpreted with great caution"
   - 方向性 FALSE 的三条证据（Steiger FALSE、N=8,735 小样本、R² 低）已在句中；
     最终裁定由 3.5 末段（WM P=0.054 不支持 + Steiger FALSE + 反向因果不成立）
     与 Methods 2.4 统一决策规则给出，"great caution" 无推导增量。
2. Discussion 4.6："; extending this association map to muscle mass is a natural next
   step rather than an outlier claim"
   - 与前一分句 "extending the same design to muscle phenotypes is the contribution
     of the present study" 语义重复；"rather than an outlier claim" 纯防御。
3. Discussion 4.5（Finally 段）：" rather than as an agnostic discovery screen"
   - "hypothesis-generating" 的定性已承担全部认识论定位；
     "no one claimed agnostic discovery screen"，对照物不存在，纯规避。

评估后保留（对推导有实质贡献）的 hedging：
- Abstract "If confirmed, screening and eradication may merit evaluation"（转化含义）
- 4.3 "A coherent, though partly hypothesised, mechanistic chain"（GroEL 降级后的
  认识论框架，标记哪些环节是假设——评审要求）
- 4.3 "If the suggestive GroEL signal is confirmed, ..."（条件式，负载荷）
- 4.3 NFKBIA 替代解释句与 GTEx nominal 证据（方向性歧义的实质处理）
- 4.4 "nulls should be read as 'no evidence of an effect', not as evidence of no
  effect"（null 的正确统计解读，直接服务特异性论证）
- 4.5 Limitations 全部条目（含 BIA 水合混杂——指明了具体替代机制）
- Methods 2.4 统一决策规则（预设定规则，是裁定依据本身）
"""
import io

P = 'D:/GEO-limma,肌少症1/70-MR/paper/'

def patch(fn, pairs):
    p = P + fn
    t = io.open(p, encoding='utf-8').read()
    for old, new in pairs:
        c = t.count(old)
        assert c == 1, 'count=%d for %r in %s' % (c, old[:60], fn)
        t = t.replace(old, new)
    io.open(p, 'w', encoding='utf-8').write(t)
    print(fn, 'patched', len(pairs))

# 1) Results 3.5 — delete pure-caution clause
patch('Results_draft.md', [
    (' and the low variance explained; the nominal signal should therefore be '
     'interpreted with great caution.',
     ' and the low variance explained.'),
])

# 2) Discussion 4.6 — delete duplicate defensive clause
patch('Discussion_draft.md', [
    ('extending the same design to muscle phenotypes is the contribution of the '
     'present study; extending this association map to muscle mass is a natural '
     'next step rather than an outlier claim.',
     'extending the same design to muscle phenotypes is the contribution of the '
     'present study.'),
])

# 3) Discussion 4.5 — trim defensive tail, keep the load-bearing designation
patch('Discussion_draft.md', [
    ('and the results at this layer should be read as hypothesis-generating '
     'rather than as an agnostic discovery screen.',
     'and the results at this layer should be read as hypothesis-generating.'),
])

print('all patches applied')
