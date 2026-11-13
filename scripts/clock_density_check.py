import sys, csv, math
from scipy import stats
inp, pre, clk, hr = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
tok = open(inp).read().split('\n')
# parse sections
lines = [l.split() for l in tok]
i = next(k for k, l in enumerate(lines) if l and l[0] == 'extant_nodes'); ne = int(lines[i][1])
ext = {lines[i + 1 + q][1]: int(lines[i + 1 + q][0]) for q in range(ne)}
i = next(k for k, l in enumerate(lines) if l and l[0] == 'branches'); K = int(lines[i][1])
br = {int(lines[i + 1 + q][0]): (int(lines[i + 1 + q][1]), int(lines[i + 1 + q][2])) for q in range(K)}
pb = {b: next((q for q, (a2, c2) in br.items() if c2 >= 0 and c2 == a), None) for b, (a, c) in br.items()}
log = list(csv.DictReader(open(pre + '.log'), delimiter='\t'))
dump = [l.rstrip('\n').split('\t') for l in open(pre + '.dump')]
byit = {r['Iteration']: r for r in log}
worst = 0; n = 0
for d in dump:
    it, m, rates, rp, sig = d[0], float(d[5]), [float(x) for x in d[6].split(',')], float(d[10]), float(d[-1])
    r = byit[it]; age = {v: float(r['age_' + nm]) for nm, v in ext.items()}
    t = lambda b: abs(age[br[b][0]] - (0.0 if br[b][1] < 0 else age[br[b][1]]))
    s = stats.expon(scale=0.1).logpdf(m)
    if clk != 'exp': s += stats.expon(scale=1 / hr).logpdf(sig)
    for b in range(K):
        x = rates[b]
        if clk == 'exp': s += stats.expon(scale=m).logpdf(x)
        elif clk == 'ln': s += stats.lognorm(s=sig, scale=math.exp(math.log(m) - sig * sig / 2)).logpdf(x)
        elif clk == 'gam': s += stats.gamma(a=sig, scale=m / sig).logpdf(x)
        elif clk == 'igr': tt = t(b) + 0.01; s += stats.gamma(a=tt / sig, scale=sig * m / tt).logpdf(x)
        else:
            tt = t(b) + 0.01; pr = m if pb[b] is None else rates[pb[b]]
            s += stats.lognorm(s=sig * math.sqrt(tt), scale=math.exp(math.log(pr) - sig * sig * tt / 2)).logpdf(x)
    worst = max(worst, abs(s - rp)); n += 1
print(f'{clk}: {n} states, {K} branches, root branches {sum(v is None for v in pb.values())}, largest |ffbd - scipy| = {worst:.2e}')
