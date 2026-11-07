import sys, csv, bisect
sys.path.insert(0,'/home/user/ressim/tools/dev')
sys.argv=[sys.argv[0]]+sys.argv[1:]
from smrydiff import load
ref, ours, keys = sys.argv[1], sys.argv[2], sys.argv[3].split(',')
kr, rr = load(ref)
rows = list(csv.DictReader(open(ours)))
t_o = [float(r['TIME']) for r in rows]
tr = [r[kr.index('TIME')] for r in rr]
for k in keys:
    j = kr.index(k)
    print('%-12s %8s %12s %12s %8s' % (k, 't', 'ref', 'ours', 'rel%'))
    for t in [x for x in t_o if x > 0][::max(1, len(t_o)//14)]:
        i = min(range(len(tr)), key=lambda q: abs(tr[q]-t))
        o = float(rows[t_o.index(t)][k]); r = rr[i][j]
        print('%-12s %8.1f %12.5g %12.5g %8.2f' % ('', t, r, o, 100*(o-r)/r if r else 0))
