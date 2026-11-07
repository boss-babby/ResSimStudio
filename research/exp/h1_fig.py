import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R = {r['tag']: r for r in json.load(open('h1_results.json'))}
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8), dpi=90)
names = ['step-up(4)', 'ramp-up-linear', 'step-down(4)', 'taper-last30%', 'cyclic(6x,50%)', 'OPTIMAL']
tags = ['Th=1Ti', 'Th=2Ti', 'baseline Th=5Ti', 'Th=20Ti']; lab = ['co-injection\nonly', '+1×Ti tail', '+4×Ti tail', '+19×Ti tail']
x = np.arange(len(tags)); wdt = 0.13
for i, n in enumerate(names):
    v = [100*(R[t]['rows'][n]['E_d']/R[t]['rows']['constant']['E_d']-1) for t in tags]
    ax[0].bar(x+(i-2.5)*wdt, v, wdt, label=n)
ax[0].axhline(0, c='k', lw=.6); ax[0].set_xticks(x, lab, fontsize=8); ax[0].set_ylabel('induced events vs constant rate [%]')
ax[0].set_title('Ranking flips with evaluation window', fontsize=10); ax[0].legend(fontsize=7)
for t, c in [('Th=1Ti', 'C3'), ('baseline Th=5Ti', 'C0')]:
    ax[1].step(R[t]['t_d'], np.array(R[t]['q_opt'])*1e3, where='mid', c=c, label=f'optimal, {t}')
ax[1].axhline(8640/(10*86400)*1e3, c='k', ls='--', lw=.8, label='constant')
ax[1].set_xlabel('time [d]'); ax[1].set_ylabel('rate [L/s]'); ax[1].set_title('Global optima (bang–singular–bang)', fontsize=10); ax[1].legend(fontsize=7)
Pi = [0.19, 0.38, 0.76, 1.53, 2.55, 3.82, 7.64]; gain = [0.6, 1.0, 2.0, 4.7, 9.8, 18.4, 46.7]; worst = [0.7, 1.0, 2.1, 6.9, 17.9, 42.2, 176.1]
ax[2].loglog(Pi, gain, 'o-', label='optimal gain vs constant'); ax[2].loglog(Pi, worst, 's-', label='worst standard schedule penalty')
ax[2].set_xlabel(r'$\Pi=\Delta S_{max}/(A\sigma)$'); ax[2].set_ylabel('%'); ax[2].set_title('Schedule matters only when Π ≳ 1', fontsize=10); ax[2].legend(fontsize=7)
fig.tight_layout(); fig.savefig('h1_fig.png')
