import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
X = json.load(open('h1b_extra.json')); B = json.load(open('h1b_results.json'))
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8), dpi=90)
f = np.array(X['T6']); ax[0].plot(f[:, 2], f[:, 1], 'o-'); ax[0].set_xlabel('peak rate R/r'); ax[0].set_ylabel('excess events [background-days]')
ax[0].set_title('Convex Pareto front: −25% peak for +3% events', fontsize=10)
r = X['T3b']; ax[1].plot([100*x['B'] for x in r], [-100*x['opt'] for x in r], 'o-', label='optimal flowback (convex)')
ax[1].plot([100*x['B'] for x in r], [-100*x['immediate'] for x in r], 'x--', label='immediate, max rate')
ax[1].set_xlabel('flowback volume [% of injected]'); ax[1].set_ylabel('event reduction [%]'); ax[1].legend(fontsize=8)
ax[1].set_title('Optimal shut-in = immediate flowback', fontsize=10)
t4 = B['T4']; ax[2].bar(['step-up better\n(co-injection)', 'step-up better\n(incl. tail)', 'ranking flips'],
       [100*t4['co_better'], 100*t4['tail_better'], 100*t4['flip']], color=['C2', 'C3', 'C7'])
ax[2].set_ylabel(f"% of {t4['n']} random sites"); ax[2].set_title('Window artefact across parameter space', fontsize=10)
fig.tight_layout(); fig.savefig('h1b_fig.png')
