import sys; sys.path.insert(0, '/home/user/research')
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from chairN import *
fr30 = {'0': (0, 1, 2), (0, 0, 0): (0, 1, 2), (0, 0, 1): (0, 2, 1), (0, 1, 0): (2, 1, 0), (0, 1, 1): (1, 2, 0),
        (1, 0, 0): (1, 0, 2), (1, 0, 1): (0, 1, 2), (1, 1, 0): (2, 0, 1)}
ch = Chair(3, fr30)
# ---- Figure: supertile 2C_3 with children (voxels)
fig = plt.figure(figsize=(11, 5))
ax = fig.add_subplot(1, 2, 1, projection='3d')
grid = np.zeros((4, 4, 4), dtype=int) - 1
keys = ['0'] + [v for v in ch.U]
cols = plt.cm.Set3(np.linspace(0, 1, 12))[[0,4,3,5,6,7,9,11]]
for idx, (h, tau) in enumerate(ch.children):
    for c in ch.cubes((h, tau)):
        u = tuple((x - 1) // 2 for x in c); grid[u] = idx
filled = grid >= 0
facecolors = np.empty(grid.shape + (4,), dtype=float)
for idx in range(8): facecolors[grid == idx] = cols[idx]
ax.voxels(filled, facecolors=facecolors, edgecolor='k', linewidth=0.4, shade=False)
ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z'); ax.set_title(r'Level-1 supertile $2C_3$: 8 children')
ax.view_init(elev=25, azim=40)
ax2 = fig.add_subplot(1, 2, 2); ax2.axis('off')
rows = [['child', 'position', 'sign flips $D_s$', r'permutation $\pi_v$', 'rotation $h_v = D_s\pi_v$']]
desc = {'0': 'identity (translation by (1,1,1))', (0,0,0): 'identity', (1,0,0): '180° about $e_2+e_3$', (0,1,0): '180° about $e_1+e_3$',
        (0,0,1): '180° about $e_1+e_2$', (1,0,1): '180° about $e_2$', (0,1,1): '120° about $(1,1,-1)$', (1,1,0): '120° about $(1,-1,-1)$'}
for idx, k in enumerate(keys):
    if k == '0': rows.append(['$T_0$', '(1,1,1)', 'none', '(0,1,2)', desc[k]])
    else:
        s = ''.join('-' if x else '+' for x in k)
        rows.append([f'$T_{{{"".join(map(str,k))}}}$', str(tuple(2*x for x in k)), s, str(fr30[k]), desc[k]])
tab = ax2.table(cellText=rows[1:], colLabels=rows[0], loc='center', cellLoc='center', colWidths=[0.11,0.13,0.14,0.16,0.3])
tab.auto_set_font_size(False); tab.set_fontsize(8); tab.scale(1.15, 1.5)
for idx in range(8):
    tab[(idx+1, 0)].set_facecolor(cols[idx])
ax2.set_title('Certified frame assignment (frames30)', fontsize=10)
plt.tight_layout(); plt.savefig('fig/supertile3.png', dpi=200); plt.close()

# ---- Figure: certificate funnel
fig, ax = plt.subplots(figsize=(10, 3.2)); ax.axis('off')
steps = [('2388', 'lattice-registered\ntouching poses'), ('44', r'admissible poses $P^*$' + '\n(30 occur hierarchically)'),
         ('33', 'complete enclosures\n(1-shells) of a tile'), ('15', 'extendable to a\nsecond shell'), ('1', 'compatible level-1\nsupertile, always\nfully present')]
x = 0.02
for i, (num, txt) in enumerate(steps):
    ax.add_patch(plt.Rectangle((x, 0.25), 0.16, 0.55, fc=plt.cm.Blues(0.25 + 0.15 * i), ec='k'))
    ax.text(x + 0.08, 0.66, num, ha='center', va='center', fontsize=20, weight='bold')
    ax.text(x + 0.08, 0.40, txt, ha='center', va='center', fontsize=8.5)
    if i < 4: ax.annotate('', xy=(x + 0.205, 0.52), xytext=(x + 0.165, 0.52), arrowprops=dict(arrowstyle='->', lw=1.5))
    x += 0.205
labels = ['facet rules $F^*$\n(135 triples), tight', 'enumerate covers\nof 24 panels', 'exhaustive 2-shell\nsearch (18 dead)', 'uniqueness (U) +\npresence (Pr)']
x = 0.02
for i, l in enumerate(labels):
    ax.text(x + 0.185, 0.12, l, ha='center', va='center', fontsize=7.5, style='italic'); x += 0.205
ax.set_xlim(0, 1.05); ax.set_ylim(0, 1)
plt.savefig('fig/funnel.png', dpi=200, bbox_inches='tight'); plt.close()

# ---- Figure: closure behaviour
fig, ax = plt.subplots(figsize=(6.5, 3.4))
series = {'N=3 certified assignment': [30, 44, 44], 'N=3 "canonical" assignment': [26, 62, 398],
          'N=4 canonical assignment': [46, 158, 3354]}
for (k, v), m in zip(series.items(), ['o', 's', '^']):
    ax.plot(range(len(v)), v, marker=m, label=k)
ax.set_yscale('log'); ax.set_xlabel('coarsening iteration'); ax.set_ylabel(r'$|P_i|$ (admissible poses)')
ax.set_xticks([0, 1, 2]); ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax.annotate('closed: $P_2 = P_1$', xy=(2, 44), xytext=(1.3, 15), fontsize=8, arrowprops=dict(arrowstyle='->'))
ax.annotate('offset (odd) supertile\nposes appear next', xy=(2, 3354), xytext=(0.9, 1200), fontsize=8, arrowprops=dict(arrowstyle='->'))
plt.tight_layout(); plt.savefig('fig/closure.png', dpi=200); plt.close()
print('ok')
