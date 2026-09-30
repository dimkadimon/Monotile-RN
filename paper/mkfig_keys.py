import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
def key(x0, y0, w=1.0, eps=0.18, kind='key', flip=False, color='#f5b942'):
    # asymmetric relief profile on a horizontal panel of width w at height y0; key protrudes up, lock is a depression down
    xs = np.array([0, 0.25, 0.30, 0.55, 0.62, 0.85, 1.0]) * w
    ys = np.array([0, 0, 1.0, 1.0, 0.4, 0.4, 0]) * eps
    if flip: xs = w - xs[::-1]; ys = ys[::-1]
    return xs + x0, ys + y0
fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
for ax, title, ok, flip in zip(axes, ['(a) panel $a$ with key, panel $b$ with lock', '(b) flush fit: $(a,b,g)\\in F^*$', '(c) wrong relative orientation: rejected'], [None, True, False], [False, False, True]):
    ax.set_aspect('equal'); ax.axis('off'); ax.set_title(title, fontsize=9)
    gap = 0.45 if ok is None else 0.0
    # lower tile A with key on top panel
    kx, ky = key(0, 1.0)
    A = Polygon(np.vstack([[0, 0], [1, 0], [1, 1.0], *zip(kx[::-1], ky[::-1])]), closed=True, fc='#a6cee3', ec='k', lw=1)
    ax.add_patch(A)
    # upper tile B with lock in bottom panel
    lx, ly = key(0, 1.0 + gap, flip=flip)
    B = Polygon(np.vstack([[0, 2.0 + gap], [0, 1.0 + gap], *zip(lx, ly), [1, 1.0 + gap], [1, 2.0 + gap]]), closed=True, fc='#b2df8a', ec='k', lw=1)
    ax.add_patch(B)
    ax.text(0.5, 0.45, 'tile $A$, panel $a$', ha='center', fontsize=8)
    ax.text(0.5, 1.55 + gap, 'tile $B$, panel $b$', ha='center', fontsize=8)
    if ok is False:
        ax.text(0.5, 1.22, 'overlap', ha='center', fontsize=8, color='red', weight='bold')
    if ok is None:
        ax.annotate('', xy=(1.12, 1.0), xytext=(1.12, 1.18), arrowprops=dict(arrowstyle='<->', lw=0.8)); ax.text(1.16, 1.06, r'$\varepsilon$', fontsize=8)
    ax.set_xlim(-0.1, 1.35); ax.set_ylim(-0.1, 2.6)
plt.tight_layout(); plt.savefig('fig/keys.png', dpi=200); print('ok')
