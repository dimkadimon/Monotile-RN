import sys; sys.path.insert(0, '/home/user/research')
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from chairN import *
fr30 = {'0': (0, 1, 2), (0, 0, 0): (0, 1, 2), (0, 0, 1): (0, 2, 1), (0, 1, 0): (2, 1, 0), (0, 1, 1): (1, 2, 0),
        (1, 0, 0): (1, 0, 2), (1, 0, 1): (0, 1, 2), (1, 1, 0): (2, 0, 1)}
ch = Chair(3, fr30)
# ---- data for the supertile figures
grid = np.zeros((4, 4, 4), dtype=int) - 1
keys = ['0'] + [v for v in ch.U]
cols = plt.cm.Set3(np.linspace(0, 1, 12))[[0,4,3,5,6,7,9,11]]
for idx, (h, tau) in enumerate(ch.children):
    for c in ch.cubes((h, tau)):
        u = tuple((x - 1) // 2 for x in c); grid[u] = idx
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
def cube_faces(u):
    x,y,z=u
    P=lambda a,b,c:(x+a,y+b,z+c)
    return {(0,-1):[P(0,0,0),P(0,1,0),P(0,1,1),P(0,0,1)], (0,1):[P(1,0,0),P(1,1,0),P(1,1,1),P(1,0,1)],
            (1,-1):[P(0,0,0),P(1,0,0),P(1,0,1),P(0,0,1)], (1,1):[P(0,1,0),P(1,1,0),P(1,1,1),P(0,1,1)],
            (2,-1):[P(0,0,0),P(1,0,0),P(1,1,0),P(0,1,0)], (2,1):[P(0,0,1),P(1,0,1),P(1,1,1),P(0,1,1)]}
rows = [['child', 'position', 'sign flips $D_s$', r'permutation $\pi_v$', 'rotation $h_v = D_s\pi_v$']]
desc = {'0': 'identity (translation by (1,1,1))', (0,0,0): 'identity', (1,0,0): '180° about $e_2+e_3$', (0,1,0): '180° about $e_1+e_3$',
        (0,0,1): '180° about $e_1+e_2$', (1,0,1): '180° about $e_2$', (0,1,1): '120° about $(1,1,-1)$', (1,1,0): '120° about $(1,-1,-1)$'}
for idx, k in enumerate(keys):
    if k == '0': rows.append(['$T_0$', '(1,1,1)', 'none', '(0,1,2)', desc[k]])
    else:
        s = ''.join('-' if x else '+' for x in k)
        rows.append([f'$T_{{{"".join(map(str,k))}}}$', str(tuple(2*x for x in k)), s, str(fr30[k]), desc[k]])

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


def tile_faces(g, val, shape):
    faces=[]
    for u in zip(*np.where(g==val)):
        for (ax_,sg),poly in cube_faces(u).items():
            nb=list(u); nb[ax_]+=sg
            if not (0<=nb[ax_]<shape) or g[tuple(nb)]!=val: faces.append(poly)
    return faces
def style(ax, L, ticks):
    ax.set_xlim(0,L); ax.set_ylim(0,L); ax.set_zlim(0,L); ax.set_box_aspect((1,1,1)); ax.grid(False)
    for a in (ax.xaxis, ax.yaxis, ax.zaxis): a.set_pane_color((1,1,1,0)); a.set_ticks(ticks)
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')

# ---- Figure: a single chair C_3 and the hat
fig = plt.figure(figsize=(9, 4))
ax = fig.add_subplot(1, 2, 1, projection='3d')
g = np.zeros((2,2,2), dtype=int); g[1,1,1] = -1
allf=[poly for u in zip(*np.where(g==0)) for (ax_,sg),poly in cube_faces(u).items()
      if not (0<=u[ax_]+sg<2) or g[tuple(np.array(u)+np.eye(3,dtype=int)[ax_]*sg)]!=0]
ax.add_collection3d(Poly3DCollection(allf, facecolors=[cols[0]]*len(allf), edgecolors='k', linewidths=0.5))
style(ax, 2, [0,1,2]); ax.view_init(elev=25, azim=-50); ax.set_axis_off()
ax.set_title('(a) the chair $C_3$ (Chair44): 7 unit cubes', fontsize=10)
ax2 = fig.add_subplot(1, 2, 2)
r = np.sqrt(3)
hat = [(0,0),(0,-1),(r/2,-1.5),(r,-2),(1.5*r,-1.5),(r,0),(1.5*r,1.5),(r,2),(r,3),(0,3),(-r/2,1.5),(-r,2),(-1.5*r,1.5),(-r,0)]
hx=[p[0] for p in hat]+[hat[0][0]]; hy=[p[1] for p in hat]+[hat[0][1]]
ax2.fill(hx, hy, color=cols[4], ec='k', lw=1.2)
# kite grid lines (hexagon tiling overlaid with its dual triangles), clipped to the hat
import matplotlib.patches as mpatches
from matplotlib.path import Path
clip = mpatches.PathPatch(Path(list(zip(hx,hy))), transform=ax2.transData, fc='none', ec='none'); ax2.add_patch(clip)
for cx in np.arange(-6, 7, 1.0):
    for cy in np.arange(-6, 7, 1.0):
        pass
R = 2.0  # hexagon circumradius so that kite short edge = 1 (centre to edge midpoint = R*sqrt3/2 = sqrt3; vertex to midpoint = 1)
centers=[(i*np.sqrt(3)*R+(j%2)*np.sqrt(3)*R/2 + np.sqrt(3), j*3*R/2) for i in range(-4,5) for j in range(-4,5)]
for (cx,cy) in centers:
    V=[(cx+R*np.cos(np.pi/3*k+np.pi/6), cy+R*np.sin(np.pi/3*k+np.pi/6)) for k in range(6)]
    for k in range(6):
        a=V[k]; b=V[(k+1)%6]; mid=((a[0]+b[0])/2,(a[1]+b[1])/2)
        for (p,q) in [(a,b),((cx,cy),mid)]:
            l,=ax2.plot([p[0],q[0]],[p[1],q[1]],color=(0,0,0,0.35),lw=0.6); l.set_clip_path(clip)
ax2.set_aspect('equal'); ax2.set_xlim(-3.2,3.2); ax2.set_ylim(-2.6,3.6); ax2.axis('off')
ax2.set_title('(b) the hat (Smith et al. 2024): 8 kites', fontsize=10)
plt.tight_layout(); plt.savefig('fig/chair3.png', dpi=200); plt.close()

# ---- Figure: exploded view, assembled supertile, and frame table (combined)
fig = plt.figure(figsize=(11, 8.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0])
for k, spread in enumerate([3.4, 0.0]):
    ax = fig.add_subplot(gs[0, k], projection='3d')
    faces_all=[]; fc_all=[]
    for idx in range(8):
        faces = tile_faces(grid, idx, 4)
        cen = np.mean([np.mean(f, axis=0) for f in faces], axis=0)
        if idx == 0: off = np.array([-0.64, 0.77, 0.0]) * spread * 0.55   # lift the central translate so T000 becomes visible
        else: off = (cen - 2.0) / np.linalg.norm(cen - 2.0) * spread
        faces_all += [[tuple(np.array(p) + off) for p in f] for f in faces]; fc_all += [cols[idx]] * len(faces)
    ax.add_collection3d(Poly3DCollection(faces_all, facecolors=fc_all, edgecolors='k', linewidths=0.7))
    style(ax, 4, [0,1,2,3,4]); lim = (-spread/2, 4 + spread/2)
    ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_zlim(-spread/2, 4 + spread)
    ax.view_init(elev=25, azim=40)
    if spread: ax.set_axis_off()
    ax.set_title(['(a) the eight children pulled apart', '(b) assembled: the level-1 supertile $2C_3$'][k], fontsize=10)
ax2 = fig.add_subplot(gs[1, :]); ax2.axis('off')
tab = ax2.table(cellText=rows[1:], colLabels=rows[0], loc='center', cellLoc='center', colWidths=[0.08,0.1,0.11,0.13,0.26])
tab.auto_set_font_size(False); tab.set_fontsize(9); tab.scale(1.15, 1.6)
for idx in range(8):
    tab[(idx+1, 0)].set_facecolor(cols[idx])
ax2.set_title('(c) certified frame assignment (frames30)', fontsize=10)
plt.tight_layout(); plt.savefig('fig/supertile3.png', dpi=200); plt.close()
