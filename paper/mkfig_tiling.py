import sys; sys.path.insert(0, '/home/user/research')
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from chairN import *
fr30 = {'0': (0, 1, 2), (0, 0, 0): (0, 1, 2), (0, 0, 1): (0, 2, 1), (0, 1, 0): (2, 1, 0), (0, 1, 1): (1, 2, 0),
        (1, 0, 0): (1, 0, 2), (1, 0, 1): (0, 1, 2), (1, 1, 0): (2, 0, 1)}
ch = Chair(3, fr30)
rng = np.random.default_rng(3)
fig = plt.figure(figsize=(12, 4.6))
# (a) level-2 supertile 4C_3, 64 tiles
ax = fig.add_subplot(1, 3, 1, projection='3d')
tiles = ch.patch(2); n = 8
grid = -np.ones((n, n, n), int)
for i, T in enumerate(tiles):
    for c in ch.cubes(T):
        grid[tuple((x - 1) // 2 for x in c)] = i
cols = plt.cm.tab20(rng.random(len(tiles)))
fc = np.empty(grid.shape + (4,)); fc[grid >= 0] = cols[grid[grid >= 0]]
ax.voxels(grid >= 0, facecolors=fc, edgecolor='k', linewidth=0.25, shade=False)
ax.view_init(elev=25, azim=40); ax.set_title('(a) level-2 supertile $4C_3$: 64 tiles', fontsize=9)
ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
# (b) same with front octant cut away to show the interior
ax = fig.add_subplot(1, 3, 2, projection='3d')
mask = grid >= 0
cut = np.zeros_like(mask); cut[4:, 4:, 4:] = True   # remove cubes with x,y,z >= 4 ... choose a corner block
cut2 = np.zeros_like(mask); cut2[:4, 4:, 4:] = True
m2 = mask & ~cut2
ax.voxels(m2, facecolors=fc, edgecolor='k', linewidth=0.25, shade=False)
ax.view_init(elev=25, azim=40); ax.set_title('(b) same, block $x<4,\\ y,z\\geq 4$ removed', fontsize=9)
ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
# (c) cross-section of level-3 patch (512 tiles), layer z in [5,6]
ax = fig.add_subplot(1, 3, 3)
tiles3 = ch.patch(3); n = 16
grid3 = -np.ones((n, n, n), int)
for i, T in enumerate(tiles3):
    for c in ch.cubes(T):
        grid3[tuple((x - 1) // 2 for x in c)] = i
layer = 5
cols3 = plt.cm.tab20(rng.random(len(tiles3)))
for x in range(n):
    for y in range(n):
        i = grid3[x, y, layer]
        if i < 0: continue
        ax.add_patch(Rectangle((x, y), 1, 1, fc=cols3[i], ec='none'))
# tile outlines: draw edges between cells of different tiles
for x in range(n):
    for y in range(n):
        i = grid3[x, y, layer]
        if i < 0: continue
        if x == n-1 or grid3[x+1, y, layer] != i: ax.plot([x+1, x+1], [y, y+1], 'k', lw=1.2)
        if x == 0 or grid3[x-1, y, layer] != i: ax.plot([x, x], [y, y+1], 'k', lw=1.2)
        if y == n-1 or grid3[x, y+1, layer] != i: ax.plot([x, x+1], [y+1, y+1], 'k', lw=1.2)
        if y == 0 or grid3[x, y-1, layer] != i: ax.plot([x, x+1], [y, y], 'k', lw=1.2)
ax.set_xlim(0, n); ax.set_ylim(0, n); ax.set_aspect('equal'); ax.set_xlabel('x'); ax.set_ylabel('y')
ax.set_title(f'(c) level-3 supertile $8C_3$ (512 tiles): slice $z\\in[{layer},{layer+1}]$', fontsize=9)
plt.tight_layout(); plt.savefig('fig/tiling3.png', dpi=200); print('ok')
