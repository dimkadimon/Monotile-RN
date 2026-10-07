import sys, time
sys.path.insert(0,'/home/user/research'); sys.path.insert(0,'/home/user/research/next')
import warnings; warnings.filterwarnings('ignore')
from rigidity import Hier
from polytile import normalize_cells, dissection_classes
from chairN import make_group
from cert3 import log
cells=normalize_cells(frozenset([(0,0,0,0),(0,0,0,1),(0,0,0,2),(0,0,1,0),(0,1,0,2),(1,0,1,0)])); G=make_group(4)
classes,ntot,_=dissection_classes(cells,G)
H=Hier(cells,classes[1]); P=H.adjacency_language(); log(f"tw6 class 1: |P|={len(P)} stabilised at level {H.Plev}")
