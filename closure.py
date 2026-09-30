import warnings; warnings.filterwarnings('ignore')
import pickle, time, sys
from cert3 import *
def closure(ch, verbose=True, maxit=8):
    P,k=language(ch,verbose=False)
    P=set(P); hist=[len(P)]
    for it in range(maxit):
        R=Rules(ch,ch.local_cubes,P)
        RS,cand=supertile_rules(ch,R)
        odd=[q for q in cand if any(x%2 for x in q[1])]
        new=set((g,tuple(x//2 for x in t)) for g,t in cand if not any(x%2 for x in t))
        if verbose: print(f"  it{it}: |P|={len(P)} pairwise supertile poses={len(cand)} offset(odd)={len(odd)} new={len(new-P)}",flush=True)
        if odd: return None, hist, 'offset'
        if new<=P: return frozenset(P), hist, 'closed'
        P|=new; hist.append(len(P))
    return None, hist, 'diverged'
if __name__=="__main__":
    ex=pickle.load(open('ex3.pkl','rb'))
    for key in sorted(ex):
        if key[1] is False: continue
        frames=ex[key]; ch=Chair(3,frames)
        print("class",key,frames)
        Pstar,hist,status=closure(ch)
        print("  ->",status,hist, None if Pstar is None else len(Pstar), flush=True)
