exec(open('curve_head.py').read())
import numpy as np
for name in ['e053_linear_May14','e053_linear_PoN_Apr9']:
    calls,items=evalvcf(VCFS[name])
    items.sort(key=lambda x:-x[0])
    ps=np.array([p for p,_,_ in items]); 
    print(name)
    for q in [0.001,0.01,0.05,0.1,0.25,0.5,0.75,0.9,1.0]:
        k=max(1,int(round(q*len(items)))); th=items[k-1][0]
        sub=[x for x in items if x[0]>=th]
        tp=len({x[1] for x in sub if x[2]=='TP'}); fp=sum(1 for x in sub if x[2]!='TP')
        P=tp/(tp+fp); R=tp/697
        print(f'  top {q:.3f} thr>={th:.4f}: TP={tp} FP={fp} P={P:.3f} R={R:.3f}')
