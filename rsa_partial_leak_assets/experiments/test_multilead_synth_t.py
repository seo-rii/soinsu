import sys,time,random
sys.path.append('/mnt/data')
from hm_bivar import linear_bivar_multilead, eval_poly
bits=1024;p=random.getrandbits(bits)|(1<<(bits-1))|1;q=random.getrandbits(bits)|(1<<(bits-1))|1;N=p*q
X=1<<155;Y=1<<230;s1=265;s2=600;x0=(p>>s1)&(X-1);y0=(p>>s2)&(Y-1);mask=((1<<bits)-1)^(((X-1)<<s1)|((Y-1)<<s2));p0=p&mask
for t in [3,4,5,6]:
 print('m8 t',t,flush=True); st=time.time(); hs,_=linear_bivar_multilead(N,p0,1<<s1,1<<s2,X,Y,m=8,t=t,lll_rows=999,verbose=True); print('time',time.time()-st)
 hs=[h for h in hs if h]; zeros=[];minb=99999;pt=p**t
 for i,h in enumerate(hs):
  v=eval_poly(h,x0,y0)
  if v==0: zeros.append(i)
  else: minb=min(minb,abs(v).bit_length())
 print('zeros',len(zeros),zeros[:10],'minb',minb,'ptbits',pt.bit_length())
