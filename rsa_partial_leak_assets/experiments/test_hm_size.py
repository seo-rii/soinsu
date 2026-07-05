import sys,time,random,sympy as sp
sys.path.append('/mnt/data')
from hm_bivar import hm_polys_bivar, eval_poly
# Generate 1024-bit primes (may take)
def rand_odd(bits): return random.getrandbits(bits) | (1<<(bits-1)) | 1
p=sp.randprime(1<<1023,1<<1024)
q=sp.randprime(1<<1023,1<<1024)
N=int(p*q); p=int(p)
# unknown x bits 155 at 265 and y 230 at 600, with arbitrary p0 known elsewhere
X=1<<155; Y=1<<230; s1=265; s2=600
x0=(p>>s1)&(X-1); y0=(p>>s2)&(Y-1)
mask=((1<<1024)-1) ^ (((X-1)<<s1)|((Y-1)<<s2))
p0=p&mask
ainv=pow(1<<s1,-1,N); c0=(p0*ainv)%N; cy=((1<<s2)*ainv)%N
print('generated',N.bit_length(),x0.bit_length(),y0.bit_length())
for m,t in [(6,3),(7,3),(8,3),(8,4),(9,4)]:
 st=time.time(); hs,_=hm_polys_bivar(N,c0,cy,X,Y,m=m,t=t,lll_rows=6,verbose=True); print('time',time.time()-st,'p^t bits',(p**t).bit_length())
 for i,h in enumerate(hs[:6]):
  v=eval_poly(h,x0,y0)
  print(i, v==0, abs(v).bit_length(), v%(p**t)==0, abs(v)<p**t)
