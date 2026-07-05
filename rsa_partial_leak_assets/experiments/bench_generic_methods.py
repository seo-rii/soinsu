import sys,time,sympy as sp, random
sys.path.append('/mnt/data')
from hm_bivar import iiter, poly_pow, poly_shift, poly_monomials_power_linear, mon_div
from fpylll import IntegerMatrix, LLL
bits=512 if len(sys.argv)<2 else int(sys.argv[1])
# use challenge constants maybe if bits=0
if bits==0:
 from solve_bits import N,LEAK
 def W(s,l): return ((1<<l)-1)<<s
 base=LEAK; p0=base & ~W(265,155) & ~W(600,230); X=1<<155; Y=1<<230; s1=265; s2=600
else:
 p=int(sp.randprime(1<<(bits-1),1<<bits)); q=int(sp.randprime(1<<(bits-1),1<<bits)); N=p*q; s1=bits//4; l1=bits*155//1024; s2=bits*600//1024; l2=bits*230//1024; X=1<<l1; Y=1<<l2; mask=((1<<bits)-1)^(((X-1)<<s1)|((Y-1)<<s2)); p0=p&mask
m,d=3,4
ainv=pow(1<<s1,-1,N); c0=(p0*ainv)%N; cy=((1<<s2)*ainv)%N
f={(1,0):1,(0,1):cy%N,(0,0):c0%N}; l=(1,0)
Msets=[]; fm_mons=poly_monomials_power_linear(m)
for k in range(m+1):
 Mk=set(); T=poly_monomials_power_linear(m-k)
 for mon in fm_mons:
  div=mon_div(mon,(k,0))
  if div is not None and div in T:
   for ex in range(d):
    for ey in range(d): Mk.add((mon[0]+ex,mon[1]+ey))
 Msets.append(Mk)
Msets.append(set()); shifts=[]
for k in range(m+1):
 fk=poly_pow(f,k)
 for mon in sorted(Msets[k]-Msets[k+1]):
  div=mon_div(mon,(k,0)); shifts.append(poly_shift(fk,div[0],div[1],pow(N,m-k)))
monoms=sorted(set().union(*[set(s.keys()) for s in shifts]),key=lambda ij:(ij[0]+ij[1],ij[0],ij[1])); factors=[X**i*Y**j for i,j in monoms]
Mat=[[0]*len(monoms) for _ in shifts]
for i,sh in enumerate(shifts):
 for j,mon in enumerate(monoms): Mat[i][j]=sh.get(mon,0)*factors[j]
print('dim',len(Mat),len(monoms),'bits',bits,'maxbits',max(abs(v).bit_length() for row in Mat for v in row if v),flush=True)
method=sys.argv[2] if len(sys.argv)>2 else 'heuristic'; delta=float(sys.argv[3]) if len(sys.argv)>3 else 0.99; prec=int(sys.argv[4]) if len(sys.argv)>4 else 256
B=IntegerMatrix.from_matrix(Mat); st=time.time(); kwargs={'method':method}
if prec: kwargs.update(float_type='mpfr',precision=prec)
else: kwargs.update(float_type='double')
try:
 LLL.reduction(B,delta=delta,eta=0.51,**kwargs); print('done',time.time()-st, 'b0bits',max(abs(int(B[0,j])).bit_length() for j in range(B.ncols)), flush=True)
except Exception as e:
 print('ERR',repr(e),time.time()-st, flush=True)
