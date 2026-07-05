import sys,time
sys.path.append('/mnt/data')
from solve_bits import N,MASK,LEAK
from sat_mul import solve_candidate
lo=int(sys.argv[1]); hi=int(sys.argv[2]); timeout=int(sys.argv[3]) if len(sys.argv)>3 else 120
p=solve_candidate(N,MASK,LEAK,lo,hi,timeout)
print('RESULT',p if p in [None,False] else hex(p), 'div', None if not isinstance(p,int) else N%p)
