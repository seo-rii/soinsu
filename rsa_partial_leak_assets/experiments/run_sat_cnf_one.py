import sys
sys.path.append('/mnt/data')
from solve_bits import N,MASK,LEAK
from sat_mul_cnf import solve
lo=int(sys.argv[1]); hi=int(sys.argv[2]); solver=sys.argv[3] if len(sys.argv)>3 else 'cadical153'; budget=int(sys.argv[4]) if len(sys.argv)>4 else 100000
p=solve(N,MASK,LEAK,lo,hi,solver,budget)
print('RESULT',p if not isinstance(p,int) else (hex(p),N%p))
