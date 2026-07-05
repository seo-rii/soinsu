#!/usr/bin/env python3
from solve7_main import N, MASK, LEAK, intervals_from_mask
print('N bits =', N.bit_length())
print('known p bits =', MASK.bit_count())
print('unknown p bits =', 1024 - MASK.bit_count())
for s,e,l in intervals_from_mask(MASK):
    print(f'bits {s}..{e} len={l}')
assert (LEAK & ~MASK) == 0
