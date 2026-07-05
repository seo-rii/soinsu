#!/usr/bin/env python3
import sys
from solve7_main import testp
if len(sys.argv) != 2:
    print('usage: python3 decrypt_with_factor.py <p-hex-or-dec>')
    raise SystemExit(2)
s=sys.argv[1].strip()
p=int(s,16) if all(c in '0123456789abcdefABCDEFxX' for c in s) else int(s)
raise SystemExit(0 if testp(p) else 1)
