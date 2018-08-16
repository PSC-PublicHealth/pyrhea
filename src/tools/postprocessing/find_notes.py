#! /usr/bin/env python

import os
import sys
dirnm = sys.argv[1]
cutoff = int(sys.argv[2])
nmL = os.listdir(dirnm)[:]
nmL.sort()
for fn in nmL:
    fullFn = os.path.join(dirnm, fn)
    if fn.endswith('_notes.pkl') and os.path.getsize(fullFn) > cutoff:
        print fn
