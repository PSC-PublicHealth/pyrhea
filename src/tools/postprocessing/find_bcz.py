#! /usr/bin/env python

import os
import sys
dirnm = sys.argv[1]
nmL = os.listdir(dirnm)[:]
nmL.sort()
for fn in nmL:
    fullFn = os.path.join(dirnm, fn)
    if fn.endswith('.bcz_0'):
        print fn
