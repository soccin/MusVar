#!/usr/bin/env python3

import os
import sys

sl=[]
for s in sys.argv[1:]:
    sl.append(s[::-1])

print(os.path.commonprefix(sl)[::-1])
