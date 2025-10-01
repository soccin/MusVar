#!/usr/bin/env python3

import os
import sys

sl=[]
for s in sys.argv[1:]:
    sl.append(s[::-1])

common_suffix = os.path.commonprefix(sl)[::-1]

# Return suffix starting from first period
if common_suffix and '.' in common_suffix:
    dot_index = common_suffix.index('.')
    print(common_suffix[dot_index:])
else:
    print("")
