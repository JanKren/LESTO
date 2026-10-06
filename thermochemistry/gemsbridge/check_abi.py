#!/usr/bin/env python3
"""Check that the declared C ABI is the complete exported interface."""
from pathlib import Path
import re
import subprocess
import sys
library=Path(sys.argv[1]);header=Path(sys.argv[2])
expected=set(re.findall(r'\b(gemsb_\w+)\s*\(',header.read_text()))
lines=subprocess.check_output(['nm','-D','--defined-only',str(library)],text=True).splitlines()
actual={line.split()[-1] for line in lines}
assert actual==expected,{'missing':expected-actual,'extra':actual-expected}
needed=re.findall(r'Shared library: \[([^]]+)\]',subprocess.check_output(['readelf','-d',str(library)],text=True))
assert all(n in ['libc.so.6','libm.so.6','ld-linux-x86-64.so.2'] for n in needed),needed
print(f'exported: {len(actual)} gemsb_*, 0 other; dependencies: '+', '.join(needed))
