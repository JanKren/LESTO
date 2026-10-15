import re, numpy as np
def _list(txt, start):
    m = re.compile(r'(nonuniform List<scalar>|uniform)\s*').search(txt, start)
    s = m.end()
    if m.group(1) == 'uniform':
        v = float(re.match(r'([^;]+);', txt[s:]).group(1)); return v
    n = int(re.match(r'\s*(\d+)', txt[s:]).group(1))
    p = txt.index('(', s); q = txt.index(')', p)
    a = np.array(txt[p+1:q].split(), float); assert len(a) == n; return a
def read(path, nx=100, ny=24):
    txt = open(path).read()
    i = txt.index('internalField')
    internal = _list(txt, i)
    if np.isscalar(internal): internal = np.full(nx*ny, internal)
    out = {'internal': internal.reshape(ny, nx)}   # [j, i]
    j = txt.index('WALL', txt.index('boundaryField'))
    k = txt.index('value', j) if 'value' in txt[j:txt.index('}', j)] else None
    out['WALL'] = _list(txt, k) if k else None
    return out
