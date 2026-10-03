#!/usr/bin/env python3
"""
make_gems3k_input.py -- write a GEMS3K (v4.x) input set without GEM-Selektor.

Chemical system: one ideal-gas mixture (phase class 'g', DC class 'G') plus any
number of pure condensed phases (phase class 's' or 'l', DC class 'O').
No aqueous phase, no charge IC ('Zz'), no sorption, no non-ideal mixing.

Output (key-value text, the format GEM_init(<lst>) reads with "-t"):
    <name>-dat.lst  <name>-dch.dat  <name>-ipm.dat  <name>-dbr-0-0000.dat
    <name>-dbr.lst
and, with --json, the same three structures as JSON documents
(<name>-dch.json ...; needed for GEM_init(dch, ipm, dbr) from strings).

Units written (GEMS3K v4 conventions, see datach_formats.cpp:61-62,147-172):
    T in K, P in Pa, molar masses in kg/mol, G0/H0 in J/mol, S0/Cp0 in J/(mol K),
    V0 in J/Pa (= m3/mol).  G0 is the standard-state (P0 = 1 bar) molar Gibbs
    energy; for gases NO ln(P) term is added here (GEMS3K adds ln(P/bar) itself,
    ipm_chemical.cpp:865-871).  Condensed G0 gets + V0*(P - 1e5 Pa).

Thermodynamic data sources accepted per species (all in one self-consistent
element reference convention -- do not mix conventions between species):
  * "nasa9": name of a record in a NASA-Glenn thermo.inp file (McBride et al.
    9-coefficient polynomials; evaluated as in nasa/cea source/fits.f90).
  * "cp": {"Tref":298.15, "H":J/mol, "S":J/(mol K),
           "intervals":[{"Tmax":K, "cp":{exponent: coeff,...}, "dHtr":J/mol}, ...]}
    Cp(T) = sum_k coeff_k * T**exponent_k  (any real exponents; e.g. HSC form
    A + B*1e-3 T + C*1e5 T^-2 + D*1e-6 T^2 -> {0:A, 1:B*1e-3, -2:C*1e5, 2:D*1e-6}).
    dHtr is a transition enthalpy added at the LOWER bound of that interval.

Everything outside the fitted ranges is extrapolated with the nearest
interval (needed: GEMS3K requires G0 of every DC at every grid T, e.g.
supercooled liquid PbI2 at 300 K and superheated crystal at 1100 K).

This is a research sketch; check every generated number before production use.
"""
import argparse
import json
import math
import os
import sys

R = 8.314462618           # J/(mol K), CODATA 2018
P0 = 1.0e5                # Pa, standard-state pressure of the data (1 bar)

# ----------------------------------------------------------------------------
# NASA-9 (McBride/Zehe/Gordon) reader and evaluator
# ----------------------------------------------------------------------------

def _f(s):
    s = s.strip().replace('D', 'E').replace('d', 'e')
    return float(s) if s else 0.0


def read_nasa9(path, wanted):
    """Return {name: record} for the names in `wanted` from a thermo.inp file."""
    wanted = set(wanted)
    out = {}
    with open(path) as fh:
        lines = fh.read().splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i]
        name = ln[:18].strip()
        if (ln and not ln[0].isspace() and not ln.startswith('!')
                and name in wanted and i + 1 < len(lines)):
            hdr = lines[i + 1]
            nint = int(hdr[0:2])
            formula = {}
            for k in range(5):
                chunk = hdr[10 + 8 * k: 18 + 8 * k]
                sym, cnt = chunk[:2].strip(), _f(chunk[2:8])
                if sym and cnt != 0.0:
                    sym = sym[0] + sym[1:].lower()     # 'PB' -> 'Pb'
                    formula[sym] = formula.get(sym, 0.0) + cnt
            phase_flag = int(hdr[50:52])
            mw = _f(hdr[52:65])                          # g/mol
            hf = _f(hdr[65:80])                          # J/mol at 298.15 K
            ints = []
            for n in range(nint):
                a = lines[i + 2 + 3 * n]
                b = lines[i + 3 + 3 * n]
                c = lines[i + 4 + 3 * n]
                tmin, tmax = _f(a[0:11]), _f(a[11:22])
                co = [_f(b[16 * k:16 * k + 16]) for k in range(5)]
                co += [_f(c[0:16]), _f(c[16:32])]
                b1, b2 = _f(c[48:64]), _f(c[64:80])
                ints.append(dict(Tmin=tmin, Tmax=tmax, a=co, b1=b1, b2=b2))
            out[name] = dict(name=name, formula=formula, condensed=phase_flag,
                             mw_g=mw, Hf298=hf, intervals=ints,
                             source='%s (line %d)' % (os.path.basename(path), i + 1))
            i += 2 + 3 * nint
            continue
        i += 1
    missing = wanted - set(out)
    if missing:
        raise KeyError('not found in %s: %s' % (path, sorted(missing)))
    return out


def _nasa_pick(rec, T):
    # same interval choice as nasa/cea source/thermo.f90 (last interval with T > Tmin)
    idx = 0
    for k, iv in enumerate(rec['intervals']):
        if T > iv['Tmin']:
            idx = k
    return rec['intervals'][idx]


def nasa_props(rec, T):
    """(Cp, H, S, G) in J/mol, J/(mol K) at T (K), P0 = 1 bar."""
    iv = _nasa_pick(rec, T)
    a1, a2, a3, a4, a5, a6, a7 = iv['a']
    b1, b2 = iv['b1'], iv['b2']
    lnT = math.log(T)
    cp = (a1 / T**2 + a2 / T + a3 + a4 * T + a5 * T**2 + a6 * T**3 + a7 * T**4)
    h = (-a1 / T + a2 * lnT + a3 * T + a4 * T**2 / 2 + a5 * T**3 / 3
         + a6 * T**4 / 4 + a7 * T**5 / 5 + b1)                 # H/R  [K]
    s = (-a1 / (2 * T**2) - a2 / T + a3 * lnT + a4 * T + a5 * T**2 / 2
         + a6 * T**3 / 3 + a7 * T**4 / 4 + b2)                 # S/R
    return R * cp, R * h, R * s, R * (h - T * s)

# ----------------------------------------------------------------------------
# Generic piecewise Cp-polynomial model (for tabulated H298, S298, Cp(T))
# ----------------------------------------------------------------------------

def _int_T(e, T1, T2):          # integral of T**e dT
    return math.log(T2 / T1) if abs(e + 1) < 1e-12 else (T2**(e + 1) - T1**(e + 1)) / (e + 1)


def _int_T_over_T(e, T1, T2):   # integral of T**(e-1) dT
    return math.log(T2 / T1) if abs(e) < 1e-12 else (T2**e - T1**e) / e


def cp_props(model, T):
    Tref, H, S = model.get('Tref', 298.15), model['H'], model['S']
    ivs = model['intervals']
    bounds = [Tref] + [iv['Tmax'] for iv in ivs]    # lower bound of iv k = bounds[k]

    def cpk(k, t):
        return sum(float(c) * t**float(e) for e, c in ivs[k]['cp'].items())

    def dH(k, t1, t2):
        return sum(float(c) * _int_T(float(e), t1, t2) for e, c in ivs[k]['cp'].items())

    def dS(k, t1, t2):
        return sum(float(c) * _int_T_over_T(float(e), t1, t2) for e, c in ivs[k]['cp'].items())

    if T <= Tref:                                   # extrapolate first interval down
        k = 0
        return cpk(0, T), H + dH(0, Tref, T), S + dS(0, Tref, T), None
    k, t = 0, Tref
    while True:
        last = (k == len(ivs) - 1)
        t_hi = T if (last or T <= bounds[k + 1]) else bounds[k + 1]
        H += dH(k, t, t_hi)
        S += dS(k, t, t_hi)
        if t_hi == T:
            return cpk(k, T), H, S, None
        t = t_hi
        k += 1
        dHtr = float(ivs[k].get('dHtr', 0.0))
        H += dHtr
        S += dHtr / t


def species_props(sp, T):
    """(Cp, H, S, G) at T and P0 in the species' own G convention (see g_shift)."""
    if 'nasa' in sp:
        return nasa_props(sp['nasa'], T)
    if 'tf' in sp:                      # (thermofun.ThermoEngine, symbol) - record as given
        eng, sym = sp['tf']
        p = eng.thermoPropertiesSubstance(T, P0, sym)
        return (p.heat_capacity_cp.val, p.enthalpy.val, p.entropy.val, p.gibbs_energy.val)
    cp, H, S, _ = cp_props(sp['cp'], T)
    return cp, H, S, H - T * S


TREF = 298.15


def g_shift(sp, system):
    """Constant (T-independent) J/mol to add to G of `sp` to express it in the
    system-wide convention.  'absolute': G = H - T*S with H(298.15) = DfH and
    third-law S (NASA/JANAF/HSC style).  'formation': Benson-Helgeson apparent
    Gibbs energy, G(298.15) = DfG (GEM-Selektor/ThermoFun/ThermoHub style).
    The two differ by TREF * sum_e nu_e * S_e(298.15) of the elements in their
    reference states - a linear function of the formula, so each convention is
    self-consistent, but species of different conventions must NOT be mixed."""
    own = sp.get('convention', 'formation' if 'tf' in sp else 'absolute')
    want = system.get('G_convention', 'absolute')
    if own == want:
        return 0.0
    Se = system['S_elements']
    corr = TREF * sum(n * Se[e] for e, n in sp['formula'].items())
    return corr if (own == 'absolute' and want == 'formation') else -corr

# ----------------------------------------------------------------------------
# Writers
# ----------------------------------------------------------------------------

def _fmt(x):
    return repr(float(x)) if isinstance(x, float) else str(x)


def _kv_array(fh, tag, values, per_line=8, quote=False):
    fh.write('<%s>\n' % tag)
    for k in range(0, len(values), per_line):
        chunk = values[k:k + per_line]
        if quote:
            fh.write(' '.join("'%s'" % v for v in chunk) + '\n')
        else:
            fh.write(' '.join(_fmt(v) for v in chunk) + '\n')


def build(system):
    """Compute all DCH/IPM/DBR content from the system description (dict)."""
    ics = system['ICs']                      # [(name, kg/mol)], order = ICNL order
    icn = [n for n, _ in ics]
    phases = system['phases']                # list of dicts, gas mixture first
    if phases[0]['class'] not in 'gfp' or len(phases[0]['species']) < 2:
        raise ValueError('first phase must be the multi-component gas phase')
    for ph in phases[1:]:
        if len(ph['species']) != 1:
            raise ValueError('this sketch supports only pure condensed phases')
    TK, Pv = system['TKval'], system['Pval']
    dcs = [(ph, sp) for ph in phases for sp in ph['species']]
    nDC, nIC, nPH = len(dcs), len(ics), len(phases)
    nTp, nPp = len(TK), len(Pv)
    A, DCmm, G0, V0, H0, S0, Cp0 = [], [], [], [], [], [], []
    for ph, sp in dcs:
        form = sp['formula']
        bad = set(form) - set(icn)
        if bad:
            raise ValueError('%s uses ICs not in ICNL: %s' % (sp['name'], bad))
        A += [float(form.get(n, 0.0)) for n in icn]
        DCmm.append(sum(form.get(n, 0.0) * m for n, m in ics))
        gas = ph['class'] in 'gfp'
        for p in Pv:
            for t in TK:
                cp, h, s, g = species_props(sp, t)
                g += g_shift(sp, system)
                if gas:
                    v = R * t / p                       # ideal gas, m3/mol = J/Pa
                else:
                    v = sp['V0']                        # m3/mol, P,T independent here
                    g += v * (p - P0)                   # condensed: V dP correction
                    h += v * (p - P0)
                G0.append(g); V0.append(v); H0.append(h); S0.append(s); Cp0.append(cp)
    return dict(icn=icn, icmm=[m for _, m in ics], dcs=dcs, A=A, DCmm=DCmm,
                G0=G0, V0=V0, H0=H0, S0=S0, Cp0=Cp0, nDC=nDC, nIC=nIC,
                nPH=nPH, nTp=nTp, nPp=nPp, TK=TK, Pv=Pv, phases=phases)


def write_kv(system, d, outdir, name, with_hs=True):
    os.makedirs(outdir, exist_ok=True)
    ph = d['phases']
    nPS = 1
    nDCs = len(ph[0]['species'])
    dcnl = [sp['name'] for _, sp in d['dcs']]
    ccdc = [('G' if p['class'] in 'gfp' else 'O') for p, _ in d['dcs']]
    grid = d['nTp'] * d['nPp']
    # ---------------- DCH ----------------
    with open(os.path.join(outdir, name + '-dch.dat'), 'w') as fh:
        fh.write('# GEMS3K DCH file written by make_gems3k_input.py (no GEM-Selektor)\n')
        fh.write('# units: T K, P Pa, masses kg/mol, G0/H0 J/mol, S0/Cp0 J/(mol K), V0 J/Pa\n')
        for tag, val in [('nIC', d['nIC']), ('nDC', d['nDC']), ('nPH', d['nPH']),
                         ('nPS', nPS), ('nDCs', nDCs), ('nICb', d['nIC']),
                         ('nDCb', d['nDC']), ('nPHb', d['nPH']), ('nPSb', nPS),
                         ('nTp', d['nTp']), ('nPp', d['nPp']), ('iGrd', 0),
                         ('fAalp', 0), ('mLook', 0)]:
            fh.write('<%s> %d\n' % (tag, val))
        fh.write('\n<END_DIM>\n\n')
        # xic/xdc/xph omitted: all ICs/DCs/phases are kept in DBR (trivial lists)
        _kv_array(fh, 'ICNL', d['icn'], quote=True)
        _kv_array(fh, 'ccIC', ['e'] * d['nIC'], quote=True)
        _kv_array(fh, 'ICmm', d['icmm'])
        fh.write('# DCs: gas species first (phase order), then pure condensed phases\n')
        _kv_array(fh, 'DCNL', dcnl, per_line=6, quote=True)
        _kv_array(fh, 'ccDC', ccdc, quote=True)
        _kv_array(fh, 'DCmm', d['DCmm'])
        _kv_array(fh, 'PHNL', [p['name'] for p in ph], per_line=6, quote=True)
        _kv_array(fh, 'ccPH', [p['class'] for p in ph], quote=True)
        _kv_array(fh, 'nDCinPH', [len(p['species']) for p in ph])
        fh.write('# A: stoichiometry, one row per DC, columns = ICNL order %s\n' % d['icn'])
        _kv_array(fh, 'A', d['A'], per_line=d['nIC'])
        fh.write('<Ttol> %s\n' % _fmt(float(system.get('Ttol', 0.5))))
        _kv_array(fh, 'TKval', d['TK'])
        fh.write('<Ptol> %s\n' % _fmt(float(system.get('Ptol', 100.0))))
        _kv_array(fh, 'Pval', d['Pv'])
        fh.write('# Arrays below: [nDC][nPp][nTp], T index fastest\n')
        for tag in (['V0', 'G0'] + (['H0', 'S0', 'Cp0'] if with_hs else [])):
            fh.write('<%s>\n' % tag)
            arr = d[tag]
            for j, nm in enumerate(dcnl):
                fh.write('# %s\n' % nm)
                blk = arr[j * grid:(j + 1) * grid]
                for p in range(d['nPp']):
                    row = blk[p * d['nTp']:(p + 1) * d['nTp']]
                    fh.write(' '.join('%.10g' % v for v in row) + '\n')
    # ---------------- IPM ----------------
    b = [float(system['bIC'].get(n, 0.0)) for n in d['icn']]
    with open(os.path.join(outdir, name + '-ipm.dat'), 'w') as fh:
        fh.write('# GEMS3K IPM file written by make_gems3k_input.py\n')
        fh.write('<ID_key> "%s"\n' % name[:60])
        fh.write('<pa_PE> 0\n<PV> 0\n<PSOL> 0\n<PAalp> \'-\'\n<PSigm> \'-\'\n')
        fh.write('<Lads> 0\n<FIa> 0\n<FIat> 0\n\n<END_DIM>\n\n')
        ctl = system.get('ipm_controls', {})
        for tag, val in ctl.items():
            fh.write('<%s> %s\n' % (tag, val))
        fh.write('# one ideal multicomponent phase (gas): TSolMod code I = ideal\n')
        _kv_array(fh, 'sMod', ['INNINNNN'], quote=True)
        _kv_array(fh, 'LsMod', [0, 0, 0])
        _kv_array(fh, 'LsMdc', [0, 0, 0])
        _kv_array(fh, 'B', b)
    # ---------------- DBR ----------------
    with open(os.path.join(outdir, name + '-dbr-0-0000.dat'), 'w') as fh:
        fh.write('# GEMS3K DBR file written by make_gems3k_input.py\n')
        fh.write('<NodeHandle> 0\n<NodeTypeHY> 0\n<NodeTypeMT> 0\n')
        fh.write('<NodeStatusFMT> -1\n<NodeStatusCH> 1\n<IterDone> 0\n')
        fh.write('<TK> %s\n<P> %s\n' % (_fmt(float(system['TK'])), _fmt(float(system['P']))))
        _kv_array(fh, 'bIC', b)
    with open(os.path.join(outdir, name + '-dat.lst'), 'w') as fh:
        fh.write('-t "%s-dch.dat" "%s-ipm.dat" "%s-dbr-0-0000.dat"\n' % (name, name, name))
    with open(os.path.join(outdir, name + '-dbr.lst'), 'w') as fh:
        fh.write('"%s-dbr-0-0000.dat"\n' % name)


def write_json(system, d, outdir, name, with_hs=True):
    """Same content as JSON documents (for GEM_init(dch_json, ipm_json, dbr_json))."""
    ph = d['phases']
    dch = dict(nIC=d['nIC'], nDC=d['nDC'], nPH=d['nPH'], nPS=1,
               nDCs=len(ph[0]['species']), nICb=d['nIC'], nDCb=d['nDC'],
               nPHb=d['nPH'], nPSb=1, nTp=d['nTp'], nPp=d['nPp'], iGrd=0,
               fAalp=0, mLook=0,
               ICNL=d['icn'], ccIC=['e'] * d['nIC'], ICmm=d['icmm'],
               DCNL=[sp['name'] for _, sp in d['dcs']],
               ccDC=[('G' if p['class'] in 'gfp' else 'O') for p, _ in d['dcs']],
               DCmm=d['DCmm'], PHNL=[p['name'] for p in ph],
               ccPH=[p['class'] for p in ph],
               nDCinPH=[len(p['species']) for p in ph], A=d['A'],
               Ttol=float(system.get('Ttol', 0.5)), TKval=d['TK'],
               Ptol=float(system.get('Ptol', 100.0)), Pval=d['Pv'],
               V0=d['V0'], G0=d['G0'])
    if with_hs:
        dch.update(H0=d['H0'], S0=d['S0'], Cp0=d['Cp0'])
    b = [float(system['bIC'].get(n, 0.0)) for n in d['icn']]
    ipm = {'ID_key': name, 'pa_PE': 0, 'PV': 0, 'PSOL': 0, 'PAalp': '-',
           'PSigm': '-', 'Lads': 0, 'FIa': 0, 'FIat': 0,
           'sMod': ['INNINNNN'], 'LsMod': [0, 0, 0], 'LsMdc': [0, 0, 0], 'B': b}
    for tag, val in system.get('ipm_controls', {}).items():
        ipm[tag] = val
    dbr = dict(NodeHandle=0, NodeTypeHY=0, NodeTypeMT=0, NodeStatusFMT=-1,
               NodeStatusCH=1, IterDone=0, TK=float(system['TK']),
               P=float(system['P']), bIC=b)
    for key, doc in (('dch', dch), ('ipm', ipm), ('dbr-0-0000', dbr)):
        tag = key.split('-')[0]
        with open(os.path.join(outdir, '%s-%s.json' % (name, key)), 'w') as fh:
            json.dump([{'_key': name + '-' + tag, 'set': name, tag: doc}], fh, indent=1)
    with open(os.path.join(outdir, name + '-dat-json.lst'), 'w') as fh:
        fh.write('-j "%s-dch.json" "%s-ipm.json" "%s-dbr-0-0000.json"\n' % (name, name, name))

def write_thermofun_json(system, d, outdir, name, element_S=None):
    """ThermoFun database (<name>-fun.json) with the same species, so that GEMS3K
    (built with USE_THERMOFUN) computes G0(T,P) at run time instead of using the
    DCH look-up grid ("-o"/"-f" modes of the lst file).  NASA-9 intervals are
    mapped 1:1 onto ThermoFun cp_ft_equation coefficients
    [a0 + a1 T + a2 T^-2 + a3 T^-0.5 + a4 T^2 + a5 T^3 + a6 T^4 + a7 T^-3 + a8 T^-1 ...]
    (ThermoFun EmpiricalCpIntegration.cpp).  Reference properties at Tst = 298.15 K
    are taken from the NASA functions, G298 = H298 - Tst*S298 (same convention as
    the DCH grid)."""
    Tst = 298.15
    subs = []
    for ph, sp in d['dcs']:
        gas = ph['class'] in 'gfp'
        if 'tf_record' in sp:            # copy a ThermoFun record (e.g. from ThermoHub) as is
            rec = json.loads(json.dumps(sp['tf_record']))
            rec['symbol'] = sp['name']
            if gas:   # GEMS3K adds ln(P/bar) for gases itself: G0 must stay the 1-bar value.
                # 'mv_constant' on a gas record makes ThermoFun add V*(P-Pst) to G
                # (found in ThermoHub HERACLES gas records) -> replace by 'mv_pvnrt'.
                rec['TPMethods'] = [m for m in rec['TPMethods']
                                    if not any(v.startswith('mv_') for v in m['method'].values())]
                rec['TPMethods'].append({'method': {'40': 'mv_pvnrt'}})
            if sp.get('convention', 'formation') != system.get('G_convention', 'absolute'):
                sh = g_shift(sp, system)
                rec['sm_gibbs_energy']['values'][0] += sh   # H and S are the same in both conventions
            subs.append(rec)
            continue
        rec = sp['nasa']
        cp, h, s, g = nasa_props(rec, Tst)
        g += g_shift(sp, system)
        tpm = []
        ivs = rec['intervals']
        for k, iv in enumerate(ivs):
            a1, a2, a3, a4, a5, a6, a7 = iv['a']
            lo = min(iv['Tmin'], 200.0) if k == 0 else iv['Tmin']
            m = {'method': {'0': 'cp_ft_equation'},
                 'limitsTP': {'range': True, 'lowerT': lo, 'upperT': iv['Tmax']},
                 'm_heat_capacity_ft_coeffs': {'values': [R * a3, R * a4, R * a1, 0.0, R * a5,
                                                          R * a6, R * a7, 0.0, R * a2, 0.0, 0.0]}}
            if k > 0:   # interval boundary without transition (NASA fits are continuous)
                m['m_phase_trans_props'] = {'values': [iv['Tmin'], 0.0, 0.0, 0.0, 0.0]}
            tpm.append(m)
        tpm.append({'method': {'40': 'mv_pvnrt'}} if gas else {'method': {'34': 'mv_constant'}})
        form = ''.join('%s%s' % (e, ('' if abs(n - 1) < 1e-12 else ('%g' % n)))
                       for e, n in sp['formula'].items())
        subs.append({
            'symbol': sp['name'], 'name': sp['name'], 'formula': form, 'formula_charge': 0,
            'aggregate_state': ({'0': 'AS_GAS'} if gas else
                                ({'1': 'AS_LIQUID'} if ph['class'] == 'l' else {'3': 'AS_CRYSTAL'})),
            'class_': ({'1': 'SC_GASFLUID'} if gas else {'0': 'SC_COMPONENT'}),
            'Tst': Tst, 'Pst': P0,
            'mass_per_mole': {'values': [sum(sp['formula'].get(n, 0) * mm for n, mm in
                                              zip(d['icn'], d['icmm'])) * 1e3], 'units': ['g/mol']},
            'sm_gibbs_energy': {'values': [g], 'units': ['J/mol']},
            'sm_enthalpy': {'values': [h], 'units': ['J/mol']},
            'sm_entropy_abs': {'values': [s], 'units': ['J/(mol*K)']},
            'sm_heat_capacity_p': {'values': [cp], 'units': ['J/(mol*K)']},
            'sm_volume': {'values': [0.0 if gas else sp['V0'] * 1e5], 'units': ['J/bar']},
            'datasources': [rec['source']],
            'TPMethods': tpm})
    element_S = element_S or {}
    elems = [{'symbol': e, 'class_': {'0': 'ELEMENT'}, 'atomic_mass': {'values': [m * 1e3]},
              'entropy': {'values': [element_S.get(e, 0.0)]}, 'datasources': ['user']}
             for e, m in zip(d['icn'], d['icmm'])]
    db = {'date': '2026-09-23', 'datasources': ['make_gems3k_input.py'],
          'thermodataset': ['user'], 'elements': elems, 'substances': subs, 'reactions': []}
    with open(os.path.join(outdir, name + '-fun.json'), 'w') as fh:
        json.dump(db, fh, indent=1)


# ----------------------------------------------------------------------------
# Example: He + PbI2 system with NASA-Glenn (Gurvich 1991) data
# ----------------------------------------------------------------------------

def pbi2_system(nasa_path, liquid_class='l', with_pbi4=False, TK=1000.0, P=101325.0,
                x_pbi2=1.0e-6):
    gas_names = ['He', 'I', 'I2', 'Pb', 'PbI', 'PbI2'] + (['PbI4'] if with_pbi4 else [])
    cond_names = [('PbI2(cr)', 's', 7.484e-5), ('PbI2(L)', liquid_class, 7.484e-5),
                  ('Pb(cr)', 's', 1.827e-5), ('Pb(L)', liquid_class, 1.827e-5)]
    # V0 (m3/mol) = M/rho with rho(PbI2) = 6.16 g/cm3 (PubChem CID 24931),
    # rho(Pb) = 11.34 g/cm3 (PubChem CID 5352425); liquids reuse the solid value
    # (placeholder: only affects reported phase volumes, <2 J/mol in G0 at 1 atm).
    recs = read_nasa9(nasa_path, gas_names + [c[0] for c in cond_names])
    mw = {'He': recs['He']['mw_g'] * 1e-3, 'I': recs['I']['mw_g'] * 1e-3,
          'Pb': recs['Pb']['mw_g'] * 1e-3}
    phases = [dict(name='gas_gen', **{'class': 'g'},
                   species=[dict(name=n + '(g)', formula=recs[n]['formula'], nasa=recs[n])
                            for n in gas_names])]
    for n, cls, v in cond_names:
        short = n.replace('(L)', '(l)')
        phases.append({'name': short, 'class': cls,
                       'species': [dict(name=short, formula=recs[n]['formula'],
                                        nasa=recs[n], V0=v)]})
    # bulk: 1 mol He carrier + x_pbi2 mol PbI2 (as elements)
    return dict(ICs=[('He', mw['He']), ('I', mw['I']), ('Pb', mw['Pb'])],
                phases=phases,
                TKval=[270.0 + 10.0 * k for k in range(99)],     # 270 ... 1250 K
                Pval=[0.8e5, 1.0e5, 1.2e5],                       # Pa
                Ttol=0.01, Ptol=1.0, TK=TK, P=P,
                bIC={'He': 1.0, 'I': 2.0 * x_pbi2, 'Pb': x_pbi2},
                # GEMS3K defaults (ms_multi_diff.cpp:34-48) except the last three:
                # pa_DHB 1e-13 -> 1e-10 and pa_DK 1e-6 -> 1e-5 cured the silent
                # mass-balance failures seen for exactly stoichiometric PbI2 bulk at
                # 300-450 K; pa_PSTALL 0 makes a stalled MBR() fail loudly (E04IPM)
                # instead of returning OK with a wrong answer (ipm_main.cpp:879-956).
                ipm_controls={'pa_DB': 1e-17, 'pa_EPS': 1e-10, 'pa_DS': 1e-20,
                              'pa_DF': 0.01, 'pa_DFM': 0.01, 'pa_IIM': 7000,
                              'pa_PC': 2, 'pa_PRD': -5, 'pa_PSM': 1, 'pa_DG': 1000,
                              'pa_DcMin': 1e-33, 'pa_PhMin': 1e-20, 'pKin': 1,
                              'pa_DHB': 1e-10, 'pa_DK': 1e-5, 'pa_PSTALL': 0})


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('--nasa', required=True, help='NASA-Glenn thermo.inp')
    ap.add_argument('--out', required=True)
    ap.add_argument('--name', default='PbI2He')
    ap.add_argument('--liquid-class', default='l', choices=['l', 's'])
    ap.add_argument('--pbi4', action='store_true')
    ap.add_argument('--TK', type=float, default=1000.0)
    ap.add_argument('--P', type=float, default=101325.0)
    ap.add_argument('--x', type=float, default=1e-6, help='mol PbI2 per mol He')
    ap.add_argument('--json', action='store_true')
    ap.add_argument('--thermofun', action='store_true',
                    help='also write <name>-fun.json and <name>-dat-fun.lst (-o mode)')
    ap.add_argument('--no-hs', action='store_true', help='omit H0/S0/Cp0 grids')
    ap.add_argument('--dump-table', action='store_true')
    a = ap.parse_args()
    sysd = pbi2_system(a.nasa, a.liquid_class, a.pbi4, a.TK, a.P, a.x)
    d = build(sysd)
    write_kv(sysd, d, a.out, a.name, with_hs=not a.no_hs)
    if a.json:
        write_json(sysd, d, a.out, a.name, with_hs=not a.no_hs)
    if a.thermofun:
        write_thermofun_json(sysd, d, a.out, a.name)
        with open(os.path.join(a.out, a.name + '-dat-fun.lst'), 'w') as fh:
            fh.write('-o "%s-dch.dat" "%s-ipm.dat" "%s-fun.json" "%s-dbr-0-0000.dat"\n'
                     % (a.name, a.name, a.name, a.name))
    if a.dump_table:
        print('%-10s %10s %12s %10s %10s %12s' % ('DC', 'M g/mol', 'H298 J/mol',
                                                  'S298', 'Cp298', 'G(1000K)'))
        for ph, sp in d['dcs']:
            cp, h, s, g = species_props(sp, 298.15)
            g1000 = species_props(sp, 1000.0)[3]
            m = sum(sp['formula'].get(n, 0) * mm for n, mm in sysd['ICs']) * 1e3
            print('%-10s %10.4f %12.1f %10.3f %10.3f %12.1f' % (sp['name'], m, h, s, cp, g1000))
    print('wrote %s/%s-{dat.lst,dch.dat,ipm.dat,dbr-0-0000.dat}' % (a.out, a.name))


if __name__ == '__main__':
    main()
