#!/usr/bin/env python3
"""
PrimerHybridizationMap.py
12-11-2026 12:22 IST
Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA

What This Script Does
The script takes a circular DNA sequence, a forward primer, a reverse primer and the reaction
conditions. It scores every gapless register of each primer against both strands of the circle,
converts the free energies into equilibrium occupancies at the stated primer and template
concentrations, and writes one self-contained interactive HTML page. The page draws the circle as
concentric annuli with one lane per primer on each strand. The reader can change the display
threshold, hide a primer or a strand, zoom to a position, and open the alignment behind every
coloured arc, extension threat and seed.

Strand and Frame Convention
The sequence as supplied is the (+) strand. Position 1 sits at 12 o'clock and the (+) strand runs
5'->3' clockwise. Nucleotide i of the (-) ring is the base paired to (+) nucleotide i, so the (-)
strand runs 5'->3' anticlockwise. A primer binds the (+) strand when it is complementary to it,
which means it reads identical to the (-) strand. For register r the primer lies antiparallel over
(+) coordinates r .. r+n-1. The energy routine takes the primer 5'->3' as seq and the base that
physically opposes each primer position as seq2:
  binding the (+) strand   primer base k faces (+) base r+n-k
  binding the (-) strand   primer base k faces the complement of (+) base r+k-1
A forward primer reads like the (+) strand, so its cognate site lies on the (-) strand. A reverse
primer is the reverse complement of the (+) strand, so its cognate site lies on the (+) strand.
With a perfect match the energy routine reduces to the nearest-neighbour sum of a fully paired duplex.

Scoring
1. Survey. A register is scored only when it holds at least SurveyRun contiguous Watson-Crick pairs.
   The survey is exhaustive over the registers of both strands. A register whose longest run is
   under 3 holds at most one nearest-neighbour stack and is never stable.
2. Best sub-duplex. Every contiguous stretch of a register that starts and ends on a Watson-Crick
   pair and holds at least one Watson-Crick stack is scored. The register reports its most stable
   stretch and the ring colours only the nucleotides that stretch spans. Unpaired primer overhangs
   carry no energy, because dangling-end terms are outside the model.
3. Energy of a stretch. Watson-Crick stacks, helix initiation keyed on the closing pair and the
   terminal A-T term are taken from SantaLucia and Hicks (2004) Table 1 (unified parameters of
   SantaLucia 1998). An isolated internal mismatch takes the Table 2 increment, with measured
   enthalpy and entropy from Allawi and SantaLucia (1997, 1998) and Peyret et al. (1999). A run of
   two or more mismatches closed by Watson-Crick pairs is one internal loop and takes Equation 12:
   the Table 4 length term for 2k unpaired nucleotides (linear interpolation between tabulated
   sizes, Jacobson-Stockmayer extrapolation with coefficient 2.44 beyond 30) plus the terminal
   mismatch increment at each closure. The loop term is purely entropic and scales as T/310.15.
   A mismatch at a helix end is a terminal mismatch and takes the terminal-mismatch increment
   instead of a Table 2 step.
4. Salt. The Owczarzy et al. (2004, 2008) correction is applied as an entropy shift, with the
   enthalpy and the GC fraction taken from the fully paired duplex of the same sequence, so that
   mismatched registers are not credited with spurious stability. Free Mg2+ follows from Mg.dNTP
   chelation with Ka = 3E4 per molar.
5. Cache. Scores are cached by primer segment and opposing bases.

Melting Temperature of a Site
A primer over a template site is bimolecular, so it has no melting temperature until a
concentration is named. Occupancy is theta = K f/(1 + K f), so theta = 1/2 is K f = 1 and
dG(Tm) = R Tm ln(f), with f the free primer. Since dG(T) is linear in T for a fixed structure, dH
and dS follow from the same duplex evaluated at Ta and Ta + dTm, and Tm = dH/(dS + R ln f). Two
temperatures are reported, one for the free-edge optimum that colours the ring and one for the
3'-reachable duplex on which the threat circles are scored.

Equilibrium and Display Cutoff
Each scored register is a site with K = exp(-dG/RT). Each primer has its own mass balance,
      [P]0 = f + sum over sites of [T] K f/(1 + K f),
summed over the strands physically present: the (+) strand alone for a circular ssDNA template,
both strands for dsDNA. The excess primer and the limiting primer of an asymmetric reaction
therefore have different free concentrations and different display cutoffs. Registers are treated
as independent sites, which is sound while little primer is bound, and the bound fraction is
reported. A site is displayed when theta >= theta_min, which is the same statement as
dG <= -RT ln[theta_min/((1 - theta_min) f)]. Both primers share one DeltaG colour axis, so a given
shade means the same free energy in every lane; the hue names the primer.

Extension Threats
Occupancy does not say whether the polymerase can extend. Taq carries no 3'->5' exonuclease and
does not extend a primer whose 3' terminus lies outside the duplex (Huang et al. 1992). A site is a
threat when the duplex reachable from the 3' end, with the right edge pinned at the outermost pair
within EndTol bases of the terminus, has theta*epsilon >= ThreatFactor*theta_min. The efficiency
epsilon is read from the last three primer positions on the table of Rejali et al. (2018), which
prices a mismatch at its own distance from the 3' end and applies the smallest value. The
cognate site is excluded. A second perfect match is not the cognate and is listed as a threat. A
proofreading polymerase (class weak or strong) removes unpaired 3' overhangs and excises a
mispaired terminus (New England Biolabs), so every epsilon becomes 1 and EndTol widens to the whole
primer. Whether trimming outruns dissociation is kinetics and is not modelled.

Exact-Complement Seeds
Two combinatorial criteria, kept apart from the energies. A 3'-terminal run of Seed3Len perfect
pairs nucleates a duplex and leaves a paired 3'-OH. An internal run of SeedIntLen pairs lies wholly
inside primer positions 1..n-1 and holds primer at the locus until an exonuclease trims it. No
free energy is quoted for a run on its own, because helix initiation charged at both ends makes a
short run read positive. The cognate site is excluded.

Amplicon
When one primer has its cognate site on the (-) strand and the other on the (+) strand, the
amplicon runs clockwise from the 5' end of the first to the 5' end of the second and is drawn as an
arc inside the (-) ring.

References
  SantaLucia J Jr, Hicks D (2004) Annu Rev Biophys Biomol Struct 33:415-440.
      doi:10.1146/annurev.biophys.32.110601.141800
  SantaLucia J Jr (1998) Proc Natl Acad Sci USA 95:1460-1465. doi:10.1073/pnas.95.4.1460
  Owczarzy R et al. (2004) Biochemistry 43:3537-3554. doi:10.1021/bi034621r
  Owczarzy R et al. (2008) Biochemistry 47:5336-5353. doi:10.1021/bi702363u
  Allawi HT, SantaLucia J Jr (1997) Biochemistry 36:10581-10594. doi:10.1021/bi962590c
  Peyret N et al. (1999) Biochemistry 38:3468-3477. doi:10.1021/bi9825091
  Huang MM, Arnheim N, Goodman MF (1992) Nucleic Acids Res 20:4567-4573.
      doi:10.1093/nar/20.17.4567
  Rejali NA, Moric E, Wittwer CT (2018) Clin Chem 64:801-809. doi:10.1373/clinchem.2017.282285
  Veneziano R et al. (2018) Sci Rep 8:6548. doi:10.1038/s41598-018-24677-5
  M13mp18 reference sequence, GenBank X02513.1. https://www.ncbi.nlm.nih.gov/nuccore/X02513.1
"""

import argparse
import csv
import json
import math
import os
import re
import sys
from pathlib import Path

# Constants

R_GAS = 8.314e-3/4.184 # kcal/(mol K), the constant used for every conversion below
KA_MG = 3e4 # 1/M, Mg.dNTP association constant
COMP = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
COMP_TRANS = str.maketrans('ACGT', 'TGCA')
MINUS = '\u2212'
DEFAULT_FWD = 'GACGCTTACTGGAACGGTAAA'
DEFAULT_REV = 'CACATTAATTGCGTTGCGCTC' # M13mp18 6092..6112, reverse complement of the (+) strand
PLASMID_DIR = r'D:\OneDrive - Indian Institute of Science\Codes\Developed_Scripts\aPCR\Plasmids'


def revcomp(s):
    return s.translate(COMP_TRANS)[::-1]


def finite(x):
    return (x is not None) and math.isfinite(x)


# Parameter Tables

def read_csv(path):
    with open(path, 'r', encoding = 'utf-8-sig') as fh:
        raw = [ln for ln in fh.read().splitlines() if ln.strip() and not ln.strip().startswith('#')]
    rows = list(csv.reader(raw))
    if len(rows) < 2:
        raise ValueError(f'{path} holds a header but no data rows.')
    hdr = [h.strip() for h in rows[0]]
    return hdr, [[c.strip() for c in r] for r in rows[1:]]


def col_index(hdr, name, path):
    norm = lambda s: re.sub('[^a-z0-9]', '', s.lower())
    hn = [norm(h) for h in hdr]
    if norm(name) not in hn:
        raise ValueError(f'Column {name} was not found in {path}. Header as read: {hdr}')
    return hn.index(norm(name))


def to_num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float('nan')


def free_mg_speciated(biv, dntp, ka = KA_MG):
    # Free Mg2+ after 1:1 chelation by dNTP, the numerically stable root of the quadratic
    b = 1 + ka*dntp - ka*biv
    disc = math.sqrt(b*b + 4*ka*biv)
    m = (-b + disc)/(2*ka)
    if (b > 0) and ((b + disc) > 0):
        m = (2*biv)/(b + disc)
    return max(m, 0.0)


class Rxn:
    # Reaction conditions at one temperature, with the Owczarzy constants resolved once
    def __init__(self, monv, biv, dntp, tann):
        self.monv = monv
        self.biv = biv
        self.dntp = dntp
        self.tann = tann
        self.tk = tann + 273.15
        self.free_mg = free_mg_speciated(biv, dntp)
        self.lmon = math.log(monv)
        self.catrat = math.sqrt(self.free_mg)/monv
        self.lmg = math.log(self.free_mg) if self.free_mg > 0 else float('-inf')
        if self.catrat > 6.0:
            a, d, g = 3.92, 1.42, 8.31
        elif self.catrat >= 0.22:
            a = 3.92*(0.843 - 0.352*math.sqrt(monv)*self.lmon)
            d = 1.42*(1.279 - 4.03e-3*self.lmon - 8.03e-3*self.lmon**2)
            g = 8.31*(0.486 - 0.258*self.lmon + 5.25e-3*self.lmon**3)
        else:
            a = d = g = 0.0
        self.a, self.d, self.g = a, d, g

    def owczarzy_shift(self, h, gc, nbp):
        # Entropy shift in cal/(K mol) for enthalpy h in kcal/mol
        if nbp < 2:
            raise ValueError('The divalent branch divides by 2*(N-1); a segment of one pair has no neighbour.')
        if self.catrat < 0.22:
            f = (4.29*gc - 3.95)*self.lmon + 0.940*self.lmon**2
        else:
            f = (self.a - 0.911*self.lmg + gc*(6.26 + self.d*self.lmg)
                 + (1/(2*(nbp - 1)))*(-48.2 + 52.5*self.lmg + self.g*self.lmg**2))
        return h*f*1e-5*1000


class Thermo:
    # Nearest-neighbour parameters and the free energy of a pre-aligned duplex segment
    def __init__(self, folder):
        folder = Path(folder)
        steps = ['AA', 'TT', 'AT', 'TA', 'CA', 'TG', 'GT', 'AC', 'CT', 'AG', 'GA', 'TC', 'CG', 'GC', 'GG', 'CC']
        dh = [-7.90, -7.90, -7.20, -7.20, -8.50, -8.50, -8.40, -8.40,
              -7.80, -7.80, -8.20, -8.20, -10.60, -9.80, -8.00, -8.00]
        ds = [-22.20, -22.20, -20.40, -21.30, -22.70, -22.70, -22.40, -22.40,
              -21.00, -21.00, -22.20, -22.20, -27.20, -24.40, -19.90, -19.90]
        self.wc_h = dict(zip(steps, dh))
        self.wc_s = dict((k, 0.001*v) for k, v in zip(steps, ds)) # kcal/(K mol)
        self.init_gc = (0.10, -2.80) # dH kcal/mol, dS cal/(K mol)
        self.init_at = (2.30, 4.10)
        self.sym = (0.00, -1.40)
        self.load_mismatch(folder/'mismatch_nn_params.csv')
        self.load_terminal(folder/'terminal_mismatch_nn_params.csv')
        self.load_loops(folder/'loop_length_params.csv')
        self._loop_cache = {}

    def load_mismatch(self, path):
        hdr, rows = read_csv(path)
        c = dict((k, col_index(hdr, k, path)) for k in
                 ['Closing', 'X', 'Y', 'dG37_kcal_per_mol', 'dH_kcal_per_mol', 'dS_cal_per_K_per_mol'])
        self.mm = {}
        for r in rows:
            key = (r[c['Closing']].upper(), r[c['X']].upper(), r[c['Y']].upper())
            h = to_num(r[c['dH_kcal_per_mol']])
            s = to_num(r[c['dS_cal_per_K_per_mol']])
            if not (math.isfinite(h) and math.isfinite(s)):
                raise ValueError(f'Mismatch entry {key} in {path} has no dH and dS. The temperature slope '
                                 'from the Vienna table is not supported here; fill dH and dS for all 48 rows.')
            self.mm[key] = (h, 0.001*s)
        if len(self.mm) != 48:
            raise ValueError(f'{path} supplied {len(self.mm)} mismatch entries where 48 were expected.')

    def load_terminal(self, path):
        hdr, rows = read_csv(path)
        names = ['i1', 'i2', 'j1', 'j2', 'dH_kcal_per_mol', 'dS_cal_per_K_per_mol', 'Suspect']
        c = dict((k, col_index(hdr, k, path)) for k in names)
        self.tm = {}
        ok = []
        for r in rows:
            key = tuple(r[c[k]].upper() for k in ['i1', 'i2', 'j1', 'j2'])
            h = to_num(r[c['dH_kcal_per_mol']])
            s = to_num(r[c['dS_cal_per_K_per_mol']])
            sus = to_num(r[c['Suspect']]) != 0
            self.tm[key] = (h, s, sus)
            if (not sus) and math.isfinite(h) and math.isfinite(s):
                ok.append(h - 310.15*s/1000)
        self.tm_cap = max(ok) if ok else 0.0

    def load_loops(self, path):
        hdr, rows = read_csv(path)
        c = dict((k, col_index(hdr, k, path)) for k in ['LoopLength', 'Internal_dG37'])
        pairs = sorted((to_num(r[c['LoopLength']]), to_num(r[c['Internal_dG37']])) for r in rows)
        self.loop_n = [int(p[0]) for p in pairs]
        self.loop_internal = [p[1] for p in pairs]

    def loop_internal_dg(self, n, tann):
        # Internal loop of n unpaired nucleotides on both strands, Table 4 at the annealing temperature
        key = (n, tann)
        if key in self._loop_cache:
            return self._loop_cache[key]
        sizes, tab = self.loop_n, self.loop_internal
        if n in sizes:
            g = tab[sizes.index(n)]
            if math.isnan(g):
                raise ValueError(f'Loop length {n} is not allowed for internal loops.')
        elif n < sizes[-1]:
            lo = max(x for x in sizes if x < n)
            hi = min(x for x in sizes if x > n)
            g_lo, g_hi = tab[sizes.index(lo)], tab[sizes.index(hi)]
            g = g_lo + (g_hi - g_lo)*(n - lo)/(hi - lo)
        else:
            g = tab[-1] + 2.44*1.9872e-3*310.15*math.log(n/sizes[-1])
        g = g*((tann + 273.15)/310.15)
        self._loop_cache[key] = g
        return g

    def term_mismatch(self, i1, i2, j1, j2, tk):
        # Terminal mismatch closing a loop or a helix end: returns dH kcal/mol and dS cal/(K mol)
        e = self.tm.get((i1, i2, j1, j2))
        if e is None:
            return 0.0, 0.0
        h, s, sus = e
        if not (math.isfinite(h) and math.isfinite(s)):
            return 0.0, 0.0
        if sus:
            g = h - tk*s/1000
            if g > self.tm_cap:
                g = self.tm_cap
                s = 1000*(h - g)/tk
        return h, s

    def dg_mismatch(self, seq, seq2, rx, apply_init = True):
        # Free energy in kcal/mol of a pre-aligned segment, seq 5'->3' over seq2 read 3'->5'
        L = len(seq)
        tk = rx.tk
        iswc = [COMP[a] == b for a, b in zip(seq, seq2)]
        h = 0.0
        s = 0.0
        dg_fixed = 0.0
        skip5 = not iswc[0]
        skip3 = not iswc[L - 1]

        # Runs of two or more mismatches closed by pairs at both ends belong to the loop term
        eq12 = [False]*L
        i = 0
        while i < L:
            if not iswc[i]:
                j = i
                while (j < L - 1) and (not iswc[j + 1]):
                    j += 1
                if ((j - i + 1) >= 2) and (i > 0) and (j < L - 1):
                    for q in range(i, j + 1):
                        eq12[q] = True
                i = j + 1
            else:
                i += 1

        for i in range(L - 1):
            at5 = skip5 and (i == 0)
            at3 = skip3 and (i == L - 2)
            if at5 or at3:
                if (at5 and (not iswc[0]) and iswc[1]) or (at3 and (not iswc[L - 1]) and iswc[L - 2]):
                    continue
            a_wc, b_wc = iswc[i], iswc[i + 1]
            if a_wc and b_wc:
                step = seq[i:i + 2]
                h += self.wc_h[step]
                s += self.wc_s[step]
            elif a_wc and (not b_wc):
                if eq12[i + 1]:
                    continue
                e = self.mm.get((seq[i], seq[i + 1], seq2[i + 1]))
                if e is None:
                    raise ValueError(f'No Table 2 parameter for closing {seq[i]}, X {seq[i + 1]}, Y {seq2[i + 1]}.')
                h += e[0]
                s += e[1]
            elif (not a_wc) and b_wc:
                if eq12[i]:
                    continue
                e = self.mm.get((seq2[i + 1], seq2[i], seq[i]))
                if e is None:
                    raise ValueError(f'No Table 2 parameter for closing {seq2[i + 1]}, X {seq2[i]}, Y {seq[i]}.')
                h += e[0]
                s += e[1]

        # Internal loops of Equation 12, one per run, with a terminal mismatch at each closure
        dg_runs = 0.0
        i = 0
        while i < L:
            if not iswc[i]:
                j = i
                while (j < L - 1) and (not iswc[j + 1]):
                    j += 1
                runlen = j - i + 1
                if (runlen >= 2) and (i > 0) and (j < L - 1):
                    dg_runs += self.loop_internal_dg(2*runlen, rx.tann)
                    h_l, s_l = self.term_mismatch(seq[i - 1], seq[i], seq2[i - 1], seq2[i], tk)
                    h_r, s_r = self.term_mismatch(seq2[j + 1], seq2[j], seq[j + 1], seq[j], tk)
                    h += h_l + h_r
                    s += 0.001*(s_l + s_r)
                i = j + 1
            else:
                i += 1

        # Helix initiation on the outermost pair of each end, and a stacking term for a terminal mismatch
        q5 = next((k for k in range(L) if iswc[k]), None)
        q3 = next((k for k in range(L - 1, -1, -1) if iswc[k]), None)
        if q5 is not None:
            init = self.init_gc if seq[q5] in 'GC' else self.init_at
            h += init[0]
            s += 0.001*init[1]
            if q5 == 1:
                h_t, s_t = self.term_mismatch(seq2[1], seq2[0], seq[1], seq[0], tk)
                h += h_t
                s += 0.001*s_t
        if q3 is not None:
            init = self.init_gc if seq[q3] in 'GC' else self.init_at
            h += init[0]
            s += 0.001*init[1]
            if q3 == L - 2:
                h_t, s_t = self.term_mismatch(seq[L - 2], seq[L - 1], seq2[L - 2], seq2[L - 1], tk)
                h += h_t
                s += 0.001*s_t

        # Symmetry correction, owed when the two strands of the complex are identical
        if apply_init and (seq == seq2[::-1]):
            s += 0.001*self.sym[1]
            h += self.sym[0]

        # Salt, enthalpy and GC fraction of the fully paired duplex of the same sequence
        h_ref = 0.0
        for i in range(L - 1):
            h_ref += self.wc_h[seq[i:i + 2]]
        h_ref += self.init_gc[0] if seq[0] in 'GC' else self.init_at[0]
        h_ref += self.init_gc[0] if seq[L - 1] in 'GC' else self.init_at[0]
        gc_ref = sum(1 for ch in seq if ch in 'GC')/L
        s += 0.001*rx.owczarzy_shift(h_ref, gc_ref, L)

        return h - tk*s + dg_fixed + dg_runs


class ExtensionEfficiency:
    # Extension efficiency of the 3' region from the last positions of the primer, Rejali et al. 2018
    def __init__(self, path):
        hdr, rows = read_csv(path)
        ip = col_index(hdr, 'PrimerBase', path)
        io = col_index(hdr, 'OpposingBase', path)
        ie = col_index(hdr, 'Epsilon', path)
        grid = []
        for k, h in enumerate(hdr):
            m = re.fullmatch(r'Eps_pos(\d+)', h)
            if m and int(m.group(1)) >= 1:
                grid.append((int(m.group(1)), k))
        grid.sort()
        self.grid = [g[0] for g in grid]
        self.logs = {}
        self.cite = ''
        icite = hdr.index('Citation') if 'Citation' in hdr else None
        for r in rows:
            key = (r[ip].upper(), r[io].upper())
            ys = [math.log10(float(r[ie]))] + [math.log10(float(r[k])) for _, k in grid]
            self.logs[key] = ys
            if (not self.cite) and (icite is not None) and r[icite]:
                self.cite = r[icite]
        self.source = str(path)
        for b in 'ACGT':
            wc = self.logs[(b, COMP[b])]
            if any(abs(v) > 1e-12 for v in wc):
                raise ValueError(f'The Watson-Crick cells of {path} must be exactly 1.')
        ct = 10**self.logs[('C', 'T')][0]
        tc = 10**self.logs[('T', 'C')][0]
        if ct < tc:
            raise ValueError(f'{path} appears transposed: primer C opposite T reads {ct:.4g} and primer T '
                             f'opposite C reads {tc:.4g}.')
        self._memo = {}

    def at_distance(self, pb, ob, d):
        ys = self.logs[(pb, ob)]
        if d <= 0:
            return 10**ys[0]
        if not self.grid:
            return 1.0
        if d >= self.grid[-1]:
            return 10**ys[-1]
        xs = [0] + self.grid
        for q in range(len(xs) - 1):
            if xs[q] <= d <= xs[q + 1]:
                t = (d - xs[q])/(xs[q + 1] - xs[q])
                return 10**(ys[q] + t*(ys[q + 1] - ys[q]))
        return 1.0

    def __call__(self, p_bases, t_bases):
        # Smallest efficiency over the mismatched positions, 3' end last
        key = (p_bases, t_bases)
        if key in self._memo:
            return self._memo[key]
        n = len(p_bases)
        eps = 1.0
        gov_pair = p_bases[-1] + '.' + t_bases[-1]
        gov_pos = 0
        for j in range(n - 1, -1, -1):
            if COMP[p_bases[j]] == t_bases[j]:
                continue
            e = self.at_distance(p_bases[j], t_bases[j], n - 1 - j)
            if e < eps:
                eps = e
                gov_pair = p_bases[j] + '.' + t_bases[j]
                gov_pos = n - 1 - j
        label = f'{gov_pair}@{gov_pos}' if gov_pos > 0 else gov_pair
        self._memo[key] = (eps, label)
        return self._memo[key]


def bimolecular_tm(g1, g2, t1, t2, ceff, r = R_GAS):
    # Tm = dH/(dS + R ln ceff) from the free energy of one structure at two temperatures
    if not (finite(g1) and finite(g2)):
        return float('nan')
    k1, k2 = t1 + 273.15, t2 + 273.15
    ds = (g1 - g2)/(k2 - k1)
    dh = g1 + k1*ds
    den = ds + r*math.log(ceff)
    if abs(den) < 1e-9:
        return float('nan')
    return dh/den - 273.15


# Register Scan

def scan_strand(plus, primer, strand, th, rx, rx2, ee, opt, cache):
    # Every register of the primer against one strand, scored as its most stable Watson-Crick stretch
    L = len(plus)
    n = len(primer)
    plus_comp = plus.translate(COMP_TRANS)
    plus_ext = plus + plus[:n]
    pcomp = [COMP[c] for c in primer]
    seed_min = min(opt['seed3'], opt['seedint'])
    keep_run = min(opt['survey_run'], seed_min)
    end_tol = opt['end_tol']
    sites = []
    for r in range(L):
        win = plus_ext[r:r + n]
        if strand == '+':
            ip = [w == c for w, c in zip(reversed(win), pcomp)]
        else:
            ip = [w == c for w, c in zip(win, primer)]
        run = 0
        mx = 0
        for k in range(n):
            run = run + 1 if ip[k] else 0
            if run > mx:
                mx = run
        if mx < keep_run:
            continue
        end3 = 0
        for k in range(n - 1, -1, -1):
            if ip[k]:
                end3 += 1
            else:
                break
        int_run = 0
        int_end = 0
        rl = 0
        for k in range(n - 1):
            rl = rl + 1 if ip[k] else 0
            if rl > int_run:
                int_run = rl
                int_end = k + 1
        if strand == '+':
            idx = [(r + n - 1 - k)%L for k in range(n)]
            opp = ''.join(plus[i] for i in idx)
        else:
            idx = [(r + k)%L for k in range(n)]
            opp = ''.join(plus_comp[i] for i in idx)

        # Prefix counts of adjacent Watson-Crick stacks, so that a stretch can be tested in one step
        cs = [0]*(n + 1)
        for k in range(n - 1):
            cs[k + 1] = cs[k] + (1 if (ip[k] and ip[k + 1]) else 0)
        cs[n] = cs[n - 1] if n > 0 else 0
        pairs = [k for k in range(n) if ip[k]]
        if opt['exact']:
            e_left = pairs
            e_right = pairs
        else:
            e_left = [k for k in range(n) if ip[k] and ((k == 0) or (not ip[k - 1]))]
            e_right = [k for k in range(n) if ip[k] and ((k == n - 1) or (not ip[k + 1]))]

        best = math.inf
        bf = -1
        bl = -1
        for f in e_left:
            for l in e_right:
                if (l < f + 1) or (cs[l] - cs[f] == 0):
                    continue
                key = primer[f:l + 1] + '/' + opp[f:l + 1]
                g = cache.get(key)
                if g is None:
                    g = th.dg_mismatch(primer[f:l + 1], opp[f:l + 1], rx, opt['apply_init'])
                    cache[key] = g
                if g < best:
                    best, bf, bl = g, f, l

        # Duplex reachable from the 3' end, right edge pinned to the outermost pair
        ga = math.inf
        af = -1
        j_r = pairs[-1] if pairs else -1
        if (j_r >= 0) and (j_r >= n - 1 - end_tol):
            for f in e_left:
                if (f > j_r - 1) or (cs[j_r] - cs[f] == 0):
                    continue
                key = primer[f:j_r + 1] + '/' + opp[f:j_r + 1]
                g = cache.get(key)
                if g is None:
                    g = th.dg_mismatch(primer[f:j_r + 1], opp[f:j_r + 1], rx, opt['apply_init'])
                    cache[key] = g
                if g < ga:
                    ga, af = g, f
        if bf < 0:
            continue

        st = {'p':0, 'd':strand, 'reg':r + 1, 'pf':bf + 1, 'pt':bl + 1, 'n':bl - bf + 1, 'g':best,
              'tp':bool(ip[n - 1]),
              'seg':primer[bf:bl + 1], 'op':opp[bf:bl + 1],
              'pr':''.join('|' if ip[k] else '.' for k in range(bf, bl + 1)),
              'perfect':(sum(ip) == n), 'e3':end3, 'ir':int_run,
              'an':False, 'ga':math.inf, 'eps':1.0, 'ep':'', 'gae':math.inf}
        st['g2'] = th.dg_mismatch(primer[bf:bl + 1], opp[bf:bl + 1], rx2, opt['apply_init'])
        if strand == '+':
            st['a'], st['b'] = idx[bl] + 1, idx[bf] + 1
        else:
            st['a'], st['b'] = idx[bf] + 1, idx[bl] + 1

        # Seed footprints on the ring, in the same frame as the start and end of the site
        if end3 > 0:
            f3 = n - end3
            st['s3a'] = (idx[n - 1] + 1) if strand == '+' else (idx[f3] + 1)
            st['s3b'] = (idx[f3] + 1) if strand == '+' else (idx[n - 1] + 1)
            st['s3m'] = (st['s3a'] - 1 + end3/2)%L
        if int_run > 0:
            li = int_end - 1
            fi = li - int_run + 1
            st['ifr'], st['ito'] = fi + 1, li + 1
            st['isa'] = (idx[li] + 1) if strand == '+' else (idx[fi] + 1)
            st['isb'] = (idx[fi] + 1) if strand == '+' else (idx[li] + 1)
            st['im'] = (st['isa'] - 1 + int_run/2)%L

        if af >= 0:
            st['an'] = True
            st['ga'] = ga
            st['ga2'] = th.dg_mismatch(primer[af:j_r + 1], opp[af:j_r + 1], rx2, opt['apply_init'])
            st['af'], st['at'] = af + 1, j_r + 1
            st['al'] = j_r - af + 1
            st['as'] = (idx[j_r] + 1) if strand == '+' else (idx[af] + 1)
            st['ae'] = (idx[af] + 1) if strand == '+' else (idx[j_r] + 1)
            st['am'] = (st['as'] - 1 + st['al']/2)%L
            st['aseg'] = primer[af:j_r + 1]
            st['aop'] = opp[af:j_r + 1]
            st['apr'] = ''.join('|' if ip[k] else '.' for k in range(af, j_r + 1))
            st['tpa'] = bool(ip[n - 1])
            k0 = max(0, n - 3)
            if opt['use_eps']:
                st['eps'], st['ep'] = ee(primer[k0:], opp[k0:])
            else:
                st['eps'], st['ep'] = 1.0, primer[-1] + '.' + opp[n - 1]
            st['gae'] = ga - R_GAS*tk_of(rx)*math.log(st['eps'])
        sites.append(st)
    return sites


def tk_of(rx):
    return rx.tk


# Equilibrium, Regions and Report

def solve_free_fraction(k_present, c_t, p0):
    # Free primer as a fraction of the total, from the mass balance over all independent sites
    ratio = c_t/p0
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5*(lo + hi)
        tot = mid
        for k in k_present:
            x = k*p0*mid
            tot += ratio*x/(1 + x)
        if tot > 1:
            hi = mid
        else:
            lo = mid
    return 0.5*(lo + hi)


def analyse_primer(plus, spec, th, rx, rx2, ee, opt, cache, tmpl, notes):
    # Scan both strands for one primer and resolve occupancy, melting temperatures and flags
    primer = spec['seq']
    n = len(primer)
    po = dict(opt)
    po['end_tol'] = (n - 1) if opt['proofread'] else opt['end_tol']
    print(f"Scanning the {spec['name']} primer, {n} nt: (+) strand ...", flush = True)
    sites = scan_strand(plus, primer, '+', th, rx, rx2, ee, po, cache)
    print(f"  (+) strand done, {len(sites)} registers kept; scanning the (-) strand ...", flush = True)
    sites += scan_strand(plus, primer, '-', th, rx, rx2, ee, po, cache)
    print(f"  {len(sites)} registers kept in all, {len(cache)} distinct duplexes evaluated so far.", flush = True)

    rt = R_GAS*rx.tk
    p0 = spec['conc_nM']*1e-9
    for s in sites:
        s['p'] = spec['index']
        s['K'] = math.exp(-s['g']/rt)
    present = [s['K'] for s in sites if (s['d'] == '+') or (tmpl['strands'] == 2)]
    u_free = solve_free_fraction(present, tmpl['cT'], p0)
    f_free = u_free*p0
    bound = 1 - u_free
    if bound > 0.01:
        notes.append(f"{spec['name']} primer: {100*bound:.2f}% of the primer is bound at equilibrium. Registers are "
                     'treated as independent sites, but overlapping registers on one template molecule exclude '
                     'each other, so the free primer is underestimated. Raise the primer concentration or lower '
                     'the template concentration to keep this below 1%.')
    for s in sites:
        s['tm'] = bimolecular_tm(s['g'], s['g2'], rx.tann, rx2.tann, f_free)
        s['tma'] = bimolecular_tm(s['ga'], s.get('ga2', math.inf), rx.tann, rx2.tann, f_free) if s['an'] else float('nan')
        s['th'] = 1/(1 + 1/(s['K']*f_free)) if s['K']*f_free > 0 else 0.0
        if s['an']:
            ka = math.exp(-s['ga']/rt)
            s['thA'] = 1/(1 + 1/(ka*f_free))
        else:
            s['thA'] = 0.0
        s['rate'] = s['thA']*s['eps']
        s['cg'] = False
    perfect = [k for k, s in enumerate(sites) if s['perfect']]
    if perfect:
        sites[perfect[0]]['cg'] = True
        if len(perfect) > 1:
            loc = ', '.join(f"({sites[k]['d']}) {sites[k]['a']}-{sites[k]['b']}" for k in perfect)
            notes.append(f"{spec['name']} primer matches the template perfectly at {len(perfect)} places: {loc}. "
                         'Every full-length perfect duplex has the same free energy, so none can be singled out as '
                         'the intended site on thermodynamics. The first is the cognate and the rest are listed as '
                         'extension threats. A primer that is not unique cannot give a clean asymmetric reaction.')
    return {'sites':sites, 'f_free':f_free, 'u_free':u_free, 'bound':bound, 'rt':rt, 'p0':p0,
            'end_tol':po['end_tol'], 'n':n}


def region_table(best, best_site, length):
    # Connected runs of displayed nucleotides on one circular strand, with the wrap from L to 1
    mask = [math.isfinite(v) for v in best]
    regs = []
    if not any(mask):
        return regs
    if all(mask):
        order = list(range(length))
        runs = [(0, length)]
    else:
        k0 = mask.index(False)
        order = [(k0 + i)%length for i in range(length)]
        runs = []
        i = 0
        while i < length:
            if mask[order[i]]:
                j = i
                while (j + 1 < length) and mask[order[j + 1]]:
                    j += 1
                runs.append((i, j + 1))
                i = j + 1
            else:
                i += 1
    for a, b in runs:
        idxs = order[a:b]
        vals = [best[x] for x in idxs]
        j = vals.index(min(vals))
        s = idxs[0] + 1
        e = idxs[-1] + 1
        ln = (e - s)%length + 1
        regs.append({'s':s, 'e':e, 'len':ln, 'g':vals[j], 'k':best_site[idxs[j]], 'x':(s - 1 + ln/2)%length})
    return regs


def fmt_signed(x, d = 2):
    s = f'{x:.{d}f}'
    return MINUS + s[1:] if s.startswith('-') else s


def report_primer(spec, res, plus_len, opt, theta_pct, n_strands):
    # Console summary of one primer, the counterpart of the figure
    sites = res['sites']
    n = res['n']
    rt = res['rt']
    f = res['f_free']
    theta_min = theta_pct/100
    cut = -rt*math.log(theta_min/((1 - theta_min)*f))
    thr_min = theta_min*opt['threat_factor']
    shown = [k for k, s in enumerate(sites) if s['th'] >= theta_min]
    is_thr = [s['an'] and (not s['cg']) and (s['rate'] >= thr_min) for s in sites]
    ord_thr = sorted([k for k in range(len(sites)) if is_thr[k]], key = lambda k:-sites[k]['rate'])
    print('\n' + '='*100)
    print(f"{spec['name'].upper()} PRIMER 5'-{spec['seq']}-3' ({n} nt), {spec['conc_nM']:g} nM, "
          f"free primer {f*1e9:.5g} nM ({100*res['bound']:.3g}% bound)")
    print(f"Display cutoff: occupancy >= {theta_pct:.3g}%, i.e. dG <= {cut:.2f} kcal/mol at {opt['ta']:.1f} C.")
    print(f"End tolerance {res['end_tol']}{' (whole primer, proofreading)' if opt['proofread'] else ''}.")
    print(f"Registers scored: {sum(1 for s in sites if s['d'] == '+')} on (+), "
          f"{sum(1 for s in sites if s['d'] == '-')} on (-).")
    print(f"Displayed sites: {len(shown)}; extension threats: {sum(is_thr)}.")
    kc = [k for k, s in enumerate(sites) if s['cg']]
    for k in kc:
        s = sites[k]
        print(f"  Cognate ({s['d']}) {s['a']}-{s['b']}  dG {s['g']:.2f} kcal/mol  Tm {s['tm']:.2f} C  theta {s['th']:.4g}")
    print('\nMost stable displayed sites (whole duplex, the 3\' end may be unpaired)')
    for k in sorted(shown, key = lambda k:sites[k]['g'])[:10]:
        s = sites[k]
        print(f"  ({s['d']}) {s['a']:5d}-{s['b']:<5d}  dG {s['g']:7.2f}  Tm {s['tm']:7.2f} C  theta {s['th']:.4g}  "
              f"3' base paired: {int(s['tp'])}  primer {s['pf']}-{s['pt']}")
        print(f"        5'-{s['seg']}-3'  primer\n           {s['pr']}\n        3'-{s['op']}-5'  template")
    print('\nExtension threats (3\'-reachable duplex), ranked by theta*epsilon')
    if not ord_thr:
        print(f'  none at theta*eps >= {100*thr_min:.3g}%.')
    for k in ord_thr[:opt['threat_max']]:
        s = sites[k]
        print(f"  ({s['d']}) {s['as']:5d}-{s['ae']:<5d}  dG {s['ga']:7.2f}  Tm {s['tma']:7.2f}  dGeff {s['gae']:7.2f}  "
              f"theta {s['thA']:9.3g}  eps {s['eps']:8.3g}  theta*eps {s['rate']:9.3g}")
        print(f"        primer {s['af']}-{s['at']} of {n}, eps from {s['ep']}, "
              f"3' base {'paired' if s['tpa'] else 'MISPAIRED'}"
              f"{'' if s['th'] >= theta_min else ', below the display cutoff'}")
        print(f"        5'-{s['aseg']}-3'  primer\n           {s['apr']}\n        3'-{s['aop']}-5'  template")
    is3 = [(s['e3'] >= opt['seed3']) and (not s['cg']) for s in sites]
    isi = [(s['ir'] >= opt['seedint']) and (not s['cg']) for s in sites]
    rank = {}
    for k in range(len(sites)):
        if isi[k]:
            rank[k] = sites[k]['ir']
        if is3[k]:
            rank[k] = 1000 + sites[k]['e3']
    ord_seed = sorted(rank, key = lambda k:-rank[k])
    print(f"\nExact-complement seeds: {len(rank)} sites ({sum(is3)} with a 3' run >= {opt['seed3']}, "
          f"{sum(isi)} with an internal run >= {opt['seedint']}); the cognate site is excluded.")
    for k in ord_seed[:opt['seed_max']]:
        s = sites[k]
        cls = 'both' if (is3[k] and isi[k]) else ("3' run" if is3[k] else 'internal')
        print(f"  ({s['d']}) {cls:9s} 3' run {s['e3']}, internal run {s['ir']}, site dG {s['g']:.2f}, "
              f"Tm {s['tm']:.2f}{', threat' if is_thr[k] else ''}")


# Input Handling

def parse_sequence_text(txt, fallback):
    # Plain text, fasta (first record) or GenBank (origin block); every non-letter is dropped
    name = fallback
    ln = re.split(r'\r\n|\n|\r', txt)
    lt = [x.strip() for x in ln]
    if any(x.startswith('LOCUS') for x in lt) and any(x.startswith('ORIGIN') for x in lt):
        tok = next(x for x in lt if x.startswith('LOCUS')).split()
        if len(tok) >= 2:
            name = tok[1]
        o = next(i for i, x in enumerate(lt) if x.startswith('ORIGIN'))
        e = next((i for i, x in enumerate(lt) if x.startswith('//') and i > o), len(lt))
        body = ''.join(ln[o + 1:e])
    elif any(x.startswith('>') for x in lt):
        h = [i for i, x in enumerate(lt) if x.startswith('>')]
        hdr = lt[h[0]][1:].strip()
        if hdr:
            name = hdr.split()[0]
        stop = h[1] if len(h) > 1 else len(lt)
        body = ''.join(ln[h[0] + 1:stop])
    else:
        body = txt
    return re.sub('[^A-Za-z]', '', body).upper().replace('U', 'T'), name


def ask(val, msg, cast = str, default = None):
    # Command-line value first, then a prompt when a terminal is attached, then the default
    if val is not None:
        return val
    if sys.stdin.isatty():
        shown = '' if default is None else f' [{default}]'
        raw = input(f'{msg}{shown}: ').strip()
        return cast(raw) if raw else default
    return default


def ask_sequence_file(given):
    if given:
        return given
    if sys.stdin.isatty():
        try:
            import tkinter
            from tkinter import filedialog
            root = tkinter.Tk()
            root.withdraw()
            path = filedialog.askopenfilename(
                title = 'Select the circular DNA sequence',
                initialdir = PLASMID_DIR if os.path.isdir(PLASMID_DIR) else os.getcwd(),
                filetypes = [('Sequence files', '*.txt *.fa *.fasta *.fna *.seq *.gb *.gbk'), ('All files', '*.*')])
            root.destroy()
            if path:
                return path
        except Exception:
            pass
        return input('Circular DNA sequence file (plain text, FASTA or GenBank): ').strip()
    raise SystemExit('A sequence file is required: pass --seq-file.')


def check_number(x, label, lo, hi):
    if (x is None) or (not math.isfinite(x)) or (x < lo) or (x > hi):
        raise SystemExit(f'{label} must be a number in [{lo:g}, {hi:g}]; received {x}.')
    return x


# Page Data and Output

def jnum(x, nd = 5):
    return round(x, nd) if finite(x) else None


def site_json(s):
    out = {'p':s['p'], 'd':s['d'], 'a':s['a'], 'b':s['b'], 'n':s['n'], 'g':jnum(s['g'], 4),
           'tm':jnum(s['tm'], 3), 'pf':s['pf'], 'pt':s['pt'], 'tp':int(s['tp']), 'seg':s['seg'], 'pr':s['pr'],
           'op':s['op'], 'cg':int(s['cg']), 'an':int(s['an']), 'e3':s['e3'], 'ir':s['ir']}
    if s['e3'] > 0:
        out.update({'s3a':s['s3a'], 's3b':s['s3b'], 's3m':jnum(s['s3m'], 2)})
    if s['ir'] > 0:
        out.update({'ifr':s['ifr'], 'ito':s['ito'], 'isa':s['isa'], 'isb':s['isb'], 'im':jnum(s['im'], 2)})
    if s['an']:
        out.update({'ga':jnum(s['ga'], 4), 'tma':jnum(s['tma'], 3), 'as':s['as'], 'ae':s['ae'], 'al':s['al'],
                    'af':s['af'], 'at':s['at'], 'am':jnum(s['am'], 2), 'eps':float(f"{s['eps']:.6g}"),
                    'ep':s['ep'], 'gae':jnum(s['gae'], 4), 'tpa':int(s['tpa']), 'aseg':s['aseg'],
                    'apr':s['apr'], 'aop':s['aop']})
    return out


def keep_site(s, primer_res, floor, opt):
    # Embed a site only when some control of the page could display it
    if s['th'] >= floor:
        return True
    if s['an'] and (s['rate'] >= floor*opt['threat_factor']):
        return True
    if (s['e3'] >= opt['seed3']) or (s['ir'] >= opt['seedint']) or s['perfect']:
        return True
    return False


def build_html(out_path, payload):
    data = json.dumps(payload, separators = (',', ':'), allow_nan = False).replace('<', '\\u003c')
    html = HTML_TEMPLATE.replace('__TITLE__', f"Primer Hybridization Map: {payload['meta']['name']}")
    html = html.replace('__PAYLOAD__', data)
    with open(out_path, 'w', encoding = 'utf-8') as fh:
        fh.write(html)


HTML_TEMPLATE = r'''<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
:root{--bg:#ffffff;--fg:#1f2937;--muted:#6b7280;--card:#f8fafc;--line:#d1d5db;--unbound:#e6e6e6;--outline:#8c8c8c;--tick:#737373;--accent:#2563eb;--amp:#7c3aed;--tip:#111827;--tipfg:#f9fafb}
:root[data-theme=dark]{--bg:#0f172a;--fg:#e5e7eb;--muted:#94a3b8;--card:#1e293b;--line:#334155;--unbound:#2f3b52;--outline:#64748b;--tick:#94a3b8;--accent:#60a5fa;--amp:#a78bfa;--tip:#f9fafb;--tipfg:#111827}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:14px/1.45 system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
header{padding:14px 20px 6px}
header h1{margin:0;font-size:1.25rem}
header .sub{color:var(--muted);font-size:.85rem;margin-top:2px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px 14px}
#controls{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:10px 16px;margin:8px 20px;align-items:start}
#controls label{display:block;font-size:.72rem;text-transform:uppercase;letter-spacing:.03em;color:var(--muted);font-weight:600;margin-bottom:3px}
#controls .row{display:flex;gap:10px;flex-wrap:wrap;align-items:center}
#controls input[type=range]{width:100%}
#controls input[type=number],#controls input[type=text]{width:64px;padding:3px 5px;background:var(--bg);color:var(--fg);border:1px solid var(--line);border-radius:5px}
#controls input[type=text]{width:110px}
button{background:var(--bg);color:var(--fg);border:1px solid var(--line);border-radius:6px;padding:4px 10px;cursor:pointer;font:inherit}
button:hover{border-color:var(--accent)}
main{display:grid;grid-template-columns:minmax(0,1.55fr) minmax(330px,1fr);gap:14px;margin:8px 20px 24px}
@media (max-width:1100px){main{grid-template-columns:1fr}}
#mapcard{padding:8px}
#map{width:100%;height:auto;display:block;cursor:grab;touch-action:none;user-select:none}
#map.drag{cursor:grabbing}
#legend{display:flex;flex-wrap:wrap;gap:8px 22px;padding:6px 8px 2px;align-items:flex-start;font-size:.78rem;color:var(--muted)}
.bar{width:210px;height:12px;border:1px solid var(--outline);border-radius:3px}
.barticks{display:flex;justify-content:space-between;width:210px;font-size:.7rem}
.key{display:flex;gap:6px;align-items:center;white-space:nowrap}
.tabs{display:flex;gap:4px;flex-wrap:wrap;margin-bottom:8px}
.tabs button{padding:3px 9px;font-size:.82rem}
.tabs button.on{background:var(--accent);border-color:var(--accent);color:#fff}
.tabp{display:none}.tabp.on{display:block}
#side{min-width:0}
.tblwrap{max-height:520px;overflow:auto;border:1px solid var(--line);border-radius:6px}
table{border-collapse:collapse;width:100%;font-size:.76rem;font-variant-numeric:tabular-nums}
th,td{padding:4px 7px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}
th{position:sticky;top:0;background:var(--card);color:var(--muted)}
tr.r:hover,tr.r.on{background:rgba(37,99,235,.15);cursor:pointer}
pre.aln{margin:6px 0;padding:8px 10px;background:var(--bg);border:1px solid var(--line);border-radius:6px;font:12px/1.35 ui-monospace,Menlo,Consolas,monospace;overflow:auto}
.kv{display:grid;grid-template-columns:auto 1fr;gap:2px 12px;font-size:.82rem}
.kv span:nth-child(odd){color:var(--muted)}
.dot{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:5px;vertical-align:-1px}
.warn{border-left:3px solid #d97706;padding:4px 10px;margin:6px 0;background:rgba(217,119,6,.10);font-size:.82rem}
#tip{position:fixed;z-index:10;pointer-events:none;background:var(--tip);color:var(--tipfg);padding:7px 10px;border-radius:6px;font-size:.76rem;max-width:430px;display:none;box-shadow:0 4px 14px rgba(0,0,0,.25)}
#tip pre{margin:4px 0 0;font:11.5px/1.3 ui-monospace,Menlo,Consolas,monospace}
.mono{font-family:ui-monospace,Menlo,Consolas,monospace}
</style></head><body>
<header>
  <h1 id="ttl"></h1>
  <div class="sub" id="sub"></div>
</header>
<section id="controls" class="card">
  <div><label>Display threshold, equilibrium occupancy</label>
    <input type="range" id="thr" min="0" max="1000" step="1">
    <div id="thrTxt" class="mono"></div></div>
  <div><label>Show</label>
    <div class="row">
      <span><input type="checkbox" id="shF" checked> <span id="shFl">Forward</span></span>
      <span><input type="checkbox" id="shR" checked> <span id="shRl">Reverse</span></span>
      <span><input type="checkbox" id="shP" checked> (+) lanes</span>
      <span><input type="checkbox" id="shM" checked> (&minus;) lanes</span>
      <span><input type="checkbox" id="shT" checked> threats</span>
    </div></div>
  <div><label>Labels per lane or primer</label>
    <div class="row">
      <span>regions <input type="number" id="nReg" min="0" max="20"></span>
      <span>threats <input type="number" id="nThr" min="0" max="15"></span>
      <span>seeds <input type="number" id="nSeed" min="0" max="15"></span>
    </div></div>
  <div><label>View</label>
    <div class="row">
      <input type="text" id="goto" placeholder="position or a-b">
      <button id="goBtn">Go</button>
      <button id="resetBtn">Reset</button>
      <button id="themeBtn">Dark</button>
    </div></div>
  <div><label>Export</label>
    <div class="row"><button id="svgBtn">Save SVG</button><button id="csvBtn">Save sites CSV</button></div></div>
</section>
<main>
  <div class="card" id="mapcard">
    <svg id="map" xmlns="http://www.w3.org/2000/svg" viewBox="-760 -620 1520 1240"></svg>
    <div id="legend"></div>
  </div>
  <div class="card" id="side">
    <div class="tabs" id="tabs">
      <button data-t="sel" class="on">Selected</button><button data-t="reg">Regions</button>
      <button data-t="thr">Threats</button><button data-t="seed">Seeds</button>
      <button data-t="cond">Conditions</button><button data-t="note">Notes</button>
    </div>
    <div class="tabp on" id="p-sel"><p style="color:var(--muted)">Click an arc, triangle, circle, diamond, label or table row to see its alignment here.</p></div>
    <div class="tabp" id="p-reg"></div>
    <div class="tabp" id="p-thr"></div>
    <div class="tabp" id="p-seed"></div>
    <div class="tabp" id="p-cond"></div>
    <div class="tabp" id="p-note"></div>
  </div>
</main>
<div id="tip"></div>
<script id="payload" type="application/json">__PAYLOAD__</script>
<script>
'use strict';
const D = JSON.parse(document.getElementById('payload').textContent);
const M = D.meta;
const L = M.L;
const NS = 'http://www.w3.org/2000/svg';
const SC = 300;
const NLEV = 256;
const MINUS = '\u2212', PRIME = '\u2032', EN = '\u2013', DEGC = '\u00b0C';
const LANE = {'0+':[0.940, 1.000], '1+':[0.875, 0.935], '0-':[0.740, 0.800], '1-':[0.675, 0.735]};
const R_THR = {'+':0.859, '-':0.816};
const R_AMP = 0.652;
const FS_LAB = 12;
const LAB_R = 1.235;
const TAGS = [D.primers[0].name, D.primers[1].name];
const $ = (id) => document.getElementById(id);

// Colour maps, the same ramps as the autumn and winter maps of the desktop figure
function mapAutumn(n) {
  const a = [];
  for (let i = 0; i < n; i++) a.push('rgb(255,' + Math.round(255*i/(n - 1)) + ',0)');
  return a;
}
function mapWinter(n) {
  const a = [];
  for (let i = 0; i < n; i++) a.push('rgb(0,' + Math.round(255*i/(n - 1)) + ',' + Math.round(255*(1 - 0.5*i/(n - 1))) + ')');
  return a;
}
const MAPS = [mapAutumn(NLEV), mapWinter(NLEV)];

function svgEl(tag, attrs, parent) {
  const e = document.createElementNS(NS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(e);
  return e;
}
function ang(x) { return -Math.PI/2 + 2*Math.PI*x/L; }
function pol(r, x) { const a = ang(x); return [r*Math.cos(a), r*Math.sin(a)]; }
function f1(p) { return p[0].toFixed(2) + ' ' + p[1].toFixed(2); }
function arcD(r, x0, x1) {
  const large = ((x1 - x0)/L > 0.5)?1:0;
  return 'M' + f1(pol(r, x0)) + 'A' + r.toFixed(2) + ' ' + r.toFixed(2) + ' 0 ' + large + ' 1 ' + f1(pol(r, x1));
}
function secD(ri, ro, x0, x1) {
  const large = ((x1 - x0)/L > 0.5)?1:0;
  return 'M' + f1(pol(ro, x0)) + 'A' + ro.toFixed(2) + ' ' + ro.toFixed(2) + ' 0 ' + large + ' 1 ' + f1(pol(ro, x1))
    + 'L' + f1(pol(ri, x1)) + 'A' + ri.toFixed(2) + ' ' + ri.toFixed(2) + ' 0 ' + large + ' 0 ' + f1(pol(ri, x0)) + 'Z';
}
function fmtSigned(x, d) { const s = x.toFixed(d); return (s[0] === '-')?MINUS + s.slice(1):s; }
function rangeTxt(a, b) { return a + EN + b; }
function esc(s) { return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;'); }
function num(x, d) { return (x === null || x === undefined || !isFinite(x))?'n/a':x.toFixed(d); }
function sci(x) { return (x === null || !isFinite(x))?'n/a':((x !== 0 && (Math.abs(x) < 1e-3 || Math.abs(x) >= 1e4))?x.toExponential(2):x.toPrecision(3)); }
function mod(a, m) { return ((a%m) + m)%m; }
function pn(p) { return TAGS[p]; }
function sdn(d) { return (d === '+')?'(+)':'(' + MINUS + ')'; }

// Per-site quantities that do not depend on the display threshold
const SITES = D.sites;
const BYP = [[], []];
SITES.forEach((s, k) => {
  s.k = k;
  BYP[s.p].push(k);
  const P = D.primers[s.p];
  const K = Math.exp(-s.g/P.RT);
  s.th = 1/(1 + 1/(K*P.f));
  if (s.an) {
    const Ka = Math.exp(-s.ga/P.RT);
    s.thA = 1/(1 + 1/(Ka*P.f));
    s.rate = s.thA*s.eps;
  } else { s.thA = 0; s.rate = 0; }
  s.sd3 = (s.e3 >= M.seed3Len) && !s.cg;
  s.sdI = (s.ir >= M.seedIntLen) && !s.cg;
  s.sdRank = s.sdI?s.ir:0;
  if (s.sd3) s.sdRank = 1000 + s.e3;
});

const ST = {thetaPct:M.thetaPct, show:[true, true], strand:{'+':true, '-':true}, threats:true,
  nReg:M.labels.reg, nThr:M.labels.thr, nSeed:M.labels.seed, sel:null, k:1, tx:0, ty:0, dark:false};
let CUR = null;
let VB = [-760, -620, 1520, 1240];

// Display model for the current threshold
function regionTable(best, bidx) {
  const regs = [];
  let k0 = -1;
  for (let i = 0; i < L; i++) { if (!isFinite(best[i])) { k0 = i; break; } }
  if (k0 < 0 && !isFinite(best[0])) return regs;
  const order = [];
  for (let i = 0; i < L; i++) order.push(mod(Math.max(k0, 0) + i, L));
  let i = 0;
  while (i < L) {
    if (isFinite(best[order[i]])) {
      let j = i;
      while (j + 1 < L && isFinite(best[order[j + 1]])) j++;
      let gmin = Infinity, kmin = -1;
      for (let q = i; q <= j; q++) {
        const v = best[order[q]];
        if (v < gmin) { gmin = v; kmin = bidx[order[q]]; }
      }
      const s = order[i] + 1, e = order[j] + 1;
      const len = j - i + 1;
      regs.push({s:s, e:e, len:len, g:gmin, k:kmin, x:mod(s - 1 + len/2, L)});
      i = j + 1;
    } else { i++; }
  }
  return regs;
}

function compute() {
  const thetaMin = ST.thetaPct/100;
  const V = {thetaMin:thetaMin, prim:[], gHi:-Infinity, gLo:Infinity, anyShown:false};
  for (let p = 0; p < 2; p++) {
    const P = D.primers[p];
    const cut = -P.RT*Math.log(thetaMin/((1 - thetaMin)*P.f));
    const best = {'+':new Float64Array(L).fill(Infinity), '-':new Float64Array(L).fill(Infinity)};
    const bidx = {'+':new Int32Array(L).fill(-1), '-':new Int32Array(L).fill(-1)};
    let nShown = 0, minG = Infinity;
    for (const k of BYP[p]) {
      const s = SITES[k];
      if (s.th < thetaMin) continue;
      nShown++;
      if (s.g < minG) minG = s.g;
      const b = best[s.d], bi = bidx[s.d];
      for (let t = 0; t < s.n; t++) {
        const pos = (s.a - 1 + t)%L;
        if (s.g < b[pos]) { b[pos] = s.g; bi[pos] = k; }
      }
    }
    const regs = {'+':regionTable(best['+'], bidx['+']), '-':regionTable(best['-'], bidx['-'])};
    const thrAll = BYP[p].filter((k) => SITES[k].an && !SITES[k].cg && SITES[k].rate >= thetaMin*M.threatFactor);
    thrAll.sort((a, b) => SITES[b].rate - SITES[a].rate);
    const seedAll = BYP[p].filter((k) => SITES[k].sd3 || SITES[k].sdI);
    seedAll.sort((a, b) => SITES[b].sdRank - SITES[a].sdRank);
    V.prim.push({cut:cut, best:best, bidx:bidx, regs:regs, nShown:nShown, minG:minG,
      thr:thrAll.slice(0, M.threatMax), nThr:thrAll.length, seeds:seedAll.slice(0, M.seedMax), nSeed:seedAll.length});
    if (nShown > 0) { V.anyShown = true; V.gLo = Math.min(V.gLo, minG); }
    V.gHi = Math.max(V.gHi, cut);
  }
  if (V.anyShown) V.gLo = Math.floor(2*V.gLo)/2; else V.gLo = V.gHi - 5;
  if (V.gHi - V.gLo < 1) V.gLo = V.gHi - 1;
  return V;
}
function cIdx(V, g) {
  const i = 1 + Math.round((g - V.gLo)/(V.gHi - V.gLo)*(NLEV - 1));
  return Math.min(Math.max(i, 1), NLEV) - 1;
}
function colour(V, p, g) { return MAPS[p][cIdx(V, g)]; }

// Labels
function siteLabel(pre, tS, tE, pf, pt, g, tm) {
  let t1 = pre + rangeTxt(tS, tE);
  if (pf !== null && pt !== null && pf !== undefined) t1 += ' | ' + rangeTxt(pf, pt);
  const lines = [t1, '\u0394G\u2095\u1d67\u2096 = ' + fmtSigned(g, 2)];
  lines[1] = '\u0394G hyb = ' + fmtSigned(g, 2);
  if (tm !== null && tm !== undefined && isFinite(tm)) lines.push('Tm = ' + tm.toFixed(1) + ' ' + DEGC);
  return lines;
}
function buildLabels(V) {
  const labs = [];
  for (let p = 0; p < 2; p++) {
    if (!ST.show[p]) continue;
    for (const d of ['+', '-']) {
      if (!ST.strand[d]) continue;
      const regs = V.prim[p].regs[d].slice().sort((a, b) => a.g - b.g).slice(0, ST.nReg);
      for (const r of regs) {
        const s = SITES[r.k];
        labs.push({cls:'reg', k:r.k, p:p, d:d, cw:360*r.x/L, col:colour(V, p, r.g), shape:'tri', edge:'#404040',
          lines:siteLabel(pn(p) + ' ' + sdn(d) + ' ', r.s, r.e, s.pf, s.pt, r.g, s.tm)});
      }
    }
    if (ST.threats) {
      for (const k of V.prim[p].thr.slice(0, ST.nThr)) {
        const s = SITES[k];
        if (!ST.strand[s.d]) continue;
        labs.push({cls:'thr', k:k, p:p, d:s.d, cw:360*s.am/L, col:colour(V, p, s.gae), shape:'circle', edge:'#000',
          lines:siteLabel(pn(p) + ' ' + sdn(s.d) + ' ', s.as, s.ae, s.af, s.at, s.ga, s.tma)});
      }
    }
    for (const k of V.prim[p].seeds.slice(0, ST.nSeed)) {
      const s = SITES[k];
      if (!ST.strand[s.d]) continue;
      const pre = pn(p) + ' ' + sdn(s.d) + ' ';
      let l1, l2, cw, shape, kind;
      if (s.sd3 && s.sdI) { shape = 'star'; cw = 360*s.s3m/L; l1 = pre + rangeTxt(s.s3a, s.s3b); l2 = '3' + PRIME + ' run ' + s.e3 + ', int ' + s.ir; kind = 'seed3'; }
      else if (s.sd3) { shape = 'square'; cw = 360*s.s3m/L; l1 = pre + rangeTxt(s.s3a, s.s3b); l2 = '3' + PRIME + ' run ' + s.e3 + ' nt'; kind = 'seed3'; }
      else { shape = 'diamond'; cw = 360*s.im/L; l1 = pre + rangeTxt(s.isa, s.isb); l2 = 'int run ' + s.ir + ' nt'; kind = 'seedI'; }
      labs.push({cls:kind, k:k, p:p, d:s.d, cw:cw, col:colour(V, p, s.g), shape:shape, edge:'#000', lines:[l1, l2]});
    }
  }
  return labs;
}
function shapeEl(shape, x, y, col, edge, parent) {
  const a = {fill:col, stroke:edge, 'stroke-width':1.1};
  if (shape === 'tri') { a.points = [x, y - 6.5, x + 6, y + 5, x - 6, y + 5].join(','); return svgEl('polygon', a, parent); }
  if (shape === 'circle') { a.cx = x; a.cy = y; a.r = 5; return svgEl('circle', a, parent); }
  if (shape === 'square') { a.x = x - 5.4; a.y = y - 5.4; a.width = 10.8; a.height = 10.8; return svgEl('rect', a, parent); }
  if (shape === 'diamond') { a.points = [x, y - 7, x + 7, y, x, y + 7, x - 7, y].join(','); return svgEl('polygon', a, parent); }
  const pts = [];
  for (let i = 0; i < 10; i++) {
    const r = (i%2 === 0)?9.5:3.9, t = -Math.PI/2 + i*Math.PI/5;
    pts.push((x + r*Math.cos(t)).toFixed(2) + ',' + (y + r*Math.sin(t)).toFixed(2));
  }
  a.points = pts.join(' ');
  return svgEl('polygon', a, parent);
}
function spreadAngles(cw, minSep) {
  const n = cw.length;
  if (n < 2) return cw.slice();
  const idx = cw.map((v, i) => i).sort((i, j) => mod(cw[i], 360) - mod(cw[j], 360));
  const s = idx.map((i) => mod(cw[i], 360));
  const gaps = s.map((v, i) => (i < n - 1)?s[i + 1] - v:s[0] + 360 - s[n - 1]);
  let g = 0;
  for (let i = 1; i < n; i++) if (gaps[i] > gaps[g]) g = i;
  const order = [];
  for (let i = g + 1; i < n; i++) order.push(i);
  for (let i = 0; i <= g; i++) order.push(i);
  const t = order.map((i) => s[i]);
  for (let i = 0; i < n; i++) if (t[i] < t[0]) t[i] += 360;
  for (let i = 1; i < n; i++) t[i] = Math.max(t[i], t[i - 1] + minSep);
  const out = new Array(n);
  order.forEach((si, q) => { out[idx[si]] = t[q]; });
  return out;
}
function sepDeg(nLines) {
  return nLines*FS_LAB*1.3/(LAB_R*SC*Math.PI/180);
}
function anyOverlap(bb) {
  for (let i = 0; i < bb.length - 1; i++) {
    for (let j = i + 1; j < bb.length; j++) {
      const a = bb[i], b = bb[j];
      if (a.x < b.x + b.width && b.x < a.x + a.width && a.y < b.y + b.height && b.y < a.y + a.height) return true;
    }
  }
  return false;
}
function placeLabel(e, cwTrue, cwLab) {
  const a0 = cwTrue*Math.PI/180, a1 = cwLab*Math.PI/180;
  const p0 = [1.05*SC*Math.sin(a0), -1.05*SC*Math.cos(a0)], p1 = [1.12*SC*Math.sin(a0), -1.12*SC*Math.cos(a0)];
  const p2 = [1.18*SC*Math.sin(a1), -1.18*SC*Math.cos(a1)];
  e.lead.setAttribute('points', [p0, p1, p2].map((q) => q[0].toFixed(1) + ',' + q[1].toFixed(1)).join(' '));
  e.mark.remove();
  e.mark = shapeEl(e.lab.shape, p2[0], p2[1], e.lab.col, e.lab.edge, e.g);
  e.mark.setAttribute('data-k', e.lab.k); e.mark.setAttribute('data-kind', e.lab.cls);
  placeText(e, cwLab);
}
function placeText(e, cwLab) {
  const a1 = cwLab*Math.PI/180;
  const sn = Math.sin(a1), cs = Math.cos(a1);
  const tx = LAB_R*SC*sn, ty = -LAB_R*SC*cs;
  const c = sn;
  const anchor = (c > 0.25)?'start':((c < -0.25)?'end':'middle');
  const lh = FS_LAB*1.25, nl = e.lab.lines.length;
  let y0;
  if (anchor === 'middle') y0 = (cs > 0)?ty - (nl - 1)*lh - 2:ty + FS_LAB;
  else y0 = ty - (nl - 1)*lh/2 + FS_LAB*0.35;
  e.text.setAttribute('text-anchor', anchor);
  e.text.setAttribute('x', tx.toFixed(1)); e.text.setAttribute('y', y0.toFixed(1));
  Array.from(e.text.children).forEach((ts, i) => { ts.setAttribute('x', tx.toFixed(1)); ts.setAttribute('dy', i === 0?0:lh); });
}
function layoutLabels(g, labs) {
  const els = [];
  if (!labs.length) return;
  const cwTrue = labs.map((l) => l.cw);
  const maxLines = Math.max(...labs.map((l) => l.lines.length));
  let minSep = Math.max(8, sepDeg(maxLines));
  let cwLab = spreadAngles(cwTrue, minSep);
  labs.forEach((lab) => {
    const eg = svgEl('g', {'data-k':lab.k, 'data-kind':lab.cls, class:'lab'}, g);
    const lead = svgEl('polyline', {fill:'none', stroke:'var(--tick)', 'stroke-width':1, 'vector-effect':'non-scaling-stroke'}, eg);
    const mark = svgEl('circle', {r:1}, eg);
    const text = svgEl('text', {'font-size':FS_LAB, class:'labtxt'}, eg);
    lab.lines.forEach((ln) => { const ts = svgEl('tspan', {}, text); ts.textContent = ln; });
    text.setAttribute('data-k', lab.k); text.setAttribute('data-kind', lab.cls);
    els.push({g:eg, lead:lead, mark:mark, text:text, lab:lab});
  });
  els.forEach((e, i) => placeText(e, cwLab[i]));

  // Relaxation: each overlapping pair is nudged apart by under a degree per pass, so crowded sides open up
  // where the text is wide and stay tight where it is not
  const n = els.length;
  for (let pass = 0; pass < 240; pass++) {
    const bb = els.map((e) => { const b = e.text.getBBox(); return {x:b.x - 2, y:b.y - 1, width:b.width + 4, height:b.height + 2}; });
    let moved = false;
    for (let i = 0; i < n - 1; i++) {
      for (let j = i + 1; j < n; j++) {
        const a = bb[i], b = bb[j];
        if (a.x < b.x + b.width && b.x < a.x + a.width && a.y < b.y + b.height && b.y < a.y + a.height) {
          const d = mod(cwLab[j] - cwLab[i] + 180, 360) - 180;
          const s = (d >= 0)?1:-1;
          cwLab[j] += s*0.8; cwLab[i] -= s*0.8;
          moved = true;
        }
      }
    }
    if (!moved) break;
    els.forEach((e, i) => placeText(e, cwLab[i]));
  }
  els.forEach((e, i) => placeLabel(e, cwTrue[i], cwLab[i]));
}

function niceStep(x) {
  const p = Math.pow(10, Math.floor(Math.log10(x)));
  for (const m of [1, 2, 2.5, 5, 10]) { if (m*p >= x) return m*p; }
  return 10*p;
}

// Drawing
const SVG_CSS = '.unbound{stroke:var(--unbound)}.outline{stroke:var(--outline);fill:none}'
  + '.labtxt{fill:var(--fg);font-family:system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}'
  + '.cap{fill:var(--muted);font-family:system-ui,Arial,sans-serif;paint-order:stroke;stroke:var(--bg);stroke-width:3px}'
  + '.ctr{fill:var(--fg);font-family:system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}'
  + '.mono{font-family:ui-monospace,Menlo,Consolas,monospace}.tick{stroke:var(--tick);fill:var(--tick)}'
  + '.lab{cursor:pointer}[data-k]{cursor:pointer}';

function render() {
  const V = compute();
  CUR = V;
  const svg = $('map');
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  svg.setAttribute('viewBox', VB.join(' '));
  const st = svgEl('style', {}, svg); st.textContent = SVG_CSS;
  const vp = svgEl('g', {id:'vp'}, svg);
  const gRing = svgEl('g', {}, vp), gArc = svgEl('g', {}, vp), gTick = svgEl('g', {}, vp);
  const gMark = svgEl('g', {}, vp), gThr = svgEl('g', {}, vp), gAmp = svgEl('g', {}, vp);
  const gHL = svgEl('g', {id:'hl'}, vp), gLab = svgEl('g', {}, vp), gCtr = svgEl('g', {}, vp);

  // Lanes
  for (let p = 0; p < 2; p++) {
    for (const d of ['+', '-']) {
      const lane = LANE[p + d];
      const rm = (lane[0] + lane[1])/2*SC, wd = (lane[1] - lane[0])*SC;
      svgEl('circle', {r:rm, fill:'none', 'stroke-width':wd, class:'unbound'}, gRing);
      svgEl('circle', {r:lane[0]*SC, class:'outline', 'stroke-width':0.4, 'vector-effect':'non-scaling-stroke'}, gRing);
      svgEl('circle', {r:lane[1]*SC, class:'outline', 'stroke-width':0.4, 'vector-effect':'non-scaling-stroke'}, gRing);
      if (!ST.show[p] || !ST.strand[d]) continue;
      const bi = V.prim[p].bidx[d], bs = V.prim[p].best[d];
      let i = 0;
      while (i < L) {
        if (bi[i] < 0) { i++; continue; }
        let j = i;
        while (j + 1 < L && bi[j + 1] === bi[i]) j++;
        const path = svgEl('path', {d:arcD(rm, i - 0.1, j + 1.1), fill:'none', stroke:colour(V, p, bs[i]), 'stroke-width':wd,
          'data-k':bi[i], 'data-kind':'reg'}, gArc);
        i = j + 1;
      }
    }
  }

  // Ticks on the (+) numbering
  const step = niceStep(L/12);
  const ticks = [1];
  for (let t = step; t <= L; t += step) ticks.push(t);
  ticks.forEach((t) => {
    if (t !== 1 && L - t < 0.2*step) return;
    const x = t - 0.5;
    const a = pol(1.000*SC, x), b = pol(1.028*SC, x), c = pol(1.085*SC, x);
    svgEl('line', {x1:a[0], y1:a[1], x2:b[0], y2:b[1], class:'tick', 'stroke-width':1, 'vector-effect':'non-scaling-stroke'}, gTick);
    const tx = svgEl('text', {x:c[0], y:c[1], 'text-anchor':'middle', 'dominant-baseline':'middle', 'font-size':11, class:'labtxt', opacity:0.7}, gTick);
    tx.textContent = t;
  });

  // Lane captions beside 12 o'clock
  for (let p = 0; p < 2; p++) {
    for (const d of ['+', '-']) {
      const lane = LANE[p + d];
      const rm = (lane[0] + lane[1])/2*SC;
      const arrow = (d === '+')?'5' + PRIME + '\u2192' + '3' + PRIME:'3' + PRIME + '\u2190' + '5' + PRIME;
      const tx = svgEl('text', {x:-6, y:-rm + 3, 'text-anchor':'end', 'font-size':8.5, class:'cap', 'pointer-events':'none'}, gTick);
      tx.textContent = pn(p) + ' ' + sdn(d) + ' ' + arrow;
    }
  }

  // Region triangles, outside the (+) lanes and inside the (-) lanes
  for (let p = 0; p < 2; p++) {
    if (!ST.show[p]) continue;
    for (const d of ['+', '-']) {
      if (!ST.strand[d]) continue;
      let regs = V.prim[p].regs[d];
      if (regs.length > M.maxMarkers) regs = regs.slice().sort((a, b) => a.g - b.g).slice(0, ST.nReg);
      const rTip = (d === '+')?1.005*SC:0.670*SC, rBase = (d === '+')?1.045*SC:0.630*SC;
      for (const r of regs) {
        const a = ang(r.x), da = 5/rBase;
        const pts = [[rTip*Math.cos(a), rTip*Math.sin(a)], [rBase*Math.cos(a - da), rBase*Math.sin(a - da)], [rBase*Math.cos(a + da), rBase*Math.sin(a + da)]];
        svgEl('polygon', {points:pts.map((q) => q[0].toFixed(1) + ',' + q[1].toFixed(1)).join(' '), fill:colour(V, p, r.g),
          stroke:'#404040', 'stroke-width':0.5, 'data-k':r.k, 'data-kind':'reg'}, gMark);
      }
    }
  }

  // Extension threat circles in the gap between the strands
  if (ST.threats) {
    for (let p = 0; p < 2; p++) {
      if (!ST.show[p]) continue;
      for (const k of V.prim[p].thr) {
        const s = SITES[k];
        if (!ST.strand[s.d]) continue;
        const c = pol(R_THR[s.d]*SC, s.am);
        svgEl('circle', {cx:c[0], cy:c[1], r:4.6, fill:colour(V, p, s.gae), stroke:'#000', 'stroke-width':0.9,
          'data-k':k, 'data-kind':'thr'}, gThr);
      }
    }
  }

  // Amplicon arc inside the (-) lanes
  if (D.amp) {
    const x0 = D.amp.a - 1, x1 = x0 + D.amp.len;
    svgEl('path', {d:arcD(R_AMP*SC, x0, Math.min(x1, x0 + L - 0.01)), fill:'none', stroke:'var(--amp)', 'stroke-width':3.2,
      'data-kind':'amp'}, gAmp);
    for (const x of [x0, x1]) {
      const a = pol(0.625*SC, x), b = pol(0.680*SC, x);
      svgEl('line', {x1:a[0], y1:a[1], x2:b[0], y2:b[1], stroke:'var(--amp)', 'stroke-width':2, 'data-kind':'amp'}, gAmp);
    }
  }

  // Labels
  layoutLabels(gLab, buildLabels(V));

  // Centre block
  const lines = [];
  lines.push({t:M.name, b:1, s:21});
  lines.push({t:L + ' nt, circular ' + M.topology, s:12.5});
  lines.push({t:'', s:6});
  for (let p = 0; p < 2; p++) {
    const P = D.primers[p];
    lines.push({t:P.name + ' primer (' + P.n + ' nt, ' + P.conc_nM + ' nM)', b:1, s:12.5});
    lines.push({t:'5' + PRIME + '-' + P.seq + '-3' + PRIME, mono:1, s:P.n > 26?10:11.5});
  }
  lines.push({t:'', s:6});
  lines.push({t:'Annealing ' + M.Ta + ' ' + DEGC + ', ' + M.monv_mM + ' mM monovalent, ' + M.biv_mM + ' mM bivalent', s:11.5});
  lines.push({t:'Total dNTP ' + M.dNTP_uM + ' \u00b5M, template ' + M.templ_nM + ' nM', s:11.5});
  if (D.amp) lines.push({t:'Amplicon ' + D.amp.len + ' nt (' + D.amp.a + EN + D.amp.b + ')', s:11.5});
  lines.push({t:'', s:6});
  for (let p = 0; p < 2; p++) {
    const P = D.primers[p];
    let best = null;
    for (const k of BYP[p]) { if (best === null || SITES[k].g < SITES[best].g) best = k; }
    if (best !== null) {
      const s = SITES[best];
      lines.push({t:P.name + ': ' + sdn(s.d) + ' ' + s.a + EN + s.b + ', \u0394G ' + fmtSigned(s.g, 2) + ', Tm ' + num(s.tm, 1) + ' ' + DEGC, s:11});
    }
  }
  lines.push({t:'Shown: \u03b8 \u2265 ' + ST.thetaPct.toPrecision(3).replace(/\.?0+$/, '') + '% (\u0394G \u2264 ' + fmtSigned(V.prim[0].cut, 2) + ' / ' + fmtSigned(V.prim[1].cut, 2) + ')', s:10.5});
  lines.push({t:'\u03b5 data: ' + M.effCite, s:10.5});
  const lh = (l) => l.s*1.32;
  const total = lines.reduce((a, l) => a + lh(l), 0);
  let y = -total/2;
  lines.forEach((l) => {
    y += lh(l);
    if (!l.t) return;
    const t = svgEl('text', {x:0, y:y - l.s*0.3, 'text-anchor':'middle', 'font-size':l.s, class:'ctr' + (l.mono?' mono':''), 'font-weight':l.b?'700':'400'}, gCtr);
    t.textContent = l.t;
  });

  // Fit the viewBox to the artwork
  const bb = vp.getBBox();
  const pad = 14;
  VB = [bb.x - pad, bb.y - pad, bb.width + 2*pad, bb.height + 2*pad];
  svg.setAttribute('viewBox', VB.join(' '));
  applyView();
  if (ST.sel) drawHL(ST.sel.k, ST.sel.kind, false);
  fillLegend(V);
  fillTables(V);
  updateHeader(V);
}

// Zoom and pan
function applyView() {
  const vp = $('vp');
  if (vp) vp.setAttribute('transform', 'translate(' + ST.tx + ' ' + ST.ty + ') scale(' + ST.k + ')');
}
function svgPoint(ev) {
  const svg = $('map');
  const pt = svg.createSVGPoint();
  pt.x = ev.clientX; pt.y = ev.clientY;
  return pt.matrixTransform(svg.getScreenCTM().inverse());
}
function zoomAt(px, py, factor) {
  const k2 = Math.min(Math.max(ST.k*factor, 1), 80);
  const f = k2/ST.k;
  ST.tx = px - (px - ST.tx)*f;
  ST.ty = py - (py - ST.ty)*f;
  ST.k = k2;
  if (ST.k === 1) { ST.tx = 0; ST.ty = 0; }
  applyView();
}
function focusX(x, r, k) {
  const c = pol(r*SC, x);
  const cx = VB[0] + VB[2]/2, cy = VB[1] + VB[3]/2;
  ST.k = k;
  ST.tx = cx - k*c[0];
  ST.ty = cy - k*c[1];
  applyView();
}
function resetView() { ST.k = 1; ST.tx = 0; ST.ty = 0; applyView(); }

// Site details
function siteKindRange(s, kind) {
  if (kind === 'thr') return [s.as, s.al];
  if (kind === 'seed3') return [s.s3a, s.e3];
  if (kind === 'seedI') return [s.isa, s.ir];
  return [s.a, s.n];
}
function drawHL(k, kind, zoom) {
  const g = $('hl');
  if (!g) return;
  while (g.firstChild) g.removeChild(g.firstChild);
  const s = SITES[k];
  const lane = LANE[s.p + s.d];
  const r = siteKindRange(s, kind);
  const x0 = r[0] - 1;
  svgEl('path', {d:secD(lane[0]*SC - 1.5, lane[1]*SC + 1.5, x0, x0 + Math.max(r[1], 1)), fill:'none', stroke:'var(--fg)', 'stroke-width':1.6,
    'vector-effect':'non-scaling-stroke', 'pointer-events':'none'}, g);
  if (zoom) focusX(x0 + r[1]/2, (lane[0] + lane[1])/2, Math.max(ST.k, 8));
}
function alnBlock(seg, pr, op) {
  return '<pre class="aln">5' + PRIME + '-' + esc(seg) + '-3' + PRIME + '  primer\n   ' + esc(pr) + '\n3' + PRIME + '-' + esc(op) + '-5' + PRIME + '  template</pre>';
}
function flagsOf(s) {
  const f = [];
  const thetaMin = ST.thetaPct/100;
  if (s.cg) f.push('cognate');
  if (s.th >= thetaMin) f.push('displayed');
  if (s.an && !s.cg && s.rate >= thetaMin*M.threatFactor) f.push('extension threat');
  if (s.sd3) f.push("3' seed");
  if (s.sdI) f.push('internal seed');
  return f.join(', ') || 'none';
}
function detailHtml(k, kind, full) {
  const s = SITES[k];
  const P = D.primers[s.p];
  let h = '<div class="kv">';
  h += '<span>Primer</span><span><b>' + esc(P.name) + '</b> ' + sdn(s.d) + ' strand</span>';
  h += '<span>Template</span><span>' + rangeTxt(s.a, s.b) + ' (' + s.n + ' nt), primer ' + rangeTxt(s.pf, s.pt) + ' of ' + P.n + '</span>';
  h += '<span>&Delta;G hyb</span><span>' + fmtSigned(s.g, 2) + ' kcal/mol, Tm ' + num(s.tm, 1) + ' ' + DEGC + '</span>';
  h += '<span>Occupancy</span><span>&theta; = ' + sci(s.th) + (s.tp?', 3&prime; terminal base paired':', 3&prime; terminal base mispaired') + '</span>';
  if (s.an) {
    h += '<span>3&prime;-reachable</span><span>' + rangeTxt(s.as, s.ae) + ', primer ' + rangeTxt(s.af, s.at) + ', &Delta;G ' + fmtSigned(s.ga, 2)
      + ', Tm ' + num(s.tma, 1) + ' ' + DEGC + '</span>';
    h += '<span>Extension</span><span>&epsilon; = ' + sci(s.eps) + ' (' + esc(s.ep) + '), &Delta;G eff ' + fmtSigned(s.gae, 2) + ', &theta;&middot;&epsilon; = ' + sci(s.rate) + '</span>';
  } else {
    h += '<span>3&prime;-reachable</span><span>none within the end tolerance</span>';
  }
  h += '<span>Seeds</span><span>3&prime; run ' + s.e3 + ' nt, internal run ' + s.ir + ' nt</span>';
  h += '<span>Flags</span><span>' + flagsOf(s) + '</span></div>';
  if (kind === 'thr' && s.an) h += alnBlock(s.aseg, s.apr, s.aop);
  else h += alnBlock(s.seg, s.pr, s.op);
  if (full && kind !== 'thr' && s.an && (s.af !== s.pf || s.at !== s.pt)) h += '<div style="font-size:.76rem;color:var(--muted)">3&prime;-reachable duplex</div>' + alnBlock(s.aseg, s.apr, s.aop);
  return h;
}
function selectSite(k, kind, zoom) {
  ST.sel = {k:k, kind:kind};
  $('p-sel').innerHTML = detailHtml(k, kind, true);
  showTab('sel');
  drawHL(k, kind, zoom);
}
function showTab(t) {
  document.querySelectorAll('#tabs button').forEach((b) => b.classList.toggle('on', b.dataset.t === t));
  document.querySelectorAll('.tabp').forEach((p) => p.classList.toggle('on', p.id === 'p-' + t));
}

// Legend, tables and header
function fillLegend(V) {
  const lg = $('legend');
  const stops = (p) => { const a = []; for (let i = 0; i <= 12; i++) a.push(MAPS[p][Math.round(i/12*(NLEV - 1))]); return a.join(','); };
  const ticks = [];
  for (let i = 0; i <= 4; i++) ticks.push(fmtSigned(V.gLo + (V.gHi - V.gLo)*i/4, 1));
  let h = '';
  for (let p = 0; p < 2; p++) {
    h += '<div><div>' + esc(pn(p)) + ' primer, &Delta;G hyb (kcal/mol)</div><div class="bar" style="background:linear-gradient(to right,' + stops(p) + ')"></div>'
      + '<div class="barticks">' + ticks.map((t) => '<span>' + t + '</span>').join('') + '</div></div>';
  }
  const sw = (svgInner) => '<svg width="16" height="16" viewBox="-8 -8 16 16">' + svgInner + '</svg>';
  h += '<div style="display:grid;gap:3px">'
    + '<span class="key">' + sw('<rect x="-7" y="-4" width="14" height="8" fill="var(--unbound)" stroke="var(--outline)" stroke-width=".5"/>') + 'No displayed duplex</span>'
    + '<span class="key">' + sw('<polygon points="0,-6 6,5 -6,5" fill="#aaa" stroke="#404040"/>') + 'Displayed region</span>'
    + '<span class="key">' + sw('<circle r="5" fill="#aaa" stroke="#000"/>') + 'Extendable off-target site</span>'
    + '<span class="key">' + sw('<rect x="-5" y="-5" width="10" height="10" fill="#aaa" stroke="#000"/>') + 'Seed: square 3&prime; run, diamond internal, star both</span>'
    + (D.amp?'<span class="key">' + sw('<line x1="-7" y1="0" x2="7" y2="0" stroke="var(--amp)" stroke-width="3"/>') + 'Amplicon</span>':'')
    + '</div>';
  lg.innerHTML = h;
}
function thead(cols) { return '<thead><tr>' + cols.map((c) => '<th>' + c + '</th>').join('') + '</tr></thead>'; }
function fillTables(V) {
  // Regions
  let rows = [];
  for (let p = 0; p < 2; p++) for (const d of ['+', '-']) for (const r of V.prim[p].regs[d]) rows.push({p:p, d:d, r:r});
  rows.sort((a, b) => a.r.g - b.r.g);
  let h = '<div class="tblwrap"><table>' + thead(['Primer', 'Strand', 'Template', 'Len', '&Delta;G', 'Tm', 'Primer']);
  rows.forEach((o) => {
    const s = SITES[o.r.k];
    h += '<tr class="r" data-k="' + o.r.k + '" data-kind="reg"><td><span class="dot" style="background:' + colour(V, o.p, o.r.g) + '"></span>' + esc(pn(o.p)) + '</td><td>' + sdn(o.d) + '</td><td>'
      + rangeTxt(o.r.s, o.r.e) + '</td><td>' + o.r.len + '</td><td>' + fmtSigned(o.r.g, 2) + '</td><td>' + num(s.tm, 1) + '</td><td>' + rangeTxt(s.pf, s.pt) + '</td></tr>';
  });
  $('p-reg').innerHTML = (rows.length?h + '</table></div>':'<p>No region at this threshold.</p>');

  // Threats
  rows = [];
  for (let p = 0; p < 2; p++) for (const k of V.prim[p].thr) rows.push(k);
  rows.sort((a, b) => SITES[b].rate - SITES[a].rate);
  h = '<div class="tblwrap"><table>' + thead(['Primer', 'Strand', 'Template', 'Primer', '&Delta;G', 'Tm', '&epsilon;', '&theta;&middot;&epsilon;', "3&prime;"]);
  rows.forEach((k) => {
    const s = SITES[k];
    h += '<tr class="r" data-k="' + k + '" data-kind="thr"><td><span class="dot" style="background:' + colour(V, s.p, s.gae) + '"></span>' + esc(pn(s.p)) + '</td><td>' + sdn(s.d) + '</td><td>'
      + rangeTxt(s.as, s.ae) + '</td><td>' + rangeTxt(s.af, s.at) + '</td><td>' + fmtSigned(s.ga, 2) + '</td><td>' + num(s.tma, 1) + '</td><td>' + sci(s.eps) + '</td><td>' + sci(s.rate) + '</td><td>' + (s.tpa?'paired':'mispaired') + '</td></tr>';
  });
  const nAll = V.prim[0].nThr + V.prim[1].nThr;
  $('p-thr').innerHTML = (rows.length?h + '</table></div>':'<p>No extension threat at this threshold.</p>')
    + '<p style="font-size:.76rem;color:var(--muted)">' + nAll + ' site(s) clear &theta;&middot;&epsilon; &ge; ' + M.threatFactor + '&times;&theta;min; the longest list shown per primer is ' + M.threatMax + '. The cognate site is excluded.</p>';

  // Seeds
  rows = [];
  for (let p = 0; p < 2; p++) for (const k of V.prim[p].seeds) rows.push(k);
  h = '<div class="tblwrap"><table>' + thead(['Primer', 'Strand', 'Class', 'Template', 'Primer', 'Run', 'Site &Delta;G', 'Tm']);
  rows.forEach((k) => {
    const s = SITES[k];
    const cls = (s.sd3 && s.sdI)?'both':(s.sd3?"3&prime; run":'internal');
    const kind = s.sd3?'seed3':'seedI';
    const tmpl = s.sd3?rangeTxt(s.s3a, s.s3b):rangeTxt(s.isa, s.isb);
    const pr = s.sd3?rangeTxt(P_N(s.p) - s.e3 + 1, P_N(s.p)):rangeTxt(s.ifr, s.ito);
    h += '<tr class="r" data-k="' + k + '" data-kind="' + kind + '"><td><span class="dot" style="background:' + colour(V, s.p, s.g) + '"></span>' + esc(pn(s.p)) + '</td><td>' + sdn(s.d) + '</td><td>' + cls + '</td><td>'
      + tmpl + '</td><td>' + pr + '</td><td>' + (s.sd3?s.e3:s.ir) + ' nt</td><td>' + fmtSigned(s.g, 2) + '</td><td>' + num(s.tm, 1) + '</td></tr>';
  });
  const nSeed = V.prim[0].nSeed + V.prim[1].nSeed;
  $('p-seed').innerHTML = (rows.length?h + '</table></div>':'<p>No seed site.</p>')
    + '<p style="font-size:.76rem;color:var(--muted)">' + nSeed + ' site(s): a 3&prime; run of at least ' + M.seed3Len + ' nt, or an internal run of at least ' + M.seedIntLen
    + ' nt. No free energy is quoted for a run alone; the site columns describe the register that carries it.</p>';
}
function P_N(p) { return D.primers[p].n; }
function updateHeader(V) {
  $('thrTxt').textContent = '\u03b8 \u2265 ' + ST.thetaPct.toPrecision(3).replace(/\.?0+$/, '') + ' %   (\u0394G \u2264 ' + fmtSigned(V.prim[0].cut, 2) + ' ' + TAGS[0] + ', '
    + fmtSigned(V.prim[1].cut, 2) + ' ' + TAGS[1] + ')';
}
function fillStatic() {
  $('ttl').textContent = M.name + ': forward and reverse primer hybridization';
  $('sub').textContent = L + ' nt circular ' + M.topology + '; ' + D.primers[0].name + ' ' + D.primers[0].n + ' nt and ' + D.primers[1].name + ' ' + D.primers[1].n
    + ' nt; annealing ' + M.Ta + ' ' + DEGC + '. Scroll to zoom, drag to pan, click an arc, circle, triangle or label for the alignment.';
  $('shFl').textContent = TAGS[0]; $('shRl').textContent = TAGS[1];
  const kv = (rows) => '<div class="kv">' + rows.map((r) => '<span>' + r[0] + '</span><span>' + r[1] + '</span>').join('') + '</div>';
  const pr = (p) => {
    const P = D.primers[p];
    const cog = P.cognate?sdn(P.cognate.d) + ' ' + rangeTxt(P.cognate.a, P.cognate.b):'none';
    return [P.name + ' primer', '<span class="mono">5' + PRIME + '-' + P.seq + '-3' + PRIME + '</span>'],
      [P.name + ' concentration', P.conc_nM + ' nM, free ' + P.f_nM.toPrecision(5) + ' nM, ' + (100*P.bound).toPrecision(3) + '% bound'],
      [P.name + ' cognate site', cog],
      [P.name + ' end tolerance', P.endTol + (M.proofread?' (whole primer)':' base(s) outside the duplex')];
  };
  const rows = [['Circular DNA', esc(M.name)], ['Length', L + ' nt'], ['Template', M.topology + ', ' + M.templ_nM + ' nM (' + M.templ_ngul + ' ng/\u00b5L)'],
    ['Annealing temperature', M.Ta + ' ' + DEGC], ['Monovalent cation', M.monv_mM + ' mM'],
    ['Bivalent cation', M.biv_mM + ' mM (' + M.freeBiv_mM + ' mM free)'], ['Total dNTP', M.dNTP_uM + ' \u00b5M']];
  for (let p = 0; p < 2; p++) {
    const P = D.primers[p];
    const cog = P.cognate?sdn(P.cognate.d) + ' ' + rangeTxt(P.cognate.a, P.cognate.b):'none';
    rows.push([P.name + ' primer', '<span class="mono">5' + PRIME + '-' + P.seq + '-3' + PRIME + '</span>']);
    rows.push([P.name + ' concentration', P.conc_nM + ' nM, free ' + P.f_nM.toPrecision(5) + ' nM, ' + (100*P.bound).toPrecision(3) + '% bound']);
    rows.push([P.name + ' cognate site', cog]);
    rows.push([P.name + ' end tolerance', P.endTol + (M.proofread?' (whole primer)':' base(s) outside the duplex')]);
  }
  rows.push(['Proofreading class', M.proofClass]);
  rows.push(['Extension efficiency', esc(M.effApplied) + '; ' + esc(M.effCite)]);
  rows.push(['Efficiency table', '<span class="mono" style="word-break:break-all">' + esc(M.effSource) + '</span>']);
  rows.push(['Survey run', M.surveyRun + ' bp']);
  rows.push(['Seed lengths', "3' " + M.seed3Len + ' nt, internal ' + M.seedIntLen + ' nt']);
  rows.push(['Symmetry correction', M.applyInit?'applied to self-complementary stretches':'not applied']);
  rows.push(['Sub-duplex edges', M.exact?'every Watson-Crick pair':'block boundaries only']);
  if (D.amp) rows.push(['Amplicon', D.amp.len + ' nt, ' + D.amp.a + EN + D.amp.b + ' clockwise']);
  $('p-cond').innerHTML = kv(rows);
  $('p-note').innerHTML = (D.notes.length?D.notes.map((n) => '<div class="warn">' + esc(n) + '</div>').join(''):'<p>No warnings.</p>')
    + '<p style="font-size:.8rem;color:var(--muted)">The map shows binding, not priming. A site can be occupied and still be catalytically dead when the 3&prime; terminus lies outside the duplex; the threat circles and the 3&prime; end column carry that distinction. '
    + 'For a single-stranded template the (' + MINUS + ') lanes show the occupancy a (' + MINUS + ') site would reach once a (' + MINUS + ') strand has been made, at the same free primer.</p>';
}

// Export
function exportSvg() {
  const svg = $('map').cloneNode(true);
  const cs = getComputedStyle(document.documentElement);
  const vars = ['--bg', '--fg', '--muted', '--unbound', '--outline', '--tick', '--amp'];
  svg.setAttribute('style', vars.map((v) => v + ':' + cs.getPropertyValue(v).trim()).join(';') + ';background:' + cs.getPropertyValue('--bg').trim());
  svg.setAttribute('width', 1600); svg.setAttribute('height', Math.round(1600*VB[3]/VB[2]));
  const blob = new Blob([new XMLSerializer().serializeToString(svg)], {type:'image/svg+xml'});
  download(blob, M.name + '_HybMap.svg');
}
function download(blob, name) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob); a.download = name;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 2000);
}
function exportCsv() {
  const thetaMin = ST.thetaPct/100;
  const cols = ['Primer', 'Strand', 'Start', 'End', 'Length', 'PrimerFrom', 'PrimerTo', 'dG_kcal_per_mol', 'Tm_site_C', 'Occupancy', 'Displayed',
    'PrimerThreePrimePaired', 'PrimerSegment_5to3', 'Pairing', 'Template_3to5', 'Anchored', 'dGanch_kcal_per_mol', 'Tm_anchored_C', 'Epsilon', 'EpsilonPair',
    'dGanchEff_kcal_per_mol', 'PrimingRate', 'AnchStart', 'AnchEnd', 'AnchPrimerFrom', 'AnchPrimerTo', 'Cognate', 'ExtensionThreat',
    'ThreePrimeRun', 'InternalRun', 'Seed3Prime', 'SeedInternal'];
  const rows = [cols.join(',')];
  const order = SITES.slice().sort((a, b) => a.g - b.g);
  order.forEach((s) => {
    const thr = s.an && !s.cg && s.rate >= thetaMin*M.threatFactor;
    const v = [pn(s.p), s.d, s.a, s.b, s.n, s.pf, s.pt, s.g, s.tm === null?'':s.tm, s.th, s.th >= thetaMin?1:0, s.tp, s.seg, s.pr, s.op, s.an,
      s.an?s.ga:'', s.an && s.tma !== null?s.tma:'', s.an?s.eps:'', s.an?s.ep:'', s.an?s.gae:'', s.an?s.rate:'', s.an?s.as:'', s.an?s.ae:'', s.an?s.af:'', s.an?s.at:'',
      s.cg, thr?1:0, s.e3, s.ir, s.sd3?1:0, s.sdI?1:0];
    rows.push(v.map((x) => (typeof x === 'string' && /[",]/.test(x))?'"' + x.replace(/"/g, '""') + '"':x).join(','));
  });
  download(new Blob([rows.join('\n')], {type:'text/csv'}), M.name + '_HybMap_Sites.csv');
}

// Tooltip and events
function tipHtml(k, kind) {
  return '<b>' + esc(D.primers[SITES[k].p].name) + '</b> ' + detailHtml(k, kind, false).replace(/class="kv"/, 'class="kv" style="color:inherit"');
}
function wire() {
  const svg = $('map'), tip = $('tip');
  let drag = null, moved = false;
  svg.addEventListener('mousemove', (ev) => {
    if (drag) return;
    const t = ev.target.closest('[data-k]');
    if (!t) { tip.style.display = 'none'; return; }
    const k = +t.getAttribute('data-k'), kind = t.getAttribute('data-kind');
    tip.innerHTML = tipHtml(k, kind);
    tip.style.display = 'block';
    const x = Math.min(ev.clientX + 14, window.innerWidth - tip.offsetWidth - 8);
    const y = Math.min(ev.clientY + 14, window.innerHeight - tip.offsetHeight - 8);
    tip.style.left = x + 'px'; tip.style.top = y + 'px';
  });
  svg.addEventListener('mouseleave', () => { tip.style.display = 'none'; });
  svg.addEventListener('wheel', (ev) => {
    ev.preventDefault();
    const p = svgPoint(ev);
    zoomAt(p.x, p.y, ev.deltaY < 0?1.25:0.8);
  }, {passive:false});
  svg.addEventListener('mousedown', (ev) => {
    drag = {x:ev.clientX, y:ev.clientY, tx:ST.tx, ty:ST.ty};
    moved = false;
  });
  window.addEventListener('mousemove', (ev) => {
    if (!drag) return;
    const dx = ev.clientX - drag.x, dy = ev.clientY - drag.y;
    if (Math.abs(dx) + Math.abs(dy) > 3) { moved = true; svg.classList.add('drag'); tip.style.display = 'none'; }
    if (!moved) return;
    const ctm = svg.getScreenCTM();
    ST.tx = drag.tx + dx/ctm.a;
    ST.ty = drag.ty + dy/ctm.d;
    applyView();
  });
  window.addEventListener('mouseup', () => { drag = null; svg.classList.remove('drag'); });
  svg.addEventListener('click', (ev) => {
    if (moved) { moved = false; return; }
    const t = ev.target.closest('[data-k]');
    if (!t) return;
    const kind = t.getAttribute('data-kind');
    if (kind === 'amp') return;
    selectSite(+t.getAttribute('data-k'), kind, false);
  });
  svg.addEventListener('dblclick', resetView);
  document.querySelectorAll('.tabp').forEach((p) => p.addEventListener('click', (ev) => {
    const tr = ev.target.closest('tr.r');
    if (!tr) return;
    document.querySelectorAll('tr.r.on').forEach((r) => r.classList.remove('on'));
    tr.classList.add('on');
    selectSite(+tr.getAttribute('data-k'), tr.getAttribute('data-kind'), true);
  }));
  document.querySelectorAll('#tabs button').forEach((b) => b.addEventListener('click', () => showTab(b.dataset.t)));
}
function sliderToPct(v) {
  const lo = Math.log10(M.thetaFloorPct), hi = Math.log10(50);
  return Math.pow(10, lo + (hi - lo)*v/1000);
}
function pctToSlider(p) {
  const lo = Math.log10(M.thetaFloorPct), hi = Math.log10(50);
  return Math.round(1000*(Math.log10(Math.max(p, M.thetaFloorPct)) - lo)/(hi - lo));
}
let pending = false;
function redraw() {
  if (pending) return;
  pending = true;
  requestAnimationFrame(() => { pending = false; render(); });
}
function init() {
  try {
    if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
      ST.dark = true;
      document.documentElement.setAttribute('data-theme', 'dark');
      $('themeBtn').textContent = 'Light';
    }
  } catch (err) { /* The page renders in the light theme when the query is unavailable */ }
  fillStatic();
  wire();
  const thr = $('thr');
  thr.value = pctToSlider(ST.thetaPct);
  thr.addEventListener('input', () => { ST.thetaPct = sliderToPct(+thr.value); ST.sel = null; $('hl') && ($('hl').innerHTML = ''); redraw(); });
  $('shF').addEventListener('change', (e) => { ST.show[0] = e.target.checked; redraw(); });
  $('shR').addEventListener('change', (e) => { ST.show[1] = e.target.checked; redraw(); });
  $('shP').addEventListener('change', (e) => { ST.strand['+'] = e.target.checked; redraw(); });
  $('shM').addEventListener('change', (e) => { ST.strand['-'] = e.target.checked; redraw(); });
  $('shT').addEventListener('change', (e) => { ST.threats = e.target.checked; redraw(); });
  [['nReg', 'nReg'], ['nThr', 'nThr'], ['nSeed', 'nSeed']].forEach((q) => {
    const el = $(q[0]);
    el.value = ST[q[1]];
    el.addEventListener('change', () => { ST[q[1]] = Math.max(0, Math.min(+el.value || 0, +el.max)); el.value = ST[q[1]]; redraw(); });
  });
  $('resetBtn').addEventListener('click', resetView);
  $('goBtn').addEventListener('click', () => {
    const m = $('goto').value.match(/(\d+)\s*[-\u2013:]?\s*(\d+)?/);
    if (!m) return;
    const a = +m[1], b = m[2]?+m[2]:a;
    const span = mod(b - a, L) + 1;
    focusX(a - 1 + span/2, 0.84, Math.min(8, Math.max(1.5, L/(span*30))));
  });
  $('themeBtn').addEventListener('click', () => {
    ST.dark = !ST.dark;
    document.documentElement.setAttribute('data-theme', ST.dark?'dark':'light');
    $('themeBtn').textContent = ST.dark?'Light':'Dark';
  });
  $('svgBtn').addEventListener('click', exportSvg);
  $('csvBtn').addEventListener('click', exportCsv);
  render();
}
init();
</script>
</body></html>
'''


# Main

def main():
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description = 'Interactive hybridization map of a forward and a reverse primer on a circular DNA.')
    ap.add_argument('--seq-file', type = str, default = None, help = 'Plain text, FASTA or GenBank file of the circular DNA')
    ap.add_argument('--name', type = str, default = None, help = 'Name shown for the circular DNA')
    ap.add_argument('--strands', type = int, default = None, choices = [1, 2], help = '1 = circular ssDNA, 2 = circular dsDNA (default 2)')
    ap.add_argument('--fwd', type = str, default = None, help = "Forward primer 5'->3', reads like the (+) strand")
    ap.add_argument('--rev', type = str, default = None, help = "Reverse primer 5'->3', the reverse complement of the (+) strand")
    ap.add_argument('--ta', type = float, default = None, help = 'Annealing temperature in C (default 60)')
    ap.add_argument('--monv', type = float, default = None, help = 'Monovalent cation in mM (default 40)')
    ap.add_argument('--biv', type = float, default = None, help = 'Bivalent cation in mM (default 2.5)')
    ap.add_argument('--dntp', type = float, default = None, help = 'Total dNTP in uM (default 800)')
    ap.add_argument('--fwd-conc-nM', type = float, default = None, help = 'Forward (excess) primer in nM (default 1000)')
    ap.add_argument('--rev-conc-nM', type = float, default = None, help = 'Reverse (limiting) primer in nM (default 20)')
    ap.add_argument('--template', type = float, default = None, help = 'Template concentration, value (default 1)')
    ap.add_argument('--template-unit', type = str, default = None, help = 'nM or ng/uL (default ng/uL)')
    ap.add_argument('--theta-pct', type = float, default = None, help = 'Display threshold, equilibrium occupancy in percent (default 0.01)')
    ap.add_argument('--proofread', type = str, default = None, help = 'Polymerase proofreading class: none, weak or strong (default none)')
    ap.add_argument('--out', type = str, default = None, help = 'Output HTML path')
    ap.add_argument('--params-dir', type = str, default = str(here), help = 'Folder holding the parameter CSV files (default: beside this script)')
    ap.add_argument('--eff-table', type = str, default = None, help = 'Extension-efficiency CSV (default: extension_efficiency_rejali2018.csv beside this script)')
    ap.add_argument('--survey-run', type = int, default = 3, help = 'Contiguous Watson-Crick pairs before a register is scored (default 3)')
    ap.add_argument('--end-tol', type = int, default = 1, help = "Primer 3' bases allowed outside the duplex for a Taq-type enzyme (default 1)")
    ap.add_argument('--no-eps', action = 'store_true', help = 'Score threats on occupancy alone (every epsilon is 1)')
    ap.add_argument('--threat-factor', type = float, default = 0.1, help = 'Threat floor as a multiple of the display threshold (default 0.1)')
    ap.add_argument('--threat-max', type = int, default = 15, help = 'Most threats listed per primer (default 15)')
    ap.add_argument('--seed3-len', type = int, default = 7, help = "Contiguous pairs ending at the primer 3' terminus (default 7)")
    ap.add_argument('--seed-int-len', type = int, default = 9, help = 'Longest internal run inside primer 1..n-1 (default 9)')
    ap.add_argument('--seed-max', type = int, default = 15, help = 'Most seed sites listed per primer (default 15)')
    ap.add_argument('--d-tm', type = float, default = 5.0, help = 'Degrees C between the two points used for dH and dS (default 5)')
    ap.add_argument('--theta-floor-pct', type = float, default = 1e-4, help = 'Lowest occupancy the threshold slider can reach, in percent (default 1e-4)')
    ap.add_argument('--region-labels', type = int, default = 3, help = 'Region labels per lane (default 3)')
    ap.add_argument('--threat-labels', type = int, default = 3, help = 'Threat labels per primer (default 3)')
    ap.add_argument('--seed-labels', type = int, default = 3, help = 'Seed labels per primer (default 3)')
    ap.add_argument('--no-exact-subduplex', action = 'store_true', help = 'Place sub-duplex edges at block boundaries only (faster, misses some optima)')
    ap.add_argument('--no-symmetry-correction', action = 'store_true', help = 'Do not apply the twofold symmetry term to self-complementary stretches')
    args = ap.parse_args()

    # Sequence
    seq_file = ask_sequence_file(args.seq_file)
    with open(seq_file, 'r', encoding = 'utf-8', errors = 'ignore') as fh:
        plus, parsed_name = parse_sequence_text(fh.read(), os.path.splitext(os.path.basename(seq_file))[0])
    bad = [i for i, c in enumerate(plus) if c not in 'ACGT']
    if bad:
        raise SystemExit(f'The circular DNA carries {len(bad)} character(s) outside {{A,C,G,T}}; the first sits at '
                         f'position {bad[0] + 1} and is {plus[bad[0]]}.')
    L = len(plus)
    if L < 20:
        raise SystemExit(f'The circular DNA is {L} nt long; at least 20 nt are needed.')

    # Conditions
    name = ask(args.name, 'Name of the circular DNA', str, parsed_name)
    n_strands = ask(args.strands, 'Template strands, 1 = circular ssDNA, 2 = circular dsDNA', int, 2)
    fwd = ask(args.fwd, "Forward primer, 5'->3'", str, DEFAULT_FWD)
    rev = ask(args.rev, "Reverse primer, 5'->3'", str, DEFAULT_REV)
    ta = ask(args.ta, 'Annealing temperature (C)', float, 60.0)
    monv_mM = ask(args.monv, 'Monovalent cation (mM)', float, 40.0)
    biv_mM = ask(args.biv, 'Bivalent cation (mM)', float, 2.5)
    dntp_uM = ask(args.dntp, 'Total dNTP (uM)', float, 800.0)
    fwd_nM = ask(args.fwd_conc_nM, 'Forward (excess) primer (nM)', float, 1000.0)
    rev_nM = ask(args.rev_conc_nM, 'Reverse (limiting) primer (nM)', float, 20.0)
    t_value = ask(args.template, 'Template concentration, value', float, 1.0)
    t_unit = ask(args.template_unit, 'Template concentration, unit (nM or ng/uL)', str, 'ng/uL')
    theta_pct = ask(args.theta_pct, 'Display threshold, equilibrium occupancy (%)', float, 0.01)
    proof_class = ask(args.proofread, 'Polymerase 3\'->5\' proofreading exonuclease (none / weak / strong)', str, 'none')
    if n_strands not in (1, 2):
        raise SystemExit('Template strands must be 1 or 2.')
    specs = []
    for idx, (nm, sq, conc) in enumerate([('Forward', fwd, fwd_nM), ('Reverse', rev, rev_nM)]):
        sq = re.sub('[^A-Za-z]', '', sq).upper().replace('U', 'T')
        if any(c not in 'ACGT' for c in sq):
            raise SystemExit(f'The {nm.lower()} primer carries characters outside {{A,C,G,T}}.')
        if (len(sq) < 6) or (len(sq) >= L):
            raise SystemExit(f'The {nm.lower()} primer must be 6 to {L - 1} nt long; it is {len(sq)} nt.')
        specs.append({'index':idx, 'name':nm, 'seq':sq, 'conc_nM':check_number(conc, f'{nm} primer concentration', 1e-12, 1e7)})
    ta = check_number(ta, 'Annealing temperature', 0, 100)
    monv_mM = check_number(monv_mM, 'Monovalent cation', 1e-12, 2000)
    biv_mM = check_number(biv_mM, 'Bivalent cation', 0, 1000)
    dntp_uM = check_number(dntp_uM, 'Total dNTP', 0, 1e5)
    t_value = check_number(t_value, 'Template concentration', 1e-12, 1e7)
    theta_pct = check_number(theta_pct, 'Display threshold', 1e-12, 100 - 1e-12)
    pc = proof_class.strip().lower()
    pc = {'y':'weak', 'yes':'weak', '1':'weak', 'true':'weak', 'n':'none', 'no':'none', '0':'none', 'false':'none'}.get(pc, pc)
    if pc not in ('none', 'weak', 'strong'):
        raise SystemExit(f"The proofreading answer must be none, weak or strong; received '{proof_class}'.")
    proofread = pc != 'none'

    opt = {'survey_run':args.survey_run, 'end_tol':args.end_tol, 'use_eps':(not args.no_eps) and (not proofread),
           'threat_factor':args.threat_factor, 'threat_max':args.threat_max, 'exact':not args.no_exact_subduplex,
           'seed3':args.seed3_len, 'seedint':args.seed_int_len, 'seed_max':args.seed_max, 'apply_init':not args.no_symmetry_correction,
           'proofread':proofread, 'ta':ta, 'd_tm':args.d_tm}
    notes = []
    if proofread:
        notes.append('A proofreading polymerase was declared, so two things changed. Every extension efficiency is set '
                     "to 1, and the end tolerance was widened from its Taq value to the whole primer, because the "
                     "exonuclease removes unpaired 3' overhanging nucleotides and primes from the outermost paired "
                     'base. Threats are ranked on occupancy of the trimmed duplex and no overhang length disqualifies '
                     'a site. Whether trimming outruns dissociation is kinetics and is not modelled, and nothing here '
                     'models the nibbling of free single-stranded primers either.')
    if pc == 'strong':
        notes.append('Class strong: Veneziano et al. report that Phusion was unable to synthesize fragments larger than '
                     '1 kb of ssDNA, which they attribute to enhanced exonuclease activity, and that Taq-based polymerases '
                     'gave the highest ssDNA production. Both are hypotheses offered by the authors. The warning concerns '
                     'yield, which nothing here models; the map is drawn exactly as for a weak exonuclease '
                     '(Sci Rep 2018;8:6548, doi:10.1038/s41598-018-24677-5).')

    # Parameters and derived quantities
    th = Thermo(args.params_dir)
    eff_path = Path(args.eff_table) if args.eff_table else Path(args.params_dir)/'extension_efficiency_rejali2018.csv'
    ee = ExtensionEfficiency(eff_path)
    monv, biv, dntp = monv_mM*1e-3, biv_mM*1e-3, dntp_uM*1e-6
    rx = Rxn(monv, biv, dntp, ta)
    rx2 = Rxn(monv, biv, dntp, ta + opt['d_tm'])
    res_mass = {'A':313.2, 'C':289.2, 'G':329.2, 'T':304.2}
    mw_plus = sum(res_mass[c] for c in plus)
    mw_minus = sum(res_mass[COMP[c]] for c in plus)
    mw = mw_plus if n_strands == 1 else mw_plus + mw_minus
    topo = 'ssDNA' if n_strands == 1 else 'dsDNA'
    unit = re.sub(r'\s', '', t_unit.lower())
    if unit == 'nm':
        c_t = t_value*1e-9
    elif unit in ('ng/ul', 'ng/\u00b5l', 'ngul'):
        c_t = t_value/(mw*1e3)
    else:
        raise SystemExit(f'Template unit must be nM or ng/uL, not {t_unit}.')
    tmpl = {'cT':c_t, 'strands':n_strands}

    cache = {}
    results = []
    for spec in specs:
        results.append(analyse_primer(plus, spec, th, rx, rx2, ee, opt, cache, tmpl, notes))
    for spec, res in zip(specs, results):
        report_primer(spec, res, L, opt, theta_pct, n_strands)

    # Amplicon from the two cognate sites
    cog = []
    for res in results:
        cog.append(next((s for s in res['sites'] if s['cg']), None))
    amp = None
    if cog[0] and cog[1]:
        sense = [c for c in cog if c['d'] == '-']
        anti = [c for c in cog if c['d'] == '+']
        if len(sense) == 1 and len(anti) == 1:
            a = sense[0]['a']
            b = anti[0]['b']
            amp = {'a':a, 'b':b, 'len':(b - a)%L + 1}
            print(f"\nAmplicon: {a}-{b} clockwise, {amp['len']} nt.")
        else:
            notes.append('The two cognate sites lie on the same strand, so no amplicon is drawn.')
    else:
        notes.append('At least one primer has no perfect-match site on this template, so no amplicon is drawn.')

    # Page payload
    floor_pct = min(args.theta_floor_pct, theta_pct)
    floor = floor_pct/100
    kept = []
    for res in results:
        kept += [s for s in res['sites'] if keep_site(s, res, floor, opt)]
    prim_json = []
    for spec, res in zip(specs, results):
        c = next((s for s in res['sites'] if s['cg']), None)
        prim_json.append({'name':spec['name'], 'seq':spec['seq'], 'n':res['n'], 'conc_nM':spec['conc_nM'],
                          'f':res['f_free'], 'f_nM':res['f_free']*1e9, 'RT':res['rt'], 'bound':res['bound'],
                          'endTol':res['end_tol'], 'cognate':({'d':c['d'], 'a':c['a'], 'b':c['b']} if c else None)})
    fb = free_mg_speciated(biv, dntp)
    meta = {'name':name, 'L':L, 'topology':topo, 'Ta':ta, 'monv_mM':monv_mM, 'biv_mM':biv_mM, 'dNTP_uM':dntp_uM,
            'freeBiv_mM':float(f'{fb*1e3:.3g}'), 'templ_nM':float(f'{c_t*1e9:.3g}'),
            'templ_ngul':float(f'{c_t*mw*1e3:.3g}'), 'proofClass':pc, 'proofread':proofread,
            'effApplied':'yes' if opt['use_eps'] else 'no', 'effSource':ee.source if opt['use_eps'] else 'not applied',
            'effCite':ee.cite if opt['use_eps'] else 'not applied', 'thetaPct':theta_pct, 'thetaFloorPct':floor_pct,
            'threatFactor':opt['threat_factor'], 'threatMax':opt['threat_max'], 'seedMax':opt['seed_max'],
            'seed3Len':opt['seed3'], 'seedIntLen':opt['seedint'], 'surveyRun':opt['survey_run'],
            'applyInit':opt['apply_init'], 'exact':opt['exact'], 'maxMarkers':25,
            'labels':{'reg':args.region_labels, 'thr':args.threat_labels, 'seed':args.seed_labels}}
    payload = {'meta':meta, 'primers':prim_json, 'sites':[site_json(s) for s in kept], 'amp':amp, 'notes':notes}

    out = args.out
    if not out:
        stem = re.sub(r'[^\w\-]', '_', f"{name}_{specs[0]['seq']}_{specs[1]['seq']}_HybMap")
        out = stem + '.html'
    build_html(out, payload)
    for nt in notes:
        print('\nNote: ' + nt)
    print(f'\nEmbedded {len(kept)} of {sum(len(r["sites"]) for r in results)} sites.')
    print(f'Saved interactive HTML map to: {os.path.abspath(out)}')


if __name__ == '__main__':
    main()
