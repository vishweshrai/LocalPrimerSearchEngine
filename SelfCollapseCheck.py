#!/usr/bin/env python3
"""
SelfCollapseCheck.py
12-11-2026 02:17 IST
Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA

What This Script Does
Predicts 3'-OH self-priming and polymerase extension on the linear ssDNA scaffolds ((+) and (-)
strands) excised from a circular ssDNA template, and writes an interactive HTML map. The 3'-OH can
fold back onto its own strand (foldback) or pair with the same stem on a second copy of the strand
(bimolecular), and two enzyme channels judge whether the paired end is extended. Expansion is
outward-only around the mandatory core [a, b], with these controls:
  1. Forward Primer 5' Index [a - fFlex .. a]
  2. Forward Primer Length   [18 .. 24]
  3. Reverse Primer 5' Index [b .. b + rFlex]
  4. Reverse Primer Length   [20 .. 28]
  5. Active ssDNA Strand     [(+) Strand (3'-OH @ rPos) | (-) Strand (3'-OH @ fPos)]
  6. Enzyme Channel          [Taq | Proofreading (Exo)]

Physics
In asymmetric (LATE) PCR the accumulating single strand ends, at its 3' terminus, where the
opposite primer's 5' end sits: rPos for the (+) strand and fPos for the (-) strand. That 3'-OH can
fold back onto the strand's own body. If the fold leaves the 3' base paired and leaves at least one
template base beyond the stem, polymerase extends the 3'-OH and the strand lengthens in an
unintended way.

Each candidate fold is a stem plus a hairpin loop (SantaLucia and Hicks 2004). The stem is scored
by nearest neighbours with internal mismatches. The loop takes the Table 4 length term, the
terminal mismatch of the first and last loop bases, and the sequence bonus of a special triloop or
tetraloop. A triloop takes no terminal mismatch and takes a penalty on an A-T closing pair. Loops
beyond 30 nt take the worm-like-chain polymer increment of Aksel et al. 2024 and Dunn et al. 2015,
anchored at the 30 nt table value (see Loop Model below). Salt enters through the equivalent
monovalent concentration Na_eq = Mono + 1.2*sqrt(max(Mg - dNTP, 0)).

Loop Model
The mean squared end-to-end distance of a worm-like chain of contour length Lc = n*Lbase is
d2 = 2*Lp*Lc*(1 - (Lp/Lc)*(1 - exp(-Lc/Lp))), with Lbase = 0.6 nm per nt and Lp = 0.9 nm
(the expression is Rivetti 1998 eq. 5, the values are those of the Aksel et al. supplement; Rivetti
measured Lp = 1.3 nm for ssDNA, and using it moves the fold Tm by about 0.1 C). The effective
concentration of one chain end around the other is j(n) = (1/N_A)*(3/(2*pi*d2))^(3/2). The loop
penalty for n > 30 is dG(n) = dG_table(30) + R*T*(ln j(30) - ln j(n)), anchored at 30 nt so the
measured table is kept where it applies. The engine form, dG_table(30) + 2.44*R*T37*ln(n/30), is
carried beside it in the fields dg_f_eng, tm_eng and pe_eng, and is used for nothing else.

Persistence and Tiers (log-only)
A hairpin is unimolecular, so its melting temperature has no concentration dependence. The fold Tm
comes from dG_fold at T and T + 10 C (two-point van't Hoff). The fold persists when
Tm >= T_anneal - HpTmMargin (HpTmMargin = 5 C, as in the engine). Tier of one state:
  REJECT     persists and eps = 1 (the last three 3' bases are paired)
  FLAG       persists and FlagEps <= eps < 1 (a 3' mismatch limits extension)
  TRANSIENT  does not persist and eps >= FlagEps (extendable, melts at T_anneal)
  NONE       otherwise
The tier of a strand is the highest tier among its states that have template beyond the stem.
For each fold the ratio j/c compares the loop's effective concentration with the strand
concentration c; j/c > 1 means the hairpin outcompetes strand-strand binding.

No dangling-end term is applied. At the polymerase end the geometry differs from the free duplex
that Bommarito et al. measured, and the suite carries the same decision in dGDangle.m.

The fold equilibrium constant is Ku = exp(-dG_fold/RT), measured against the unfolded coil. The
extension efficiency eps comes from the 3'-terminal three bases (Rejali 2018) and is a kinetic
factor on the folded end. dG_eff = dG_fold - RT*ln(eps) is the free energy shown to rank folds.

Folds of one 3' tail compete for the same bases and exclude each other. Their partition function
is Z = 1 + sum(Ku_i), and the probability that the tail sits in fold i is Ku_i/Z. The headline
probability of an extendable foldback collapse, P_ext = sum(Ku_i*eps_i)/Z, is taken over the sites
that the selected enzyme channel can extend.

A state counts only when template lies beyond its stem. On the active strand the stem's 5'-most
paired base must lie inside the active amplicon with at least one base upstream of it.

Scan range. The window start r0 runs to the last position that keeps the window inside the strand,
and a tail base takes part only if pairing it leaves a loop of at least 3 nt, so the scan reaches
the tight hairpins at the very 3' end, whose loop is shorter than about 2(W - stem) + 3 nt. Every
physical state is kept, one per window, and a state is dropped only when dG_fold exceeds 13.5
kcal/mol.

Bimolecular Channel
The 3'-OH can also pair with the same stem on a second copy of the strand. That complex has no loop
and owes helix initiation, so its free energy dG_bi is the stem scored with the initiation and the
terminal A-T terms. With [S] the concentration of partner strands, referred to the 1 M standard
state,

      Kb = [S]*exp(-dG_bi/RT)          Ku = exp(-dG_fold/RT)

The states of one 3' tail exclude each other: Zall = 1 + sum(Ku) + sum(Kb). The extendable fraction
is sum(eps*(Ku + Kb))/Zall and the share carried by the bimolecular states is
sum(eps*Kb)/sum(eps*(Ku + Kb)). Partner strands are taken at their ceiling concentration and
binding does not deplete them, so the share is an upper bound.

A bimolecular complex persists when its melting temperature reaches T_anneal - HpTmMargin. The
temperature comes from Tm = dH/(dS + R*ln(ceff)), with dH and dS recovered from dG_bi at T and
T + 10 C. Tails and sites are equimolar, so the half-bound point is K*[S]/2 = 1, with
K = exp(-dG_bi/RT), and ceff = [S]/2. Polymerase extends the 3'-OH along the partner's body, so
the template rule is the same as for foldback.

[S] is the ceiling of the free strand that carries the stem. For the (+) strand it is the excess
primer concentration (--fwd-conc-nM). For the (-) strand it is the limiting primer concentration
(--rev-conc-nM), because the (-) strand cannot exceed the amount of limiting primer extended.

Enzyme Channels
Two channels are evaluated for every strand, because the end an enzyme will extend depends on it.
  Taq   has no 3' to 5' exonuclease. A site counts when at most 1 base of the tail hangs unpaired
        beyond the stem, and its efficiency is eps from the Rejali table.
  Exo   a proofreading exonuclease trims an unpaired 3' flap back to the outermost paired base and
        excises a mispaired terminus. Efficiency is 1 and a flap of up to --exo-end-tol bases
        counts.
P_ext(Taq) and P_ext(Exo) are the extendable fractions in each channel, over the same partition
function. A blend in which a fraction phi of the 3' end engagements is served by the proofreading
enzyme gives

      P_ext(blend) = (1 - phi)*P_ext(Taq) + phi*P_ext(Exo)

which is exact when every engagement is served by one enzyme molecule. NEB LongAmp Taq is a blend of
Taq and Deep Vent DNA polymerases whose proportions are not published, and NEB reports a fidelity
two-fold higher than Taq alone, which does not fix phi. The blend is therefore shown only when
--exo-blend is given. The tier of all channels takes the highest rank over both channels and both
topologies, so a design that fails under either enzyme is visible. The selector in the page sets
the channel of the headline numbers and of the table.

References
  SantaLucia J Jr, Hicks D. The thermodynamics of DNA structural motifs. Annu Rev Biophys Biomol
      Struct 2004;33:415-40. doi:10.1146/annurev.biophys.32.110601.141800
  Rejali NA, Moric E, Wittwer CT. The effect of single mismatches on primer extension. Clin Chem
      2018;64(5):801-809. doi:10.1373/clinchem.2017.282285
  Bommarito S, Peyret N, SantaLucia J Jr. Thermodynamic parameters for DNA sequences with dangling
      ends. Nucleic Acids Res 2000;28(9):1929-34. doi:10.1093/nar/28.9.1929
  Aksel T, Navarro EJ, Fong N, Douglas SM. Design principles for accurate folding of DNA origami.
      Proc Natl Acad Sci USA 2024;121(48):e2406769121. doi:10.1073/pnas.2406769121
  Dunn KE, Dannenberg F, Ouldridge TE, Kwiatkowska M, Turberfield AJ, Bath J. Guiding the folding
      pathway of DNA origami. Nature 2015;525:82-86. doi:10.1038/nature14860
  Rivetti C, Walker C, Bustamante C. Polymer chain statistics and conformational analysis of DNA
      molecules with bends or sections of different flexibility. J Mol Biol 1998;280:41-59.
      doi:10.1006/jmbi.1998.1830
  New England Biolabs. LongAmp Taq DNA Polymerase, M0323. Product page, accessed 12-11-2026:
      https://www.neb.com/en-us/products/m0323-longamp-taq-dna-polymerase
"""

import argparse, json, math, os, re
from typing import Dict, List, Tuple

R_GAS = 8.314e-3 / 4.184  # kcal/(mol*K)
T_REF = 310.15            # 37 C in K
N_AVO = 6.02214076e23     # 1/mol
L_BASE_NM = 0.6           # contour length per nt of ssDNA (nm)
L_PERS_NM = 0.9           # persistence length of ssDNA (nm)
MIN_RUN = 3 # Shortest Watson-Crick run that counts as a stem (nt)
TAQ_END_TOL = 1 # Unpaired 3' bases that a polymerase without exonuclease tolerates
FOLD_DG_CUT = 13.5 # A state whose dG_fold exceeds this (kcal/mol) carries no weight and is dropped
COMP =str.maketrans("ACGTacgt", "TGCAtgca")
WC_PAIRS = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}

def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]

# 1. Thermodynamic Parameter Tables (SantaLucia & Hicks 2004 / aPCR Suite)
WC_NN = {
    "AA": (-7.6, -21.3), "TT": (-7.6, -21.3), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7), "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0), "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2), "GC": (-9.8, -24.4), "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
INIT_PARAMS = (0.2, -5.7)
AT_TERM_PENALTY = (2.2, 6.9)

INTERNAL_MM = {
    "AG/TT": (1.0, 0.9), "AT/TG": (-2.5, -8.3), "CG/GT": (-4.1, -11.7), "CT/GG": (-2.8, -8.0),
    "GG/CT": (3.3, 10.4), "GT/CG": (-4.4, -12.3), "TG/AT": (-0.1, -1.7), "TT/AG": (-1.3, -5.3),
    "AA/TG": (-0.6, -2.3), "AG/TA": (-0.7, -2.3), "CA/GG": (-0.7, -2.3), "CG/GA": (-4.0, -13.2),
    "GA/CG": (-0.6, -1.0), "GG/CA": (0.5, 3.2), "TA/AG": (0.7, 0.7), "TG/AA": (3.0, 7.4),
    "AC/TT": (0.7, 0.2), "AT/TC": (-1.2, -6.2), "CC/GT": (-0.8, -4.5), "CT/GC": (-1.5, -6.1),
    "GC/CT": (2.3, 5.4), "GT/CC": (5.2, 13.5), "TC/AT": (1.2, 0.7), "TT/AC": (1.0, 0.7),
    "AA/TC": (2.3, 4.6), "AC/TA": (5.3, 14.6), "CA/GC": (1.9, 3.7), "CC/GA": (0.6, -0.6),
    "GA/CC": (5.2, 14.2), "GC/CA": (-0.7, -3.8), "TA/AC": (3.4, 8.0), "TC/AA": (7.6, 20.2),
    "AA/TA": (1.2, 1.7), "CA/GA": (-0.9, -4.2), "GA/CA": (-2.9, -9.8), "TA/AA": (4.7, 12.9),
    "AC/TC": (0.0, -4.4), "CC/GC": (-1.5, -7.2), "GC/CC": (3.6, 8.9), "TC/AC": (6.1, 16.4),
    "AG/TG": (-3.1, -9.5), "CG/GG": (-4.9, -15.3), "GG/CG": (-6.0, -15.8), "TG/AG": (1.6, 3.6),
    "AT/TT": (-2.7, -10.8), "CT/GT": (-5.0, -15.8), "GT/CT": (-2.2, -8.4), "TT/AT": (0.2, -1.5),
}

TERMINAL_MM = {
    "AA/TA": (-3.1, -7.8), "AA/TC": (-1.6, -4.0), "AA/TG": (-1.9, -4.4), "AC/TA": (-1.8, -3.8),
    "AC/TC": (-0.1, -0.5), "AC/TT": (-0.9, -1.7), "AG/TA": (-2.5, -5.9), "AG/TG": (-1.1, -2.1),
    "AG/TT": (-3.2, -8.7), "AT/TC": (-2.3, -6.3), "AT/TG": (-3.5, -9.4), "AT/TT": (-2.4, -6.5),
    "CA/GA": (-4.3, -10.7), "CA/GC": (-2.6, -5.9), "CA/GG": (-3.9, -9.6), "CC/GA": (-2.7, -6.0),
    "CC/GC": (-2.1, -5.1), "CC/GT": (-3.2, -8.0), "CG/GA": (-6.0, -15.5), "CG/GG": (-3.8, -9.5),
    "CG/GT": (-3.8, -9.0), "CT/GC": (-3.9, -10.6), "CT/GG": (-6.6, -18.7), "CT/GT": (-6.1, -16.9),
    "GA/CA": (-8.0, -22.5), "GA/CC": (-5.0, -13.8), "GA/CG": (-4.3, -11.1), "GC/CA": (-3.2, -7.1),
    "GC/CC": (-3.9, -10.6), "GC/CT": (-4.9, -13.5), "GG/CA": (-4.6, -11.4), "GG/CG": (-0.7, -19.2),
    "GG/CT": (-5.7, -15.9), "GT/CC": (-3.0, -7.8), "GT/CG": (-5.9, -16.1), "GT/CT": (-7.4, -21.2),
    "TA/AA": (-2.5, -6.3), "TA/AC": (-2.3, -5.9), "TA/AG": (-2.0, -4.7), "TC/AA": (-2.7, -7.0),
    "TC/AC": (-0.7, -1.3), "TC/AT": (-2.5, -6.3), "TG/AA": (-2.4, -5.8), "TG/AG": (-1.1, -2.7),
    "TG/AT": (-3.9, -10.5), "TT/AC": (-0.7, -1.2), "TT/AG": (-3.6, -9.8), "TT/AT": (-3.2, -8.9),
}

HAIRPIN_LOOP_37 = {3: 3.5, 4: 3.5, 5: 3.3, 6: 4.0, 7: 4.2, 8: 4.3, 9: 4.5, 10: 4.6, 12: 5.0, 14: 5.1, 16: 5.3, 18: 5.5, 20: 5.7, 25: 6.1, 30: 6.3}
INTERNAL_LOOP_37 = {3: 3.2, 4: 3.6, 5: 4.0, 6: 4.4, 7: 4.6, 8: 4.8, 9: 4.9, 10: 4.9, 12: 5.2, 14: 5.4, 16: 5.6, 18: 5.8, 20: 5.9, 25: 6.3, 30: 6.6}

SPECIAL_LOOPS = {
    "AGAAT": (-1.5, 0.0), "AGCAT": (-1.5, 0.0), "AGGAT": (-1.5, 0.0), "AGTAT": (-1.5, 0.0),
    "CGAAG": (-2.0, 0.0), "CGCAG": (-2.0, 0.0), "CGGAG": (-2.0, 0.0), "CGTAG": (-2.0, 0.0),
    "GGAAC": (-2.0, 0.0), "GGCAC": (-2.0, 0.0), "GGGAC": (-2.0, 0.0), "GGTAC": (-2.0, 0.0),
    "TGAAA": (-1.5, 0.0), "TGCAA": (-1.5, 0.0), "TGGAA": (-1.5, 0.0), "TGTAA": (-1.5, 0.0),
    "AAAAAT": (0.5, -0.65), "AAAACT": (0.7, 1.61), "AAACAT": (1.0, 1.61), "ACTTGT": (0.0, 4.19),
    "AGAAAT": (-1.1, 1.61), "AGAGAT": (-1.1, 1.61), "AGATAT": (-1.5, 1.61), "AGCAAT": (-1.6, 1.61),
    "AGCGAT": (-1.1, 1.61), "AGCTTT": (0.2, 1.61), "AGGAAT": (-1.1, 1.61), "AGGGAT": (-1.1, 1.61),
    "AGGGGT": (0.5, 0.64), "AGTAAT": (-1.6, 1.61), "AGTGAT": (-1.1, 1.61), "AGTTCT": (0.8, 1.61),
    "ATTCGT": (-0.2, 1.61), "ATTTGT": (0.0, 1.61), "ATTTTT": (-0.5, 1.61), "CAAAAG": (0.5, -1.29),
    "CAAACG": (0.7, 0.0), "CAACAG": (1.0, 0.0), "CAACCG": (0.0, 0.0), "CCTTGG": (0.0, 2.57),
    "CGAAAG": (-1.1, 0.0), "CGAGAG": (-1.1, 0.0), "CGATAG": (-1.5, 0.0), "CGCAAG": (-1.6, 0.0),
    "CGCGAG": (-1.1, 0.0), "CGCTTG": (0.2, 0.0), "CGGAAG": (-1.1, 0.0), "CGGGAG": (-1.0, 0.0),
    "CGGGGG": (0.5, -0.97), "CGTAAG": (-1.6, 0.0), "CGTGAG": (-1.1, 0.0), "CGTTCG": (0.8, 0.0),
    "CTTCGG": (-0.2, 0.0), "CTTTGG": (0.0, 0.0), "CTTTTG": (-0.5, 0.0), "GAAAAC": (0.5, -3.23),
    "GAAACC": (0.7, 0.0), "GAACAC": (1.0, 0.0), "GCTTGC": (0.0, 2.57), "GGAAAC": (-1.1, 0.0),
    "GGAGAC": (-1.1, 0.0), "GGATAC": (-1.6, 0.0), "GGCAAC": (-1.6, 0.0), "GGCGAC": (-1.1, 0.0),
    "GGCTTC": (0.2, 0.0), "GGGAAC": (-1.1, 0.0), "GGGGAC": (-1.1, 0.0), "GGGGGC": (0.5, -0.97),
    "GGTAAC": (-1.6, 0.0), "GGTGAC": (-1.1, 0.0), "GGTTCC": (0.8, 0.0), "GTTCGC": (-0.2, 0.0),
    "GTTTGC": (0.0, 0.0), "GTTTTC": (-0.5, 0.0), "TAAAAA": (0.5, 0.32), "TAAACA": (0.7, 1.61),
    "TAACAA": (1.0, 1.61), "TCTTGA": (0.0, 4.19), "TGAAAA": (-1.1, 1.61), "TGAGAA": (-1.1, 1.61),
    "TGATAA": (-1.6, 1.61), "TGCAAA": (-1.6, 1.61), "TGCGAA": (-1.1, 1.61), "TGCTTA": (0.2, 1.61),
    "TGGAAA": (-1.1, 1.61), "TGGGAA": (-1.1, 1.61), "TGGGGA": (0.5, 0.64), "TGTAAA": (-1.6, 1.61),
    "TGTGAA": (-1.1, 1.61), "TGTTCA": (0.8, 1.61), "TTTCGA": (-0.2, 1.61), "TTTTGA": (0.0, 1.61),
    "TTTTTA": (-0.5, 1.61),
}

REJALI_EPS = {
    ("A", "A"): (2.12e-4, 2.40e-3, 1.21e-1), ("A", "C"): (4.25e-3, 2.02e-2, 1.44e-1),
    ("A", "G"): (6.73e-5, 1.14e-3, 2.79e-2), ("A", "T"): (1.0, 1.0, 1.0),
    ("C", "A"): (7.82e-3, 4.36e-2, 2.24e-1), ("C", "C"): (1.98e-4, 6.62e-3, 2.57e-2),
    ("C", "G"): (1.0, 1.0, 1.0),             ("C", "T"): (1.25e-2, 6.34e-2, 1.81e-1),
    ("G", "A"): (4.71e-5, 4.51e-4, 9.60e-3), ("G", "C"): (1.0, 1.0, 1.0),
    ("G", "G"): (8.85e-4, 1.21e-2, 3.45e-1), ("G", "T"): (1.55e-2, 3.98e-2, 1.61e-1),
    ("T", "A"): (1.0, 1.0, 1.0),             ("T", "C"): (3.47e-3, 9.00e-3, 9.21e-2),
    ("T", "G"): (1.30e-2, 6.56e-2, 2.82e-1), ("T", "T"): (1.07e-3, 9.43e-3, 1.03e-1),
}

# 2. Thermodynamic Helpers
def monv_equiv(monv_M: float, biv_M: float, dntp_M: float) -> float:
    return max(monv_M + 1.2 * math.sqrt(max(biv_M - dntp_M, 0.0)), 1e-5)

def primer_thermo(seq: str, temp_c: float, na_eq: float, conc_M: float) -> Tuple[float, float, float]:
    n, tk = len(seq), temp_c + 273.15
    dh, ds = INIT_PARAMS[0], INIT_PARAMS[1]
    for end_b in (seq[0], seq[-1]):
        if end_b in "AT":
            dh += AT_TERM_PENALTY[0]; ds += AT_TERM_PENALTY[1]
    for k in range(n - 1):
        s_dh, s_ds = WC_NN[seq[k:k+2]]
        dh += s_dh; ds += s_ds
    ds_salt = ds + 0.368 * (n - 1) * math.log(na_eq)
    dg = dh - tk * (ds_salt / 1000.0)
    tm = (dh * 1000.0) / (ds_salt + R_GAS * 1000.0 * math.log(conc_M / 4.0)) - 273.15
    w = seq[-min(5, n):]
    dh3, ds3 = INIT_PARAMS[0], INIT_PARAMS[1]
    for end_b in (w[0], w[-1]):
        if end_b in "AT":
            dh3 += AT_TERM_PENALTY[0]; ds3 += AT_TERM_PENALTY[1]
    for k in range(len(w) - 1):
        s_dh, s_ds = WC_NN[w[k:k+2]]
        dh3 += s_dh; ds3 += s_ds
    return round(dg, 2), round(tm, 1), round(dh3 - T_REF * (ds3 / 1000.0), 2)

def ln_j_eff(n: int) -> float:
    # Natural log of the effective concentration (mol/L) of a worm-like-chain end around its partner
    lc = n * L_BASE_NM
    d2 = 2.0 * L_PERS_NM * lc * (1.0 - (L_PERS_NM / lc) * (1.0 - math.exp(-lc / L_PERS_NM)))
    j_per_m3 = (3.0 / (2.0 * math.pi * d2 * 1e-18)) ** 1.5
    return math.log(j_per_m3 / N_AVO / 1000.0)

def dg_loop(loop_type: str, n: int, temp_c: float, model: str = "polymer") -> float:
    tk = temp_c + 273.15
    tbl = HAIRPIN_LOOP_37 if loop_type == "hairpin" else INTERNAL_LOOP_37
    if n > 30 and loop_type == "hairpin" and model == "polymer":
        return tbl[30] * (tk / T_REF) + R_GAS * tk * (ln_j_eff(30) - ln_j_eff(n))
    if n in tbl:
        dg37 = tbl[n]
    elif n > 30:
        dg37 = tbl[30] + 2.44 * R_GAS * T_REF * math.log(n / 30.0)
    else:
        keys = sorted(tbl.keys())
        lo = max(k for k in keys if k < n)
        hi = min(k for k in keys if k > n)
        dg37 = tbl[lo] + (tbl[hi] - tbl[lo]) * ((n - lo) / (hi - lo))
    return dg37 * (tk / T_REF)

def dg_term_mm(i1: str, i2: str, j1: str, j2: str, temp_c: float) -> float:
    tk = temp_c + 273.15
    for k in (f"{i1}{i2}/{j1}{j2}", f"{j1}{j2}/{i1}{i2}"):
        if k in TERMINAL_MM:
            dh, ds = TERMINAL_MM[k]
            return dh - tk * (ds / 1000.0)
    return 0.0

def score_stem_duplex(s1: str, s2: str, temp_c: float, na_eq: float, init_terms: bool = False) -> float:
    # init_terms adds helix initiation and the terminal A-T penalty, which a duplex between two molecules owes
    n, tk = len(s1), temp_c + 273.15
    is_wc = [(s1[k], s2[k]) in WC_PAIRS for k in range(n)]
    dh, ds, dg_lp, k = 0.0, 0.0, 0.0, 0
    if init_terms:
        dh += INIT_PARAMS[0]; ds += INIT_PARAMS[1]
        for end_b in (s1[0], s1[-1]):
            if end_b in "AT":
                dh += AT_TERM_PENALTY[0]; ds += AT_TERM_PENALTY[1]
    while k < n - 1:
        if is_wc[k] and is_wc[k + 1]:
            s_dh, s_ds = WC_NN[s1[k:k+2]]
            dh += s_dh; ds += s_ds; k += 1
        elif is_wc[k] and not is_wc[k + 1]:
            m_end = k + 1
            while m_end < n and not is_wc[m_end]:
                m_end += 1
            if m_end >= n:
                break
            mm_len = m_end - k - 1
            if mm_len == 1:
                dh1, ds1 = INTERNAL_MM.get(f"{s1[k]}{s1[k+1]}/{s2[k]}{s2[k+1]}", (1.0, 0.0))
                dh2, ds2 = INTERNAL_MM.get(f"{s2[k+2]}{s2[k+1]}/{s1[k+2]}{s1[k+1]}", (1.0, 0.0))
                dh += dh1 + dh2; ds += ds1 + ds2
            else:
                dg_lp += dg_loop("internal", 2 * mm_len, temp_c)
                dg_lp += dg_term_mm(s1[k], s1[k+1], s2[k], s2[k+1], temp_c)
                dg_lp += dg_term_mm(s2[m_end], s2[m_end-1], s1[m_end], s1[m_end-1], temp_c)
            k = m_end
        else:
            k += 1
    ds_salt = ds + 0.368 * (n - 1) * math.log(na_eq)
    return dh - tk * (ds_salt / 1000.0) + dg_lp

def stem_energy(scaf: str, tail: str, opp: str, r0: int, i_s: int, j_s: int, temp_c: float, na_eq: float, init_terms: bool = False) -> float:
    # Stem of the fold, or of the bimolecular complex of the same stem when init_terms is true
    L_max, W = len(scaf), len(tail)
    t_end = L_max - W + j_s
    b_start = r0 + W - 1 - j_s
    dg_st = score_stem_duplex(tail[i_s:j_s+1], opp[i_s:j_s+1], temp_c, na_eq, init_terms)
    if (j_s < W - 1) and (b_start > 0): # Flap end: terminal mismatch of the unpaired 3' base
        dg_st += dg_term_mm(scaf[t_end], scaf[t_end+1], scaf[b_start], scaf[b_start-1], temp_c)
    return dg_st

def eval_fold(scaf: str, tail: str, opp: str, r0: int, i_s: int, j_s: int, temp_c: float, na_eq: float, model: str = "polymer") -> Tuple[float, float, float]:
    L_max, W, tk = len(scaf), len(tail), temp_c + 273.15
    t_start = L_max - W + i_s
    b_end = r0 + W - 1 - i_s
    loop_len = t_start - b_end - 1
    if loop_len < 3:
        return 1e6, 1e6, 1e6
    dg_st = stem_energy(scaf, tail, opp, r0, i_s, j_s, temp_c, na_eq)
    dg_lp = dg_loop("hairpin", loop_len, temp_c, model)
    c5, c3 = scaf[b_end], scaf[t_start]
    if loop_len == 3:
        if (c5, c3) in {("A", "T"), ("T", "A")}:
            dg_lp += 0.5
        sp = SPECIAL_LOOPS.get(c5 + scaf[b_end+1:t_start] + c3)
        if sp:
            dg_lp += sp[0] - tk * (sp[1] / 1000.0)
    else:
        dg_lp += dg_term_mm(c5, scaf[b_end+1], c3, scaf[t_start-1], temp_c)
        if loop_len == 4:
            sp = SPECIAL_LOOPS.get(c5 + scaf[b_end+1:t_start] + c3)
            if sp:
                dg_lp += sp[0] - tk * (sp[1] / 1000.0)
    return dg_st + dg_lp, dg_st, dg_lp

def ext_efficiency(tail: str, opp: str, i_s: int) -> Tuple[float, str]:
    # Rejali efficiency of a polymerase without exonuclease, from the 3'-terminal three bases of the tail
    w = len(tail)
    eps, gov = 1.0, f"{tail[-1]}.{opp[-1]}"
    for d in (0, 1, 2):
        idx = (w - 1) - d
        if idx < i_s:
            break
        bp, bt = tail[idx], opp[idx]
        if (bp, bt) not in WC_PAIRS:
            d_eps = REJALI_EPS.get((bp, bt), (1e-4, 1e-4, 1e-4))[d]
            if d_eps < eps:
                eps, gov = d_eps, (f"{bp}.{bt}" if d == 0 else f"{bp}.{bt}@{d}")
    return max(eps, 1e-7), gov

TIER_RANK = {"NONE": 0, "TRANSIENT": 1, "FLAG": 2, "REJECT": 3}

def classify_fold(eps: float, tm_c: float, temp_c: float, hp_tm_margin: float, flag_eps: float) -> str:
    persists = tm_c >= temp_c - hp_tm_margin
    if persists and eps >= 0.999:
        return "REJECT"
    if persists and eps >= flag_eps:
        return "FLAG"
    if (not persists) and eps >= flag_eps:
        return "TRANSIENT"
    return "NONE"

def fold_tm(dg_t1: float, dg_t2: float, temp_c: float, dt: float = 10.0) -> float:
    # Two-point melting temperature of a unimolecular fold from dG at T and T + dt (no concentration term)
    ds = (dg_t1 - dg_t2) / dt
    if abs(ds) <= 1e-7:
        return -99.0
    return (dg_t1 + (temp_c + 273.15) * ds) / ds - 273.15

def bimol_tm(dg_t1: float, dg_t2: float, t1_c: float, t2_c: float, ceff_M: float) -> float:
    # Melting temperature of one bimolecular complex: Tm = dH/(dS + R*ln(ceff)), dH and dS from dG at two temperatures
    k1, k2 = t1_c + 273.15, t2_c + 273.15
    ds = (dg_t1 - dg_t2)/(k2 - k1)
    dh = dg_t1 + k1*ds
    den = ds + R_GAS*math.log(ceff_M)
    if abs(den) < 1e-9:
        return float("nan")
    return dh/den - 273.15

# 3. Sequence Parsing & Circular Slicing
def parse_seq_file(path: str) -> Tuple[str, str]:
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = [ln.strip() for ln in f.read().splitlines()]
    name = os.path.splitext(os.path.basename(path))[0]
    if any(ln.startswith("LOCUS") for ln in lines) and any(ln.startswith("ORIGIN") for ln in lines):
        o_i = next(i for i, ln in enumerate(lines) if ln.startswith("ORIGIN"))
        e_i = next((i for i, ln in enumerate(lines) if ln.startswith("//") and i > o_i), len(lines))
        body = "".join(lines[o_i+1:e_i])
    elif any(ln.startswith(">") for ln in lines):
        h_i = [i for i, ln in enumerate(lines) if ln.startswith(">")]
        name = lines[h_i[0]][1:].split()[0] if lines[h_i[0]][1:].strip() else name
        body = "".join(lines[h_i[0]+1 : (h_i[1] if len(h_i) > 1 else len(lines))])
    else:
        body = "".join(lines)
    seq = re.sub(r"[^A-Za-z]", "", body).upper().replace("U", "T")
    return seq, name

def wrap_pos(idx_1based: int, L: int) -> int:
    return ((idx_1based - 1) % L) + 1

def extract_circular(seq: str, s1: int, e1: int) -> Tuple[str, List[int]]:
    L, s0, e0 = len(seq), (s1 - 1) % len(seq), (e1 - 1) % len(seq)
    coords, chars, idx = [], [], s0
    while True:
        coords.append(idx + 1); chars.append(seq[idx])
        if idx == e0 or len(coords) >= L:
            break
        idx = (idx + 1) % L
    return "".join(chars), coords

# 4. Scan 3' Terminus Across Maximal Admissible Scaffold (5' -> 3' along either strand)
def scan_rpos(max_scaf: str, max_coords: List[int], tail_win: int, temp_c: float, na_eq: float,
              strand_conc_M: float = 1e-6, exo_end_tol: int = 15, bimolecular: bool = True) -> Dict:
    L_max = len(max_scaf)
    W = min(tail_win, L_max // 3)
    tail = max_scaf[L_max - W:]
    rt = R_GAS*(temp_c + 273.15)
    # The scan runs at the wider of the two flap tolerances. Each site records its own flap, so the Taq
    # channel can drop the sites that only the exonuclease would extend.
    scan_tol = max(TAQ_END_TOL, min(exo_end_tol, W - MIN_RUN))
    sites: List[Dict] = []

    # Window start r0 runs to the last position that keeps the window inside the strand. A tail base k
    # can pair only if it leaves a loop of at least 3 nt: loop = L_max - 2W - r0 + 2k >= 3.
    for r0 in range(0, L_max - W + 1):
        opp = "".join(max_scaf[r0 + W - 1 - k] for k in range(W))
        k_min = max(0, -((L_max - 3 - 2*W - r0)//2))
        is_p = [(k >= k_min) and ((tail[k], opp[k]) in WC_PAIRS) for k in range(W)]
        p_idx = [k for k, p in enumerate(is_p) if p]
        if not p_idx or p_idx[-1] < (W - 1 - scan_tol):
            continue
        j_r = p_idx[-1]
        best_dg, best_sub = 1e6, None
        for i_s in p_idx:
            if i_s >= j_r:
                break
            run_m, cr = 0, 0
            for k in range(i_s, j_r + 1):
                cr = (cr + 1) if is_p[k] else 0
                run_m = max(run_m, cr)
            if run_m < MIN_RUN:
                continue
            dg_tot, dg_st, dg_lp = eval_fold(max_scaf, tail, opp, r0, i_s, j_r, temp_c, na_eq)
            if dg_tot < best_dg:
                best_dg, best_sub = dg_tot, (i_s, j_r, run_m, dg_tot, dg_st, dg_lp)
        if not best_sub or best_dg > FOLD_DG_CUT:
            continue
        i_s, j_r, run_m, dg_tot, dg_st, dg_lp = best_sub
        dg_t2, _, _ = eval_fold(max_scaf, tail, opp, r0, i_s, j_r, temp_c + 10.0, na_eq)
        tm_f = fold_tm(dg_tot, dg_t2, temp_c)
        dg_eng, _, _ = eval_fold(max_scaf, tail, opp, r0, i_s, j_r, temp_c, na_eq, "engine")
        dg_eng2, _, _ = eval_fold(max_scaf, tail, opp, r0, i_s, j_r, temp_c + 10.0, na_eq, "engine")
        tm_eng = fold_tm(dg_eng, dg_eng2, temp_c)
        end3_r = 0
        for k in range(W - 1, -1, -1):
            if is_p[k]: end3_r += 1
            else: break
        flap = (W - 1) - j_r
        eps_taq, term_p = ext_efficiency(tail, opp, i_s)
        term_x = f"{tail[j_r]}.{opp[j_r]} (exo+)"
        dg_eff = round(dg_tot - rt*math.log(eps_taq), 2)
        # Bimolecular state of the same stem on a second strand: initiation and terminal terms, no loop
        dg_bi = stem_energy(max_scaf, tail, opp, r0, i_s, j_r, temp_c, na_eq, True)
        tm_bi = None
        if bimolecular:
            dg_bi10 = stem_energy(max_scaf, tail, opp, r0, i_s, j_r, temp_c + 10.0, na_eq, True)
            tm_b = bimol_tm(dg_bi, dg_bi10, temp_c, temp_c + 10.0, strand_conc_M/2.0)
            tm_bi = None if math.isnan(tm_b) else round(tm_b, 1)
        b_s0, b_e0 = r0 + W - 1 - j_r, r0 + W - 1 - i_s
        loop_n = (L_max - W + i_s) - b_e0 - 1
        j_over_c = math.exp(ln_j_eff(loop_n))/strand_conc_M
        al_t, al_b = tail[i_s:], opp[i_s:]
        al_m = "".join("|" if (al_t[k], al_b[k]) in WC_PAIRS else ("." if k <= j_r - i_s else " ") for k in range(len(al_t)))
        sites.append({
            "b_s0": b_s0, "b_e0": b_e0, "p_s": max_coords[b_s0], "p_e": max_coords[b_e0],
            "t_from": i_s + 1, "t_to": j_r + 1, "stem": j_r - i_s + 1,
            "loop": loop_n, "run_m": run_m, "end3_r": end3_r, "flap": flap,
            "dg_st": round(dg_st, 2), "dg_lp": round(dg_lp, 2), "dg_f": round(dg_tot, 2),
            "tm": round(tm_f, 1), "eps_taq": eps_taq, "dg_eff": dg_eff, "term_p": term_p, "term_x": term_x,
            "dg_f_eng": round(dg_eng, 2), "tm_eng": round(tm_eng, 1), "j_over_c": float(f"{j_over_c:.4g}"),
            "dg_bi": round(dg_bi, 2), "tm_bi": tm_bi,
            "al_t": f"5'-{al_t}-3'", "al_m": f"   {al_m}   ", "al_b": f"3'-{al_b}-5'",
        })

    return {
        "tail_coords": [max_coords[L_max - W], max_coords[-1]],
        "strand_conc_M": strand_conc_M,
        "sites": sorted(sites, key=lambda x: x["dg_eff"]),
    }

# 5. Build Search Space & HTML Output
def build_html(out_path: str, payload: Dict) -> None:
    data_json = json.dumps(payload)
    html = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8"><title>3'-OH Collapse Map ((+) & (-) Strands)</title>
<style>
  body {{ font-family: system-ui, sans-serif; background: #0b1120; color: #f8fafc; margin: 0; padding: 18px; }}
  .bar, .card {{ background: #1e293b; border: 1px solid #334155; border-radius: 10px; padding: 14px 18px; margin-bottom: 14px; }}
  .controls {{ display: grid; grid-template-columns: repeat(6, 1fr) auto auto; gap: 10px; align-items: end; }}
  label {{ display: block; font-size: 0.72rem; color: #94a3b8; text-transform: uppercase; margin-bottom: 5px; font-weight: 600; }}
  select {{ width: 100%; padding: 8px; background: #0f172a; color: #fff; border: 1px solid #475569; border-radius: 6px; font-family: monospace; }}
  .risk-box {{ background: #0f172a; padding: 6px 12px; border-radius: 8px; border: 1px solid #334155; text-align: center; min-width: 135px; cursor: pointer; transition: border-color 0.15s; }}
  .risk-box.active-strand {{ border: 2px solid #38bdf8; background: #172554; }}
  .p-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 12px; margin-bottom: 14px; font-family: monospace; font-size: 0.82rem; }}
  .main {{ display: grid; grid-template-columns: 560px minmax(0, 1fr); gap: 16px; }}
  .stats {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 8px; margin-bottom: 12px; }}
  .tier {{ font-weight: 800; }} .t-REJECT {{ color: #ef4444; }} .t-FLAG {{ color: #f59e0b; }} .t-TRANSIENT {{ color: #38bdf8; }} .t-NONE {{ color: #10b981; }}
  .tile {{ background: #0f172a; padding: 8px 10px; border-radius: 6px; border: 1px solid #334155; font-family: monospace; }}
  .tile .k {{ font-size: 0.68rem; color: #94a3b8; }} .tile .v {{ font-size: 0.95rem; font-weight: 700; margin-top: 2px; }}
  .tbl-wrap {{ max-height: 320px; overflow: auto; border: 1px solid #334155; border-radius: 6px; }}
  .stats2 {{ margin-top: 12px; }} .sec {{ font-size: 0.68rem; color: #94a3b8; text-transform: uppercase; margin: 10px 0 6px; font-weight: 600; }}
  table {{ width: 100%; border-collapse: collapse; font-family: monospace; font-size: 0.78rem; }}
  th, td {{ padding: 6px 8px; border-bottom: 1px solid #1e293b; text-align: left; }}
  th {{ background: #0f172a; color: #94a3b8; position: sticky; top: 0; }}
  tr.row:hover, tr.row.active {{ background: #1e3a5f; cursor: pointer; }}
  .align {{ margin-top: 12px; background: #070b14; padding: 12px; border-radius: 6px; border: 1px solid #334155; font-family: monospace; font-size: 0.8rem; white-space: pre; overflow-x: auto; }}
  .leg {{ display: flex; flex-wrap: wrap; gap: 10px; font-size: 0.75rem; color: #94a3b8; margin-top: 8px; }}
  .sw {{ width: 11px; height: 11px; display: inline-block; border-radius: 2px; margin-right: 4px; vertical-align: middle; }}
</style></head><body>
<div class="bar" style="display:flex;justify-content:space-between;align-items:center;">
  <div><strong style="font-size:1.15rem;">Intramolecular 3'-OH Collapse & Extension Visualizer ((+) & (-) Strands)</strong><br>
  <span style="font-size:0.82rem;color:#94a3b8;">Template: <b>{payload['dnaName']}</b> ({payload['plasmidLen']} nt) | Mandatory Core [a, b]: <b>[{payload['sStart']}..{payload['sEnd']}]</b> ({payload['coreLen']} nt) | fFlex={payload['fFlex']}, rFlex={payload['rFlex']}</span></div>
  <div style="font-size:0.78rem;color:#94a3b8;text-align:right;">T={payload['params']['temp_c']}&deg;C | Na={payload['params']['monv_mM']}mM | Mg={payload['params']['biv_mM']}mM | dNTP={payload['params']['dntp_uM']}&mu;M</div>
</div>
<div class="bar controls">
  <div><label>1. Forward Primer 5' Index</label><select id="fpos"></select></div>
  <div><label>2. Forward Primer Length [18..24]</label><select id="flen"></select></div>
  <div><label>3. Reverse Primer 5' Index</label><select id="rpos"></select></div>
  <div><label>4. Reverse Primer Length [20..28]</label><select id="rlen"></select></div>
  <div><label>5. Active ssDNA Strand</label>
    <select id="strand">
      <option value="+" selected>(+) Strand (3'-OH @ rPos)</option>
      <option value="-">(-) Strand (3'-OH @ fPos)</option>
    </select>
  </div>
  <div><label>6. Enzyme Channel</label>
    <select id="channel">
      <option value="taq">Taq (no exonuclease)</option>
      <option value="exo">Proofreading (Exo)</option>
    </select>
  </div>
  <div class="risk-box active-strand" id="box-plus" title="Click to inspect (+) strand 3'-OH collapse">
    <div style="font-size:0.66rem;color:#94a3b8;">(+) 3'-COLLAPSE RISK</div>
    <div id="prob-plus" style="font-size:1.25rem;font-weight:800;">--%</div>
    <div id="sub-plus" style="font-size:0.66rem;color:#94a3b8;">--</div>
    <div id="all-plus" style="font-size:0.66rem;color:#94a3b8;">--</div>
  </div>
  <div class="risk-box" id="box-minus" title="Click to inspect (-) strand 3'-OH collapse">
    <div style="font-size:0.66rem;color:#94a3b8;">(-) 3'-COLLAPSE RISK</div>
    <div id="prob-minus" style="font-size:1.25rem;font-weight:800;">--%</div>
    <div id="sub-minus" style="font-size:0.66rem;color:#94a3b8;">--</div>
    <div id="all-minus" style="font-size:0.66rem;color:#94a3b8;">--</div>
  </div>
</div>
<div class="p-grid">
  <div class="card" style="margin:0;"><div style="color:#38bdf8;font-weight:700;" id="f-hdr">Forward Primer</div><div id="f-seq"></div><div style="color:#94a3b8;font-size:0.75rem;" id="f-meta"></div></div>
  <div class="card" style="margin:0;"><div style="color:#ec4899;font-weight:700;" id="r-hdr">Reverse Primer</div><div id="r-seq"></div><div style="color:#94a3b8;font-size:0.75rem;" id="r-meta"></div></div>
</div>
<div class="main">
  <div class="card" style="display:flex;flex-direction:column;align-items:center;">
    <svg id="map" viewBox="-290 -290 580 580" style="width:100%;max-width:530px;background:#070b14;border-radius:8px;"></svg>
    <div class="leg">
      <span><i class="sw" style="background:#0284c7;"></i>Core [a, b]</span>
      <span><i class="sw" style="background:#38bdf8;"></i>Outward Flex</span>
      <span><i class="sw" style="background:#a855f7;"></i>Fwd Primer</span>
      <span><i class="sw" style="background:#ec4899;"></i>Rev Primer</span>
      <span><i class="sw" style="background:#f97316;"></i>Active 3'-OH Tail</span>
    </div>
  </div>
  <div class="card">
    <div class="stats">
      <div class="tile"><div class="k">Scaffold Len</div><div class="v" id="st-len">--</div></div>
      <div class="tile"><div class="k">Active 3'-OH End</div><div class="v" id="st-flk">--</div></div>
      <div class="tile"><div class="k">Best &Delta;G_fold</div><div class="v" id="st-dg">--</div></div>
      <div class="tile"><div class="k">Best &Delta;G_eff</div><div class="v" id="st-dge">--</div></div>
      <div class="tile"><div class="k">Hairpin Tm</div><div class="v" id="st-tm">--</div></div>
      <div class="tile"><div class="k">Strand Tier (log-only)</div><div class="v tier" id="st-tier">--</div></div>
    </div>
    <div class="sec">Channels and Bimolecular State (log-only)</div>
    <div class="stats" style="margin-bottom:12px;">
      <div class="tile"><div class="k">P_ext, Taq</div><div class="v" id="st-ptaq">--</div></div>
      <div class="tile"><div class="k">P_ext, Exo</div><div class="v" id="st-pexo">--</div></div>
      <div class="tile"><div class="k" id="k-pbl">P_ext, Blend</div><div class="v" id="st-pbl">--</div></div>
      <div class="tile"><div class="k">P_ext, Foldback + Bimol.</div><div class="v" id="st-ptot">--</div></div>
      <div class="tile"><div class="k">Bimolecular Share</div><div class="v" id="st-bi">--</div></div>
      <div class="tile"><div class="k">Strongest Bimol. Tm</div><div class="v" id="st-tmbi">--</div></div>
      <div class="tile"><div class="k">Tier, All Channels</div><div class="v tier" id="st-tierall">--</div></div>
    </div>
    <div class="tbl-wrap"><table><thead><tr><th>#</th><th>Plasmid</th><th id="th-scaf">Scaf Pos (5'&rarr;3')</th><th>Loop</th><th>Flap</th><th>&Delta;G_st</th><th>&Delta;G_lp</th><th>&Delta;G_fold</th><th>&epsilon;</th><th>&Delta;G_eff</th><th>Tm</th><th>j/c</th><th>&Delta;G_bi</th><th>Tm_bi</th><th>Tier</th><th>P_ext</th></tr></thead><tbody id="tbody"></tbody></table></div>
    <div class="align" id="inspector">Click a row or chord to inspect the 3'-OH stem-loop alignment...</div>
  </div>
</div>
<script>
const D = {data_json};
const RT = (8.314e-3 / 4.184) * (D.params.temp_c + 273.15);
const fPosEl = document.getElementById('fpos'), fLenEl = document.getElementById('flen');
const rPosEl = document.getElementById('rpos'), rLenEl = document.getElementById('rlen');
const strandEl = document.getElementById('strand'), chEl = document.getElementById('channel');
chEl.value = D.params.channel;
const boxPlus = document.getElementById('box-plus'), boxMinus = document.getElementById('box-minus');
const svg = document.getElementById('map'), tbody = document.getElementById('tbody');

D.fwdPositions.forEach(x => {{
  const o = new Option(x.label, x.fPos); o.dataset.fidx = x.f_idx; o.dataset.out = x.out5;
  if (x.out5 === 0) o.selected = true; fPosEl.add(o);
}});
for (let l = 18; l <= 24; l++) fLenEl.add(new Option(l + ' nt', l, l === 20, l === 20));
D.revPositions.forEach(x => {{
  const o = new Option(x.label, x.rPos); o.dataset.out = x.out3;
  if (x.out3 === 0) o.selected = true; rPosEl.add(o);
}});
for (let l = 20; l <= 28; l++) rLenEl.add(new Option(l + ' nt', l, l === 22, l === 22));

boxPlus.onclick = () => {{ strandEl.value = '+'; update(); }};
boxMinus.onclick = () => {{ strandEl.value = '-'; update(); }};

function ang(p, L) {{ return -Math.PI/2 + ((p - 1) / L) * 2 * Math.PI; }}
function xy(r, a) {{ return [r * Math.cos(a), r * Math.sin(a)]; }}
function arc(r1, r2, s, e, L) {{
  let sp = e - s; if (sp < 0) sp += L; if (sp === 0) sp = 1; if (sp >= L - 1) return '';
  const a1 = ang(s, L), a2 = ang(s + sp, L), lg = (sp / L) > 0.5 ? 1 : 0;
  const [x1o, y1o] = xy(r2, a1), [x2o, y2o] = xy(r2, a2), [x2i, y2i] = xy(r1, a2), [x1i, y1i] = xy(r1, a1);
  return `M ${{x1o.toFixed(1)}} ${{y1o.toFixed(1)}} A ${{r2}} ${{r2}} 0 ${{lg}} 1 ${{x2o.toFixed(1)}} ${{y2o.toFixed(1)}} L ${{x2i.toFixed(1)}} ${{y2i.toFixed(1)}} A ${{r1}} ${{r1}} 0 ${{lg}} 0 ${{x1i.toFixed(1)}} ${{y1i.toFixed(1)}} Z`;
}}
function col(dge) {{
  const t = (Math.max(-5, Math.min(5, dge)) + 5) / 10;
  return `rgb(${{Math.round(239*(1-t)+34*t)}},${{Math.round(68*(1-t)+197*t)}},${{Math.round(68*(1-t)+94*t)}})`;
}}
function colorProb(el, p) {{
  el.textContent = p.toFixed(2) + '%';
  el.style.color = p >= 25 ? '#ef4444' : (p >= 5 ? '#f59e0b' : '#10b981');
}}

// States of one 3' tail exclude each other. The headline foldback probability uses Z = 1 + sum(Ku) over the
// sites of the selected channel: P_ext = sum(Ku*eps)/Z. The channel totals use every physical state of the
// tail, Zall = 1 + sum(Ku) + sum(Kb), with Kb = [S]*exp(-dG_bi/RT) for the same stem on a second strand.
// A state needs template beyond its stem, so its 5'-most paired base must lie past the amplicon start.
// Tier of one state: persistence by Tm >= T_anneal - HpTmMargin, extension by eps (log-only).
const TIER_NAMES = ['NONE', 'TRANSIENT', 'FLAG', 'REJECT'];
function tierRank(eps, tm) {{
  const persists = (tm !== null) && (tm >= D.params.temp_c - D.params.hp_tm_margin);
  if (persists && eps >= 0.999) return 3;
  if (persists && eps >= D.params.flag_eps) return 2;
  if ((!persists) && eps >= D.params.flag_eps) return 1;
  return 0;
}}
function collapse(allSites, start0, conc, ch) {{
  const P = D.params, taq = (ch === 'taq');
  const kept = allSites.filter(s => s.b_s0 > start0).map(s => {{
    const isTaq = s.flap <= P.taq_end_tol;
    return {{ ...s, isTaq, epsT: isTaq ? s.eps_taq : 0, epsC: taq ? s.eps_taq : 1, epsD: taq ? (isTaq ? s.eps_taq : 0) : 1,
      Ku: Math.exp(-s.dg_f/RT), Ke: Math.exp(-s.dg_f_eng/RT), Kb: P.bimolecular ? conc*Math.exp(-s.dg_bi/RT) : 0 }};
  }});
  const sum = (arr, f) => arr.reduce((a, s) => a + f(s), 0);

  // Headline: foldback in the selected channel
  const head = taq ? kept.filter(s => s.isTaq) : kept;
  const Z = 1 + sum(head, s => s.Ku), Ze = 1 + sum(head, s => s.Ke);
  let pExt = 0, pExtEng = 0, rk = 0;
  const sites = head.map(s => {{
    const pe = s.Ku*s.epsC/Z, t = tierRank(s.epsC, s.tm);
    pExt += pe; pExtEng += s.Ke*s.epsC/Ze; rk = Math.max(rk, t);
    return {{ ...s, sc_s: s.b_s0 - start0 + 1, sc_e: s.b_e0 - start0 + 1, pe, tier: TIER_NAMES[t],
      dgEff: taq ? s.dg_eff : s.dg_f, kbku: s.Kb/s.Ku }};
  }});
  sites.sort((a, b) => a.dgEff - b.dgEff);

  // Channel totals over every physical state of the tail: foldback, bimolecular and unfolded
  const Zall = 1 + sum(kept, s => s.Ku + s.Kb);
  const pExtTaq = sum(kept, s => (s.Ku + s.Kb)*s.epsT)/Zall;
  const pExtExo = sum(kept, s => s.Ku + s.Kb)/Zall;
  const phi = P.exo_blend;
  const pExtBlend = (phi === null) ? null : (1 - phi)*pExtTaq + phi*pExtExo;
  const tot = sum(kept, s => (s.Ku + s.Kb)*s.epsD);
  const pExtTot = tot/Zall;
  const biShare = tot > 0 ? sum(kept, s => s.Kb*s.epsD)/tot : 0;
  let tmBi = null, wMax = -1;
  if (P.bimolecular) kept.forEach(s => {{ const w = s.Kb*s.epsD; if (w > wMax) {{ wMax = w; tmBi = s.tm_bi; }} }});

  // Worst tier over foldback and bimolecular states, in the Taq and the Exo channel
  let rkAll = 0;
  kept.forEach(s => {{
    rkAll = Math.max(rkAll, tierRank(s.epsT, s.tm), tierRank(1, s.tm));
    if (P.bimolecular) rkAll = Math.max(rkAll, tierRank(s.epsT, s.tm_bi), tierRank(1, s.tm_bi));
  }});
  return {{ sites, pExt, pExtEng, tier: TIER_NAMES[rk], pExtTaq, pExtExo, pExtBlend, pExtTot, biShare, tmBi, tierAll: TIER_NAMES[rkAll] }};
}}
function fmtJ(x) {{ return x >= 100 ? x.toExponential(0) : x.toPrecision(2); }}
function fmtP(x) {{ return x === null ? '--' : (100*x).toPrecision(3) + '%'; }}

function update() {{
  const fp = fPosEl.value, fl = fLenEl.value, rp = rPosEl.value, rl = rLenEl.value, st = strandEl.value, ch = chEl.value;
  const fOpt = fPosEl.options[fPosEl.selectedIndex], rOpt = rPosEl.options[rPosEl.selectedIndex];
  const fIdx = parseInt(fOpt.dataset.fidx), out5 = parseInt(fOpt.dataset.out), out3 = parseInt(rOpt.dataset.out);
  const mOff = D.rFlex - out3;
  const fI = D.fwdPrimers[fp][fl], rI = D.revPrimers[rp][rl];
  const scanPlus = D.scansByRPos[rp], scanMinus = D.scansByFPos[fp];

  boxPlus.classList.toggle('active-strand', st === '+');
  boxMinus.classList.toggle('active-strand', st === '-');

  document.getElementById('f-hdr').textContent = `Forward Primer — plus(${{fI.fPos}}..${{fI.fEnd}}) [${{fl}} nt]`;
  document.getElementById('f-seq').textContent = `5'-${{fI.seq}}-3'`;
  document.getElementById('f-meta').textContent = `GC: ${{fI.gc}}% | Tm: ${{fI.Tm}}°C | dG: ${{fI.dG}} | dG3end: ${{fI.dG3}} kcal/mol`;
  document.getElementById('r-hdr').textContent = `Reverse Primer — revcomp(plus(${{rI.rLow}}..${{rI.rPos}})) [${{rl}} nt]`;
  document.getElementById('r-seq').textContent = `5'-${{rI.seq}}-3'`;
  document.getElementById('r-meta').textContent = `GC: ${{rI.gc}}% | Tm: ${{rI.Tm}}°C | dG: ${{rI.dG}} | dG3end: ${{rI.dG3}} kcal/mol`;

  const resPlus = collapse(scanPlus.sites, fIdx, scanPlus.strand_conc_M, ch), resMinus = collapse(scanMinus.sites, mOff, scanMinus.strand_conc_M, ch);
  const sitesPlus = resPlus.sites, sitesMinus = resMinus.sites;
  colorProb(document.getElementById('prob-plus'), resPlus.pExt * 100);
  colorProb(document.getElementById('prob-minus'), resMinus.pExt * 100);
  document.getElementById('sub-plus').textContent = `${{resPlus.tier}} | eng. loop ${{(resPlus.pExtEng*100).toFixed(2)}}%`;
  document.getElementById('sub-minus').textContent = `${{resMinus.tier}} | eng. loop ${{(resMinus.pExtEng*100).toFixed(2)}}%`;
  document.getElementById('all-plus').textContent = `all channels ${{resPlus.tierAll}}`;
  document.getElementById('all-minus').textContent = `all channels ${{resMinus.tierAll}}`;

  const scan = st === '+' ? scanPlus : scanMinus;
  const sites = st === '+' ? sitesPlus : sitesMinus;
  const termCoord = st === '+' ? rI.rPos : fI.fPos;

  document.getElementById('st-len').textContent = (D.coreLen + out5 + out3) + ' nt';
  document.getElementById('st-flk').textContent = `(${{st}}) @ ${{termCoord}}`;
  document.getElementById('st-dg').textContent = sites.length ? sites[0].dg_f + ' kcal/mol' : 'N/A';
  document.getElementById('st-dge').textContent = sites.length ? sites[0].dgEff.toFixed(2) + ' kcal/mol' : 'N/A';
  document.getElementById('st-tm').textContent = sites.length ? sites[0].tm + ' °C' : 'N/A';
  const stRes = st === '+' ? resPlus : resMinus, tEl = document.getElementById('st-tier');
  tEl.textContent = stRes.tier; tEl.className = 'v tier t-' + stRes.tier;
  document.getElementById('st-ptaq').textContent = fmtP(stRes.pExtTaq);
  document.getElementById('st-pexo').textContent = fmtP(stRes.pExtExo);
  document.getElementById('k-pbl').textContent = D.params.exo_blend === null ? 'P_ext, Blend (set --exo-blend)' : `P_ext, Blend (phi = ${{D.params.exo_blend}})`;
  document.getElementById('st-pbl').textContent = fmtP(stRes.pExtBlend);
  document.getElementById('st-ptot').textContent = fmtP(stRes.pExtTot);
  document.getElementById('st-bi').textContent = D.params.bimolecular ? fmtP(stRes.biShare) : 'off';
  document.getElementById('st-tmbi').textContent = (D.params.bimolecular && stRes.tmBi !== null) ? stRes.tmBi.toFixed(1) + ' °C' : '--';
  const taEl = document.getElementById('st-tierall');
  taEl.textContent = stRes.tierAll; taEl.className = 'v tier t-' + stRes.tierAll;
  document.getElementById('th-scaf').textContent = `(${{st}}) Pos (5'→3')`;

  const L = D.plasmidLen, rI_r = 195, rO_r = 222, rC = 186;
  svg.innerHTML = `<circle cx="0" cy="0" r="208.5" fill="none" stroke="#1e293b" stroke-width="27"/>`;
  const addP = (d, f, op=1) => {{ const p = document.createElementNS('http://www.w3.org/2000/svg','path'); p.setAttribute('d',d); p.setAttribute('fill',f); p.setAttribute('opacity',op); svg.appendChild(p); }};
  addP(arc(rI_r, rO_r, fI.fPos, rI.rPos, L), '#38bdf8', 0.85);
  addP(arc(rI_r, rO_r, D.sStart, D.sEnd, L), '#0284c7', 1.0);
  addP(arc(rO_r+2, rO_r+9, fI.fPos, fI.fEnd, L), '#a855f7', 1.0);
  addP(arc(rO_r+2, rO_r+9, rI.rLow, rI.rPos, L), '#ec4899', 1.0);
  if (st === '+') {{
    addP(arc(rI_r-6, rI_r, scan.tail_coords[0], scan.tail_coords[1], L), '#f97316', 1.0);
  }} else {{
    addP(arc(rI_r-6, rI_r, scan.tail_coords[1], scan.tail_coords[0], L), '#f97316', 1.0);
  }}

  const [x3, y3] = xy(rC, ang(scan.tail_coords[1], L));
  tbody.innerHTML = '';
  [...sites].reverse().forEach((s, ri) => {{
    const rk = sites.length - ri, [xh, yh] = xy(rC, ang(s.p_s, L)), c = col(s.dgEff);
    const ch = document.createElementNS('http://www.w3.org/2000/svg','path');
    ch.setAttribute('d', `M ${{x3.toFixed(1)}} ${{y3.toFixed(1)}} Q 0 0 ${{xh.toFixed(1)}} ${{yh.toFixed(1)}}`);
    ch.setAttribute('fill','none'); ch.setAttribute('stroke', c);
    ch.setAttribute('stroke-width', Math.max(1.2, Math.min(4.2, (5 - s.dgEff)*0.55)));
    ch.setAttribute('opacity', rk <= 5 ? '0.9' : '0.45'); ch.style.cursor = 'pointer';
    ch.onclick = () => pick(s, rk, termCoord, st); svg.appendChild(ch);
  }});
  sites.forEach((s, i) => {{
    const tr = document.createElement('tr'); tr.className = 'row'; tr.id = 'r-' + (i+1);
    tr.innerHTML = `<td>${{i+1}}</td><td>${{s.p_s}}..${{s.p_e}}</td><td>${{s.sc_s}}..${{s.sc_e}}</td><td>${{s.loop}}</td><td>${{s.flap}}</td><td>${{s.dg_st}}</td><td>+${{s.dg_lp}}</td><td>${{s.dg_f}}</td><td>${{s.epsC.toExponential(1)}}</td><td style="color:${{col(s.dgEff)}};font-weight:700;">${{s.dgEff.toFixed(2)}}</td><td>${{s.tm}}</td><td>${{fmtJ(s.j_over_c)}}</td><td>${{s.dg_bi}}</td><td>${{s.tm_bi === null ? '--' : s.tm_bi}}</td><td class="tier t-${{s.tier}}">${{s.tier}}</td><td>${{(s.pe*100).toFixed(2)}}%</td>`;
    tr.onclick = () => pick(s, i+1, termCoord, st); tbody.appendChild(tr);
  }});
  if (sites.length) {{
    pick(sites[0], 1, termCoord, st);
  }} else {{
    document.getElementById('inspector').textContent = `No qualifying (${{st}}) strand 3'-OH collapse sites found under the current cutoff.`;
  }}
}}

function pick(s, rk, termCoord, st) {{
  document.querySelectorAll('.row').forEach(r => r.classList.remove('active'));
  const rw = document.getElementById('r-' + rk); if (rw) rw.classList.add('active');
  const taq = (chEl.value === 'taq'), cNM = st === '+' ? D.params.strand_conc_plus_nM : D.params.strand_conc_minus_nM;
  const biTxt = D.params.bimolecular ? `dG_bi = ${{s.dg_bi}} kcal/mol | Tm_bi = ${{s.tm_bi === null ? '--' : s.tm_bi}} °C at ceff = ${{cNM/2}} nM | Kb/Ku = ${{fmtJ(s.kbku)}}` : 'off';
  document.getElementById('inspector').textContent =
`Rank #${{rk}} [(${{st}}) Strand, ${{taq ? 'Taq' : 'Exo'}} channel] | (${{st}}) Body Pos ${{s.sc_s}}..${{s.sc_e}} (Plasmid ${{s.p_s}}..${{s.p_e}}) | Loop = ${{s.loop}} nt | Stem = ${{s.stem}} bp (Max WC run = ${{s.run_m}}, 3' WC run = ${{s.end3_r}}) | Flap = ${{s.flap}} nt
dG_stem = ${{s.dg_st}} | dG_loop = +${{s.dg_lp}} | dG_fold = ${{s.dg_f}} kcal/mol | Tm = ${{s.tm}} °C | 3' Pair = ${{taq ? s.term_p : s.term_x}} (eps = ${{s.epsC.toExponential(2)}}) => dG_eff = ${{s.dgEff.toFixed(2)}} kcal/mol
Tier = ${{s.tier}} (persists if Tm >= ${{(D.params.temp_c - D.params.hp_tm_margin).toFixed(1)}} °C) | j/c = ${{fmtJ(s.j_over_c)}} at c = ${{cNM}} nM (above 1 the hairpin outcompetes strand-strand binding)
Bimolecular state of the same stem on a second strand (log only): ${{biTxt}}
Engine loop form (log only): dG_fold = ${{s.dg_f_eng}} kcal/mol, Tm = ${{s.tm_eng}} °C

(${{st}}) 3'-Tail Window (pos ${{s.t_from}}..${{s.t_to}}, 3'-OH @ ${{termCoord}}):  ${{s.al_t}}
                                                     ${{s.al_m}}
(${{st}}) Scaffold Body  (pos ${{s.sc_e}}..${{s.sc_s}}):                ${{s.al_b}}`;
}}
[fPosEl, fLenEl, rPosEl, rLenEl, strandEl, chEl].forEach(el => el.addEventListener('change', update));
update();
</script></body></html>"""
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)

def ask(val, msg: str, cast=str, default=None):
    if val is not None:
        return val
    raw = input(f"{msg}" + (f" [{default}]: " if default is not None else ": ")).strip()
    return cast(raw) if raw else default

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seq-file", type=str, default=None)
    ap.add_argument("--start", type=int, default=None, help="Mandatory scaffold start a (1-based)")
    ap.add_argument("--end", type=int, default=None, help="Mandatory scaffold end b (1-based)")
    ap.add_argument("--fflex", type=int, default=None, help="Outward 5' flexibility fFlex (nt)")
    ap.add_argument("--rflex", type=int, default=None, help="Outward 3' flexibility rFlex (nt)")
    ap.add_argument("--temp", type=float, default=58.0)
    ap.add_argument("--monv", type=float, default=40.0)
    ap.add_argument("--biv", type=float, default=2.5)
    ap.add_argument("--dntp", type=float, default=1200.0)
    ap.add_argument("--tail-win", type=int, default=18)
    ap.add_argument("--proof-exo", action="store_true", help="Open the page in the proofreading (Exo) channel; both channels are always computed")
    ap.add_argument("--hp-tm-margin", type=float, default=5.0, help="Fold persists when Tm >= temp - this (C)")
    ap.add_argument("--flag-eps", type=float, default=1e-3, help="Lowest eps that still earns a FLAG or TRANSIENT tier")
    ap.add_argument("--fwd-conc-nM", type=float, default=1000.0, help="Excess (forward) primer concentration (nM), the ceiling of the (+) strand")
    ap.add_argument("--rev-conc-nM", type=float, default=20.0, help="Limiting (reverse) primer concentration (nM), the ceiling of the (-) strand")
    ap.add_argument("--no-bimolecular", action="store_true", help="Leave the bimolecular state out of the channel totals")
    ap.add_argument("--exo-blend", type=float, default=None, help="Fraction phi of 3' end engagements served by the proofreading enzyme (0 to 1); the blend is shown only when set")
    ap.add_argument("--exo-end-tol", type=int, default=15, help="Unpaired 3' bases that the proofreading channel counts (clamped to window - 3)")
    ap.add_argument("--out", type=str, default="scaffold_3prime_collapse_map.html")
    args = ap.parse_args()
    if (args.exo_blend is not None) and not (0.0 <= args.exo_blend <= 1.0):
        ap.error("--exo-blend must lie between 0 and 1")
    bimol = not args.no_bimolecular

    seq_file = ask(args.seq_file, "Circular ssDNA sequence file (.fasta/.gb/.txt)", str)
    seq, dna_name = parse_seq_file(seq_file)
    L = len(seq)
    a = ask(args.start, f"Scaffold start index a (1..{L})", int, 2173)
    b = ask(args.end, f"Scaffold end index b (1..{L})", int, min(L, 6112))
    f_flex = ask(args.fflex, "Outward 5' flexibility fFlex (nt)", int, 20)
    r_flex = ask(args.rflex, "Outward 3' flexibility rFlex (nt)", int, 20)

    na_eq = monv_equiv(args.monv * 1e-3, args.biv * 1e-3, args.dntp * 1e-6)
    core_seq, _ = extract_circular(seq, a, b)

    fwd_pos, fwd_primers, scans_by_fpos = [], {}, {}
    r_latest = wrap_pos(b + r_flex, L)
    for f_idx, i_unw in enumerate(range(a - f_flex, a + 1)):
        fp = wrap_pos(i_unw, L)
        out5 = a - i_unw
        fwd_pos.append({"f_idx": f_idx, "fPos": fp, "out5": out5, "label": f"{fp} (a - {out5} nt)" if out5 else f"{fp} (exact a)"})
        fwd_primers[str(fp)] = {}
        for fl in range(18, 25):
            fe = wrap_pos(fp + fl - 1, L)
            f_s, _ = extract_circular(seq, fp, fe)
            dg, tm, dg3 = primer_thermo(f_s, args.temp, na_eq, args.fwd_conc_nM*1e-9)
            fwd_primers[str(fp)][str(fl)] = {"fPos": fp, "fEnd": fe, "seq": f_s, "gc": round(100*(f_s.count("G")+f_s.count("C"))/fl, 1), "dG": dg, "Tm": tm, "dG3": dg3}
        # (-) strand runs 5' -> 3' from r_latest down to fp (3'-OH at fp)
        plus_scaf, plus_coords = extract_circular(seq, fp, r_latest)
        minus_scaf = revcomp(plus_scaf)
        minus_coords = plus_coords[::-1]
        scans_by_fpos[str(fp)] = scan_rpos(minus_scaf, minus_coords, args.tail_win, args.temp, na_eq, args.rev_conc_nM*1e-9, args.exo_end_tol, bimol)

    rev_pos, rev_primers, scans_by_rpos = [], {}, {}
    f_earliest = wrap_pos(a - f_flex, L)
    for r_idx, k_unw in enumerate(range(b, b + r_flex + 1)):
        rp = wrap_pos(k_unw, L)
        out3 = k_unw - b
        rev_pos.append({"r_idx": r_idx, "rPos": rp, "out3": out3, "label": f"{rp} (b + {out3} nt)" if out3 else f"{rp} (exact b)"})
        rev_primers[str(rp)] = {}
        for rl in range(20, 29):
            r_low = wrap_pos(rp - rl + 1, L)
            p_foot, _ = extract_circular(seq, r_low, rp)
            r_s = revcomp(p_foot)
            dg, tm, dg3 = primer_thermo(r_s, args.temp, na_eq, args.rev_conc_nM*1e-9)
            rev_primers[str(rp)][str(rl)] = {"rPos": rp, "rLow": r_low, "seq": r_s, "gc": round(100*(r_s.count("G")+r_s.count("C"))/rl, 1), "dG": dg, "Tm": tm, "dG3": dg3}
        # (+) strand runs 5' -> 3' from f_earliest up to rp (3'-OH at rp)
        max_scaf, max_coords = extract_circular(seq, f_earliest, rp)
        scans_by_rpos[str(rp)] = scan_rpos(max_scaf, max_coords, args.tail_win, args.temp, na_eq, args.fwd_conc_nM*1e-9, args.exo_end_tol, bimol)

    build_html(args.out, {
        "dnaName": dna_name, "plasmidLen": L, "sStart": a, "sEnd": b,
        "coreLen": len(core_seq), "fFlex": f_flex, "rFlex": r_flex,
        "params": {"temp_c": args.temp, "monv_mM": args.monv, "biv_mM": args.biv, "dntp_uM": args.dntp,
                   "hp_tm_margin": args.hp_tm_margin, "flag_eps": args.flag_eps,
                   "strand_conc_plus_nM": args.fwd_conc_nM, "strand_conc_minus_nM": args.rev_conc_nM,
                   "bimolecular": bimol, "exo_blend": args.exo_blend, "exo_end_tol": args.exo_end_tol,
                   "taq_end_tol": TAQ_END_TOL, "channel": "exo" if args.proof_exo else "taq"},
        "fwdPositions": fwd_pos, "revPositions": rev_pos,
        "fwdPrimers": fwd_primers, "revPrimers": rev_primers,
        "scansByRPos": scans_by_rpos, "scansByFPos": scans_by_fpos,
    })
    print(f"Saved interactive HTML map to: {os.path.abspath(args.out)}")

if __name__ == "__main__":
    main()