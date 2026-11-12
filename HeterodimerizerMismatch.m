function [dGfrmin, fAlignr, dGfr, check, dGfrend] = HeterodimerizerMismatch(primj, primk, f0, r0, phi, plascon, mon, biv, dNTP, tmf, dTsub, tAnneal, tmplmin)
%% HETERODIMERIZERMISMATCH  Cross-dimer screen for a forward/reverse primer pair
% Cross-dimer screen for a forward/reverse primer pair.
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 12-11-2026 01:53 IST

% dGfrmin: threshold DeltaG at which the heterodimer would sit at phi times the reverse-primer-
%          template duplex concentration. See dGHetero.m.
% fAlignr: 3-row alignment block for the most stable cross-dimer register.
% dGfr: DeltaG of that register. 0 when no register carries three consecutive Watson-Crick pairs.
% check: false when the heterodimer is more stable than the threshold.
% dGfrend: DeltaG of the most stable register in which the FORWARD primer 3' terminus is paired,
%          i.e. the extendable subset. Call the function a second time with primj and primk swapped
%          to get the reverse-primer-anchored value, and take the more stable of the two. The optional
%          thirteenth argument tmplmin is the least number of partner nucleotides required 5' of the base
%          paired to the forward terminus, forwarded to DimerScan, which explains it. It defaults to 0.

% CHANGES
%  1. swalign has been replaced by DimerScan for the same reasons as in HomodimerizerMismatch: with
%     'GapOpen', 1000 the search space is just the set of ungapped registers, swalign picks from that
%     set by NUC44 identity score rather than by DeltaG, and it returns only one alignment, so a
%     short extendable 3' dimer is hidden behind any longer internal complementarity.
%  2. dGmaxMismatch was called even when the alignment came back empty. DimerScan returns a clean 0
%     sentinel instead.
%  3. The temperature offset is the argument dTsub rather than a hardcoded 5, defaulting to 5.
%  4. The rejection now also requires the dimer to be stabilising. dGfr = 0 means no register was
%     found, and if plascon is low enough for dGfrmin to go positive, a pair with no cross-
%     complementarity at all would otherwise be rejected.

%  5. BOTH PRIMER-TEMPLATE FREE ENERGIES ARE EVALUATED HERE, AT ONE TEMPERATURE.            (G2,H9)
%     dGHetero solves a single coupled equilibrium at one temperature and forms
%     Kr = exp(-dGrT/(R*tann)), so a free energy computed at any other temperature does not belong
%     in it. The forward and reverse scans each evaluate their own primer at its own reference,
%     which is the right thing for comparing candidates and the wrong thing to feed a joint
%     equilibrium, so both are recomputed here.
%
%     WHY IT MATTERS, measured on the 3940 nt run of 16-10-2026 with reverse primer
%     GTGAGCTAACTCACATTAATTGCGTTGC: dS = -0.631 kcal/(mol K) and tmr - tmf = 2.27 K, so dGr is
%     -14.48 at tmr - dTsub against -15.91 at the joint temperature, a gap of 1.43 kcal/mol. Kr
%     came out 4.0e9 where it should be 3.6e10, UNDERSTATED BY A FACTOR OF NINE, leaving the
%     heterodimer threshold about 1.4 kcal/mol too permissive. The RT denominator accounts for only
%     1.16x of that; the stale free energy is the rest.
%
%     THE dGf AND dGr ARGUMENTS ARE GONE.                                                     (H9)
%     They were advisory: recomputed here and then discarded, with a warning when the supplied
%     values drifted beyond 0.05 kcal/mol. That warning could not not fire. The caller stores dGf
%     at the forward primer's own reference and dGr at the reverse primer's, and LATE-PCR sets the
%     two Tm deliberately apart, so the drift is structural rather than accidental: measured at
%     0.571 kcal/mol per degree of offset, which over the engine's gated offset window is between
%     2 and 4 kcal/mol on every pair it will ever build, against a threshold of 0.05. A warning
%     that fires on every run teaches the reader to skim past warnings, so both the arguments and
%     the check are removed rather than retuned.

if (nargin < 11)||isempty(dTsub),  dTsub   = 5;  end
if (nargin < 12),                  tAnneal = []; end
if (nargin < 13)||isempty(tmplmin), tmplmin = 0; end

% EVALUATION TEMPERATURE.                                                                      (H9)
% A heterodimer involves both primers, which sit in one tube at one temperature, so both free
% energies must be taken at the same temperature and that temperature should be the one the
% reaction reaches. The older code used the forward primer's Tm minus dTsub, which approximates it
% but is the excess primer's convention imposed on a structure that belongs to neither primer.
% Falls back to that convention when the annealing temperature is not resolvable.
if ~isempty(tAnneal)&&isfinite(tAnneal)
    tann = tAnneal;
else
    tann = tmf - dTsub;
end

primj = char(primj); primk = char(primk);

dGfAt = dGPrim(primj, mon, biv, dNTP, tann);
dGrAt = dGPrim(primk, mon, biv, dNTP, tann);

dGfrmin = dGHetero(dGfAt, dGrAt, f0, r0, phi, plascon, tann);
[dGfr, dGfrend, fAlignr, ~] = DimerScan(primj, primk, mon, biv, dNTP, tann, [], tmplmin);

if (dGfr < 0)&&(dGfr < dGfrmin)
    check = false;
else
    check = true;
end
end
