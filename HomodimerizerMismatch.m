function [dgp, dgppmin, pAlignp, dgpp, check, dgppend] = HomodimerizerMismatch(sequence, phi, pcon, plascon, monv, biv, dNTP, tmf, dTsub, tAnneal, tmplmin)
%% HOMODIMERIZERMISMATCH  Homodimer screen for a single primer
% Homodimer screen for a single primer.
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 12-11-2026 01:53 IST

% dgp: DeltaG of the primer against its perfect Watson-Crick complement, i.e. the primer-template duplex
%      that dGHomo needs as its specific-binding reference.
% dgppmin: threshold DeltaG at which the homodimer would sit at exactly phi times the primer-template
%      duplex concentration. See dGHomo.m for the derivation and for what phi means numerically.
% pAlignp: 3-row alignment block for the most stable homodimer register.
% dgpp: DeltaG of that register. 0 when no register carries three consecutive Watson-Crick pairs.
% check: false when the homodimer is more stable than the threshold.
% dgppend: DeltaG of the most stable register in which the primer 3' terminus is paired. This is the
%      extendable subset and is returned so the caller can apply a separate 3'-dimer cutoff without a
%      second alignment pass. Ignore it if the caller only wants the original five outputs. The optional
%      eleventh argument tmplmin is the least number of partner nucleotides required 5' of the base paired to
%      the terminus, forwarded to DimerScan, which explains it. It defaults to 0.

% CHANGES
%  1. The single swalign call has been replaced by DimerScan, which enumerates every ungapped register
%     and returns the thermodynamic minimum. swalign with 'GapOpen', 1000 was searching that same space
%     but selecting from it by NUC44 identity score rather than by DeltaG, and returning one alignment
%     only, which hid 3'-anchored dimers behind longer internal complementarities.
%  2. The evaluation temperature offset is now the argument dTsub rather than a hardcoded 5, so it stays
%     in step with the dTsub constant in LocalPrimerSearchEngine.m. It defaults to 5 when omitted, so
%     existing eight-argument call sites are unaffected.
%  3. The "no dimer detected" sentinel is 0 on every path. The previous version returned 0 when swalign
%     came back empty but Inf when dGmaxMismatch found no run of three matches, and the Inf propagated
%     into the candidate logs.
%  4. The rejection now also requires the dimer to be stabilising. dgpp = 0 means no register was found,
%     and dGHomo's threshold goes positive at low enough template concentration, at which point a primer
%     with no self-complementarity at all would otherwise be rejected.

if (nargin < 9)||isempty(dTsub),   dTsub   = 5;  end
if (nargin < 10),                  tAnneal = []; end
if (nargin < 11)||isempty(tmplmin), tmplmin = 0; end

% EVALUATION TEMPERATURE.                                                                      (H9)
% A self-dimer has to be judged at the temperature the reaction actually reaches, which is one
% temperature for every primer in the tube. The older convention here was the primer's OWN Tm minus
% dTsub, which puts each candidate on equal footing in a symmetric reaction and breaks down in an
% asymmetric one. LATE-PCR sets the limiting primer's Tm deliberately ABOVE the excess primer's, by
% dTrfmin to dTrfmax in the engine, so that convention screened the limiting primer several degrees
% hotter than the annealing step it has to survive and made its dimers look weaker than they will
% be. Measured on real dimerising sequences the error runs 0.15 to 0.25 kcal/mol per degree, so a
% six-degree offset is worth 1 to 1.6 kcal/mol, which is the range the dimer thresholds themselves
% occupy.
% tAnneal is the reaction's own lowest annealing temperature. When it is not supplied or not
% resolvable, which happens when the thermal protocol cannot be parsed, the old per-primer
% convention is used so existing callers keep working rather than failing.
if ~isempty(tAnneal)&&isfinite(tAnneal)
    tann = tAnneal;
else
    tann = tmf - dTsub;
end

sequence = char(sequence);
dgp = dGPrim(sequence, monv, biv, dNTP, tann);
dgppmin = dGHomo(dgp, phi, pcon, plascon, tann);

[dgpp, dgppend, pAlignp, ~] = DimerScan(sequence, sequence, monv, biv, dNTP, tann, [], tmplmin);

if (dgpp < 0)&&(dgpp < dgppmin)
    check = false;
else
    check = true;
end
end
