function [dG, dH, dS, nPhos] = dGDangle(paired, unpaired, side, tann, N)
%% DGDANGLE  Free energy of one unpaired base stacked on the end of a duplex
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 05-11-2026 22:10 IST
% ____________________________________________________________________________________________________________
% STATUS: COMPUTED, VALIDATED, AND DELIBERATELY NOT APPLIED
%
% Nothing in the suite calls this function. NNParams loads the table, asserts its enthalpies
% against an independently transcribed Published_dH column and its entropies through the published
% dG37, and the sigreq row keeps a stale copy from going unnoticed. No free energy in the screen
% carries a dangling-end term. That is a decision taken on 05-11-2026, not an oversight, and these
% are the three reasons.
%
% 1. AT THE END THAT MATTERS, THE GEOMETRY IS NOT THE ONE THAT WAS MEASURED. Bommarito's values
%    are for a free duplex melting at equilibrium, where the unpaired base stacks on the helix
%    terminus with nothing else bound. At a site being scored for PRIMING the 3' terminus sits in
%    the polymerase active site. Whatever the first unpaired template base is doing there, it is
%    not the free stack that was measured, and no measurement of that geometry exists.
%
% 2. THE CORRECTION IS LEAST USEFUL WHERE THE SCREEN DOES ITS WORK. See the cancellation note
%    below: a dangling end stabilises the intended site and an off-target site alike, so it largely
%    cancels from a difference and survives mainly in absolute free energies and melting
%    temperatures. The mispriming screen thresholds margins, which is the quantity the note says
%    loses least by omitting it.
%
% 3. SWITCHING IT ON WOULD ADD HALF TO ONE KCAL/MOL PER END, unvalidated, to every off-target free
%    energy in the run, in a direction set by a table measured in a different geometry. That is
%    larger than several of the terms the screen is trying to resolve against each other.
%
% THE RULE IF IT IS EVER SWITCHED ON: 5' SIDE ONLY.
% The 5' end of a primer-template duplex has template running past it with nothing bound, which is
% exactly the geometry Bommarito measured, so a 5' dangle there is defensible. The 3' end is the
% one in the enzyme. Applying both sides would mean applying a measured number at one end and a
% guess at the other while calling them one correction. Melting-temperature work on a duplex that
% is not being extended is a different case and may use both ends.
%
% ____________________________________________________________________________________________________________
% PHYSICS
%
% A duplex does not end in a vacuum. When a primer anneals to a long template, template continues
% past both ends of the paired region, and the first unpaired base on each side stacks on the
% terminal pair. The stack is real and it is worth having: across the 32 measured combinations the
% contribution runs from +0.48 to -0.96 kcal/mol, averaging -0.45 for a base hanging off the 5' end
% and -0.29 off the 3' end. At the stabilising extreme an unpaired neighbour is worth as much as a
% whole A-T pair, so a duplex free energy computed without it is biased, not merely imprecise.
%
% The stabilisation comes from base stacking rather than hydrogen bonding, so it depends on the
% identity of the unpaired base, on the pair it stacks against, and on which side it sits. Two of
% the DNA values are destabilising, which happens in no RNA dangling end, and B-form DNA shows no
% systematic preference between the two ends while A-form RNA strongly prefers 3'. Those two facts
% are why an RNA table cannot stand in for this one.
%
% WHERE THE BIAS CANCELS AND WHERE IT DOES NOT
%
% A dangling end stabilises the intended site and a mispriming site alike, so a quantity computed as
% a DIFFERENCE between the two loses less than either absolute number does. The cancellation is
% partial, because the flanking bases differ between sites and the spread across the table is
% 1.44 kcal/mol. Treat this as a correction that matters most to absolute free energies and to
% melting temperatures, and least to margins.
%
% SALT
%
% The table is quoted at 1 M NaCl and is returned at 1 M NaCl. The motif adds one phosphate, which
% is reported as nPhos so the caller can fold it into whatever per-phosphate correction it already
% applies to the whole structure. Applying a correction here as well would count it twice.
%
% ____________________________________________________________________________________________________________
% INPUTS
%   paired    the terminal PAIRED base of the duplex, on the strand that carries the unpaired base
%   unpaired  the unpaired base itself
%   side      3 or 5, which end of the paired strand the unpaired base sits on
%   tann      temperature, degrees C
%   N         the NNParams struct. Fetched if omitted
%
% OUTPUTS
%   dG        kcal/mol at tann and 1 M NaCl. Zero when the combination is not in the measured set
%   dH, dS    kcal/mol and cal/mol/K
%   nPhos     phosphates this motif adds, for the caller's salt accounting. 1 when a value was
%             found, 0 when it was not
%
% ____________________________________________________________________________________________________________
% REFERENCES
%   Bommarito S, Peyret N, SantaLucia J Jr. Thermodynamic parameters for DNA sequences with
%       dangling ends. Nucleic Acids Res 2000;28:1929-34.  doi:10.1093/nar/28.9.1929
%   SantaLucia J Jr, Hicks D. The thermodynamics of DNA structural motifs. Table 3.
%       Annu Rev Biophys Biomol Struct 2004;33:415-40.  doi:10.1146/annurev.biophys.32.110601.141800
% ____________________________________________________________________________________________________________

if (nargin < 5)||isempty(N), N = NNParams(); end
dG = 0; dH = 0; dS = 0; nPhos = 0;

bases = 'ACGT';
ip = strfind(bases, upper(paired(1)));
ix = strfind(bases, upper(unpaired(1)));
if isempty(ip)||isempty(ix), return, end

% The measured set covers a dangling base against a WATSON-CRICK terminal pair only. A dangle
% stacked on a terminal mismatch was never measured, so the partner is taken as the complement of
% the paired base and nothing is invented for other cases.
partner = [4 3 2 1];
ic = partner(ip);

if ~ismember(side, [3 5])
    error('dGDangle:BadSide', 'side must be 3 or 5, not %s.', mat2str(side));
end
slot = 1 + (side == 5);

if ~isfield(N, 'Dangle')||isempty(N.Dangle.dH), return, end
h = N.Dangle.dH(ip, ix, ic, slot);
s = N.Dangle.dS(ip, ix, ic, slot);
if ~isfinite(h)||~isfinite(s), return, end

dH = h; dS = s; nPhos = 1;
dG = dH - (tann + 273.15)*dS/1000;
end
