function [dG, dH, dS, info] = dGTerminalMismatch(i1, i2, j1, j2, tann, N)
%% DGTERMINALMISMATCH  Free energy of the unpaired pair that closes a hairpin or internal loop
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 26-10-2026 11:55 UTC
% ____________________________________________________________________________________________________________
% PHYSICS
%
% A loop does not begin where the pairing stops. The first unpaired pair beyond the closing
% Watson-Crick pair still stacks on it, and that stack carries most of the sequence dependence of
% loop stability. This is not a refinement bolted onto the hairpin model; it is part of it.
% Equations 9 and 10 of the SantaLucia and Hicks review define hairpin stability as a loop-length
% term PLUS this increment, so a hairpin free energy computed without it is incomplete.
%
% The sign matters for a screen. Of the 48 measured combinations with a Watson-Crick closing pair,
% 46 are stabilising, spanning -1.19 kcal/mol at CG/GA to roughly zero at AC/TC. Substituting a
% conservative constant therefore does not merely add noise: it removes stabilisation that is almost
% always present, which makes every hairpin look weaker than it is. For a gate whose whole purpose
% is catching hairpins that survive to the annealing step, that error runs the wrong way.
%
% GEOMETRY AND NOTATION
%
% The review writes the motif as i1 i2 / j1 j2, meaning
%
%       5'- i1 i2 -3'
%       3'- j1 j2 -5'
%
% where i1.j1 is the closing Watson-Crick pair and i2.j2 is the mismatch. For a hairpin, i2 and j2
% are the first and last bases of the loop. For an internal loop, Equation 12 applies this twice,
% once at each end, each oriented so that its closing pair faces outward.
%
% THE CLAMP
%
% One entry, GG/CG, carries dH -0.70 kcal/mol with dS -19.2 cal/mol/K, giving dG37 +5.26. That is
% inconsistent with every neighbouring entry and lies far outside the range the review reports, whose
% upper end it gives as AC/TC. Rather than import a number that would distort any hairpin unlucky
% enough to close that way, the value is clamped to the largest physically supported entry in the
% table and info.Clamped is set. The entry is flagged in the CSV, not silently rewritten, so the
% defect stays visible to anyone reading the data.
%
% SALT
%
% Returned at 1 M NaCl. The review takes the terminal mismatch to carry the same salt dependence as
% a base pair while the remaining loop nucleotides are salt independent, so the motif contributes one
% pair to the caller's per-phosphate accounting. That count is reported in info.Phosphates rather
% than applied here, to avoid correcting the same phosphate twice.
%
% ____________________________________________________________________________________________________________
% INPUTS   i1, i2, j1, j2 as single characters; tann in degrees C; N the NNParams struct, fetched
%          if omitted
% OUTPUTS  dG kcal/mol at tann; dH, dS; info with fields Found, Clamped, Suspect, Phosphates
%
% ____________________________________________________________________________________________________________
% REFERENCES
%   SantaLucia J Jr, Hicks D. The thermodynamics of DNA structural motifs. Equations 9, 10 and 12.
%       Annu Rev Biophys Biomol Struct 2004;33:415-40.  doi:10.1146/annurev.biophys.32.110601.141800
%   The complete terminal-mismatch measurements are attributed there to S. Varma and
%       J. SantaLucia, manuscript in preparation, so no primary paper reproduces the table.
% ____________________________________________________________________________________________________________

if (nargin < 6)||isempty(N), N = NNParams(); end
dG = 0; dH = 0; dS = 0;
info = struct('Found', false, 'Clamped', false, 'Suspect', false, 'Phosphates', 0);

bases = 'ACGT';
q = zeros(1, 4);
in = {i1, i2, j1, j2};
for t = 1:4
    hit = strfind(bases, upper(in{t}(1)));
    if isempty(hit), return, end
    q(t) = hit(1);
end

if ~isfield(N, 'TM')||isempty(N.TM.dH), return, end
h = N.TM.dH(q(1), q(2), q(3), q(4));
s = N.TM.dS(q(1), q(2), q(3), q(4));
if ~isfinite(h)||~isfinite(s), return, end

info.Found = true;
info.Suspect = N.TM.Suspect(q(1), q(2), q(3), q(4));
info.Phosphates = 1;
dH = h; dS = s;
dG = dH - (tann + 273.15)*dS/1000;

if info.Suspect
    % Clamp in free energy at the working temperature, then hold dH and report the entropy that
    % reproduces the clamped value, so a caller that recombines dH and dS gets the same answer.
    cap = N.TM.MaxPhysical;
    if isfinite(cap) && (dG > cap)
        dG = cap;
        dS = 1000*(dH - dG)/(tann + 273.15);
        info.Clamped = true;
    end
end
end
