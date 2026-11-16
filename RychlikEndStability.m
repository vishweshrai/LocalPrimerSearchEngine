function dG = RychlikEndStability(prim, endwin)
%% Rychlik 3' End Stability
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 14-11-2026 20:43 IST
%
% Revision history is not kept here. See RychlikEndStability_<YYYYMMDD>.txt beside this file.
% ____________________________________________________________________________________________________________
% What This Function Does
%
% Returns the stability of the 3'-terminal window of a primer on the scale Rychlik (1995) uses for
% his primer selection rule: oligonucleotides whose 3'-terminal pentamer is less stable than
% -9 kcal/mol are more likely to be specific primers. The value is the sum of the nearest-neighbour
% steps of his Table 1, the simplified Breslauer et al. (1986) free energies at 25 C, with no helix
% initiation and no terminal term. For a pentamer that is four steps. The Primer3 manual quotes the
% same quantity as Max 3' Stability on its Breslauer table (GCGCG 13.4, TATAC 4.6), and Primer3Plus
% ships the limit at 9.0, so a pentamer more negative than -9.0 here is what Primer3Plus rejects.
%
% The scale is not the SantaLucia scale that EndStability in LocalPrimerSearchEngine returns, and
% the two cannot be compared by threshold. On all 1024 pentamers they correlate at r = 0.90 with a
% residual SD of 0.48 kcal/mol, which is too loose to translate a limit from one to the other.
%
% Rychlik's Table 1 rounds TA to -1.0 where Breslauer's original has -0.9, so TATAC reads -4.8 here
% against the 4.6 Primer3 quotes. The table is used as printed in Rychlik, since that is the source
% of the -9 limit.
% ____________________________________________________________________________________________________________
% Usage
%   dG = RychlikEndStability(prim, endwin)
%
%   prim     primer sequence, 5' to 3', char or string, letters A C G T
%   endwin   number of 3'-terminal bases in the window, default 5
%
% Output
%   dG       kcal/mol, negative for a stabilising window. NaN when the window holds a letter other
%            than A C G T or fewer than two bases
% ____________________________________________________________________________________________________________
% References
%   Rychlik W. Selection of primers for polymerase chain reaction. Mol Biotechnol 1995;3:129-34.
%       doi:10.1007/BF02789108. Table 1 and section 2.2.1.4.
%   Rychlik W. Selection of primers for polymerase chain reaction. In: White BA, editor. Methods in
%       Molecular Biology, vol 15, PCR Protocols. Humana Press; 1993. p 31-40.
%   Breslauer KJ, Frank R, Blocker H, Marky LA. Predicting DNA duplex stability from the base
%       sequence. Proc Natl Acad Sci USA 1986;83:3746-50.  doi:10.1073/pnas.83.11.3746
%   Primer3 manual, PRIMER_MAX_END_STABILITY. https://primer3.org/manual.html, read 14-11-2026
% ____________________________________________________________________________________________________________

if nargin < 2||isempty(endwin), endwin = 5; end
tab = [-1.9 -1.3 -1.6 -1.5
       -1.9 -3.1 -3.6 -1.6
       -1.6 -3.1 -3.1 -1.3
       -1.0 -1.6 -1.9 -1.9];     % rows are the 5' base and columns the 3' base, order A C G T
prim = upper(char(prim));
n = numel(prim);
w = prim(max(1, n - endwin + 1):n);
dG = NaN;
if numel(w) < 2, return; end
idx = zeros(1, numel(w));
for q = 1:numel(w)
    k = find('ACGT' == w(q), 1);
    if isempty(k), return; end
    idx(q) = k;
end
dG = 0;
for q = 1:(numel(w) - 1)
    dG = dG + tab(idx(q), idx(q + 1));
end
end
