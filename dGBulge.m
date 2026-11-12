function [dG, info] = dGBulge(bulgeLen, pair5, pair3, tann, N)
%% DGBULGE  Free energy of a bulge loop, one strand carrying unpaired bases the other does not
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 26-10-2026 11:55 UTC
% ____________________________________________________________________________________________________________
% PHYSICS
%
% A bulge is an insertion on one strand with nothing opposite it. The helix must accommodate the
% extra nucleotide or nucleotides without losing register, which costs entropy and, for longer
% bulges, breaks the stack across the insertion.
%
% THE SINGLE-BASE BULGE IS A SPECIAL CASE AND THE REASON MATTERS
%
% For a bulge of one, the bulged base flips out and the two flanking base pairs stack on each other
% as though it were absent. The nearest-neighbour term between those flanking pairs is therefore
% RETAINED and added to the bulge initiation. For a bulge of two or more the strand cannot bridge
% the gap that way, the stack is lost, and only the initiation and the closure penalties remain.
% This is why bulge initiation is NOT monotonic in length: a bulge of 1 costs +4.0 kcal/mol and a
% bulge of 2 costs +2.9, which looks wrong until the retained stack of roughly -1 to -2 kcal/mol is
% restored to the first and makes it the more stable of the two.
%
% Following the review, the calculation is
%
%   bulge of 1      dG = dG(bulge of 1) + dG(intervening nearest neighbour) + A-T closure penalties
%   bulge of n > 1  dG = dG(bulge of n) + A-T closure penalties
%
% The review works one example by hand: 5'-CAT-3' over 5'-AG-3', bulging the A, gives
% +4.0 - 1.28 + 0.5 = +3.22 kcal/mol. NNParams asserts that this file's tables reproduce it.
%
% TEMPERATURE
%
% Every loop enthalpy in the database is zero by construction, so a bulge penalty is purely
% entropic and scales as dG(T) = dG37 * T / 310.15. It grows with temperature rather than staying
% fixed, which is the opposite of the intuition that a destabilising term matters less when hot.
% The retained nearest-neighbour stack of a single-base bulge carries its own enthalpy and is scaled
% properly from dH and dS.
%
% WHAT IS NOT MODELLED
%
% Sequence dependence of the bulged base itself. The review reports that A bulges are appreciably
% more stable than C bulges with G and T between, attributing it to unpublished work, and gives no
% parameters. Bulge lengths beyond 30 use Jacobson-Stockmayer.
%
% ____________________________________________________________________________________________________________
% INPUTS
%   bulgeLen  number of unpaired bases on the bulged strand, 1 or more
%   pair5     the base pair closing the bulge on its 5' side, two characters, top then bottom
%   pair3     the base pair closing it on the 3' side, two characters, top then bottom
%   tann      temperature, degrees C
%   N         the NNParams struct, fetched if omitted
%
% OUTPUTS
%   dG        kcal/mol at tann, 1 M NaCl
%   info      struct: Init, Stack, Closure, Extrapolated, Phosphates
%
% ____________________________________________________________________________________________________________
% REFERENCES
%   SantaLucia J Jr, Hicks D. The thermodynamics of DNA structural motifs. Table 4 and Equation 13;
%       Jacobson-Stockmayer as Equation 7.
%       Annu Rev Biophys Biomol Struct 2004;33:415-40.  doi:10.1146/annurev.biophys.32.110601.141800
%   Turner DH. Bulges in nucleic acids. Curr Opin Struct Biol 1992;2:334-37.
%       doi:10.1016/0959-440X(92)90222-S
% ____________________________________________________________________________________________________________

if (nargin < 5)||isempty(N), N = NNParams(); end
info = struct('Init', NaN, 'Stack', 0, 'Closure', 0, 'Extrapolated', false, 'Phosphates', 0);
dG = NaN;

if ~isfinite(bulgeLen)||(bulgeLen < 1)||(bulgeLen ~= fix(bulgeLen))
    error('dGBulge:BadLength', 'bulgeLen must be a positive whole number, not %s.', mat2str(bulgeLen));
end
if ~isfield(N, 'Loop')||isempty(N.Loop.N)
    return
end

TK = tann + 273.15;

%% ---- initiation, interpolated inside the table and Jacobson-Stockmayer beyond it --------------
nTab = N.Loop.N; gTab = N.Loop.Bulge;
ok = isfinite(gTab); nTab = nTab(ok); gTab = gTab(ok);
if bulgeLen <= max(nTab)
    g37 = interp1(nTab, gTab, bulgeLen, 'linear');
    info.Extrapolated = ~any(nTab == bulgeLen);
else
    % dG37(loop n) = dG37(loop x) + 2.44 * R * 310.15 * ln(n/x), x the longest measured length
    g37 = gTab(end) + N.Loop.JSCoeff*N.Rgas*310.15*log(bulgeLen/nTab(end));
    info.Extrapolated = true;
end
% Loop enthalpy is zero, so the whole initiation is entropic and scales linearly in temperature.
info.Init = g37*TK/310.15;

%% ---- retained nearest-neighbour stack, single-base bulge only ---------------------------------
if bulgeLen == 1
    if (numel(pair5) < 2)||(numel(pair3) < 2)
        error('dGBulge:BadPair', ...
            'A bulge of 1 needs both flanking pairs as two characters each to form the retained stack.');
    end
    step = [upper(pair5(1)) upper(pair3(1))];
    idx = find(strcmp(N.WC.Step, step), 1);
    if isempty(idx)
        error('dGBulge:NoStep', ...
            'The intervening step ''%s'' is not in the Watson-Crick table.', step);
    end
    info.Stack = N.WC.dH(idx) - TK*N.WC.dS(idx)/1000;
    info.Phosphates = 1;
end

%% ---- A-T closure penalty, once per closing pair that is A-T ------------------------------------
% Applied to each closing pair independently. The review's worked example carries one A-T closure
% and one G-C closure and takes the penalty once.
info.Closure = N.Loop.ATClosureG*(IsAT(pair5) + IsAT(pair3));

dG = info.Init + info.Stack + info.Closure;
end

%% ____________________________________________________________________________________________________________
function t = IsAT(pr)
% True when the closing pair is A-T in either orientation.
if numel(pr) < 2
    t = false; return
end
p = upper(pr(1:2));
t = strcmp(p, 'AT')||strcmp(p, 'TA');
end
