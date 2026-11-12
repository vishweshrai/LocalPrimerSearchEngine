function dGhomo = dGHomo(dGprim, phi, p0, target, Ta)
%% DGHOMO  Threshold free energy for a primer self-dimer, as a fraction of the primer pool
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 29-10-2026 07:05 UTC
% ____________________________________________________________________________________________________________
% WHAT phi MEANS HERE, AND WHY IT IS NOT WHAT IT MEANS AT THE MISPRIMING STAGES
%
%       phi = 2*[p2] / [p]0
%
% the fraction of the PRIMER POOL tied up in self-dimer. The factor of two is there because a
% homodimer consumes two primer molecules.
%
% This is a different denominator from dGMispTarget and dGMispOffTarget, where phi is the ratio of
% unwanted to wanted complex AT ONE SITE, and the difference is physical rather than conventional.
% A mispriming event competes for TEMPLATE: it occupies a site that the real primer wanted, so the
% quantity that matters is spurious product against intended product and the primer-template duplex
% is the right reference. A dimer competes for PRIMER: what it costs the reaction is primer removed
% from circulation, and that has nothing to do with how much template was loaded.
%
% THE OLD REFERENCE COUPLED THIS GATE TO TEMPLATE LOADING, WHICH IS THE WRONG DEPENDENCE.
% Referenced to the primer-template duplex, the threshold moved by 4.5 kcal/mol across three
% decades of template concentration. Diluting the template cannot make a primer dimer more harmful,
% because the template is not a reactant in it. Under the present reference the dependence is
% second order and the threshold is template-independent IN THE PRIMER-EXCESS LIMIT, which is the
% regime this suite runs in. It is not independent identically: [T]0 enters the quadratic through
% B = 1 - Kp*([p]0*(1 - phi) - target), so the free primer it solves for, and therefore Kh, do
% move with template. What the change bought is the size of that movement. Referenced to the
% primer-template duplex the threshold swung 4.5 kcal/mol over three decades of template;
% referenced to the primer pool it moves by less than 0.001 kcal/mol over the same range.
%
% ____________________________________________________________________________________________________________
% THE ALGEBRA
%
% Species: free primer p, primer-template duplex pT*, self-dimer p2.
%       [p]0  = [p] + [pT*] + 2[p2]
%       [pT*] = Kp[p][T]0 / (1 + Kp[p])
%       [p2]  = Kh[p]^2
% Pinning 2[p2] = phi*[p]0 leaves a single quadratic in the free primer,
%
%       Kp[p]^2 + ( 1 - Kp( [p]0(1-phi) - [T]0 ) )[p] - [p]0(1-phi) = 0
%
% and the dimer constant that satisfies the constraint follows,
%
%       Kh = phi*[p]0 / (2[p]^2),        dGhomo = -R*Ta*ln(Kh)
%
% With primer in large excess over template the free primer tends to [p]0 and this reduces to
%
%       dGhomo -> -R*Ta*ln( phi / (2[p]0) )
%
% which says the threshold is set by the primer concentration alone. That is the useful reading: a
% dilute primer meets itself less often, so it takes a MORE stable dimer to sequester the same
% fraction of it, and its threshold is correspondingly more permissive. At 1 uM and 50 nM with
% phi = 1E-2 the two thresholds are -5.60 and -7.58 kcal/mol at 58 C.
%
% WHAT THIS GATE DOES NOT COVER. Sequestration only. Whether a dimer can be EXTENDED, which is the
% failure that actually ruins a reaction because an extended dimer amplifies as a competing
% template, is tested separately against dGEndMin in LocalPrimerSearchEngine.m, on a quantity that
% depends on whether the polymerase proofreads.
%
% ____________________________________________________________________________________________________________
% INPUTS   dGprim  primer-template duplex free energy at Ta, kcal/mol
%          phi     tolerated fraction of the primer pool held as dimer
%          p0      total primer, M
%          target  total template, M
%          Ta      temperature, degrees C
% OUTPUT   dGhomo  kcal/mol. A candidate is rejected when its self-dimer dG falls BELOW this value,
%                  since more negative means more stable means more dimer.
%
% ____________________________________________________________________________________________________________
% REFERENCES
%   SantaLucia J Jr, Hicks D. The thermodynamics of DNA structural motifs.
%       Annu Rev Biophys Biomol Struct 2004;33:415-40.  doi:10.1146/annurev.biophys.32.110601.141800
% ____________________________________________________________________________________________________________

R  = 8.314E-3/4.184;
TK = Ta + 273.15;
Kp = exp(-dGprim./(R.*TK));

if any(phi(:) <= 0)||any(phi(:) >= 1)
    error('dGHomo:BadPhi', 'phi is a fraction of the primer pool and must lie strictly in (0,1).');
end

A = Kp;
B = 1 - Kp.*(p0.*(1 - phi) - target);
C = -p0.*(1 - phi);
p = StableRoot(A, B, C);

Kh = phi.*p0./(2.*p.^2);
dGhomo = -R.*TK.*log(Kh);
end

%% ____________________________________________________________________________________________________________
function x = StableRoot(A, B, C)
% Positive root of A*x^2 + B*x + C = 0 with A > 0 and C < 0, in the form that does not lose
% precision to cancellation when B > 0 and 4AC is small against B^2.
disc = B.^2 - 4.*A.*C;
if any(disc(:) < 0)
    error('dGHomo:ComplexRoot', 'Negative discriminant solving the primer mass balance.');
end
signB = ones(size(B)); signB(B < 0) = -1;
q = -0.5.*(B + signB.*sqrt(disc));
x = max(q./A, C./q);
if any(x(:) <= 0)
    error('dGHomo:NonPositiveRoot', 'Non-positive free primer concentration from the mass balance.');
end
end
