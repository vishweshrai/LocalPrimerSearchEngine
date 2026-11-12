function dGfr = dGHetero(dGfTc, dGrT, f0, r0, phi, template, Ta)
%% DGHETERO  Threshold free energy for a forward-reverse cross-dimer, as a fraction of the
%% limiting primer
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 26-10-2026 18:24 UTC
% ____________________________________________________________________________________________________________
% WHAT phi MEANS HERE
%
%       phi = [fr] / [r]0
%
% the fraction of the LIMITING primer tied up in cross-dimer. The limiting primer is the one that
% matters: it sets how much duplex forms before the reaction switches to linear synthesis, so
% primer lost from it is product not made. The excess primer loses the same absolute amount and
% barely notices.
%
% This is a different denominator from the mispriming stages, where phi is the ratio of unwanted to
% wanted complex at one site. A misprime competes for TEMPLATE and is rightly judged against the
% primer-template duplex. A dimer competes for PRIMER. Judging a dimer against the duplex pool
% makes the gate tighten as template is diluted, which is backwards: removing a species that is not
% a reactant cannot make the dimer more harmful.
%
% ____________________________________________________________________________________________________________
% THE ALGEBRA
%
% Species: free forward f, free reverse r, forward-template fT*, reverse-template rT, dimer fr.
%       [f]0 = [f] + [fT*] + [fr],      [fT*] = Kf[f][T*]0 / (1 + Kf[f])
%       [r]0 = [r] + [rT]  + [fr],      [rT]  = Kr[r][T]0  / (1 + Kr[r])
%       [fr] = Kh[f][r]
%
% Pinning [fr] = phi*[r]0 DECOUPLES the two balances, because the dimer term is then a known
% constant in each. Each becomes one quadratic in one unknown:
%
%       Kr[r]^2 + ( 1 - Kr( [r]0(1-phi) - [T]0  ) )[r] - [r]0(1-phi)      = 0
%       Kf[f]^2 + ( 1 - Kf( [f]0 - phi[r]0 - [T*]0 ) )[f] - ( [f]0 - phi[r]0 ) = 0
%
% and then
%       Kh = phi*[r]0 / ([f][r]),       dGfr = -R*Ta*ln(Kh)
%
% This is simpler than the coupled pair the previous reference required, where phi multiplied the
% reverse-template term inside the forward balance and the two had to be solved in sequence.
%
% THE LIMIT IS WORTH KNOWING. With primer in large excess over template, [f] -> [f]0 and
% [r] -> [r]0, and the limiting primer CANCELS:
%
%       Kh -> phi / [f]0,        dGfr -> -R*Ta*ln( phi / [f]0 )
%
% The threshold is set by the EXCESS primer concentration alone: not by template, not by the
% limiting primer, not by either primer-template free energy. That is the right dependence, because
% the fraction of limiting primer consumed depends on how much excess primer there is to collide
% with, and [r]0 scales the budget and the collision rate equally. Checked against the exact
% solution at 58 C with 1 uM forward and 50 nM reverse: -4.549 against -4.545 at phi = 1E-3, and
% -6.071 against -6.061 at phi = 1E-2.
%
% WHAT THIS GATE DOES NOT COVER. Sequestration only. Extensibility is tested separately against
% dGEndMin, on a quantity that depends on whether the polymerase carries a 3' -> 5' exonuclease.
%
% ____________________________________________________________________________________________________________
% INPUTS   dGfTc, dGrT   primer-template duplex free energies at Ta, kcal/mol
%          f0, r0        total excess and limiting primer, M
%          phi           tolerated fraction of the LIMITING primer held as cross-dimer
%          template      total template, M
%          Ta            temperature, degrees C
% OUTPUT   dGfr          kcal/mol; reject when the measured cross-dimer dG falls below it
%
% ____________________________________________________________________________________________________________
% REFERENCES
%   Pierce KE, Sanchez JA, Rice JE, Wangh LJ. LATE-PCR primer design criteria.
%       Proc Natl Acad Sci USA 2005;102:8609-14.  doi:10.1073/pnas.0501946102
%   SantaLucia J Jr, Hicks D. The thermodynamics of DNA structural motifs.
%       Annu Rev Biophys Biomol Struct 2004;33:415-40.  doi:10.1146/annurev.biophys.32.110601.141800
% ____________________________________________________________________________________________________________

R  = 8.314E-3/4.184;
TK = Ta + 273.15;
Kf = exp(-dGfTc./(R.*TK));
Kr = exp(-dGrT ./(R.*TK));

if any(phi(:) <= 0)||any(phi(:) >= 1)
    error('dGHetero:BadPhi', 'phi is a fraction of the limiting primer and must lie strictly in (0,1).');
end
if any((f0(:) - phi(:).*r0(:)) <= 0)
    error('dGHetero:PhiExceedsExcess', ...
        ['phi*[r]0 = %.3e M exceeds the excess primer at %.3e M, so the dimer would consume more ' ...
        'forward primer than exists. Lower phi or raise the excess primer.'], phi.*r0, f0);
end

rsol = StableRoot(Kr, 1 - Kr.*(r0.*(1 - phi) - template), -r0.*(1 - phi), 'reverse primer');
fsol = StableRoot(Kf, 1 - Kf.*(f0 - phi.*r0 - template), -(f0 - phi.*r0), 'forward primer');

Kh   = phi.*r0./(fsol.*rsol);
dGfr = -R.*TK.*log(Kh);
end

%% ____________________________________________________________________________________________________________
function x = StableRoot(A, B, C, label)
% Positive root of A*x^2 + B*x + C = 0 for A > 0 and C < 0. The textbook
% (-B + sqrt(B^2-4AC))/(2A) loses precision to cancellation whenever B > 0 and 4AC is small against
% B^2, which is the regime here whenever the association constant is modest.
disc = B.^2 - 4.*A.*C;
if any(disc(:) < 0)
    error('dGHetero:ComplexRoot', ...
        'Negative discriminant solving the %s mass balance.', label);
end
signB = ones(size(B)); signB(B < 0) = -1;
q = -0.5.*(B + signB.*sqrt(disc));
x = max(q./A, C./q);
if any(x(:) <= 0)
    error('dGHetero:NonPositiveRoot', 'Non-positive free %s concentration.', label);
end
end
