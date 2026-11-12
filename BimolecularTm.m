function [tm, dH, dS] = BimolecularTm(g1, g2, T1, T2, ceff, R)
%% BIMOLECULARTM  Melting temperature of one bimolecular complex from its free energy at two points
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 06-11-2026 13:20 IST
% ____________________________________________________________________________________________________________
% PHYSICS
%
% A unimolecular structure melts where dG = 0, because its folded fraction is K/(1 + K) and K = 1
% there whatever the concentration. An association has no melting temperature until a concentration
% is named, so the condition carries a concentration term and the answer moves when the
% concentration moves.
%
% Write the complex as A + B <-> AB and let ceff be the concentration group that the equilibrium
% constant equals the reciprocal of at half conversion. Then
%
%     K(Tm) = 1/ceff      and      K = exp(-dG/RT)
%
%     dG(Tm) = R*Tm*ln(ceff)      ->      Tm = dH/(dS + R*ln(ceff))
%
% WHAT TO PASS AS ceff. The group depends on the stoichiometry and on which species the half point
% is referred to.
%
%   homodimer, A + A <-> A2                     ceff = CT, the total strand concentration
%     At Tm half the strands are paired, so [A2] = CT/4 and [A] = CT/2, giving K = 1/CT.
%
%   heterodimer, a0 >= b0, half of B bound       ceff = a0 - b0/2
%     [AB] = b0/2, [B] = b0/2, [A] = a0 - b0/2, giving K = 1/(a0 - b0/2). At a0 = b0 = CT/2 this
%     collapses to CT/4, the familiar non-self-complementary form.
%
%   primer over a template site, primer in       ceff = f, the FREE primer concentration
%   large excess, half the sites occupied
%     Occupancy of a site is theta = K*f/(1 + K*f), so theta = 1/2 is K*f = 1 and K = 1/f. This is
%     the a0 >> b0 limit of the line above, where a0 - b0/2 tends to a0.
%
% HOW dH AND dS ARE RECOVERED. dG(T) is linear in T for a FIXED structure, so two evaluations give
% the pair exactly:
%
%     dS = (g1 - g2)/(T2 - T1)    in kelvin      dH = g1 + T1*dS
%
% The caller is responsible for the two free energies describing the SAME structure. Scoring a
% different register or a different sub-duplex at the two temperatures puts a straight line through
% two unrelated points and the number that comes back is meaningless. DimerMeltTm compares the
% winning register between its probes for exactly this reason.
%
% RANGE. For a well-chosen primer the dimer answers land tens of degrees below zero, which is a
% linear extrapolation of dG(T) far outside the range the nearest-neighbour parameters were fitted
% over. Such a value orders designs and states distance from the annealing step. It is not a
% temperature anyone could measure.
% ____________________________________________________________________________________________________________
% INPUTS
%   g1     free energy of the structure at T1, kcal/mol. An array scores many structures at once,
%          one per element, which is how PrimerHybridizationMap converts a whole site list in one
%          call. Non-finite entries, such as the Inf marking a site with no anchored duplex, come
%          back as NaN rather than raising
%   g2     free energy of the SAME structure at T2, kcal/mol. Same size as g1
%   T1     first temperature, degrees C
%   T2     second temperature, degrees C. Must differ from T1
%   ceff   concentration group entering the logarithm, M. See WHAT TO PASS AS ceff above
%   R      gas constant, kcal/(mol K). Default 1.98720425E-3. Pass the same constant the caller
%          used to produce g1 and g2: dGMismatch and the threshold functions carry
%          8.314E-3/4.184, which differs from the default in the sixth figure
%
% OUTPUTS
%   tm     melting temperature, degrees C, the size of g1. NaN where the recovered entropy plus
%          the concentration term is too close to zero to divide by, and NaN where either input
%          free energy was not finite
%   dH     recovered enthalpy, kcal/mol
%   dS     recovered entropy, kcal/mol/K
% ____________________________________________________________________________________________________________
% REFERENCES
%
% SantaLucia J Jr, Hicks D (2004) The thermodynamics of DNA structural motifs. Annu Rev Biophys
%   Biomol Struct 33:415-440. doi:10.1146/annurev.biophys.32.110601.141800. Source of the
%   nearest-neighbour parameters and of the concentration term in the melting relation.
% ____________________________________________________________________________________________________________

if (nargin < 6)||isempty(R), R = 1.98720425E-3; end

if ~isnumeric(g1)||~isnumeric(g2)
    error('BimolecularTm:BadEnergy', 'g1 and g2 must be numeric, in kcal/mol.');
end
if ~isequal(size(g1), size(g2))
    error('BimolecularTm:SizeMismatch', ...
        'g1 and g2 must be the same size; each pair describes one structure.');
end
if ~isscalar(T1)||~isscalar(T2)||~isfinite(T1)||~isfinite(T2)
    error('BimolecularTm:BadTemperature', 'T1 and T2 must be finite scalars in degrees C.');
end
if T1 == T2
    error('BimolecularTm:EqualTemperatures', ...
        'T1 and T2 must differ; two points at one temperature fix no slope.');
end
if ~isscalar(ceff)||~isfinite(ceff)||(ceff <= 0)
    error('BimolecularTm:BadConcentration', 'ceff must be a finite positive scalar in M.');
end
if ~isscalar(R)||~isfinite(R)||(R <= 0)
    error('BimolecularTm:BadGasConstant', 'R must be a finite positive scalar in kcal/(mol K).');
end

K1 = T1 + 273.15;
K2 = T2 + 273.15;
dS = (g1 - g2)/(K2 - K1);           % kcal/mol/K
dH = g1 + K1*dS;                    % kcal/mol

den = dS + R*log(ceff);
tm  = dH./den - 273.15;
tm(~isfinite(den)|(abs(den) < 1E-9)) = NaN;
tm(~isfinite(tm)) = NaN;
end
