function [tm, info] = DimerMeltTm(seqa, seqb, ca, cb, monv, biv, dNTP, tann, dT)
%% DIMERMELTTM  Melting temperature of the most stable dimer register, at the working concentrations
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 06-11-2026 13:20 IST
% ____________________________________________________________________________________________________________
% PHYSICS
%
% A dimer is bimolecular, so unlike a hairpin it has no melting temperature until a concentration
% is named. The hairpin case solves dG(Tm) = 0 because the folded fraction of a unimolecular
% species is K/(1 + K) and K = 1 at dG = 0, which is half folded whatever the concentration. That
% argument does not carry over. For an association the bound fraction depends on how much material
% is present, so the melting condition carries a concentration term and the answer moves when the
% primer concentration moves.
%
% HOMODIMER.  A + A <-> A2, with CT the total strand concentration. At the melting temperature half
% the strands sit in duplex, so [A2] = CT/4 and [A] = CT/2, giving
%
%     K(Tm) = [A2]/[A]^2 = (CT/4)/(CT/2)^2 = 1/CT
%
% and since K = exp(-dG/RT),
%
%     dG(Tm) = R*Tm*ln(CT)        ->      Tm = dH/(dS + R*ln(CT))
%
% HETERODIMER.  A + B <-> AB with a0 >= b0. The melting point is referred to the LIMITING strand,
% because that is the one whose titration the temperature sweep reports. At Tm half of B is bound,
% so [AB] = b0/2, [B] = b0/2 and [A] = a0 - b0/2, giving
%
%     K(Tm) = [AB]/([A][B]) = 1/(a0 - b0/2)
%
%     Tm = dH/(dS + R*ln(a0 - b0/2))
%
% At a0 = b0 = CT/2 the group a0 - b0/2 collapses to CT/4, which is the familiar non-self-
% complementary form, so the two expressions above are one formula written for two stoichiometries
% rather than two separate rules.
%
% WHY THIS MATTERS IN A LATE-PCR DESIGN. The excess and limiting primers are deliberately held at
% different concentrations, commonly twenty-fold apart. Through the logarithm that difference is
% worth R*ln(20) = 5.95 cal/mol/K of effective entropy, which for a dimer of dH near -60 kcal/mol
% displaces the melting temperature by roughly
%
%     dTm  ~  Tm^2 * R * ln(20) / abs(dH)  ~  11 K
%
% So the homodimer of the excess primer and the homodimer of the limiting primer cannot be read
% against one temperature threshold. Each is reported here at its own concentration.
%
% HOW dH AND dS ARE OBTAINED. dG(T) is linear in T for a FIXED register, so two evaluations recover
% the pair exactly. The relation itself lives in BimolecularTm, which both this function and
% PrimerHybridizationMap call, so there is one implementation of it rather than two that can drift.
% DimerScan is called at tann and at tann + dT and the winning register is
% compared between them. A register that changes across the interval breaks the linearity the fit
% assumes, and in that case the function returns NaN rather than a number obtained by fitting a
% straight line through two different structures.
%
% RANGE. A well-chosen primer gives a dimer melting temperature far below anything the thermal
% programme visits, often tens of degrees below zero. Such a value is a linear extrapolation of
% dG(T) well outside the range the nearest-neighbour parameters were fitted over, so it should be
% read as headroom on an ordinal scale and not as a temperature anyone could measure. What it
% states reliably is the distance from the annealing step, and that distance is the point.
%
% The returned temperature describes the register DimerScan judges most stable at tann. It is a
% diagnostic of how much headroom the design has, not a gate. Dimer admission is decided on free
% energy against a tolerance that already carries the concentration, which is the same physics
% expressed in the currency the gate needs.
% ____________________________________________________________________________________________________________
% INPUTS
%   seqa   first sequence, char. For a homodimer pass the same sequence twice
%   seqb   second sequence, char
%   ca     working concentration of seqa, M
%   cb     working concentration of seqb, M. Ignored for a homodimer, where CT = ca
%   monv   total monovalent cation, M
%   biv    total divalent cation, M
%   dNTP   total dNTP, M
%   tann   temperature at which the register is chosen, degrees C. The engine passes Tm - dTsub,
%          the same probe DimerScan is given elsewhere, so the register reported here is the one
%          the free-energy screen acted on
%   dT     offset of the second probe, degrees C. Default 5
%
% OUTPUTS
%   tm     melting temperature of that register, degrees C. NaN when no register carries the
%          minimum run DimerScan requires, when the register moves between the two probes, or when
%          the recovered entropy is too small to divide by
%   info   struct with
%            .Kind            'homodimer' or 'heterodimer'
%            .CEff            the concentration group entering the logarithm, M
%            .dH              recovered enthalpy, kcal/mol
%            .dS              recovered entropy, kcal/mol/K
%            .dGAtProbe       dG of the register at tann, kcal/mol
%            .RegisterStable  false when the winning register moved between the probes
%            .Align           3-row alignment block for the register at tann
% ____________________________________________________________________________________________________________
% REFERENCES
%
% SantaLucia J Jr, Hicks D (2004) The thermodynamics of DNA structural motifs. Annu Rev Biophys
%   Biomol Struct 33:415-440. doi:10.1146/annurev.biophys.32.110601.141800. Source of the
%   nearest-neighbour parameters DimerScan sums and of the concentration term in the melting
%   relation used above.
% ____________________________________________________________________________________________________________

RGAS = 1.98720425E-3;               % kcal/mol/K

if (nargin < 9)||isempty(dT), dT = 5; end

if ~ischar(seqa)&&~isstring(seqa), error('DimerMeltTm:NonChar', 'seqa must be char.'); end
if ~ischar(seqb)&&~isstring(seqb), error('DimerMeltTm:NonChar', 'seqb must be char.'); end
seqa = upper(char(seqa));
seqb = upper(char(seqb));
if isempty(seqa)||isempty(seqb)
    error('DimerMeltTm:EmptySequence', 'Both sequences must be non-empty.');
end
if ~isscalar(ca)||~isfinite(ca)||(ca <= 0)
    error('DimerMeltTm:BadConcentration', 'ca must be a finite positive scalar in M.');
end
if ~isscalar(dT)||~isfinite(dT)||(dT == 0)
    error('DimerMeltTm:BadOffset', 'dT must be a finite non-zero scalar in degrees C.');
end

isHomo = strcmp(seqa, seqb);
if isHomo
    kind = 'homodimer';
    ceff = ca;
else
    if ~isscalar(cb)||~isfinite(cb)||(cb <= 0)
        error('DimerMeltTm:BadConcentration', 'cb must be a finite positive scalar in M.');
    end
    kind = 'heterodimer';
    hi   = max(ca, cb);
    lo   = min(ca, cb);
    ceff = hi - lo/2;
end
if ceff <= 0
    error('DimerMeltTm:BadConcentration', ...
        'The concentration group entering the logarithm came out non-positive.');
end

info = struct('Kind', kind, 'CEff', ceff, 'dH', NaN, 'dS', NaN, ...
    'dGAtProbe', NaN, 'RegisterStable', false, 'Align', '');
tm = NaN;

[g1, ~, align1, ~] = DimerScan(seqa, seqb, monv, biv, dNTP, tann);
info.dGAtProbe = g1;
info.Align     = align1;
if g1 == 0
    % DimerScan reports exactly zero when no register carries its minimum run, which is the
    % absence of a structure rather than a structure of zero free energy.
    return;
end

[g2, ~, align2, ~] = DimerScan(seqa, seqb, monv, biv, dNTP, tann + dT);
if g2 == 0, return; end

info.RegisterStable = isequal(align1, align2);
if ~info.RegisterStable
    % Two different structures cannot be joined by one straight line in T.
    return;
end

[tm, dH, dS] = BimolecularTm(g1, g2, tann, tann + dT, ceff, RGAS);
info.dH = dH;
info.dS = dS;
end
