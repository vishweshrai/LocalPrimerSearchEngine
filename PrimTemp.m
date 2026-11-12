function [tm, h, s] = PrimTemp(seq, exstrand, limstrand, monv, biv, dNTP, N)
%% PRIMTEMP  Melting temperature of a fully paired primer-template duplex
% Written by Vishwesh Rai, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 30-10-2026 15:48 IST
% ____________________________________________________________________________________________________________
% WHAT THIS FUNCTION COMPUTES
%
% The temperature at which half the limiting strand of a bimolecular duplex is melted, from the
% unified nearest-neighbour model. For two NON-IDENTICAL strands at total concentrations CT and CL
% with CT the strand in excess,
%
%       Tm = dH / ( dS + R * ln(CT - CL/2) )
%
% where dH sums the sixteen dinucleotide steps plus one helix initiation at each end, keyed on
% whether that end is a G.C or an A.T pair. The enthalpy and entropy returned alongside Tm are the
% 1 M NaCl values, before any salt correction.
%
% EVERY PARAMETER COMES FROM NNParams.m
%
% The nearest-neighbour table, the two initiation rows, the magnesium-dNTP speciation constant and
% the Owczarzy coefficients are all read from the shared module. Nothing is written down here. A
% table duplicated across files is free to drift, and a drifting copy inside this function is worse
% than most: the Tm it returns sets the annealing temperature at which every other gate in the
% suite is then evaluated, so a discrepancy here moves the whole screen without leaving a trace.
%
% MAGNESIUM SPECIATION
%
% dNTPs chelate Mg2+ with association constant Ka = 3E4 /M, so the magnesium available to screen the
% backbone is the positive root of
%
%       Ka*m^2 + m*(1 + Ka*[dNTP] - Ka*[Mg]tot) - [Mg]tot = 0
%
% rather than the difference [Mg]tot - [dNTP]. The two agree within about one percent at ordinary
% PCR compositions and diverge once dNTP approaches magnesium. N.Salt.FreeMg solves it.
%
% SALT CORRECTION
%
% N.Salt.Owczarzy returns the entropy shift equivalent to the Owczarzy correction on 1/Tm. The two
% forms are algebraically the same quantity: the correction adds f*1E-5 to 1/Tm, and since
% Tm = dH/(dS + R ln CT'), that is the same as adding dH*f*1E-5 to dS. The monovalent branch is used
% when sqrt([Mg]free)/[Na+] falls below 0.22 and the divalent branch above it, which is the
% switching rule of the 2008 paper.
%
% NO SYMMETRY CORRECTION IS APPLIED
%
% The -1.4 cal/(K mol) symmetry entropy of Equation 3 exists because two IDENTICAL strands forming a
% homodimer are rotationally indistinguishable. A primer annealing to a plasmid template is a
% heteroduplex of two distinct species, the concentration term above already uses the non-identical
% strand form, and charging the symmetry term as well would penalise a self-complementary primer by
% about 0.46 kcal/mol at 58 C for binding its own target.
%
% REFERENCES
%   SantaLucia J Jr, Hicks D. The thermodynamics of DNA structural motifs.
%       Annu Rev Biophys Biomol Struct 2004;33:415-40. doi:10.1146/annurev.biophys.32.110601.141800
%   Owczarzy R, Moreira BG, You Y, Behlke MA, Walder JA. Predicting stability of DNA duplexes in
%       solutions containing magnesium and monovalent cations.
%       Biochemistry 2008;47:5336-53. doi:10.1021/bi702363u
%   von Ahsen N, Wittwer CT, Schutz E. Oligonucleotide melting temperatures under PCR conditions.
%       Clin Chem 2001;47:1956-61. doi:10.1093/clinchem/47.11.1956
% ____________________________________________________________________________________________________________
% INPUTS   seq        char row, uppercase A C G T
%          exstrand   total concentration of the strand in excess, M
%          limstrand  total concentration of the limiting strand, M
%          monv       monovalent cation concentration, M
%          biv        TOTAL divalent cation concentration, M
%          dNTP       total dNTP concentration, M
%          N          optional NNParams struct; built and cached here when omitted
% OUTPUTS  tm         degrees C, NaN when the concentrations do not define a bimolecular melt
%          h          duplex enthalpy at 1 M NaCl, kcal/mol
%          s          duplex entropy at 1 M NaCl, kcal/(K mol)
% ____________________________________________________________________________________________________________

persistent Ncache
if (nargin < 7)||isempty(N)
    if isempty(Ncache), Ncache = NNParams(); end
    N = Ncache;
end

seq = upper(char(seq));
L = numel(seq);
if L < 2
    error('PrimTemp:SegmentTooShort', ...
        'A duplex of %d base(s) carries no nearest-neighbour step.', L);
end
if any(~ismember(seq, 'ACGT'))
    error('PrimTemp:BadBase', 'seq must contain only A, C, G and T.');
end

%% ---- nearest-neighbour sum at 1 M NaCl ------------------------------------------------------------------
h = 0; s = 0;
for i = 1:numel(N.WC.Step)
    nn = numel(strfind(seq, N.WC.Step{i}));
    if nn == 0, continue; end
    h = h + nn*N.WC.dH(i);
    s = s + 0.001*nn*N.WC.dS(i);
end

%% ---- one helix initiation at each end -------------------------------------------------------------------
for b = [seq(1), seq(L)]
    if (b == 'G')||(b == 'C')
        h = h + N.WC.InitGC(1); s = s + 0.001*N.WC.InitGC(2);
    else
        h = h + N.WC.InitAT(1); s = s + 0.001*N.WC.InitAT(2);
    end
end

%% ---- concentration term ----------------------------------------------------------------------------------
if exstrand < limstrand
    [exstrand, limstrand] = deal(limstrand, exstrand);
elseif limstrand == 0
    tm = NaN; return
end
R  = N.Rgas;
tm = h/(s + R*log(exstrand - limstrand/2));

%% ---- salt --------------------------------------------------------------------------------------------------
gc    = nnz((seq == 'G')|(seq == 'C'))/L;
dsAdd = N.Salt.Owczarzy(h, gc, L, monv, biv, dNTP);     % cal/(K mol)
tm    = h/(s + 0.001*dsAdd + R*log(exstrand - limstrand/2));
tm    = tm - 273.15;
end
