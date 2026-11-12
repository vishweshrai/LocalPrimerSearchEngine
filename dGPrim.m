function [dG, info] = dGPrim(seq, monv, biv, dNTP, tann, ctx)
%% DGPRIM  Watson-Crick duplex free energy on the shared parameter module
% Written by Vishwesh Rai, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 30-10-2026 15:48 IST
% ____________________________________________________________________________________________________________
% WHAT THIS FUNCTION COMPUTES
%
% The standard free energy of a fully Watson-Crick paired DNA duplex at an arbitrary temperature and
% salt, by the unified nearest-neighbour model. The sum runs over the sixteen dinucleotide steps,
% plus one helix initiation at each end keyed on whether that end is a G.C or an A.T pair:
%
%       dG(T) = sum_steps dH - T * sum_steps dS + dG(init 5') + dG(init 3')
%
% The nearest-neighbour table, the initiation rows and both salt corrections come from NNParams.m.
% No parameter is written down in this file.
%
% SALT
%
% Two corrections are available and both are computed on every call, so info reports the gap between
% them. ctx.SaltModel selects which one is returned.
%
%   'owczarzy_ref' (default) and 'owczarzy'   The empirical correction of Owczarzy et al., which
%       scales with duplex enthalpy and carries a GC-fraction dependence. On a fully paired duplex
%       the reference enthalpy IS the actual enthalpy, so the two names give the same number here by
%       construction. Both are kept so that this function and dGMismatch can be driven from a single
%       setting, where the distinction does matter.
%
%   'perphosphate'   The counterion-condensation form of SantaLucia 1998. It has no GC dependence,
%       which is a real physical disagreement with the measurement rather than a simplification.
%
% The two differ by about a kcal/mol across the GC range, so a conclusion that survives only one of
% them is a conclusion about the salt model. Run both.
%
% NO SYMMETRY CORRECTION IS APPLIED
%
% The -1.4 cal/(K mol) symmetry entropy of Equation 3 arises from the twofold rotational
% indistinguishability of two IDENTICAL strands associating, A + A -> A2. A primer binding a plasmid
% template is A + B -> AB, where the two species are distinct and the rotational degeneracy does not
% exist. Charging it here would penalise a self-complementary primer for binding its own target by
% about 0.46 kcal/mol at 58 C. The correction belongs in the homodimer path, where two copies of one
% primer genuinely do associate, and nowhere else.
%
% WHY THE STEP COUNTING IS NOT SHARED WITH dGMismatch
%
% This function counts steps with strfind over the sixteen dinucleotide patterns. dGMismatch walks
% positions. The two must return the same number for a fully paired duplex, and keeping the
% implementations independent makes that agreement a test rather than a tautology.
%
% DANGLING ENDS ARE NOT INCLUDED
%
% A primer terminates; the template does not. Every cognate duplex here is scored as a blunt
% oligomer on a continuing strand, which omits roughly 0.5 to 1.0 kcal/mol of stacking at each end.
% dGMismatch scores off-target registers as blunt too, so the omission largely cancels in the margin
% the engine actually gates on. It stops cancelling once registers are no longer trimmed, and it
% needs the strand orientation of the dangle tables settled first.
%
% REFERENCES
%   SantaLucia J Jr, Hicks D. The thermodynamics of DNA structural motifs.
%       Annu Rev Biophys Biomol Struct 2004;33:415-40. doi:10.1146/annurev.biophys.32.110601.141800
%   SantaLucia J Jr. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor
%       thermodynamics. Proc Natl Acad Sci USA 1998;95:1460-5. doi:10.1073/pnas.95.4.1460
%   Owczarzy R, Moreira BG, You Y, Behlke MA, Walder JA. Predicting stability of DNA duplexes in
%       solutions containing magnesium and monovalent cations.
%       Biochemistry 2008;47:5336-53. doi:10.1021/bi702363u
%   von Ahsen N, Wittwer CT, Schutz E. Oligonucleotide melting temperatures under PCR conditions.
%       Clin Chem 2001;47:1956-61. doi:10.1093/clinchem/47.11.1956
% ____________________________________________________________________________________________________________
% INPUTS   seq, monv, biv, dNTP, tann   as before, concentrations in M, temperature in C
%          ctx   optional: .SaltModel ('owczarzy_ref' default, 'owczarzy' or 'perphosphate')
%                          .MonovalentOnly (default false; bare monv rather than the von Ahsen
%                           equivalent in the per-phosphate term)
% OUTPUTS  dG    kcal/mol at tann
%          info  dH, dS, fGC, SaltRoute, dS_salt_used, dS_salt_alternative, dG_salt_model_gap
%
% REFERENCES, none of which I have been able to open; verify before citing
%   SantaLucia J, Hicks D. Annu Rev Biophys Biomol Struct 2004;33:415.
%                                                doi:10.1146/annurev.biophys.32.110601.141800
%   SantaLucia J. PNAS 1998;95:1460.             doi:10.1073/pnas.95.4.1460
%   Owczarzy R et al. Biochemistry 2008;47:5336. doi:10.1021/bi702363u
% ____________________________________________________________________________________________________________

if (nargin < 6)||isempty(ctx), ctx = struct(); end
% Default matches dGMismatch. On a fully paired duplex the reference enthalpy IS the actual
% enthalpy, so 'owczarzy_ref' and 'owczarzy' are the same number here by construction; the two names
% are kept so that both functions can be driven from one setting.
if ~isfield(ctx, 'SaltModel')||isempty(ctx.SaltModel), ctx.SaltModel = 'owczarzy_ref'; end
if ~isfield(ctx, 'MonovalentOnly')||isempty(ctx.MonovalentOnly), ctx.MonovalentOnly = false; end

seq = char(seq);
L = numel(seq);
if L < 2
    error('dGPrim_new:SequenceTooShort', ...
        'A duplex of %d base(s) has no nearest-neighbour step.', L);
end
if any(~ismember(seq, 'ACGT'))
    error('dGPrim_new:NonACGT', 'seq carries a base outside {A,C,G,T}: ''%s''.', seq);
end

N = NNParams(struct('Verbose', false));

%% ---- nearest-neighbour steps, counted by pattern ---------------------------------------------------
h = 0; s = 0; nsteps = 0;
for i = 1:16
    pat = N.WC.Step{i};
    nn = numel(strfind(seq, pat));
    nsteps = nsteps + nn;
    h = h + nn*N.WC.dH(i);
    s = s + 0.001*nn*N.WC.dS(i);
end
if nsteps ~= (L - 1)
    error('dGPrim_new:StepCount', ...
        ['The sixteen dinucleotide patterns matched %d steps in a %d-mer, where %d are expected. ' ...
        'The pattern list in NNParams.WC.Step is no longer a complete cover of the alphabet.'], ...
        nsteps, L, L - 1);
end

%% ---- initiation at both ends ----------------------------------------------------------------------------
for b = [seq(1), seq(end)]
    if (b == 'G')||(b == 'C')
        h = h + N.WC.InitGC(1); s = s + 0.001*N.WC.InitGC(2);
    else
        h = h + N.WC.InitAT(1); s = s + 0.001*N.WC.InitAT(2);
    end
end
% The symmetry correction of Equation 3 is deliberately absent. See the header: this function scores
% a primer against a template, which is a heteroduplex, and the rotational degeneracy that the term
% accounts for does not exist for two distinct strands.

%% ---- salt -----------------------------------------------------------------------------------------------
% Every position of a fully paired duplex is a base pair, so the GC fraction over the sequence and
% over the paired positions are the same number here. dGMismatch has to count over paired positions
% only; this function does not, and the two agree on any input this function accepts.
gc = nnz((seq == 'G')|(seq == 'C'))/L;

if ctx.MonovalentOnly
    monEq = monv;
else
    monEq = N.Salt.MonovalentEquivalent(monv, biv, dNTP);
end
dsOwcz = N.Salt.Owczarzy(h, gc, L, monv, biv, dNTP);
if monEq > 0
    dsPer = N.Salt.PerPhosphate(L, monEq);
else
    dsPer = NaN;
end

switch lower(ctx.SaltModel)
    case 'perphosphate'
        if ~isfinite(dsPer)
            error('dGPrim_new:NonPositiveSalt', ...
                'The per-phosphate term needs a positive monovalent concentration; got %.3g M.', monEq);
        end
        dsUsed = dsPer; dsAlt = dsOwcz; route = 'perphosphate';
    case {'owczarzy', 'owczarzy_ref'}
        dsUsed = dsOwcz; dsAlt = dsPer; route = lower(ctx.SaltModel);
    otherwise
        error('dGPrim_new:BadSaltModel', ...
            ['ctx.SaltModel must be ''owczarzy_ref'' (default), ''owczarzy'' or ''perphosphate''; ' ...
            'got ''%s''. There is no ''auto'' here: a fully paired duplex has nothing to switch on.'], ...
            ctx.SaltModel);
end
s = s + 0.001*dsUsed;

%% ---- total ------------------------------------------------------------------------------------------------
TK = tann + 273.15;
dG = h - TK*s;

if nargout > 1
    info = struct('dH', h, 'dS', s, 'fGC', gc, 'Length', L, 'SaltRoute', route, ...
        'dS_salt_used', dsUsed, 'dS_salt_alternative', dsAlt, ...
        'dG_salt_model_gap', TK*abs(dsUsed - dsAlt)/1000, 'Temperature_C', tann);
end
end
