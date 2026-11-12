function [tmA, info] = TmAmplicon(ampseq, monv, biv, dNTP, useEquiv)
%% TMAMPLICON  Melting temperature of a long double-stranded amplicon
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 05-11-2026 18:24 IST
% ____________________________________________________________________________________________________________
% PHYSICS
%
% A nearest-neighbour calculation is wrong for an amplicon of this length, and not because it is
% expensive. For a duplex of a few thousand base pairs the helix growth steps dominate the helix
% initiation step by three orders of magnitude, so the bimolecular initiation term that makes
% oligonucleotide Tm depend on strand concentration becomes negligible. The transition is
% pseudo-first-order and no concentration effect survives. Feeding such a duplex to a two-state
% oligonucleotide formula therefore adds a concentration dependence that is not physically there.
%
% The empirical %GC form is used instead, which is what Pierce et al. apply to the same quantity:
%
%       TmA = 81.5 + 16.6*log10( [M] / (1 + 0.7*[M]) ) + 0.41*(%G + %C) - 500/length
%
% The three sequence-dependent terms carry distinct physics. The 0.41*(%G+%C) term is the extra
% stacking and hydrogen bonding of G-C relative to A-T, averaged over the duplex. The -500/length
% term is the end effect: a finite duplex frays at both ends, and the fraying costs a fixed amount
% of stability that is diluted as the duplex lengthens. The logarithmic salt term is counterion
% condensation, written in the Wetmur form where the denominator (1 + 0.7*[M]) bends the response at
% high salt where the polyelectrolyte is already screened.
%
% WHAT THIS QUANTITY IS FOR
%
% In an asymmetric reaction the accumulating single strand is itself a competitor. Once the limiting
% primer is spent, every excess primer molecule has to find its site on a template strand that the
% previously made excess-primer strand is also trying to re-anneal to. The rate and extent of linear
% synthesis are therefore set by the relative affinities of the excess primer and the accumulated
% product for the same template, which is the difference TmA - TmX. Pierce et al. measured the
% consequence directly: the early slope of single-strand accumulation fell from 7.6 to 4.2 to 2.3 as
% TmA - TmX widened through 17, 20 and 23 degrees, with the best rates at a difference near 13.
%
% SALT CONVENTION
%
% The published formula takes [M] as the TOTAL MONOVALENT cation concentration, and that is the
% default here, so the number reproduces the figures in the source. Passing useEquiv = true
% substitutes the von Ahsen monovalent equivalent instead,
%
%       [Na+]eq = [Na+] + 120*sqrt([Mg2+]free)             concentrations in MILLIMOLAR
%
% which accounts for the magnesium the reaction actually contains. The two differ by several degrees
% at ordinary PCR magnesium, so the flag must be set deliberately and the same convention used on
% both sides of any TmA - TmX comparison. The coefficient 120 is calibrated for millimolar and the
% conversion is written out below rather than folded into a constant, because folding it in is how
% it goes missing.
%
% WHERE THE FREE MAGNESIUM COMES FROM, AND A DEPARTURE FROM THE SOURCE
%
% von Ahsen writes the bracket as [Mg2+] - [dNTP], a plain difference, and regressed the coefficient
% 120 against that definition. This file passes the speciated root from FreeMgSpeciated.m, which
% solves the 1:1 chelation equilibrium rather than taking its strong-binding limit. That is a
% departure from the published form and is recorded here as one.
%
% The reason is the other consumer. NNParams.m has used the speciated root throughout, so the suite
% held two different definitions of one quantity, and the amplicon melting temperature reported to
% the user was computed from the definition that clamps. MEASURED through the formula above on a
% 3940 nt amplicon at 49.85 percent GC in 40 mM monovalent, the two definitions differ by
%
%       2.5 mM Mg / 0.4 mM dNTP     0.008 degrees
%       2.5 mM Mg / 0.8 mM dNTP     0.023 degrees
%       3.0 mM Mg / 1.2 mM dNTP     0.030 degrees
%       2.5 mM Mg / 1.6 mM dNTP     0.141 degrees
%       1.5 mM Mg / 1.2 mM dNTP     0.580 degrees
%       1.0 mM Mg / 0.8 mM dNTP     0.703 degrees
%       2.0 mM Mg / 2.0 mM dNTP     6.252 degrees
%
% Everything above magnesium-in-excess sits inside the 1.76 degree standard error von Ahsen reports
% for his own fit on 217 matched duplexes, by one to two orders of magnitude at the compositions
% this suite runs at. The last row is the difference form clamping to zero, which drops the whole
% magnesium contribution in one step at a composition the chemistry passes through smoothly. The
% substitution is therefore invisible wherever his data constrain the answer, and it keeps the
% number physical where they do not.
%
% The price is that the speciated root carries the chelation constant into the far field. Over
% Ka = 1E4 to 1E5 /M the melting temperature moves 0.057 degrees at 2.5 mM Mg against 0.8 mM dNTP
% and 2.14 degrees at 2.0 mM against 2.0 mM. The body below therefore evaluates the root at both
% constants, forms the salt term at each, and warns when the spread exceeds von Ahsen's own standard
% error. info.KaSpread carries that number on every call. The buffer announces when it has entered
% the regime where an association constant is setting the answer, rather than a threshold being
% guessed in advance.
%
% ____________________________________________________________________________________________________________
% INPUTS
%   ampseq    amplicon sequence, char or string, or a scalar GC FRACTION in [0,1] when the length is
%             supplied through info.Length by the caller (see the two-argument form below)
%   monv      monovalent cation, M
%   biv       divalent cation, M.      Used only when useEquiv is true
%   dNTP      total dNTP, M.           Used only when useEquiv is true
%   useEquiv  false (default) for the published monovalent form, true for the von Ahsen equivalent
%
% OUTPUTS
%   tmA       amplicon melting temperature, degrees C
%   info      struct: .Length, .GC (fraction), .NaEff (M), .Formula,
%             .FreeMg (M, NaN when useEquiv is false), .KaSpread (degrees C, 0 when useEquiv
%             is false) the sensitivity of tmA to the chelation constant over Ka = 1E4 to 1E5
%
% ____________________________________________________________________________________________________________
% REFERENCES
%   Pierce KE, Sanchez JA, Rice JE, Wangh LJ. Linear-After-The-Exponential (LATE)-PCR: primer design
%       criteria for high yields of specific single-stranded DNA and improved real-time detection.
%       Proc Natl Acad Sci USA 2005;102:8609-14.  doi:10.1073/pnas.0501946102
%   Wetmur JG. DNA probes: applications of the principles of nucleic acid hybridization.
%       Crit Rev Biochem Mol Biol 1991;26:227-59.  doi:10.3109/10409239109114069
%   von Ahsen N, Wittwer CT, Schutz E. Oligonucleotide melting temperatures under PCR conditions:
%       nearest-neighbor corrections for Mg2+, deoxynucleotide triphosphate, and dimethyl sulfoxide
%       concentrations with comparison to alternative empirical formulas.
%       Clin Chem 2001;47:1956-61.  PMID 11673362
%   Owczarzy R, Moreira BG, You Y, Behlke MA, Walder JA. Predicting stability of DNA duplexes in
%       solutions containing magnesium and monovalent cations.
%       Biochemistry 2008;47:5336-53.  doi:10.1021/bi702363u
% ____________________________________________________________________________________________________________

% THE SALT ARGUMENTS ARE REQUIRED, BECAUSE THE WRONG ANSWER HERE IS A PLAUSIBLE ONE.
% These three used to default to false, 0 and 0, so TmAmplicon(ampSeq, monv) returned a
% monovalent-only melting temperature with nothing printed to say so. On a 40 mM monovalent buffer
% with 2.5 mM magnesium against 0.8 mM dNTP that reading is 10.12 degrees below the corrected one,
% and it is perfectly believable on its face, which is what makes a silent default the worst kind
% here. Both conventions are wanted by real callers: the engine takes the corrected value as the
% amplicon's real melting temperature and the monovalent one for Pierce's gap, so neither can be
% the default. Pass all five, with useEquiv saying which convention you are asking for; pass zeros
% explicitly for the monovalent reading.
if nargin < 5
    error('TmAmplicon:SaltNotStated', ...
        ['TmAmplicon needs all five arguments: (ampseq, monv, biv, dNTP, useEquiv). Pass the ' ...
        'reaction magnesium and dNTP with useEquiv true for the corrected melting temperature, ' ...
        'or (ampseq, monv, 0, 0, false) for the monovalent reading Pierce''s criteria are ' ...
        'defined in. The two differ by about 10 degrees in a typical PCR buffer, so there is no ' ...
        'safe default.']);
end
if isempty(useEquiv), useEquiv = false; end
if isempty(biv),      biv  = 0; end
if isempty(dNTP),     dNTP = 0; end

seq = upper(char(ampseq));
n   = numel(seq);
if n < 100
    error('TmAmplicon:TooShort', ...
        ['The amplicon is %d nt. The %%GC form assumes helix growth dominates initiation, which ' ...
        'fails below a few hundred base pairs. Use the nearest-neighbour route for short duplexes.'], n);
end
gc = nnz((seq == 'G')|(seq == 'C'))/n;

if ~isnumeric(monv)||~isscalar(monv)||~isfinite(monv)||(monv <= 0)
    error('TmAmplicon:NoSalt', ...
        ['Monovalent cation concentration must be a positive finite scalar, in MOLAR. The salt ' ...
        'term is a logarithm, so zero has no reading and a millimolar value passed here is off ' ...
        'by three orders of magnitude.']);
end

% von Ahsen's reported standard error on his own fit, 1.76 degrees over 217 matched duplexes. It is
% used as the threshold for the sensitivity warning below rather than a number chosen here, because
% a departure smaller than the scatter of the fit the correction comes from is not a departure that
% the fit can resolve.
KaSpreadLimit = 1.76;

kaSpread = 0;
freemg   = NaN;
if useEquiv
    if ~isscalar(biv)||~isscalar(dNTP)
        error('TmAmplicon:SaltNotScalar', ...
            'biv and dNTP must be scalar concentrations in molar when useEquiv is true.');
    end
    % Millimolar inside the bracket, as the coefficient requires, then back to molar. The free
    % magnesium is the speciated root, not the [Mg] - [dNTP] difference; see the header.
    [freemg, mgInfo] = FreeMgSpeciated(biv, dNTP);
    naEff = NaEquivalent(monv, freemg);
    % How much of this number the chelation constant is setting. The GC and length terms are common
    % to both evaluations and cancel exactly, so the spread is the salt term alone.
    kaSpread = 16.6*abs(SaltTerm(NaEquivalent(monv, mgInfo.FreeAtKaLo)) ...
                      - SaltTerm(NaEquivalent(monv, mgInfo.FreeAtKaHi)));
    if kaSpread > KaSpreadLimit
        warning('TmAmplicon:ChelationConstantSensitive', ...
            ['At %.2f mM magnesium against %.2f mM dNTP the free magnesium is %.4f mM, which is ' ...
            '%.0f percent of total, and the melting temperature moves %.2f degrees when the ' ...
            'chelation constant is varied over Ka = 1E4 to 1E5 /M. That exceeds the %.2f degree ' ...
            'standard error of the von Ahsen fit this conversion comes from, so Ka is setting ' ...
            'this answer rather than the buffer. Either move dNTP further below magnesium or ' ...
            'read this TmA as indicative.'], ...
            1E3*biv, 1E3*dNTP, 1E3*freemg, 100*mgInfo.FractionFree, kaSpread, KaSpreadLimit);
    end
else
    naEff = monv;
end

tmA = 81.5 + 16.6*SaltTerm(naEff) + 0.41*100*gc - 500/n;

info = struct('Length', n, 'GC', gc, 'NaEff', naEff, ...
    'Formula', "81.5 + 16.6 log10([M]/(1+0.7[M])) + 0.41(%GC) - 500/L", ...
    'FreeMg', freemg, 'KaSpread', kaSpread);
end

%% ____________________________________________________________________________________________________________
function e = NaEquivalent(monv, freemg)
% von Ahsen's monovalent equivalent. The coefficient 120 is calibrated for MILLIMOLAR, so the
% conversion in and back out is written rather than folded into a constant.
e = monv + 120*sqrt(max(freemg, 0)*1E3)*1E-3;
end

%% ____________________________________________________________________________________________________________
function t = SaltTerm(m)
% The counterion-condensation term of the Wetmur form, without its 16.6 coefficient. Factored out so
% that the sensitivity spread above is computed from the same expression the melting temperature is.
t = log10(m./(1 + 0.7*m));
end
