function [mfree, info] = FreeMgSpeciated(biv, dNTP, Ka)
%% FREEMGSPECIATED  Free magnesium in the presence of chelating dNTP
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 05-11-2026 18:24 IST
% ____________________________________________________________________________________________________________
% PHYSICS
%
% Nucleoside triphosphates carry three ionised phosphates and bind divalent magnesium tightly. In a
% PCR buffer the dNTPs are therefore a magnesium sink, and the quantity that screens the DNA
% backbone is what is left over. Writing the 1:1 complex as
%
%       Mg2+ + dNTP <-> Mg.dNTP        Ka = [Mg.dNTP] / ([Mg2+][dNTP])
%
% and imposing mass balance on both magnesium and dNTP gives, with m the free magnesium,
%
%       Ka*m^2 + m*(1 + Ka*[dNTP]tot - Ka*[Mg]tot) - [Mg]tot = 0
%
% whose positive root is returned. Ka = 3E4 /M is the default, the value used by Owczarzy et al.
%
% WHY NOT THE DIFFERENCE
%
% The stand-in max([Mg]tot - [dNTP]tot, 0) is the Ka -> infinity limit of the same equation, which
% is where every dNTP molecule carries a magnesium and none dissociates. At 3E4 /M and 0.8 mM dNTP
% the product Ka*[dNTP] is about 24, so binding is strong and finite, and the two agree to under one
% percent at ordinary PCR compositions. The difference form fails in one specific way: it clamps at
% exactly zero once dNTP reaches magnesium, and the clamp is a step rather than a taper. Two
% consequences follow. In the von Ahsen monovalent equivalent the sodium-equivalent collapses to the
% bare monovalent concentration, worth up to 4.6 degrees on a 21-mer at 40 mM monovalent. In the
% Owczarzy correction the ratio sqrt([Mg]free)/[Mon+] collapses to zero, which is below the 0.22
% competition threshold, so the function drops into its monovalent-only branch and the magnesium
% correction vanishes outright rather than shrinking. Owczarzy et al. place typical PCR buffers of
% 1.5 to 5 mM magnesium with 20 to 100 mM monovalent above that threshold.
%
% SENSITIVITY TO Ka
%
% The root is insensitive to Ka while magnesium is in comfortable excess and becomes governed by it
% near stoichiometry. MEASURED through the nearest-neighbour route on the 21-mer
% GCCAACCAAGGCGGCGCTAAT at 40 mM monovalent and 1 uM primer, varying Ka over 5E3 to 1E5 /M moves the
% melting temperature by
%
%       2.5 mM Mg / 0.8 mM dNTP     0.09 degrees
%       2.0 mM Mg / 1.2 mM dNTP     0.39 degrees
%       1.5 mM Mg / 1.6 mM dNTP     2.37 degrees
%       2.0 mM Mg / 2.0 mM dNTP     2.14 degrees
%
% Through the Wetmur %GC route in TmAmplicon.m the same Ka range over 1E4 to 1E5 /M gives 0.057
% degrees at 2.5 mM against 0.8 mM and 2.14 degrees at 2.0 mM against 2.0 mM.
%
% The second output therefore carries the root evaluated at Ka = 1E4 and 1E5 /M, so a caller that
% forms a temperature can report its own Ka sensitivity and warn when the buffer has entered the
% regime where an uncited constant is setting the answer. TmAmplicon does exactly this.
%
% ____________________________________________________________________________________________________________
% INPUTS
%   biv       total divalent cation (magnesium), M. Scalar or array.
%   dNTP      total dNTP across all four species, M. Scalar or array, default 0.
%   Ka        association constant for the 1:1 complex, 1/M. Default 3E4.
%
% OUTPUTS
%   mfree     free magnesium, M, same size as the expanded inputs
%   info      struct, computed only when a second output is requested:
%               .Ka            the constant used
%               .Crude         max(biv - dNTP, 0), the Ka -> infinity limit, for comparison
%               .FreeAtKaLo    the root at Ka = 1E4 /M
%               .FreeAtKaHi    the root at Ka = 1E5 /M
%               .FractionFree  mfree ./ biv, set to 1 where biv is zero
%
% ____________________________________________________________________________________________________________
% REFERENCES
%   Owczarzy R, Moreira BG, You Y, Behlke MA, Walder JA. Predicting stability of DNA duplexes in
%       solutions containing magnesium and monovalent cations.
%       Biochemistry 2008;47:5336-53.  doi:10.1021/bi702363u
%   von Ahsen N, Wittwer CT, Schutz E. Oligonucleotide melting temperatures under PCR conditions:
%       nearest-neighbor corrections for Mg2+, deoxynucleotide triphosphate, and dimethyl sulfoxide
%       concentrations with comparison to alternative empirical formulas.
%       Clin Chem 2001;47:1956-61.  PMID 11673362
% ____________________________________________________________________________________________________________

if (nargin < 2)||isempty(dNTP), dNTP = 0;   end
if (nargin < 3)||isempty(Ka),   Ka   = 3E4; end

if ~isnumeric(biv)||~isnumeric(dNTP)||~isnumeric(Ka)
    error('FreeMgSpeciated:NonNumeric', ...
        ['biv, dNTP and Ka must be numeric. A char or string concentration reaches here when a ' ...
        'caller has passed a table cell through without converting it.']);
end
if any(~isfinite(biv(:)))||any(~isfinite(dNTP(:)))
    error('FreeMgSpeciated:NonFiniteConcentration', ...
        'biv and dNTP must be finite. A NaN here propagates silently into every free energy.');
end
if any(biv(:) < 0)||any(dNTP(:) < 0)
    error('FreeMgSpeciated:NegativeConcentration', ...
        ['Concentrations must be non-negative, and they are in MOLAR. A negative value almost ' ...
        'always means a millimolar quantity has been subtracted from a molar one somewhere ' ...
        'upstream.']);
end
if ~isscalar(Ka)||(Ka <= 0)
    error('FreeMgSpeciated:BadKa', ...
        'Ka must be a positive scalar in 1/M. Received a value of %g.', Ka);
end

mfree = PositiveRoot(biv, dNTP, Ka);

if nargout > 1
    info.Ka         = Ka;
    info.Crude      = max(biv - dNTP, 0);
    info.FreeAtKaLo = PositiveRoot(biv, dNTP, 1E4);
    info.FreeAtKaHi = PositiveRoot(biv, dNTP, 1E5);
    ff              = mfree./max(biv, realmin);
    ff(biv == 0)    = 1;
    info.FractionFree = ff;
end
end

%% ____________________________________________________________________________________________________________
function m = PositiveRoot(biv, dNTP, Ka)
% The positive root of Ka*m^2 + b*m - biv = 0, written in whichever of the two algebraically
% identical forms is numerically stable. When b is large and positive, which is the dNTP-rich
% corner, (-b + sqrt(b^2 + 4*Ka*biv)) subtracts two nearly equal numbers and sheds significant
% figures. The conjugate form 2*biv/(b + sqrt(...)) adds instead of subtracting and keeps them.
b    = 1 + Ka.*dNTP - Ka.*biv;
disc = sqrt(b.^2 + 4.*Ka.*biv);
m    = (-b + disc)./(2.*Ka);
pos  = (b > 0) & ((b + disc) > 0);
if any(pos(:))
    alt      = (2.*biv)./(b + disc);
    bigm     = m + zeros(size(pos));      % force full size before indexed assignment
    bigalt   = alt + zeros(size(pos));
    bigm(pos)= bigalt(pos);
    m        = bigm;
end
m = max(m, 0);
end
