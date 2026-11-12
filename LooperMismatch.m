function [stemalign, check, dGhp, hpTm, lockTm, foldFrac, best] = LooperMismatch(primer, monv, biv, dNTP, tann, dGhpmin, tAnneal, TmMargin, BulgeMax, N)
%% LOOPERMISMATCH  Hairpin screen for a single primer, judged on melting temperature
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 04-11-2026 23:06 IST
% ____________________________________________________________________________________________________________
% PHYSICS
%
% A hairpin costs the reaction in two unrelated ways, and only one of them is a specificity problem.
%
% The first is availability. A folded primer cannot find its template, so the reaction sees a lower
% effective primer concentration than the one pipetted. Folding is an equilibrium, so what matters is
% the FOLDED FRACTION at the annealing temperature, f = K/(1+K) with K = exp(-dG/RT) for a
% unimolecular transition. That fraction is one half when the annealing temperature equals the
% hairpin melting temperature, which is why the melting temperature and not the free energy is the
% quantity to compare against the reaction.
%
% The second is blockage. If the stem pairs the primer's own 3'-terminal base, the polymerase has no
% free 3'-OH to extend even when the primer is correctly bound to template. That is a hard stop
% rather than a dilution, and it is the case this function gates on.
%
% WHY FREE ENERGY AT THE ANNEALING TEMPERATURE IS THE WRONG GATE
%
% Free energy passes through zero at a structure's own melting point by definition. Evaluating dG
% within a few degrees of the hairpin melting temperature therefore returns a number near zero
% whatever the hairpin is doing, and a dG cutoff applied there cannot fire for exactly the hairpins
% that melt near the annealing step. Measured on this project's own primers at the reaction
% conditions: hairpin dG came to -0.104 and -0.088 kcal/mol against a -3.0 cutoff, while the two
% hairpins melt at 56.3 and 60.9 degrees, so roughly half of each primer is folded while annealing.
% The dG rule passed both. The melting-temperature rule rejects the one that also blocks the 3' end.
%
% HOW THE MELTING TEMPERATURE IS OBTAINED
%
% dG(T) = dH - T*dS. Both contributions here are linear in temperature by construction: the stem is
% a sum of nearest-neighbour stacks with temperature-independent dH and dS, and the loop increment
% is purely entropic, since every loop dH is taken as zero and dGLoop scales its 37-degree value by
% T/310.15. Scoring the SAME structure at two temperatures therefore recovers both coefficients
% exactly,
%
%       dS = (dG(T1) - dG(T2)) / (T2 - T1),      dH = dG(T1) + T1*dS
%
% and the unimolecular melting temperature follows from dG(Tm) = 0 as Tm = dH/dS. No concentration
% term enters, because one molecule folds on itself.
%
% STRUCTURE ENUMERATION
%
% For a start m, arm length n and loop length o the register is fully determined: stem53(i) pairs
% with stem35(n+1-i). No alignment search is run, because a shifted register implies a different
% loop length than the o being scored, which would make the stem and loop terms describe different
% structures. The innermost and outermost pairs must both be Watson-Crick: an open innermost pair
% means the loop is really longer than o, and an open outermost pair means those bases are dangling
% ends, and both cases are reached by other (m, n, o) triples in the same enumeration.
%
% Stem thermodynamics use stacks only, salt corrected, with no bimolecular initiation, because a
% hairpin is unimolecular and Equation 10 of the review asks for the base-pair nearest-neighbour
% contributions rather than the initiation rows. Loops of three follow Equation 8, taking the
% hairpin-of-three increment plus a closing A-T penalty of 0.5 kcal/mol and no terminal mismatch.
% Longer loops follow Equations 9 and 10 and carry the MEASURED terminal mismatch at the closing
% pair, read from terminal_mismatch_nn_params.csv through dGTerminalMismatch. The single constant
% of -1.23 kcal/mol that stood here until H8 was the most stabilising entry in the whole table,
% handed to every hairpin regardless of sequence; the measured value is weaker for almost every
% sequence and the hairpins found here are correspondingly less stable.
%
% WHAT THE ENUMERATION DOES NOT REACH
%
% BULGES ARE IMPLEMENTED. H8 added them: one arm may carry b unpaired bases between an outer helix
% of n1 and an inner helix of n2, with b running to BulgeMax and either arm eligible. Setting
% b = 0 recovers the symmetric register exactly, so the earlier bulge-free search is a strict
% subset of this one. The under-detection figure that stood here, +2.30 kcal/mol and a melting
% temperature of 27.2 degrees on GACGCTTACTGGAACGGTAAA against Primer3's 56.3, was measured BEFORE
% that change and before the triloop bonus, and must be re-measured before it is quoted anywhere.
%
% TRILOOP AND TETRALOOP BONUSES ARE IMPLEMENTED, through SpecialLoop at loop lengths 3 and 4,
% backed by special_loop_params.csv. The review places the values in supplementary material and
% they come from the Primer3 distribution.
%
% STILL ABSENT: internal loops with unpaired bases on BOTH arms at once, multibranch structures,
% and pseudoknots. A stem interrupted on one side only is reachable; one interrupted on both sides
% at the same position is not.
%
% ____________________________________________________________________________________________________________
% INPUTS
%   primer    primer sequence, 5' -> 3', char
%   monv      monovalent cation, M
%   biv       divalent cation, M
%   dNTP      total dNTP, M
%   tann      temperature at which dG is reported, degrees C
%   dGhpmin   free-energy fallback cutoff, used only when tAnneal is not supplied
%   tAnneal   reaction annealing temperature, degrees C. Supplying it selects the melting-temperature
%             gate, which is the intended route
%   TmMargin  degrees of headroom. A 3'-blocking hairpin rejects when its melting temperature
%             reaches tAnneal - TmMargin. Default 0
%   BulgeMax  largest bulge admitted on either arm, nt. Default 3
%   N         NNParams struct. Default a fresh NNParams() call, which is cached internally
%
% OUTPUTS
%   stemalign 3-row alignment block for the most stable hairpin; empty when none qualifies
%   check     false when the candidate is rejected
%   dGhp      dG of the most stable enumerated hairpin at tann, kcal/mol. A POSITIVE value means no
%             hairpin forms in the bulge-free family; it is reported rather than zeroed, because the
%             margin to zero is itself informative
%   hpTm      melting temperature of that hairpin, degrees C; NaN when none was found
%   lockTm    melting temperature of the most stable hairpin ENDING at the primer's 3' terminus,
%             so that its outermost pair buries the 3'-OH; NaN when no such hairpin exists
%   foldFrac  folded fraction of the most stable hairpin at tAnneal, 0 to 1
%   best      struct describing the winning structure: S, dG, Bulge, BulgeSide, Span and LockS
%
% ____________________________________________________________________________________________________________
% REFERENCES
%   SantaLucia J, Hicks D. The thermodynamics of DNA structural motifs. Annu Rev Biophys Biomol
%       Struct 2004;33:415-40.  doi:10.1146/annurev.biophys.32.110601.141800
% ____________________________________________________________________________________________________________

if (nargin < 6)||isempty(dGhpmin),  dGhpmin  = -3.0; end
if (nargin < 9)||isempty(BulgeMax), BulgeMax = 3;    end
if (nargin < 10)||isempty(N),       N = NNParams();  end
if (nargin < 7),                    tAnneal  = [];   end
if (nargin < 8)||isempty(TmMargin), TmMargin = 0;    end

RGAS   = 1.9872E-3;   % kcal/mol/K
DT     = 10;          % degrees between the two scoring temperatures used to recover dH and dS

check     = true;
stemalign = repmat(' ', 3, 0);
best      = struct('S', [], 'dG', Inf, 'Bulge', 0, 'BulgeSide', 0, 'Span', 0, 'LockS', []);
dGhp      = 0;
hpTm      = NaN;
lockTm    = NaN;
foldFrac  = 0;

primer  = char(primer);
plen    = strlength(primer);
loopmin = 3;          % DNA hairpin loops shorter than three nucleotides are sterically unfavourable
stemmin = 3;          % shortest stem arm scored

bestdG = Inf;  bestS = [];       % most stable hairpin anywhere in the primer
lockTmBest = -Inf;  lockS = [];  % HOTTEST hairpin that pairs the 3'-terminal base
% THE BLOCKING HAIRPIN IS CHOSEN BY MELTING TEMPERATURE, BECAUSE THAT IS WHAT THE GATE ASKS.
% The gate below fires when a structure pairing the 3' terminus is still half folded or more at
% the annealing temperature, which is a statement about Tm. Choosing the blocking structure by
% minimum free energy answers a different question. Writing dG(T) = -|dS|*(Tm - T) shows why:
% the most negative dG at one temperature maximises |dS|*(Tm - T), so a long stem with a large
% entropy and a modest Tm beats a short stem with a small entropy and a higher Tm. The short one
% is the one that is still folded while annealing, and under a minimum-dG rule it is never
% tested. Over primers carrying two or more blocking structures the gap between the two rules
% reaches 20 degrees.
% THE LOCK TEST IS ON WHERE THE STRUCTURE ENDS, NOT ON HOW LONG IT IS.
% A structure blocks extension when its last paired base IS the primer's 3' terminus, which is
% m + SpanOf(S) - 1 == plen. Testing SpanOf(S) == plen instead is not equivalent: since
% ScoreHairpinGeneral already rejects m + SpanOf(S) - 1 > plen, that form can only be satisfied at
% m == 1, which silently drops every blocking hairpin that leaves a 5' overhang from the lock set
% and leaves lockTm at NaN.

% A STRUCTURE IS SIX NUMBERS: [m n1 b side n2 o].
%   m     start of the structure in the primer
%   n1    length of the OUTER helix, the one carrying the free end of the stem
%   b     bulge size, 0 for a continuous stem
%   side  5 or 3, which arm carries the bulge; meaningless when b is 0
%   n2    length of the INNER helix, the one closing on the loop; 0 when b is 0
%   o     loop length
% With b = 0 this reduces exactly to the symmetric register the earlier version enumerated, so the
% old search is a strict subset of this one and nothing it used to find is lost.
for m = 1:plen
    nmax = floor(0.5*(plen - m + 1 - loopmin));
    for n1 = stemmin:nmax
        % ---- continuous stem, b = 0 -----------------------------------------------------------
        for o = loopmin:(plen - m - 2*n1 + 1)
            S = [m n1 0 0 0 o];
            [ok, dG, align] = ScoreHairpinGeneral(primer, S, monv, biv, dNTP, tann, stemmin, N);
            if ~ok, continue; end
            if dG < bestdG, bestdG = dG; bestS = S; stemalign = align; end
            if m + SpanOf(S) - 1 == plen
                tmL = TmFromdG(primer, S, dG, monv, biv, dNTP, tann, DT, stemmin, N);
                if isfinite(tmL)&&(tmL > lockTmBest), lockTmBest = tmL; lockS = S; end
            end
        end
        % ---- bulged stem ----------------------------------------------------------------------
        % One arm carries b unpaired bases between an outer helix of n1 and an inner helix of n2.
        % Both helices must reach stemmin, otherwise the structure is a short stem with a long
        % overhang rather than a bulged hairpin, and the loop database does not describe it.
        for b = 1:BulgeMax
            for side = [5 3]
                for n2 = stemmin:nmax
                    rem = plen - m + 1 - 2*n1 - 2*n2 - b;
                    for o = loopmin:rem
                        S = [m n1 b side n2 o];
                        [ok, dG, align] = ScoreHairpinGeneral(primer, S, monv, biv, dNTP, tann, stemmin, N);
                        if ~ok, continue; end
                        if dG < bestdG, bestdG = dG; bestS = S; stemalign = align; end
                        if m + SpanOf(S) - 1 == plen
                            tmL = TmFromdG(primer, S, dG, monv, biv, dNTP, tann, DT, stemmin, N);
                            if isfinite(tmL)&&(tmL > lockTmBest), lockTmBest = tmL; lockS = S; end
                        end
                    end
                end
            end
        end
    end
end

if isempty(bestS), return; end

dGhp = bestdG;
hpTm = MeltFromTwoPoints(primer, bestS, monv, biv, dNTP, tann, DT, stemmin, N);
if ~isempty(lockS)
    % Already measured during the enumeration, on the same two-point construction
    % MeltFromTwoPoints uses, so there is nothing to recompute here.
    lockTm = lockTmBest;
end
best = struct('S', bestS, 'dG', bestdG, 'Bulge', bestS(3), 'BulgeSide', bestS(4), ...
    'Span', SpanOf(bestS), 'LockS', lockS);

if ~isempty(tAnneal)&&isfinite(tAnneal)
    % Availability, reported for ranking: folded fraction of the best hairpin while annealing.
    [okA, dGatAnneal] = ScoreHairpinGeneral(primer, bestS, monv, biv, dNTP, tAnneal, stemmin, N);
    if okA
        K = exp(-dGatAnneal/(RGAS*(tAnneal + 273.15)));
        foldFrac = K/(1 + K);
    end
    % Blockage, the gate: a hairpin that pairs the 3' terminus and is half folded or more while
    % annealing leaves the polymerase no free 3'-OH on that fraction of the primer.
    if isfinite(lockTm)&&(lockTm >= tAnneal - TmMargin)
        check = false;
    end
else
    % No annealing temperature supplied, so fall back to the free-energy rule with its known blind
    % spot near the hairpin melting temperature.
    if dGhp < dGhpmin
        check = false;
    end
end
end

%% ____________________________________________________________________________________________________________
%% SPANOF  How many bases of the primer a structure occupies
function L = SpanOf(S)
L = 2*S(2) + 2*S(5) + S(3) + S(6);
end

%% ____________________________________________________________________________________________________________
%% SCOREHAIRPINGENERAL  One hairpin, continuous or bulged, at one temperature
%
% GEOMETRY. Reading the primer 5' to 3' from position m, the structure lays out as
%
%   b = 0            [ A5 n1 ][ loop o ][ A3 n1 ]
%   b > 0, side 5    [ A5 n1 ][ bulge b ][ B5 n2 ][ loop o ][ B3 n2 ][ A3 n1 ]
%   b > 0, side 3    [ A5 n1 ][ B5 n2 ][ loop o ][ B3 n2 ][ bulge b ][ A3 n1 ]
%
% A5 pairs with A3 reversed and B5 with B3 reversed. The bulged bases pair with nothing.
%
% WHY NO HELIX INITIATION IS CHARGED. The unified nearest-neighbour model charges an initiation at
% each end of a BIMOLECULAR duplex, because two separate strands have to find one another and
% nucleate a helix before any stack can form. A hairpin is unimolecular. The strand never parts, so
% there is no nucleation of that kind to pay for, and Equations 8 to 10 of SantaLucia and Hicks
% carry the whole cost of closing the structure in the loop initiation term instead.
%
% Both arms are therefore scored with applyTerm5 and applyTerm3 set false. Inside dGMismatch those
% two flags, and not applyInit, are what gate AddInit; applyInit gates the symmetry correction
% alone. With both false no initiation is charged at any of the four arm ends, and the loop block
% below supplies the only initiation the structure owes. Nothing is subtracted afterwards.
function [ok, dG, align] = ScoreHairpinGeneral(primer, S, monv, biv, dNTP, T, stemmin, N)
ok = false; dG = Inf; align = repmat(' ', 3, 0);
m = S(1); n1 = S(2); b = S(3); side = S(4); n2 = S(5); o = S(6);
plen = strlength(primer);
if (m < 1)||(m + SpanOf(S) - 1 > plen), return; end
if (o < 3)||(n1 < stemmin), return; end
if (b > 0)&&(n2 < stemmin), return; end

% ---- cut the primer into the pieces the layout defines ----
k = m;
A5 = primer(k:k+n1-1); k = k + n1;
% The bulged bases are kept as they are stepped over, so the alignment can show them on the arm
% that actually carries them instead of marking both arms with a gap.
bulge5 = ''; bulge3 = '';
if (b > 0)&&(side == 5), bulge5 = primer(k:k+b-1); k = k + b; end
if b > 0
    B5 = primer(k:k+n2-1); k = k + n2;
else
    B5 = '';
end
loop = primer(k:k+o-1); k = k + o;
if b > 0
    B3 = primer(k:k+n2-1); k = k + n2;
else
    B3 = '';
end
if (b > 0)&&(side == 3), bulge3 = primer(k:k+b-1); k = k + b; end
A3 = primer(k:k+n1-1);

% ---- the outer helix must be Watson-Crick at both of its own ends ----
pairA = (char(seqcomplement(A5)) == fliplr(A3));
if ~pairA(1), return; end
if (b == 0)&&~pairA(end), return; end
if MaxRun(pairA) < stemmin, return; end

if b > 0
    pairB = (char(seqcomplement(B5)) == fliplr(B3));
    if ~pairB(end), return; end                 % the pair that closes the loop
    if ~pairA(end)||~pairB(1), return; end      % the two pairs that close the bulge
    if MaxRun(pairB) < stemmin, return; end
else
    pairB = [];
end

% ---- stem free energy ----
dG = dGMismatch(A5, fliplr(A3), monv, biv, dNTP, T, false, false, false);
if b > 0
    dG = dG + dGMismatch(B5, fliplr(B3), monv, biv, dNTP, T, false, false, false);
    % The bulge itself. Its flanking pairs are the last pair of A and the first pair of B, each
    % written top base then bottom base.
    %
    % ARGUMENT ORDER ACROSS A SINGLE-BASE BULGE. dGBulge builds the retained nearest-neighbour step
    % from the first character of each pair argument, so the order in which the two pairs arrive
    % decides which step is looked up. The retained stack is the one read 5' to 3' along whichever
    % strand runs through the bulge unbroken. For a bulge on the 5' arm that strand is the 3' arm,
    % B3(end) to A3(1), and [A5(end) B5(1)] is its reverse complement. For a bulge on the 3' arm the
    % unbroken strand is the 5' arm itself, A5(end) to B5(1). Both cases therefore want pr5 first.
    % A nearest-neighbour table is symmetric under reverse complementation and NOT under plain
    % reversal, so passing pr3 first would look up the reverse of the correct step: AT would be read
    % as TA, which differs by 0.30 kcal/mol at 37 C, the widest such gap in Table 1.
    pr5 = [A5(end) A3(1)];
    pr3 = [B5(1) B3(end)];
    [dgb, ~] = dGBulge(b, pr5, pr3, T, N);
    dG = dG + dgb;
    closeTop = B5(end); closeBot = B3(1);
else
    closeTop = A5(end); closeBot = A3(1);
end

% ---- loop ----
% Equation 8 for a triloop, Equation 10 for longer loops. The terminal mismatch is the measured
% stack of the first and last loop bases on the closing pair, in place of the single constant the
% earlier version applied to every hairpin. That constant was -1.23, the MOST stabilising entry in
% the whole table, so it handed every hairpin the best closure that exists; the measured value is
% weaker for almost every sequence and the hairpins found here are correspondingly less stable.
if o == 3
    % Equation 8. A triloop takes no terminal mismatch; it takes a sequence bonus and, when the
    % closing pair is A-T, a penalty.
    dGloop = dGLoop('hairpin', 3, T, N) + SpecialLoop(closeTop, loop, closeBot, T, N);
    if any(closeTop == 'AT')
        dGloop = dGloop + N.Loop.ATClosureG;
    end
else
    % Equations 9 and 10. Length term plus the measured terminal mismatch, plus, at length 4, the
    % tetraloop bonus where the sequence has one.
    dGloop = dGLoop('hairpin', o, T, N) ...
        + dGTerminalMismatch(closeTop, loop(1), closeBot, loop(end), T, N);
    if o == 4
        dGloop = dGloop + SpecialLoop(closeTop, loop, closeBot, T, N);
    end
end
dG = dG + dGloop;
ok = isfinite(dG);
if ~ok, return; end

% ---- alignment block for reporting ----
if b == 0
    % THE BOTTOM ROW CARRIES THE 3' ARM'S OWN BASES, REVERSED.
    % A5(i) pairs with A3(end-i+1), so the partner of the top row read left to right is fliplr(A3),
    % which is the same frame the pairing test and dGMismatch use a few lines above. Taking the
    % reverse COMPLEMENT instead prints the complement of the real bases, and on a fully paired
    % stem that reproduces the top row, which is how the error hid.
    align = [A5; PairSymbols(pairA); fliplr(A3)];
else
    % THE BULGE IS DRAWN ON THE ARM THAT CARRIES IT.
    % b unpaired bases sit on ONE arm; the other arm runs straight past them. So the arm carrying
    % the bulge shows its bases and the opposite arm shows b gap characters. Marking both arms
    % would depict a symmetric internal loop, which is a different motif with different energetics.
    % The 3'-arm case is drawn reversed, like the rest of the bottom row, because the bulged bases
    % lie between B3 and A3 in the primer and the bottom row reads the arm from its 3' end.
    gap = repmat('-', 1, b);
    pad = repmat(' ', 1, b);
    if side == 5
        top = [A5 bulge5 B5];
        mid = [PairSymbols(pairA) pad PairSymbols(pairB)];
        bot = [fliplr(A3) gap fliplr(B3)];
    else
        top = [A5 gap B5];
        mid = [PairSymbols(pairA) pad PairSymbols(pairB)];
        bot = [fliplr(A3) fliplr(bulge3) fliplr(B3)];
    end
    align = [top; mid; bot];
end
end

%% ____________________________________________________________________________________________________________
%% SPECIALLOOP  Sequence bonus for an unusually stable triloop or tetraloop
% Keyed by the closing pair's top base, the loop 5' to 3', then the closing pair's bottom base.
% Returns zero when the loop is not one of the measured sequences, which is the common case.
function g = SpecialLoop(closeTop, loop, closeBot, T, N)
g = 0;
if ~isfield(N, 'Loop')||~isfield(N.Loop, 'SpecialKey')||isempty(N.Loop.SpecialKey)
    return
end
key = upper([closeTop loop closeBot]);
idx = find(strcmp(N.Loop.SpecialKey, key), 1);
if isempty(idx), return; end
g = N.Loop.SpecialdH(idx) - (T + 273.15)*N.Loop.SpecialdS(idx)/1000;
end

%% ____________________________________________________________________________________________________________
function tm = TmFromdG(primer, S, g1, monv, biv, dNTP, T1, dT, stemmin, N)
% Melting temperature of one structure when its free energy at T1 is already in hand.
% dG(T) is linear in T for a hairpin, so a single further evaluation recovers dH and dS exactly
% and Tm follows from dG(Tm) = 0. Used inside the enumeration, where every candidate has just
% been scored at T1 and a second full call would repeat work already done. Identical in value to
% MeltFromTwoPoints on the same structure.
tm = NaN;
[ok2, g2] = ScoreHairpinGeneral(primer, S, monv, biv, dNTP, T1 + dT, stemmin, N);
if ~ok2, return; end
K1 = T1 + 273.15;  K2 = K1 + dT;
dS = (g1 - g2)/(K2 - K1);          % kcal/mol/K
if ~isfinite(dS)||(abs(dS) < 1E-9), return; end
tm = (g1 + K1*dS)/dS - 273.15;
end

%% ____________________________________________________________________________________________________________
function tm = MeltFromTwoPoints(primer, S, monv, biv, dNTP, T1, dT, stemmin, N)
% dG(T) is linear in T for this structure, so two evaluations recover dH and dS exactly and the
% unimolecular melting temperature follows from dG(Tm) = 0.
tm = NaN;
T2 = T1 + dT;
[ok1, g1] = ScoreHairpinGeneral(primer, S, monv, biv, dNTP, T1, stemmin, N);
[ok2, g2] = ScoreHairpinGeneral(primer, S, monv, biv, dNTP, T2, stemmin, N);
if ~ok1||~ok2, return; end
K1 = T1 + 273.15;  K2 = T2 + 273.15;
dS = (g1 - g2)/(K2 - K1);          % kcal/mol/K
if ~isfinite(dS)||(abs(dS) < 1E-9), return; end
dH = g1 + K1*dS;                   % kcal/mol
tm = dH/dS - 273.15;
end

%% ____________________________________________________________________________________________________________
function runmax = MaxRun(mask)
runmax = 0; runlen = 0;
for q = 1:numel(mask)
    if mask(q)
        runlen = runlen + 1;
        if runlen > runmax, runmax = runlen; end
    else
        runlen = 0;
    end
end
end

%% ____________________________________________________________________________________________________________
function sym = PairSymbols(mask)
sym = repmat(' ', 1, numel(mask));
sym(mask) = '|';
end
