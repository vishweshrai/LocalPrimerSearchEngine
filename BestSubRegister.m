function s = BestSubRegister(basa, basb, ispair, minrun, monv, biv, dNTP, tann, applySym, gaterun)
%% BESTSUBREGISTER  Most stable duplex within one alignment register
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 07-11-2026 23:14 IST
% ____________________________________________________________________________________________________________
% WHAT THIS REPLACES, AND WHY
%
% DimerScan and MisprimerMismatch/SeedScan both scored a register by trimming once to its outermost
% Watson-Crick identity at each end:
%
%       first = find(ispair, 1, 'first');  last = find(ispair, 1, 'last');
%       basa = basa(first:last); ...
%
% That assumes the duplex runs between the outermost identities. It does not. A single accidental
% identity upstream of the real duplex drags the register out to meet it, and every mismatch in
% between then becomes INTERIOR and is charged as an internal loop rather than left as the frayed
% tail it is.
%
% MEASURED on the reverse primer GTGAGCTAACTCACATTAATTGCGTTGC against the forward primer of the
% same run, 17-10-2026, at 56.43 C:
%
%     what the trimming scores            what is physically there
%       TTGCGTTGC                           CGTTGC
%       |  |||| |                           |||| |
%       TACCGTTCC                           CGTTCC
%       +4.544 kcal/mol                     +0.245 kcal/mol
%
% One stray T.A identity three bases upstream makes the two intervening mismatches interior, which
% costs an internal loop of four, +3.84 at this temperature. The reported number is 4.3 kcal/mol
% away from the truth, which is larger than the whole dGEndMin threshold of -3.0. It did not flip a
% verdict on that pair because both values are positive. It would flip one on a pair whose real 3'
% duplex sits near -4, which would be reported as comfortably positive and would pass.
%
% In SeedScan it is worse, because the partner is a 7249 nt template rather than a 19 nt primer. The
% probability of at least one accidental identity in the first eleven primer positions at a random
% register is 1 - (3/4)^11 = 0.958, so the artifact fires on essentially every register, and an
% internal loop of 26 costs +6.79. It also corrupted the anchoring test, since that reads the same
% ballooned bounds.
%
% ____________________________________________________________________________________________________________
% WHAT THIS DOES INSTEAD
%
% Minimise over every admissible sub-register rather than taking one trim. A sub-register [i, j] is
% admissible when both its ends are base-paired, it is at least two positions long, and it carries a
% contiguous run of at least minrun pairs. Ends must be paired because a duplex terminates at a base
% pair; anything beyond it is a dangling end or a frayed tail and must not be charged as a loop.
%
% Two answers come out of the one sweep:
%   s.dGfree   the minimum over all admissible [i, j]
%   s.dGgate   the same minimum restricted to sub-registers whose longest run reaches gaterun,
%              which is the population a rejection rule may act on. Equal to s.dGright whenever
%              gaterun equals minrun. iGate, jGate and runGate describe that sub-register.
%   s.dGright  the minimum over admissible [i, jRight] with the right edge pinned at the OUTERMOST
%              3' pair. This is the quantity an extension test needs: the duplex has to reach the
%              3' terminus to be extendable, so the right edge is not free to move inward.
%
% COST. The sweep is O(p^2) calls to dGMismatch per register, p being the number of paired
% positions, against one call before. Sized against the 17-10-2026 run: 100 to 170 surviving
% registers per primer per strand, but only 23 primers ever reach the mispriming stage, so this is
% roughly a hundredfold on a stage that was taking seconds. Efficiency was explicitly not a
% constraint for this project; correctness of this number is the whole point of the stage.
%
% RELATION TO THE PRINCIPLED VERSION. The physical quantity is an ensemble, so the right answer is
% -RT*ln(SUM over sub-registers of exp(-dG/RT)) rather than the minimum. The minimum sits within
% about RT*ln(number of comparable modes) of that sum, a few tenths of a kcal/mol, and the sum is a
% one-line change to the accumulator below once anyone wants it. See Stage 2.3 of the revision plan.
%
% ____________________________________________________________________________________________________________
% INPUTS
%   basa     query bases over the UNTRIMMED register, char
%   basb     opposing bases over the same register, char, in the frame where IDENTITY means a
%            Watson-Crick pair (that is, basb is the partner read so that basa(k) == basb(k) when
%            paired). Both callers already build it that way.
%   ispair   logical, (basa == basb)
%   minrun   contiguous paired run a sub-register must carry
%   monv, biv, dNTP   concentrations, M
%   tann     temperature, C
%
% OUTPUT  s, struct
%   .dGfree  .iFree  .jFree  .runFree     best anywhere in the register; Inf and NaN when none
%   .dGright .iRight .jRight .runRight    best with the right edge at the outermost 3' pair
%   .trunRight .iTrunRight               length and start of the contiguous paired run TOUCHING the
%                                        right edge, which runRight does not give; see note below (G9)
%   .nTested number of sub-registers scored, for sizing
% ____________________________________________________________________________________________________________

s = struct('dGfree', Inf, 'iFree', NaN, 'jFree', NaN, 'runFree', 0, ...
           'dGright', Inf, 'iRight', NaN, 'jRight', NaN, 'runRight', 0, ...
           'dGgate', Inf, 'iGate', NaN, 'jGate', NaN, 'runGate', 0, ...
           'trunRight', 0, 'iTrunRight', NaN, 'nTested', 0);

basa = char(basa); basb = char(basb);
if (nargin < 9)||isempty(applySym), applySym = false; end   % heteroduplex unless told otherwise
% THE GATE FLOOR, WHICH IS NOT ALWAYS THE SURVEY FLOOR.
% minrun decides which sub-registers are SCORED. gaterun decides which of the scored ones are
% eligible to be the single worst site a rejection rule acts on. They are the same number unless
% the caller deliberately widens the survey below the gate, and when they differ the two winners
% can sit at the same locus and be different sub-registers: the survey admits a shorter run, that
% shorter run can be the more stable duplex, and it wins dGright while being ineligible to reject.
% Tracking the gate winner separately is what stops a widened survey from hiding the locus's
% gate-eligible register, which is the one the margin and the rejection both need.
if (nargin < 10)||isempty(gaterun), gaterun = minrun; end
ispair = logical(ispair(:)).';
p = find(ispair);
if numel(p) < 1, return; end
jRight = p(end);

% Contiguous run of pairs carrying the primer 3' terminus, which is not the same object as
% runRight. runRight is the longest run anywhere inside the winning anchored sub-register; this is
% the run that touches position n, and it is the quantity Rychlik (1995, Biotechniques 18:84-90)
% measured when he reported priming becoming detectable once the continuous 3'-terminal duplex
% beat -11 kcal/mol. His primers carried mismatched 5' termini while their 3' portion formed a
% continuous duplex with the target, so the duplex he measured reaches the 3' terminus and the
% count has to start there.
% Counting back from jRight instead returns an interior run whenever position n is mispaired, and
% an interior run scored against a terminal threshold rejects a primer that cannot prime from that
% site at all: a mispaired 3' terminus is not extendable, so no run behind it can make it prime.
% Anchored at n, a mispaired terminus gives trunRight = 0 and iTrunRight = n + 1, which
% TerminalRundG reads as fewer than two positions and scores as zero.
% Computed here because this function already owns ispair; it is scored by the caller, at the
% caller's choice of reference state.
n = numel(ispair);
q = n; s.trunRight = 0;
while (q >= 1) && ispair(q), s.trunRight = s.trunRight + 1; q = q - 1; end
s.iTrunRight = n - s.trunRight + 1;

for a = 1:numel(p)
    i = p(a);
    for b = numel(p):-1:a
        j = p(b);
        if (j - i + 1) < 2, continue; end          % dGMismatch needs a nearest-neighbour step
        if MaxRunLocal(ispair(i:j)) < minrun, continue; end
        % Argument 7 of dGMismatch gates ONLY the symmetry correction of Equation 3, and it has to
        % follow the caller rather than be fixed here, because this function serves two kinds of
        % duplex. MisprimerMismatch passes a primer against a TEMPLATE footprint, where the two
        % strands are distinct molecules, the twofold rotational degeneracy does not exist and the
        % term must be off. DimerScan passes a primer against ITSELF, where two identical strands
        % associate and a self-complementary register genuinely does carry that degeneracy, so the
        % term must be on. Hardcoding either value is wrong for the other caller by about
        % 0.46 kcal/mol at 58 C. Arguments 8 and 9 stay true: these are real duplex ends and do owe
        % helix initiation.
        g = dGMismatch(basa(i:j), char(seqcomplement(basb(i:j))), monv, biv, dNTP, tann, ...
            applySym, true, true);
        s.nTested = s.nTested + 1;
        if ~isfinite(g), continue; end
        if g < s.dGfree
            s.dGfree = g; s.iFree = i; s.jFree = j; s.runFree = MaxRunLocal(ispair(i:j));
        end
        runij = MaxRunLocal(ispair(i:j));
        if (j == jRight) && (g < s.dGright)
            s.dGright = g; s.iRight = i; s.jRight = j; s.runRight = runij;
        end
        if (j == jRight) && (runij >= gaterun) && (g < s.dGgate)
            s.dGgate = g; s.iGate = i; s.jGate = j; s.runGate = runij;
        end
    end
end
end

%% ____________________________________________________________________________________________________________
function r = MaxRunLocal(mask)
r = 0; c = 0;
for q = 1:numel(mask)
    if mask(q)
        c = c + 1; if c > r, r = c; end
    else
        c = 0;
    end
end
end
