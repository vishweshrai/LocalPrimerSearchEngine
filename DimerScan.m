function [dGdimer, dGend, dimeralign, endalign] = DimerScan(seqa, seqb, monv, biv, dNTP, tann, minrun, tmplmin)
%% DIMERSCAN  Most stable antiparallel duplex between two sequences
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 12-11-2026 01:53 IST
% Most stable antiparallel duplex between seqa and seqb, located by enumerating every ungapped register
% rather than by a single Smith-Waterman call.

% Rationale: swalign is invoked throughout this suite with 'GapOpen', 1000, which forbids gaps, so the
% space it searches is exactly the set of ungapped registers. For two primers of 18 to 26 nt that set
% holds fewer than 50 members and can be enumerated in full. Doing so removes two defects of the swalign
% route:
%   (1) swalign returns the register that maximises the NUC44 identity score, not the register that
%       minimises DeltaG. Four contiguous G-C pairs score 20 and sit near -6 kcal/mol; six pairs across
%       an AT-rich stretch score 30 and may be less stable. The two optima are frequently different
%       registers, and it is the thermodynamic one that governs whether the dimer forms.
%   (2) swalign returns one alignment only, so a short 3'-anchored dimer, which is the only kind the
%       polymerase can extend, is discarded whenever a longer internal complementarity outscores it.

% Geometry: seqa is aligned against seqbrc = seqrcomplement(seqb). Identity in that frame is Watson-Crick
% pairing between seqa and seqb, because seqbrc(q) = comp(seqb(nb+1-q)), so seqa(p) == seqbrc(q) implies
% seqa(p) pairs with seqb(nb+1-q), and the index on seqb runs backwards as p advances. The base physically
% opposing seqa(p) is therefore comp(seqbrc(q)), which is what dGMismatch expects as its seq2 argument.

% CHANGE of 17-10-2026 (G3): each register is no longer trimmed once to its outermost Watson-Crick pair.
% That assumed the duplex runs between the outermost identities, and a single accidental identity
% upstream of the real duplex dragged the register out to meet it, turning the intervening mismatches
% into an interior loop. BestSubRegister.m now minimises over every admissible sub-register and
% returns two answers: the best anywhere, and the best whose right edge sits at the outermost 3' pair,
% which is the one an extension test needs. The old behaviour over-penalised by up to 4.3 kcal/mol on
% a register from this project's own primers. MaxRun below is retained for callers and is no longer
% used by the register loop.

% Template beyond the paired 3' end: a polymerase extends the 3' terminus of seqa by copying the partner,
% and the partner strand supplies a template only past the base that seqa's terminus is paired to. In the
% register frame that base is seqb(m) with m = nb + 1 - na + off, and the template available to copy is
% seqb(1:m-1), read 5' to 3' away from the pair. If the partner ends at the paired base the duplex has a
% blunt or recessed 3' end of seqa and nothing is left to copy, so the dimer is not extendable. tmplmin is
% the least number of partner nucleotides required 5' of the paired base. It defaults to 0, which keeps the
% earlier rule that any register pairing the 3' terminal base counts, and 1 is the physical minimum.
% The requirement acts on dGend and endalign only. dGdimer is an occupancy quantity and ignores it.

% dGdimer: lowest DeltaG over all registers, kcal/mol. Returns 0 when no register carries a run of minrun
%          consecutive Watson-Crick pairs, matching the "no dimer detected" sentinel used elsewhere.
% dGend: lowest DeltaG over those registers whose duplex includes and pairs the 3'-terminal base of seqa
%        and leaves at least tmplmin partner nucleotides 5' of the paired base. This is the extendable
%        subset, and it is the quantity to compare against a 3'-dimer cutoff.
% dimeralign, endalign: 3-row char arrays in swalign layout (query, symbol row, partner) for the two
%        winning registers, so callers and logs that expect an Alignment matrix keep working unchanged.

% WARNING STATE: the three warnings dGMismatch can raise are silenced for the duration of this
% function and the previous state is restored on exit. Enumerating registers means deliberately
% scoring structures that will be discarded, and most of them are badly mismatched, so a warning
% that is informative in a scored structure is noise here.
%   dGMismatch:NoPairedPosition is the one this enumeration actually provokes. A register with no
%     Watson-Crick pair at all initiates no helix, which dGMismatch reports because it is unexpected
%     in a structure someone meant to score. In an exhaustive sweep it arrives at rate (3/4)^L and
%     carries no information; dGMismatch's own message names this silencing as the remedy.
%   dGMismatch:MixedMismatchRoute fires when an isolated mismatch has one Table 2 step with measured
%     enthalpy and one without. With mismatch_nn_params.csv complete at 48 of 48 measured entries it
%     cannot fire, and it is silenced here only so that an incomplete table does not turn a register
%     sweep into a wall of text.
%   dGMismatch:UnverifiedTerminalOrientation fires only when ctx.UseViennaTerminal is set, which
%     this function never sets. Silenced for the same reason.
% All three are once-per-session warnings, so none of them can flood on its own. What they would do
% is fire on the first discarded register of a run and be attributed to a structure nobody scored.
% They stay active in LooperMismatch, MisprimerMismatch and InternalLoopBridge, where every scored
% structure is one the caller meant to score and any of the three is worth seeing.

if (nargin < 7)||isempty(minrun), minrun = 3; end
if (nargin < 8)||isempty(tmplmin), tmplmin = 0; end

warnstate = warning;
warnrestore = onCleanup(@() warning(warnstate)); %#ok<NASGU>
warning('off', 'dGMismatch:NoPairedPosition');
warning('off', 'dGMismatch:MixedMismatchRoute');
warning('off', 'dGMismatch:UnverifiedTerminalOrientation');

seqa = char(seqa); seqb = char(seqb);

% HOMODIMER OR HETERODIMER, DECIDED HERE AND FORWARDED TO BestSubRegister.
% The symmetry correction of Equation 3 applies when two IDENTICAL strands associate and not
% otherwise. BestSubRegister cannot decide it for itself, because the same function also scores
% primer-against-template registers for MisprimerMismatch, where the two strands are distinct
% molecules and the twofold rotational degeneracy does not exist. This is the only place that knows
% which case it is. A self-complementary register inside a self-dimer genuinely carries that
% degeneracy, so the term is owed there and is worth about 0.46 kcal/mol at 58 C.
isHomo = isequal(upper(seqa), upper(seqb));
na = numel(seqa); nb = numel(seqb);
seqbrc = char(seqrcomplement(seqb));

dGdimer = 0; dGend = 0;
dimeralign = repmat(' ', 3, 0); endalign = repmat(' ', 3, 0);   % char.empty is MATLAB-only; zeros takes numeric classes only (G3)
bestdG = Inf; bestenddG = Inf;

% Register offset: seqa index p is aligned with seqbrc index p - off. Both strands are scanned from
% complete 5' overhang to complete 3' overhang. For a homodimer, off and -off are genuinely different
% duplexes (the pairing map is self-symmetric under molecule swap for a fixed off), so the full range
% is required in that case too.
for off = -(nb - 1):(na - 1)
    lo = max(1, 1 + off); hi = min(na, nb + off);
    if (hi - lo + 1) < minrun, continue; end
    idxa = lo:hi;
    basa = seqa(idxa);
    basb = seqbrc((lo - off):(hi - off));
    ispair = (basa == basb);

    % Minimise over sub-registers instead of trimming once to the outermost identity. See the note
    % at the top of BestSubRegister.m for the measurement that forced this: on the reverse primer of
    % the 17-10-2026 run the single trim reported +4.544 kcal/mol for a duplex worth +0.245.   (G3)
    sub = BestSubRegister(basa, basb, ispair, minrun, monv, biv, dNTP, tann, isHomo);

    if isfinite(sub.dGfree)&&(sub.dGfree < bestdG)
        k = sub.iFree:sub.jFree;
        bestdG = sub.dGfree;
        dimeralign = [basa(k); PairSymbols(ispair(k)); basb(k)];
    end
    % Extendability needs the duplex to REACH the 3' terminus of seqa, so the right edge stays
    % pinned at the outermost 3' pair rather than being free to move inward. idxa(sub.jRight) is
    % that pair's coordinate in seqa, so the terminus is paired exactly when it equals na.
    % Partner nucleotides 5' of the base paired to the terminus of seqa: m - 1 with m = nb + 1 - na + off
    tmplbeyond = nb - na + off;
    if isfinite(sub.dGright)&&(idxa(sub.jRight) == na)&&(tmplbeyond >= tmplmin)&&(sub.dGright < bestenddG)
        k = sub.iRight:sub.jRight;
        bestenddG = sub.dGright;
        endalign = [basa(k); PairSymbols(ispair(k)); basb(k)];
    end
end

if isfinite(bestdG), dGdimer = bestdG; end
if isfinite(bestenddG), dGend = bestenddG; end
end

%% ____________________________________________________________________________________________________________
function runmax = MaxRun(mask)
% Longest run of true values in a logical vector.
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
% Middle row of a swalign-style alignment block.
sym = repmat(' ', 1, numel(mask));
sym(mask) = '|';
end
