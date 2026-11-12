function dGm = dGmaxMismatch(Alignment, monv, biv, dNTP, tann)
%% DGMAXMISMATCH  Mismatch- and bulge-aware duplex free energy from a Smith-Waterman alignment
% Mismatch- and bulge-aware version of dGmax. Prior versions either
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 04-11-2026 23:10 IST
%       (a) extracted only bare '|' match runs, discarding mismatches, or
%       (b) scored mismatches but still split the alignment into separate, independently-scored
%           spans at every gap ('-'), discarding any single-strand bulge as if it were two
%           unrelated duplexes.

% This version treats a gap run in ONE row (the other row having real bases) as a BULGE
% (SantaLucia & Hicks 2004, Eq. 13, Table 4) and scores the ENTIRE bulge-bridged region as ONE
% physical duplex:
%               dG(total) = dG(left arm) + dG(bulge) + dG(right arm) + ...,
% with initiation counted once and terminal-AT penalties applied only at the true 5'/3' ends of the
% whole candidate, not at the internal cuts adjacent to a bulge (via dGMismatch's
% applyInit/applyTerm flags). A run of >=3 consecutive matches is still required SOMEWHERE in the
% merged candidate. True INTERNAL loops (unpaired bases on BOTH strands simultaneously) cannot occur
% in a standard pairwise local alignment's gap structure, so they are not auto-detected here;
% InternalLoopBridge.m handles that case.

% FRAME CONVENTION, stated once because three separate defects came from losing track of it: the
% caller aligns a query against the REVERSE COMPLEMENT of its partner, so identity between row 1 and
% row 3 means Watson-Crick pairing, and the base physically opposing alignbas(i) is
% comp(alignref(i)). Every seq2 handed to dGMismatch must therefore be seqcomplement(alignref), and
% every "is this position paired" test must be alignbas(i) == alignref(i), NOT
% comp(alignbas(i)) == alignref(i).

% CHANGES
%  1. The "no candidate found" fallback (dGm = 0) sat inside the numel(validSegIdx) >= 2 block. Every
%     caller invokes swalign with 'GapOpen', 1000, so gaps never appear, so there is always exactly
%     one segment and the fallback was unreachable. Any alignment without a run of three identities
%     returned Inf, which then propagated into dGff, dGrr and dGfr in the candidate logs. The reset
%     now runs unconditionally.
%  2. The single-base-bulge flanking stack passed alignref raw as seq2 while every other call site
%     passed seqcomplement(alignref). dGMismatch therefore computed the pairing test inverted and
%     scored the flanking stack as a mismatch. Now complemented.
%  3. closingAT_left/right tested comp.(leftBase) against leftPair, which in this frame is true
%     precisely when the closing bases are NOT paired. Now tests identity.
%  4. Segments shorter than two bases are skipped rather than handed to dGMismatch, where they would
%     divide by zero in the divalent salt term and apply both initiation terms to the same base.
%  5. The warning string was missing a space across its concatenation and had no identifier.
%
% The +0.5 AT-closing penalty is correctly sourced: Equation 13 reads "dG37(Loop total) =
% dG37(Bulge Loop of 1) + dG37(intervening NN) + closing AT penalty", and the worked example gives
% +4.0 - 1.28 + 0.5 = +3.22 for 5'-CAT-3'/5'-AG-3'. It applies per AT closing pair, which is why the
% two ends are tested independently.

AT_CLOSING_PENALTY = 0.5; % kcal/mol per AT pair closing a bulge, SantaLucia & Hicks Eq. 13

comp = struct('A','T','T','A','C','G','G','C');
alignsym = Alignment(2, :);
alignbas = char(Alignment(1, :)); % Primer/query bases, with any gaps
alignref = char(Alignment(3, :)); % Aligned reference bases, with any gaps
symchar = replace(alignsym, '|', 'X');
L = length(alignbas);
gapmask = (alignbas == '-')|(alignref == '-');
dGm = Inf;

% Candidate 1: Each gap-free span scored alone (fallback)
breaks = [0, find(gapmask), L + 1];
nSeg = numel(breaks) - 1;
segStart = zeros(1, nSeg); segEnd = zeros(1, nSeg); segValid = false(1, nSeg);
for k = 1:nSeg
    spanIdx = (breaks(k) + 1):(breaks(k + 1) - 1);
    if numel(spanIdx) < 2, continue; end
    segStart(k) = spanIdx(1); segEnd(k) = spanIdx(end);
    matches = regexp(symchar(spanIdx), 'X+', 'match');
    numocc = cellfun(@length, matches);
    if ~isempty(numocc) && max(numocc) >= 3
        segValid(k) = true;
        seq = alignbas(spanIdx);
        seq2 = char(seqcomplement(alignref(spanIdx)));
        dG = dGMismatch(seq, seq2, monv, biv, dNTP, tann);
        if dG < dGm, dGm = dG; end
    end
end

% Candidate 2: Bridge across gap runs as bulges, score as one duplex
validSegIdx = find(segEnd > 0); % Non-empty segments of at least two bases, in order
if numel(validSegIdx) >= 2
    dGbridge = 0;
    anyMatchAnchor = false;
    for j = 1:numel(validSegIdx)
        k = validSegIdx(j);
        spanIdx = segStart(k):segEnd(k);
        seq = alignbas(spanIdx);
        seq2 = char(seqcomplement(alignref(spanIdx)));
        applyInit = (j == 1);
        applyTerm5 = (j == 1);
        applyTerm3 = (j == numel(validSegIdx));
        dGbridge = dGbridge + dGMismatch(seq, seq2, monv, biv, dNTP, tann, applyInit, applyTerm5, applyTerm3);
        if segValid(k), anyMatchAnchor = true; end
        if j < numel(validSegIdx)
            nextk = validSegIdx(j+1);
            gapIdx = (segEnd(k)+1):(segStart(nextk)-1);
            bulgelen = numel(gapIdx);
            isPrimerBulge = all(alignref(gapIdx) == '-'); % Extra bases in alignbas = primer-side bulge
            isRefBulge = all(alignbas(gapIdx) == '-'); % Extra bases in alignref = target-side bulge
            if ~(isPrimerBulge||isRefBulge)
                warning('dGmaxMismatch:MixedGapPattern', ...
                    ['Mixed gap pattern at positions %d-%d is not a simple single-strand bulge; ' ...
                    'skipping the bulge-bridged candidate.'], gapIdx(1), gapIdx(end));
                dGbridge = NaN; break;
            end
            leftBase = alignbas(segEnd(k)); leftPair = alignref(segEnd(k));
            rightBase = alignbas(segStart(nextk)); rightPair = alignref(segStart(nextk));
            % Identity, not complementarity: see the frame convention note in the header.
            closingAT_left = (leftBase == leftPair) && any(leftBase == 'AT');
            closingAT_right = (rightBase == rightPair) && any(rightBase == 'AT');
            dGbulge = dGLoop('bulge', bulgelen, tann);
            if bulgelen == 1
                % A bulge of one retains the flanking nearest-neighbour stack (Eq. 13).
                flankTop = [leftBase rightBase];
                flankBot = char(seqcomplement([leftPair rightPair]));
                dGbulge = dGbulge + dGMismatch(flankTop, flankBot, monv, biv, dNTP, tann, false, false, false);
            end
            if closingAT_left, dGbulge = dGbulge + AT_CLOSING_PENALTY; end
            if closingAT_right, dGbulge = dGbulge + AT_CLOSING_PENALTY; end
            dGbridge = dGbridge + dGbulge;
        end
    end
    if ~isnan(dGbridge) && anyMatchAnchor && dGbridge < dGm
        dGm = dGbridge;
    end
end

% Fallback to 0 when no candidate with a run of three matches was found at all. This sits outside the
% bridging block, which is the whole point: with gaps forbidden there is never more than one segment.
if isinf(dGm)
    dGm = 0;
end
end
