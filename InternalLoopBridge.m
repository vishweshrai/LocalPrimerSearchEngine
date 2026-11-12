function [dGil, found] = InternalLoopBridge(alstruct, primer, foot, monv, biv, dNTP, tann, endTol, N)
%% INTERNALLOOPBRIDGE  Two-sided internal loops bridging two colinear alignment arms
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 09-11-2026 00:48 IST
%
% Revision history is NOT kept here. See InternalLoopBridge_<YYYYMMDD>.txt beside this file.
% ____________________________________________________________________________________________________________
% WHAT THIS FUNCTION DETECTS
%
% A genuine DNA internal loop carries unpaired bases on BOTH strands at once, between two anchored
% arms. A single pairwise alignment cannot represent one: a gap column always has a real base in
% only one row, so it encodes a bulge and never a two-sided loop. The structure is recovered here by
% comparing pairs of registers already held in alstruct. Two registers that are colinear and
% non-overlapping on both the primer and the target coordinate, with unpaired bases on both sides in
% between, are exactly a bridged internal loop.
%
% HOW IT IS SCORED
%
% Equation 12 of SantaLucia and Hicks, with N = n1 + n2 the total unpaired nucleotides:
%
%   dG(loop total) = dG(internal loop of N) + dG(asymmetry)
%                  + dG(left terminal mismatch) + dG(right terminal mismatch)
%
% with dG(asymmetry) = |n1 - n2| * 0.3 kcal/mol. The two arms are scored by dGMismatch with the
% interior-facing ends suppressed, so neither arm charges a helix initiation into the loop.
%
% THE TWO CLOSURES ARE SEQUENCE SPECIFIC
%
% Each closure is the measured stack of the first unpaired pair on the Watson-Crick pair that closes
% that side, taken from dGTerminalMismatch. The 48 measured entries span -1.19 kcal/mol to roughly
% zero, so a sequence-averaged stand-in misstates any individual loop by up to about 0.55 kcal/mol
% per closure. The bases are read from the primer and the footprint at the coordinates alstruct
% carries, which is why both sequences are arguments. Should either base fall outside A C G T the
% averaged value N.TM.MeanG37 stands in for that one closure and nothing else.
%
% ORIENTATION. The motif table is indexed i1 i2 / j1 j2, meaning 5'-i1 i2-3' over 3'-j1 j2-5', with
% i1.j1 the closing pair. The primer runs 5' to 3' with increasing index and the base physically
% opposing primer position k is the COMPLEMENT of foot(k). At the LEFT closure the loop opens toward
% increasing index, so the primer supplies i1 i2. At the RIGHT closure it opens toward decreasing
% index, which along the primer is 3' to 5', so the two strands exchange roles and the complemented
% footprint supplies i1 i2 instead.
%
% THE RIGHT ARM MUST REACH THE PRIMER 3' TERMINUS
%
% A bridge is a reason to reject a candidate only when the structure it describes could actually be
% extended. A loop whose right arm stops short of the primer 3' end leaves an unpaired flap carrying
% the 3'-OH, which no polymerase without a 3' exonuclease can remove, so the site is catalytically
% dead however tightly it binds. The arm is therefore required to reach within endTol of the
% terminus, which is the same test the anchored register path applies and the same enzyme property:
% endTol is 1 for Taq and is widened to the primer length only when a proofreading exonuclease is
% declared.
%
% TARGET COORDINATES ARE UNWRAPPED, AND TOTAL LOOP SIZE IS BOUNDED AT TABLE 4 (Ntot <= 30)
%
% The caller leaves target coordinates unwrapped so that loop size stays a plain coordinate
% difference across the origin of a circular template. Indices into foot are therefore wrapped here
% before use. Because alstruct holds the top registers across the whole circle, two colinear hits
% can sit kilobases apart on the template. Bounding Ntot = n1 + n2 at 30 nucleotides (the domain of
% Table 4 in SantaLucia and Hicks 2004) prevents unrelated distant hits from being bridged into
% phantom loops that would otherwise increment the anchored register counts in MisprimerMismatch.
%
% REFERENCES
%   SantaLucia J Jr, Hicks D. The thermodynamics of DNA structural motifs.
%       Annu Rev Biophys Biomol Struct 2004;33:415-40. doi:10.1146/annurev.biophys.32.110601.141800
%   Peyret N, Seneviratne PA, Allawi HT, SantaLucia J Jr. Nearest-neighbor thermodynamics and NMR of
%       DNA sequences with internal A.A, C.C, G.G and T.T mismatches.
%       Biochemistry 1999;38:3468-77. doi:10.1021/bi9825091
% ____________________________________________________________________________________________________________
% INPUTS   alstruct  localalign-shaped struct: .Alignment cell array, .Start and .Stop Nx2 with
%                    column 1 the primer coordinate and column 2 the target coordinate
%          primer    char row, the primer, 5' to 3'
%          foot      char row, the target strand that was scanned, in the same frame as alstruct
%          monv, biv, dNTP, tann   as in dGPrim.m, concentrations in M and temperature in C
%          endTol    bases of the primer 3' end allowed outside the bridged duplex
%          N         optional NNParams struct
% OUTPUTS  dGil      dG of the most stable bridged loop, kcal/mol; 0 when none was found
%          found     true only when a loop was actually scored. Callers MUST test this rather than
%                    testing dGil against a threshold, because the 0 sentinel would otherwise fail
%                    every primer whenever the acceptance threshold is positive
% ____________________________________________________________________________________________________________

if (nargin < 9)||isempty(N), N = NNParams(); end
if (nargin < 8)||isempty(endTol), endTol = 1; end

asym_coeff = 0.3;              % kcal/mol per nt of loop asymmetry, Equation 12
dGil  = Inf;
found = false;

if ~isfield(alstruct, 'Alignment')||~iscell(alstruct.Alignment)
    dGil = 0; return
end
nA = numel(alstruct.Alignment);
if nA < 2
    dGil = 0; return                     % an internal loop is physically impossible
end

haveSeq = (nargin >= 3) && ~isempty(primer) && ~isempty(foot);
if haveSeq
    primer = char(primer); foot = char(foot);
    plen = numel(primer); flen = numel(foot);
else
    plen = NaN;
end

for i = 1:nA
    for j = 1:nA
        if i == j, continue; end
        % Hit i must precede hit j on BOTH the primer and the target coordinate
        if ~(alstruct.Stop(i, 1) < alstruct.Start(j, 1) && alstruct.Stop(i, 2) < alstruct.Start(j, 2))
            continue
        end
        n1 = alstruct.Start(j,1) - alstruct.Stop(i,1) - 1;      % unpaired primer bases
        n2 = alstruct.Start(j,2) - alstruct.Stop(i,2) - 1;      % unpaired target bases
        if (n1 <= 0)||(n2 <= 0)
            continue        % one-sided, so a bulge rather than a loop, and handled elsewhere
        end
        Ntot = n1 + n2;
        if (Ntot == 2)||(Ntot > 30)
            continue        % Ntot == 2 is a single mismatch (Table 2); Ntot > 30 exceeds Table 4
        end

        % THE RIGHT ARM MUST REACH THE 3' TERMINUS. Without this the bridge can reject a candidate
        % on a structure that leaves the 3'-OH on an unpaired flap, which cannot prime.
        if haveSeq && (alstruct.Stop(j, 1) < plen - endTol)
            continue
        end

        leftAlign  = alstruct.Alignment{i};
        rightAlign = alstruct.Alignment{j};
        if any(leftAlign(1,:) == '-')||any(leftAlign(3,:) == '-')|| ...
                any(rightAlign(1,:) == '-')||any(rightAlign(3,:) == '-')
            continue        % gapped arm; seqcomplement is undefined on '-'
        end
        if (size(leftAlign, 2) < 2)||(size(rightAlign, 2) < 2)
            continue        % too short for a nearest-neighbour step
        end
        leftSym  = strrep(leftAlign(2,:),  '|', 'X');
        rightSym = strrep(rightAlign(2,:), '|', 'X');
        leftMatches  = regexp(leftSym,  'X+', 'match');
        rightMatches = regexp(rightSym, 'X+', 'match');
        if isempty(leftMatches)||isempty(rightMatches)|| ...
                max(cellfun(@length, leftMatches))  < 3|| ...
                max(cellfun(@length, rightMatches)) < 3
            continue        % neither arm carries a three-base seed anchor
        end

        seqL  = char(leftAlign(1,:));
        seq2L = char(seqcomplement(leftAlign(3,:)));
        seqR  = char(rightAlign(1,:));
        seq2R = char(seqcomplement(rightAlign(3,:)));

        dGarmL = dGMismatch(seqL, seq2L, monv, biv, dNTP, tann, false, true,  false);
        dGarmR = dGMismatch(seqR, seq2R, monv, biv, dNTP, tann, false, false, true);

        % ---- the two closures ----
        if haveSeq
            pL = alstruct.Stop(i,1);  tL = alstruct.Stop(i,2);
            pR = alstruct.Start(j,1); tR = alstruct.Start(j,2);
            dGtmL = ClosureDG(primer, foot, flen, pL, tL, pL + 1, tL + 1, false, tann, N);
            dGtmR = ClosureDG(primer, foot, flen, pR, tR, pR - 1, tR - 1, true,  tann, N);
        else
            dGtmL = N.TM.MeanG37;
            dGtmR = N.TM.MeanG37;
        end

        dGtotal = dGarmL + dGarmR + dGLoop('internal', Ntot, tann, N) ...
                + abs(n1 - n2)*asym_coeff + dGtmL + dGtmR;

        if dGtotal < dGil
            dGil  = dGtotal;
            found = true;
        end
    end
end

if isinf(dGil), dGil = 0; end
end

%% ____________________________________________________________________________________________________________
function g = ClosureDG(primer, foot, flen, pPair, tPair, pMism, tMism, flipStrands, tann, N)
%% CLOSUREDG  One terminal-mismatch closure of a bridged internal loop
% pPair and tPair index the closing Watson-Crick pair, pMism and tMism the first unpaired pair
% beyond it. flipStrands is true at the right-hand closure, where the loop opens toward decreasing
% index and the two strands therefore exchange the i and j roles of the motif table.
g = N.TM.MeanG37;
if (pPair < 1)||(pMism < 1)||(pPair > numel(primer))||(pMism > numel(primer)), return; end
i1 = primer(pPair);  i2 = primer(pMism);
j1 = Opposite(foot, flen, tPair);
j2 = Opposite(foot, flen, tMism);
if any(~ismember([i1 i2 j1 j2], 'ACGT')), return; end
if flipStrands
    [i1, j1] = deal(j1, i1);
    [i2, j2] = deal(j2, i2);
end
[dg, ~, ~, info] = dGTerminalMismatch(i1, i2, j1, j2, tann, N);
if info.Found, g = dg; end
end

%% ____________________________________________________________________________________________________________
function b = Opposite(foot, flen, k)
%% OPPOSITE  The base physically facing primer position k, with the target coordinate wrapped
% The caller leaves target coordinates unwrapped across the origin of a circular template, and the
% primer is identical to the strand being scanned, so the base opposing it is the COMPLEMENT of the
% footprint base rather than the footprint base itself.
idx = mod(k - 1, flen) + 1;
c = upper(foot(idx));
switch c
    case 'A', b = 'T';
    case 'T', b = 'A';
    case 'G', b = 'C';
    case 'C', b = 'G';
    otherwise, b = 'N';
end
end
