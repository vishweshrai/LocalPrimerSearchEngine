function out = EnzymeBracket(runNone, runExo, opts)
%% Enzyme Blend Bracket for the aPCR Primer Screening Engine
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 12-11-2026 19:53 IST
%
% Revision history is not kept here. See EnzymeBracket_<YYYYMMDD>.txt beside this file.
% ____________________________________________________________________________________________________________
% What This Script Does
%
% NEB LongAmp Taq is a blend of Taq and Deep Vent DNA polymerases whose proportions are not published,
% so the fraction phi of 3' end engagements served by the proofreading enzyme is unknown. The engine
% runs in one channel at a time: ProofreadingClass 'none' judges every gate and the collapse stage
% with Taq, and 'weak' (or 'strong') judges them with the proofreading enzyme. This script brackets
% the unknown phi from two completed runs, one in each class.
%
% A pair belongs to the robust set when it survives both runs, because a pair that survives the
% Taq channel and the exonuclease channel separately has no reject-tier collapse in either. The
% collapse term of the Score is then recomputed on a grid of phi from 0 (Taq) to 1 (exonuclease),
% with the extendable fraction of each strand taken as
%
%       pColl(phi) = max over the two strands of ((1 - phi)*pExtTaq + phi*pExtExo)
%
% which is exact when each engagement is served by one enzyme molecule. A strand in the flag tier in
% either run adds the flag increment for 0 < phi < 1. The other Score terms are taken from the run
% in the class 'none', less its own collapse term. The output gives, for every robust pair, its Score
% and rank at each phi, the best and the worst rank over the grid and whether it stays inside the top
% set at every phi. A pair that does so is a choice that does not depend on the unknown blend.
%
% Both runs need the settings CollapseLog, CollapseGate, CollapseTol, CollapseFlagDecades, CollapseCap
% and WtCollapse to be the same, since only the proofreading class should differ between them.
% With one run only (runExo empty) the survival test is skipped and the bracket uses the Taq and Exo
% fractions that the single run logs, which are computed in both channels whatever the class.
% ____________________________________________________________________________________________________________
% Usage
%   out = EnzymeBracket(runNone, runExo, opts)
%
%   runNone  Pair_Log table of the run in class 'none', or a path to its xlsx (sheet Pair_Log) or csv
%   runExo   the same for the run in class 'weak' or 'strong'; empty for a single-run bracket
%   opts     optional struct. Fields and defaults: Tol 1E-3, FlagDecades 1.0, Cap 3.0, WtCollapse 5/3,
%            PhiGrid 0:0.1:1, Top 10
%
% Output
%   out.Pairs      n by 4: forward position, forward length, reverse position, reverse length
%   out.Phi        the grid
%   out.Score      n by numel(Phi), Score at each phi
%   out.Rank       n by numel(Phi), rank at each phi, 1 is best
%   out.RankBest, out.RankWorst, out.InTop  per pair over the grid
%   out.Dropped    number of pairs of the first run that did not survive the second
% ____________________________________________________________________________________________________________
% References
%   New England Biolabs. LongAmp Taq DNA Polymerase, M0323. Product page, accessed 12-11-2026:
%       https://www.neb.com/en-us/products/m0323-longamp-taq-dna-polymerase
% ____________________________________________________________________________________________________________

if nargin < 2, runExo = []; end
if nargin < 3||isempty(opts), opts = struct(); end
opts = FillOptions(opts);

A = LoadRun(runNone);
selA = strcmp(cellstr(string(A.Status)), 'Selected');
keyA = PairKeys(A);
keyA = keyA(selA);
dropped = 0;
keep = true(size(keyA));
flagB = false(size(keyA));
if ~isempty(runExo)
    B = LoadRun(runExo);
    selB = strcmp(cellstr(string(B.Status)), 'Selected');
    keyB = PairKeys(B);
    [inB, loc] = ismember(keyA, keyB(selB));
    dropped = nnz(~inB);
    keep = inB;
    idxB = find(selB);
    flagB = false(size(keyA));
    flagB(inB) = HasFlag(B, idxB(loc(inB)));
end
idxA = find(selA);
idxA = idxA(keep);
flagA = HasFlag(A, idxA);
flagAny = flagA|flagB(keep);

pTaqP = Col(A, 'cpPExtTaq', idxA); pExoP = Col(A, 'cpPExtExo', idxA);
pTaqM = Col(A, 'cmPExtTaq', idxA); pExoM = Col(A, 'cmPExtExo', idxA);
base = Col(A, 'Score', idxA);
if ismember('CollapseDecades', A.Properties.VariableNames)
    base = base - opts.WtCollapse*Col(A, 'CollapseDecades', idxA);
end

phi = opts.PhiGrid(:).';
n = numel(idxA);
score = zeros(n, numel(phi));
for g = 1:numel(phi)
    pc = max((1 - phi(g))*pTaqP + phi(g)*pExoP, (1 - phi(g))*pTaqM + phi(g)*pExoM);
    if phi(g) <= 0, fl = flagA(:); elseif phi(g) >= 1, fl = flagB(keep).'; fl = fl(:); else, fl = flagAny(:); end
    if isempty(runExo), fl = flagA(:); end
    for k = 1:n
        dec = ScaffoldCollapseCheck('penalty', pc(k), fl(k), opts.Tol, opts.FlagDecades, opts.Cap);
        score(k, g) = base(k) + opts.WtCollapse*dec;
    end
end
rank = zeros(size(score));
for g = 1:numel(phi)
    [~, ord] = sort(score(:, g));
    rank(ord, g) = (1:n).';
end
out = struct('Pairs', [Col(A, 'fPos', idxA) Col(A, 'flen', idxA) Col(A, 'rPos', idxA) Col(A, 'rlen', idxA)], 'Phi', phi, 'Score', score, 'Rank', rank, ...
    'RankBest', min(rank, [], 2), 'RankWorst', max(rank, [], 2), 'InTop', (max(rank, [], 2) <= opts.Top), 'Dropped', dropped);

%% Print
fprintf('\nEnzymeBracket: %d pairs in the robust set, %d dropped by the second run\n', n, dropped);
[~, ord] = sort(out.RankWorst);
shown = min(n, opts.Top);
fprintf('%-10s %-10s %10s %10s %10s %6s %6s\n', 'forward', 'reverse', 'Score 0', 'Score mid', 'Score 1', 'best', 'worst');
mid = ceil(numel(phi)/2);
for q = 1:shown
    k = ord(q);
    fprintf('%-10s %-10s %10.3f %10.3f %10.3f %6d %6d\n', sprintf('%d/%d', out.Pairs(k, 1), out.Pairs(k, 2)), sprintf('%d/%d', out.Pairs(k, 3), out.Pairs(k, 4)), score(k, 1), score(k, mid), score(k, end), ...
        out.RankBest(k), out.RankWorst(k));
end
fprintf('%d of %d pairs stay in the top %d at every phi\n', nnz(out.InTop), n, opts.Top);
end

%% Fill Options
function o = FillOptions(o)
d = struct('Tol', 1E-3, 'FlagDecades', 1.0, 'Cap', 3.0, 'WtCollapse', 5.0/3.0, 'PhiGrid', 0:0.1:1, 'Top', 10);
f = fieldnames(d);
for q = 1:numel(f)
    if ~isfield(o, f{q}), o.(f{q}) = d.(f{q}); end
end
end

%% Load Run
function T = LoadRun(src)
if ischar(src)||(isstring(src)&&isscalar(src))
    src = char(src);
    if endsWith(lower(src), '.xlsx'), T = readtable(src, 'Sheet', 'Pair_Log'); else, T = readtable(src); end
else
    T = src;
end
end

%% Pair Keys
function k = PairKeys(T)
k = arrayfun(@(a, b, c, d) sprintf('%d_%d_%d_%d', a, b, c, d), Col(T, 'fPos', 'all'), Col(T, 'flen', 'all'), ...
    Col(T, 'rPos', 'all'), Col(T, 'rlen', 'all'), 'UniformOutput', false);
end

%% Column Access
function v = Col(T, name, idx)
v = T.(name);
if iscell(v), v = str2double(v); end
v = double(v(:));
if ~ischar(idx), v = v(idx); end
end

%% Flag Test
function f = HasFlag(T, idx)
f = false(numel(idx), 1);
vn = T.Properties.VariableNames;
if all(ismember({'cpTierDecl', 'cmTierDecl'}, vn))
    tp = cellstr(string(T.cpTierDecl)); tm = cellstr(string(T.cmTierDecl));
    f = (strcmp(tp(idx), 'FLAG')|strcmp(tm(idx), 'FLAG'));
    f = f(:);
end
end
