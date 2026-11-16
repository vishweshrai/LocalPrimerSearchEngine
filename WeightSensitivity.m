function out = WeightSensitivity(run, opts)
%% Score Weight Sensitivity for the aPCR Primer Screening Engine
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 12-11-2026 19:53 IST
%
% Revision history is not kept here. See WeightSensitivity_<YYYYMMDD>.txt beside this file.
% ____________________________________________________________________________________________________________
% What This Script Does
%
% Shows how far the ranking of the accepted pairs depends on the Score weights, using a finished run
% and no new search. The engine logs every term as it enters the Score (columns tDelta, tTmA,
% tMargin, tMisp, tEns, tCollapse, tAccess, t3End, tProf, tFlank, tHairpin and t3Base of Pair_Log),
% and Score is their exact sum, so scaling a weight is scaling a column.
%
% Three analyses
%   1. One weight at a time. Each term is multiplied by each factor in opts.Factors (default 0, 0.5
%      and 2). For each case the script reports the Spearman correlation of the new ranking with the
%      base ranking, how many of the base top-K pairs remain in the top K, and whether the first pair
%      changes. A term whose removal or doubling leaves these unchanged does not steer the choice.
%   2. The collapse tolerance. CollapseTol is not a published number, so its column is recomputed from
%      the logged extendable fractions (Collapse_Plus_pExtTot and Collapse_Minus_pExtTot, here
%      cpPExtTot and cmPExtTot) at each tolerance in opts.Tols, and the same measures are reported.
%   3. All weights at once. Every term is multiplied by a log-normal factor, ln-scale SD opts.Sigma
%      (default ln 2/2, so a factor of 2 sits at two SD), for opts.Draws draws. The script reports how
%      often each pair is first and how often it sits in the top K. A pair that is first in most
%      draws is a choice that does not hinge on any particular set of weights.
%
% Usage
%   out = WeightSensitivity(run, opts)
%
%   run   Pair_Log table of one run, or a path to its xlsx (sheet Pair_Log) or csv
%   opts  optional struct. Fields and defaults: Factors [0 0.5 2], Top 5, Tols [1E-4 1E-3 1E-2 1E-1],
%         Tol 1E-3, FlagDecades 1.0, Cap 3.0, WtCollapse 5/3, Sigma log(2)/2, Draws 1000, Seed 1
%
% Output
%   out.Terms, out.Factors, out.Rho, out.TopKept, out.FirstChanged   one weight at a time
%   out.Tols, out.TolRho, out.TolTopKept, out.TolFirstChanged         collapse tolerance
%   out.FirstShare, out.TopShare                                      all weights, per pair
%   out.Pairs   n by 4: forward position, forward length, reverse position, reverse length
% ____________________________________________________________________________________________________________

if nargin < 2||isempty(opts), opts = struct(); end
d = struct('Factors', [0 0.5 2], 'Top', 5, 'Tols', [1E-4 1E-3 1E-2 1E-1], 'Tol', 1E-3, 'FlagDecades', 1.0, 'Cap', 3.0, ...
    'WtCollapse', 5.0/3.0, 'Sigma', log(2)/2, 'Draws', 1000, 'Seed', 1);
f = fieldnames(d);
for q = 1:numel(f)
    if ~isfield(opts, f{q}), opts.(f{q}) = d.(f{q}); end
end

T = LoadRun(run);
vn = T.Properties.VariableNames;
terms = {'tDelta', 'tTmA', 'tMargin', 'tMisp', 'tEns', 'tCollapse', 'tAccess', 't3End', 'tProf', 'tFlank', 'tHairpin', 't3Base'};
if ~all(ismember(terms, vn))
    error('WeightSensitivity:NoTerms', 'Pair_Log has no per-term columns. Re-run the search with the current engine.');
end
sel = find(strcmp(cellstr(string(T.Status)), 'Selected'));
n = numel(sel);
M = zeros(n, numel(terms));
for q = 1:numel(terms)
    M(:, q) = Col(T, terms{q}, sel);
end
base = sum(M, 2);
top = min(opts.Top, n);
r0 = RankOf(base);
[~, first0] = min(base);
set0 = find(r0 <= top);

%% One Weight at a Time
nf = numel(opts.Factors);
rho = zeros(numel(terms), nf); kept = rho; chg = rho;
for q = 1:numel(terms)
    for g = 1:nf
        s = base + (opts.Factors(g) - 1)*M(:, q);
        [rho(q, g), kept(q, g), chg(q, g)] = Compare(s, r0, first0, set0, top);
    end
end

%% Collapse Tolerance
pColl = max(Col(T, 'cpPExtTot', sel), Col(T, 'cmPExtTot', sel));
flag = false(n, 1);
if all(ismember({'cpTierDecl', 'cmTierDecl'}, vn))
    flag = (strcmp(cellstr(string(T.cpTierDecl(sel))), 'FLAG')|strcmp(cellstr(string(T.cmTierDecl(sel))), 'FLAG'));
end
iC = find(strcmp(terms, 'tCollapse'));
nt = numel(opts.Tols);
tolRho = zeros(nt, 1); tolKept = tolRho; tolChg = tolRho;
for g = 1:nt
    s = base - M(:, iC);
    for k = 1:n
        s(k) = s(k) + opts.WtCollapse*ScaffoldCollapseCheck('penalty', pColl(k), flag(k), opts.Tols(g), opts.FlagDecades, opts.Cap);
    end
    [tolRho(g), tolKept(g), tolChg(g)] = Compare(s, r0, first0, set0, top);
end

%% All Weights at Once
rand('seed', opts.Seed); randn('seed', opts.Seed);
firstCount = zeros(n, 1); topCount = zeros(n, 1);
for k = 1:opts.Draws
    w = exp(opts.Sigma*randn(1, numel(terms)));
    s = M*w(:);
    [~, i1] = min(s);
    firstCount(i1) = firstCount(i1) + 1;
    r = RankOf(s);
    topCount(r <= top) = topCount(r <= top) + 1;
end

out = struct('Terms', {terms}, 'Factors', opts.Factors, 'Rho', rho, 'TopKept', kept, 'FirstChanged', chg, ...
    'Tols', opts.Tols, 'TolRho', tolRho, 'TolTopKept', tolKept, 'TolFirstChanged', tolChg, ...
    'FirstShare', firstCount/opts.Draws, 'TopShare', topCount/opts.Draws, ...
    'Pairs', [Col(T, 'fPos', sel) Col(T, 'flen', sel) Col(T, 'rPos', sel) Col(T, 'rlen', sel)]);

%% Print
fprintf('\nWeightSensitivity: %d accepted pairs, top %d\n', n, top);
fprintf('One weight at a time (Spearman rho / top-%d kept / first pair changed)\n', top);
fprintf('%-10s', 'term');
for g = 1:nf, fprintf('  x%-14g', opts.Factors(g)); end
fprintf('\n');
for q = 1:numel(terms)
    fprintf('%-10s', terms{q});
    for g = 1:nf, fprintf('  %5.3f/%d/%d      ', rho(q, g), kept(q, g), chg(q, g)); end
    fprintf('\n');
end
fprintf('Collapse tolerance (Spearman rho / top-%d kept / first pair changed)\n', top);
for g = 1:nt, fprintf('  Tol %-8.0e %5.3f/%d/%d\n', opts.Tols(g), tolRho(g), tolKept(g), tolChg(g)); end
[~, ord] = sort(out.FirstShare, 'descend');
fprintf('All weights at once, %d draws: pairs most often first\n', opts.Draws);
for q = 1:min(5, n)
    k = ord(q);
    fprintf('  %d/%d  %d/%d  first %5.1f%%  top-%d %5.1f%%  base rank %d\n', out.Pairs(k, 1), out.Pairs(k, 2), out.Pairs(k, 3), out.Pairs(k, 4), ...
        100*out.FirstShare(k), top, 100*out.TopShare(k), r0(k));
end
end

%% Compare Rankings
function [rho, kept, chg] = Compare(s, r0, first0, set0, top)
r = RankOf(s);
c = corrcoef(r(:), r0(:));
rho = c(1, 2);
kept = nnz(r(set0) <= top);
[~, i1] = min(s);
chg = double(i1 ~= first0);
end

%% Rank of Each Element
function r = RankOf(s)
[~, ord] = sort(s(:));
r = zeros(numel(s), 1);
r(ord) = (1:numel(s)).';
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

%% Column Access
function v = Col(T, name, idx)
v = T.(name);
if iscell(v), v = str2double(v); end
v = double(v(:));
v = v(idx);
end
