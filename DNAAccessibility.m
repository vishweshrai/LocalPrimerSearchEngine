function [siteTable, lunp, posMap] = DNAAccessibility(seq, sites, opts)
%% DNAAccessibility  RNAplfold accessibility of named sites on a DNA strand
% Scripted for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore (INDIA).
% 04-11-2026 23:10 IST
% ____________________________________________________________________________________________________________
% Companion to ViennaRNAhandler.m, with three differences that matter for DNA rather than mRNA.
%
%  1. --paramFile IS MANDATORY. Without it RNAplfold uses the RNA energy model, which allows G-U
%     pairs. Those do not exist in DNA, so every probability would be wrong in a direction that
%     makes structure look more stable than it is. dna_mathews2004.par ships with the package.
%
%  2. Circular templates are handled by triplication. RNAplfold has no circular mode, so a site near
%     position 1 of the M13 circle would be scored with no 5' context at all. The sequence is run as
%     three tandem copies and the middle copy is read, which gives every position a full window of
%     context on both sides. posMap carries the offset.
%
%  3. Two accessibility numbers are returned per site, not one.
%       pFull  probability that the ENTIRE footprint is unpaired. The right question for a short
%              probe like a molecular beacon loop, which has to displace whatever is there using
%              only the ~16 kcal/mol its own duplex is worth.
%       pNuc   the best probability over any nucleation-length sub-window inside the footprint. The
%              right question for a long oligo. A 44-mer guide forms a duplex worth roughly
%              -60 kcal/mol, so it does not need the whole site open, it needs a toehold to
%              nucleate on and then it zips. Asking for a 44 nt stretch to be wholly unpaired would
%              reject every site on a structured template.
%
% Free energy of opening is dGopen = -R*T*ln(p), which can be added directly to a binding free
% energy computed by dGMismatch to give the effective free energy of binding a structured target.
%
% SALT CAVEAT, and it is not small: ViennaRNA's --salt models monovalent cations only. There is no
% magnesium term. Magnesium stabilises secondary structure, so predictions run optimistic about
% accessibility. opts.Salt can be set to a monovalent-equivalent to compensate; see the note on
% opts.Salt below.
% ____________________________________________________________________________________________________________
%
% seq    char row vector, uppercase ACGT.
% sites  struct array or table with fields/variables Label, Start, Stop (1-based, inclusive,
%        positions on seq; Start may exceed Stop for a circular site that wraps).
% opts   struct, all fields optional:
%        .W          sliding window, default min(240, numel(seq))
%        .L          maximum pairing span, default round(0.7*W)
%        .uMax       longest stretch to tabulate, default max site length capped at 60
%        .Temp       degrees C, default 55
%        .Salt       molar monovalent, default 0.04. For a PCR buffer the von Ahsen rule of thumb
%                    gives a monovalent equivalent of [Na+] + 120*sqrt([Mg2+] - [dNTP]) in mM, so
%                    40 mM Na+ with 2.5 mM Mg2+ and 0.8 mM dNTP comes to about 0.196 M. That rule
%                    was derived for duplex melting, not for folding, so treat it as a sensitivity
%                    test rather than a correction: run the analysis at both 0.04 and 0.196 and see
%                    whether any conclusion moves.
%        .Circular   logical, default false
%        .NucLen     nucleation window for pNuc, default 12
%        .ParamFile  full path to dna_mathews2004.par; searched for if omitted
%        .OutDir     working directory for the RNAplfold output, default a tempdir folder
%        .Header     FASTA header, default 'DNAstrand'
%        .KeepFiles  logical, default true
%
% siteTable  one row per site: Label, Start, Stop, Len, pFull, dGopenFull, pNuc, dGopenNuc, NucStop
% lunp       raw RNAplfold table, rows are positions in the RUN sequence, column u+1 is
%            P(stretch of length u ending at that row is unpaired)
% posMap     function handle mapping a position on seq to its row in lunp

R = 8.314E-3/4.184;   % kcal/mol/K, same constant the rest of the suite uses

%% ---- input handling ----------------------------------------------------------------------------
seq = upper(char(seq));
if any(~ismember(seq, 'ACGT'))
    error('DNAAccessibility:NonACGT', 'seq carries characters outside {A,C,G,T}.');
end
N = numel(seq);
if istable(sites), sites = table2struct(sites); end

if nargin < 3 || isempty(opts), opts = struct(); end
def = struct('W', min(240, N), 'L', [], 'uMax', [], 'Temp', 55, 'Salt', 0.04, ...
    'Circular', false, 'NucLen', 12, 'ParamFile', '', 'OutDir', '', ...
    'Header', 'DNAstrand', 'KeepFiles', true);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts, fn{k}) || isempty(opts.(fn{k})), opts.(fn{k}) = def.(fn{k}); end
end
if isempty(opts.L), opts.L = round(0.7*opts.W); end
if opts.L > opts.W
    error('DNAAccessibility:BadWindow', 'L (%d) must not exceed W (%d).', opts.L, opts.W);
end

% Site lengths, wrap-aware.
nS = numel(sites);
sLen = zeros(nS, 1);
for k = 1:nS
    sLen(k) = mod(sites(k).Stop - sites(k).Start, N) + 1;
end
if isempty(opts.uMax), opts.uMax = min(60, max([sLen; opts.NucLen])); end
if opts.uMax > opts.W
    error('DNAAccessibility:uTooLong', 'uMax (%d) exceeds W (%d).', opts.uMax, opts.W);
end

%% ---- locate the DNA parameter file --------------------------------------------------------------
if isempty(opts.ParamFile)
    opts.ParamFile = FindDNAParams();
end
if isempty(opts.ParamFile) || ~isfile(opts.ParamFile)
    error('DNAAccessibility:NoParamFile', ...
        ['dna_mathews2004.par was not found. RNAplfold would fall back to the RNA energy model, ' ...
        'which permits G-U pairs and is simply wrong for DNA, so this is a hard stop rather than ' ...
        'a warning. Pass opts.ParamFile with the full path. On a Windows install it normally sits ' ...
        'under the ViennaRNA Package folder in Misc\\ or share\\ViennaRNA\\.']);
end

%% ---- build the run sequence ---------------------------------------------------------------------
if opts.Circular
    runSeq = [seq seq seq];    % read the middle copy so every position has full context both sides
    off = N;
else
    runSeq = seq;
    off = 0;
end
posMap = @(p) mod(p - 1, N) + 1 + off;

%% ---- run RNAplfold -------------------------------------------------------------------------------
if isempty(opts.OutDir)
    opts.OutDir = fullfile(tempdir, ['DNAaccess_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
end
if ~exist(opts.OutDir, 'dir'), mkdir(opts.OutDir); end
fastaFile = [opts.Header '.fas'];
oldDir = pwd; cd(opts.OutDir);
cleanup = onCleanup(@() cd(oldDir)); %#ok<NASGU>
fid = fopen(fastaFile, 'w'); fprintf(fid, ['>' opts.Header '\n%s\n'], runSeq); fclose(fid);

cmd = sprintf('RNAplfold -d2 -W %d -L %d -u %d --temp=%f --salt=%f --paramFile="%s" < %s', ...
    opts.W, opts.L, opts.uMax, opts.Temp, opts.Salt, opts.ParamFile, fastaFile);
[status, cmdout] = system(cmd);
cd(oldDir);
if status ~= 0
    error('DNAAccessibility:RNAplfoldFailed', 'RNAplfold failed: %s', cmdout);
end

lunpFile = fullfile(opts.OutDir, [opts.Header '_lunp']);
if ~isfile(lunpFile)
    error('DNAAccessibility:NoOutput', 'RNAplfold produced no _lunp file in %s.', opts.OutDir);
end
lunp = readmatrix(lunpFile, 'FileType', 'text', 'NumHeaderLines', 2);

%% ---- read the sites ------------------------------------------------------------------------------
% Column u+1 of lunp holds P(the stretch of length u ENDING at that row is unpaired), so a footprint
% [a..b] of length u is read at row posMap(b), column u+1.
Label = strings(nS, 1); Start = zeros(nS, 1); Stop = zeros(nS, 1); Len = zeros(nS, 1);
pFull = nan(nS, 1); dGopenFull = nan(nS, 1);
pNuc = nan(nS, 1); dGopenNuc = nan(nS, 1); NucStop = nan(nS, 1);
TK = opts.Temp + 273.15;

for k = 1:nS
    Label(k) = string(sites(k).Label);
    Start(k) = sites(k).Start; Stop(k) = sites(k).Stop; Len(k) = sLen(k);

    if Len(k) <= opts.uMax
        pFull(k) = lunp(posMap(Stop(k)), Len(k) + 1);
    end

    % Best nucleation window anywhere inside the footprint.
    uN = min(opts.NucLen, Len(k));
    bestp = 0; bestEnd = NaN;
    for e = 0:(Len(k) - uN)
        endPos = mod(Stop(k) - 1 - e, N) + 1;
        p = lunp(posMap(endPos), uN + 1);
        if ~isnan(p) && p > bestp, bestp = p; bestEnd = endPos; end
    end
    pNuc(k) = bestp; NucStop(k) = bestEnd;
end
dGopenFull = -R*TK*log(pFull);
dGopenNuc  = -R*TK*log(pNuc);

siteTable = table(Label, Start, Stop, Len, pFull, dGopenFull, pNuc, dGopenNuc, NucStop);
siteTable.Properties.Description = sprintf(...
    'RNAplfold DNA accessibility, W=%d L=%d u=%d, %.1f C, %.3f M monovalent, params %s', ...
    opts.W, opts.L, opts.uMax, opts.Temp, opts.Salt, opts.ParamFile);

if ~opts.KeepFiles
    try, rmdir(opts.OutDir, 's'); catch, end %#ok<CTCH>
end
end

%% ____________________________________________________________________________________________________________
function p = FindDNAParams()
% Look in the usual places for the DNA energy parameters before giving up.
cands = { ...
    'C:\Program Files (x86)\ViennaRNA Package\Misc\dna_mathews2004.par', ...
    'C:\Program Files\ViennaRNA Package\Misc\dna_mathews2004.par', ...
    'C:\Program Files (x86)\ViennaRNA Package\share\ViennaRNA\dna_mathews2004.par', ...
    'C:\Program Files\ViennaRNA Package\share\ViennaRNA\dna_mathews2004.par', ...
    '/usr/share/ViennaRNA/dna_mathews2004.par', ...
    '/usr/local/share/ViennaRNA/dna_mathews2004.par', ...
    '/opt/homebrew/share/ViennaRNA/dna_mathews2004.par'};
p = '';
for k = 1:numel(cands)
    if isfile(cands{k}), p = cands{k}; return; end
end
% Last resort: ask.
[f, d] = uigetfile({'*.par', 'ViennaRNA parameter file'}, ...
    'Locate dna_mathews2004.par (DNA energy parameters)');
if ischar(f), p = fullfile(d, f); end
end
