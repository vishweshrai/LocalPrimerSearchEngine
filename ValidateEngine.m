function report = ValidateEngine(mode, folder)
%% Validation Harness for the aPCR Primer Screening Engine
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 12-11-2026 21:26 IST
%
% Revision history is not kept here. See ValidateEngine_<YYYYMMDD>.txt beside this file.
% ____________________________________________________________________________________________________________
% What This Script Does
%
% Runs the checks that can be made without a bench result and reports each one against a stated
% tolerance. A pass means the engine still reproduces its recorded behaviour and stays within the
% accuracy of the published model it is built on. A pass is not a measurement of the reaction.
%
% Tests in mode 'quick' (about 30 s)
%   T1  Hairpin panel. LooperMismatch against Primer3 on the 59 hairpins of ValidationHairpinsPrimer3.csv,
%       at 58 C in 40 mM monovalent, 2.5 mM divalent, 0.8 mM dNTP. Reports the mean and SD of the dG and
%       Tm differences and the correlation. Tolerances: |mean Tm error| <= 4.0 C and SD <= 3.9 C, the
%       accuracy SantaLucia and Hicks 2004 report for their own model on 61 hairpins and 320 beacons.
%       The dG tolerances (|mean| <= 0.5 kcal/mol, r >= 0.9) are manual values with no published source.
%   T2  Collapse tiers on two fixed strands. A strand that ends in the complement of a 14 bp GC stem
%       upstream must read as the reject tier in the declared channel. A strand without a designed stem
%       must not.
%   T3  Collapse penalty arithmetic. Exact values of the Score input at the tolerance, one and two
%       decades above it, with the flag increment and at the cap.
%   T4  Channel invariants. The Taq-channel extendable fraction never exceeds the exonuclease
%       channel's, every fraction lies in [0 1], and the declared-channel tier never outranks the tier
%       over both channels.
%   T5  Parameter consistency. PrimerSearchParameters runs on a reference input, which also runs its
%       own check that the sheet and the cfg agree, and the collapse weight stays tied to MarginCap.
%
% Test in mode 'full' (adds about 15 min in Octave)
%   T6  Golden regression. The 47 pairs of the 12-11-2026 02:27 MATLAB run are recomputed from the
%       plasmid and compared with the extendable fractions that run logged
%       (ValidationGoldenCollapse.csv). Tolerance: relative difference <= 1E-6 on pExtTot.
%
% Usage
%   report = ValidateEngine()                 % quick, files read from the folder of this script
%   report = ValidateEngine('full')
%   report = ValidateEngine('quick', folder)  % folder holding the csv files and ValidationM13mp18.txt
%
% Output
%   report  struct array with fields name, pass, value and detail, one element per test.
%   The table is also printed. The function errors out only when a data file is missing.
%
% Regenerating the reference files
%   ValidationHairpinsPrimer3.csv   primer3-py 2.4.0, calc_hairpin(seq, mv_conc = 40, dv_conc = 2.5,
%                                   dntp_conc = 0.8, temp_c = 58), structures that Primer3 found
%   ValidationGoldenCollapse.csv    Collapse_Plus_pExtTot and Collapse_Minus_pExtTot of the
%                                   Candidate_Primers sheet of the logged run
%   ValidationM13mp18.txt           the M13mp18 sequence the logged run used
% ____________________________________________________________________________________________________________
% References
%   SantaLucia J Jr, Hicks D. The thermodynamics of DNA structural motifs. Annu Rev Biophys Biomol
%       Struct 2004;33:415-40. doi:10.1146/annurev.biophys.32.110601.141800
%   Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M, Rozen SG. Primer3 - new
%       capabilities and interfaces. Nucleic Acids Res 2012;40:e115. doi:10.1093/nar/gks596
% ____________________________________________________________________________________________________________

if nargin < 1||isempty(mode), mode = 'quick'; end
if nargin < 2||isempty(folder), folder = fileparts(mfilename('fullpath')); end
mode = lower(char(mode));
if ~any(strcmp(mode, {'quick', 'full'}))
    error('ValidateEngine:BadMode', 'Mode must be quick or full.');
end
addpath(folder);

%% Tolerances
TolTmMean = 4.0;    % C, mean Tm error of the hairpin model, SantaLucia and Hicks 2004
TolTmSd = 3.9;      % C, SD of the Tm error on the beacon set, SantaLucia and Hicks 2004
TolDGMean = 0.5;    % kcal/mol, manual value
TolCorr = 0.9;      % correlation with Primer3, manual value
TolGolden = 1E-6;   % relative, MATLAB against the recorded run

report = struct('name', {}, 'pass', {}, 'value', {}, 'detail', {});

%% T1 Hairpin Panel
[seqs, p3dG, p3Tm] = ReadPanel(fullfile(folder, 'ValidationHairpinsPrimer3.csv'));
N = NNParams();
n = numel(seqs);
dG = zeros(n, 1); Tm = zeros(n, 1);
for k = 1:n
    [~, ~, dGhp, hpTm] = LooperMismatch(seqs{k}, 0.04, 0.0025, 0.0008, 58, -3.0, 58, 5, 3, N);
    dG(k) = dGhp; Tm(k) = hpTm;
end
ok = isfinite(Tm);
ddG = dG(ok) - p3dG(ok); dTm = Tm(ok) - p3Tm(ok);
rG = corrcoef(dG(ok), p3dG(ok)); rG = rG(1, 2);
rT = corrcoef(Tm(ok), p3Tm(ok)); rT = rT(1, 2);
pass = ((abs(mean(dTm)) <= TolTmMean)&&(std(dTm) <= TolTmSd)&&(abs(mean(ddG)) <= TolDGMean)&&(rG >= TolCorr)&&(rT >= TolCorr));
detail = sprintf('%d hairpins, dG error mean %+.2f SD %.2f kcal/mol (r %.3f), Tm error mean %+.2f SD %.2f C (r %.3f)', ...
    nnz(ok), mean(ddG), std(ddG), rG, mean(dTm), std(dTm), rT);
report(end + 1) = struct('name', 'T1 hairpin panel against Primer3', 'pass', pass, 'value', mean(dTm), 'detail', detail);

%% T2 Collapse Tiers
body = 'CCGTAATGCCTTTCCCTAACAGAGTTTTTCGAACTCGTGTTGTCGAGCGACGGAATTAGA';
ctlTail = 'TCAGTTAAATGGCAGAAAACTGGCAGGGCTTTTAGTCGTGGGATGATCAGTGGGTAAAGGTGGCGCGGGGTAAC';
stem = 'GGCCGCGCCGGCGC';
designed = [body stem repmat('AT', 1, 15) char(seqrcomplement(stem))];
control = [body ctlTail];
opts = struct('StrandConc', 1E-6, 'LBase', 0.6, 'LPers', 0.9, 'HpTmMargin', 5, 'FlagEps', 1E-3, ...
    'Bimolecular', true, 'ExoBlend', NaN, 'ExoEndTol', 15);
evD = zeros(1, 2); ctl = struct();
for q = 1:2
    sr = ScaffoldCollapseCheck('scan', designed, 18, 3, 1, 0.04, 0.0025, 0.0008, 58, (q == 1), opts);
    ev = ScaffoldCollapseCheck('eval', sr, 1, 0);
    evD(q) = ev.tierDeclRank;
end
sr = ScaffoldCollapseCheck('scan', control, 18, 3, 1, 0.04, 0.0025, 0.0008, 58, true, opts);
ctl = ScaffoldCollapseCheck('eval', sr, 1, 0);
pass = ((evD(1) == 3)&&(evD(2) == 3)&&(ctl.tierDeclRank < 3));
detail = sprintf('designed stem: tier rank %d in the Taq channel and %d in the exonuclease channel (3 is reject); control: rank %d, tier %s', ...
    evD(1), evD(2), ctl.tierDeclRank, char(ctl.tierDecl));
report(end + 1) = struct('name', 'T2 collapse tiers on fixed strands', 'pass', pass, 'value', evD(1), 'detail', detail);

%% T3 Penalty Arithmetic
tol = 1E-3;
got = [ScaffoldCollapseCheck('penalty', tol, 0, tol, 1.0, 3.0), ...
    ScaffoldCollapseCheck('penalty', 10*tol, 0, tol, 1.0, 3.0), ...
    ScaffoldCollapseCheck('penalty', 100*tol, 0, tol, 1.0, 3.0), ...
    ScaffoldCollapseCheck('penalty', 10*tol, 1, tol, 1.0, 3.0), ...
    ScaffoldCollapseCheck('penalty', 0.5, 1, tol, 1.0, 3.0), ...
    ScaffoldCollapseCheck('penalty', 0, 1, tol, 1.0, 3.0), ...
    ScaffoldCollapseCheck('penalty', NaN, 0, tol, 1.0, 3.0)];
want = [0 1 2 2 3 0 0];
pass = all(abs(got - want) < 1E-9);
detail = sprintf('penalty at (1, 10, 100) x tolerance, flag, cap, zero, NaN = %s', mat2str(got, 4));
report(end + 1) = struct('name', 'T3 collapse penalty arithmetic', 'pass', pass, 'value', max(abs(got - want)), 'detail', detail);

%% T4 Channel Invariants
okAll = true; worst = 0;
for q = 1:2
    if q == 1, sr = ScaffoldCollapseCheck('scan', designed, 18, 3, 1, 0.04, 0.0025, 0.0008, 58, true, opts); ev = ScaffoldCollapseCheck('eval', sr, 1, 0);
    else, sr = ScaffoldCollapseCheck('scan', control, 18, 3, 1, 0.04, 0.0025, 0.0008, 58, true, opts); ev = ScaffoldCollapseCheck('eval', sr, 1, 0); end
    fr = [ev.pExtTaq ev.pExtExo ev.pExtTot];
    okAll = okAll&&all((fr >= 0)&(fr <= 1 + 1E-12))&&(ev.pExtTaq <= ev.pExtExo + 1E-12)&&(ev.tierDeclRank <= ev.tierAllRank);
    worst = max(worst, ev.pExtTaq - ev.pExtExo);
end
report(end + 1) = struct('name', 'T4 channel invariants', 'pass', okAll, 'value', worst, ...
    'detail', 'Taq fraction <= exonuclease fraction, all fractions in [0 1], declared tier <= tier over both channels');

%% T5 Parameter Consistency
plasFile = fullfile(folder, 'ValidationM13mp18.txt');
if exist(plasFile, 'file')
    plas = upper(regexprep(fileread(plasFile), '[^ACGTacgt]', ''));
    input = struct('MonovalentCationConcentration', 0.04, 'BivalentCationConcentration', 0.0025, 'dNTPsConcentration', 0.0008, ...
        'ScaffoldStart', 2173, 'ScaffoldEnd', 6112, 'ScaffoldLength', 3940, 'ForwardPrimerConcentration', 1E-6, ...
        'ReversePrimerConcentration', 2E-8, 'AnnealingTemperature', 58, 'ExtensionTemperature', 65, 'ExtensionTime', 60, ...
        'FivePrimeExo', false, 'ForwardPrimer', '', 'ReversePrimer', '', 'MaxLaps', 1, 'NumAmplificationCycles', 40, ...
        'PlasmidConcentration', 1E-10, 'PlasmidFileName', 'M13mp18.txt', 'PlasmidSequence', plas, 'PlasmidStranded', 'ssDNA', ...
        'PolymeraseActivity', 1, 'PolymeraseEfficiency', 0.9, 'PolymeraseName', 'LongAmp Taq', 'ProofreadingClass', 'none', ...
        'ProofreadingExo', false, 'UnitDNTP', 'mM', 'UnitTime', 's');
    cfg = PrimerSearchParameters(input);
    tie = abs(cfg.WtCollapse*cfg.CollapseCap - cfg.MarginCap);
    pass = ((tie < 1E-9)&&(cfg.CollapseTol > 0)&&(cfg.CollapseCap > 0)&&(cfg.CollapseFlagDecades >= 0)&&islogical(cfg.CollapseGate));
    detail = sprintf('WtCollapse*CollapseCap - MarginCap = %.1e, CollapseTol %.1e, CollapseGate %d', tie, cfg.CollapseTol, cfg.CollapseGate);
    report(end + 1) = struct('name', 'T5 parameter consistency', 'pass', pass, 'value', tie, 'detail', detail);
else
    report(end + 1) = struct('name', 'T5 parameter consistency', 'pass', false, 'value', NaN, 'detail', 'ValidationM13mp18.txt not found beside the script');
end

%% T6 Golden Regression
if strcmp(mode, 'full')
    plasFile = fullfile(folder, 'ValidationM13mp18.txt');
    plas = upper(regexprep(fileread(plasFile), '[^ACGTacgt]', ''));
    G = ReadGolden(fullfile(folder, 'ValidationGoldenCollapse.csv'));
    PlasLen = numel(plas); PlasDbl = [plas plas];
    cacheP = containers.Map('KeyType', 'double', 'ValueType', 'any');
    cacheM = containers.Map('KeyType', 'double', 'ValueType', 'any');
    optsP = opts; optsM = opts; optsM.StrandConc = 2E-8;
    relP = zeros(size(G, 1), 1); relM = relP;
    for k = 1:size(G, 1)
        [pP, pM] = PairCollapse(G(k, 1), G(k, 2), PlasDbl, PlasLen, cacheP, cacheM, optsP, optsM);
        relP(k) = abs(pP - G(k, 3))/max(G(k, 3), realmin);
        relM(k) = abs(pM - G(k, 4))/max(G(k, 4), realmin);
    end
    worst = max([relP; relM]);
    report(end + 1) = struct('name', 'T6 golden regression', 'pass', (worst <= TolGolden), 'value', worst, ...
        'detail', sprintf('worst relative difference in pExtTot %.2e over %d pairs and two strands', worst, size(G, 1)));
end

%% Print
fprintf('\nValidateEngine, mode %s\n', mode);
fprintf('%s\n', repmat('-', 1, 110));
for k = 1:numel(report)
    if report(k).pass, tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('%-4s  %-38s  %s\n', tag, report(k).name, report(k).detail);
end
fprintf('%s\n%d of %d passed\n', repmat('-', 1, 110), nnz([report.pass]), numel(report));
end

%% Read Panel
function [seqs, dG, Tm] = ReadPanel(path)
fid = fopen(path, 'r');
if fid < 0, error('ValidateEngine:NoFile', 'Cannot open %s', path); end
fgetl(fid); seqs = {}; dG = []; Tm = [];
while true
    l = fgetl(fid);
    if ~ischar(l), break; end
    c = strsplit(strtrim(l), ',');
    seqs{end + 1} = c{1}; dG(end + 1, 1) = str2double(c{2}); Tm(end + 1, 1) = str2double(c{3});
end
fclose(fid);
end

%% Read Golden
function G = ReadGolden(path)
fid = fopen(path, 'r');
if fid < 0, error('ValidateEngine:NoFile', 'Cannot open %s', path); end
fgetl(fid); G = [];
while true
    l = fgetl(fid);
    if ~ischar(l), break; end
    c = strsplit(strtrim(l), ',');
    G(end + 1, :) = [str2double(c{1}) str2double(c{2}) str2double(c{3}) str2double(c{4})];
end
fclose(fid);
end

%% Pair Collapse
% Mirrors CollapseForPair in LocalPrimerSearchEngine for the reference run: 20 nt flexibility on both
% ends, a window of 18 nt, stems of 3 nt or more and the annealing temperature of 58 C.
function [pPlus, pMinus] = PairCollapse(fPos, rPos, PlasDbl, PlasLen, cacheP, cacheM, optsP, optsM)
monv = 0.04; biv = 0.0025; dNTP = 0.0008; tann = 58;
scaffStart = 2173; scaffEnd = 6112; fFlex = 20; rFlex = 20;
fEarliest = mod(scaffStart - fFlex - 1, PlasLen) + 1;
rLatest = mod(scaffEnd + rFlex - 1, PlasLen) + 1;
if ~isKey(cacheP, rPos)
    len = mod(rPos - fEarliest, PlasLen) + 1;
    cacheP(rPos) = ScaffoldCollapseCheck('scan', PlasDbl(fEarliest:fEarliest + len - 1), 18, 3, 1, monv, biv, dNTP, tann, true, optsP);
end
if ~isKey(cacheM, fPos)
    len = mod(rLatest - fPos, PlasLen) + 1;
    cacheM(fPos) = ScaffoldCollapseCheck('scan', char(seqrcomplement(PlasDbl(fPos:fPos + len - 1))), 18, 3, 1, monv, biv, dNTP, tann, true, optsM);
end
evP = ScaffoldCollapseCheck('eval', cacheP(rPos), mod(fPos - fEarliest, PlasLen) + 1, 0);
evM = ScaffoldCollapseCheck('eval', cacheM(fPos), mod(rLatest - rPos, PlasLen) + 1, 0);
pPlus = evP.pExtTot; pMinus = evM.pExtTot;
end
