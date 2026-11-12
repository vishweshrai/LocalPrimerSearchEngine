function T = CalibrateAggregateFlux(results, opts)
%% CALIBRATEAGGREGATEFLUX  How far the off-target ensemble sits from the single worst site
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 05-11-2026 22:46 IST
% ____________________________________________________________________________________________________________
% WHAT THIS ANSWERS
%
% The mispriming gate asks about one register: the most stable site the primer can prime from. The
% quantity that actually competes with the cognate site is the sum over every such site, because a
% candidate carrying many comparable off-target registers initiates more spurious product per cycle
% than one carrying a single slightly worse register. This script measures the difference on real
% candidates before anything is allowed to gate on it.
%
% THE TWO QUANTITIES BEING COMPARED
%
% The engine already records, per candidate and per strand:
%
%   dG_eff_min   the most stable anchored register, extension-weighted        (…AnchEff)
%   dG_ens       -RT*log( sum_k exp(-dG_eff_k / RT) ) over ALL anchored       (…AnchEns)
%   dG_ens_run   the same sum over registers carrying MispMinRun identities   (…AnchEnsRun)
%
% The shift dG_ens - dG_eff_min is never positive. It is zero when one site dominates and equals
% -RT*log(M) when M sites are degenerate, which at 58 C is -0.46 kcal/mol for two sites,
% -1.52 for ten and -3.03 for a hundred. Because the weights are exponential, a site 3 kcal/mol
% weaker than the best contributes about a hundredth as much, so the shift tracks the multiplicity
% of COMPARABLE sites rather than the raw count of anchored registers.
%
% WHY THE ENSEMBLE IS EXACT RATHER THAN APPROXIMATE
%
% At phi = 1E-4 every off-target site sits at an occupancy of order 1E-4, far inside the linear
% regime where theta_i = K_i*[p]. There the ensemble association constant is literally the sum of
% the individual ones, so adding Boltzmann factors introduces no approximation.
%
% WHAT THE OUTPUT DECIDES
%
% Three things, none of which should be settled by argument:
%
%   1. Whether a gate on the ensemble would reject candidates the present gate accepts, and how
%      many. If the shift is a few hundredths of a kcal/mol the change is cosmetic.
%   2. Whether to sum over all anchored registers or only over run-qualified ones. The two columns
%      diverge exactly when sub-minRun registers carry appreciable weight.
%   3. Whether adjacent overlapping registers inflate the sum. Overlapping duplexes on one stretch
%      of template are mutually exclusive and summing them over-counts. A shift that tracks
%      RT*log(nAnchRun) rather than saturating at a small multiplicity is the signature of that.
%
% ____________________________________________________________________________________________________________
% INPUTS   results  the second output of LocalPrimerSearchEngine, or a .mat file holding it
%          opts     optional struct
%                   .Tann      annealing temperature in C. OMIT IT: the temperature the run
%                              actually annealed at is read from results(1).Parameters, and a
%                              value passed here overrides that. Used for the RT*log(M) column
%                   .Selected  true to restrict to Status == "Selected", default true
% OUTPUT   T        one row per candidate per strand, sorted by the size of the shift
% ____________________________________________________________________________________________________________

if (nargin < 2)||isempty(opts), opts = struct(); end
if ~isfield(opts, 'Selected')||isempty(opts.Selected), opts.Selected = true; end
if ischar(results)||isstring(results)
    S = load(char(results));
    f = fieldnames(S);
    results = S.(f{1});
end

% THE TEMPERATURE IS READ FROM THE RUN, NOT ASSUMED HERE.
% RT scales every column in this table, so assuming a temperature means reporting a shift in
% kcal/mol that the run never saw. This used to default to 58 while the reaction dialog pre-fills
% 55, which is the kind of disagreement that produces a table nobody can reconcile with the run it
% came from. The resolution happens AFTER the .mat load, because a file handed in by name carries
% its own Parameters sheet and that is the authority.
if ~isfield(opts, 'Tann')||isempty(opts.Tann)
    [opts.Tann, tannSrc] = TannFromResults(results);
else
    tannSrc = 'caller-supplied opts.Tann';
end

RT = 1.987E-3*(273.15 + opts.Tann);

rows = {};
rows = Harvest(rows, results(1).Forward_Log, 'forward', 'fseq', ...
    {'dGfTcAnchEff','dGfTcAnchEns','dGfTcAnchEnsRun','fAnchC','fAnchRunC'}, 'cognate',  opts);
rows = Harvest(rows, results(1).Forward_Log, 'forward', 'fseq', ...
    {'dGfTAnchEff','dGfTAnchEns','dGfTAnchEnsRun','fAnchT','fAnchRunT'},   'opposite', opts);
rows = Harvest(rows, results(1).Reverse_Log, 'reverse', 'rseq', ...
    {'dGrTcAnchEff','dGrTcAnchEns','dGrTcAnchEnsRun','rAnchC','rAnchRunC'}, 'cognate',  opts);
rows = Harvest(rows, results(1).Reverse_Log, 'reverse', 'rseq', ...
    {'dGrTAnchEff','dGrTAnchEns','dGrTAnchEnsRun','rAnchT','rAnchRunT'},   'opposite', opts);

if isempty(rows)
    warning('CalibrateAggregateFlux:NoRows', ...
        ['No candidate carried the flux columns. They are written only by a build that includes ' ...
        'the Phase 2a instrumentation, and only for candidates that reached the mispriming ' ...
        'stage. Check that the log came from this build.']);
    T = table(); return
end

T = vertcat(rows{:});
T.ShiftAll = T.dG_ens     - T.dG_eff_min;
T.ShiftRun = T.dG_ens_run - T.dG_eff_min;
T.RTlogNrun = -RT*log(max(T.nAnchRun, 1));
T = sortrows(T, 'ShiftAll', 'ascend');

fprintf('\n');
fprintf('AGGREGATE OFF-TARGET FLUX, %d candidate-strand rows at %.1f C (RT = %.4f kcal/mol)\n', ...
    height(T), opts.Tann, RT);
fprintf('  annealing temperature taken from: %s\n', tannSrc);
fprintf('%s\n', repmat('-', 1, 96));
Report('shift, all anchored registers', T.ShiftAll);
Report('shift, run-qualified only',     T.ShiftRun);
fprintf('\n');
fprintf('  anchored registers per row : median %.0f, max %.0f\n', ...
    median(T.nAnch), max(T.nAnch));
fprintf('  run-qualified per row      : median %.0f, max %.0f\n', ...
    median(T.nAnchRun), max(T.nAnchRun));
fprintf('\n');
fprintf('  If ShiftAll clusters near -RT*log(2) = %.3f the sum is dominated by a couple of\n', -RT*log(2));
fprintf('  comparable sites. If it tracks the RTlogNrun column instead, overlapping registers are\n');
fprintf('  being counted as independent and the population needs narrowing before anything gates.\n');
fprintf('%s\n', repmat('-', 1, 96));
fprintf('\nTen largest shifts:\n');
disp(head(T(:, {'primer','strand','side','dG_eff_min','dG_ens','dG_ens_run','nAnch','nAnchRun','ShiftAll'}), 10));
end

%% ____________________________________________________________________________________________________________
function rows = Harvest(rows, L, strand, seqcol, cols, side, opts)
%% HARVEST  Pull one strand's flux columns out of one log, skipping rows that never reached it
if isempty(L)||(height(L) == 0), return; end
have = all(ismember([cols, {seqcol}], L.Properties.VariableNames));
if ~have, return; end
keep = true(height(L), 1);
if opts.Selected && ismember('Status', L.Properties.VariableNames)
    keep = strcmp(string(L.Status), "Selected");
end
% A candidate rejected before the mispriming stage carries NaN in these columns by construction.
keep = keep & isfinite(L.(cols{2})) & isfinite(L.(cols{1}));
if ~any(keep), return; end
rows{end+1} = table(string(L.(seqcol)(keep)), repmat(string(strand), nnz(keep), 1), ...
    repmat(string(side), nnz(keep), 1), ...
    L.(cols{1})(keep), L.(cols{2})(keep), L.(cols{3})(keep), ...
    L.(cols{4})(keep), L.(cols{5})(keep), ...
    'VariableNames', {'primer','strand','side','dG_eff_min','dG_ens','dG_ens_run','nAnch','nAnchRun'});
end

%% ____________________________________________________________________________________________________________
function Report(label, v)
%% REPORT  One distribution line, in the units the gate is written in
v = v(isfinite(v));
if isempty(v), fprintf('  %-32s no finite values\n', label); return; end
fprintf('  %-32s median %+7.4f   mean %+7.4f   worst %+7.4f   |shift| > 0.5 : %d of %d\n', ...
    label, median(v), mean(v), min(v), nnz(abs(v) > 0.5), numel(v));
end

%% ____________________________________________________________________________________________________________
function [tann, src] = TannFromResults(results)
%% TANNFROMRESULTS  The annealing temperature the run actually used, from its own Parameters sheet
% LocalPrimerSearchEngine writes a resolved numeric TAnneal row. Older result files predate that
% row and carry only input.AnnealingTemperature, which is rendered as text because it may hold a
% mode flag or a touchdown recipe; a plain number there is still usable, anything else is not.
% When neither yields a number this returns 58 and says so, because a silent assumption is the
% defect this function exists to remove.
tann = 58; src = 'FALLBACK 58 C, no usable temperature found on the Parameters sheet';
if ~isstruct(results)||~isfield(results, 'Parameters')||isempty(results(1).Parameters)
    warning('CalibrateAggregateFlux:NoParameterSheet', ...
        ['The results carry no Parameters sheet, so the annealing temperature could not be read ' ...
        'and 58 C was assumed. Every kcal/mol in this table scales with RT, so check that 58 is ' ...
        'what the run used before quoting any of it.']);
    return
end
P = results(1).Parameters;
if ~ismember('Parameter', P.Properties.VariableNames)||~ismember('Value', P.Properties.VariableNames)
    return
end
names = string(P.Parameter);
for want = ["TAnneal", "input.AnnealingTemperature"]
    k = find(names == want, 1);
    if isempty(k), continue; end
    v = str2double(string(P.Value(k)));
    if isfinite(v)&&(v > 0)&&(v < 120)
        tann = v;
        src = sprintf('%s = %.4g C on the run''s Parameters sheet', want, v);
        return
    end
end
warning('CalibrateAggregateFlux:TemperatureNotOnSheet', ...
    ['Neither TAnneal nor a numeric input.AnnealingTemperature was found on the Parameters ' ...
    'sheet, so 58 C was assumed. A touchdown recipe or a mode flag renders as text and cannot ' ...
    'be read as a number; re-run the engine to get the resolved TAnneal row, or pass opts.Tann.']);
end
