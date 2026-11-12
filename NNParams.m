function N = NNParams(opts)
%% NNPARAMS  Nearest-neighbour tables and salt corrections for the aPCR suite
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 05-11-2026 18:49 IST
%
% Revision history is NOT kept here. See NNParams_<YYYYMMDD>.txt beside this file.
% ____________________________________________________________________________________________________________
% WHAT THIS MODULE IS FOR
%
% Four files in this suite each carry their own copy of the Watson-Crick nearest-neighbour table and
% their own copy of the Owczarzy salt block: PrimTemp.m, dGPrim.m, dGMismatch.m and, transitively,
% everything that calls them. The tables are identical by intention and divergent by accident, which
% is how dGPrim.m came to apply a correction that dGMismatch.m applies differently to the same
% duplex. One module, read once, removes that class of defect outright.
%
% It also carries the two things Stage 1 adds: mismatch enthalpies and entropies, so a duplex can be
% scored at the temperature the reaction runs at instead of at 37 C, and a per-phosphate salt
% correction for duplexes that carry mismatches.
%
% ____________________________________________________________________________________________________________
% THE MISMATCH PARAMETERS ARE LOADED, NOT HARD-CODED, AND THIS IS DELIBERATE
%
% dGMismatch.m currently carries 48 Watson-Crick-free entries of Table 2 as DeltaG37 ONLY, with a
% header note conceding that they are added to the Watson-Crick term with no temperature scaling.
% At an annealing temperature twenty degrees above the reference that is the single largest
% systematic error in the specificity calculation, because discrimination is a difference of
% differences and a mismatch-dependent offset does not cancel out of it.
%
% The fix needs DeltaH and DeltaS for all 48. Those exist, in four papers. They are NOT reproduced
% here from memory, because a wrong entry in a hard-coded table of 96 numbers is undetectable by any
% assertion this suite could carry and would silently corrupt the one quantity the whole exercise
% turns on. mismatch_nn_params.csv ships beside this file carrying, for each of the 48 entries, the
% DeltaG37 value from the block you verified against the paper yourself, the mismatch identity, the
% nearest-neighbour step in the notation the papers use, and the paper and DOI it comes from. Two
% columns to fill. Rows are grouped so each source paper is contiguous.
%
% Until they are filled, N.MM.HasHS is false for those entries and the caller is told which route it
% is on. See the three routes below.
%
% ____________________________________________________________________________________________________________
% THREE ROUTES TO MISMATCH TEMPERATURE DEPENDENCE, IN DESCENDING ORDER OF DEFENSIBILITY
%
%  (a) MEASURED. DeltaH and DeltaS from the four papers, filled into the CSV. The mismatch term then
%      carries its own temperature dependence through the same h and s accumulators the
%      Watson-Crick steps use, and nothing is extrapolated. This is what you should end up on.
%
%  (b) VIENNA SLOPE. Keep the verified SantaLucia DeltaG37 magnitudes and take only the temperature
%      SLOPE from the int11 table of dna_mathews2004.par, which ships DeltaG37 and DeltaH together:
%
%          dG_mm(T) = dG_mm_SantaLucia(37) + [ int11_Vienna(T) - int11_Vienna(37) ]
%
%      An isolated single mismatch flanked by Watson-Crick pairs IS a 1x1 internal loop, so int11 is
%      the direct counterpart of the two Table 2 steps SantaLucia uses for it, and the bracketed
%      difference is a pure temperature correction with the 37 C magnitudes divided out. Mixing two
%      parameterisations is a real cost; taking only the slope is the cheapest form of it. Note that
%      this route applies only to an ISOLATED mismatch, since int11 describes nothing else, and
%      tandem mismatches fall back to route (c).
%
%  (c) NONE. DeltaG37 with no temperature scaling, which is the status quo and the defect.
%
% The resolution order is (a) where the CSV supplies both columns for that entry, then (b) where
% ViennaParams is available and the mismatch is isolated, then (c) with a warning naming the count.
% Route (c) is never taken silently. N.MM.Route records what actually happened per entry.
%
% ____________________________________________________________________________________________________________
% SALT: TWO CORRECTIONS, AND WHICH DUPLEX EACH ONE BELONGS TO
%
% N.Salt.Owczarzy is the existing treatment, the one every file in the suite already applies. It
% shifts DeltaS by DeltaH times a function of salt, which is the algebraic image of a correction
% defined on 1/Tm. It is valid for a FULLY PAIRED duplex, where a two-state Tm exists and DeltaH
% tracks length closely enough to stand in for phosphate count.
%
% N.Salt.PerPhosphate is new and belongs to any duplex carrying a mismatch:
%
%       DeltaS(Na+) = DeltaS(1 M) + 0.368*(N-1)*ln[Na+]      cal/(K mol)
%
% The reason the distinction matters is not cosmetic. A trimmed misprime register of twelve
% positions with five mismatches has a small DeltaH and an unchanged phosphate count, so the
% DeltaH-scaled form gives it almost none of the low-salt destabilisation it should receive. The
% bias runs one way: mismatched registers come out looking more stable than they are, and the engine
% over-predicts mispriming. Every off-target free energy in the suite currently carries it.
%
% The per-phosphate form is monovalent-only, so N.Salt.MonovalentEquivalent supplies the von Ahsen
% rule, [Na+] + 120*sqrt([Mg2+]free), with free magnesium from the same speciation the rest of the
% suite uses. That rule was derived for duplex melting rather than for this, so run anything that
% turns on it at both the bare monovalent and the equivalent value and keep what survives both.
%
% ____________________________________________________________________________________________________________
% INPUT  opts, struct, all optional
%   .MismatchFile  path to mismatch_nn_params.csv, default beside this file
%   .ViennaFile    path to dna_mathews2004.par, default the location ViennaParams resolves
%   .UseVienna     allow route (b), default true
%   .Verbose       print the load summary, default true
%   .Force         bypass the persistent cache, default false
%
% OUTPUT  N, struct
%   .WC.Step, .WC.dH, .WC.dS            16 Watson-Crick steps, kcal/mol and cal/(K mol)
%   .WC.InitGC, .WC.InitAT, .WC.Sym     [dH dS] for the three correction rows
%   .MM.dG37, .MM.dH, .MM.dS            4x4x4 over closing base, X, Y. NaN at Watson-Crick slots.
%   .MM.HasHS                           4x4x4 logical, true where measured dH and dS were supplied
%   .MM.Route                           4x4x4 uint8: 1 measured, 2 Vienna slope available, 3 none
%   .MM.Missing                         table of entries still lacking dH and dS
%   .MM.ViennaSlope                     @(W, X, Y, V, TdegC) kcal/mol, route (b); NaN when absent
%   .Salt.FreeMg, .Salt.Owczarzy, .Salt.PerPhosphate, .Salt.MonovalentEquivalent
%   .Index.Base                         struct A C G T to 1 2 3 4
%   .Meta
%
% REFERENCES, none of which I have been able to open; verify before citing
%   SantaLucia J, Hicks D. Annu Rev Biophys Biomol Struct 2004;33:415.
%                                               doi:10.1146/annurev.biophys.32.110601.141800
%   SantaLucia J. PNAS 1998;95:1460.            doi:10.1073/pnas.95.4.1460   (per-phosphate salt)
%   Owczarzy R et al. Biochemistry 2004;43:3537. doi:10.1021/bi034621r
%   Owczarzy R et al. Biochemistry 2008;47:5336. doi:10.1021/bi702363u
%   von Ahsen N et al. Clin Chem 2001;47:1956.  doi:10.1093/clinchem/47.11.1956
% ____________________________________________________________________________________________________________

persistent CACHE

if (nargin < 1)||isempty(opts), opts = struct(); end
here = fileparts(mfilename('fullpath'));
def = struct('MismatchFile', fullfile(here, 'mismatch_nn_params.csv'), ...
    'ViennaFile', '', 'UseVienna', true, 'Verbose', true, 'Force', false);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts, fn{k})||isempty(opts.(fn{k})), opts.(fn{k}) = def.(fn{k}); end
end

% EVERY TABLE THIS FUNCTION READS IS RESOLVED HERE, BEFORE THE KEY IS BUILT.
% Each path is wanted twice, once to stamp the cache key and once to read the file further down,
% so all four are constructed in one place and referenced thereafter.
danglingFile = fullfile(here, 'dangling_nn_params.csv');
termMismFile = fullfile(here, 'terminal_mismatch_nn_params.csv');
loopFile     = fullfile(here, 'loop_length_params.csv');
specialFile  = fullfile(here, 'special_loop_params.csv');

% THE KEY STAMPS ALL FIVE TABLES, NOT THE MISMATCH ONE ALONE.
% A key covering a single file lets an edit to any of the other four return the structure that was
% assembled before the edit, for as long as the session lives, with nothing printed to say so. A
% parameter sweep that walks the dangling-end or loop table then reports the same numbers at every
% step. Size and modification time are both carried because a sweep that rewrites a table in place
% usually leaves its byte count untouched.
key = sprintf('UseVienna=%d', opts.UseVienna);
stampFiles = {opts.MismatchFile, danglingFile, termMismFile, loopFile, specialFile};
for k = 1:numel(stampFiles)
    if isfile(stampFiles{k})
        dd = dir(stampFiles{k});
        key = sprintf('%s|%s|%d|%.10f', key, stampFiles{k}, dd.bytes, dd.datenum);
    else
        key = sprintf('%s|%s|absent', key, stampFiles{k});
    end
end
if ~opts.Force && ~isempty(CACHE) && strcmp(CACHE.Meta.Key, key) ...
        && strcmp(ViennaStamp(CACHE.Meta.ViennaPath), CACHE.Meta.ViennaStamp)
    N = CACHE; return
end

bases = 'ACGT';
N = struct();
N.Index.Base = struct('A', 1, 'C', 2, 'G', 3, 'T', 4);

%% ---- Watson-Crick table, identical to the block in dGPrim.m and PrimTemp.m --------------------------
N.WC.Step = {'AA'; 'TT'; 'AT'; 'TA'; 'CA'; 'TG'; 'GT'; 'AC'; 'CT'; 'AG'; 'GA'; 'TC'; 'CG'; 'GC'; 'GG'; 'CC'};
N.WC.dH = [-7.90; -7.90; -7.20; -7.20; -8.50; -8.50; -8.40; -8.40; ...
           -7.80; -7.80; -8.20; -8.20; -10.60; -9.80; -8.00; -8.00];
N.WC.dS = [-22.20; -22.20; -20.40; -21.30; -22.70; -22.70; -22.40; -22.40; ...
           -21.00; -21.00; -22.20; -22.20; -27.20; -24.40; -19.90; -19.90];
N.WC.InitGC = [0.10, -2.80];      % [dH kcal/mol, dS cal/(K mol)]
N.WC.InitAT = [2.30,  4.10];
N.WC.Sym    = [0.00, -1.40];

%% ---- mismatch table --------------------------------------------------------------------------------
% 4x4x4 over closing base, X (top strand), Y (bottom strand). The Watson-Crick (X,Y) slots stay NaN
% and are never reached, because the two lookup branches in dGMismatch only fire on genuine
% mismatches.
N.MM.dG37 = nan(4, 4, 4);
N.MM.dH   = nan(4, 4, 4);
N.MM.dS   = nan(4, 4, 4);
N.MM.HasHS = false(4, 4, 4);

if ~isfile(opts.MismatchFile)
    error('NNParams:NoMismatchFile', ...
        ['mismatch_nn_params.csv was not found at %s. It ships beside NNParams.m and carries the ' ...
        '48 verified DeltaG37 entries; without it there is no mismatch table at all.'], opts.MismatchFile);
end

% READ THE FILE DIRECTLY RATHER THAN THROUGH readtable.
% readtable auto-detects the delimiter, renames headers that it considers invalid identifiers, and
% behaves differently across MATLAB releases and locales. On the first run against this very file it
% renamed the header row and the lookup for 'Closing' then failed, with an error that named what was
% missing but not what had been found. For a 49-line table none of that machinery buys anything.
% ReadSimpleCsv below strips a byte-order mark, picks the delimiter by consistency between the
% header and the first data row, and reports the header row verbatim when a column cannot be
% matched. Column names are matched case-insensitively after punctuation is stripped, so a
% round trip through Excel that turns dG37_kcal_per_mol into dG37 kcal per mol still resolves.
[hdr, cells] = ReadSimpleCsv(opts.MismatchFile);
col = MatchColumns(hdr, ...
    {'Closing', 'X', 'Y', 'dG37_kcal_per_mol', 'dH_kcal_per_mol', 'dS_cal_per_K_per_mol'}, ...
    opts.MismatchFile);

for r = 1:size(cells, 1)
    ic = strfind(bases, upper(strtrim(cells{r, col.Closing})));
    ix = strfind(bases, upper(strtrim(cells{r, col.X})));
    iy = strfind(bases, upper(strtrim(cells{r, col.Y})));
    if isempty(ic)||isempty(ix)||isempty(iy)
        error('NNParams:BadBase', ...
            ['Row %d of %s carries a base outside {A,C,G,T}: Closing ''%s'', X ''%s'', Y ''%s''.'], ...
            r + 1, opts.MismatchFile, cells{r, col.Closing}, cells{r, col.X}, cells{r, col.Y});
    end
    N.MM.dG37(ic, ix, iy) = ToNum(cells{r, col.dG37_kcal_per_mol});
    h = ToNum(cells{r, col.dH_kcal_per_mol});
    s = ToNum(cells{r, col.dS_cal_per_K_per_mol});
    N.MM.dH(ic, ix, iy) = h;
    N.MM.dS(ic, ix, iy) = s;
    N.MM.HasHS(ic, ix, iy) = isfinite(h) && isfinite(s);
end

nEntries = nnz(isfinite(N.MM.dG37));
if nEntries ~= 48
    error('NNParams:WrongEntryCount', ...
        ['The mismatch file supplied %d DeltaG37 entries where 48 were expected. Table 2 holds ' ...
        '48 non-Watson-Crick combinations across the four closing bases.'], nEntries);
end

% Consistency of a filled row: dG37 must equal dH - 310.15*dS/1000 to within rounding.
badrow = {};
for ic = 1:4
    for ix = 1:4
        for iy = 1:4
            if ~N.MM.HasHS(ic, ix, iy), continue; end
            recon = N.MM.dH(ic, ix, iy) - 310.15*N.MM.dS(ic, ix, iy)/1000;
            % 0.10, not 0.06. DeltaH is published to 0.1 kcal/mol and DeltaS to 0.1 cal/(mol K),
            % so reconstructing DeltaG37 from them carries a rounding floor of about
            % 0.05 + 310.15*0.05/1000 = 0.065 kcal/mol before anything is wrong. Measured across the
            % full 48 on 17-10-2026: mean residual +0.0000, sd 0.026, worst 0.063, with two rows
            % above 0.06. A tolerance of 0.06 therefore flagged correct data. 0.10 still catches a
            % sign error (which moves DeltaG37 by twice the T*dS term, several kcal/mol) and DeltaS
            % entered in kcal rather than cal (which moves it by a factor of a thousand).
            if abs(recon - N.MM.dG37(ic, ix, iy)) > 0.10
                badrow{end+1} = sprintf('%s%s/%s%s  dG37 %.2f but dH - T dS gives %.2f', ...
                    bases(ic), bases(ix), Comp(bases(ic)), bases(iy), ...
                    N.MM.dG37(ic, ix, iy), recon); %#ok<AGROW>
            end
        end
    end
end
if ~isempty(badrow)
    warning('NNParams:InconsistentMismatchRow', ...
        ['%d filled mismatch row(s) are internally inconsistent at 37 C by more than 0.10 kcal/mol. ' ...
        'Either a sign or a unit is wrong, or DeltaS was entered in kcal rather than cal. ' ...
        'DeltaS belongs in cal/(K mol), matching the Watson-Crick column.\n   %s'], ...
        numel(badrow), strjoin(badrow, sprintf('\n   ')));
end

%% ---- route assignment -------------------------------------------------------------------------------
N.MM.Route = uint8(zeros(4, 4, 4));
N.MM.Route(isfinite(N.MM.dG37)) = 3;
N.MM.Route(N.MM.HasHS) = 1;

% Only reach for ViennaParams when something actually needs route (b). Once every entry carries
% measured DeltaH and DeltaS the int11 slope has no job, and loading the parameter file anyway
% emits a warning about an asymmetry in a stack table this module never reads. Added 17-10-2026,
% when the CSV was filled and that warning started firing for no reason.
nNeedSlope = nnz(N.MM.Route == 3);
V = [];
if opts.UseVienna && (nNeedSlope > 0)
    try
        if isempty(opts.ViennaFile)
            V = ViennaParams([], struct('Verbose', false));
        else
            V = ViennaParams(opts.ViennaFile, struct('Verbose', false));
        end
    catch ME
        warning('NNParams:NoVienna', ...
            ['ViennaParams failed (%s), so route (b) is unavailable and every unfilled mismatch ' ...
            'falls to route (c), DeltaG37 with no temperature scaling.'], ME.message);
    end
end
if ~isempty(V)
    N.MM.Route(N.MM.Route == 3) = 2;
    N.MM.ViennaSlope = MakeViennaSlope(V);
else
    N.MM.ViennaSlope = @(W, X, Y, Vb, TdegC) NaN;
end

miss = {}; mstep = {}; mmm = {};
for ic = 1:4
    for ix = 1:4
        for iy = 1:4
            if ~isfinite(N.MM.dG37(ic, ix, iy)), continue; end
            if N.MM.HasHS(ic, ix, iy), continue; end
            miss{end+1} = sprintf('%s,%s,%s', bases(ic), bases(ix), bases(iy)); %#ok<AGROW>
            mstep{end+1} = sprintf('%s%s/%s%s', bases(ic), bases(ix), Comp(bases(ic)), bases(iy)); %#ok<AGROW>
            mmm{end+1} = sprintf('%s.%s', bases(ix), bases(iy)); %#ok<AGROW>
        end
    end
end
N.MM.Missing = table(string(miss(:)), string(mstep(:)), string(mmm(:)), ...
    'VariableNames', {'Index', 'Step', 'Mismatch'});

%% ---- dangling ends, terminal mismatches and loop initiation ------------------------------------------
% THE MOTIFS A PRIMER ACTUALLY FORMS, BEYOND A BLUNT FULLY PAIRED DUPLEX.
%
% A primer annealed to a long template does not end in a vacuum: unpaired template hangs off both
% ends of the duplex and stacks on it. A hairpin does not close on nothing: the first unpaired pair
% stacks on the stem. A mispriming register may carry an unpaired base on one strand. Each motif
% carries measured free energy that the Watson-Crick table cannot express.
%
% DANGLING ENDS are complete. Bommarito, Peyret and SantaLucia measured all 32 single-base dangling
% ends; nothing in that table is extrapolated. The contribution reaches -0.96 kcal/mol, which is
% comparable to a whole A-T pair, so omitting it biases every duplex free energy the same way.
% B-form DNA shows no systematic preference between the 3' and 5' end and two of the values are
% destabilising, neither of which is true of A-form RNA, so an RNA table is not a substitute.
%
% TERMINAL MISMATCHES are required by Equations 9 and 10 of the review for hairpin stability, which
% means a hairpin free energy computed without them is incomplete rather than approximate. 46 of
% the 48 are stabilising, so replacing the table with a conservative constant makes hairpins look
% weaker than they are. The full measured set is attributed in the review to a manuscript in
% preparation, so there is no primary paper to transcribe and the values come from Primer3.
%
% LOOP INITIATION is Table 4. All loop enthalpies are zero by construction, so a loop penalty is
% purely entropic and does not move with temperature. Lengths off the table are reached by
% Jacobson-Stockmayer with the DNA kinetics coefficient 2.44.
%
% WHAT IS DELIBERATELY ABSENT. DNA tandem mismatches, the 2x2 internal loops that RNA folding calls
% int22, are not here because they have not been measured. The review states that DNA parameters
% exist only for tandem G.T and that other 2x2 sequences remain undetermined. RNA has a complete
% periodic table of them (Wu, McDowell & Turner 1995) and transferring it would be a modelling
% choice presented as a measurement, which this module does not do.
%
% EVERY NUMBER IS CHECKED RATHER THAN TRUSTED. The dangling-end file carries an independently
% transcribed Published_dH column from Table 3 and the loader asserts agreement on all 32 rows. The
% loop file is asserted against the review's own worked bulge example. The terminal-mismatch file is
% asserted against the two sequences the review names as the extremes of its range.
% MOTIFCAP IS A UNIT TEST WEARING A PHYSICS COAT.
% No single measured motif increment in this suite reaches three kcal/mol at 37 C: the dangling
% ends span +0.48 to -0.96, the terminal mismatches -1.19 to +0.06, and the special loops -2.0 to
% +1.5. Entering dS in the Primer3 files' scaled units instead of cal/(K mol) multiplies the
% T*dS term by a thousand and throws dG37 into the hundreds, so any bound between about three and
% a hundred separates the two cases. Three is used because it is the tightest value that passes
% every row of all three tables as shipped, which makes the assertion informative rather than
% decorative. It is a SCALE test and nothing finer; the range and provenance tests below do the
% rest.
MOTIFCAP = 3.0;                 % kcal/mol, largest |dG37| any one motif may reconstruct to

N.Dangle.dH = nan(4, 4, 4, 2);   % (paired, unpaired neighbour, partner of paired, 1 = 3' side, 2 = 5')
N.Dangle.dS = nan(4, 4, 4, 2);
if isfile(danglingFile)
    [hdrD, cellsD] = ReadSimpleCsv(danglingFile);
    cD = MatchColumns(hdrD, {'Side', 'Paired', 'Dangling', 'Complement', ...
        'dH_kcal_per_mol', 'dS_cal_per_K_per_mol', 'Published_dH_kcal_per_mol', ...
        'Published_dG37_kcal_per_mol'}, danglingFile);
    % TOLERANCE ON THE ENTROPY CHECK, DERIVED RATHER THAN CHOSEN.
    % Table 3 prints no dS column. Its footnote states that dS is to be recovered from dH and
    % dG37, so the only way to test the stored dS against the source is to reconstruct
    % dG37 = dH - 310.15*dS/1000 and compare. The residual is carried in ENTROPY UNITS, which is
    % where the quantity under test actually lives, and the bound follows from the three printed
    % precisions rather than from taste:
    %       dH   printed to 0.1  kcal/mol  ->  0.050  / 0.31015  =  0.1612 e.u.
    %       dG37 printed to 0.01 kcal/mol  ->  0.005  / 0.31015  =  0.0161 e.u.
    %       dS   stored  to 0.1  e.u.      ->                       0.0500 e.u.
    %                                                   worst case  0.2273 e.u.
    % No correctly transcribed row can exceed that. The worst of the 32 as shipped is 0.1967 e.u.
    % at 5' T/A with a T dangling end, which is 87 percent of the bound, so the test has almost no
    % slack and is informative rather than decorative. A tolerance of 0.02 kcal/mol, which looks
    % like the natural choice against a table printed to 0.01, rejects 20 of the 32 correct rows.
    % In this form the check catches any slip of 0.45 e.u. or more wherever the row started, which
    % covers every transcription error in the integer place and most in the first decimal.
    DANGLEDSTOL = 0.25;             % entropy units
    worstD = 0; nD = 0; worstS = 0; nS = 0; worstSrow = '';
    for r = 1:size(cellsD, 1)
        side = ToNum(cellsD{r, cD.Side});
        ip = strfind(bases, upper(strtrim(cellsD{r, cD.Paired})));
        ix = strfind(bases, upper(strtrim(cellsD{r, cD.Dangling})));
        ic = strfind(bases, upper(strtrim(cellsD{r, cD.Complement})));
        if isempty(ip)||isempty(ix)||isempty(ic)||~ismember(side, [3 5])
            error('NNParams:BadDangle', 'Row %d of %s is malformed.', r + 1, danglingFile);
        end
        slot = 1 + (side == 5);
        h = ToNum(cellsD{r, cD.dH_kcal_per_mol});
        N.Dangle.dH(ip, ix, ic, slot) = h;
        N.Dangle.dS(ip, ix, ic, slot) = ToNum(cellsD{r, cD.dS_cal_per_K_per_mol});
        pubh = ToNum(cellsD{r, cD.Published_dH_kcal_per_mol});
        if isfinite(pubh)
            worstD = max(worstD, abs(h - pubh)); nD = nD + 1;
        end
        sVal = ToNum(cellsD{r, cD.dS_cal_per_K_per_mol});
        pubg = ToNum(cellsD{r, cD.Published_dG37_kcal_per_mol});
        if isfinite(pubg)&&isfinite(sVal)&&isfinite(h)
            resS = abs((h - 310.15*sVal/1000 - pubg)*1000/310.15);
            if resS > worstS
                worstS = resS;
                worstSrow = sprintf('%d-prime %s/%s with %s dangling', side, ...
                    upper(strtrim(cellsD{r, cD.Paired})), upper(strtrim(cellsD{r, cD.Complement})), ...
                    upper(strtrim(cellsD{r, cD.Dangling})));
            end
            nS = nS + 1;
        end
    end
    % The two columns were produced independently, one from the Primer3 distribution and one by
    % transcription from Table 3. If they disagree, one of them was mistyped and the whole table is
    % untrustworthy, so this stops rather than warns.
    if (nD > 0) && (worstD > 0.01)
        error('NNParams:DangleMismatch', ...
            ['The dangling-end enthalpies in %s disagree with the Published_dH column by up to ' ...
            '%.3f kcal/mol over %d rows. One of the two transcriptions is wrong.'], ...
            danglingFile, worstD, nD);
    end
    % THE ENTROPY COLUMN, CHECKED THROUGH THE PUBLISHED FREE ENERGY.
    % dH is already pinned by the assertion above, so a reconstruction of dG37 that disagrees with
    % the printed value by more than the rounding bound can only be the dS column. This stops for
    % the same reason the enthalpy check does: one mistyped entropy makes every dangling-end
    % contribution at every temperature other than 37 C wrong, by an amount that grows with the
    % distance from 37, and nothing downstream would notice.
    if (nS > 0)&&(worstS > DANGLEDSTOL)
        error('NNParams:DanglingEntropyMismatch', ...
            ['The dangling-end entropies in %s are inconsistent with the Published_dG37 column ' ...
            'by up to %.4f entropy units over %d rows, against a tolerance of %.2f derived from ' ...
            'the printed precisions of the source. The worst row is %s. Since the enthalpy ' ...
            'column has already been asserted against Published_dH, the discrepancy is in the ' ...
            'dS column.'], danglingFile, worstS, nS, DANGLEDSTOL, worstSrow);
    end
    AssertMotifScale(N.Dangle.dH, N.Dangle.dS, MOTIFCAP, danglingFile, 'dangling end', {});
    N.Dangle.Checked = nD;
    N.Dangle.CheckedS = nS;
    N.Dangle.WorstS = worstS;
else
    N.Dangle.Checked = 0;
    N.Dangle.CheckedS = 0;
    N.Dangle.WorstS = NaN;
    warning('NNParams:NoDanglingFile', ...
        ['dangling_nn_params.csv was not found at %s; dangling ends will read as NaN and ' ...
        'dGDangle returns zero.'], danglingFile);
end

N.TM.dH = nan(4, 4, 4, 4);       % (i1, i2, j1, j2): 5'-i1 i2-3' over 3'-j1 j2-5'
N.TM.dS = nan(4, 4, 4, 4);
N.TM.Suspect = false(4, 4, 4, 4);
if isfile(termMismFile)
    [hdrT, cellsT] = ReadSimpleCsv(termMismFile);
    cT = MatchColumns(hdrT, {'i1', 'i2', 'j1', 'j2', ...
        'dH_kcal_per_mol', 'dS_cal_per_K_per_mol', 'Suspect'}, termMismFile);
    for r = 1:size(cellsT, 1)
        q = zeros(1, 4);
        nm = {'i1', 'i2', 'j1', 'j2'};
        for t = 1:4
            hit = strfind(bases, upper(strtrim(cellsT{r, cT.(nm{t})})));
            if isempty(hit)
                error('NNParams:BadTermMism', ...
                    'Row %d of %s carries ''%s'' in column %s, which is not one of A, C, G, T.', ...
                    r + 1, termMismFile, strtrim(cellsT{r, cT.(nm{t})}), nm{t});
            end
            q(t) = hit(1);
        end
        N.TM.dH(q(1), q(2), q(3), q(4)) = ToNum(cellsT{r, cT.dH_kcal_per_mol});
        N.TM.dS(q(1), q(2), q(3), q(4)) = ToNum(cellsT{r, cT.dS_cal_per_K_per_mol});
        N.TM.Suspect(q(1), q(2), q(3), q(4)) = ToNum(cellsT{r, cT.Suspect}) ~= 0;
    end
    % The review names CG/GA and AC/TC as the extremes of the measured range. If this file is the
    % set the review describes, those two must come out as the extremes here.
    g = N.TM.dH - 310.15*N.TM.dS/1000;
    ok = g(~isnan(g) & ~N.TM.Suspect);
    if ~isempty(ok)
        lo = g(2, 3, 3, 1);          % CG/GA
        hi = g(1, 2, 4, 2);          % AC/TC
        % THE NaN GUARD IS THE POINT OF THIS CHECK, NOT AN AFTERTHOUGHT.
        % The obvious corruption of a four-index table is a transposed column pair, and that
        % permutes which key holds which value while leaving the SET of values untouched. A test on
        % min and max alone is therefore blind to it. What a transposition does do is send CG/GA to
        % a slot the file never fills, because its mirror CG/AG has no Watson-Crick closing pair and
        % is not shipped. The slot then reads NaN, and every inequality against NaN is false, so a
        % bare range test passes in silence. Testing for NaN first is what catches the transposition.
        if ~isfinite(lo)||~isfinite(hi)
            error('NNParams:TermMismIndex', ...
                ['%s has no entry at CG/GA or at AC/TC. The review names both as the extremes of ' ...
                'the measured range, so both must be present. The most likely cause is that two ' ...
                'of the four base columns are transposed.'], termMismFile);
        end
        if (abs(lo - min(ok)) > 1E-9)||(abs(hi - max(ok)) > 1E-9)
            error('NNParams:TermMismRange', ...
                ['The terminal-mismatch range in %s runs %.2f to %.2f, but the review gives its ' ...
                'extremes as CG/GA (%.2f here) and AC/TC (%.2f here). Check the index order.'], ...
                termMismFile, min(ok), max(ok), lo, hi);
        end
        % THE RANGE TEST ABOVE IS SCALE-FREE AND THIS ONE IS NOT. Multiplying every dS by a
        % thousand leaves CG/GA the minimum and AC/TC the maximum, so the extremes test passes a
        % unit error untouched. Suspect rows are excluded here as they are there, because the one
        % flagged entry is already known to be unphysical and is clamped at lookup.
        if max(abs(ok)) > MOTIFCAP
            [~, kworst] = max(abs(ok));
            error('NNParams:TermMismScale', ...
                ['A terminal mismatch in %s reconstructs to %.1f kcal/mol at 37 C, past the ' ...
                '%.1f kcal/mol any single motif can reach. dH must be kcal/mol and dS ' ...
                'cal/(K mol).'], termMismFile, ok(kworst), MOTIFCAP);
        end
        N.TM.MaxPhysical = max(ok);  % the clamp dGTerminalMismatch applies to a suspect entry
        % MEAN OVER THE MEASURED SET, for the one caller that needs a terminal mismatch without
        % being able to say which one. InternalLoopBridge closes a two-sided internal loop whose
        % flanking loop bases lie between two alignment blocks and are therefore not in its scope,
        % so it cannot index the table; it uses this instead. Measured over the 47 non-suspect
        % entries on 29-10-2026: -0.646 kcal/mol, against a range of -1.193 to +0.055.
        N.TM.MeanG37 = mean(ok);
    else
        N.TM.MaxPhysical = 0;
        N.TM.MeanG37 = 0;
    end
else
    N.TM.MaxPhysical = 0;
    N.TM.MeanG37 = 0;
    warning('NNParams:NoTermMismFile', ...
        ['terminal_mismatch_nn_params.csv was not found at %s; hairpin and internal-loop closures ' ...
        'fall back to the conservative constant.'], termMismFile);
end

N.Loop.N = []; N.Loop.Internal = []; N.Loop.Bulge = []; N.Loop.Hairpin = [];
if isfile(loopFile)
    [hdrL, cellsL] = ReadSimpleCsv(loopFile);
    cL = MatchColumns(hdrL, {'LoopLength', 'Internal_dG37', 'Bulge_dG37', 'Hairpin_dG37'}, loopFile);
    nRow = size(cellsL, 1);
    N.Loop.N = zeros(1, nRow); N.Loop.Internal = nan(1, nRow);
    N.Loop.Bulge = nan(1, nRow); N.Loop.Hairpin = nan(1, nRow);
    for r = 1:nRow
        N.Loop.N(r)        = ToNum(cellsL{r, cL.LoopLength});
        N.Loop.Internal(r) = ToNum(cellsL{r, cL.Internal_dG37});
        N.Loop.Bulge(r)    = ToNum(cellsL{r, cL.Bulge_dG37});
        N.Loop.Hairpin(r)  = ToNum(cellsL{r, cL.Hairpin_dG37});
    end
    [N.Loop.N, ord] = sort(N.Loop.N);
    N.Loop.Internal = N.Loop.Internal(ord);
    N.Loop.Bulge    = N.Loop.Bulge(ord);
    N.Loop.Hairpin  = N.Loop.Hairpin(ord);
    % The review works one bulge by hand: 5'-CAT-3' over 5'-AG-3' gives +4.0 - 1.28 + 0.5 = +3.22,
    % where -1.28 is the intervening CT/GA stack and +0.5 the A-T closure. Reproducing that number
    % tests the bulge column, the closure penalty and the sign convention in one line.
    b1 = N.Loop.Bulge(N.Loop.N == 1);
    if ~isempty(b1) && (abs((b1 - 1.28 + 0.5) - 3.22) > 0.005)
        error('NNParams:LoopTableCheck', ...
            ['The worked bulge example in SantaLucia and Hicks gives +3.22 kcal/mol; this table ' ...
            'gives %.2f. The bulge column or the closure penalty is wrong.'], b1 - 1.28 + 0.5);
    end
else
    warning('NNParams:NoLoopFile', ...
        'loop_length_params.csv was not found at %s; bulge and hairpin loops have no length term.', ...
        loopFile);
end
% SPECIAL SHORT LOOPS. Certain triloop and tetraloop sequences are far more stable than their
% length alone implies, and the review gives them a sequence-specific bonus in Equations 8 and 9.
% Every measured DNA triloop bonus is a G-N-A loop. These are not a refinement: GGAAC is worth
% -2.0 kcal/mol, which is larger than the entire terminal-mismatch range, and such loops arise by
% chance often enough in primers that ignoring them mis-ranks real candidates. The review places
% the values in supplementary material, so they are taken from the Primer3 distribution.
% Keyed by the closing pair's top base, the loop bases 5' to 3', then the closing pair's bottom base.
N.Loop.SpecialKey = {};
N.Loop.SpecialdH  = [];
N.Loop.SpecialdS  = [];
if isfile(specialFile)
    [hdrS, cellsS] = ReadSimpleCsv(specialFile);
    cS = MatchColumns(hdrS, {'Loop', 'Key', 'dH_kcal_per_mol', 'dS_cal_per_K_per_mol'}, specialFile);
    nS = size(cellsS, 1);
    N.Loop.SpecialKey = cell(nS, 1);
    N.Loop.SpecialdH  = zeros(nS, 1);
    N.Loop.SpecialdS  = zeros(nS, 1);
    for r = 1:nS
        N.Loop.SpecialKey{r} = upper(strtrim(cellsS{r, cS.Key}));
        N.Loop.SpecialdH(r)  = ToNum(cellsS{r, cS.dH_kcal_per_mol});
        N.Loop.SpecialdS(r)  = ToNum(cellsS{r, cS.dS_cal_per_K_per_mol});
    end
    % Every DNA triloop bonus measured is a G-N-A loop, so the second and fourth characters of a
    % five-character key are fixed. If that stops being true the table has been edited or the key
    % convention has shifted, and the hairpin energies would silently change.
    for r = 1:nS
        k = N.Loop.SpecialKey{r};
        if (numel(k) == 5) && ~(k(2) == 'G' && k(4) == 'A')
            error('NNParams:TriloopKey', ...
                ['Triloop key ''%s'' in %s is not a G-N-A loop. Every measured DNA triloop bonus ' ...
                'is, so either the key convention or the table is wrong.'], k, specialFile);
        end
    end
    AssertMotifScale(N.Loop.SpecialdH, N.Loop.SpecialdS, MOTIFCAP, specialFile, ...
        'triloop or tetraloop bonus', N.Loop.SpecialKey);
else
    warning('NNParams:NoSpecialLoopFile', ...
        ['special_loop_params.csv was not found at %s; stable triloops and tetraloops will be ' ...
        'scored on length alone, which understates their stability badly.'], specialFile);
end

N.Loop.ATClosureG = 0.5;      % review Equations 8 and 13, loop closed by an A-T pair
N.Loop.JSCoeff    = 2.44;     % Jacobson-Stockmayer, DNA kinetics value preferred over 1.75
N.Rgas            = 1.9872E-3;  % kcal/mol/K, the value the rest of the suite uses

%% ---- salt ---------------------------------------------------------------------------------------------
Ka = 3E4;    % 1/M, the 1:1 Mg.dNTP association constant, after Owczarzy et al. 2008
% THE SPECIATION ITSELF LIVES IN FreeMgSpeciated.m, NOT HERE.
% It used to be a local function at the foot of this file, which put it out of reach of any caller
% that does not load the parameter struct. TmAmplicon.m is one, and it had been carrying its own
% max([Mg] - [dNTP], 0) stand-in instead. The two agree to under one percent at ordinary PCR
% compositions and part company once dNTP approaches magnesium, where the stand-in clamps at zero
% and the quadratic does not, which was worth up to 4.6 degrees in a reported amplicon melting
% temperature. One file now holds the root and every caller takes it from there. FreeMgSpeciated
% also returns the root at Ka = 1E4 and 1E5 /M, so a caller that forms a temperature can report how
% much of its answer this constant is setting.
N.Salt.FreeMg = @(biv, dNTP) FreeMgSpeciated(biv, dNTP, Ka);
N.Salt.Owczarzy = @(h, gc, nbp, monv, biv, dNTP) ...
    OwczarzyShift(h, gc, nbp, monv, FreeMgSpeciated(biv, dNTP, Ka));
N.Salt.PerPhosphate = @(nbp, monvEq) 0.368*(nbp - 1).*log(monvEq);
% UNITS. The von Ahsen coefficient 120 is calibrated for concentrations in MILLIMOLAR, so it
% cannot be applied to molar quantities. The conversion is written out rather than folded into a
% single constant, because folding it in is how the error below happened in the first place: the
% first version of this line read monv + 120*sqrt(FreeMgSpeciated(...)) with both in molar and
% returned 3.529 M where 0.150 M is correct, a factor of 23. Through the per-phosphate correction
% that is worth about 7 kcal/mol on a 20-mer, in the direction that makes every mismatched duplex
% look more stable than it is.
N.Salt.MonovalentEquivalent = @(monv, biv, dNTP) ...
    1E-3*( 1E3*monv + 120*sqrt(1E3*max(FreeMgSpeciated(biv, dNTP, Ka), 0)) );

% Unit self-check against the worked example already in the header of DNAAccessibility.m: 40 mM
% sodium with 2.5 mM magnesium and 0.8 mM dNTP comes to about 196 mM. That number was derived there
% with the crude [Mg] - [dNTP] stand-in for free magnesium, which at this composition happens to sit
% within 1 percent of the speciation used here, so the two routes should agree closely. Any future
% edit that reintroduces a unit slip fails here rather than in a free energy three files away.
eqRef = N.Salt.MonovalentEquivalent(40E-3, 2.5E-3, 0.8E-3);
if ~isfinite(eqRef)||(abs(eqRef - 0.197) > 0.010)
    error('NNParams:MonovalentEquivalentUnits', ...
        ['The monovalent equivalent self-check returned %.4f M at 40 mM Na, 2.5 mM Mg and ' ...
        '0.8 mM dNTP, where 0.197 M is expected from the worked example in DNAAccessibility.m. ' ...
        'This is almost always a millimolar against molar slip in the von Ahsen term.'], eqRef);
end

%% ---- metadata and summary ------------------------------------------------------------------------------
if isempty(V)
    viennaPath = '';
else
    viennaPath = V.Meta.Path;
end
N.Meta = struct('Key', key, 'MismatchFile', opts.MismatchFile, ...
    'ViennaLoaded', ~isempty(V), 'ViennaPath', viennaPath, ...
    'ViennaStamp', ViennaStamp(viennaPath), ...
    'Units', 'dH kcal/mol, dS cal/(K mol)');
CACHE = N;

if opts.Verbose
    nA = nnz(N.MM.Route == 1); nB = nnz(N.MM.Route == 2); nC = nnz(N.MM.Route == 3);
    fprintf('\nNNParams loaded.\n');
    fprintf('  Watson-Crick: 16 steps plus 3 correction rows, from the verified block.\n');
    fprintf('  Mismatch entries: %d of 48 measured (route a), %d on the Vienna slope (route b), %d unscaled (route c).\n', ...
        nA, nB, nC);
    if nC > 0
        warning('NNParams:UnscaledMismatches', ...
            ['%d mismatch entries have no temperature dependence at all and will be applied at ' ...
            'their 37 C value. This is the defect Stage 1 exists to remove. Fill dH and dS in %s, ' ...
            'or make ViennaParams available for the slope route.'], nC, opts.MismatchFile);
    elseif nB > 0
        fprintf('  Route (b) is in use for %d entries. Keep the Vienna dependency documented in any write-up.\n', nB);
    end
    fprintf('  Free Mg  (2.0 mM total, 1.2 mM dNTP) : %6.3f mM   [crude Mg-dNTP gives %.3f mM]\n', ...
        1E3*N.Salt.FreeMg(2E-3, 1.2E-3), 1E3*(2E-3 - 1.2E-3));
    fprintf('  Monovalent equivalent (40 mM Na)      : %6.1f mM   [von Ahsen, mM-calibrated]\n', ...
        1E3*N.Salt.MonovalentEquivalent(40E-3, 2E-3, 1.2E-3));
    fprintf('  Unit self-check at 2.5 mM Mg, 0.8 mM dNTP: %6.1f mM   [expect 197, from DNAAccessibility.m]\n\n', ...
        1E3*eqRef);
end
end

%% ____________________________________________________________________________________________________________
function v = ToNum(x)
if isnumeric(x)
    v = double(x);
    if isempty(v), v = NaN; end
    return
end
% char() rather than string(): it accepts a char row, a string scalar and a 1x1 cellstr alike,
% and str2double already returns NaN for anything it cannot parse. The previous form tested
% s == "" against a string scalar, which assumed one specific input type for no benefit.
s = strtrim(char(x));
if isempty(s)
    v = NaN;
else
    v = str2double(s);
end
end

%% ____________________________________________________________________________________________________________
function AssertMotifScale(dH, dS, cap, fileName, what, labels)
%% ASSERTMOTIFSCALE  Refuse a motif table whose dH and dS are not in the stated units
% dG37 = dH - 310.15*dS/1000 is formed for every filled entry and compared against cap. The test
% exists because a column scaled by a thousand is invisible to every other check in this file: a
% range test, an extremes test and a key-convention test are all scale-free, and the consumer
% happily returns a hairpin at minus five hundred kcal/mol. Entries that are NaN are unfilled
% slots and are skipped. labels is a cell array naming each entry, or {} to report linear indices.
g = dH - 310.15*dS/1000;
bad = find(isfinite(g) & (abs(g) > cap));
if isempty(bad)
    return
end
lines = cell(1, min(numel(bad), 6));
for q = 1:numel(lines)
    k = bad(q);
    if isempty(labels)
        nm = sprintf('entry %d', k);
    else
        nm = labels{k};
    end
    lines{q} = sprintf('%s: dH %.4g, dS %.4g, dG37 %.1f', nm, dH(k), dS(k), g(k));
end
more = '';
if numel(bad) > numel(lines)
    more = sprintf(', and %d further entries', numel(bad) - numel(lines));
end
error('NNParams:MotifScale', ...
    ['%d %s entries in %s reconstruct past %.1f kcal/mol at 37 C. dH must be kcal/mol and ' ...
    'dS cal/(K mol); a column carrying the Primer3 files'' scaled entropy looks exactly like ' ...
    'this. Offenders: %s%s'], ...
    numel(bad), what, fileName, cap, strjoin(lines, ' | '), more);
end

%% ____________________________________________________________________________________________________________
function c = Comp(b)
switch b
    case 'A', c = 'T';
    case 'T', c = 'A';
    case 'C', c = 'G';
    case 'G', c = 'C';
    otherwise, c = 'N';
end
end

%% ____________________________________________________________________________________________________________
function ds = OwczarzyShift(h, gc, nbp, monv, bivFree)
% Entropy shift equivalent to the Owczarzy correction on 1/Tm, in cal/(K mol). FOR FULLY PAIRED
% DUPLEXES ONLY. Use N.Salt.PerPhosphate for anything carrying a mismatch, and read the salt section
% of the header for why the distinction is not cosmetic.
catrat = sqrt(bivFree)./monv;
if catrat > 6.0
    a = 3.92; d = 1.42; g = 8.31;
elseif catrat >= 0.22
    a = 3.92*(0.843 - 0.352*sqrt(monv)*log(monv));
    d = 1.42*(1.279 - 4.03E-3*log(monv) - 8.03E-3*(log(monv))^2);
    g = 8.31*(0.486 - 0.258*log(monv) + 5.25E-3*(log(monv))^3);
else
    a = 0; d = 0; g = 0;
end
if nbp < 2
    error('NNParams:SegmentTooShort', ...
        'The divalent branch divides by 2*(N-1); a segment of %d base pair(s) has no nearest neighbour.', nbp);
end
if catrat < 0.22
    f = (4.29*gc - 3.95)*log(monv) + 0.940*(log(monv))^2;
else
    f = a - 0.911*log(bivFree) + gc*(6.26 + d*log(bivFree)) ...
        + (1/(2*(nbp - 1)))*(-48.2 + 52.5*log(bivFree) + g*(log(bivFree))^2);
end
ds = h*f*1E-5*1000;     % h is kcal/mol, the correction is on 1/Tm, output in cal/(K mol)
end

%% ____________________________________________________________________________________________________________
function fh = MakeViennaSlope(V)
% Route (b). Temperature correction for an ISOLATED single mismatch, taken as the change in the
% Vienna 1x1 internal loop energy between 37 C and the working temperature, with the 37 C magnitude
% divided out so only the slope transfers.
%
%   W is the closing base on the top strand 5' of the mismatch, V the closing base on the top strand
%   3' of it, X and Y the mismatch bases on the top and bottom strands. The loop is then closed by
%   pair (W, comp(W)) on one side and by (V, comp(V)) read from the other side.
%
% Returns NaN whenever any index is outside the Vienna alphabet, which the caller must treat as
% route (c) rather than as zero.
fh = @slope;
    function dg = slope(W, X, Y, Vb, TdegC)
        dg = NaN;
        p1 = [W, CompL(W)];
        p2 = [CompL(Vb), Vb];
        if ~isfield(V.PairIndex, p1)||~isfield(V.PairIndex, p2), return; end
        if ~isfield(V.BaseIndex, X)||~isfield(V.BaseIndex, Y), return; end
        i1 = V.PairIndex.(p1); i2 = V.PairIndex.(p2);
        b1 = V.BaseIndex.(X);  b2 = V.BaseIndex.(Y);
        g37 = V.int11.dG37(i1, i2, b1, b2);
        hh  = V.int11.dH(  i1, i2, b1, b2);
        if ~isfinite(g37)||~isfinite(hh), return; end
        dg = V.AtTemp(g37, hh, TdegC) - g37;
    end
    function c = CompL(b)
        switch b
            case 'A', c = 'T';
            case 'T', c = 'A';
            case 'C', c = 'G';
            case 'G', c = 'C';
            otherwise, c = 'N';
        end
    end
end

%% ____________________________________________________________________________________________________________
function [hdr, cells] = ReadSimpleCsv(file)
% Minimal delimited-text reader for a small, flat table. No auto-detected import options, no header
% renaming, no release-dependent behaviour. Handles a byte-order mark, any line ending, optional
% surrounding double quotes, and a delimiter chosen by consistency rather than by guesswork.
txt = fileread(file);
if ~isempty(txt) && (double(txt(1)) == 65279)
    txt = txt(2:end);                       % UTF-8 byte-order mark, which Excel adds on save
end
raw = regexp(txt, '\r\n|\n|\r', 'split');
raw = raw(~cellfun(@(x) isempty(strtrim(x)), raw));
% A leading # marks provenance, which belongs beside the numbers rather than in a separate note.
% Dropped before the header is read, so a commented file and a bare one parse identically.
raw = raw(~cellfun(@(x) strncmp(strtrim(x), '#', 1), raw));
if numel(raw) < 2
    error('NNParams:EmptyCsv', '%s holds a header but no data rows.', file);
end

% Pick the delimiter that splits the header and the first data row into the same number of fields,
% and into more than one. Comma first, since that is what ships.
cand = {',', ';', sprintf('\t'), '|'};
delim = '';
for k = 1:numel(cand)
    nh = numel(SplitQuoted(raw{1}, cand{k}));
    nd = numel(SplitQuoted(raw{2}, cand{k}));
    if (nh > 1) && (nh == nd)
        delim = cand{k}; break
    end
end
if isempty(delim)
    error('NNParams:NoDelimiter', ...
        ['Could not find a delimiter that splits the header and the first data row of %s into the ' ...
        'same number of fields.\n  header: %s\n  row 1 : %s'], file, raw{1}, raw{2});
end

hdr = CleanFields(SplitQuoted(raw{1}, delim));
nCol = numel(hdr);
cells = cell(numel(raw) - 1, nCol);
for i = 2:numel(raw)
    f = CleanFields(SplitQuoted(raw{i}, delim));
    if numel(f) ~= nCol
        error('NNParams:RaggedCsv', ...
            ['Row %d of %s holds %d fields where the header declares %d. An unbalanced double ' ...
            'quote would do this, since a delimiter inside quotes is not a separator.\n  %s'], ...
            i, file, numel(f), nCol, raw{i});
    end
    cells(i - 1, :) = f;
end
end

%% ____________________________________________________________________________________________________________
function f = CleanFields(f)
for k = 1:numel(f)
    s = strtrim(f{k});
    if (numel(s) >= 2) && (s(1) == '"') && (s(end) == '"')
        s = strtrim(s(2:end-1));
    end
    f{k} = s;
end
end

%% ____________________________________________________________________________________________________________
function col = MatchColumns(hdr, want, file)
% Case-insensitive match after stripping everything that is not a letter or a digit, so a header
% mangled into "dG37 kcal per mol" or "dG37_kcal_per_mol_1" still resolves. The error names the
% header row as read, which is the one piece of information a missing-column message must carry.
norm = @(s) lower(regexprep(s, '[^a-zA-Z0-9]', ''));
hn = cellfun(norm, hdr, 'UniformOutput', false);
col = struct();
for k = 1:numel(want)
    ix = find(strcmp(hn, norm(want{k})), 1);
    if isempty(ix)
        error('NNParams:MissingColumn', ...
            ['Column "%s" was not found in %s.\nHeader row as read: %s\nMatching ignores case ' ...
            'and punctuation, so this means the column is genuinely absent or the file was split ' ...
            'on the wrong delimiter.'], want{k}, file, strjoin(hdr, ' | '));
    end
    col.(want{k}) = ix;
end
end

%% ____________________________________________________________________________________________________________
function f = SplitQuoted(line, delim)
% Split one line on delim, honouring RFC 4180 double quotes: a delimiter inside a quoted field is
% not a separator, and a doubled quote inside a quoted field is one literal quote.
%
% Added 17-10-2026. The previous version used strsplit, which has no notion of quoting, so a field
% containing the delimiter split into two and the row came back with more fields than the header.
% That is not hypothetical: Excel quotes any field holding a comma, so the first time this file was
% edited in a spreadsheet and saved it would have failed to load, with an error naming the
% delimiter rather than the quoting. Found by running NNParams under Octave.
f = {}; cur = ''; inq = false; i = 1; n = numel(line);
while i <= n
    c = line(i);
    if inq
        if c == '"'
            if (i < n) && (line(i+1) == '"')
                cur(end+1) = '"'; i = i + 2; continue    %#ok<AGROW>
            end
            inq = false; i = i + 1; continue
        end
        cur(end+1) = c; i = i + 1;                       %#ok<AGROW>
    else
        if c == '"'
            inq = true; i = i + 1; continue
        end
        if strncmp(line(i:end), delim, numel(delim))
            f{end+1} = cur; cur = ''; i = i + numel(delim); continue   %#ok<AGROW>
        end
        cur(end+1) = c; i = i + 1;                       %#ok<AGROW>
    end
end
f{end+1} = cur;
end

%% ____________________________________________________________________________________________________________
function s = ViennaStamp(p)
%% VIENNASTAMP  Path, size and modification time of the Vienna parameter file, as one string
% The Vienna file cannot enter the cache key the way the five CSV tables do. Its path is resolved
% inside ViennaParams when the caller leaves opts.ViennaFile empty, and that call is reached only
% after the mismatch table has been parsed and the slope routes counted, which is well past the
% point where the key is formed. Stamping the file after the fact and re-testing that stamp on a
% cache hit closes the same hole from the other side: a swapped or re-edited parameter file fails
% the comparison and the structure is rebuilt.
s = '';
if isempty(p), return; end
if ~isfile(p), s = 'absent'; return; end
d = dir(p);
s = sprintf('%s|%d|%.10f', p, d.bytes, d.datenum);
end
