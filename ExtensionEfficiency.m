function [eps, info] = ExtensionEfficiency(pBases, tBases, opts)
%% EXTENSIONEFFICIENCY  Relative rate of extension from a given primer 3' terminus
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 07-11-2026 00:24 IST
%
% PROVENANCE OF THE TABLE: the per-mispair values are digitised from Figure 2 of Huang et al.
% 1992, which plots them as bars on a log axis rather than tabulating them, and which was supplied
% here as a character-recognised copy. They replace the range bounds quoted in that paper's
% abstract.
% ____________________________________________________________________________________________________________
% WHAT THIS FUNCTION DOES, AND WHAT IT REPLACES
%
% Occupancy and extensibility are separate physics. aPCREquilibrium tells you how often a primer
% sits on a site; it says nothing about whether the polymerase will extend from there. A duplex can
% be thermodynamically respectable and catalytically dead if its 3'-terminal base is mispaired,
% because Taq carries no 3'->5' proofreading exonuclease and will not trim the terminus
% (Huang et al. 1992).
%
% MisprimerMismatch.m currently answers this with a binary test, anchored = (l >= n - endTol) with
% endTol = 2, which treats a perfectly matched terminus and a terminus carrying two mismatches in
% the last two positions as the same object. They differ by three to five orders of magnitude in
% extension rate. This function replaces that test, and with it four heuristic constants:
% dG3min, MispEndTol, MispMinRun and MispMargin. All four were standing in for the quantity
% computed here.
%
% The priming rate at a site becomes
%       rate  proportional to  theta(site) * epsilon(terminus)
% and the expected number of extension events per cycle is
%       N = theta * [strand]_ss * V * N_A * epsilon
%
% ____________________________________________________________________________________________________________
% THE TABLE IS NOW SOURCED. Huang MM, Arnheim N, Goodman MF. Nucleic Acids Res 1992;20:4567.
% doi:10.1093/nar/20.17.4567. Retrieved via PubMed on 20-10-2026.
%
% Those authors measured extension efficiency for every configuration of matched and mismatched
% base pair at a Taq primer-template 3' terminus, relative to a correct pair. Two structural results
% from that paper shape how this function is used:
%
%   1. Taq binds matched and mispaired 3' ends with SIMILAR AFFINITY. Their Table II gives the
%      relative binding constant Kp(rel) for all sixteen termini: every value lies between 0.29 and
%      1.00, a spread of about threefold, against the 10^3 to beyond 10^6 spread in extension. The
%      failure to extend a mismatch is therefore a property of the chemistry, not of binding, so
%      epsilon multiplies the catalytic step and must NOT be folded into an equilibrium constant or
%      a binding free energy.
%   2. The ratios are INDEPENDENT OF TEMPERATURE across 45, 55 and 70 C. One table serves every
%      annealing temperature, and no temperature argument is needed or accepted.
%
% WHAT IS IN THE TABLE AND HOW EXACT IT IS. The per-mispair values live in that paper's Figure 2 as
% bars on a log axis, not in a table, so they are read off the figure. The axis runs from 10^0 at the
% frame top to 10^-6 at the baseline over 366 pixels in the 220 dpi rendering used, 61.00 px per
% decade, and each bar top was located by the first row whose ink density over the next eight rows
% exceeds a third of the bar width. The digitised values and the values adopted:
%
%   mispair          digitised   adopted    -RT ln(eps) at 59 C
%   C(primer).T      1.35E-2     2.0E-2     2.58 kcal/mol      <- adopted from the TEXT, see below
%   T.G              8.93E-4     9.0E-4     4.63
%   A.C              5.07E-4     5.0E-4     5.02
%   G.T              2.77E-4     2.8E-4     5.40
%   C.A              1.51E-4     1.5E-4     5.81
%   T.C              1.08E-4     1.1E-4     6.02
%   T.T              3.10E-5     3.1E-5     6.85
%   A.A              2.67E-6     2.7E-6     8.46
%   G.A, A.G, G.G, C.C   no bar  1.0E-6     9.12          <- below the assay's detection limit
%
% C(primer).T is not digitised. The authors state it twice in the body text, as "roughly 1/50 the
% rate for extension of the corresponding correct A.T basepair" and as "f0_ext = 2 x 10^-2", so
% their own number is used in place of the 1.35E-2 the bar reads. The four pairs at the right of the
% figure carry no bar at all; the caption gives f0_ext < 10^-6 and the text estimates "at least a
% million-fold lower efficiency", so 1E-6 is a BOUND standing in for a value, and it is the
% conservative choice because a larger epsilon keeps a site in the threat set.
%
% THREE STATEMENTS IN THE TEXT CHECK THE DIGITISATION, and all three hold:
%   "The values of the transition mispairs, A(primer).C, C.A, G.T, and T.G are in the range of
%    10^-3 to 10^-4"  -> adopted values span 1.5E-4 to 9.0E-4.
%   "the transversion mispairs, T(primer).C, T.T, and A.A were from 5-fold to at least 100-fold
%    less efficient than the transition mismatches"  -> 1 to 8, 5 to 29, and 56 to 333 fold below
%    the four transitions respectively.
%   "C.C, G.G, A.G and G.A were poorly extended, even after an hour incubation"  -> no bar.
%
% C(primer).T is the outlier the authors flag: it extends roughly a hundred times better than the
% other transversions and ten times better than the transitions. A primer ending in C sitting
% opposite a T is therefore the single most dangerous 3'-terminal mismatch there is, and the old
% placeholder table in this file had it at 1E-3, a decade too strict.
%
% ____________________________________________________________________________________________________________
% WHICH TABLE A CALL ACTUALLY GETS, IN ORDER
%
%   1. opts.Table, a 4-by-4 matrix supplied by the caller
%   2. opts.File, a CSV named by the caller
%   3. the session choice held in ExtensionEfficiencyTable, which PrimerSearchParameters sets
%   4. extension_efficiency_rejali2018.csv, found beside this file
%   5. the built-in Huang table, with a warning
%
% STEP 4 EXISTS BECAUSE STEP 3 DOES NOT SURVIVE A FRESH SESSION. The session choice is a persistent
% inside ExtensionEfficiencyTable and only an engine run sets it. A standalone caller in a new
% MATLAB session, PrimerHybridizationMap above all, therefore used to fall straight past it to the
% built-in table and score threats on Huang's four censored bounds while telling the operator that
% Rejali was "the shipped default". It was not a default at all; nothing had set it, and nothing
% had cleared it either. The shipped table is now resolved from this file's own folder, so the
% phrase is true and a map drawn without a search in front of it prices sites the way the search
% would have.
%
% TO ASK FOR THE BUILT-IN TABLE DELIBERATELY, call ExtensionEfficiencyTable('set', 'builtin').
% An empty path there means "resolve the shipped table", which is what 'clear' now restores.
%
% TO OVERRIDE, supply opts.File pointing at a CSV with the columns
%       PrimerBase,OpposingBase,Epsilon
% holding all sixteen combinations, or set opts.Table directly as a 4-by-4 matrix in ACGT order
% with rows indexed by the primer base and columns by the base physically opposing it. Watson-Crick
% entries must be 1. A file may carry further columns named Eps_pos<d>, read as described next; a
% matrix supplied through opts.Table cannot, and scores the terminus alone.
%
% ____________________________________________________________________________________________________________
% INTERIOR POSITIONS ARE NOW MEASURED, NOT ASSUMED
%
% Huang et al. measured the 3' terminus and nothing else, so for years the two positions behind it
% were charged nothing. A mismatch one base in was priced identically to a perfect duplex. Rejali
% NA, Moric E, Wittwer CT, Clin Chem 2018;64(5):801-809, doi:10.1373/clinchem.2017.282285, close
% that gap: their Supplemental Table 2 gives extension rate against mismatch type AND distance from
% the 3' end, at 0, 1, 2, 4, 6, 8 and 10 bases, on two templates, for nine of the twelve mismatch
% types. Relative to the matched rate on the same template, the geometric mean over the templates:
%
%   distance from 3' end      0        1        2        4        6        8       10
%   geometric mean eps      7.8E-4   9.8E-3   9.5E-2   4.4E-1   5.4E-1   8.2E-1   9.6E-1
%
% The identity of the mismatch stops mattering as it moves inside. Regressing log10 eps(d) on
% log10 eps(0) over the nine fully profiled types gives a slope of 0.756 at d = 1, 0.432 at d = 2,
% and 0.073 at d = 4; by d = 6 the slope is 0.006 and the correlation is gone. What a mismatch is
% governs extension at the terminus, and where it sits governs extension behind it.
%
% A FILE SUPPLIES THESE AS COLUMNS Eps_pos1, Eps_pos2 AND SO ON, one per measured distance, with
% the column name carrying the distance. A table with none of them scores the terminus alone and
% behaves exactly as this function did before, which is what the built-in Huang table does.
%
% WHAT IS CHARGED WHEN A REGISTER CARRIES MORE THAN ONE MISMATCH. Rejali measured single
% mismatches, so there is no measured value for a pair of them. Each mismatch present is priced at
% its own distance and the SMALLEST of those efficiencies is returned. That is the single-mismatch
% reading applied conservatively: the register is charged for the worst thing in it and nothing is
% multiplied. Note that the worst is not always the most 3'-proximal, since eps at distance 0
% spans 4.7E-5 to 1.6E-2 while at distance 1 it spans 4.5E-4 to 6.3E-2, so a G.A one base in is
% more inhibitory than a G.T at the terminus.
%
% THE CALLER DECIDES HOW DEEP TO LOOK, by how many positions it passes. MisprimerMismatch passes
% three, so distances 0, 1 and 2 are scored and no further. That is not a limit of this function.
% It is that an anchored register is only guaranteed to reach three bases back from the primer
% terminus; past that the register may have ended, and a mismatch outside the helix is not the
% object Rejali measured. The shipped CSV carries distances 4 through 10 for whoever scores them
% with the register bounds in hand.
%
% ____________________________________________________________________________________________________________
% FRAME CONVENTION, stated once because the suite has lost track of it before
%
% tBases holds the base PHYSICALLY OPPOSING each primer position, the same convention dGMismatch
% takes as its seq2 argument. A Watson-Crick pair means tBases(k) is the complement of pBases(k).
% If you are reading a register out of SeedScan, the opposing base is seqcomplement of the footprint
% at that index, NOT the footprint base itself. Passing the raw footprint inverts every pairing test
% and returns epsilon = 1 for every mismatch, which looks like a working function.
%
% ____________________________________________________________________________________________________________
% INPUTS
%   pBases   char row, primer positions ending at the 3' terminus, 5' to 3', so pBases(end) is
%            position n. One to eleven positions; positions not passed are not examined, which for
%            a matched-terminus register is the same thing as treating them as matched.
%   tBases   char row, the bases opposing those positions, same ordering and length.
%   opts     struct, optional:
%            .Table    4-by-4 numeric, ACGT order, rows primer base, columns opposing base
%            .File     path to a CSV as described above
%
% OUTPUTS
%   eps      relative extension efficiency, 1 for a fully matched duplex, in (0, 1].
%   info     struct:
%            .TerminalPair     the 3'-terminal pair as primer.opposing, always reported
%            .TerminalEps      efficiency of that terminal pair alone
%            .TerminalIsBound  the terminal entry is a detection-limit bound, not a measurement
%            .GoverningPair    the mismatch that set eps, equal to TerminalPair when matched
%            .GoverningPos     its distance from the 3' terminus, 0 at the terminus
%            .nMismatch        mismatches found across the positions passed
%            .PositionsScored  the active table carries interior distances
%            .Citation         the table's own short citation, for a figure caption
%            .Source           where the table came from, a full path for a file
% ____________________________________________________________________________________________________________

% A TABLE SHOULD CARRY ITS OWN CITATION, AND A FIGURE SHOULD NOT GUESS ONE FROM A FILENAME.
% Source names where the table came from and is a full path once a file is in force, which is the
% right answer for a log and the wrong one for a caption. Citation is the short human form, read
% from the file's own Citation column when it has one, so a caller drawing 'epsilon Data: ...' on
% a figure never has to know which file is in force or how to spell its authors.
persistent cachedTable cachedBound cachedPos cachedGrid cachedCite cachedSource placeholderWarned

if (nargin < 3)||isempty(opts), opts = struct(); end
% PenFac and AnteFac are gone. They were multiplicative stand-ins for the interior positions, both
% defaulting to 1 because nothing had measured them. Rejali's distance columns measure them, so a
% caller still passing either is reading from a parameterisation this file no longer has, and is
% told so rather than being silently ignored.
gone = intersect({'PenFac', 'AnteFac'}, fieldnames(opts));
if ~isempty(gone)
    error('ExtensionEfficiency:RetiredOption', ...
        ['opts.%s is retired. Interior positions are no longer multiplicative factors on the ' ...
        'terminal value; they are read from the Eps_pos<d> columns of the efficiency table, ' ...
        'where Rejali et al. 2018 measured them. Remove the field and supply a table that ' ...
        'carries those columns, such as extension_efficiency_rejali2018.csv.'], gone{1});
end
def = struct('Table', [], 'File', '');
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts, fn{k}), opts.(fn{k}) = def.(fn{k}); end
end

pBases = upper(char(pBases));
tBases = upper(char(tBases));
if numel(pBases) ~= numel(tBases)
    error('ExtensionEfficiency:LengthMismatch', ...
        'pBases is %d long and tBases is %d; they must be aligned position by position.', ...
        numel(pBases), numel(tBases));
end
% ELEVEN, because that is how far Rejali measured. The cap is on this function, not on the
% caller's judgement: MisprimerMismatch passes three because an anchored register is guaranteed to
% reach only that far back, and the header says why.
if isempty(pBases)||(numel(pBases) > 11)
    error('ExtensionEfficiency:BadLength', ...
        ['Pass between one and eleven positions, 3'' end last. Received %d. The upper limit is ' ...
        'the depth to which extension efficiency has been measured.'], numel(pBases));
end
if any(~ismember(pBases, 'ACGT'))||any(~ismember(tBases, 'ACGT'))
    error('ExtensionEfficiency:NonACGT', ...
        'Both arguments must be unambiguous DNA. Got primer ''%s'' against ''%s''.', pBases, tBases);
end

%% ---- resolve the table -----------------------------------------------------------------------------
% A caller naming neither a matrix nor a file falls back to THIS SESSION'S CHOICE, which
% PrimerSearchParameters sets once per run through ExtensionEfficiencyTable. That file's header
% says why the choice lives there rather than in an argument, and where it is announced.
useFile = opts.File;
askedBuiltin = false;
if isempty(opts.Table)&&isempty(useFile)
    useFile = ExtensionEfficiencyTable('get');
    if strcmpi(strtrim(useFile), 'builtin')
        askedBuiltin = true; useFile = '';
    elseif isempty(useFile)
        useFile = ShippedTable();
    end
end

if ~isempty(opts.Table)
    tab = opts.Table; bnd = false(4); src = 'caller-supplied matrix';
    ptab = []; pgrid = []; cite = 'caller-supplied matrix';
    ValidateTable(tab);
elseif ~isempty(useFile)
    if ~strcmp(cachedSource, useFile)
        [cachedTable, cachedBound, cachedPos, cachedGrid, cachedCite] = LoadTable(useFile);
        cachedSource = useFile;
    end
    tab = cachedTable; bnd = cachedBound; ptab = cachedPos; pgrid = cachedGrid; src = useFile;
    cite = cachedCite;
    if isempty(cite), cite = Basename(useFile); end
else
    tab = HuangTable(); bnd = HuangBounds(); src = 'Huang et al. 1992, Figure 2 digitised';
    ptab = []; pgrid = []; cite = 'Huang et al. (1992)';
    if isempty(placeholderWarned), placeholderWarned = false; end
    if ~placeholderWarned
        placeholderWarned = true;
        warning('ExtensionEfficiency:DigitisedFromFigure', ...
            ['Extension efficiencies are read off Figure 2 of Huang et al. 1992 (NAR 20:4567), ' ...
            'which plots them as bars on a log axis rather than tabulating them; C.T is taken ' ...
            'from the body text instead, where the authors state it twice. Four entries, G.A, ' ...
            'A.G, G.G and C.C, fall below that assay''s detection limit and are used AS 1E-6, so ' ...
            'they are bounds standing in for values. Everything else is good to roughly 15 per ' ...
            'cent, which is finer than the decade-level resolution the screen needs. See the ' ...
            'header for the digitisation and the three text statements that check it.\n' ...
            'G.G IS THE ONE TO WATCH. Two later studies put it far above that bound rather than ' ...
            'merely unresolved: Rejali et al. 2018 measure 8.8E-4 by stopped-flow, 880 times the ' ...
            'bound and mid-pack among the twelve, and Chen et al. 2015 place it outside the ' ...
            'mismatches showing observable inhibition, and Ayyadevara et al. 2000 get 1.8 per ' ...
            'cent of matched product from G.G in 35-cycle PCR on Taq itself, which no bound of ' ...
            '1E-6 can produce. On this table G.G is charged 9.09 kcal/mol at 58 C against the ' ...
            '4.63 that measurement gives, a difference of 4.47. The shipped ' ...
            'extension_efficiency_rejali2018.csv scores on the measured values instead.\n' ...
            '%s\n' ...
            'Reported once per session.'], ...
            Tern3(askedBuiltin, ...
            ['You are seeing this because ExtensionEfficiencyTable was set to ''builtin'', so ' ...
            'the built-in table was asked for deliberately.'], ...
            ['You are seeing this because extension_efficiency_rejali2018.csv was not found ' ...
            'beside ExtensionEfficiency.m and nothing named a table. Put it back in that folder, ' ...
            'or name one with opts.File, and this call will score on measured values instead.']));
    end
end

%% ---- score every mismatch at its own distance from the terminus ------------------------------------
% comp is index-aligned with bases, so comp(k) is the Watson-Crick partner of bases(k) and a pair
% is matched exactly when tBases(j) equals comp of pBases(j). The old struct form did the same
% thing; this one cannot be indexed with a base that is not ACGT, which is already checked above.
bases = 'ACGT';
comp  = 'TGCA';
nPos  = numel(pBases);

ipT = strfind(bases, pBases(end));
itT = strfind(bases, tBases(end));
epsTerm = tab(ipT, itT);

eps = 1;
govPair = [pBases(end) '.' tBases(end)];
govPos  = 0;
nMis    = 0;
for j = nPos:-1:1
    ip = strfind(bases, pBases(j));
    if strcmp(comp(ip), tBases(j)), continue; end
    it = strfind(bases, tBases(j));
    nMis = nMis + 1;
    e = EpsAtDistance(tab, ptab, pgrid, ip, it, nPos - j);
    if e < eps
        eps = e;
        govPair = [pBases(j) '.' tBases(j)];
        govPos  = nPos - j;
    end
end

info = struct('TerminalIsBound', bnd(ipT, itT), ...
    'TerminalPair', [pBases(end) '.' tBases(end)], 'TerminalEps', epsTerm, ...
    'GoverningPair', govPair, 'GoverningPos', govPos, 'nMismatch', nMis, ...
    'PositionsScored', ~isempty(pgrid), 'Citation', cite, 'Source', src);
end

%% ____________________________________________________________________________________________________________
function s = Basename(p)
%% BASENAME  Last element of a path, split on BOTH separators
% fileparts follows the separator of the platform it runs on, so a Windows path handed to it under
% Linux comes back whole. This does not care.
s = char(p);
if isempty(s), return; end
k = max([strfind(s, '\'), strfind(s, '/'), 0]);
if (k > 0)&&(k < numel(s)), s = s(k + 1:end); end
end

%% ____________________________________________________________________________________________________________
function e = EpsAtDistance(tab, ptab, pgrid, ip, it, d)
%% EPSATDISTANCE  Efficiency for one mismatch sitting d bases from the primer 3' terminus
% d = 0 is the terminus itself and reads the ordinary table. Beyond it, the value comes from the
% distance columns, interpolated LINEARLY IN log10(eps) because the quantity spans four decades and
% a linear interpolation on eps itself would sit an order of magnitude high between measured
% distances. Past the largest measured distance the last measured value is held rather than
% extrapolated, since nothing in the source speaks to what happens further in.
%
% A TABLE WITH NO DISTANCE COLUMNS RETURNS 1 FOR ANY d >= 1. That is not a claim that an interior
% mismatch is free. It is the statement that this table has not measured it, and it reproduces
% exactly what this function did when the interior factors defaulted to 1.
if d <= 0, e = tab(ip, it); return; end
if isempty(pgrid), e = 1; return; end
if d >= pgrid(end)
    e = ptab(ip, it, numel(pgrid));
    return
end
x = [0, pgrid(:)'];
y = log10([tab(ip, it), reshape(ptab(ip, it, :), 1, [])]);
e = 10^interp1(x, y, d, 'linear');
end

%% ____________________________________________________________________________________________________________
function f = ShippedTable()
%% SHIPPEDTABLE  The efficiency table that ships with this suite, or '' when it is not there
% RESOLVED FROM THIS FILE'S OWN FOLDER FIRST, not from the working directory and not from the
% MATLAB path. A bare filename resolved against the working directory identifies a different file
% depending on where MATLAB happened to be, which is the same failure PrimerSearchParameters guards
% against when it expands EffTableFile to a full path before recording it. A copy elsewhere on the
% path is accepted afterwards, for an installation that keeps data apart from code.
f = '';
cand = fullfile(fileparts(mfilename('fullpath')), 'extension_efficiency_rejali2018.csv');
if isfile(cand), f = cand; return; end
w = which('extension_efficiency_rejali2018.csv');
if ~isempty(w)&&isfile(w), f = w; end
end

%% ____________________________________________________________________________________________________________
function s = Tern3(c, a, b)
% Local ternary, so the warning can carry an accurate closing sentence without a branch around the
% whole call.
if c, s = a; else, s = b; end
end

%% ____________________________________________________________________________________________________________
function tab = HuangTable()
% Rows: primer 3' base, in ACGT order. Columns: the base physically opposing it, ACGT order.
% Watson-Crick entries are exactly 1 by definition. Every other entry is digitised from Figure 2 of
% Huang et al. 1992, except C.T, which comes from that paper's body text, and the four entries at
% 1E-6, which are the assay's detection limit standing in for values. See the header.
%
%              opposing A   opposing C   opposing G   opposing T
tab = [ ...
    2.7E-6,      5.0E-4,      1.0E-6,      1.0      ; ...   % primer A:  A.A  A.C  A.G  A.T(WC)
    1.5E-4,      1.0E-6,      1.0,         2.0E-2   ; ...   % primer C:  C.A  C.C  C.G(WC)  C.T
    1.0E-6,      1.0,         1.0E-6,      2.8E-4   ; ...   % primer G:  G.A  G.C(WC)  G.G  G.T
    1.0,         1.1E-4,      9.0E-4,      3.1E-5   ];      % primer T:  T.A(WC)  T.C  T.G  T.T
end

%% ____________________________________________________________________________________________________________
function [tab, bnd, ptab, pgrid, cite] = LoadTable(file)
% CSV with columns PrimerBase, OpposingBase, Epsilon. All sixteen combinations required. An
% optional IsBound column marks entries that are detection-limit bounds rather than measurements;
% absent, nothing is marked.
%
% ANY FURTHER COLUMN NAMED Eps_pos<d> IS READ AS THE EFFICIENCY AT DISTANCE d from the 3' terminus,
% with d a positive integer taken from the column name itself. The columns need not be contiguous
% or evenly spaced, and the shipped Rejali table is not: it carries 1, 2, 4, 6, 8 and 10, which is
% where that paper measured. Every one of the sixteen rows must fill every distance column, the
% four Watson-Crick rows with 1. A file carrying none of these columns loads as before and scores
% the terminus alone. Columns the loader does not recognise, such as the per-template values and
% the provenance strings in the shipped file, are ignored, with one exception: an optional
% Citation column is read and its first non-empty value becomes info.Citation, which is what a
% figure caption should carry. It must hold no comma, since a caption is not worth a quoted-CSV
% parser. Absent, info.Citation falls back to the file's basename.
if ~isfile(file)
    error('ExtensionEfficiency:NoFile', 'Extension-efficiency file not found: %s', file);
end
T = readtable(file, 'TextType', 'string');
need = {'PrimerBase', 'OpposingBase', 'Epsilon'};
for k = 1:numel(need)
    if ~ismember(need{k}, T.Properties.VariableNames)
        error('ExtensionEfficiency:BadFile', ...
            'Column %s is missing from %s. Required columns: %s.', need{k}, file, strjoin(need, ', '));
    end
end
% Distance columns, discovered from the header rather than named here, so that a table measured at
% other distances needs no change to this file.
vn = T.Properties.VariableNames;
pgrid = []; pcol = {};
for k = 1:numel(vn)
    nm = char(vn{k});
    if (numel(nm) > 7)&&strcmp(nm(1:7), 'Eps_pos')
        d = str2double(nm(8:end));
        if isfinite(d)&&(d >= 1)&&(d == fix(d))
            pgrid(end+1) = d;   %#ok<AGROW>
            pcol{end+1}  = nm;  %#ok<AGROW>
        end
    end
end
[pgrid, ord] = sort(pgrid);
pcol = pcol(ord);

bases = 'ACGT';
tab = nan(4, 4);
bnd = false(4, 4);
ptab = nan(4, 4, numel(pgrid));
hasBound = ismember('IsBound', T.Properties.VariableNames);
cite = '';
if ismember('Citation', vn)
    for r = 1:height(T)
        v = T.Citation(r);
        if iscell(v), v = v{1}; end
        v = strtrim(char(v));
        if ~isempty(v), cite = v; break; end
    end
end
for r = 1:height(T)
    ip = strfind(bases, upper(char(T.PrimerBase(r))));
    it = strfind(bases, upper(char(T.OpposingBase(r))));
    if isempty(ip)||isempty(it)
        error('ExtensionEfficiency:BadBase', 'Row %d of %s carries a base outside {A,C,G,T}.', r, file);
    end
    tab(ip, it) = T.Epsilon(r);
    for q = 1:numel(pgrid)
        v = T.(pcol{q})(r);
        if iscell(v), v = str2double(v{1}); end
        ptab(ip, it, q) = v;
    end
    if hasBound
        v = T.IsBound(r);
        if isnumeric(v)||islogical(v)
            bnd(ip, it) = logical(v);
        else
            bnd(ip, it) = ismember(lower(strtrim(char(v))), {'1', 'true', 'y', 'yes'});
        end
    end
end
ValidateTable(tab);
ValidateDistances(ptab, pgrid, file);
end

%% ____________________________________________________________________________________________________________
function ValidateDistances(ptab, pgrid, file)
%% VALIDATEDISTANCES  The distance columns must be as complete and as bounded as the terminal one
% MONOTONICITY IS NOT CHECKED AND MUST NOT BE. Extension efficiency rises with distance overall,
% but Rejali's template 1 has C.A, C.T and C.C each slightly lower at 6 bases than at 4, every one
% of those pairs inside a single standard deviation of the other. Forcing the sequence upward would
% reject a faithful transcription of the measurement to satisfy an expectation the measurement does
% not quite meet.
if isempty(pgrid), return; end
for q = 1:numel(pgrid)
    v = ptab(:, :, q);
    if any(~isfinite(v(:)))
        error('ExtensionEfficiency:IncompleteDistance', ...
            ['Column Eps_pos%d of %s has a missing entry. Every distance column must carry all ' ...
            'sixteen primer-base by opposing-base combinations, the Watson-Crick four as 1.'], ...
            pgrid(q), file);
    end
    if any(v(:) <= 0)||any(v(:) > 1)
        error('ExtensionEfficiency:DistanceOutOfRange', ...
            'Column Eps_pos%d of %s has an efficiency outside (0, 1].', pgrid(q), file);
    end
    wc = [v(1,4), v(2,3), v(3,2), v(4,1)];
    if any(abs(wc - 1) > 1E-12)
        error('ExtensionEfficiency:DistanceWatsonCrickNotUnity', ...
            ['The four Watson-Crick cells of column Eps_pos%d of %s must be exactly 1. A matched ' ...
            'pair does not become less extendable for sitting inside the duplex, and a value off ' ...
            'unity there would charge every register with a matched base at that distance. Got ' ...
            '%.4g %.4g %.4g %.4g for A.T, C.G, G.C and T.A.'], ...
            pgrid(q), file, wc(1), wc(2), wc(3), wc(4));
    end
end
end

%% ____________________________________________________________________________________________________________
function bnd = HuangBounds()
%% HUANGBOUNDS  Which entries of the built-in table are detection-limit bounds, not measurements
% Huang et al. report G.A, A.G, G.G and C.C as "less than 1E-6" rather than as values, so the 1E-6
% stored for them is an upper bound standing in for an unmeasured number. Everything else in that
% table is a measurement. Marking them lets a caller say so when a rejection turns on one.
%
% THE FOUR ARE NOT EQUALLY SUSPECT. Rejali et al. 2018, measuring extension rates by stopped-flow
% on KlenTaq at 65 C, resolve all twelve mismatches with nothing left at a detection limit. As
% epsilon relative to a matched terminus, geometric mean over their three templates:
%       G.A  4.7E-5   rank 1 of 12, least extendable. Huang's ordering holds.
%       A.G  6.7E-5   rank 2 of 12. Huang's ordering holds.
%       C.C  2.0E-4   rank 3 of 12, with A.A alongside it.
%       G.G  8.8E-4   rank 5 of 12. MID-PACK, above T.T, and 880 times the stored bound.
% Three of the four are therefore too strict by one to two decades while keeping their rank, and
% G.G is wrong in rank as well as magnitude. Chen et al. 2015 reach the same conclusion about G.G
% independently, placing it outside the mismatches that show observable inhibition.
bnd = false(4, 4);
bases = 'ACGT';
pairs = {'G', 'A'; 'A', 'G'; 'G', 'G'; 'C', 'C'};
for k = 1:size(pairs, 1)
    bnd(strfind(bases, pairs{k, 1}), strfind(bases, pairs{k, 2})) = true;
end
end

%% ____________________________________________________________________________________________________________
function ValidateTable(tab)
if ~isequal(size(tab), [4 4])
    error('ExtensionEfficiency:BadTableSize', 'The table must be 4-by-4 in ACGT order.');
end
if any(~isfinite(tab(:)))
    error('ExtensionEfficiency:IncompleteTable', ...
        'The table has a missing entry. All sixteen primer-base by opposing-base combinations are required.');
end
if any(tab(:) <= 0)||any(tab(:) > 1)
    error('ExtensionEfficiency:OutOfRange', ...
        'Every efficiency must lie in (0, 1], with 1 for a matched terminus.');
end
% Watson-Crick entries, ACGT order: A-T is (1,4), C-G is (2,3), G-C is (3,2), T-A is (4,1).
wc = [tab(1,4), tab(2,3), tab(3,2), tab(4,1)];
if any(abs(wc - 1) > 1E-12)
    error('ExtensionEfficiency:WatsonCrickNotUnity', ...
        ['The four Watson-Crick entries must be exactly 1, since epsilon is defined relative to a ' ...
        'matched terminus. Got %.4g %.4g %.4g %.4g for A.T, C.G, G.C and T.A. A table with these ' ...
        'off unity rescales every site by a constant and changes nothing except the absolute ' ...
        'event counts, which is almost certainly not what was meant.'], wc(1), wc(2), wc(3), wc(4));
end

% ORIENTATION. THE TEST ABOVE CANNOT SEE A TRANSPOSED TABLE.
% The four Watson-Crick cells sit on the anti-diagonal, and the anti-diagonal maps onto itself
% under transposition. A table entered with rows and columns swapped therefore passes the size,
% range and Watson-Crick tests untouched, and every subsequent lookup reads epsilon out of the
% wrong cell. The failure is silent and systematic: it exchanges each mismatch with its mirror, so
% a primer ending in C opposite a T would be scored with the efficiency measured for a primer
% ending in T opposite a C.
%
% The table's own asymmetry settles the orientation. Huang et al. measure primer C opposite T at
% 2.0E-2 and primer T opposite C at 1.1E-4, a factor of 182, which is the largest asymmetry in the
% set and the one the header singles out as the most dangerous 3' terminus there is. Transposition
% reverses that inequality, so one comparison decides it.
ct = tab(2, 4);     % primer C, opposing T
tc = tab(4, 2);     % primer T, opposing C
if ct < tc
    error('ExtensionEfficiency:TableTransposed', ...
        ['The table appears to be TRANSPOSED. Primer C opposite T reads %.4g and primer T ' ...
        'opposite C reads %.4g, but the measured set has the first larger than the second by ' ...
        'about two orders of magnitude. Rows must be indexed by the PRIMER base and columns by ' ...
        'the base physically OPPOSING it, in ACGT order. The Watson-Crick test cannot catch ' ...
        'this, because the four matched cells lie on the anti-diagonal and transposing the ' ...
        'table leaves them in place.'], ct, tc);
elseif ct < 2*tc
    warning('ExtensionEfficiency:OrientationUnverifiable', ...
        ['Primer C opposite T reads %.4g against primer T opposite C at %.4g, a ratio of %.2f. ' ...
        'The orientation test needs these two to differ to tell a table from its transpose, and ' ...
        'at this ratio it cannot. Check by hand that rows are the primer base and columns the ' ...
        'opposing base before trusting any efficiency from this table.'], ct, tc, ct/tc);
end
end
