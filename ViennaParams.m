function P = ViennaParams(parFile, opts)
%% VIENNAPARAMS  Load dna_mathews2004.par into the tables the rebuilt dGMismatch.m reads
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 12-10-2026 13:40 UTC
% ____________________________________________________________________________________________________________
% WHAT THIS FUNCTION SUPPLIES, AND WHAT IT DELIBERATELY DOES NOT
%
% The suite keeps the SantaLucia and Hicks engine, because that is where the Owczarzy monovalent and
% divalent salt treatment lives and this file has no equivalent for it. What SantaLucia does not
% publish is taken from here:
%
%   TERMINAL MISMATCHES   from mismatch_exterior. This replaces the constant term_mismatch_DG =
%                         -1.23 that currently appears in dGMismatch.m, LooperMismatch.m and
%                         InternalLoopBridge.m with three comments telling you to keep them in sync.
%                         A constant applied to every candidate shifts the absolute scale and
%                         discriminates nothing; the table here runs from 0 to about -1.0 kcal/mol
%                         and depends on both mismatch bases and the closing pair.
%   DANGLING ENDS         from dangle5 and dangle3. SeedScan currently trims every register to its
%                         outermost identity and scores it as a blunt oligomer, which on a 7249 nt
%                         circle is never the physical situation.
%   SMALL INTERNAL LOOPS  from int11, int21 and int22, replacing the generic Table 4 length term for
%                         runs of one and two non-Watson-Crick positions.
%   LOOP LENGTH TABLES    from hairpin, bulge and internal, as a cross-check against the Table 4
%                         values hard-coded in dGLoop.m rather than as a replacement for them.
%
% NOT taken from here: the stack table, which comes from SantaLucia, and ML_params and NINIO, which
% describe multiloops that cannot occur in a primer-template duplex.
%
% ____________________________________________________________________________________________________________
% THE ONE BOOKKEEPING RULE THAT MATTERS
%
% Vienna and SantaLucia partition the same physical energy differently. In the Turner decomposition
% a helix carries a separate duplex initiation plus a terminal AU or GU penalty, and dangles and
% terminal mismatches are additive bonuses on top. In the SantaLucia unified model the bimolecular
% initiation and the terminal A.T penalty are already folded into rows 17 and 18 of the
% nearest-neighbour table that dGPrim.m and dGMismatch.m use.
%
% So: dangles and terminal mismatches from this file ARE additive increments on a helix end in both
% frameworks and transplant cleanly. The Vienna terminal AU/GU penalty MUST NOT be applied on top of
% a SantaLucia duplex, because row 18 already carries it. That penalty lives in the Misc block and
% is loaded below for completeness, flagged, and never used.
%
% ____________________________________________________________________________________________________________
% WHY EVERY TABLE IS BUILT FROM ITS ROW LABELS
%
% The file labels its own rows: "/* CG,GC,N */" on a line of int11 names the first pair, the second
% pair and the first mismatch base, and the five values across that line run over the second
% mismatch base. Reshaping by position instead would require guessing the nesting order, and a
% guess that is wrong by one level produces a table of entirely plausible numbers in the wrong
% slots, which nothing downstream would catch.
%
% So each labelled block is filled BY NAME and then checked for completeness. Positional reshape is
% used only for the unlabelled blocks, which are the three loop-length tables and the three scalar
% parameter rows.
%
% ____________________________________________________________________________________________________________
% CONVENTIONS, READ OFF THE FILE RATHER THAN ASSUMED
%
%   Pair index       CG=1 GC=2 GT=3 TG=4 AT=5 TA=6 NN=7.  DNA, so T and not U.
%   Base index       N=1 A=2 C=3 G=4 T=5, from the column header "/*  N  A  C  G  T  */" above
%                    dangle3. INDEX ONE IS A PLACEHOLDER, NOT ADENINE. A parser that assumes
%                    A,C,G,T shifts every dangle and every terminal mismatch by one base.
%   int22 ONLY       six pair types, NN absent, and four bases A=1 C=2 G=3 T=4 with no placeholder.
%                    The two conventions coexist; P.BaseIndex and P.BaseIndex4 keep them apart.
%   Units            the file stores 10 cal/mol, so every value is divided by 100 on load and every
%                    table this function returns is in kcal/mol.
%   Temperature      dG(T) = dH - (dH - dG37)*T/310.15, T in Kelvin, the Turner relation. Use
%                    P.AtTemp.
%
% ____________________________________________________________________________________________________________
% KNOWN DEFECTS IN THE SHIPPED FILE, both confirmed against ViennaRNA 2.7.2
%
%   Line 3 reads "/* DNA stacking enthalpies not symmetric in SOURCE!!! /*", opening a comment and
%   then opening it again instead of closing it. Any parser that searches across newlines for the
%   next "*/" will blank the "# stack" header on line 5 and lose the whole stack table. Comments are
%   therefore stripped ONE LINE AT A TIME here. Do not replace StripLineComments with a dotall
%   pattern, however much tidier it looks.
%
%   The warning in that comment is accurate. The stack ENERGY matrix is symmetric; the stack
%   ENTHALPY matrix is not, by 0.3 kcal/mol at AT against TG. We do not use the stack table, but the
%   enthalpy companions of the dangle and mismatch blocks share its provenance, so symmetry is
%   audited on load and reported.
%
%   Hexaloops, Tetraloops and Triloops ship EMPTY. The triloop and tetraloop bonus gap flagged in
%   LooperMismatch.m cannot be closed from this file. GNRA-class tetraloops stay under-detected by
%   up to about 1.6 kcal/mol; source those from the SantaLucia and Hicks supplementary material or
%   the Peyret thesis if you need them.
%
% ____________________________________________________________________________________________________________
% INPUT
%   parFile  optional path. Defaults to the location DNAAccessibility.m resolves.
%   opts     .Verbose  print the load summary and the three hand-checkable values, default true
%            .Force    bypass the persistent cache, default false
%
% OUTPUT  P, struct
%   .stack .dangle5 .dangle3 .mismatch_exterior .mismatch_hairpin .mismatch_internal
%   .mismatch_internal_1n .mismatch_internal_23 .mismatch_multi .int11 .int21 .int22
%   .hairpin .bulge .internal        each with .dG37 and .dH in kcal/mol
%   .PairIndex .BaseIndex .BaseIndex4    lookup structs
%   .AtTemp                              @(dG37, dH, TdegC) -> dG at that temperature
%   .Dangle3 .Dangle5 .TerminalMismatch .Int11    convenience accessors taking letters
%   .Audit                               symmetry report
%   .Meta                                path, bytes, sha, block line numbers
%
% REFERENCES
%   Mathews DH et al. PNAS 2004;101:7287.           doi:10.1073/pnas.0401799101
%   SantaLucia J, Hicks D. Annu Rev Biophys Biomol Struct 2004;33:415.
%                                                   doi:10.1146/annurev.biophys.32.110601.141800
%   Lorenz R et al. Algorithms Mol Biol 2011;6:26.  doi:10.1186/1748-7188-6-26   (ViennaRNA 2.0)
%   I have not been able to open these; verify before citing.
% ____________________________________________________________________________________________________________

persistent CACHE

if (nargin < 1)||isempty(parFile)
    parFile = 'C:\Program Files (x86)\ViennaRNA Package\Misc\dna_mathews2004.par';
end
if (nargin < 2)||isempty(opts), opts = struct(); end
if ~isfield(opts, 'Verbose')||isempty(opts.Verbose), opts.Verbose = true; end
if ~isfield(opts, 'Force')||isempty(opts.Force), opts.Force = false; end

if ~isfile(parFile)
    error('ViennaParams:NoFile', ...
        ['Parameter file not found: %s\nRun FindViennaParams to locate it. Check you are pointing ' ...
        'at the 2004 file and not dna_mathews1999.par, which ships alongside and uses a different ' ...
        'parameterisation.'], parFile);
end
d = dir(parFile);
key = sprintf('%s|%d|%.10f', parFile, d.bytes, d.datenum);
if ~opts.Force && ~isempty(CACHE) && strcmp(CACHE.Meta.Key, key)
    P = CACHE; return
end

%% ---- read, strip per line, index blocks -------------------------------------------------------------
rawlines = regexp(fileread(parFile), '\r\n|\n|\r', 'split');
nL = numel(rawlines);
clean = cell(1, nL);
for i = 1:nL
    clean{i} = StripLineComments(rawlines{i});
end

blockLine = struct();
order = strings(0,1);
for i = 1:nL
    s = strtrim(clean{i});
    if isempty(s), continue; end
    tok = regexp(s, '^#(?!#)\s*(\S+)', 'tokens', 'once');
    if isempty(tok), continue; end
    nm = matlab.lang.makeValidName(tok{1});
    blockLine.(nm) = i;
    order(end+1,1) = string(tok{1}); %#ok<AGROW>
end

%% ---- alphabets --------------------------------------------------------------------------------------
pairs  = {'CG', 'GC', 'GT', 'TG', 'AT', 'TA', 'NN'};
bases5 = {'N', 'A', 'C', 'G', 'T'};
bases4 = {'A', 'C', 'G', 'T'};
PairIndex = BuildIndex(pairs);
BaseIndex = BuildIndex(bases5);
BaseIndex4 = BuildIndex(bases4);

P = struct();
P.PairIndex = PairIndex; P.BaseIndex = BaseIndex; P.BaseIndex4 = BaseIndex4;
P.PairNames = pairs; P.BaseNames = bases5; P.BaseNames4 = bases4;

%% ---- labelled blocks ----------------------------------------------------------------------------------
% spec: {fieldname, filename, dims, kindsOfLabelFields, alphabetForBases}
% The trailing dimension is filled by the values across each line; the leading dimensions come from
% the row label.
lab = { ...
    'stack',                'stack',                [7 7],         {'pair'},                      BaseIndex; ...
    'dangle5',              'dangle5',              [7 5],         {'pair'},                      BaseIndex; ...
    'dangle3',              'dangle3',              [7 5],         {'pair'},                      BaseIndex; ...
    'mismatch_hairpin',     'mismatch_hairpin',     [7 5 5],       {'pair', 'base'},              BaseIndex; ...
    'mismatch_internal',    'mismatch_internal',    [7 5 5],       {'pair', 'base'},              BaseIndex; ...
    'mismatch_internal_1n', 'mismatch_internal_1n', [7 5 5],       {'pair', 'base'},              BaseIndex; ...
    'mismatch_internal_23', 'mismatch_internal_23', [7 5 5],       {'pair', 'base'},              BaseIndex; ...
    'mismatch_multi',       'mismatch_multi',       [7 5 5],       {'pair', 'base'},              BaseIndex; ...
    'mismatch_exterior',    'mismatch_exterior',    [7 5 5],       {'pair', 'base'},              BaseIndex; ...
    'int11',                'int11',                [7 7 5 5],     {'pair', 'pair', 'base'},      BaseIndex; ...
    'int21',                'int21',                [7 7 5 5 5],   {'pair', 'pair', 'base', 'base'}, BaseIndex; ...
    'int22',                'int22',                [6 6 4 4 4 4], {'pair', 'pair', 'base', 'base', 'base'}, BaseIndex4};

for r = 1:size(lab, 1)
    fld = lab{r,1}; nm = lab{r,2}; dims = lab{r,3}; kinds = lab{r,4}; bidx = lab{r,5};
    P.(fld).dG37 = ReadLabelled(rawlines, clean, blockLine, nm,                dims, kinds, PairIndex, bidx)/100;
    P.(fld).dH   = ReadLabelled(rawlines, clean, blockLine, [nm '_enthalpies'], dims, kinds, PairIndex, bidx)/100;
end

%% ---- unlabelled blocks ---------------------------------------------------------------------------------
% Loop-length tables run from size 0 to 30, so index k holds size k-1.
for nm = {'hairpin', 'bulge', 'internal'}
    f = nm{1};
    P.(f).dG37 = ReadPositional(clean, blockLine, f,                31)/100;
    P.(f).dH   = ReadPositional(clean, blockLine, [f '_enthalpies'], 31)/100;
end
P.ML_params = ReadPositional(clean, blockLine, 'ML_params', 6)/100;
P.NINIO     = ReadPositional(clean, blockLine, 'NINIO',     3)/100;
P.Misc      = ReadPositional(clean, blockLine, 'Misc',      4)/100;

P.EmptyBlocks = {'Hexaloops', 'Tetraloops', 'Triloops'};

%% ---- temperature scaling and accessors ------------------------------------------------------------------
P.AtTemp = @(dG37, dH, TdegC) dH - (dH - dG37).*((TdegC + 273.15)/310.15);

d3 = P.dangle3; d5 = P.dangle5; mx = P.mismatch_exterior; i11 = P.int11;
At = P.AtTemp;
P.Dangle3 = @(pair, base, TdegC) ...
    At(d3.dG37(PairIndex.(pair), BaseIndex.(base)), d3.dH(PairIndex.(pair), BaseIndex.(base)), TdegC);
P.Dangle5 = @(pair, base, TdegC) ...
    At(d5.dG37(PairIndex.(pair), BaseIndex.(base)), d5.dH(PairIndex.(pair), BaseIndex.(base)), TdegC);
P.TerminalMismatch = @(pair, b1, b2, TdegC) ...
    At(mx.dG37(PairIndex.(pair), BaseIndex.(b1), BaseIndex.(b2)), ...
       mx.dH(  PairIndex.(pair), BaseIndex.(b1), BaseIndex.(b2)), TdegC);
P.Int11 = @(p1, p2, b1, b2, TdegC) ...
    At(i11.dG37(PairIndex.(p1), PairIndex.(p2), BaseIndex.(b1), BaseIndex.(b2)), ...
       i11.dH(  PairIndex.(p1), PairIndex.(p2), BaseIndex.(b1), BaseIndex.(b2)), TdegC);

%% ---- audits -----------------------------------------------------------------------------------------------
A = struct();
A.stack_dG37 = max(max(abs(P.stack.dG37 - P.stack.dG37')));
A.stack_dH   = max(max(abs(P.stack.dH   - P.stack.dH')));
% int11 is invariant under rotating the loop end for end: swapping the two closing pairs and the two
% mismatch bases together must leave the value unchanged.
e = 0;
for t1 = 1:7
    for t2 = 1:7
        for b1 = 1:5
            for b2 = 1:5
                e = max(e, abs(P.int11.dG37(t1,t2,b1,b2) - P.int11.dG37(t2,t1,b2,b1)));
            end
        end
    end
end
A.int11_dG37 = e;
P.Audit = A;

if A.stack_dH > 1E-9
    warning('ViennaParams:AsymmetricEnthalpies', ...
        ['The stack enthalpy matrix is asymmetric by %.2f kcal/mol, which the shipped file warns ' ...
        'about on line 3. Stacks are taken from SantaLucia here so this does not bite directly, ' ...
        'but the dangle and mismatch enthalpies share its provenance. Treat any conclusion that ' ...
        'turns on %.2f kcal/mol of enthalpy as unsupported.'], A.stack_dH, A.stack_dH);
end
if A.int11_dG37 > 1E-9
    warning('ViennaParams:Int11NotSymmetric', ...
        ['int11 is not invariant under swapping the two closing pairs together with the two ' ...
        'mismatch bases, maximum deviation %.2f kcal/mol. Either the file carries a genuine ' ...
        'asymmetry or the index order assumed from the row labels is wrong. Check a labelled row ' ...
        'by hand before using int11.'], A.int11_dG37);
end

%% ---- metadata and cache ------------------------------------------------------------------------------------
P.Meta = struct('Path', parFile, 'Bytes', d.bytes, 'Key', key, ...
    'Modified', datestr(d.datenum, 'yyyy-mm-dd HH:MM:SS'), ...
    'BlockLines', blockLine, 'BlockOrder', order, 'Units', 'kcal/mol');
CACHE = P;

%% ---- load summary ---------------------------------------------------------------------------------------
if opts.Verbose
    fprintf('\nViennaParams loaded %s\n', parFile);
    fprintf('  %d labelled tables, %d loop-length tables, all values kcal/mol.\n', size(lab,1), 3);
    fprintf('  Symmetry: stack dG37 %.3g, stack dH %.3g, int11 dG37 %.3g  (zero is the expectation)\n', ...
        A.stack_dG37, A.stack_dH, A.int11_dG37);
    fprintf('  Empty in this file: %s\n', strjoin(P.EmptyBlocks, ', '));
    fprintf('  Hand-checkable against the raw dump:\n');
    fprintf('    dangle3(CG, T)           = %+6.2f   expect -0.30  (file line 491, fifth value -30)\n', ...
        P.dangle3.dG37(PairIndex.CG, BaseIndex.T));
    fprintf('    mismatch_exterior(CG,A,A)= %+6.2f   expect -1.00  (file line 397, second value -100)\n', ...
        P.mismatch_exterior.dG37(PairIndex.CG, BaseIndex.A, BaseIndex.A));
    fprintf('    int11(CG,CG,A,A)         = %+6.2f   expect +0.90  (file line 511, second value 90)\n', ...
        P.int11.dG37(PairIndex.CG, PairIndex.CG, BaseIndex.A, BaseIndex.A));
    fprintf('    stack(CG,GC)             = %+6.2f   expect -1.80  (file line 7, second value -180)\n\n', ...
        P.stack.dG37(PairIndex.CG, PairIndex.GC));
end
end

%% ____________________________________________________________________________________________________________
function s = StripLineComments(s)
% Remove C-style comments from a SINGLE line. Well-formed pairs first, then an unterminated opener
% truncates the rest of its own line, then an orphan closer is dropped. Nothing crosses a newline,
% which is what makes the malformed comment on line 3 of the shipped file harmless.
if isempty(s), return; end
s = regexprep(s, '/\*.*?\*/', ' ');
ix = strfind(s, '/*');
if ~isempty(ix), s = s(1:ix(1) - 1); end
s = strrep(s, '*/', ' ');
end

%% ____________________________________________________________________________________________________________
function lbl = LineLabel(rawline)
% Last /* ... */ on a raw line, which is where this file puts its row labels.
lbl = '';
m = regexp(rawline, '/\*(.*?)\*/', 'tokens');
if isempty(m), return; end
lbl = strtrim(m{end}{1});
end

%% ____________________________________________________________________________________________________________
function idx = BuildIndex(names)
idx = struct();
for k = 1:numel(names)
    idx.(names{k}) = k;
end
end

%% ____________________________________________________________________________________________________________
function v = TokenValues(s)
% Numeric tokens of a stripped line, with INF mapped to +Inf.
t = regexp(s, '[-+]?\d+|INF', 'match');
v = zeros(1, numel(t));
for k = 1:numel(t)
    if strcmp(t{k}, 'INF')
        v(k) = Inf;
    else
        v(k) = str2double(t{k});
    end
end
end

%% ____________________________________________________________________________________________________________
function M = ReadLabelled(rawlines, clean, blockLine, name, dims, kinds, pairIdx, baseIdx)
% Fill an array from a labelled block. Each data line carries a label naming its leading indices and
% a row of values spanning the trailing dimension.
key = matlab.lang.makeValidName(name);
if ~isfield(blockLine, key)
    error('ViennaParams:MissingBlock', ...
        'Block "%s" is not present in the file. Run FindViennaParams to list what is.', name);
end
i0 = blockLine.(key);
M = nan(dims);
filled = false(dims);
nTrail = dims(end);
nLead = numel(dims) - 1;
if numel(kinds) ~= nLead
    error('ViennaParams:SpecMismatch', ...
        'Block "%s": %d label fields declared for %d leading dimensions.', name, numel(kinds), nLead);
end
posRow = 0;

for i = (i0 + 1):numel(clean)
    s = strtrim(clean{i});
    if ~isempty(s) && (s(1) == '#'), break; end
    v = TokenValues(s);
    if isempty(v), continue; end
    if numel(v) ~= nTrail
        error('ViennaParams:RowWidth', ...
            'Block "%s" line %d holds %d values; %d were expected across the trailing dimension.', ...
            name, i, numel(v), nTrail);
    end

    lbl = LineLabel(rawlines{i});
    sub = cell(1, nLead);
    if ~isempty(lbl)
        f = strtrim(strsplit(lbl, ','));
        if numel(f) ~= nLead
            error('ViennaParams:LabelFields', ...
                ['Block "%s" line %d is labelled "%s", which names %d indices where %d were ' ...
                'expected. The nesting order assumed for this block is wrong.'], ...
                name, i, lbl, numel(f), nLead);
        end
        for k = 1:nLead
            if strcmp(kinds{k}, 'pair')
                if ~isfield(pairIdx, f{k})
                    error('ViennaParams:BadPairLabel', ...
                        'Block "%s" line %d: "%s" is not a known pair type.', name, i, f{k});
                end
                sub{k} = pairIdx.(f{k});
            else
                if ~isfield(baseIdx, f{k})
                    error('ViennaParams:BadBaseLabel', ...
                        ['Block "%s" line %d: "%s" is not a known base for this block. Note that ' ...
                        'int22 uses A,C,G,T with no placeholder while every other block uses ' ...
                        'N,A,C,G,T.'], name, i, f{k});
                end
                sub{k} = baseIdx.(f{k});
            end
        end
    else
        % Unlabelled row inside a labelled block: fall back to position, which must then be exact.
        rem = posRow;
        for k = nLead:-1:1
            sub{k} = mod(rem, dims(k)) + 1;
            rem = floor(rem/dims(k));
        end
    end
    posRow = posRow + 1;

    for t = 1:nTrail
        M(sub{:}, t) = v(t);
        filled(sub{:}, t) = true;
    end
end

if ~all(filled(:))
    error('ViennaParams:Incomplete', ...
        ['Block "%s" left %d of %d entries unfilled. Either the block is truncated or the label ' ...
        'convention differs from the one declared for it.'], name, nnz(~filled), numel(filled));
end
end

%% ____________________________________________________________________________________________________________
function v = ReadPositional(clean, blockLine, name, nWant)
% Unlabelled block: read tokens forward until nWant have been collected.
key = matlab.lang.makeValidName(name);
if ~isfield(blockLine, key)
    error('ViennaParams:MissingBlock', 'Block "%s" is not present in the file.', name);
end
i0 = blockLine.(key);
v = [];
for i = (i0 + 1):numel(clean)
    s = strtrim(clean{i});
    if ~isempty(s) && (s(1) == '#'), break; end
    v = [v, TokenValues(s)]; %#ok<AGROW>
    if numel(v) >= nWant, break; end
end
if numel(v) < nWant
    error('ViennaParams:ShortBlock', ...
        'Block "%s" yielded %d values where %d were expected.', name, numel(v), nWant);
end
v = v(1:nWant);
end
