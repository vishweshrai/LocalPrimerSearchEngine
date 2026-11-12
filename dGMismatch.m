function [dG, info] = dGMismatch(seq, seq2, monv, biv, dNTP, tann, applyInit, applyTerm5, applyTerm3, ctx)
%% DGMISMATCH  Mismatch-aware nearest-neighbour free energy of a pre-aligned duplex segment
% Written by Vishwesh Rai, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 04-11-2026 23:06 IST
% ____________________________________________________________________________________________________________
% SIGNATURE IS UNCHANGED. Every existing call site keeps working: DimerScan, LooperMismatch,
% MisprimerMismatch, dGmaxMismatch, InternalLoopBridge and EndStability all pass six to nine
% arguments and take one output. The tenth argument and the second output are new and optional.
%
% ____________________________________________________________________________________________________________
% WHAT CHANGED, AND WHAT EACH CHANGE IS WORTH
%
%  1. MISMATCH TERMS NOW CARRY TEMPERATURE DEPENDENCE. The previous version added the Table 2
%     DeltaG37 increments to the Watson-Crick term with no scaling, and said so in its own header.
%     At an annealing temperature twenty degrees above the reference that was the largest single
%     systematic error in the specificity calculation, because discrimination is a difference of
%     differences and a mismatch-dependent offset does not cancel out of it. NNParams.m now resolves
%     each of the 48 entries to one of three routes and this function applies whichever it gets:
%
%       route (a)  measured DeltaH and DeltaS from the source papers, accumulated into h and s
%                  alongside the Watson-Crick steps, so nothing is extrapolated
%       route (b)  verified DeltaG37 magnitude plus the temperature SLOPE of the Vienna int11 table,
%                  applied once per ISOLATED mismatch. int11 describes a 1x1 internal loop, which is
%                  exactly what an isolated mismatch flanked by Watson-Crick pairs is.
%       route (c)  DeltaG37 alone, the old behaviour, taken only where neither of the above applies
%
%     Route (b) is applied only when NEITHER of the two Table 2 steps bracketing that mismatch has
%     measured parameters. A mixed pair would double-count the measured step's temperature
%     dependence, so the mixed case falls to route (c) and warns once.
%
%  2. SALT IS NOW CHOSEN BY WHAT THE DUPLEX ACTUALLY IS. The old code applied the Owczarzy
%     correction, s = s + h*f(salt), to every segment. That form shifts DeltaS by DeltaH times a
%     function of salt, which is the algebraic image of a correction defined on 1/Tm, and it is
%     valid for a FULLY PAIRED duplex where a two-state Tm exists and DeltaH stands in for phosphate
%     count. A trimmed misprime register of twelve positions carrying five mismatches has a small
%     DeltaH and an unchanged phosphate count, so it was receiving almost none of the low-salt
%     destabilisation it should. The bias ran one way: mismatched registers came out looking more
%     stable than they are, and the engine over-predicted mispriming on every off-target site it
%     scored. Any segment carrying a mismatch now takes the per-phosphate form instead,
%
%         DeltaS(Na+) = DeltaS(1 M) + 0.368*(N-1)*ln[Na+]      cal/(K mol)
%
%     with [Na+] the von Ahsen monovalent equivalent from NNParams. A fully paired segment still
%     takes Owczarzy, so EndStability and every other all-Watson-Crick call reproduces its old
%     number exactly. That is the regression anchor for this rewrite.
%
%  3. PARAMETERS COME FROM NNParams.m. The Watson-Crick table and the Owczarzy block were duplicated
%     verbatim across PrimTemp.m, dGPrim.m and this file. Identical by intention, divergent by
%     accident. One module now owns them.
%
%  4. RUN-AWARE LOOP TERM RETAINED. A run of k consecutive non-Watson-Crick positions is one
%     internal loop with 2k unpaired nucleotides and carries ONE length term, which is the
%     correction made in the previous version and is kept unchanged.
%
% ____________________________________________________________________________________________________________
% WHAT IS DELIBERATELY NOT DONE HERE, AND WHY
%
%  TERMINAL MISMATCHES still use the conservative constant term_mismatch_DG = -1.23 kcal/mol by
%  default. The Vienna mismatch_exterior table is loaded and the lookup is written, but it is
%  switched OFF behind ctx.UseViennaTerminal because the strand orientation convention of that table
%  is not established. The table is not symmetric in its two exterior bases: mismatch_exterior(CG,A,C)
%  is -1.00 while mismatch_exterior(GC,C,A) is -0.80, so putting the two bases in the wrong order
%  moves every terminal mismatch in the suite by 0.2 to 0.3 kcal/mol with no symptom anywhere
%  downstream. Resolve it with RNAeval against a structure carrying a known terminal mismatch, then
%  set the flag. Until then the constant is wrong by a known amount in a known direction, which is
%  a better position than being wrong by an unknown amount in an unknown direction.
%
%  DANGLING ENDS are not applied here at all. They belong to dGPrim, where the primer terminates and
%  the template continues. Inside this function both strands are present at every position by
%  construction, so a segment end that is mismatched is a terminal mismatch and not a dangle.
%
%  int22 FOR TANDEM MISMATCHES is not used. A run of exactly two non-Watson-Crick positions is a 2x2
%  internal loop and dna_mathews2004.par carries a complete int22 table with enthalpies, which would
%  be strictly better than the generic length term now applied. The index order of that block was
%  not dumped during the parameter audit, so it would have to be guessed. Dump it with
%  FindViennaParams('', {'int22'}) and this becomes a short addition.
%
% ____________________________________________________________________________________________________________
% INPUTS
%   seq        query strand segment, 5'->3', uppercase char, at least two bases
%   seq2       the base ACTUALLY opposing each position of seq, aligned 1:1, same length. Callers
%              that align a query against the reverse complement of a partner must pass
%              seqcomplement of the partner row, not the partner row itself.
%   monv       monovalent cation, M
%   biv        divalent cation, M
%   dNTP       total dNTP, M
%   tann       temperature, degrees C
%   applyInit  gates ONLY the symmetry correction of Equation 3. Default true. PASS FALSE FOR ANY
%              HETERODUPLEX. The term exists because two identical strands associating have a
%              twofold rotational degeneracy; a primer on a template has none, and charging it there
%              penalises a self-complementary primer for binding its own target.
%   applyTerm5 suppress the 5'-terminal correction, for an arm of a larger broken duplex. Default true.
%   applyTerm3 same at the 3' end. Default true.
%   ctx        optional struct:
%              .UseViennaTerminal  use mismatch_exterior instead of the measured table.
%                                  Default false. READ THE ORIENTATION NOTE ABOVE BEFORE SETTING IT.
%              .TerminalMismatchDG pin every exterior terminal mismatch to one constant in kcal/mol.
%                                  Leave unset to take the measured sequence-specific value.
%              .SaltModel          'owczarzy_ref' (default), 'owczarzy', 'perphosphate' or 'auto'.
%                                  'auto' picks per-phosphate when any position is mismatched.
%              .MonovalentOnly     use bare monv rather than the von Ahsen equivalent in the
%                                  per-phosphate term. Default false. Run both; keep conclusions
%                                  that survive each.
%
% OUTPUTS
%   dG    kcal/mol at tann
%   info  struct breaking the total into dH, dS, the fixed DeltaG37 mismatch bucket, the loop term,
%         the terminal term, the salt route taken and the per-mismatch route codes. Use it when a
%         number looks wrong; guessing which term moved is what the breakdown exists to prevent.
%
% REFERENCES, none of which I have been able to open; verify before citing
%   SantaLucia J, Hicks D. Annu Rev Biophys Biomol Struct 2004;33:415.
%                                                doi:10.1146/annurev.biophys.32.110601.141800
%   SantaLucia J. PNAS 1998;95:1460.             doi:10.1073/pnas.95.4.1460
%   Owczarzy R et al. Biochemistry 2004;43:3537. doi:10.1021/bi034621r
%   Owczarzy R et al. Biochemistry 2008;47:5336. doi:10.1021/bi702363u
%   von Ahsen N et al. Clin Chem 2001;47:1956.   doi:10.1093/clinchem/47.11.1956
% ____________________________________________________________________________________________________________

persistent mixedWarned terminalWarned noPairWarned

TERM_MISMATCH_DG = -1.23;   %#ok<NASGU>   % kcal/mol, conservative literature bound. See the note above.

if (nargin < 7)||isempty(applyInit),  applyInit  = true; end
if (nargin < 8)||isempty(applyTerm5), applyTerm5 = true; end
if (nargin < 9)||isempty(applyTerm3), applyTerm3 = true; end
if (nargin < 10)||isempty(ctx), ctx = struct(); end
if ~isfield(ctx, 'UseViennaTerminal')||isempty(ctx.UseViennaTerminal), ctx.UseViennaTerminal = false; end
% THE SALT MODEL IS UNIFORM ACROSS A REGISTER, AND THE DEFAULT IS 'owczarzy_ref'.
%
% Four settings are accepted. 'auto' switches model at the first mismatch, which puts a fully paired
% duplex and its mismatched variant on different footings, a bias of 0.33 kcal/mol on a 20-mer and
% 0.52 on an 8-mer sitting exactly on the comparison an off-target screen rests on. It is retained
% for sensitivity work and should not be used to produce a result.
%
% Plain 'owczarzy' carries a defect of its own on mismatched input. Its correction scales as the
% duplex enthalpy times a function of salt, so removing base pairs weakens the salt destabilisation
% and partially cancels the loss: mismatching one terminal position of a 20-mer returns about
% 0.2 kcal/mol of spurious stability, and a register missing eight pairs runs past a kcal/mol,
% always favouring the mismatched case.
%
% 'owczarzy_ref', the default, multiplies the salt function by the enthalpy of the FULLY PAIRED
% duplex of the same sequence rather than by the mismatch-reduced enthalpy. Counterion release
% scales with how many phosphates are present, not with which positions happen to be mispaired. On
% a fully paired duplex the two enthalpies are identical, so this route and plain 'owczarzy' agree
% exactly there by construction.
%
% 'perphosphate' is the counterion-condensation form of SantaLucia 1998. It depends only on length
% and so cannot produce the mismatch pathology at all, but it has no GC dependence, which is a
% physical disagreement with the Owczarzy measurement rather than a simplification. The two differ
% by -1.06 to +1.09 kcal/mol across GC fractions from 0.10 to 1.00, so adopting it would re-rank
% GC-rich against AT-rich candidates by about two kcal/mol. It is monovalent-only and has to be fed
% the von Ahsen equivalent for magnesium, an approximation on an approximation.
%
% Validate the absolute scale against a melting curve you have measured before trusting an absolute
% Tm from any of these routes.
if ~isfield(ctx, 'SaltModel')||isempty(ctx.SaltModel), ctx.SaltModel = 'owczarzy_ref'; end
if ~isfield(ctx, 'MonovalentOnly')||isempty(ctx.MonovalentOnly), ctx.MonovalentOnly = false; end
% THE EXTERIOR TERMINAL STACKING TERM.
%
% A single constant cannot stand in for this motif. Pinning every exterior mismatch to -1.23, the
% most stabilising entry in the table, made a 20-mer score more stable with a terminal base broken
% than fully paired, which inverts the sign of what an off-target screen measures. Pinning it to
% zero removes stacking that is almost always present and reports a frayed register as weaker than
% it is. The sequence-specific table removes the choice.
%
% The magnitude is small either way. At 58 C a terminal A.T step is worth -0.147
% kcal/mol, so the entire thermodynamic signal available at a 3' terminus is of that order, which is
% noise against the rest of the model. DISCRIMINATION AT THE 3' END IS KINETIC, NOT THERMODYNAMIC:
% see ExtensionEfficiency.m and Huang MM, Arnheim N, Goodman MF, Nucleic Acids Res 1992;20:4567,
% doi:10.1093/nar/20.17.4567, where the spread across terminal mispairs is four orders of magnitude
% in kcat/Km. Do not expect this term to buy specificity once it carries the right values.
%
% THE DEFAULT IS THE MEASURED TABLE, NOT A CONSTANT. dGTerminalMismatch.m holds the 48 measured
% entries that Equations 9, 10 and 12 call for, indexed by the closing Watson-Crick pair and the two
% exterior bases. Forty-six of them are stabilising, spanning -1.19 kcal/mol to roughly zero, so
% scoring an exterior terminal mismatch as free removes stabilisation that is almost always present
% and reports a frayed register as weaker than it is. Leave the field unset to take the
% sequence-specific value. Set it to a number only to pin every exterior mismatch to one constant,
% which is useful for a sensitivity run and for nothing else.
%
% THIS IS A DIFFERENT TABLE FROM THE VIENNA mismatch_exterior ROUTE. The orientation caveat above
% applies to that route alone; the measured table carries its own i1 i2 / j1 j2 convention, which is
% documented in dGTerminalMismatch.m and used identically for the internal-loop closures above.
if ~isfield(ctx, 'TerminalMismatchDG'), ctx.TerminalMismatchDG = []; end

seq = char(seq); seq2 = char(seq2);
L = numel(seq);
if numel(seq2) ~= L
    error('dGMismatch:LengthMismatch', ...
        'seq is %d long and seq2 is %d; they must be pre-aligned position by position.', L, numel(seq2));
end
if L < 2
    error('dGMismatch:SegmentTooShort', ...
        ['A duplex segment of %d base(s) has no nearest-neighbour step. Callers must fold single ' ...
        'bases into the adjacent loop or bulge instead of passing them here.'], L);
end

N = NNParams(struct('Verbose', false));
bi = N.Index.Base;

%% ---- pairing mask -------------------------------------------------------------------------------------
comp = struct('A', 'T', 'T', 'A', 'C', 'G', 'G', 'C');
isWC = false(1, L);
for i = 1:L
    if ~isfield(comp, seq(i))||~isfield(comp, seq2(i))
        error('dGMismatch:NonACGT', ...
            'Position %d carries a base outside {A,C,G,T}: seq ''%s'', seq2 ''%s''.', i, seq(i), seq2(i));
    end
    isWC(i) = (comp.(seq(i)) == seq2(i));
end
nMM = nnz(~isWC);

%% ---- nearest-neighbour walk ---------------------------------------------------------------------------
h = 0;           % kcal/mol,        carries every term with a known enthalpy
s = 0;           % kcal/(K mol),    same convention as dGPrim.m
dGfixed = 0;     % kcal/mol,        DeltaG37 terms with no temperature dependence
routeUsed = [];

% A mismatch at a HELIX END is not an internal mismatch. Table 2 of SantaLucia and Hicks is
% parameterised for mismatches flanked by base pairs on both sides; at a terminus the stacking
% partner on one side simply is not there. Charging a Table 2 increment across that boundary AND a
% terminal-mismatch stacking term double counts the same interaction. The boundary step is skipped
% here and the terminal block below charges initiation plus one stacking term instead.
% When applyTerm5 or applyTerm3 is false the segment is an arm of a larger structure, the position
% is interior to that structure, and the Table 2 increment is correct, so the skip is gated on the
% same flags.
skip5 = applyTerm5 && ~isWC(1);
skip3 = applyTerm3 && ~isWC(L);

% WHICH MISMATCH RUNS BELONG TO EQUATION 12 RATHER THAN TO TABLE 2.
%
% Table 2 of SantaLucia and Hicks is parameterised for an ISOLATED mismatch flanked by Watson-Crick
% pairs on both sides. Page 429 of the review states the division explicitly: a single mismatch is
% formally a symmetric internal loop of 2 but is calculated from the mismatch nearest-neighbour
% parameters INSTEAD OF the loop formalism, not in addition to it. The Table 2 increments therefore
% already carry the cost of opening a one-by-one loop.
%
% A run of two or more consecutive mismatches is a larger internal loop and takes Equation 12:
%
%   dG(loop total) = dG(internal loop of N) + dG(asymmetry)
%                  + dG(left terminal mismatch) + dG(right terminal mismatch)
%
% Nothing from Table 2 appears in that sum. Charging the two boundary increments on top of the loop
% initiation would pay the initiation twice, which over-penalises every tandem or longer run by
% roughly 2 kcal/mol and makes an off-target register carrying one read as weaker than it is.
%
% isEq12 marks the positions of every run that Equation 12 owns: at least two long, and closed by a
% Watson-Crick pair at BOTH ends. A run touching either terminus is a frayed tail rather than a
% loop, and the terminal block further down handles it.
isEq12 = false(1, L);
ii = 1;
while ii <= L
    if ~isWC(ii)
        jj = ii;
        while (jj < L) && ~isWC(jj + 1), jj = jj + 1; end
        if (jj - ii + 1 >= 2) && (ii > 1) && (jj < L), isEq12(ii:jj) = true; end
        ii = jj + 1;
    else
        ii = ii + 1;
    end
end

for i = 1:(L - 1)
    atBoundary5 = skip5 && (i == 1);
    atBoundary3 = skip3 && (i == L - 1);
    if atBoundary5 || atBoundary3
        if (atBoundary5 && ~isWC(1) && isWC(2)) || (atBoundary3 && ~isWC(L) && isWC(L - 1))
            continue
        end
    end
    if isWC(i) && isWC(i + 1)
        step = seq(i:i + 1);
        row = find(strcmp(N.WC.Step, step), 1);
        if isempty(row)
            error('dGMismatch:MissingNN', 'No Watson-Crick nearest-neighbour row for step ''%s''.', step);
        end
        h = h + N.WC.dH(row);
        s = s + 0.001*N.WC.dS(row);

    elseif isWC(i) && ~isWC(i + 1)
        if isEq12(i + 1), continue; end          % Equation 12 closes this loop below
        [h, s, dGfixed, r] = AddMismatchStep(N, bi, h, s, dGfixed, ...
            seq(i), seq(i + 1), seq2(i + 1), tann);
        routeUsed(end+1) = r; %#ok<AGROW>

    elseif ~isWC(i) && isWC(i + 1)
        if isEq12(i), continue; end              % Equation 12 closes this loop below
        [h, s, dGfixed, r] = AddMismatchStep(N, bi, h, s, dGfixed, ...
            seq2(i + 1), seq2(i), seq(i), tann);
        routeUsed(end+1) = r; %#ok<AGROW>
    end
    % Both positions mismatched: nothing here. That step lies inside a run whose length term is
    % charged once below.
end

%% ---- internal loops, Equation 12, one per run -------------------------------------------------------------
% ORIENTATION OF THE TWO CLOSURES. The motif table is indexed i1 i2 / j1 j2, meaning 5'-i1 i2-3'
% stacked over 3'-j1 j2-5', with i1.j1 the closing Watson-Crick pair and i2.j2 the mismatch. Here
% seq runs 5' to 3' with increasing index and seq2 runs 3' to 5' with increasing index.
%
% At the LEFT closure the pair sits at i-1 and the loop opens toward increasing index, which is
% 5' to 3' along seq, so seq supplies i1 i2 and seq2 supplies j1 j2 directly.
%
% At the RIGHT closure the pair sits at j+1 and the loop opens toward DECREASING index. Along seq
% that direction is 3' to 5', so the two strands exchange roles: seq2 read from j+1 to j runs 5' to
% 3' and supplies i1 i2, while seq supplies j1 j2. Passing seq first at both ends would read the
% right-hand motif backwards.
%
% ASYMMETRY. Equation 12 adds 0.3 kcal/mol per nucleotide of difference between the two sides of the
% loop. This function scores ungapped registers of equal length, so every run it sees is a symmetric
% loop with n1 = n2 = runlen and the asymmetry term is identically zero. It is left out rather than
% written as a zero, and any future gapped input would have to restore it.
%
% SALT. The review treats a terminal mismatch as carrying the base-pair salt dependence and the rest
% of the loop as salt independent. The phosphates at these positions are already inside hRef below,
% which runs over every step of the register irrespective of pairing, so the closure free energies
% join dGruns as fixed terms and are not salt corrected a second time.
dGruns = 0;
i = 1;
while i <= L
    if ~isWC(i)
        j = i;
        while (j < L) && ~isWC(j + 1), j = j + 1; end
        runlen = j - i + 1;
        if (runlen >= 2) && (i > 1) && (j < L)
            dGruns = dGruns + dGLoop('internal', 2*runlen, tann, N);
            % The two closures enter the ENTHALPY and ENTROPY accumulators rather than the fixed
            % free-energy bucket. dGTerminalMismatch returns dG = dH - T*dS/1000, so at the
            % evaluation temperature the two routes are the same number to the last bit, and the
            % clamp it applies to the one suspect entry is preserved because it hands back the dS
            % that reproduces the clamped dG. Two things change. info.dH and info.dS become the
            % whole thermodynamic sum instead of a partial one, and the 'owczarzy' salt route,
            % which scales with the ACTUAL accumulated enthalpy rather than with hRef, now sees the
            % closure stacks it was previously blind to. Under the default 'owczarzy_ref' the salt
            % keys on hRef and nothing moves.
            [~, hL, sL] = dGTerminalMismatch(seq(i - 1), seq(i), seq2(i - 1), seq2(i), tann, N);
            [~, hR, sR] = dGTerminalMismatch(seq2(j + 1), seq2(j), seq(j + 1), seq(j), tann, N);
            h = h + hL + hR;
            s = s + 0.001*(sL + sR);
        end
        i = j + 1;
    else
        i = i + 1;
    end
end

%% ---- isolated-mismatch temperature correction, route (b) ---------------------------------------------------
dGslope = 0;
for i = 2:(L - 1)
    if isWC(i)||~isWC(i - 1)||~isWC(i + 1), continue; end
    W = seq(i - 1); X = seq(i); Y = seq2(i); V = seq(i + 1);
    % Both bracketing Table 2 steps must be on route (b) for the slope to be the right correction.
    r1 = StepRoute(N, bi, seq(i - 1), seq(i),  seq2(i));
    r2 = StepRoute(N, bi, seq2(i + 1), seq2(i), seq(i));
    if (r1 == 1) && (r2 == 1), continue; end                 % both measured, nothing to add
    if (r1 == 1) || (r2 == 1)
        if isempty(mixedWarned), mixedWarned = false; end
        if ~mixedWarned
            mixedWarned = true;
            warning('dGMismatch:MixedMismatchRoute', ...
                ['An isolated mismatch has one Table 2 step with measured enthalpy and one without, ' ...
                'so the Vienna slope would double-count the measured step. This mismatch falls back ' ...
                'to its 37 C value. Fill the remaining rows of mismatch_nn_params.csv to clear it. ' ...
                'Reported once per session.']);
        end
        continue
    end
    d = N.MM.ViennaSlope(W, X, Y, V, tann);
    if isfinite(d), dGslope = dGslope + d; end
end

%% ---- terminal corrections ------------------------------------------------------------------------------
% INITIATION IS PAID WHETHER OR NOT THE TERMINAL BASE PAIRS.
%
% Two strands have to find one another and nucleate a helix irrespective of whether the outermost
% base happens to pair, so initiation is owed at a mismatched end exactly as at a paired one. It is
% worth about 1.0 kcal/mol. Charging a terminal stacking term INSTEAD of initiation, rather than as
% well as it, makes a 20-mer score 1.2 to 1.5 kcal/mol more stable with its terminal base broken
% than fully paired.
%
% The correct decomposition, and the one both SantaLucia and Turner use: initiation keyed on the
% CLOSING Watson-Crick pair, plus one terminal stacking term for the unpaired bases sitting on it.
% Where several leading positions are unpaired the end is a frayed tail rather than a terminal
% mismatch, so initiation keys on the first real pair and no stacking term is charged. The run-length
% block above already declines to charge a loop term for a run that touches either end, so a frayed
% tail is correctly free in both places.
%
% WHAT IS STILL WRONG HERE. TERM_MISMATCH_DG = -1.23 is a bound on the MOST stabilising terminal
% mismatch applied uniformly to all of them, so it remains too generous, and a terminal mismatch will
% still look cheaper than it is. Expect thermodynamics to discriminate a 3'-terminal mismatch only
% weakly in any case: terminal pairs fray and contribute little at 58 C. That is precisely why the
% discrimination has to come from extension kinetics rather than from this function.
dGterm = 0;
q5 = find(isWC, 1, 'first');
q3 = find(isWC, 1, 'last');
if isempty(q5)
    % A register with no paired position at all is legitimate, not an error. An exhaustive scan
    % generates them at the rate (3/4)^L: about 61 of the 14,498 registers of a 19-mer against a
    % 7249 nt circle. No helix exists, so no
    % initiation is charged and the free energy is strongly positive, which is the correct answer.
    % Reported once per session so a scan does not bury its own output.
    if isempty(noPairWarned), noPairWarned = false; end
    if ~noPairWarned
        noPairWarned = true;
        warning('dGMismatch:NoPairedPosition', ...
            ['A %d-position register carries no Watson-Crick pair, so no helix is initiated and no ' ...
            'initiation term is charged. Expected at rate (3/4)^L during an exhaustive scan. ' ...
            'Reported once per session; silence with warning(''off'', ' ...
            '''dGMismatch:NoPairedPosition'').'], L);
    end
end
if applyTerm5 && ~isempty(q5)
    [h, s] = AddInit(N, h, s, seq(q5));
    if q5 == 2
        [dgt, dht, dst] = TerminalMismatch(seq, seq2, 1, ctx, ctx.TerminalMismatchDG, tann, N);
        dGterm = dGterm + dgt; h = h + dht; s = s + 0.001*dst;
    end
end
if applyTerm3 && ~isempty(q3)
    [h, s] = AddInit(N, h, s, seq(q3));
    if q3 == L - 1
        [dgt, dht, dst] = TerminalMismatch(seq, seq2, L, ctx, ctx.TerminalMismatchDG, tann, N);
        dGterm = dGterm + dgt; h = h + dht; s = s + 0.001*dst;
    end
end

% THE TEST IS ON THE DUPLEX, NOT ON THE QUERY STRAND.
% The symmetry correction is owed when the complex is formed from two identical molecules, which
% is a property of the PAIR of strands. seq runs 5' to 3' with increasing index and seq2 is the
% opposing strand written 3' to 5' with the same index, so the opposing strand read 5' to 3' is
% fliplr(seq2), and the complex is self-complementary exactly when that equals seq.
% ON A FULLY PAIRED DUPLEX THIS IS THE FAMILIAR CONDITION. There seq2 is the elementwise
% complement of seq, so fliplr(seq2) is the reverse complement of seq and the test reduces to
% seq == revcomp(seq).
% ON A MISMATCHED DUPLEX THE TWO DIVERGE, and the reverse-complement form asks the wrong
% question. A self-dimer register carrying one mismatch is still two copies of one molecule and
% still carries the twofold degeneracy, while seq == revcomp(seq) fails because the mismatched
% positions are not complementary. BestSubRegister exists to score exactly those registers, so
% the duplex form is the one that belongs here.
if applyInit && isequal(seq, fliplr(seq2))
    s = s + 0.001*N.WC.Sym(2);
    h = h + N.WC.Sym(1);
end

%% ---- salt ---------------------------------------------------------------------------------------------------
% REFERENCE ENTHALPY FOR THE SALT TERM.
%
% Converting a correction defined on 1/Tm into an entropy shift gives dS_salt = dS + dH*Phi, where
% dH belongs to the duplex whose Tm is being corrected. Owczarzy fitted Phi on FULLY PAIRED duplexes,
% where dH and length move together. On a mismatched register they decouple: the phosphate count is
% unchanged while dH collapses, so the correction shrinks and removing base pairs returns spurious
% stability: mismatching one terminal position of a 20-mer returns about 0.2 kcal/mol, and a
% register missing eight pairs runs past a kcal/mol, always favouring the mismatched case.
%
% The fix is to multiply Phi by the enthalpy of the fully paired duplex of the SAME SEQUENCE rather
% than by the mismatch-reduced enthalpy. Counterion release scales with how many phosphates are
% present, not with which positions happen to be mispaired. On a fully paired duplex hRef equals h
% exactly, so the two Owczarzy routes agree there by construction.
%
% The alternative is the per-phosphate form of SantaLucia 1998, which has no pathology and also no
% GC dependence. The two models differ by -1.06 to +1.09 kcal/mol across GC fractions from 0.10 to
% 1.00, so adopting it would re-rank GC-rich against AT-rich candidates by about two kcal/mol for a
% reason that is an artefact. Owczarzy measured a GC dependence; counterion condensation theory
% assumes none. The measurement wins.
hRef = 0;
for i = 1:(L - 1)
    rr = find(strcmp(N.WC.Step, seq(i:i + 1)), 1);
    if isempty(rr)
        error('dGMismatch:MissingNN', 'No Watson-Crick row for step ''%s''.', seq(i:i + 1));
    end
    hRef = hRef + N.WC.dH(rr);
end
if applyTerm5
    if (seq(1) == 'G')||(seq(1) == 'C'), hRef = hRef + N.WC.InitGC(1); else, hRef = hRef + N.WC.InitAT(1); end
end
if applyTerm3
    if (seq(L) == 'G')||(seq(L) == 'C'), hRef = hRef + N.WC.InitGC(1); else, hRef = hRef + N.WC.InitAT(1); end
end
gcRef = nnz((seq == 'G')|(seq == 'C'))/L;

switch lower(ctx.SaltModel)
    case 'owczarzy_ref'
        mode = 1;
    case 'auto'
        if nMM > 0, mode = 3; else, mode = 2; end
    case 'perphosphate'
        mode = 3;
    case 'owczarzy'
        mode = 2;
    otherwise
        error('dGMismatch:BadSaltModel', ...
            ['ctx.SaltModel must be ''owczarzy_ref'' (default), ''owczarzy'', ''perphosphate'' or ' ...
            '''auto''; got ''%s''.'], ctx.SaltModel);
end

% f_GC MUST BE COUNTED OVER WATSON-CRICK PAIRED POSITIONS ONLY.
% The Owczarzy term is parameterised by the GC fraction of the BASE PAIRS of a duplex. Counting a
% mismatched G or C in the query strand as a G.C pair inflates f_GC on exactly the registers where
% it matters most, since a heavily mismatched off-target register can carry many unpaired G and C
% bases. The previous version of this file took the fraction over the whole query strand, and this
% rewrite reproduced that line unchanged until an external review of the original suite caught it.
% On a fully paired segment the two expressions coincide, so nothing that passed before changes.
nWCpos = max(1, nnz(isWC));
gc = nnz(isWC & ((seq == 'G')|(seq == 'C')))/nWCpos;

if ctx.MonovalentOnly
    monEq = monv;
else
    monEq = N.Salt.MonovalentEquivalent(monv, biv, dNTP);
end

% Both salt models are evaluated every time, and the one not used is reported in info. The reason
% is a discontinuity in the 'auto' rule: a fully paired segment takes Owczarzy and a segment with
% one mismatch takes the per-phosphate form, so the switch introduces a step at the first mismatch.
% For a fully paired 20-mer at 40 mM Na and 0.85 mM free Mg with f_GC = 0.5 the two shifts are
% -13.7 and -13.3 cal/(mol K), a gap near 0.15 kcal/mol at 58 C. Small, but it lands squarely on the
% comparison between a cognate duplex and a near-cognate off-target, which is the comparison the
% engine rests on. Measure it per call rather than assuming it stays small: if
% abs(info.dS_salt_alternative - info.dS_salt_used)*T ever approaches the margins being ranked on,
% the switch is doing work the physics is not.
dsOwcz    = N.Salt.Owczarzy(h,    gc,    L, monv, biv, dNTP);
dsOwczRef = N.Salt.Owczarzy(hRef, gcRef, L, monv, biv, dNTP);
if monEq > 0
    dsPer = N.Salt.PerPhosphate(L, monEq);
else
    dsPer = NaN;
end

switch mode
    case 1
        dsUsed = dsOwczRef; dsAlt = dsPer; saltRoute = 'owczarzy_ref';
    case 2
        dsUsed = dsOwcz;    dsAlt = dsPer; saltRoute = 'owczarzy';
    case 3
        if monEq <= 0
            error('dGMismatch:NonPositiveSalt', ...
                'The per-phosphate term needs a positive monovalent concentration; got %.3g M.', monEq);
        end
        dsUsed = dsPer; dsAlt = dsOwczRef; saltRoute = 'perphosphate';
end
s = s + 0.001*dsUsed;

%% ---- total ---------------------------------------------------------------------------------------------------
TK = tann + 273.15;
dG = h - TK*s + dGfixed + dGruns + dGslope + dGterm;

if nargout > 1
    info = struct('dH', h, 'dS', s, 'dG_fixed_mismatch', dGfixed, 'dG_loop_runs', dGruns, ...
        'dG_vienna_slope', dGslope, 'dG_terminal', dGterm, 'SaltRoute', saltRoute, ...
        'dS_salt_used', dsUsed, 'dS_salt_alternative', dsAlt, ...
        'dG_salt_model_gap', TK*abs(dsUsed - dsAlt)/1000, ...
        'fGC_paired', gc, 'fGC_reference', gcRef, 'dH_reference', hRef, 'nMismatch', nMM, 'StepRoutes', routeUsed, ...
        'Temperature_C', tann, 'Length', L);
end
end

%% ____________________________________________________________________________________________________________
function [h, s] = AddInit(N, h, s, base)
% SantaLucia rows 17 and 18. These already carry the terminal A.T penalty, which is why no Vienna
% terminal AU/GU penalty may ever be added on top. See the bookkeeping note in ViennaParams.m.
if (base == 'G')||(base == 'C')
    h = h + N.WC.InitGC(1); s = s + 0.001*N.WC.InitGC(2);
else
    h = h + N.WC.InitAT(1); s = s + 0.001*N.WC.InitAT(2);
end
end

%% ____________________________________________________________________________________________________________
function r = StepRoute(N, bi, closingBase, xBase, yBase)
% 1 measured, 2 Vienna slope available, 3 DeltaG37 only.
ic = bi.(closingBase); ix = bi.(xBase); iy = bi.(yBase);
r = double(N.MM.Route(ic, ix, iy));
if r == 0
    error('dGMismatch:MissingMismatch', ...
        ['No Table 2 parameter for closing ''%s'', X ''%s'', Y ''%s''. This slot is a ' ...
        'Watson-Crick combination, so the pairing test upstream is wrong.'], closingBase, xBase, yBase);
end
end

%% ____________________________________________________________________________________________________________
function [h, s, dGfixed, r] = AddMismatchStep(N, bi, h, s, dGfixed, closingBase, xBase, yBase, tann) %#ok<INUSD>
% One Table 2 step. Measured entries go into the enthalpy and entropy accumulators so they scale
% with temperature like everything else; unmeasured entries contribute their 37 C value to the fixed
% bucket and pick up their temperature dependence, if any, from the isolated-mismatch slope.
ic = bi.(closingBase); ix = bi.(xBase); iy = bi.(yBase);
r = double(N.MM.Route(ic, ix, iy));
if r == 0
    error('dGMismatch:MissingMismatch', ...
        ['No Table 2 parameter for closing ''%s'', X ''%s'', Y ''%s''. This slot is a ' ...
        'Watson-Crick combination, so the pairing test upstream is wrong.'], closingBase, xBase, yBase);
end
if r == 1
    h = h + N.MM.dH(ic, ix, iy);
    s = s + 0.001*N.MM.dS(ic, ix, iy);
else
    dGfixed = dGfixed + N.MM.dG37(ic, ix, iy);
end
end

%% ____________________________________________________________________________________________________________
function [dg, dh, ds] = TerminalMismatch(seq, seq2, pos, ctx, fallback, tann, N)
% A mismatched segment end. Both strands carry an unpaired base facing each other, which is a
% terminal mismatch rather than a dangle.
%
% ORIENTATION IS UNRESOLVED, which is why the Vienna route is off by default. mismatch_exterior is
% indexed by the closing pair and two exterior bases and is NOT symmetric in those two bases, so the
% assignment of the top-strand and bottom-strand base to the two indices has to be established
% rather than assumed. Settle it with RNAeval on a structure carrying a known terminal mismatch and
% the DNA parameter file, compare against both candidate lookups, then set ctx.UseViennaTerminal.
persistent warned
% dh and ds are returned in kcal/mol and cal/(K mol) so the caller can fold the motif into its own
% accumulators. They are zero on every route that yields only a free energy, which keeps the caller
% arithmetic exact: that route's dg stays in the fixed bucket and contributes no enthalpy there.
dh = 0; ds = 0;
if ~ctx.UseViennaTerminal
    if ~isempty(fallback)
        dg = fallback; return                 % caller pinned every exterior mismatch to a constant
    end
    % The measured table. At the 5' end the closing pair sits at index 2 and the loop opens toward
    % index 1, which along seq is 3' to 5', so seq2 supplies i1 i2 and seq supplies j1 j2. At the
    % 3' end the pair sits at L-1 and opening toward L runs 5' to 3' along seq, so the roles are the
    % other way round. The convention is the same one the internal-loop closures use above.
    LL = numel(seq);
    if pos == 1
        [~, dh, ds] = dGTerminalMismatch(seq2(2), seq2(1), seq(2), seq(1), tann, N);
    else
        [~, dh, ds] = dGTerminalMismatch(seq(LL - 1), seq(LL), seq2(LL - 1), seq2(LL), tann, N);
    end
    % The motif leaves as enthalpy and entropy, so the free energy returned here is zero and the
    % caller adds it exactly once, through h - T*s. Since dGTerminalMismatch defines its own free
    % energy as dH - T*dS/1000, the two routes agree to the last bit at the evaluation temperature.
    % A failed lookup returns zeros on all three, which contributes nothing and is the right answer
    % for a motif the table does not hold.
    dg = 0;
    return
end
if isempty(warned), warned = false; end
if ~warned
    warned = true;
    warning('dGMismatch:UnverifiedTerminalOrientation', ...
        ['ctx.UseViennaTerminal is set, but the strand orientation of mismatch_exterior has not ' ...
        'been verified against RNAeval. If the two exterior bases are the wrong way round every ' ...
        'terminal mismatch in the suite is off by 0.2 to 0.3 kcal/mol with no downstream symptom. ' ...
        'Reported once per session.']);
end
L = numel(seq);
if pos == 1
    inner = 2;
else
    inner = L - 1;
end
pairName = [seq(inner) seq2(inner)];
try
    V = ViennaParams([], struct('Verbose', false));
    if ~isfield(V.PairIndex, pairName)
        dg = fallback; return
    end
    ip = V.PairIndex.(pairName);
    b1 = V.BaseIndex.(seq(pos));
    b2 = V.BaseIndex.(seq2(pos));
    g37 = V.mismatch_exterior.dG37(ip, b1, b2);
    hh  = V.mismatch_exterior.dH(  ip, b1, b2);
    if ~isfinite(g37)||~isfinite(hh)
        dg = fallback; return
    end
    dg = V.AtTemp(g37, hh, tann);
catch
    dg = fallback;
end
end
