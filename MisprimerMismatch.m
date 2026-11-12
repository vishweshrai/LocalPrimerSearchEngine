function [dGpTcmin, dGpTc, dGpTmin, dGpT, check, misp] = MisprimerMismatch(lintemp, primer, phi, pcon, T, monv, biv, dNTP, tmp, dTsub, minRun, endTol, endWin, useEps, ryc, tAnneal, surveyRun)
%% MISPRIMERMISMATCH  Whole-template mispriming screen by exhaustive register enumeration
% Written by Vishwesh Rai, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 09-11-2026 00:48 IST
%
% Revision history is NOT kept here. See MisprimerMismatch_<YYYYMMDD>.txt beside this file.
% ____________________________________________________________________________________________________________
% Two scans run per primer. Both look for IDENTITY between the primer and the strand being scanned,
% because a primer identical to a stretch of one strand anneals to the other:
%   scan C: the primer against lintemp, the strand carrying its cognate site, with that site
%           excluded. Identity here means the primer would anneal to the COMPLEMENTARY strand at
%           that position, i.e. mispriming on the strand it is supposed to bind. The cognate site
%           sits on the same molecule, which is why dGMispTarget solves a cubic over four occupancy
%           states.
%   scan T: the primer against seqrcomplement(lintemp). Identity there means the primer would anneal
%           to lintemp itself, which it should never do. That molecule carries no cognate site for
%           this primer, which is why dGMispOffTarget solves a quadratic.
%
% ____________________________________________________________________________________________________________
% DESIGN NOTES: THE ANCHORED ACCUMULATOR, THE EXTENSION WEIGHTING AND THE SURVEY FLOOR
%
%  F1. THE ANCHORED ACCUMULATOR. dGpTc and dGpT were the most stable register found anywhere in the
%      survey, anchored or not, and the engine turned exactly that number into the discrimination
%      margin of cascade stage 13. A register whose duplex stops two bases short of the primer 3'
%      terminus therefore set the margin, and a candidate could be eliminated on the strength of a
%      site the polymerase cannot extend. Taq carries no 3'->5' proofreading exonuclease, so an
%      unpaired primer terminus is not trimmed and the insertion site will not accept it
%      (Huang et al. 1992 NAR 20:4567). Worked example on GACGCTTACTGGAACGGTAAA against M13mp18:
%      the register at (+) 5384-5395 scores -5.05 kcal/mol with f = 8, l = 19, n = 21, so it carries
%      a two-base 3' flap and cannot prime, yet it governed every margin and every Score term
%      derived from the forward primer.
%      SeedScan now carries dGanch in parallel with dGbest. dGanch admits a register only when
%      l >= n - endTol, which is the same terminus test the rejection rule already applied, so the
%      number that feeds the margin and the number that decides pass or fail are finally the same
%      quantity. dGbest survives unchanged in the sixth output as a diagnostic.
%      dGanch stays at 0 when no anchored register exists. The margin is then -threshold, which is
%      large and positive, and that is the correct reading: no extendable off-target site at all.
%
%  F2. PER-LOCUS TERMINAL-WINDOW dG. For the most stable anchored register of each scan, the
%      stacking-only free energy of the last endWin primer bases is reported separately from the
%      whole-register total. Initiation, terminal and symmetry terms are suppressed through the
%      applyInit/applyTerm5/applyTerm3 flags of dGMismatch, so this is the plain nearest-neighbour
%      stack sum in the same currency as the Rychlik (1995) 3'-end cutoff that cascade stage 6
%      applies to the primer in isolation. This is the quantity that separates a 9 bp duplex ending
%      on a matched 3'-OH from a 12 bp duplex ending in a flap; total dG does not.
%      endWin defaults to 5, matching EndWin in the engine.
%
%  F3. WRAP-SPANNING REGISTERS ARE NO LONGER DISCARDED. SeedScan searched the doubled string and
%      then executed "if q > L, continue; end". Consider a register starting near L whose only
%      minRun-identity block sits at primer offset o with st + o - 1 > L. The seed occurrence lies
%      past the seam, q > L, and the register was dropped. No other offset recovers it, because no
%      other qualifying seed exists in that register. Wrap-around handling is advertised at the top
%      of this file as the specific advantage over localalign, and for this class of register it was
%      not happening. The guard was also redundant: seen(st) already prevents a circular site from
%      being scored twice, once per copy of the doubled string. Deleted.
%
%  F4. SURVEYRUN. CORRECTED AT H9; THE EARLIER CLAIM HERE WAS WRONG. This note previously said the
%      test "has never once evaluated true". That is false, and the error was in reading surveyRun
%      as a filter on whole REGISTERS. It is not: it is passed to BestSubRegister as its minrun and
%      gates each SUB-RANGE of the footprint. Plenty of sub-ranges carry no run of three. Simulated
%      over 4000 realistic registers on a 21-mer with one exact 5-mer seed and background identity
%      at the expected one-in-four rate, it skips 38 per cent of the sweep at surveyRun = 3, 75 per
%      cent at 5 and 91 per cent at 6.
%      What is true is narrower and worth stating precisely: it cannot change the RESULT, because
%      the sub-range that wins contains the seed and so carries a run of SEEDLEN. It is a live cost
%      filter that happens to be inert with respect to the answer, which is not the same as a test
%      that never fires.
%      RAISING IT ABOVE minRun IS UNSAFE, and the earlier note recommended exactly that. The
%      rejection rule acts on sub.runRight >= minRun, and sub.runRight can only come from a
%      sub-range that survived this filter, so a surveyRun above minRun silently discards the very
%      registers the rule exists to catch. SeedScan now errors on that setting instead of
%      suggesting it.
%
%  F5. SEEDLEN was described as held BELOW minRun. With the engine passing minRun = 5 the expression
%      min(5, minRun) returns 5, so it is held EQUAL. Completeness is untouched, since the proof
%      needs only SEEDLEN <= minRun. The stated widening of the survey is not happening at this
%      setting, and the comment below now says so. Lowering SEEDLEN to 4 would restore it at roughly
%      four times the register count; that is left as a deliberate decision rather than a silent one.
%
%  F6. The misp output struct is new. Every field is a scalar or a string scalar, so the engine can
%      drop it straight into struct2table without unpacking.
% ____________________________________________________________________________________________________________
%
% WHY THE SEARCH IS AN ENUMERATION RATHER THAN localalign
%
% localalign ranks by the NUC44 substitution matrix, +5 for an identity and -4 for a mismatch. That
% score is a linear function of how many bases pair. Free energy is not a function of base pairs at
% all, it is a function of ADJACENT pairs, and NUC44 carries no stacking term, no initiation, no
% terminal mismatch, no salt and no temperature. Measured on GACGCTTACTGGAACGGTAAA against the
% M13mp18 minus strand, 1551 scorable registers at 55 C:
%
%     rank   --- top 10 by NUC44 ---        --- top 10 by dG ---
%             pos    score      dG           pos      dG    score
%       1    5067      105   -13.66         5067   -13.66     105
%       2    1315       47    -0.52         1651    -2.00      31
%       3    1938       44   +12.66          199    -1.38      21
%
%     overlap between the two top-10 sets: 2 of 10
%     most stable off-target register: dG -2.00 at 1651, NUC44 rank 42 of 1551
%
% COMPLETENESS
%
% The acceptance rule demands a contiguous run of minRun identities. Any register that can trigger a
% rejection therefore contains a run of minRun identities, hence contains an exact SEEDLEN-mer of the
% primer whenever SEEDLEN <= minRun. Enumerating every exact SEEDLEN-mer match with strfind and
% scoring each implied register with dGMismatch is therefore EXHAUSTIVE with respect to the rule
% rather than a top-N sample of it, and cheaper, since string search is linear where localalign fills
% a 7249 by 21 dynamic programming matrix twice.
%
% The cognate site is excluded by coordinate rather than masked with 'N'. The old mask depended on
% 'N' behaving sensibly inside NUC44 and inside seqcomplement, and the
% seqrcomplement(seqrcomplement(...)) round trip wrapped around it was a no-op. Note that the
% exclusion drops any register OVERLAPPING [exclLo, exclHi], not merely the cognate register itself,
% which costs 2*plen-1 start positions and is the conservative choice.
%
% InternalLoopBridge still detects two-sided internal loops, which a gapless enumeration cannot
% represent, and still takes a localalign-shaped struct. That struct is built from the NUMALN most
% STABLE registers rather than the NUMALN highest-scoring alignments. Its Start and Stop columns for
% the TARGET are left unwrapped on purpose: InternalLoopBridge measures loop size as a coordinate
% difference, and folding those coordinates modulo L would manufacture small, spurious loop sizes
% between arms that sit at unrelated loci.
%
% REFERENCES
%   SantaLucia J, Hicks D. Annu Rev Biophys Biomol Struct 2004;33:415.  doi:10.1146/annurev.biophys.32.110601.141800
%   Huang MM, Arnheim N, Goodman MF. Nucleic Acids Res 1992;20:4567.    doi:10.1093/nar/20.17.4567
%   Rychlik W. Mol Biotechnol 1995;3:129.                               doi:10.1007/BF02789108
% ____________________________________________________________________________________________________________

if (nargin < 10)||isempty(dTsub),  dTsub  = 5;  end
if (nargin < 11)||isempty(minRun), minRun = 7;  end
if (nargin < 12)||isempty(endTol), endTol = 3;  end
if (nargin < 13)||isempty(endWin), endWin = 5;  end
if (nargin < 14)||isempty(useEps), useEps = true; end                                        % (G7)
% RYCHLIK CONTINUOUS-TERMINAL-DUPLEX GATE.                                                     (G9)
% Rychlik W, Biotechniques 1995;18(1):84-90, PMID 7702859, measured Taq and Pfu
% priming from primers mismatched over progressively longer 5' termini and reported that, at
% standard PCR conditions, priming became DETECTABLE once the 3'-terminal portion of the primer
% formed a continuous duplex more stable than -11 kcal/mol with the target. That is a measured
% detection threshold on the continuous run of Watson-Crick pairs carrying the 3' terminus, which
% is a different object from the whole anchored register and from the fixed pentamer that stage 6
% of the engine screens. BestSubRegister returns the run (trunRight, iTrunRight) as of G9 and the
% run is scored here.
%
% OPEN QUESTION, recorded rather than hidden: the paper's reference state for that -11 is not
% established from its abstract, and OLIGO-era figures were not always quoted at 37 C / 1 M Na+.
% The reference state is therefore a named field rather than a literal, so correcting it is one
% assignment. Until the methods section is checked, ryc.T = 37 and ryc.Na = 1 are an assumption,
% and the gate is worth perhaps 1 to 2 kcal/mol of uncertainty in consequence.
if (nargin < 15)||isempty(ryc), ryc = struct(); end
if ~isfield(ryc, 'On')  ||isempty(ryc.On),  ryc.On  = true;  end
if ~isfield(ryc, 'Max') ||isempty(ryc.Max), ryc.Max = -11.0; end
if ~isfield(ryc, 'T')   ||isempty(ryc.T),   ryc.T   = 37;    end
if ~isfield(ryc, 'Na')  ||isempty(ryc.Na),  ryc.Na  = 1;     end
% EVALUATION TEMPERATURE.                                                                      (H9)
% Mispriming is judged at the temperature the reaction reaches, not at an offset below this
% particular primer's Tm. In an asymmetric reaction the limiting primer's Tm is deliberately higher
% than the excess primer's, so the per-primer convention screened it several degrees hotter than
% the annealing step and understated how stable its off-target duplexes will be. Falls back to the
% old convention when the annealing temperature is not resolvable.
if (nargin < 16), tAnneal = []; end
if ~isempty(tAnneal)&&isfinite(tAnneal)
    tann = tAnneal;
else
    tann = tmp - dTsub;
end

SEEDLEN   = min(5, minRun);  % exact identity seeding the register search; must be <= minRun.  (F5)
                             % At minRun = 5 this is EQUAL to minRun, not below it, so the survey
                             % is no wider than the rejection-eligible set. Lower it to 4 to widen.
% SURVEY FLOOR: the contiguous paired run a SUB-REGISTER must carry to be COSTED.
% Set by the caller, and defaulting to minRun when it is not supplied, which is the setting that
% cannot change any reported number.
%
% WHY IT IS A SETTING AND NOT A CONSTANT. The two ensemble sums exist to be compared: fluxAnch
% runs over every anchored register the survey finds, fluxAnchRun over the subset the rejection
% acts on. BestSubRegister refuses to score any sub-register whose longest run falls below the
% floor it is handed, so with the floor AT minRun every register that reaches the accumulators
% already satisfies the run condition and the two sums are the same number by construction. The
% calibration they were built for then returns a difference of zero for a reason that has nothing
% to do with the physics. A floor below minRun is what gives the comparison something to measure.
%
% WHAT A LOWER FLOOR COSTS. Each nucleotide removed multiplies the admitted sub-register count by
% roughly four, which is the same scaling the MaxRunInt note in PrimerSearchParameters records.
%
% WHAT A LOWER FLOOR DOES NOT DO. SEEDLEN still decides which WINDOWS are examined at all, so
% lowering this floor widens the sub-register sweep inside the windows the seed already found
% rather than finding new ones. If a calibration run still shows the two sums close, SEEDLEN is
% the next lever and not this one.
if (nargin < 17)||isempty(surveyRun), surveyRun = minRun; end
SURVEYRUN = surveyRun;
NUMALN    = 10;              % registers handed to InternalLoopBridge, most stable first

check = true;
primer = upper(char(primer)); lintemp = upper(char(lintemp));

% THE PRIMER-TEMPLATE FREE ENERGY IS COMPUTED HERE, NOT SUPPLIED.                             (H10)
% dGp is the free energy of this primer against its perfect Watson-Crick complement, and it is the
% REFERENCE both mispriming thresholds are built from: dGMispTarget and dGMispOffTarget ask how
% stable an off-target duplex would have to be to steal a fraction phi of the primer away from its
% real site, and that question has no meaning without knowing how tightly the real site is held.
%
% It used to arrive as an argument, carrying whatever temperature the caller had evaluated it at.
% That was consistent while tann was tmp - dTsub, because the engine computed it at the same
% offset below the same Tm. It stopped being consistent the moment tann became the reaction's own
% annealing temperature, since the caller still held a value taken at the primer's own reference.
% For a limiting primer whose Tm sits several degrees above the excess primer's, those two
% temperatures differ by the whole LATE-PCR offset: measured at 3.573 kcal/mol on
% CGGAAGCATAAAGTGTAAAGCCTGGG, which through exp(-dG/RT) is a factor of about 230 in the binding
% constant and moves both thresholds with it.
% Computing it here removes the possibility outright. There is no advisory argument to drift, and
% no warning needed to police one.
dGp = dGPrim(primer, monv, biv, dNTP, tann);
plen = strlength(primer);

if strlength(lintemp) < plen
    error('MisprimerMismatch:TemplateTooShort', ...
        'Template is %d nt and the primer is %d nt; the circular register frame is undefined.', ...
        strlength(lintemp), plen);
end

dGpTcmin = NaN; dGpTc = 0; dGpTmin = NaN; dGpT = NaN;
misp = EmptyMisp();

%% ---- scan C: the cognate strand, with the cognate site excluded --------------------------------
% The engine rotates the template so the primer's own site occupies positions 1..plen.
dGpTcmin = dGMispTarget(dGp, phi, pcon, T, tann);
RTm = 1.987E-3*(273.15 + tann);
[scC, alstructC] = SeedScan(primer, lintemp, monv, biv, dNTP, tann, dGpTcmin, ...
    SEEDLEN, minRun, SURVEYRUN, endTol, endWin, NUMALN, 1, plen, useEps, ryc);
dGpTc = scC.dGbest;
misp.ScannedC      = true;
misp.dGpTcAnch     = scC.dGanch;
misp.dGpTcAnchEff  = scC.dGanchEff;
% The ensemble sums are formed AFTER the internal-loop bridge below, not here. See the note at
% that block for why the order matters.
misp.dGpTcEps      = scC.eps;
misp.dGpTcEpsPair  = scC.epsPair;
misp.dGpTcEnd   = scC.dGend;
misp.dGpTcPos   = scC.pos;
misp.dGpTcEndRn = scC.endrun;
misp.dGpTcAlign = scC.align;
misp.nSurvC     = scC.nsurv;
misp.dGpTcTrun   = scC.dGtrun;                                                                % (G9)
misp.dGpTcTrunRn = scC.trun;                                                                  % (G9)
misp.dGpTcTrunMin = scC.dGtrunMin;                                                            % (G9)
if scC.offend
    check = false;
end
if scC.rych                                                                                   % (G9)
    check = false;
end
if check == true
    [dGpTcil, ilfound] = InternalLoopBridge(alstructC, primer, lintemp, monv, biv, dNTP, ...
        tann, endTol);
    if ilfound
        dGpTc = min(dGpTc, dGpTcil);
        % THE ACCUMULATORS MOVE WHENEVER A BRIDGE IS FOUND, NOT ONLY WHEN ONE IS REJECTED.
        % A bridge that is stable but still weaker than the rejection threshold is a real structure
        % competing for the same primer, and the reported margin is the distance from the WORST
        % structure to the threshold. Updating these only inside the rejection branch left such a
        % bridge out of the margin entirely, so a candidate could report headroom it did not have.
        % The rejection itself still keys on the threshold, below.
        misp.dGpTcAnch    = min(misp.dGpTcAnch, dGpTcil);
        % A bridge is two-sided and its terminal pairing is not resolved here, so it is carried at
        % epsilon = 1, which is the conservative reading.                                     (G7)
        misp.dGpTcAnchEff = min(misp.dGpTcAnchEff, dGpTcil);
        % AND IT ENTERS THE ENSEMBLE SUM, which it previously did not.
        % EnsembleDG used to run above this block, so a bridge that beat every ungapped register
        % updated dGpTcAnchEff and never reached fluxAnch. The ensemble then came out LESS
        % negative than the per-site minimum it is supposed to sit at or below, which is
        % thermodynamically impossible for a Boltzmann sum that contains that site, and it showed
        % up in CalibrateAggregateFlux as a positive ShiftAll.
        % epsilon = 1 for a bridge, so the effective free energy is the free energy and the weight
        % is exp(-dG/RT) with no extension term. It goes into BOTH sums because the rejection acts
        % on it: the test two lines below is exactly the rejection fluxAnchRun is defined over.
        % THE COUNTS MOVE WITH THE SUMS. EnsembleDG returns 0 both for an empty set and for a set
        % whose free energy really is zero, and its docstring says to read it with the
        % accompanying count. Adding weight without adding the count would make a bridge-only
        % ensemble read as empty.
        scC.fluxAnch    = scC.fluxAnch    + exp(-dGpTcil/RTm);
        scC.fluxAnchRun = scC.fluxAnchRun + exp(-dGpTcil/RTm);
        scC.nanch       = scC.nanch    + 1;
        scC.nanchRun    = scC.nanchRun + 1;
        if dGpTcil < dGpTcmin
            check = false;
        end
    end
end
misp.dGpTcAnchEns    = EnsembleDG(scC.fluxAnch,    RTm);
misp.dGpTcAnchEnsRun = EnsembleDG(scC.fluxAnchRun, RTm);
misp.nAnchC          = scC.nanch;
misp.nAnchRunC       = scC.nanchRun;

%% ---- scan T: the opposite strand ---------------------------------------------------------------
if check == true
    dGpTmin = dGMispOffTarget(dGp, phi, pcon, T, tann);
    [scT, alstructT] = SeedScan(primer, seqrcomplement(lintemp), monv, biv, dNTP, tann, ...
        dGpTmin, SEEDLEN, minRun, SURVEYRUN, endTol, endWin, NUMALN, [], [], useEps, ryc);
    dGpT = scT.dGbest;
    misp.ScannedT     = true;
    misp.dGpTAnch     = scT.dGanch;
    misp.dGpTAnchEff  = scT.dGanchEff;
    % As on the cognate strand, the ensemble sums are formed after the bridge block below.
    misp.dGpTEps      = scT.eps;
    misp.dGpTEpsPair  = scT.epsPair;
    misp.dGpTEnd   = scT.dGend;
    misp.dGpTPos   = scT.pos;
    misp.dGpTEndRn = scT.endrun;
    misp.dGpTAlign = scT.align;
    misp.nSurvT    = scT.nsurv;
    misp.dGpTTrun   = scT.dGtrun;                                                             % (G9)
    misp.dGpTTrunRn = scT.trun;                                                               % (G9)
    misp.dGpTTrunMin = scT.dGtrunMin;                                                         % (G9)
    if scT.offend
        check = false;
    end
    if scT.rych                                                                               % (G9)
        check = false;
    end
    if check == true
        [dGpTil, ilfound] = InternalLoopBridge(alstructT, primer, seqrcomplement(lintemp), ...
            monv, biv, dNTP, tann, endTol);
        if ilfound
            dGpT = min(dGpT, dGpTil);
            % As on the cognate strand: the accumulators carry every bridge that was found, so the
            % reported margin describes the worst structure rather than only the rejected ones.
            misp.dGpTAnch    = min(misp.dGpTAnch, dGpTil);
            misp.dGpTAnchEff = min(misp.dGpTAnchEff, dGpTil);
            scT.fluxAnch    = scT.fluxAnch    + exp(-dGpTil/RTm);
            scT.fluxAnchRun = scT.fluxAnchRun + exp(-dGpTil/RTm);
            scT.nanch       = scT.nanch    + 1;
            scT.nanchRun    = scT.nanchRun + 1;
            if dGpTil < dGpTmin
                check = false;
            end
        end
    end
    misp.dGpTAnchEns    = EnsembleDG(scT.fluxAnch,    RTm);
    misp.dGpTAnchEnsRun = EnsembleDG(scT.fluxAnchRun, RTm);
    misp.nAnchT         = scT.nanch;
    misp.nAnchRunT      = scT.nanchRun;
end
end

%% ____________________________________________________________________________________________________________
function g = EnsembleDG(flux, RT)
%% ENSEMBLEDG  Ensemble free energy of a set of sites from their summed Boltzmann weights
% g = -RT*log(sum_k exp(-dG_k/RT)). A single dominant site returns that site's free energy; M sites
% degenerate at dG return dG - RT*log(M). An empty set returns 0, which is the same sentinel the
% anchored accumulators use for "no qualifying register was found" and must be read with the
% accompanying count rather than as a free energy.
if ~isfinite(flux)||(flux <= 0)
    g = 0; return
end
g = -RT*log(flux);
end

%% ____________________________________________________________________________________________________________
function [sc, alstruct] = SeedScan(primer, template, monv, biv, dNTP, tann, thr, ...
    seedLen, minRun, surveyRun, endTol, endWin, numaln, exclLo, exclHi, useEps, ryc)
% Every register of the primer against one strand that carries an exact seedLen-mer identity, scored
% with dGMismatch.
%
% sc.dGbest  most stable register found, kcal/mol; 0 when none is scored
% sc.dGanch  most stable register whose duplex reaches within endTol of the primer 3' terminus;
%            0 when none exists, which reads as full headroom                                  (F1)
% sc.dGend   stacking-only dG of the last endWin primer bases of THAT register                 (F2)
% sc.pos     start coordinate of that register on the strand scanned
% sc.endrun  contiguous identities ending at the duplex 3' end of that register
% sc.align   three rows of that register joined with slashes, as a string scalar
% sc.nsurv   registers scored
% sc.nanch   distinct 3'-anchored template spans, counted once per physical duplex
% sc.offend  true when some register is more stable than thr AND carries a run of minRun identities
%            AND places the primer 3' terminus within endTol of the duplex end
% sc.trun    length of the contiguous paired run carrying the 3' terminus of the winning register (G9)
% sc.dGtrun  free energy of that run at the Rychlik reference state                            (G9)
% sc.dGtrunMin  most stable such run over ALL anchored registers of this scan                  (G9)
% sc.rych    true when any anchored register's terminal run beats ryc.Max                      (G9)
% alstruct   localalign-shaped struct holding the numaln most stable registers, for InternalLoopBridge
%
% Frame: the primer is IDENTICAL to the template over the footprint, so the base opposing primer(i)
% is the complement of the footprint at i, and dGMismatch takes seq2 = seqcomplement(footprint) with
% no reversal. This is the same convention dGmaxMismatch applies to a localalign result.

% Declared first, before any other statement, which is where MATLAB expects it.
persistent surveywarned

n = strlength(primer); L = strlength(template);
D = [template template];
seen = false(1, L);
% Distinct template spans, not distinct full-primer registers.
% seen() stops one register being scored twice, which is a different question from whether two
% registers describe the same duplex. BestSubRegister trims both edges, so registers starting one
% or more bases apart can trim to the same physical span on the template whenever an accidental
% upstream identity or a terminal repeat shifts the window. Keying the additive accumulators on
% the untrimmed start therefore put one duplex into the count and into the Boltzmann sum k times,
% which shifts an ensemble by -RT*log(k) for a multiplicity that does not exist. These two maps
% key on the trimmed span in template coordinates, so one physical duplex enters once.
% The anchored winner and the gate winner can trim to different left edges, so they are tracked
% separately rather than through one key.
% Each map stores the most negative effective free energy (double) seen at that template span.
% If a later seed offset trims to the same template span with a more stable effective free energy
% (for instance, lAnch = n with a matched terminus after an earlier register trimmed to the same
% template coordinates with lAnch = n - 1), the weaker Boltzmann weight is replaced by the stronger
% one while keeping the distinct-span count at one.
seenAnchSpan = containers.Map('KeyType', 'char', 'ValueType', 'double');
seenGateSpan = containers.Map('KeyType', 'char', 'ValueType', 'double');

% One shot per MATLAB session. The engine calls this function a few thousand times, so a plain
% warning would bury the console.                                                              (F4)
% SURVEYRUN ABOVE minRun BREAKS THE SCREEN, SO IT STOPS RATHER THAN WARNS.                    (H9)
% surveyRun is handed to BestSubRegister as its minrun and skips any sub-range whose longest
% contiguous paired run falls short of it. The rejection rule downstream fires on
% sub.runRight >= minRun, and sub.runRight can only ever come from a sub-range that already
% survived this filter. Set surveyRun above minRun and every sub-range whose longest run is
% between the two is discarded before it is costed, so a genuinely rejection-eligible mispriming
% site returns no qualifying run and passes. That is a silent loss of specificity, which is why it
% is an error and not advice.
if surveyRun > minRun
    error('MisprimerMismatch:UnsafeSurveyRun', ...
        ['SURVEYRUN = %d exceeds minRun = %d. Sub-registers carrying a run between %d and %d ' ...
        'would be discarded before being scored, yet those are exactly the registers the ' ...
        'rejection rule acts on, so off-target sites would pass unseen. Keep SURVEYRUN <= ' ...
        'minRun.'], surveyRun, minRun, minRun, surveyRun);
end
% The remaining notice is about cost, not correctness. Every register reaching the sweep carries an
% exact seedLen-mer by construction, so the sub-range that WINS always clears a threshold of
% seedLen or less and the filter cannot change the result. It does still fire: on a 21-mer with one
% exact 5-mer seed it skips about 38 per cent of sub-ranges at surveyRun = 3 and about 75 per cent
% at surveyRun = 5. So it is a live cost filter whose only inert property is its effect on the
% answer, and raising it to minRun is free speed.
% A SURVEY FLOOR BELOW SEEDLEN IS A DELIBERATE SETTING, SO THIS REPORTS RATHER THAN SCOLDS.
% The earlier text here said a floor below SEEDLEN "cannot change the result" and advised raising
% it back to minRun. That was written before the two ensemble sums existed. It is still true of
% dGanchEff, whose winning sub-register always carries at least SEEDLEN identities, and it is
% false of fluxAnch, which is a sum over every scored sub-register and gains a term for each one
% the lower floor admits. Widening the survey is the whole purpose of the setting.
if isempty(surveywarned), surveywarned = false; end
if (surveyRun < seedLen) && ~surveywarned
    surveywarned = true;
    warning('MisprimerMismatch:WideSurvey', ...
        ['Survey floor %d is below SEEDLEN %d, so the ensemble sums run over sub-registers that ' ...
        'cannot win the single-site minimum and cannot trigger the rejection. That is what makes ' ...
        'fluxAnch and fluxAnchRun differ, and it is what a calibration run needs. The sweep is ' ...
        'roughly four times wider per nucleotide below minRun = %d. Reported once per session.'], ...
        surveyRun, seedLen, minRun);
end

sc = struct('dGbest', 0, 'dGanch', 0, 'dGanchEff', 0, 'eps', 1, 'epsPair', "", ...
    'dGend', NaN, 'pos', NaN, 'endrun', NaN, ...
    'align', "", 'nsurv', 0, 'nanch', 0, 'offend', false, ...
    'trun', 0, 'dGtrun', 0, 'dGtrunMin', 0, 'rych', false, ...
    'fluxAnch', 0, 'fluxAnchRun', 0, 'nanchRun', 0);                                          % (G9)
RT = 1.987E-3*(273.15 + tann);                                                              % (G7)
cdg = []; cst = []; cf = []; cl = []; ca1 = {}; ca2 = {}; ca3 = {};

for o = 1:(n - seedLen + 1)
    idx = strfind(D, primer(o:o + seedLen - 1));
    for q = idx
        % No "q > L" guard. A seed occurrence lying past the seam still implies a real circular
        % register, mod folds its start back into [1, L], and seen() stops the double count.  (F3)
        st = mod(q - o, L) + 1;
        if seen(st), continue; end
        seen(st) = true;

        if ~isempty(exclLo)
            cov = mod((st:st + n - 1) - 1, L) + 1;   % circular footprint
            if any(cov >= exclLo & cov <= exclHi), continue; end
        end

        foot = D(st:st + n - 1);
        ispair = (primer == foot);
        if ~any(ispair), continue; end

        % SUB-REGISTER SWEEP, replacing a single trim to the outermost identity.               (G3)
        % The old code set f and l to the outermost identities and scored everything between them
        % as one duplex. On a 7249 nt template the probability of at least one accidental identity
        % in the first eleven primer positions is 1 - (3/4)^11 = 0.958, so a stray pair dragged the
        % register out on essentially every hit and the intervening mismatches were charged as an
        % interior loop; an internal loop of 26 costs +6.79 at annealing temperature. It also
        % corrupted the anchoring test, which read l off those same ballooned bounds. See the note
        % at the top of BestSubRegister.m for the measured example.
        sub = BestSubRegister(primer, foot, ispair, surveyRun, monv, biv, dNTP, tann, false, minRun);
        if ~isfinite(sub.dGfree), continue; end

        dg  = sub.dGfree;
        f   = sub.iFree;   l = sub.jFree;
        seg = primer(f:l); fseg = foot(f:l); ip = ispair(f:l);

        sc.nsurv = sc.nsurv + 1;
        if dg < sc.dGbest, sc.dGbest = dg; end

        % The duplex must REACH the primer 3' terminus to be extendable, so the anchoring test uses
        % the outermost 3' pair, sub.jRight, and NOT the free-edge optimum above, whose right edge
        % may have moved inward. lAnch is that pair's coordinate in the primer, so lAnch >= n -
        % endTol is the terminus test. Testing the coordinate in the TARGET instead is a no-op: on
        % a multi-kilobase template essentially every hit satisfies it.
        lAnch = sub.jRight;
        anchored = isfinite(sub.dGright)&&(lAnch >= n - endTol);
        dGeffA = Inf;
        if anchored
            ipA = ispair(sub.iRight:lAnch);
            anchKey = sprintf('%d:%d', mod(st + sub.iRight - 2, L) + 1, mod(st + lAnch - 2, L) + 1);
            % EXTENSION EFFICIENCY.                                                            (G7)
            % Occupancy says how often a primer sits at this site. It says nothing about whether Taq
            % will add a nucleotide once it is there, and for a mispaired 3' terminus the answer is
            % usually no: 1E-3 to below 1E-6 of a matched terminus, measured in Huang MM, Arnheim N,
            % Goodman MF, Nucleic Acids Res 1992;20:4567, doi:10.1093/nar/20.17.4567. Those authors
            % also showed Taq BINDS matched and mispaired 3' ends with similar affinity, so epsilon
            % belongs on the catalytic step and must not be folded into a binding constant.
            %
            % The rate of productive priming goes as theta*epsilon, and thr is the free energy at
            % which theta alone reaches phi. Demanding theta*epsilon <= phi is therefore the same as
            % demanding dG >= thr + RT*log(epsilon), which is written here as a penalty on the site
            % instead: dGeff = dG - RT*log(epsilon). With a matched terminus epsilon is 1 and
            % nothing moves; with a C(primer).T terminus, epsilon = 2.0E-2 in the shipped table, it
            % costs 2.58 kcal/mol at 59 C; with a G.G terminus, epsilon = 1.0E-6, it costs 9.12.
            % Continuous 3'-terminal duplex, scored at the Rychlik reference state rather than at
            % tann, because the -11 figure is a property of that reference state and not of this
            % reaction. Evaluated over every anchored register, not just the winner, because the
            % gate asks whether any extendable site carries a run that stable.
            % The run ends at position n rather than at lAnch. Rychlik (1995, Biotechniques
            % 18:84-90) measured primers whose 5' termini were mismatched and whose 3' portion
            % formed a continuous duplex with the target, so the -11 kcal/mol figure describes a
            % run that carries the 3' terminus. Ending it at lAnch scores an interior run against
            % a threshold defined on a terminal one, and eliminates a primer whose 3' base is
            % mispaired and therefore cannot prime from this site whatever sits behind it.
            % A mispaired position n leaves TerminalRundG fewer than two positions, so it returns
            % zero and the gate stays silent at that locus.
            dGtrA = TerminalRundG(primer, foot, sub.iTrunRight, n, ryc);
            if dGtrA < sc.dGtrunMin, sc.dGtrunMin = dGtrA; end
            if ryc.On&&(dGtrA <= ryc.Max), sc.rych = true; end

            % epsA and dGtrA are properties of the locus, not of the winning sub-register:
            % TerminalEps reads the primer 3' terminus against the footprint, and TerminalRundG
            % reads the contiguous run carrying that same terminus. Both winners therefore share
            % them, and only the duplex free energy differs between the two.
            [epsA, pairA] = TerminalEps(primer, foot, n, useEps);
            dGeffA = sub.dGright - RT*log(epsA);
            % THE GATE'S OWN WINNER AT THIS LOCUS.
            % With the survey floor below minRun a locus can have its most stable sub-register
            % carry a run too short to reject on. BestSubRegister reports that one as dGright and
            % the best REJECTION-ELIGIBLE one as dGgate, and the two are the same object whenever
            % the two floors agree. Everything the rejection and the margin consume is taken from
            % the gate winner, so widening the survey cannot hide a locus from the gate.
            gated  = isfinite(sub.dGgate);
            dGeffG = Inf;
            if gated, dGeffG = sub.dGgate - RT*log(epsA); end

            % AGGREGATE OFF-TARGET FLUX. MEASUREMENT ONLY: NOTHING GATES ON THIS.
            %
            % The rejection below asks about the single most stable anchored register. The quantity
            % that actually competes with the cognate site is the SUM over every site a primer can
            % prime from, because a candidate carrying many comparable off-target registers
            % initiates more spurious product per cycle than one carrying a single slightly worse
            % register, and the present gate cannot tell them apart.
            %
            % At phi = 1E-4 every off-target site sits at an occupancy of order 1E-4, far inside the
            % linear regime where theta_i = K_i*[p]. There the ensemble association constant is
            % literally the sum of the individual ones, so adding Boltzmann factors is exact rather
            % than an approximation, and -RT*log of the sum is the ensemble free energy. One
            % dominant site reproduces the minimum; M degenerate sites shift it by RT*log(M).
            %
            % TWO SUMS, BECAUSE THE POPULATION IS A CHOICE AND NOT AN OBVIOUS ONE. fluxAnch runs
            % over every anchored register, which is the physics. fluxAnchRun runs only over those
            % carrying minRun identities, which is the population the current rejection acts on.
            % They answer different questions and the calibration run decides which one a gate
            % should use.
            %
            % Overflow is not a concern: dGeffA would have to reach about -460 kcal/mol before
            % exp(-dGeffA/RT) left double range, so the sum is taken directly rather than through a
            % log-sum-exp.
            % One physical duplex contributes one weight. newAnch was decided with the count
            % above, so the sum and the count it has to be read with stay in step.
            if ~isKey(seenAnchSpan, anchKey)
                seenAnchSpan(anchKey) = dGeffA;
                sc.nanch    = sc.nanch + 1;
                sc.fluxAnch = sc.fluxAnch + exp(-dGeffA/RT);
            elseif dGeffA < seenAnchSpan(anchKey)
                sc.fluxAnch = sc.fluxAnch - exp(-seenAnchSpan(anchKey)/RT) + exp(-dGeffA/RT);
                seenAnchSpan(anchKey) = dGeffA;
            end
            if gated
                gateKey = sprintf('%d:%d', mod(st + sub.iGate - 2, L) + 1, mod(st + lAnch - 2, L) + 1);
                if ~isKey(seenGateSpan, gateKey)
                    seenGateSpan(gateKey) = dGeffG;
                    sc.nanchRun    = sc.nanchRun + 1;
                    sc.fluxAnchRun = sc.fluxAnchRun + exp(-dGeffG/RT);
                elseif dGeffG < seenGateSpan(gateKey)
                    sc.fluxAnchRun = sc.fluxAnchRun - exp(-seenGateSpan(gateKey)/RT) + exp(-dGeffG/RT);
                    seenGateSpan(gateKey) = dGeffG;
                end
            end

            % THE WORST REGISTER IS THE GATE'S WORST REGISTER.
            % dGanchEff is what the engine turns into fTcAnchMargin and hands to Score, and the
            % rejection below acts on the same quantity. Both read the gate winner, so the margin
            % and the gate describe one object whatever the survey floor is set to. The alignment
            % and the trailing-run figure describe that same sub-register.
            if gated && (dGeffG < sc.dGanchEff)
                ipG = ispair(sub.iGate:lAnch);
                sc.dGanchEff = dGeffG;
                sc.dGanch = sub.dGgate;
                sc.eps    = epsA;
                sc.epsPair = pairA;
                sc.pos    = st;
                sc.endrun = TrailRun(ipG);
                sc.trun   = sub.trunRight;                                                  % (G9)
                sc.dGtrun = dGtrA;                                                          % (G9)
                sc.dGend  = EndWindowdG(primer, foot, lAnch, n, endWin, monv, biv, dNTP, tann);
                sc.align  = string([primer(sub.iGate:lAnch), '/', PairSymbols(ipG), '/', ...
                    foot(sub.iGate:lAnch)]);
            end
        end

        % Rejection needs all three conditions together, and all three must describe the SAME
        % structure. Under the single trim the register dG and the anchored dG were necessarily the
        % same number, so this read dg; with a free left edge they diverge, and the anchored one is
        % what an extension event would see.                                                   (G3)
        % The rejection now reads the EXTENSION-WEIGHTED free energy, so a site that is occupied
        % often enough to matter but ends on a terminus Taq will not extend no longer rejects the
        % candidate. Set useEps false to recover the occupancy-only test exactly.              (G7)
        if anchored&&gated&&(dGeffG < thr)
            sc.offend = true;
        end

        cdg(end+1) = dg;              %#ok<AGROW>
        cst(end+1) = st;              %#ok<AGROW>
        cf(end+1)  = f;               %#ok<AGROW>
        cl(end+1)  = l;               %#ok<AGROW>
        ca1{end+1} = seg;             %#ok<AGROW>
        ca2{end+1} = PairSymbols(ip); %#ok<AGROW>
        ca3{end+1} = fseg;            %#ok<AGROW>
    end
end

alstruct = struct();
alstruct.Alignment = {};
alstruct.Start = zeros(0, 2);
alstruct.Stop  = zeros(0, 2);
if ~isempty(cdg)
    [~, ord] = sort(cdg);
    keep = ord(1:min(numaln, numel(ord)));
    alstruct.Alignment = arrayfun(@(k) [ca1{k}; ca2{k}; ca3{k}], keep, 'UniformOutput', false);
    alstruct.Start = [cf(keep)', (cst(keep)' + cf(keep)' - 1)];
    alstruct.Stop  = [cl(keep)', (cst(keep)' + cl(keep)' - 1)];
end
end

%% ____________________________________________________________________________________________________________
function [e, pairTxt] = TerminalEps(primer, foot, n, useEps)
% Relative extension efficiency from this register's 3' terminus.                              (G7)
%
% FRAME. The primer is IDENTICAL to the strand being scanned over the footprint, so the base
% PHYSICALLY OPPOSING primer position k is the COMPLEMENT of foot(k), not foot(k) itself. Passing
% the raw footprint inverts every pairing test and returns 1 for every mismatch, which looks like a
% working function and silently disables the whole correction.
if ~useEps
    e = 1; pairTxt = ""; return
end
% THREE POSITIONS, AND THE REASON IS THE REGISTER, NOT THE TABLE. Rejali et al. 2018 measured
% extension out to ten bases from the 3' terminus and ExtensionEfficiency will score every one of
% them, but an anchored register here is only guaranteed to reach n - MispEndTol, and a run of
% MispMinRun contiguous identities need not sit at the 3' end. Past the third position the primer
% base passed may lie outside the helix at this locus, and a mismatch outside the helix is not the
% object that paper measured. Widening k0 is therefore a change to what counts as a register, not
% a change to this line.
k0 = max(1, n - 2);
pB = primer(k0:n);
tB = char(seqcomplement(foot(k0:n)));
[e, einfo] = ExtensionEfficiency(pB, tB);
% The string names the mismatch that SET epsilon, with its distance from the terminus appended
% when the terminus is not where it sits. A matched terminus carrying a mismatch one base in used
% to report the matched pair against an epsilon of exactly 1; it now reports that mismatch and the
% measured value, and the suffix is what separates the two cases in a log.
if einfo.GoverningPos > 0
    pairTxt = string(sprintf('%s@%d', einfo.GoverningPair, einfo.GoverningPos));
else
    pairTxt = string(einfo.GoverningPair);
end
end

%% ____________________________________________________________________________________________________________
function dgT = TerminalRundG(primer, foot, iTrun, l, ryc)
% Free energy of the CONTIGUOUS run of Watson-Crick pairs carrying the duplex 3' terminus, scored
% at the Rychlik reference state (ryc.T, ryc.Na) rather than at the reaction temperature.       (G9)
%
% Terminal and initiation rows are KEPT here, unlike EndWindowdG below. What Rychlik measured is
% the stability of a real, separable duplex between two molecules, so it carries helix initiation;
% the stage 6 window is a stacking profile and does not. Mixing those two conventions is exactly
% the error that put dG3min on the wrong scale, so the two are written out separately and neither
% is reused for the other.
dgT = 0;
if ~isfinite(iTrun)||~isfinite(l), return; end
if (l - iTrun + 1) < 2, return; end          % one pair is not a nearest-neighbour step
segT  = primer(iTrun:l);
footT = foot(iTrun:l);
try
    dgT = dGMismatch(segT, char(seqcomplement(footT)), ryc.Na, 0, 0, ryc.T, false, true, true);
catch
    dgT = 0; return                          % a failure must never fabricate a rejection
end
if ~isfinite(dgT), dgT = 0; end
end

%% ____________________________________________________________________________________________________________
function dgE = EndWindowdG(primer, foot, l, n, endWin, monv, biv, dNTP, tann)
% Stacking-only dG of the 3'-terminal window of one register, initiation and terminal terms
% suppressed, so it sits in the same currency as the Rychlik cutoff used on the free primer.   (F2)
w = max(1, n - endWin + 1);
if l <= w
    dgE = 0; return          % fewer than two positions inside the window, so no stack to sum
end
segE  = primer(w:l);
footE = foot(w:l);
% This value is reported, never gated on, so a failure here must not take the scan down with it.
try
    dgE = dGMismatch(segE, char(seqcomplement(footE)), monv, biv, dNTP, tann, false, false, false);
catch
    dgE = NaN; return
end
if ~isfinite(dgE), dgE = 0; end
end

%% ____________________________________________________________________________________________________________
function runmax = MaxRun(mask)
runmax = 0; runlen = 0;
for q = 1:numel(mask)
    if mask(q)
        runlen = runlen + 1;
        if runlen > runmax, runmax = runlen; end
    else
        runlen = 0;
    end
end
end

%% ____________________________________________________________________________________________________________
function r = TrailRun(mask)
% Contiguous identities ending at the last position of the mask, which by construction of f and l is
% the 3' end of the duplex. A five-base run sitting at positions 3 to 7 with one lone identity at the
% terminus gives r = 1 and is a very different object from a five-base run at the terminus itself.
r = 0;
for q = numel(mask):-1:1
    if mask(q), r = r + 1; else, break; end
end
end

%% ____________________________________________________________________________________________________________
function sym = PairSymbols(mask)
sym = repmat(' ', 1, numel(mask));
sym(mask) = '|';
end

%% ____________________________________________________________________________________________________________
function misp = EmptyMisp()
% Every field scalar or string scalar, so struct2table takes it without unpacking.             (F6)
misp = struct('dGpTcAnch', 0, 'dGpTcAnchEff', 0, 'dGpTcEps', 1, 'dGpTcEpsPair', "", ...
    'dGpTcEnd', NaN, 'dGpTcPos', NaN, 'dGpTcEndRn', NaN, ...
    'dGpTcAlign', "", 'nSurvC', 0, 'nAnchC', 0, ...
    'dGpTcTrun', 0, 'dGpTcTrunRn', 0, 'dGpTcTrunMin', 0, ...                                  % (G9)
    'dGpTAnch', 0, 'dGpTAnchEff', 0, 'dGpTEps', 1, 'dGpTEpsPair', "", ...
    'dGpTEnd', NaN, 'dGpTPos', NaN, 'dGpTEndRn', NaN, ...
    'dGpTAlign', "", 'nSurvT', 0, 'nAnchT', 0, ...
    'dGpTTrun', 0, 'dGpTTrunRn', 0, 'dGpTTrunMin', 0, ...                                     % (G9)
    ... % Aggregate off-target flux, MEASUREMENT ONLY. Nothing in this file or in the engine gates
    ... % on these. Ens is the ensemble free energy over every anchored register, EnsRun the same
    ... % sum restricted to registers carrying minRun identities. Compare either against the
    ... % corresponding AnchEff to see how far the ensemble sits from the single worst site.
    'dGpTcAnchEns', 0, 'dGpTcAnchEnsRun', 0, 'nAnchRunC', 0, ...
    'dGpTAnchEns', 0, 'dGpTAnchEnsRun', 0, 'nAnchRunT', 0, ...
    ... % WHICH STRANDS WERE ACTUALLY LOOKED AT. Every count and every ensemble figure above uses
    ... % 0 for "no qualifying register was found", and scan T is skipped whenever scan C already
    ... % rejected the primer. Without these two flags a skipped scan is indistinguishable in the
    ... % log from a scan that ran and came back clean, which is the stronger of the two claims.
    'ScannedC', false, 'ScannedT', false);
end
