function [cfg, Parameters] = PrimerSearchParameters(input)
%% PRIMERSEARCHPARAMETERS  Every constant, flag and threshold the aPCR primer search engines run on
% Written by Vishwesh Rai, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 12-11-2026 02:16 IST
%
% Revision history is NOT kept here. See PrimerSearchParameters_<YYYYMMDD>.txt beside this file.
% ____________________________________________________________________________________________________________
% WHAT THIS FILE IS
%
% One place where every number the search depends on is set, and the only place where the reason
% for each number is written down. The engines hold the algorithm. This file holds the settings
% the algorithm is run with, so a reader chasing a result opens one file to find out what was
% asked for and a different file to find out how it was computed.
%
% Three kinds of thing live here.
%   CONSTANTS, thresholds and windows that a run is configured with.
%   FLAGS that switch a stage between gating and logging, or switch it off.
%   DERIVED SETTINGS that cannot be written as a literal because they depend on the reaction: the
%   fold salt, and the four settings the polymerase rewrites.
%
% HOW A CALLER USES IT
%
%       [cfg, Parameters] = PrimerSearchParameters(input);
%
%   cfg         a struct with one field per setting, named exactly as the setting is named in the
%               engines. An engine unpacks the fields it consumes into local variables in a single
%               block at the top, so the engine body reads in the names a user of the suite knows.
%   Parameters  the human-readable parameter sheet, section by section, with units and meanings.
%               It carries the user's inputs alongside the engine's configuration and it rides in
%               the results struct, so a run that found nothing still says what it was asked to
%               find.
%
% WHY THE ENZYME CONDITIONING AND THE FOLD SALT ARE HERE RATHER THAN IN THE ENGINE
%
% Both decide the VALUE of a setting, which makes them parameter determination and not algorithm.
% A proofreading exonuclease rewrites four settings and changes the shape of three gates. The fold
% salt is derived from the buffer unless it is overridden. Keeping both here means the parameter
% sheet is assembled after everything that can move a setting has moved it, so the sheet describes
% the configuration that actually executed rather than the one that was typed.
%
% THE SHEET CHECKS ITSELF, TWICE
%
% AssertParamSheetComplete compares the workspace at the moment the sheet is built against the
% names on the sheet, so a setting added here and never written into ParamSpec is caught. At that
% point the workspace holds the settings and nothing else, which is what makes the comparison
% meaningful. AssertConfigMatchesSheet then compares the returned struct against the same sheet in
% both directions, so a setting on the sheet but missing from cfg, or in cfg but missing from the
% sheet, is caught as well. A table of settings that quietly omits a setting is worse than no
% table, because it invites the reader to trust it.
%
% ____________________________________________________________________________________________________________
% REFERENCES FOR THE VALUES SET BELOW
%
%   SantaLucia J Jr, Hicks D. The thermodynamics of DNA structural motifs. Annu Rev Biophys Biomol
%       Struct 2004;33:415-40.  doi:10.1146/annurev.biophys.32.110601.141800
%   Owczarzy R, You Y, Moreira BG, et al. Effects of sodium ions on DNA duplex oligomers.
%       Biochemistry 2004;43:3537-54.  doi:10.1021/bi034621r
%   von Ahsen N, Wittwer CT, Schutz E. Oligonucleotide melting temperatures under PCR conditions.
%       Clin Chem 2001;47:1956-61.  PMID 11673362
%   Sanchez JA, Pierce KE, Rice JE, Wangh LJ. Linear-after-the-exponential (LATE)-PCR.
%       Proc Natl Acad Sci USA 2004;101:1933-8.  doi:10.1073/pnas.0305476101
%   Pierce KE, Sanchez JA, Rice JE, Wangh LJ. LATE-PCR: primer design criteria for high yields of
%       specific single-stranded DNA. Proc Natl Acad Sci USA 2005;102:8609-14.
%       doi:10.1073/pnas.0501946102
%   Huang MM, Arnheim N, Goodman MF. Extension of base mispairs by Taq DNA polymerase.
%       Nucleic Acids Res 1992;20:4567-73.  doi:10.1093/nar/20.17.4567
%   Kwok S, Kellogg DE, McKinney N, et al. Effects of primer-template mismatches on the polymerase
%       chain reaction. Nucleic Acids Res 1990;18:999-1005.  doi:10.1093/nar/18.4.999
%   Rychlik W. Priming efficiency in PCR. Biotechniques 1995;18:84-90.  PMID 7702859
%   Rychlik W. Selection of primers for polymerase chain reaction. Mol Biotechnol 1995;3:129-34.
%       doi:10.1007/BF02789108
%   Ye J, Coulouris G, Zaretskaya I, Cutcutache I, Rozen S, Madden TL. Primer-BLAST: a tool to
%       design target-specific primers for polymerase chain reaction. BMC Bioinformatics
%       2012;13:134.  doi:10.1186/1471-2105-13-134
%   Wallace RB, Shaffer J, Murphy RF, et al. Hybridization of synthetic oligodeoxyribonucleotides
%       to phi chi 174 DNA. Nucleic Acids Res 1979;6:3543-57.  doi:10.1093/nar/6.11.3543
% ____________________________________________________________________________________________________________

%% ---- SEARCH GEOMETRY --------------------------------------------------------------------------
% Both windows only ever grow the amplicon. The requested scaffold is always fully contained, so a
% candidate can start earlier than the scaffold and end later than it, never inside it.
fFlex = 20;                     % Scaffold extension allowed at the 5' end, nt. Sets the forward
                                % search width: fFlex + 1 start positions.
rFlex = 20;                     % Same at the 3' end, for the reverse primer.
fLmin = 18; fLmax = 24;         % Forward primer length
rLmin = 20; rLmax = 28;         % Reverse primer length

%% ---- MELTING TEMPERATURE ----------------------------------------------------------------------
Tmin = 58.5; Tmax = 70;         % Primer Tm window, degrees C.
                                % TMAX IS NOT A PIERCE NUMBER. Converting Pierce's own Series C
                                % primers into the engine's magnesium-corrected units, with a
                                % magnesium lift of 5.08 degrees at their buffer of 70 mM
                                % monovalent, 3.5 mM Mg and 1.6 mM dNTP:
                                %     gap 23   TmX 52.1   TmL 58.1   early slope 2.3
                                %     gap 20   TmX 55.1   TmL 61.1   early slope 4.2
                                %     gap 17   TmX 59.1   TmL 64.1   early slope 7.6
                                %     gap 13   TmX 64.1   TmL 69.1   early slope 11.3
                                % Their best design needs a limiting primer at 69.1 in these units.
                                % A ceiling of 65 rejects it outright. 72 spans Pierce's range with
                                % room for the limiting primer to sit dTrfmin above a qualifying
                                % excess primer.
                                % WATCH rLmax TOO. The limiting primer has to reach TmX + dTrfmin,
                                % so if the reverse window cannot get there, lengthen it before
                                % loosening TmAmaxGap. Loosening the gate buys linear-phase rate
                                % that Pierce measured you losing; a longer primer does not.
                                % Stated as an explicit ceiling rather than as a floor plus a
                                % width, so that raising the lower bound does not silently raise
                                % the upper one by the same amount.
dTrfmin = 5.0; dTrfmax = 8;     % LATE-PCR offset, Tm(limiting) - Tm(excess), degrees
                                % Pierce et al. report the highest reaction efficiency AND
                                % specificity at TmL - TmX of about +5 C, measuring 0.85, 0.89,
                                % 0.89, 0.97 and 1.00 at offsets of -3, 0, +3, +5 and +7. READ THAT
                                % SERIES CAREFULLY, BECAUSE IT DOES NOT SAY WHAT A WINDOW ENDING AT
                                % 6.5 WOULD LIKE IT TO SAY. Efficiency is flat from 0 to +3 and then
                                % RISES monotonically, 0.89 to 0.97 to 1.00, with +7 the highest
                                % offset they measured and the best one. There is no measured
                                % plateau and no measured turning point. On efficiency alone the
                                % data argue for a ceiling at or above +7.
                                % THE FLOOR IS THE PUBLISHED CONSTRAINT. THE CEILING IS OURS AND IT
                                % COSTS SOMETHING. Pierce et al. require the limiting primer to melt
                                % ABOVE the excess primer, and the floor here sits just below the
                                % point where their efficiency starts climbing. The ceiling trades
                                % a few per cent of amplification efficiency for a tighter spread of
                                % offsets, and the trade is a judgement rather than a reading of
                                % their figure. Anyone widening it towards +7 is moving TOWARDS the
                                % published optimum.
                                % Their measurements were made at a 20:1 primer ratio and the
                                % optimum shifts with that ratio, which is the honest reason not to
                                % treat +7 as a target at a different ratio.
                                % A NARROW WINDOW COSTS CANDIDATES, AND THAT IS THE POINT.
                                % Every degree of width here admits pairs whose offset sits further
                                % from the measured optimum. The ranking term already pulls towards
                                % the window centre, so widening the gate and letting Score sort it
                                % out is the alternative; narrowing it instead makes the engine
                                % refuse rather than merely rank down.
dTsub = 5;                      % Degrees below Tm at which structure is evaluated
% AMPLICON COMPETITION IS RANKED, NOT GATED.
%
% TmA - TmX decides how the accumulating single strand competes with the excess primer for the
% template once the limiting primer is spent. Pierce measured the consequence on their own
% amplicons and reported early linear-phase slopes at four gaps. Those four points are written
% down below as the measurement they are, and the ranking reads them directly.
%
% THE REASON THIS IS NOT A GATE. The measurement was made on Pierce's constructs, and whether the
% same relationship holds for a multi-kilobase amplicon is unestablished. Wetmur and Davidson give
% the renaturation rate constant as proportional to the square root of length at fixed complexity
% and inversely proportional to complexity, so for an amplicon, whose complexity is its own length,
% the per-molecule rate falls as L^-0.5. A 3940 nt product therefore re-anneals at roughly 0.28 of
% the rate of a 300 nt one at equal molar concentration, which makes the competition this quantity
% stands for WEAKER at long length and Pierce's thresholds correspondingly conservative. That is a
% direction, not a correction, and nothing here converts it into degrees. Excluding a pair on a
% threshold whose transferability is unknown throws away candidates for a reason that cannot be
% checked; ranking them keeps every candidate visible and lets the evidence decide later.
TmAPierceGap   = [13.0 17.0 20.0 23.0];    % TmA - TmX, degrees, Pierce's four tested values
TmAPierceSlope = [11.3  7.6  4.2  2.3];    % early linear-phase slope measured at each
                                % MONOVALENT CONVENTION on both sides. These apply to dTmAXP and
                                % never to dTmAX, which carries magnesium on both sides and runs
                                % about ten degrees hotter on an ordinary PCR buffer.
TmAoptGap = 13.0;               % The gap at which the penalty is zero. Pierce's best TESTED value
                                % rather than a fitted optimum: their curve is still climbing at the
                                % narrow end, so nothing is known about gaps below this and nothing
                                % is rewarded for reaching them.
TmAmaxGap = 20.0;               % ADVISORY ONLY. No longer gates anything. Retained because the
                                % startup feasibility check reports against it, and because it is
                                % the number a reader of Pierce would expect to find here.
UseTmAEquiv = true;             % false uses the published monovalent-only salt form in TmAmplicon,
                                % true substitutes the von Ahsen monovalent equivalent.
                                % THIS FLAG DOES NOT REACH THE GATE. The gate and the ranking both
                                % run on the monovalent pair TmAP and TmXP. The flag sets
                                % record.TmA, the number to quote as the amplicon's actual melting
                                % temperature in the reaction buffer, and true is the honest answer
                                % there because the buffer does contain magnesium.

%% ---- COMPOSITION ------------------------------------------------------------------------------
UseGCWindow = true;             % Enforce the global GC fraction window below.
                                % The window carries three arguments and they do not age alike.
                                % FIRST, before nearest-neighbour thermodynamics the GC fraction WAS
                                % the melting temperature: Tm = 2(A+T) + 4(G+C), the Wallace rule.
                                % At fixed length a GC window and a Tm window are one constraint
                                % written twice, and this suite computes Tm from the SantaLucia
                                % parameters at the true strand concentration, so that argument is
                                % spent. Of the candidates this window rejects on a 7 kb target,
                                % only about three per cent would have cleared the Tm window.
                                % SECOND, at FIXED Tm a GC-poor primer must be longer, and a longer
                                % duplex dilutes the cost of a single mismatch. Ensembles of 120
                                % primers each, all melting at 60.0 C, every interior position
                                % mismatched in turn, give an equilibrium selectivity K_perf/K_mm of
                                % 19.6 at 30 % GC (29.1 nt mean) rising monotonically to 52.2 at
                                % 70 % GC (15.9 nt). This is real content the Tm gate cannot see,
                                % and the mispriming scan covers it better by scoring mismatched
                                % duplexes against the actual template rather than against a
                                % composition proxy for a hypothetical one.
                                % THIRD, and in the opposite direction, GC-rich primers at matched
                                % Tm fold: in the same ensembles the fraction with a hairpin above
                                % 50 C runs 14 % at 40 % GC and 32 % at 70 % GC, and the fraction
                                % with a homodimer below -5 kcal/mol runs 27 % and 55 %. The dimer
                                % and hairpin stages measure both directly.
                                % Left ON by default because the floor's content is only partly
                                % covered. Set false to let Tm and the structure stages decide.
fGCmin = 0.40; fGCmax = 0.65;   % Global GC fraction window
MonoLim = 4; DiRepLim = 4;      % Longest tolerated homopolymer (nt) and dinucleotide repeat (units).
                                % NucleoRun rejects when the longest run EXCEEDS MonoLim, so 4
                                % admits AAAA and rejects AAAAA. Primer3 expresses the same rule as
                                % PRIMER_MAX_POLY_X, the maximum allowable mononucleotide repeat,
                                % and ships it at 5, so this is still one nucleotide stricter than
                                % the most widely used reference implementation.
                                % WHAT THE PREVIOUS VALUE OF 3 COST. On M13mp18 it removed 144 of
                                % 707 forward candidates and 205 of 1010 reverse candidates, and
                                % every one of those 349 carried a homopolymer of exactly 4. Not a
                                % single candidate on either strand was rejected for a run of 5 or
                                % longer, and not one was rejected on the dinucleotide repeat. The
                                % gate was doing one thing: removing any primer containing AAAA,
                                % CCCC, GGGG or TTTT, which is 20.3 per cent of the whole pool on a
                                % natural template.
                                % Reference: Untergasser A et al. Primer3, new capabilities and
                                % interfaces. Nucleic Acids Res 2012;40:e115.
                                % doi:10.1093/nar/gks596

%% ---- TERMINAL BASE AND SKEW RULES -------------------------------------------------------------
Force5primeGC        = false;   % Demand G or C at the 5' terminus of both primers
ForwardForce3primeAT = false;    % Demand A or T at the forward (excess) primer 3' terminus
                                % ON. THE CLAMP SURVIVES THE EXTENSION-EFFICIENCY WEIGHTING, AND
                                % THE REASON IS A DISTINCTION BETWEEN TWO DIFFERENT 3'-END
                                % FAILURES. Huang et al. measured how fast Taq extends a MISPAIRED
                                % terminus, and that measurement enters the mispriming stage as
                                % epsilon. What the clamp protects against is a terminus that is
                                % correctly paired and weakly stacked, which frays. Extension from
                                % a transiently unpaired 3'-OH is a rate, and every screen in the
                                % cascade integrates over bound states and is blind to it. The two
                                % are not the same quantity and the second does not subsume the
                                % first. The clamp also generalises to templates the mispriming
                                % scan has never seen.
                                % WHAT THE CLAMP COSTS, so the trade stays visible: it removes
                                % roughly two forward candidates in five on a mixed-composition
                                % target. By the geometric mean of epsilon over the three mispairs
                                % each base can make the ranking is G (6.5e-6), then A (1.1e-5),
                                % then C (1.4e-4), then T (1.5e-4), so a rule reading "A or T" does
                                % pair the second-best discriminating base with the worst. That is
                                % an argument for replacing the rule with a better one rather than
                                % for deleting it. The base is logged and scored under Wt3Base
                                % whether or not the clamp is enforced. If you turn it off, log it.
ForwardAvoid3primeT  = true;    % Narrow that clamp to A alone, following Kwok et al.
ReverseForce3primeGC = false;   % Demand G or C at the reverse (limiting) primer 3' terminus
                                % OFF, for the same reason, and because "C or G" pairs the best
                                % base with one carrying a 2e-2 exception at C.T that Huang
                                % singles out.
ForwardGCSkew        = false;   % Demand a GC-poorer 3' half than 5' half, forward primer
                                % OFF. This is a two-bin composition threshold standing in for
                                % Rychlik's internal stability profile, which asks that a primer be
                                % less tightly bound at its 3' end than at its 5' end so the
                                % polymerase sees a 3' terminus that breathes. That profile is
                                % computed directly in free-energy terms as dG3stack - dG5stack
                                % over the terminal windows and carried in the Score under WtProf,
                                % so the continuous quantity is already in the objective function.
                                % Reducing it to a hard split on GC counts discards information and
                                % then rejects on the remainder. Of the forward candidates the flag
                                % removes on a mixed-composition target, roughly two thirds have a
                                % melting temperature inside the search window, against three per
                                % cent for the GC window above. Logged always, enforced only if set
                                % true.
ReverseGCSkew        = false;   % Same on the reverse primer. Turning it on rejects a reverse
                                % primer whose 3' half is GC-richer than its 5' half, which for a
                                % LIMITING primer is desirable rather than a fault, since that
                                % primer has to hold its site while the excess primer competes for
                                % the same strand. drlocgc35 is computed and logged either way.

%% ---- 3'-END STABILITY AND HAIRPINS ------------------------------------------------------------
EndWin   = 5;                   % 3'-terminal window, both ends. Primer3 defines Max 3' Stability
                                % on the last five 3' bases, which is the scale dG3min is quoted on.
UseEndStabilityGate = true;     % Gate on dG3min rather than only scoring the 3'-terminal window
dG3min   = -5.08;               % kcal/mol on the PUBLISHED scale, used only when UseEndStabilityGate is on
                                %
                                % WHICH SCALE THIS IS. Primer3's Max 3' Stability INCLUDES helix
                                % initiation. EndStability returns the published-scale number
                                % first and the stacking-only number second, and the two differ by
                                % the initiation: 1.94 kcal/mol for a GC-terminated pentamer and
                                % 2.06 for an AT one. Measured against Primer3's own worked values,
                                % GCGCG is -8.793 stacking-only and -6.856 with the terminal rows
                                % restored, against their published 6.86; TATAT is -2.933 and
                                % -0.877 against their 0.86. Compare dG3min against the first of
                                % the two, never the second.
                                %
                                % WHY IT DOES NOT GATE BY DEFAULT. On a mixed-composition target
                                % the rule rejects almost nothing while a cluster of primers sits
                                % within a few thousandths of a kcal/mol of the line. A gate that
                                % decides nothing while standing one rounding error from
                                % eliminating most of the reverse family is not a gate. Rychlik's
                                % own framing is relative priming efficiency, not pass or fail, so
                                % the quantity is carried in Score. Set UseEndStabilityGate = true to
                                % restore the hard cut.
dGEndMin = -3.0;                % kcal/mol, least stable EXTENDABLE dimer allowed. Which registers
                                % count as extendable depends on the enzyme: 3'-anchored ones
                                % without proofreading, every register with it, because the
                                % exonuclease trims a dangling 3' tail. See the ProofExo block
                                % below. This constant carries most of the dimer protection, since
                                % phiDimer covers sequestration alone, so it is the number to
                                % revisit first if dimers appear.
dGhpmin  = -3.0;                % kcal/mol, free-energy hairpin fallback, used ONLY when the
                                % annealing temperature is unavailable
DimerTemplateBeyond = 1;        % nt, least number of partner nucleotides that must lie 5' of the
                                % base paired to a primer's 3' terminus before that dimer counts as
                                % extendable (stages 12 and P4, 3'-anchored minimum dGffend, dGrrend
                                % and dGfrend). A polymerase extends the 3' terminus by copying the
                                % partner, and the partner supplies a template only past the paired
                                % base. A duplex whose partner ends at that base leaves nothing to
                                % copy. 1 is the physical minimum; 0 restores the rule that any
                                % register pairing the terminal base counts. It does not act under a
                                % proofreading class, where the gates read the whole-molecule minimum.
HpTmMargin = 5.0;               % Degrees of headroom on the hairpin gate. A hairpin that pairs the
                                % primer 3'-terminal base rejects the candidate once its melting
                                % temperature reaches AnnealingTemperature - HpTmMargin, that is
                                % once half or more of the primer is folded with its 3'-OH buried
                                % while annealing.
                                %
                                % WHY THE HEADROOM IS NOT ZERO. The gate is unimolecular, so the
                                % melting temperature is a complete statement of occupancy: at the
                                % hairpin Tm exactly half the primer population carries a buried
                                % 3'-OH, with no concentration dependence to soften it. A zero
                                % margin therefore accepts a primer that is 49 per cent blocked
                                % while annealing. The headroom covers block-to-well temperature
                                % spread, which runs a few tenths of a degree on a typical
                                % instrument, and the residual error on a two-point melting
                                % temperature. A structure melting within a few degrees of the
                                % annealing step is treated as present rather than absent.

%% ---- AMPLICON 3'-OH COLLAPSE (LOG ONLY) ----------------------------------------------------------
% The accumulating single strand ends in a free 3'-OH. If that end folds back onto the strand's own
% body, leaves its terminal base paired and has template beyond the stem, the polymerase extends it
% and the strand lengthens in an unintended way. ScaffoldCollapseCheck screens both strands of every
% accepted pair for this. The result is WRITTEN to Pair_Log and Candidate_Primers and is never used
% to reject a pair or to move Score. Gating waits until the logged numbers have been read against
% real runs.
%
% The hairpin is unimolecular, so its melting temperature has no concentration dependence and is a
% complete statement of occupancy. A fold persists when its Tm reaches AnnealingTemperature less
% HpTmMargin, the same test stage 13 applies to primers. Tiers of one fold:
%   REJECT     persists and the last three 3' bases are paired (eps = 1)
%   FLAG       persists and a 3' mismatch limits extension, with eps >= CollapseFlagEps
%   TRANSIENT  does not persist but is extendable, eps >= CollapseFlagEps
%   NONE       otherwise
% Loops longer than 30 nt take the worm-like-chain increment of Aksel et al. 2024 and Dunn et al.
% 2015, anchored at the 30 nt table value. The engine's Jacobson-Stockmayer form is logged beside it.
%
% The same 3'-OH can pair with the same stem on a second copy of the strand. That bimolecular state
% owes helix initiation and has no loop, and its weight is [S]*exp(-dGbi/RT) against exp(-dGfold/RT)
% for foldback. [S] is the ceiling of the free strand: the excess primer concentration for the plus
% strand and the limiting primer concentration for the minus strand. At [S] = 1 uM it overtakes
% foldback for loops longer than about 1800 to 2500 nt, which only the longest loops of a
% multi-kilobase strand reach.
%
% Two enzyme channels are scored for every strand. Taq extends an end with at most one unpaired
% 3' base and weights it by the Rejali efficiency. A proofreading exonuclease trims an unpaired flap
% back to the outermost paired base, so efficiency is 1 and a flap of up to CollapseExoEndTol bases
% counts. A blend in which a fraction phi of the 3' end engagements is served by the proofreading
% enzyme gives (1 - phi)*pExtTaq + phi*pExtExo. No published value of phi exists for NEB LongAmp Taq,
% so CollapseExoBlend is NaN and pExtBlend stays empty until the user supplies one. The tier named
% TierAll is the highest rank over both channels and both topologies.
CollapseLog = true;             % Screen accepted pairs for 3'-OH collapse and log the result. Costs
                                % one scan per distinct terminus position, kept for the whole run.
CollapseWindow = 18;            % nt, length of the 3'-terminal window searched for stems
CollapseMinRun = 3;             % nt, shortest contiguous Watson-Crick run that counts as a stem
CollapseFlagEps = 1E-3;         % Lowest extension efficiency that still earns a FLAG or TRANSIENT
                                % tier. A fold limited below this by a 3' mismatch is NONE. A manual
                                % threshold, and the first number to revisit once the log is read.
CollapseLBase = 0.6;            % nm per nt, contour length of ssDNA in the loop polymer model
CollapseLPers = 0.9;            % nm, persistence length of ssDNA, Rivetti et al. 1998 equation 5
CollapseBimolecular = true;     % Include the bimolecular channel, in which the 3' tail pairs with
                                % the same stem on a second copy of the strand. False restores the
                                % foldback-only totals.
CollapseExoBlend = NaN;         % Fraction phi of 3' end engagements served by a proofreading
                                % enzyme, between 0 and 1. NaN leaves the blend column empty. The
                                % Taq and Exo channels are always reported separately.
CollapseExoEndTol = 15;         % nt, longest unpaired 3' flap the Exo channel counts. Clamped to
                                % CollapseWindow - CollapseMinRun. The Taq channel always uses one.

%% ---- SPECIFICITY -------------------------------------------------------------------------------
phi        = 1E-4;              % MISPRIMING tolerance: the ratio of unwanted to wanted complex AT
                                % ONE SITE, which is spurious product against intended product. A
                                % misprime competes for TEMPLATE, so the primer-template duplex is
                                % the right denominator and dGMispTarget and dGMispOffTarget keep
                                % it. Per-site: the scan scores of order 122 anchored registers
                                % per primer on a 7 kb circle, so the aggregate is nearer 122*phi
                                % than phi.
phiDimer   = 1E-3;              % DIMER tolerance: the fraction of the PRIMER POOL held as dimer,
                                % which is a different denominator and a deliberate one.
                                % A dimer competes for PRIMER, not for template. Referenced to the
                                % primer-template duplex the threshold moves 4.5 kcal/mol across
                                % three decades of template loading, which is backwards: removing a
                                % species that is not a reactant cannot make a dimer more harmful.
                                % Referenced to the primer pool it does not move with template at
                                % all. See the derivations in dGHomo.m and dGHetero.m.
                                % WHY IT IS LOOSER THAN phi. This gate covers SEQUESTRATION only,
                                % and sequestration is the weaker of the two dimer concerns: 1 per
                                % cent of the limiting primer lost is 1 per cent less duplex before
                                % the switch to linear synthesis, which is below the run-to-run
                                % scatter of everything else. The concern that actually ruins a
                                % reaction is EXTENSIBILITY, because an extended dimer amplifies as
                                % a competing template, and that is carried by dGEndMin above.
                                % At 1 uM and 50 nM primer this gives homodimer thresholds of
                                % -5.62 and -7.59 and a cross-dimer threshold of -6.07 kcal/mol at
                                % 58 C, computed at a representative primer-template dG of -22.
                                % Each sits within 0.02 kcal/mol of its dilute-limit closed form
                                % and moves by less than 0.001 kcal/mol across three decades of
                                % template, which is the point of the choice.
                                % THE VALUE IS CONDITIONED ON THE ENZYME. With a proofreading
                                % exonuclease declared it is tightened to phi below, because
                                % dGEndMin can no longer carry extensibility on its own.
EndKmer    = 7;                 % Length of the 3'-terminal k-mer tested for uniqueness
EndKmerMax = 1;                 % Maximum occurrences of that k-mer on both strands. One means the
                                % primer's own site and nowhere else.

% ---- CONTIGUOUS NON-COGNATE RUNS ----
% A primer binds its own site perfectly. Everywhere else on the circle it binds through some number
% of contiguous Watson-Crick pairs interrupted by mismatches, and the longest uninterrupted run is
% what decides whether that register nucleates at all. The mispriming scan scores those registers
% thermodynamically, which is the right currency and the expensive one. These two gates ask the
% combinatorial question first and answer it in O(1) per candidate from maps built once.
%
% TWO THRESHOLDS, BECAUSE THERE ARE TWO PHYSICAL SITUATIONS AND THEY ARE NOT THE SAME.
%
% A run that CARRIES THE PRIMER 3' TERMINUS is a nucleus the polymerase can work from. Taq has no
% 3' -> 5' proofreading exonuclease, so a paired 3'-OH is the precondition for adding a nucleotide,
% and Kwok et al. showed that a mismatch at the terminus is what stops PCR rather than merely
% slowing it. A terminal run is therefore both a binding site and a catalytically competent one,
% and it is held to the stricter standard.
%
% A run that lies ENTIRELY INSIDE the primer, with at least the 3'-terminal base outside it, cannot
% be extended without proofreading however long it is. What it still does is nucleate a duplex that
% competes for primer and that can zip outwards into a register the scan will then have to score.
% That is a weaker concern, so the internal threshold is looser.
%
% WHERE SEVEN COMES FROM. The expected number of occurrences of one specific k-mer elsewhere on a
% circular duplex of N bases, counting both strands, is 2N/4^k. On a 7249 nt circle that is 3.54 at
% k = 6, 0.885 at k = 7 and 0.221 at k = 8, so a unique 7-mer is already an unusual word at that
% template length and a unique 6-mer essentially does not exist. Seven is also the floor of
% detectability for the standard off-target search: Primer-BLAST runs at a default BLAST word size
% of seven, so an unintended target sharing six or fewer consecutive bases with a primer is not
% found at all. Below seven this gate rejects everything and carries no information.
%
% WHERE NINE COMES FROM. The internal gate is two bases looser, which is a factor of sixteen in the
% chance expectation, because the register it describes cannot prime on its own. Nine leaves a
% workable pool on a kilobase template while still removing the primers that carry a long internal
% homology, which are the ones whose off-target registers zip furthest once nucleated.
%
% THE THRESHOLDS AND THE TEMPLATE LENGTH ARE NOT INDEPENDENT. The expectation scales as 2N/4^k, so
% the same constant is free on a 500 nt construct, where a 7-mer is expected elsewhere 0.061 times,
% and punishing on a 7 kb plasmid. AssertRunGatesSatisfiable measures the admitted pool on whatever
% template is actually supplied, before the fold and before the scan, and refuses a run whose gates
% admit nobody.
%
% RELATION TO THE 3'-TERMINAL K-MER TEST ABOVE. At MaxRun3 = EndKmer and EndKmerMax = 1 the two are
% the same test written twice: the terminal EndKmer-mer occurring more than once on the duplex is
% exactly a 3'-terminal non-cognate run of EndKmer or longer. The engine runs the k-mer test first,
% so at these settings the terminal run gate shows no rejections of its own and the k-mer stage
% holds them. The two part company as soon as either constant is moved, and both quantities are
% logged for every candidate, so a log can be re-read at a different threshold without rescanning.
Use3PrimeRunGate  = true;       % Gate the 3'-terminal contiguous non-cognate run at MaxRun3
MaxRun3    = 7;                 % Reject when the longest contiguous run ENDING AT the primer 3'
                                % terminus that also occurs somewhere else on either strand reaches
                                % this length. Rejection is pRun3 >= MaxRun3.
UseInternalRunGate = true;      % Gate the longest non-3'-terminal contiguous non-cognate run
MaxRunInt  = 9;                 % Reject when the longest contiguous run lying wholly within the
                                % primer's first L-1 bases, and occurring somewhere else on either
                                % strand, reaches this length. Rejection is pRunInt >= MaxRunInt.
                                % Raising it by one nucleotide multiplies the admitted pool by
                                % roughly four, since that is what each extra base does to the
                                % chance of a chance match.
RunMinReport = 4;               % Shortest run the maps resolve. Below this the quantity is reported
                                % as RunMinReport - 1 rather than computed, because every primer
                                % shares 3-mers with a kilobase of template and the number carries
                                % no information. It also bounds the cost of building the maps.
RunPoolFloor = 1;               % Candidates per side the two run gates must still admit, measured
                                % on the supplied template with no other screen in front of them,
                                % before the scan is allowed to start. A JUDGEMENT, stated as one.
                                % One means refuse only when a side is empty, which is the point at
                                % which no pair can exist by arithmetic. Raise it to refuse early on
                                % a pool too thin to survive the thermodynamic cascade; of order a
                                % tenth of combinatorially clean candidates do survive it, and
                                % pairing needs both members.

MispMinRun = 5;                 % Contiguous identities a misprime REGISTER must carry to be scored.
                                % This is a sensitivity setting on the scan and not a gate. It is
                                % unrelated to MaxRun3 and MaxRunInt above, which are gates on the
                                % primer's own words. A 3'-anchored five-base run on a stable
                                % upstream duplex primes perfectly well, and pairing a tight end
                                % test with a looser seed is the right shape.
MispSurveyRun = 5;              % SURVEY floor: the contiguous run a sub-register must carry to be
                                % PRODUCTION VALUE. Equal to MispMinRun, which makes fluxAnch and
                                % fluxAnchRun identical by construction and the ensemble diagnostic
                                % uninformative. Set it to 3 to run the calibration again; the gate
                                % and the margin are invariant to this floor by construction, since
                                % BestSubRegister reports the best REJECTION-ELIGIBLE sub-register
                                % separately, so selection does not change either way.
                                % WHAT THE 04-11 CALIBRATION FOUND, at a floor of 3: the ensemble
                                % over every anchored register sits a median 0.79 kcal/mol below
                                % the single worst site on the forward strand and 0.88 on the
                                % reverse, reaching 2.19. At RT = 0.652 that median is a factor of
                                % 3.4 in total off-target priming rate. The two sums differed on 56
                                % of 57 forward and 202 of 207 reverse primers scored.
                                % WHAT IT COST: 1876 s on the reverse scan alone, against 2289 s
                                % for the whole run.
                                % COSTED, against MispMinRun, which is the run it must carry to
                                % REJECT. Handed to BestSubRegister as its minrun.
                                % WHY THE TWO DIFFER. MisprimerMismatch keeps two ensemble sums,
                                % fluxAnch over every anchored register the survey finds and
                                % fluxAnchRun over the subset the rejection acts on. With this
                                % value equal to MispMinRun every register reaching the
                                % accumulators already satisfies the run condition, so the two
                                % sums are identical by construction and a calibration comparing
                                % them returns zero for a reason unrelated to the physics.
                                % WHY 3. The survey has to reach below MispMinRun = 5 to see the
                                % registers the gate discards. BestSubRegister needs at least two
                                % paired positions for a nearest-neighbour step, and the run maps
                                % resolve nothing shorter than RunMinReport = 4. Three admits the
                                % short runs the gate ignores, keeps the enumeration finite, and
                                % leaves a clear nucleotide of separation from the gate so the
                                % comparison has something to measure.
                                % COST. Each nucleotide below MispMinRun multiplies the admitted
                                % sub-register count by roughly four.
                                % LIMIT. SEEDLEN = min(5, MispMinRun) still decides which WINDOWS
                                % are examined, so this widens the sweep inside the windows the
                                % seed already found rather than finding new ones. If a
                                % calibration still shows the two sums close, SEEDLEN is the next
                                % lever. Must not exceed MispMinRun; SeedScan raises if it does.
MispEndTol = 1;                 % Bases of the primer 3' end allowed outside the misprime duplex.
                                % THIS IS AN ENZYME PROPERTY, not a tuning knob. The value 1 is a
                                % Taq number: Taq cannot remove a 3' overhang, so a register that
                                % leaves the terminus dangling is catalytically dead however tightly
                                % it binds. The ProofExo block below widens it.
MispMargin = 0.0;               % EXTRA kcal/mol of headroom over the threshold. This is a SECOND
                                % phi written in kcal/mol: m kcal/mol multiplies phi by
                                % exp(-m/RT), 1.02E-2 at 3.0 and 56 C. Keep it at zero and set
                                % the specificity target with phi alone.
MispEndWin = 5;                 % Window for the per-locus 3'-terminal stacking dG
Rychlik = struct('On', true, 'Max', -11.0, 'T', 37, 'Na', 1);
                                % Rychlik measured priming becoming DETECTABLE once the
                                % 3'-terminal CONTINUOUS duplex of a 5'-mismatched primer beat
                                % -11 kcal/mol with the target. That is a measured threshold on the
                                % run of Watson-Crick pairs carrying the 3' terminus, which is
                                % neither the whole anchored register nor the fixed pentamer of the
                                % 3'-stability stage, so MisprimerMismatch computes that object and
                                % gates on it. Reach at this reference state: 7 bp of alternating
                                % GC gets there, a mixed-composition run needs about 10, and an
                                % AT-rich run never does inside 14. Rarely tripped, which is what a
                                % detection threshold should look like.
                                % The reference state is an ASSUMPTION until the paper's methods are
                                % read; see the note in MisprimerMismatch.m. Set On = false to drop
                                % the gate while keeping the quantity logged.
                                % THE CAP IS NEITHER LIFTED NOR RELAXED UNDER PROOFREADING, AND
                                % THAT IS DELIBERATE RATHER THAN AN OVERSIGHT.
                                % Huang et al. report that Taq "was observed to bind with similar
                                % affinities to each of the correctly paired and mispaired
                                % primer-template 3'-ends", and conclude that its failure to
                                % extend mismatches "appears to be an intrinsic property of the
                                % enzyme and not due to an inability to bind". The discrimination
                                % that table encodes is therefore CATALYTIC, not thermodynamic. A
                                % 3' -> 5' exonuclease removes the catalytic block and leaves the
                                % occupancy exactly where it was, so the same 3' anchor primes
                                % MORE readily under proofreading, not less. Rychlik measured -11
                                % with that block in place; without it the true threshold moves
                                % toward less negative, which makes -11 already the permissive end
                                % of the range. Lifting it would be wrong in direction as well as
                                % in magnitude, so the cap stays and the extension weighting goes
                                % instead. See the ENZYME CONDITIONING block below.

%% ---- ENZYME CONDITIONING ------------------------------------------------------------------------
% THE SETTING IS THREE-VALUED, BECAUSE EXONUCLEASE STRENGTH SPANS AN ORDER OF MAGNITUDE.
%   'none'    no 3' -> 5' activity. Taq, AccuStart, Stoffel Fragment, Deep Vent (exo-).
%   'weak'    a Taq and Deep Vent blend such as NEB LongAmp Taq.
%   'strong'  Pfu, Vent exo+, Phusion, Q5, KOD.
% weak AND strong CONDITION THE GATES IDENTICALLY, and that is the honest position rather than a
% missing feature. Nothing published maps an exonuclease rate onto either the Huang extension
% efficiencies or Rychlik's threshold, so a cap graded by enzyme would be an invented number
% wearing a citation. What strong adds is a warning about single-stranded YIELD, which this suite
% does not model and which is the property that actually separates the two in this application.
ProofExoClass = 'none';
ProofExo = false;               % 3' -> 5' proofreading exonuclease declared. Derived from
                                % ProofExoClass below. It changes the SHAPE of three gates, not
                                % just a threshold.
UseExtensionEfficiency = true;  % Weight each off-target site by the rate at which the polymerase
                                % extends from its 3' terminus, following Huang et al. false
                                % recovers the occupancy-only screen exactly.
% The weighting describes an enzyme that CANNOT remove a mispaired 3' terminus. A 3' -> 5'
% proofreading exonuclease excises it and carries on, so for such an enzyme every epsilon collapses
% towards 1 and the screen has to fall back to occupancy. input.ProofreadingClass carries the
% answer from CollectaPCRData_Claude.m; an input struct built without it leaves the weighting on,
% which is correct for every enzyme the preset list offers.
if isfield(input, 'ProofreadingClass')&&~isempty(input.ProofreadingClass)&& ...
        (ischar(input.ProofreadingClass)||isstring(input.ProofreadingClass))
    ProofExoClass = lower(strtrim(char(input.ProofreadingClass)));
elseif isfield(input, 'ProofreadingExo')&&~isempty(input.ProofreadingExo)&& ...
        ~isnan(input.ProofreadingExo)&&input.ProofreadingExo
    % An input struct built before the class existed carries only the flag. A bare true says that
    % proofreading was declared and nothing about its strength, and weak conditions the gates in
    % exactly the way that flag used to, so an old struct keeps its old behaviour.
    ProofExoClass = 'weak';
end
if ~ismember(ProofExoClass, {'none', 'weak', 'strong'})
    error('PrimerSearchParameters:BadProofExoClass', ...
        ['input.ProofreadingClass must be ''none'', ''weak'' or ''strong''. Received ''%s''. ' ...
        'The older ''Y'' and ''N'' spellings are translated in CollectaPCRData_Claude.m and do ' ...
        'not reach here.'], ProofExoClass);
end
ProofExo = ~strcmp(ProofExoClass, 'none');
if ProofExo
    UseExtensionEfficiency = false;
    warning('PrimerSearchParameters:ProofreadingPolymerase', ...
        ['A 3'' -> 5'' proofreading polymerase was declared, so the mispriming margin is ' ...
        'thresholding OCCUPANCY and not occupancy times extension efficiency. Every off-target ' ...
        'site with a mispaired 3'' terminus is now treated as fully extendable, which is the ' ...
        'conservative reading for an enzyme that excises the mispair. Expect far fewer ' ...
        'candidates to survive. Declared class: %s.'], ProofExoClass);
    % The same excision changes the DIMER gates, and it changes their SHAPE and not only their
    % level. Without proofreading, a dimer is a threat only if the 3' terminus is PAIRED, since a
    % terminus hanging over as a single-stranded tail cannot be extended; dGEndMin on the anchored
    % minimum therefore carries extensibility by itself, and phiDimer only has to cover how much
    % primer is tied up. A proofreading exonuclease removes that tail. Any register stable enough
    % to persist at the annealing temperature is trimmed back to a flush 3' end and then extended,
    % so the distinction between anchored and internal registers stops meaning anything, and the
    % gate has to look at the WHOLE-MOLECULE minimum. Two moves follow, and neither introduces a
    % new judgement number:
    %   (i)  dGEndMin is applied to dGff, dGrr and dGfr rather than to dGffend, dGrrend and
    %        dGfrend. The whole-molecule minimum sits of order 1.35 kcal/mol below the anchored one
    %        on an accepted pair, so this is a real tightening and not a relabelling.
    %   (ii) phiDimer is tightened to phi, because the anchored gate no longer does the
    %        extensibility work on its own and the occupancy gate has to pick up the slack.
    % NOT MODELLED: the exonuclease also degrades unbound primer, 3' end first, which lowers the
    % effective primer concentration over the run and shifts every threshold that depends on it.
    % Pfu chews a free oligo at a rate comparable with its polymerisation rate under some buffer
    % conditions. Nothing here accounts for that; phosphorothioate 3' linkages are the usual
    % protection and the suite has no way to know whether they were ordered.
    % THE MISPRIMING STAGES CARRY THE SAME ASSUMPTION, and it has to go too.
    % MisprimerMismatch calls a register anchored when the duplex reaches within endTol of the
    % primer 3' terminus, and both the mispriming rejection and the margin read ONLY anchored
    % registers. That test is the same Taq premise the dimer gates carried: an overhang cannot be
    % extended. A proofreading exonuclease removes unpaired 3' overhanging nucleotides to create
    % blunt ends, so no overhang length disqualifies a register and the tolerance is the whole
    % primer. Widening endTol is the right instrument rather than a second code path, because the
    % right edge inside MisprimerMismatch is already PINNED at the outermost Watson-Crick pair.
    % That is exactly where the exonuclease stops; it does not chew into the helix. So the object
    % the scan builds under a wide tolerance IS the trimmed duplex, and not the free-edge optimum,
    % which would be a different and wrong duplex.
    % MEASURED COST on a 3940 nt target: the free-edge minimum, which is the loosest bound the
    % widening can reach, sits 1.1 to 2.6 kcal/mol below the anchored one on average and 3.8 at
    % worst, and every candidate keeps a positive margin against it, the tightest being
    % +1.34 kcal/mol. On that target the widening rejects nothing; on a different target it will.
    % THE RUN GATES ARE UNAFFECTED. A contiguous perfect run is a property of the two sequences and
    % of nothing else, so no enzyme moves it.
    MispEndTol = max(fLmax, rLmax);
    phiDimer = phi;
    warning('PrimerSearchParameters:ProofreadingDimerGate', ...
        ['Proofreading declared, so the dimer gates changed SHAPE: dGEndMin = %.2f kcal/mol is ' ...
        'now applied to the WHOLE-MOLECULE dimer minimum (dGff, dGrr, dGfr) instead of the ' ...
        '3''-anchored minimum, because the exonuclease trims an unpaired 3'' tail and makes ' ...
        'an internal register extendable. phiDimer tightened to %.0e. MispEndTol widened from 1 ' ...
        'to %d, the longest primer, so the mispriming stages stop requiring the 3'' terminus to ' ...
        'sit inside the misprime duplex. Primer degradation by the same exonuclease is NOT ' ...
        'modelled.'], dGEndMin, phiDimer, MispEndTol);
    % THE strong TIER ADDS A YIELD WARNING AND NOTHING ELSE.
    % The gates above are already conditioned identically for weak and strong. What separates the
    % two in THIS application is single-stranded yield, which this suite does not model at any
    % point, so the only honest thing to do with the distinction is to say what the literature
    % reports and let the reader decide.
    if strcmp(ProofExoClass, 'strong')
        warning('PrimerSearchParameters:StrongProofreadingYield', ...
            ['A STRONGLY proofreading polymerase was declared for a %d nt single-stranded ' ...
            'target. Veneziano et al. report that Phusion "was unable to synthesize fragment ' ...
            'larger than 1 kb ssDNA", which they attribute to "enhanced exonuclease activity", ' ...
            'and that "Taq polymerases gave the highest ssDNA production", "possibly due to ' ...
            'lower exonuclease activity". Both are the authors'' own hypotheses rather than ' ...
            'measurements, and they are worded that way in the paper. They used a Taq-based ' ...
            'enzyme at 1000 and 3281 nt and reached for LongAmp only at 10 kb and above. This ' ...
            'warning is about YIELD, which nothing in this suite models; the gates are ' ...
            'conditioned exactly as for a weak exonuclease. Veneziano R et al. Sci Rep ' ...
            '2018;8:6548, doi:10.1038/s41598-018-24677-5.'], round(input.ScaffoldLength));
    end
end

EffTableFile = 'extension_efficiency_rejali2018.csv';
                                % Extension-efficiency table. A bare filename is resolved beside
                                % this file. EMPTY NO LONGER MEANS HUANG: it means nothing is
                                % named here, and ExtensionEfficiency then resolves this same
                                % shipped CSV from its own folder. Put the word 'builtin' here to
                                % ask for the built-in Huang et al. 1992 values on purpose, which
                                % the warning in ExtensionEfficiency.m then explains. It passes
                                % through the resolver below untouched, since it names no file.
                                %
                                % WHY REJALI IS THE DEFAULT FOR EVERY PRESET, NOT JUST THE STOFFEL
                                % ONE. This block used to argue that Rejali was the right table for
                                % preset 3 alone, because KlenTaq 1 shares the Stoffel Fragment
                                % architecture, and that Huang stayed right for full-length Taq.
                                % That reasoning is retired. Two things overturn it.
                                %
                                % First, Huang leaves G.A, A.G, G.G and C.C at a detection limit of
                                % 1E-6, and those four are USED AS VALUES. A floor charged as
                                % -RT*ln(eps) becomes 9.09 kcal/mol at 58 C, which reads as "this
                                % site cannot prime" when the experiment only said the gel could
                                % not see it. G.G is the case that bites: Rejali measure 8.8E-4 by
                                % stopped-flow, Chen et al. 2015 class it outside the mismatches
                                % showing observable inhibition, and Ayyadevara et al. 2000 get
                                % 1.8 per cent of matched product from it in 35-cycle PCR on Taq
                                % itself. Three studies, two enzymes, three assays, against one
                                % censored bound.
                                %
                                % Second, the enzyme difference does not reorder the discrimination
                                % pattern. Over the eight mismatches Huang actually resolved, the
                                % two tables rank-correlate at +0.91 and differ by a near-constant
                                % factor of 19.3. A truncation that removed discrimination would
                                % scramble the order, not shift the scale. Huang's own result says
                                % why: Taq binds matched and mispaired 3' ends with similar
                                % affinity, so the selectivity is catalytic and lives in the
                                % polymerase domain the deletion leaves intact.
                                %
                                % WHAT THIS COSTS. Rejali's values are uniformly higher, so every
                                % off-target site scores as more extendable and the mispriming
                                % screen tightens. That is the conservative direction for a design
                                % tool, where a false flag costs one candidate and a missed flag
                                % costs an experiment. Base3Spread moves from 1.347 to 0.967
                                % decades and Wt3Base self-calibrates from 1.000 to 1.394, so the
                                % 3'-base ranking term keeps its weight in Score without
                                % intervention.
                                %
                                % WHAT IT DOES NOT FIX. Rejali measured at 65 C in their own
                                % buffer, and report that the optimum for mismatch extension sits
                                % at 55 to 60 C while KCl suppresses it. At a 58 C anneal these
                                % values may understate eps. Their three templates also spread a
                                % single mismatch by up to 7.7-fold, so any one number here is good
                                % to about half a decade. That spread, not the choice of table, is
                                % the resolution floor.
                                %
                                % Rejali NA, Moric E, Wittwer CT. Clin Chem 2018;64:801-9.
                                %   doi:10.1373/clinchem.2017.282285
                                % Chen F, Zhao Y, Fan C, Zhao Y. Anal Chem 2015;87:8718-23.
                                %   doi:10.1021/acs.analchem.5b01545
                                % Ayyadevara S, Thaden JJ, Shmookler Reis RJ. Anal Biochem
                                %   2000;284:11-18. doi:10.1006/abio.2000.4635
% THE CHOICE IS HELD IN ExtensionEfficiencyTable, NOT PASSED AS AN ARGUMENT.
% The three call sites are LocalPrimerSearchEngine.m:1877, MisprimerMismatch.m:597 and
% PrimerHybridizationMap.m:1374, all inside helpers with no access to this configuration. That
% file's header says why a persistent is the proportionate answer, and where the choice is
% announced so that it cannot be invisible.
% RESOLVED TO A FULL PATH BEFORE IT IS RECORDED. A bare filename on the Parameters sheet does not
% identify a file: it depends on the working directory MATLAB happened to be in. The sheet has to
% name the table a reader could go and open.
if ~isempty(EffTableFile)
    if ~isfile(EffTableFile)
        cand = fullfile(fileparts(mfilename('fullpath')), EffTableFile);
        if isfile(cand), EffTableFile = cand; end
    end
    dEff = dir(EffTableFile);
    if ~isempty(dEff)&&isfield(dEff, 'folder')&&~isempty(dEff(1).folder)
        EffTableFile = fullfile(dEff(1).folder, dEff(1).name);
    end
end
ExtensionEfficiencyTable('set', EffTableFile);
clear cand dEff

%% ---- TEMPLATE ACCESSIBILITY --------------------------------------------------------------------
UseAccessibility = true;        % Needs DNAAccessibility.m, RNAplfold and dna_mathews2004.par
AccW       = 240;               % RNAplfold sliding window, nt
AccL       = 168;               % RNAplfold maximum pairing span, nt
AccSaltOverride = [];           % Molar monovalent for the fold. EMPTY DERIVES IT FROM THE
                                % REACTION, which is the only defensible default. A hardcoded fold
                                % salt produces no symptom while it happens to match the buffer,
                                % which is the worst way for a constant to be wrong: change the
                                % salt in the dialog and every Tm and every free energy in the
                                % suite moves while the FOLD stays where it was, so the
                                % accessibility map describes a different reaction from the one
                                % being designed, silently. The fold is also cached on its
                                % parameters, so a stale salt would be preserved run after run.
                                % Set it to a number only to answer a what-if, and expect the
                                % console line and the Parameters sheet to say you did.
UseAccSaltEquiv = true;         % RNAplfold --salt carries NO magnesium term.
                                % false passes the monovalent concentration alone, which is what
                                % the option was parameterised on. It UNDERSTATES the real ionic
                                % screening, so structure comes out less stable than it is and the
                                % template reads as more accessible than it is. Optimistic, in the
                                % direction that lets a bad primer through rather than rejecting a
                                % good one.
                                % true substitutes the von Ahsen equivalent, [Na+] + 120*sqrt([Mg]
                                % - [dNTP]) with concentrations in mM, which is the right ionic
                                % strength but feeds a monovalent correction a number no monovalent
                                % buffer would produce. On an ordinary PCR buffer the two settings
                                % differ by more than a factor of four, so they are not a small
                                % perturbation of each other. Run both and keep only conclusions
                                % that survive both.
if ~isempty(AccSaltOverride)
    AccSalt = AccSaltOverride;
    AccSaltSrc = 'explicit override';
elseif UseAccSaltEquiv
    % Through the shared speciation rather than a local [Mg] - [dNTP] stand-in. The difference is
    % under one percent at ordinary PCR compositions, because the chelation constant Ka = 3E4 /M is
    % large enough that dNTP binds magnesium almost stoichiometrically, and it grows once dNTP
    % approaches magnesium, where the stand-in clamps hard at zero while the quadratic does not.
    % Routing both through NNParams also guarantees that the salt the fold is computed at is the
    % same number every other free energy in the suite is computed at.
    AccSalt = AccSaltEquivalent(input.MonovalentCationConcentration, ...
        input.BivalentCationConcentration, input.dNTPsConcentration);
    AccSaltSrc = 'von Ahsen equivalent of the reaction buffer, on the NNParams speciation';
else
    AccSalt = input.MonovalentCationConcentration;
    AccSaltSrc = 'input.MonovalentCationConcentration';
end
AccMaxOpen = 8.0;               % kcal/mol; reject a site costing more than this to open
AccParamFile = 'C:\Program Files (x86)\ViennaRNA Package\Misc\dna_mathews2004.par';
                                % Full path to dna_mathews2004.par, set explicitly rather than left
                                % to DNAAccessibility's own search. The search itself is sound and
                                % hard-stops when nothing is found, so the fold never silently runs
                                % on RNA parameters. What the empty setting costs is the CACHE KEY:
                                % it records 'auto' and cannot see the energy file change under it,
                                % so swapping the .par would leave a stale fold looking valid. With
                                % the path here the key carries the file's path, byte count and
                                % modification time, which closes that.
                                % VERIFIED 05-11-2026 on NANOBIOLOGY: 381156 bytes, modified
                                % 21-12-2025, header "## RNAfold parameter file v2.0" with the DNA
                                % note on line 3, and a stacking alphabet of CG GC GT TG AT TA. The
                                % T rather than U is what separates it from rna_turner2004.par,
                                % which sits in the same folder at a near-identical 381610 bytes.
if ~isempty(AccParamFile)&&~isfile(AccParamFile)
    % Loud rather than silent, and then back to the previous behaviour. A hardcoded absolute path
    % is correct for one machine and wrong for the next, and the failure it would otherwise cause
    % is a hard stop inside DNAAccessibility several steps later.
    warning('PrimerSearchParameters:AccParamFileMissing', ...
        ['AccParamFile is set to %s, which does not exist on this machine. DNAAccessibility will ' ...
        'fall back to its own search and the accessibility cache key will record ''auto'' instead ' ...
        'of the file identity, so a change of energy parameters will not invalidate the cache. ' ...
        'Correct the path for this machine, or clear the constant deliberately to silence this.'], ...
        AccParamFile);
    AccParamFile = '';
end

%% ---- ACCESSIBILITY CACHE -----------------------------------------------------------------------
% The fold is the slow step. It is also a PURE FUNCTION of the sequence and seven fold parameters,
% so a run that repeats those parameters is recomputing a result it already had. The cache keeps
% lunp for both strands on disk under the plasmid's own folder, keyed on everything that can change
% the answer.
% WHAT THE KEY COVERS. The full plus-strand sequence, its length, and W, L, uMax, Temp, Salt,
% Circular and NucLen. uMax is max(fLmax, rLmax) + 2 and Temp is Tmin - dTsub, so changing a primer
% length window or the Tm floor invalidates the cache, which is correct and not obvious.
% WHAT IT CANNOT COVER unless AccParamFile is set. When that constant is empty DNAAccessibility
% resolves the energy file itself and the engine never learns which one it used, so swapping
% dna_mathews2004.par would leave a stale cache looking valid. The key records 'auto' in that case
% and the loader says so once. Set AccParamFile and the key carries the path, byte count and
% modification time, which closes the hole.
% THE FILENAME CARRIES THE PARAMETERS so the folder is readable at a glance, and ends in eight hex
% digits of a hash over the whole key. The filename is a convenience, never the test: the loader
% compares the full key stored INSIDE the file, sequence and all, and treats any mismatch as a miss.
AccCacheOn   = true;            % false recomputes every run and writes nothing
AccCacheRoot = '';              % Empty picks F:\aPCR_Computations\DNA_Accessibility when F: exists
                                % and D:\ otherwise, matching how aPCR_Claude_v2.m chooses resFold.
                                % Set it to any absolute path to override.
AccCacheWrite = true;           % false reads an existing cache but never creates one, which is what
                                % you want on a shared drive you are not supposed to fill

%% ---- RANKING WEIGHTS ---------------------------------------------------------------------------
WtMargin = 1.0; MarginCap = 5.0;   % Per kcal/mol of dimer headroom, and its ceiling
WtMisp = 0.2;   MispCap   = 10.0;  % Per kcal/mol of mispriming headroom, and its ceiling
% OFF-TARGET MULTIPLICITY, RANKED AND NOT GATED.
% WHAT IS SCORED IS THE ENSEMBLE SHIFT, dG_ens - dG_eff_min, and NOT the ensemble margin.
% The shift is -RT*ln(M_eff): it is what the Boltzmann sum over every anchored register adds on
% top of the single worst one, and it is zero for a primer whose off-target burden is one site.
% Scoring the ensemble MARGIN instead, as dG_ens - dG_min, would carry the per-site margin inside
% it, and WtMisp already weights that margin. The two terms would then double-count the per-site
% part and differ only by this shift, which is the part that is actually new.
% WHY IT RANKS RATHER THAN GATES. A gate on the ensemble at phi = 1E-4 eliminates 8 of 9 forward
% candidates on this template and leaves no valid pair at all. Beyond that, the ensemble sum runs
% over whatever registers the survey found, so it moves with SEEDLEN and MispSurveyRun, which are
% sensitivity settings on the scan rather than properties of the primer. A threshold on a quantity
% that shifts when the search is widened means something different at every setting. The per-site
% minimum does not have that problem: BestSubRegister reports the best rejection-eligible
% sub-register separately, so selection is invariant to the survey floor.
EnsCap = 2.5;                   % kcal/mol, ceiling on the ensemble shift before it is weighted.
                                % The 04-11 calibration at a survey floor of 3 found a median
                                % shift of 0.79 kcal/mol on the forward strand and 0.88 on the
                                % reverse, reaching 2.19 at worst. 2.5 clears the measured range
                                % with a little room and bounds a term that would otherwise grow
                                % with the log of a register count.
WtEns = MarginCap/EnsCap;       % TIED TO MarginCap BY CONSTRUCTION, exactly as WtAccess is.
                                % The term's full span is then MarginCap Score points, so
                                % off-target multiplicity can move a candidate by as much as dimer
                                % headroom can and no further. Changing EnsCap or MarginCap means
                                % this number follows automatically.
WtAccess = 0.3125;              % Per kcal/mol of summed opening cost, forward site plus reverse.
                                % TIED TO AccMaxOpen BY CONSTRUCTION: MarginCap/(2*AccMaxOpen).
                                % The gate admits a site costing up to AccMaxOpen to open, so a
                                % pair can carry twice that before it is rejected. At a weight of
                                % 1.0 the term's full span was 16 Score points, larger than every
                                % other term's ceiling put together, and larger than the whole
                                % spread of Score across the pairs that survive. At this weight the
                                % span is exactly MarginCap, so accessibility can move a candidate
                                % by as much as dimer headroom can and no further. Changing
                                % AccMaxOpen or MarginCap means recomputing this number.
WtTmA    = 2.0;                 % Score points per E-FOLD loss of early linear-phase slope,
                                % read off Pierce's four measured points. The old form was per
                                % degree of distance from a target, which penalised a gap narrower
                                % than the target as heavily as a wider one and said nothing about
                                % how much yield the distance costs. In these units the number is a
                                % decision: at 2.0 a pair predicted to yield 1/e of the best tested
                                % rate gives up two Score points, which puts amplicon competition on
                                % a par with the hairpin term. Across the gaps a long amplicon can
                                % actually reach, roughly 19 to 26 degrees, the whole span is worth
                                % about 2.6 points at this weight.
TmACap   = 2.5;                 % Ceiling on the amplicon yield penalty, in E-FOLDS of lost slope,
                                % applied before WtTmA. Pierce's widest tested gap is 23 degrees,
                                % where the penalty is 1.59. Beyond that the log-slope trend is
                                % extended rather than measured, and an extension that runs to
                                % arbitrary width would let this one term outweigh every other.
                                % At 2.5 the term is bounded by WtTmA*TmACap = 5.0 points, equal to
                                % MarginCap. The cap first binds near a 28 degree gap, which is
                                % wider than any pair the search has produced, so ordering survives
                                % across the whole realistic range and only runaway gaps flatten.
WtFlank  = 1.0;                 % Score points across the WHOLE admissible amplicon overhang.
                                % The term is WtFlank*flank/(fFlex + rFlex), so it is bounded by
                                % WtFlank by construction and the number reads directly against the
                                % other ceilings: at 1.0 it is the smallest term in Score, below
                                % even the Tm-offset window at 1.5.
                                % WHAT IT EXPRESSES. flank = ampLen - ScaffoldLength, the bases the
                                % amplicon carries beyond the scaffold. It is non-negative, since
                                % both search windows only ever grow outward. A shorter overhang is
                                % preferred, weakly: this states a direction without letting
                                % geometry overrule thermodynamics.
                                % WHAT IT WILL NOT DO AT THIS WEIGHT. On the 04-11 run the best
                                % pair overhangs by 57 nt and leads the best pair within 20 nt by
                                % 4.571 points, a gap no flank weight below about 20 can close, and
                                % 20 would make this term twice the dimer ceiling. At 1.0 the best
                                % compliant pair moves from rank 222 to 205 of 315. Raise it if the
                                % overhang is a hard design constraint rather than a preference; in
                                % that case a gate on flank is the honest instrument, not a weight.
WtHairpin= 2.0;                 % Per unit folded fraction at the annealing step, 0 to 1
Wt3End   = 0.5;                 % Per kcal/mol of 3'-end stability, published scale
WtProf   = 0.5;                 % Per kcal/mol of internal-stability profile.
                                % Rychlik's recommendation is that the 3' end be LOOSER than the
                                % 5' end. In the stacking-only currency the initiation cancels
                                % between the two windows, so dG3stack - dG5stack is clean and
                                % positive is good. On a GC-rich target most surviving reverse
                                % primers fail it, every one ending in a GC-rich stretch, which is
                                % the case Rychlik singles out as prone to false priming. As a hard
                                % gate it would empty the reverse set, so it is scored.
% Wt3Base IS CALIBRATED AGAINST WHICHEVER TABLE IS ACTIVE, rather than being a fixed number.
% The term is -log10 of the geometric mean of epsilon over the three mispairs a 3'-terminal base
% can form as the PRIMER base. A LARGER value is better discrimination, so it enters Score with a
% MINUS sign. What the weight has to do is make the SPREAD across the four bases worth a fixed
% number of score points, because that spread is the whole discriminating power of the term.
%       Huang    G 5.18  A 4.96  C 3.84  T 3.84     spread 1.347 decades
%       Rejali   G 3.06  A 3.41  C 2.57  T 2.44     spread 0.967 decades
% At Wt3Base = 1.0 the Huang spread is worth 1.35 points, which the original comment called about
% as much as 1.35 kcal/mol of dimer headroom. Holding that fixed gives 1.393 on the Rejali table.
% Computing it here rather than storing a number means the term carries the same weight in Score
% whichever table is chosen, and a future table needs no hand edit.
BASE3REF = 1.347;               % decades, the Huang spread, the reference the weight preserves
base3dec = zeros(1, 4);
wc3 = struct('A', 'T', 'C', 'G', 'G', 'C', 'T', 'A');
allb3 = 'ACGT';
for k3 = 1:4
    b3 = allb3(k3);
    o3 = setdiff(allb3, wc3.(b3));
    acc3 = 0;
    for q3 = 1:numel(o3), acc3 = acc3 + log10(ExtensionEfficiency(b3, o3(q3))); end
    base3dec(k3) = -acc3/numel(o3);
end
Base3Spread = max(base3dec) - min(base3dec);
if ~isfinite(Base3Spread)||(Base3Spread <= 0)
    error('PrimerSearchParameters:Base3SpreadDegenerate', ...
        ['The 3''-terminal base discrimination spread came out %.4g. Four bases scoring alike ' ...
        'makes the Wt3Base term carry no information, so a table that does that is almost ' ...
        'certainly wrong.'], Base3Spread);
end
Wt3Base  = BASE3REF/Base3Spread;
clear BASE3REF base3dec wc3 allb3 k3 b3 o3 acc3 q3

%% ---- PROGRESS REPORTING ------------------------------------------------------------------------
Verbose       = true;           % Narrate the run on the console.
                                % An engine spends most of its wall time in two places that
                                % produce no output of their own: the RNAplfold template fold,
                                % which is a single external call over the whole circle, and the
                                % pair loop, which is O(forward x reverse). Without narration a run
                                % is indistinguishable from a hang, and the waitbar does not appear
                                % at all under -nodisplay or -batch. Every line goes to stdout, so
                                % it survives redirection into a log file.
ProgressEvery = 50;             % Console progress is emitted every this many candidates scanned.
                                % The waitbar updates on every iteration regardless.

%% ---- PARAMETER SHEET ---------------------------------------------------------------------------
% EVERY SETTING A RUN DEPENDS ON, in one place: what the user supplied through
% CollectaPCRData_Claude, and what the search was configured with. The two are kept visibly apart
% by the Section column and by the input. prefix on the names, because they are changed in
% different places and a reader chasing a result needs to know which file to open.
ParamSpec = {
% Section                Name                      Value                   Units        Meaning
'Template (input)'     , 'input.PlasmidFileName' , input.PlasmidFileName , 'file'     , 'Sequence file the run was launched on'
'Template (input)'     , 'input.PlasmidSequence' , SeqText(input.PlasmidSequence)   , 'summary'  , 'Length, GC fraction and an FNV-1a checksum. The checksum is what makes a run reproducible, since two files can share a name'
'Template (input)'     , 'input.PlasmidStranded' , input.PlasmidStranded , '1 or 2'   , 'Strandedness of the template as supplied'
'Template (input)'     , 'input.PlasmidConcentration'  , ConcText(input.PlasmidConcentration)  , 'M'  , 'Template molarity, converted from the ng/uL the user entered'
'Template (input)'     , 'input.ScaffoldStart'   , input.ScaffoldStart   , 'nt'       , 'First base of the scaffold on the plus strand'
'Template (input)'     , 'input.ScaffoldEnd'     , input.ScaffoldEnd     , 'nt'       , 'Last base of the scaffold on the plus strand'
'Template (input)'     , 'input.ScaffoldLength'  , input.ScaffoldLength  , 'nt'       , 'Scaffold length requested'
'Primers (input)'      , 'input.ForwardPrimer'   , input.ForwardPrimer   , 'sequence' , 'Chosen forward primer. Empty during a search, filled for a yield test'
'Primers (input)'      , 'input.ReversePrimer'   , input.ReversePrimer   , 'sequence' , 'Chosen reverse primer. Empty during a search, filled for a yield test'
'Primers (input)'      , 'input.ForwardPrimerConcentration' , ConcText(input.ForwardPrimerConcentration) , 'M' , 'EXCESS primer starting concentration. Sets its Tm and the dimer thresholds'
'Primers (input)'      , 'input.ReversePrimerConcentration' , ConcText(input.ReversePrimerConcentration) , 'M' , 'LIMITING primer starting concentration, entered as a ratio to the forward'
'Buffer (input)'       , 'input.MonovalentCationConcentration' , ConcText(input.MonovalentCationConcentration) , 'M' , 'Total monovalent cation. Enters every Tm and every nearest-neighbour free energy'
'Buffer (input)'       , 'input.BivalentCationConcentration'   , ConcText(input.BivalentCationConcentration)   , 'M' , 'Total magnesium. Owczarzy corrects every primer Tm for the FREE fraction after dNTP chelation'
'Buffer (input)'       , 'input.dNTPsConcentration'            , ConcText(input.dNTPsConcentration)            , 'M' , 'Total dNTP. Chelates magnesium, so raising it lowers every primer Tm'
'Thermal programme (input)', 'input.AnnealingTemperature' , AnnealText(input.AnnealingTemperature) , 'C or mode' , 'Annealing step. A number is a fixed temperature; 0, 1 and 2 select programmed, variable and touchdown; a vector is the touchdown recipe'
'Thermal programme (input)', 'input.ExtensionTemperature' , input.ExtensionTemperature , 'C'    , 'Extension step, set by the polymerase'
'Thermal programme (input)', 'input.NumAmplificationCycles' , input.NumAmplificationCycles , 'cycles' , 'Cycles requested'
'Polymerase (input)'   , 'input.PolymeraseName'  , input.PolymeraseName  , 'name'     , 'Polymerase selected'
'Polymerase (input)'   , 'input.PolymeraseActivity' , input.PolymeraseActivity , 'U/uL' , 'Activity, converted from the mU/uL the user entered'
'Polymerase (input)'   , 'input.ExtensionTime'   , input.ExtensionTime   , 'min'      , 'Extension time per cycle'
'Polymerase (input)'   , 'input.UnitDNTP'        , input.UnitDNTP        , 'mol'      , 'dNTP incorporated per unit, converted from nmol'
'Polymerase (input)'   , 'input.UnitTime'        , input.UnitTime        , 'min'      , 'Window the unit definition is measured over'
'Polymerase (input)'   , 'input.PolymeraseEfficiency' , input.PolymeraseEfficiency , 'fraction' , 'Relative activity at the extension temperature, 0 to 1'
'Polymerase (input)'   , 'input.MaxLaps'         , input.MaxLaps         , 'laps'     , 'Maximum laps on a circular template'
'Polymerase (input)'   , 'input.ProofreadingClass', ClassText(input)      , 'class'    , 'none, weak or strong. The value ProofExoClass is read from; absent on a struct built before the class existed'
'Polymerase (input)'   , 'input.ProofreadingExo' , FlagText(input.ProofreadingExo) , 'yes/no' , 'A 3''->5'' exonuclease. Rewrites ProofExo, UseExtensionEfficiency, MispEndTol and phiDimer below'
'Polymerase (input)'   , 'input.FivePrimeExo'    , FlagText(input.FivePrimeExo)    , 'yes/no' , 'A 5''->3'' exonuclease. Governs what happens to a primer bound downstream of an extending polymerase'
'Search window'        , 'fFlex'                 , fFlex                 , 'nt'       , 'Forward 5'' start positions scanned upstream of ScaffoldStart'
'Search window'        , 'rFlex'                 , rFlex                 , 'nt'       , 'Reverse 5'' end positions scanned downstream of ScaffoldEnd'
'Search window'        , 'fLmin'                 , fLmin                 , 'nt'       , 'Shortest forward primer'
'Search window'        , 'fLmax'                 , fLmax                 , 'nt'       , 'Longest forward primer. Also sets the fold uMax, so changing it invalidates the accessibility cache'
'Search window'        , 'rLmin'                 , rLmin                 , 'nt'       , 'Shortest reverse primer'
'Search window'        , 'rLmax'                 , rLmax                 , 'nt'       , 'Longest reverse primer. Also sets the fold uMax'
'Melting temperature'  , 'Tmin'                  , Tmin                  , 'C'        , 'Lowest primer Tm accepted, at that primer''s own concentration'
'Melting temperature'  , 'Tmax'                  , Tmax                  , 'C'        , 'Highest primer Tm accepted. Not a Pierce number; see the note beside it'
'Melting temperature'  , 'dTrfmin'               , dTrfmin               , 'C'        , 'Smallest limiting-minus-excess Tm offset. Pierce measured efficiency rising to +7'
'Melting temperature'  , 'dTrfmax'               , dTrfmax               , 'C'        , 'Largest limiting-minus-excess Tm offset'
'Melting temperature'  , 'dTsub'                 , dTsub                 , 'C'        , 'Degrees below Tm at which secondary structure is evaluated'
'Amplicon competition' , 'TmAmaxGap'             , TmAmaxGap             , 'C'        , 'ADVISORY ONLY, no longer a gate. Reported by the startup feasibility check'
'Amplicon competition' , 'TmAoptGap'             , TmAoptGap             , 'C'        , 'Gap at which the yield penalty is zero. Pierce best tested value'
'Amplicon competition' , 'TmAPierceGap'          , mat2str(TmAPierceGap)   , 'C'      , 'Pierce TmA - TmX values, the ranking interpolates between them'
'Amplicon competition' , 'TmAPierceSlope'        , mat2str(TmAPierceSlope) , 'slope'  , 'Early linear-phase slope Pierce measured at each of those gaps'
'Amplicon competition' , 'UseTmAEquiv'           , UseTmAEquiv           , 'flag'     , 'true reports TmA at the von Ahsen monovalent equivalent. Does not reach the gate'
'Composition'          , 'UseGCWindow'           , UseGCWindow           , 'flag'     , 'Enforce the global GC fraction window'
'Composition'          , 'fGCmin'                , fGCmin                , 'fraction' , 'Lowest GC fraction accepted'
'Composition'          , 'fGCmax'                , fGCmax                , 'fraction' , 'Highest GC fraction accepted'
'Composition'          , 'MonoLim'               , MonoLim               , 'nt'       , 'Longest mononucleotide run allowed'
'Composition'          , 'DiRepLim'              , DiRepLim              , 'repeats'  , 'Longest dinucleotide repeat allowed'
'Terminal bases'       , 'Force5primeGC'         , Force5primeGC         , 'flag'     , 'Require a G or C at the 5'' end of both primers'
'Terminal bases'       , 'ForwardForce3primeAT'  , ForwardForce3primeAT  , 'flag'     , 'Require an A or T at the forward 3'' end'
'Terminal bases'       , 'ForwardAvoid3primeT'   , ForwardAvoid3primeT   , 'flag'     , 'Reject a forward primer ending in T'
'Terminal bases'       , 'ReverseForce3primeGC'  , ReverseForce3primeGC  , 'flag'     , 'Require a G or C at the reverse 3'' end'
'Terminal bases'       , 'ForwardGCSkew'         , ForwardGCSkew         , 'flag'     , 'Gate the forward primer on local 3''-versus-5'' GC skew'
'Terminal bases'       , 'ReverseGCSkew'         , ReverseGCSkew         , 'flag'     , 'Gate the reverse primer on local 3''-versus-5'' GC skew'
'3'' end stability'     , 'EndWin'                , EndWin                , 'nt'       , 'Window for the 3''-terminal stacking free energy'
'3'' end stability'     , 'UseEndStabilityGate'            , UseEndStabilityGate            , 'flag'     , 'Gate the 3''-terminal window on dG3min rather than only scoring it'
'3'' end stability'     , 'dG3min'                , dG3min                , 'kcal/mol' , 'Least stable 3'' terminal window allowed when the gate is on'
'3'' end stability'     , 'Rychlik.On'            , Rychlik.On            , 'flag'     , 'Apply the Rychlik cap on per-locus 3'' terminal stability'
'3'' end stability'     , 'Rychlik.Max'           , Rychlik.Max           , 'kcal/mol' , 'Rychlik cap value'
'3'' end stability'     , 'Rychlik.T'             , Rychlik.T             , 'C'        , 'Temperature the Rychlik cap is evaluated at'
'3'' end stability'     , 'Rychlik.Na'            , Rychlik.Na            , 'M'        , 'Monovalent salt the Rychlik cap is evaluated at'
'Dimers and hairpins'  , 'dGEndMin'              , dGEndMin              , 'kcal/mol' , 'Least stable EXTENDABLE dimer allowed. Carries most of the dimer protection'
'Dimers and hairpins'  , 'dGhpmin'               , dGhpmin               , 'kcal/mol' , 'Hairpin free-energy fallback, used only when the annealing temperature is unavailable'
'Dimers and hairpins'  , 'DimerTemplateBeyond'   , DimerTemplateBeyond   , 'nt'       , 'Partner nucleotides required 5'' of the base paired to a primer 3'' terminus for the dimer to count as extendable. 0 restores the earlier rule'
'Dimers and hairpins'  , 'HpTmMargin'            , HpTmMargin            , 'C'        , 'Headroom on the hairpin melting-temperature gate'
'Amplicon collapse'    , 'CollapseLog'           , CollapseLog           , 'flag'     , 'Screen accepted pairs for amplicon 3''-OH collapse. Log only, never gates or scores'
'Amplicon collapse'    , 'CollapseWindow'        , CollapseWindow        , 'nt'       , 'Length of the 3''-terminal window searched for stems'
'Amplicon collapse'    , 'CollapseMinRun'        , CollapseMinRun        , 'nt'       , 'Shortest contiguous Watson-Crick run that counts as a stem'
'Amplicon collapse'    , 'CollapseFlagEps'       , CollapseFlagEps       , 'ratio'    , 'Lowest extension efficiency that still earns a FLAG or TRANSIENT tier'
'Amplicon collapse'    , 'CollapseLBase'         , CollapseLBase         , 'nm'       , 'Contour length per nt of ssDNA in the loop polymer model beyond 30 nt'
'Amplicon collapse'    , 'CollapseLPers'         , CollapseLPers         , 'nm'       , 'Persistence length of ssDNA in the loop polymer model beyond 30 nt'
'Amplicon collapse'    , 'CollapseBimolecular'   , CollapseBimolecular   , 'flag'     , 'Include the bimolecular channel: the 3'' tail pairing with the same stem on a second strand'
'Amplicon collapse'    , 'CollapseExoBlend'      , CollapseExoBlend      , 'fraction' , 'Fraction of 3'' end engagements served by a proofreading enzyme. NaN leaves the blend column empty'
'Amplicon collapse'    , 'CollapseExoEndTol'     , CollapseExoEndTol     , 'nt'       , 'Longest unpaired 3'' flap the Exo channel counts. The Taq channel counts one'
'Specificity'          , 'phi'                   , phi                   , 'ratio'    , 'Mispriming tolerance, off-target over cognate duplex AT ONE SITE'
'Specificity'          , 'phiDimer'              , phiDimer              , 'fraction' , 'Dimer tolerance, fraction of the PRIMER POOL held as dimer'
'Specificity'          , 'EndKmer'               , EndKmer               , 'nt'       , 'Length of the 3''-terminal k-mer tested for uniqueness'
'Specificity'          , 'EndKmerMax'            , EndKmerMax            , 'count'    , 'Most occurrences of that k-mer allowed on the circle'
'Specificity'          , 'Use3PrimeRunGate'      , Use3PrimeRunGate      , 'flag'     , 'Gate the 3''-terminal contiguous non-cognate run. Logged either way'
'Specificity'          , 'MaxRun3'               , MaxRun3               , 'nt'       , 'Reject when the 3''-terminal non-cognate contiguous run reaches this length'
'Specificity'          , 'UseInternalRunGate'    , UseInternalRunGate    , 'flag'     , 'Gate the longest non-3''-terminal contiguous non-cognate run. Logged either way'
'Specificity'          , 'MaxRunInt'             , MaxRunInt             , 'nt'       , 'Reject when a non-cognate run inside the primer, clear of the 3'' terminus, reaches this length'
'Specificity'          , 'RunMinReport'          , RunMinReport          , 'nt'       , 'Shortest contiguous run the maps resolve. Below it the value is reported as RunMinReport - 1'
'Specificity'          , 'RunPoolFloor'          , RunPoolFloor          , 'count'    , 'Candidates per side the run gates alone must admit before the scan is allowed to start'
'Specificity'          , 'MispMinRun'            , MispMinRun            , 'nt'       , 'Contiguous identities a misprime REGISTER must carry to be scored. Unrelated to MaxRun3 and MaxRunInt'
'Specificity'          , 'MispSurveyRun'         , MispSurveyRun         , 'nt'       , 'Run a sub-register must carry to be COSTED. Below MispMinRun so the two ensemble sums differ'
'Specificity'          , 'MispEndTol'            , MispEndTol            , 'nt'       , 'Primer 3'' bases allowed outside the misprime duplex. Widened automatically under proofreading'
'Specificity'          , 'MispMargin'            , MispMargin            , 'kcal/mol' , 'Extra headroom over the mispriming threshold. Keep at 0 and set the target with phi'
'Specificity'          , 'MispEndWin'            , MispEndWin            , 'nt'       , 'Window for the per-locus 3''-terminal stacking free energy'
'Enzyme'               , 'ProofExoClass'         , ProofExoClass         , 'class'    , 'none, weak or strong. weak and strong condition the gates identically; strong adds a yield warning'
'Enzyme'               , 'ProofExo'              , ProofExo              , 'flag'     , 'True for weak and strong. Derived from ProofExoClass, which is set from input.ProofreadingClass'
'Enzyme'               , 'UseExtensionEfficiency', UseExtensionEfficiency, 'flag'     , 'Weight off-target sites by extension efficiency. Forced off under proofreading'
'Accessibility'        , 'UseAccessibility'      , UseAccessibility      , 'flag'     , 'Fold the template and charge each site its opening cost'
'Accessibility'        , 'AccW'                  , AccW                  , 'nt'       , 'RNAplfold sliding window'
'Accessibility'        , 'AccL'                  , AccL                  , 'nt'       , 'RNAplfold maximum pairing span'
'Accessibility'        , 'AccSaltOverride'       , ConcText2(AccSaltOverride) , 'M'   , 'Forces the fold salt. Empty derives it from the reaction, which is the default'
'Accessibility'        , 'UseAccSaltEquiv'       , UseAccSaltEquiv       , 'flag'     , 'true feeds the von Ahsen equivalent to RNAplfold instead of the monovalent concentration alone'
'Accessibility'        , 'AccSalt'               , ConcText(AccSalt)     , 'M'        , 'The salt the fold ACTUALLY ran at, after the two settings above'
'Accessibility'        , 'AccSaltSrc'            , AccSaltSrc            , 'provenance' , 'Where AccSalt came from. Anything other than input.MonovalentCationConcentration means the fold and the rest of the run describe different buffers'
'Accessibility'        , 'AccMaxOpen'            , AccMaxOpen            , 'kcal/mol' , 'Reject a site costing more than this to open'
'Accessibility'        , 'AccParamFile'          , AccParamFile          , 'path'     , 'dna_mathews2004.par. Empty lets DNAAccessibility find it and leaves the cache unable to see a swap'
'Accessibility cache'  , 'AccCacheOn'            , AccCacheOn            , 'flag'     , 'Reuse a stored fold when the sequence and fold parameters match'
'Accessibility cache'  , 'AccCacheRoot'          , AccCacheRoot          , 'path'     , 'Cache root. Empty picks the F: or D: default'
'Accessibility cache'  , 'AccCacheWrite'         , AccCacheWrite         , 'flag'     , 'Create new cache entries. false reads an existing cache but never writes'
'Ranking weights'      , 'WtMargin'              , WtMargin              , 'per kcal' , 'Weight on the worst dimer headroom'
'Ranking weights'      , 'MarginCap'             , MarginCap             , 'kcal/mol' , 'Ceiling on that headroom'
'Ranking weights'      , 'WtMisp'                , WtMisp                , 'per kcal' , 'Weight on the worst mispriming headroom'
'Ranking weights'      , 'MispCap'               , MispCap               , 'kcal/mol' , 'Ceiling on that headroom'
'Ranking weights'      , 'WtAccess'              , WtAccess              , 'per kcal' , 'Weight on the summed opening cost, set to MarginCap/(2*AccMaxOpen)'
'Ranking weights'      , 'TmACap'                , TmACap                , 'e-folds'  , 'Ceiling on the amplicon yield penalty before WtTmA is applied'
'Ranking weights'      , 'WtFlank'               , WtFlank               , 'points'   , 'Score points across the whole admissible amplicon overhang, fFlex + rFlex'
'Ranking weights'      , 'WtTmA'                 , WtTmA                 , 'per e-fold', 'Score points per e-fold loss of linear-phase slope, from Pierce'
'Ranking weights'      , 'WtHairpin'             , WtHairpin             , 'per unit' , 'Weight on the worst folded fraction at the annealing step'
'Ranking weights'      , 'Wt3End'                , Wt3End                , 'per kcal' , 'Weight on the 3''-terminal stacking free energy'
'Ranking weights'      , 'WtProf'                , WtProf                , 'per kcal' , 'Weight on the internal-stability profile'
'Enzyme'               , 'EffTableFile'          , EffTableFile          , 'path'     , 'Extension-efficiency table. Empty uses the built-in Huang et al. 1992 values'
'Ranking weights'      , 'EnsCap'                , EnsCap                , 'kcal/mol' , 'Ceiling on the off-target ensemble shift before WtEns is applied'
'Ranking weights'      , 'WtEns'                 , WtEns                 , 'per kcal' , 'Weight on the off-target ensemble shift, set to MarginCap/EnsCap so its full span equals dimer headroom'
'Ranking weights'      , 'Base3Spread'           , Base3Spread           , 'decades'  , 'Spread of 3''-base discrimination across A C G T on the ACTIVE table. Wt3Base holds its Score contribution at 1.347'
'Ranking weights'      , 'Wt3Base'               , Wt3Base               , 'per unit' , 'Weight on 3''-terminal base discrimination, calibrated against Base3Spread'
'Reporting'            , 'Verbose'               , Verbose               , 'flag'     , 'Narrate the run on the console'
'Reporting'            , 'ProgressEvery'         , ProgressEvery         , 'count'    , 'Console progress emitted every this many candidates'
};
Parameters = table(string(ParamSpec(:, 1)), string(ParamSpec(:, 2)), ...
    string(cellfun(@ParamText, ParamSpec(:, 3), 'UniformOutput', false)), ...
    string(ParamSpec(:, 4)), string(ParamSpec(:, 5)), ...
    'VariableNames', {'Section', 'Parameter', 'Value', 'Units', 'Meaning'});
AssertParamSheetComplete(who, fieldnames(input), ParamSpec);

%% ---- THE RETURNED CONFIGURATION ------------------------------------------------------------------
% One field per setting, named exactly as the setting is named above and in the engines. The list
% below is a manifest and it is checked against the sheet immediately afterwards, in both
% directions, so it cannot drift from the settings it is supposed to carry.
cfg = struct( ...
    'fFlex', fFlex, 'rFlex', rFlex, 'fLmin', fLmin, 'fLmax', fLmax, 'rLmin', rLmin, 'rLmax', rLmax, ...
    'Tmin', Tmin, 'Tmax', Tmax, 'dTrfmin', dTrfmin, 'dTrfmax', dTrfmax, 'dTsub', dTsub, ...
    'TmAmaxGap', TmAmaxGap, 'TmAoptGap', TmAoptGap, 'UseTmAEquiv', UseTmAEquiv, ...
    'TmAPierceGap', TmAPierceGap, 'TmAPierceSlope', TmAPierceSlope, ...
    'UseGCWindow', UseGCWindow, 'fGCmin', fGCmin, 'fGCmax', fGCmax, ...
    'MonoLim', MonoLim, 'DiRepLim', DiRepLim, ...
    'Force5primeGC', Force5primeGC, 'ForwardForce3primeAT', ForwardForce3primeAT, ...
    'ForwardAvoid3primeT', ForwardAvoid3primeT, 'ReverseForce3primeGC', ReverseForce3primeGC, ...
    'ForwardGCSkew', ForwardGCSkew, 'ReverseGCSkew', ReverseGCSkew, ...
    'EndWin', EndWin, 'UseEndStabilityGate', UseEndStabilityGate, 'dG3min', dG3min, ...
    'dGEndMin', dGEndMin, 'dGhpmin', dGhpmin, 'HpTmMargin', HpTmMargin, ...
    'DimerTemplateBeyond', DimerTemplateBeyond, ...
    'CollapseLog', CollapseLog, 'CollapseWindow', CollapseWindow, 'CollapseMinRun', CollapseMinRun, ...
    'CollapseFlagEps', CollapseFlagEps, 'CollapseLBase', CollapseLBase, 'CollapseLPers', CollapseLPers, ...
    'CollapseBimolecular', CollapseBimolecular, 'CollapseExoBlend', CollapseExoBlend, ...
    'CollapseExoEndTol', CollapseExoEndTol, ...
    'phi', phi, 'phiDimer', phiDimer, 'EndKmer', EndKmer, 'EndKmerMax', EndKmerMax, ...
    'Use3PrimeRunGate', Use3PrimeRunGate, 'MaxRun3', MaxRun3, ...
    'UseInternalRunGate', UseInternalRunGate, 'MaxRunInt', MaxRunInt, ...
    'RunMinReport', RunMinReport, 'RunPoolFloor', RunPoolFloor, ...
    'MispMinRun', MispMinRun, 'MispSurveyRun', MispSurveyRun, ...
    'MispEndTol', MispEndTol, 'MispMargin', MispMargin, ...
    'MispEndWin', MispEndWin, 'Rychlik', Rychlik, ...
    'ProofExo', ProofExo, 'ProofExoClass', ProofExoClass, ...
    'UseExtensionEfficiency', UseExtensionEfficiency, ...
    'UseAccessibility', UseAccessibility, 'AccW', AccW, 'AccL', AccL, ...
    'AccSaltOverride', AccSaltOverride, 'UseAccSaltEquiv', UseAccSaltEquiv, ...
    'AccSalt', AccSalt, 'AccSaltSrc', AccSaltSrc, 'AccMaxOpen', AccMaxOpen, ...
    'AccParamFile', AccParamFile, ...
    'AccCacheOn', AccCacheOn, 'AccCacheRoot', AccCacheRoot, 'AccCacheWrite', AccCacheWrite, ...
    'WtMargin', WtMargin, 'MarginCap', MarginCap, 'WtMisp', WtMisp, 'MispCap', MispCap, ...
    'WtEns', WtEns, 'EnsCap', EnsCap, ...
    'WtAccess', WtAccess, 'WtTmA', WtTmA, 'TmACap', TmACap, 'WtFlank', WtFlank, ...
    'WtHairpin', WtHairpin, ...
    'Wt3End', Wt3End, 'WtProf', WtProf, 'Wt3Base', Wt3Base, ...
    'EffTableFile', EffTableFile, 'Base3Spread', Base3Spread, ...
    'Verbose', Verbose, 'ProgressEvery', ProgressEvery);
AssertConfigMatchesSheet(cfg, ParamSpec);
end

%% ____________________________________________________________________________________________________________
function t = ConcText2(v)
%% CONCTEXT2  ConcText that renders an EMPTY override as what empty means
if isempty(v), t = '(empty: derived from the reaction)'; return; end
t = ConcText(v);
end

%% ____________________________________________________________________________________________________________
function t = ConcText(v)
%% CONCTEXT  A molar concentration in both the stored unit and the one it was typed in
% Everything in the suite is molar, and nobody types molar. Showing both is what lets a reader
% check the conversion without recomputing it.
if isempty(v)||~isnumeric(v)||~isscalar(v)||~isfinite(v), t = '-'; return; end
if v == 0
    t = '0';
elseif v >= 1E-3
    t = sprintf('%.4g M  (%.4g mM)', v, v*1E3);
elseif v >= 1E-6
    t = sprintf('%.4g M  (%.4g uM)', v, v*1E6);
else
    t = sprintf('%.4g M  (%.4g nM)', v, v*1E9);
end
end

%% ____________________________________________________________________________________________________________
function t = AnnealText(v)
%% ANNEALTEXT  The annealing field, which carries a temperature OR a mode OR a recipe
if isempty(v), t = '-'; return; end
if ~isscalar(v)
    t = sprintf('touchdown recipe, %d steps: %s', numel(v), ...
        strjoin(arrayfun(@(x) sprintf('%.1f', x), v(:)', 'UniformOutput', false), ', '));
    return
end
switch v
    case 0, t = '0  (programmed: each cycle set by the engine)';
    case 1, t = '1  (variable)';
    case 2, t = '2  (touchdown)';
    otherwise, t = sprintf('%.1f C  (fixed)', v);
end
end

%% ____________________________________________________________________________________________________________
function t = FlagText(v)
%% FLAGTEXT  A yes/no field that CollectaPCRData may leave as NaN
if isempty(v), t = 'not stated'; return; end
if isnumeric(v)&&isscalar(v)&&isnan(v), t = 'not stated'; return; end
if logical(v), t = 'yes'; else, t = 'no'; end
end

%% ____________________________________________________________________________________________________________
function t = SeqText(seq)
%% SEQTEXT  The template in one cell: length, composition and a checksum
s = upper(char(seq));
n = numel(s);
if n == 0, t = '(empty)'; return; end
gc = 100*nnz((s == 'G')|(s == 'C'))/n;
t = sprintf('%d nt, GC %.2f%%, FNV-1a %s', n, gc, Fnv1a32(s));
end

%% ____________________________________________________________________________________________________________
function t = ParamText(v)
%% PARAMTEXT  One cell of the parameter sheet, rendered for a human
if ischar(v), t = v; return; end
if isstring(v), t = char(v); return; end
if isempty(v), t = '-'; return; end
if islogical(v)
    if v, t = 'true'; else, t = 'false'; end
    return
end
if isnumeric(v)&&isscalar(v)
    if v == fix(v)&&(abs(v) < 1E6)
        t = sprintf('%d', v);
    elseif (abs(v) < 1E-3)||(abs(v) >= 1E5)
        t = sprintf('%.3g', v);
    else
        t = sprintf('%.4g', v);
    end
    return
end
if isnumeric(v)
    t = ['[' strjoin(arrayfun(@(x) sprintf('%.4g', x), v(:)', 'UniformOutput', false), ' ') ']'];
    return
end
t = '(unrenderable)';
end

%% ____________________________________________________________________________________________________________
function t = ClassText(input)
%% CLASSTEXT  The declared exonuclease class, or what stood in for it on an older input struct
if isfield(input, 'ProofreadingClass')&&~isempty(input.ProofreadingClass)
    t = lower(strtrim(char(input.ProofreadingClass)));
elseif isfield(input, 'ProofreadingExo')&&~isempty(input.ProofreadingExo)&& ...
        ~isnan(input.ProofreadingExo)&&input.ProofreadingExo
    t = 'weak (inferred from the ProofreadingExo flag)';
else
    t = 'none (not declared)';
end
end

%% ____________________________________________________________________________________________________________
function eq = AccSaltEquivalent(monv, biv, dNTP)
%% ACCSALTEQUIVALENT  von Ahsen monovalent equivalent on the shared magnesium speciation
% A LOCAL FUNCTION RATHER THAN A WORKSPACE VARIABLE, DELIBERATELY. AssertParamSheetComplete is
% handed who() and warns about anything in this function's workspace that is not on the Parameters
% sheet, because such a thing is a setting the run depends on and the reader will never see. An
% NNParams struct held in a variable here would trip that guard while being no kind of setting.
% Holding it in a persistent inside a local function keeps the workspace clean and still builds the
% parameter module only once per session.
persistent N
if isempty(N), N = NNParams(); end
eq = N.Salt.MonovalentEquivalent(monv, biv, dNTP);
end

%% ____________________________________________________________________________________________________________
function hex = Fnv1a32(str)
%% FNV1A32  FNV-1a 32-bit over the bytes of a char row, returned as 8 lowercase hex digits
% THE MULTIPLY IS SPLIT INTO 16-BIT HALVES, AND IT HAS TO BE. The state reaches 2^32 and the prime
% is 16777619, so the product reaches 7.2E16. A double holds integers exactly only to 2^53, which is
% 9.0E15, so a plain h*p rounds to the nearest representable value for every state above 2^29. That
% is seven eighths of the state space, and because the state feeds the next iteration the result
% diverges from FNV-1a at the first byte that trips it rather than drifting slowly.
%
% Splitting keeps every partial product exact:
%
%       h*p mod 2^32  =  ( (hiH*p mod 2^16)*2^16 + loH*p )  mod 2^32
%
% with hiH and loH below 2^16 and p below 2^24, so each product stays under 2^40.
%
% The arithmetic is done in doubles rather than uint32 or uint64 on purpose: MATLAB integer types
% SATURATE at their maximum instead of wrapping, and a hash needs wrapping.
%
% This is byte-for-byte the routine AccKeyHash uses in LocalPrimerSearchEngine.m, so the checksum
% printed on the Parameters sheet and the cache key computed during a run are the same function of
% the same bytes. Verified against the canonical FNV-1a 32 test vectors: '' gives 811c9dc5,
% 'a' gives e40c292c and 'foobar' gives bf9cf968.
b   = double(uint8(char(str)));
hv  = 2166136261;                 % offset basis
p   = 16777619;                   % prime
M32 = 4294967296;
for i = 1:numel(b)
    hv  = bitxor(hv, b(i));
    hiH = floor(hv/65536);
    loH = hv - hiH*65536;
    hv  = mod(mod(hiH*p, 65536)*65536 + loH*p, M32);
end
hex = lower(dec2hex(hv, 8));
end

%% ____________________________________________________________________________________________________________
function AssertParamSheetComplete(wsNames, inFields, spec)
%% ASSERTPARAMSHEETCOMPLETE  Warn when a setting exists but is not on the sheet
% Called at the moment the sheet is built, when the workspace holds the settings and nothing else.
% Anything in the workspace that is not on the sheet is a setting a run depends on and a reader
% will never see, which is the failure this guards against.
listed = spec(:, 2);
listedEng = listed(~startsWith(listed, 'input.'));
% A struct-valued setting appears on the sheet one row per field, as Name.Field, so the struct
% itself counts as listed when any of its fields is.
listedRoot = unique(cellfun(@(s) strtok(s, '.'), listedEng, 'UniformOutput', false));
skip = {'input', 'ParamSpec', 'Parameters', 'cfg', 'ans'};
missing = setdiff(wsNames(:)', [listedRoot(:)', skip]);
if ~isempty(missing)
    warning('PrimerSearchParameters:SheetIncomplete', ...
        ['These settings exist but are not on the Parameters sheet: %s.\n' ...
        'Add them to ParamSpec, or the sheet will tell the reader a run used defaults it did ' ...
        'not use.'], strjoin(missing, ', '));
end
% The same test for the user's side, so a field added to the input struct is noticed too.
listedIn = strrep(listed(startsWith(listed, 'input.')), 'input.', '');
missingIn = setdiff(inFields(:)', listedIn(:)');
if ~isempty(missingIn)
    warning('PrimerSearchParameters:InputSheetIncomplete', ...
        'These input fields are not on the Parameters sheet: %s.', strjoin(missingIn, ', '));
end
end

%% ____________________________________________________________________________________________________________
function AssertConfigMatchesSheet(cfg, spec)
%% ASSERTCONFIGMATCHESSHEET  The returned struct and the sheet must describe the same settings
% Both directions, because the two failures are different and both are silent. A setting on the
% sheet but absent from cfg is one an engine cannot reach, and it will fail at the unpack line with
% an unhelpful message. A field in cfg but absent from the sheet is one a reader of the run record
% will never see.
listed = spec(:, 2);
listedEng = listed(~startsWith(listed, 'input.'));
listedRoot = unique(cellfun(@(s) strtok(s, '.'), listedEng, 'UniformOutput', false));
have = fieldnames(cfg);
onSheetOnly = setdiff(listedRoot(:)', have(:)');
inCfgOnly   = setdiff(have(:)', listedRoot(:)');
if ~isempty(onSheetOnly)||~isempty(inCfgOnly)
    error('PrimerSearchParameters:ConfigSheetMismatch', ...
        ['The returned configuration and the Parameters sheet disagree.\n' ...
        '  On the sheet but not in cfg: %s\n' ...
        '  In cfg but not on the sheet: %s\n' ...
        'Both lists must be empty. Add the setting to whichever of the two is missing it.'], ...
        JoinOrNone(onSheetOnly), JoinOrNone(inCfgOnly));
end
end

%% ____________________________________________________________________________________________________________
function t = JoinOrNone(c)
%% JOINORNONE  Comma-separated list, or the word none, so an error message reads as a sentence
if isempty(c), t = '(none)'; else, t = strjoin(c, ', '); end
end
