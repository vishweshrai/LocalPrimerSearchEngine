function [tab, results] = LocalPrimerSearchEngine(input)
%% LOCALPRIMERSEARCHENGINE  Fixed-scaffold primer pair search for asymmetric (LATE) PCR
% Written by Vishwesh Rai, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 12-11-2026 01:53 IST
%
% Revision history is NOT kept here. See LocalPrimerSearchEngine_<YYYYMMDD>.txt beside this file.
% Every constant, flag and threshold is set in PrimerSearchParameters.m, which is also where the
% reason for each value is written down. This file carries the algorithm and nothing else.
% ____________________________________________________________________________________________________________
% WHAT THIS FUNCTION DOES, END TO END
%
% It scans a narrow window around a requested scaffold, builds every forward and reverse primer that
% fits inside it, pushes each one through a conjunctive filter cascade, then pairs the survivors
% under the LATE-PCR melting-temperature relations and ranks the accepted pairs.
%
% SEARCH GEOMETRY
%
% Forward candidates start at every position from ScaffoldStart - fFlex to ScaffoldStart, at every
% length from fLmin to fLmax. Reverse candidates end at every position from ScaffoldEnd to
% ScaffoldEnd + rFlex, at every length from rLmin to rLmax. Both windows only ever grow the
% amplicon, never shrink it, so the requested scaffold is always fully contained. The plasmid is
% treated as circular throughout: circshift rotates it so the candidate sits at index 1, register
% scans run over a doubled string so sites crossing the origin are found, and every reported
% coordinate is folded back with mod.
%
% A forward candidate of length L occupies plus(fPos : fPos+L-1), with its 3' terminus at the HIGH
% end of that interval. A reverse candidate of length L occupies plus(rPos-L+1 : rPos) and is the
% reverse complement of it, so its 3' terminus sits at the LOW end. Every position-indexed lookup
% below depends on that asymmetry and states which end it is reading.
%
% ____________________________________________________________________________________________________________
% THE PHYSICS, STAGE BY STAGE
%
% Stages are ordered cheapest first. Every filter is conjunctive, so the accepted set does not
% depend on the order, but the run time does: the whole-plasmid mispriming scan is the most
% expensive stage and is reached only by candidates that have cleared everything else. A candidate
% is charged to the FIRST gate it fails, which is what makes the reject tally readable.
%
%  1. TERMINAL BASES.  Composition rules on the 3'-terminal base and optionally on the 5' one.
%     The distinction the rules rest on is between two different 3'-end failures. Extension
%     efficiency, which enters stage 14 as epsilon, measures how fast the polymerase extends a
%     MISPAIRED terminus. A terminal-base clamp guards against a terminus that is correctly paired
%     and weakly stacked, which frays; extension from a transiently unpaired 3'-OH is a rate, and
%     every screen in this cascade integrates over bound states and is blind to it. The two are not
%     the same quantity and the second does not subsume the first. The base is logged and scored
%     under Wt3Base whether or not any clamp is enforced.
%
%  2. GLOBAL GC FRACTION.  A composition window. Historically this WAS the melting temperature:
%     Tm = 2(A+T) + 4(G+C) in the era of the Wallace rule, so at fixed length the window and a Tm
%     band are one constraint stated twice, and stage 10 computes Tm properly from the SantaLucia
%     parameters at the true strand concentration. What survives the redundancy appears only at
%     FIXED Tm, where a GC-poor primer must be longer and the longer duplex dilutes the free-energy
%     cost of a single mismatch. Stage 14 addresses the same quantity with better information,
%     scoring mismatched duplexes against the actual template. The upper bound answers the opposite
%     effect, that GC-rich primers at matched Tm fold more readily, which stages 12 and 13 measure
%     directly.
%
%  3. HALF-PRIMER GC SKEW.  A GC-poorer 3' half than 5' half lowers the free energy available to a
%     3'-anchored mispriming event while leaving the overall melting temperature alone, which is
%     the internal stability profile Rychlik recommends, in composition rather than free-energy
%     form. The free-energy form is computed at stage 8 as dG3stack - dG5stack and enters the Score
%     under WtProf, so enforcing a two-bin threshold on GC counts rejects candidates on a coarsened
%     copy of a quantity already being optimised. Logged always, gated only if asked.
%
%  4. LOW COMPLEXITY.  Longest homopolymer run and longest perfect dinucleotide repeat. A repeat
%     offers many registers that differ by one repeat unit, so it multiplies the number of nearly
%     equivalent binding sites without changing the free energy of any one of them.
%
%  5. 3'-TERMINAL K-MER UNIQUENESS.  The 3'-terminal EndKmer-mer must occur no more than EndKmerMax
%     times across BOTH strands of the circular plasmid. This is combinatorial rather than
%     thermodynamic, and it expresses mispriming directly: a primer whose 3' end has one home on the
%     template cannot nucleate an extendable duplex anywhere else. Two strfind calls per candidate.
%
%  6. 3'-TERMINAL CONTIGUOUS NON-COGNATE RUN.  The longest contiguous word ENDING AT the primer's
%     3' terminus that also occurs somewhere on the duplex other than at the primer's own site.
%     Rejected once it reaches MaxRun3. A run carrying the 3' terminus is the nucleus an off-target
%     product can actually be made from, because without a 3' -> 5' proofreading exonuclease a
%     paired 3'-OH is the precondition for adding a nucleotide. The quantity generalises stage 5
%     from one fixed word length to the longest length that still matches, so the log can be
%     re-read at a different threshold without rescanning.
%
%  7. INTERNAL CONTIGUOUS NON-COGNATE RUN.  The same question asked of the primer's first L-1
%     bases, so every run counted here is clear of the 3' terminus by at least one base. Rejected
%     once it reaches MaxRunInt. Such a register cannot be extended by an enzyme that carries no
%     proofreading exonuclease, however long it is, and that is why the threshold is looser than
%     the one above rather than equal to it. What it still does is nucleate a duplex that competes
%     for primer and that can zip outwards into a register stage 14 then has to score.
%
%     HOW BOTH ARE COMPUTED, since the cost would otherwise be prohibitive. Two maps are built once
%     for the whole circle, before the fold and before the scan, and each candidate then reads its
%     two numbers out of an array.
%       NonUniq{k}(p)   true when the k-mer beginning at plus-strand position p occurs MORE THAN
%                       ONCE across the two strands, counting its own occurrence. False means the
%                       primer sharing that word has exactly one home on the duplex.
%       RunMap          per window length L and per plus-strand coordinate of the window's LOW end,
%                       the longest such k inside the window, at its high end, and at its low end.
%     COUNTING BOTH STRANDS INTO ONE NUMBER IS WHAT MAKES ONE TABLE SERVE BOTH PRIMERS. A primer
%     pairs its template through Watson-Crick complementarity, so a contiguous match of the primer
%     against the minus strand IS an occurrence of that word on the plus strand, and a match against
%     the plus strand is an occurrence of its reverse complement. The total over the two strands is
%     therefore invariant under reverse complementation. A forward primer reads its words off the
%     plus strand and a reverse primer reads the reverse complements of plus-strand words, so both
%     consult the same table at the same plus-strand coordinates.
%     THE COUNT IS EXACT ON A CIRCLE. Every start position on each strand is enumerated from a
%     doubled string, so a word straddling the origin is counted once and only once.
%     MONOTONICITY IS WHAT MAKES A LONGEST-RUN SEARCH CHEAP. If a (k+1)-mer occurs twice then both
%     of its k-mer halves occur twice, so the maps can be filled in increasing k and let the larger
%     value overwrite the smaller.
%     WHAT THE TWO GATES DO NOT KNOW is the free energy of anything. They count identities. Stage 14
%     scores the same registers thermodynamically and applies the catalytic test through MispEndTol
%     and the extension-efficiency weighting; these two are deliberately blunter and deliberately
%     cheaper, and they stand in front of it.
%     Both quantities are logged for every candidate whether or not either gate is enforced.
%
%  8. 3'-TERMINAL STABILITY.  Free energy of the last EndWin bases at 37 C in 1 M Na+, WITH the
%     terminal and initiation rows, which is the scale Primer3 quotes Max 3' Stability on. An
%     over-stabilised 3' end lets a partial match hold long enough to prime. EndStability also
%     returns the stacking-only value, the 5' window and their difference, which is the internal
%     stability profile in free-energy form.
%
%  9. TEMPLATE ACCESSIBILITY.  The binding site must not be buried in template secondary structure.
%     The reverse primer has to invade a structured circle in cycle one and the forward primer has
%     to invade the nascent minus strand. dGopen = -R*T*ln(P_unpaired) is the free energy cost of
%     opening the site, from one RNAplfold fold computed before the search, so each candidate is an
%     O(1) lookup. Optional; the search runs with dGopen = 0 when ViennaRNA is absent.
%
% 10. MELTING TEMPERATURE, at the primer's OWN starting concentration. This is not optional in
%     LATE-PCR. Lowering a primer's concentration lowers its melting temperature, and designing an
%     asymmetric pair on nominal concentrations is precisely the failure mode LATE-PCR exists to
%     correct. PrimTemp builds dH and dS exactly as dGPrim does, forms the two-state melting
%     temperature for unequal strand concentrations, then applies the Owczarzy shift to 1/Tm.
%
% 11. HOMODIMER.  Priming competes with self-association. The threshold is not a constant: it is the
%     free energy at which the equilibrium homodimer would hold exactly phiDimer of the PRIMER
%     POOL, obtained by solving the mass balance at the stated primer and template concentrations.
%     The denominator is the primer pool and not the primer-template complex, which is where this
%     stage parts company with the mispriming stages. A homodimer competes for PRIMER. Referenced
%     to the primer pool the threshold does not move with template at all, and the only template
%     dependence left is the second-order one, that sequestering primer into template lowers the
%     free primer available to dimerise. Template STILL enters the mispriming thresholds, where it
%     belongs, because a misprime does compete for template.
%
% 12. 3'-ANCHORED SELF-DIMER.  Computed separately from stage 11 because a whole-primer search
%     returns the most stable complementarity, which is not the one the polymerase can extend. Only
%     a duplex reaching the primer's own 3' terminus can be extended, since Taq carries no 3'->5'
%     proofreading exonuclease to trim a dangling tail. Stage 11 therefore measures SEQUESTRATION
%     and this stage measures EXTENSIBILITY, and the second is the one that ruins a reaction,
%     because an extended dimer becomes a competing template with two primer sites and amplifies
%     exponentially from cycle one. Extensibility also needs a template: the polymerase copies the
%     partner strand, which supplies one only past the base paired to the terminus. DimerTemplateBeyond
%     is the number of partner nucleotides required 5' of that base, so a duplex whose partner ends at
%     the paired base, leaving nothing to copy, is not counted as extendable.
%     THE GATE CHANGES SHAPE WITH THE ENZYME. Declare a proofreading exonuclease and the premise
%     above fails: the enzyme trims the dangling 3' tail, so any register stable enough to persist
%     at the annealing temperature is extendable. dGEndMin then reads the whole-molecule minimum
%     from stage 11 rather than the anchored minimum from this one. Same constant, applied to the
%     quantity the enzyme can act on.
%
% 13. HAIRPIN, judged on MELTING TEMPERATURE.  A hairpin costs the reaction two unrelated ways.
%     It sequesters primer, lowering the effective concentration, which is an equilibrium governed
%     by the folded fraction at the annealing temperature. And if the stem pairs the primer's own
%     3'-terminal base it blocks extension outright. Free energy is the wrong currency for either,
%     because dG passes through zero at the structure's own melting point, so a dG cutoff evaluated
%     near the annealing temperature cannot fire for exactly the hairpins that matter. This stage
%     gates on a 3'-blocking hairpin reaching the annealing temperature and scores the folded
%     fraction otherwise.
%     The enumeration admits BULGED stems as well as continuous ones: a structure is an outer
%     helix, an optional bulge on either arm, an inner helix and a loop, which reduces to the
%     symmetric register when the bulge is zero. A bulge of one keeps the stack across the
%     insertion and a longer bulge does not, which is why bulge initiation is not monotonic in
%     length.
%     Loop stability uses the measured terminal mismatch at the closing pair and the sequence
%     bonuses for exceptionally stable triloops and tetraloops. Both matter more than their size
%     suggests: applying one constant terminal mismatch to every hairpin regardless of sequence,
%     and carrying no triloop bonus, misplaces a hairpin melting temperature by tens of degrees on
%     a stem that closes on a stable triloop. The register is found either way; only the energy is
%     wrong.
%
% 14. WHOLE-PLASMID MISPRIMING.  Every ungapped register of the primer against both strands that
%     carries a seed match is scored with mismatch nearest-neighbour parameters, with internal-loop
%     bridging for two-sided structures, and the cognate site excluded. Three things are then asked
%     of each anchored register. Its extension-weighted free energy, dGeff = dG - RT*ln(epsilon),
%     because occupancy says how often a primer sits at a site and says nothing about whether the
%     polymerase will add a nucleotide once it is there. Its contiguous run of identities. And the
%     stability of the continuous Watson-Crick run carrying the 3' terminus, against the measured
%     figure at which Rychlik found priming becoming detectable.
%
%     WHICH REGISTERS COUNT depends on the enzyme, through MispEndTol. Without proofreading only a
%     register reaching the primer 3' terminus to within one base can prime, and the stage sees
%     those alone. With proofreading the exonuclease trims any overhang back to the outermost
%     paired base, so every register with a pair is reachable and MispEndTol is widened to the
%     longest primer. The right edge of the scan is pinned at that outermost pair either way, which
%     is what makes the wide case the TRIMMED duplex rather than the free-edge optimum.
%
% 15. DISCRIMINATION MARGIN.  The smallest finite anchored margin must reach MispMargin. At zero
%     the stage asks only that the threshold be met and the mispriming target lives entirely in
%     phi. Raising it reinstates a second, hidden phi, because the threshold enters as
%     -RT*ln(phi*K/(1 + K*f)) and demanding m kcal/mol of headroom above it multiplies phi by
%     exp(-m/RT). AssertMispSatisfiable raises before that can silently reject every candidate.
%
% THE PAIRING STAGES
%
%  P1. LATE-PCR OFFSET.  The limiting primer must melt above the excess primer by dTrfmin to
%      dTrfmax degrees, each at its own starting concentration. Pierce et al. measured relative
%      amplification efficiency across offsets of -3, 0, +3, +5 and +7 as 0.85, 0.89, 0.89, 0.97 and
%      1.00, with 0 and +3 statistically indistinguishable. Efficiency then rises monotonically to
%      the highest offset they measured. The optimum shifts with primer ratio, and those
%      measurements were made at 20:1.
%
%  P2. AMPLICON MELTING TEMPERATURE.  Once the limiting primer is spent, the accumulating single
%      strand competes with the excess primer for the same template strand, so the rate and extent
%      of linear synthesis are set by TmA - TmX. Pierce et al. measured early slopes of 11.3, 7.6,
%      4.2 and 2.3 at differences of 13, 17, 20 and 23 degrees. Gated at TmAmaxGap and scored
%      towards TmAoptGap. TmAmplicon uses the %GC form rather than a nearest-neighbour sum, because
%      helix growth dominates initiation in a duplex of this length and no strand-concentration
%      effect survives.
%      BOTH SIDES OF THE GATE ARE EVALUATED IN PIERCE'S OWN CONVENTION, monovalent cation only,
%      because that is the convention TmAmaxGap and TmAoptGap were measured in and a threshold
%      means nothing apart from its convention. The engine's working temperatures carry magnesium
%      through Owczarzy and through the equivalent-salt TmA, and those stay in the log as TmA and
%      dTmAX; the gated quantity is dTmAXP, built from TmAP and TmXP.
%      A LONG AMPLICON IS INTRINSICALLY HARD HERE. Past a kilobase the -500/length term is dead and
%      TmA is set by GC and salt alone, so no primer choice moves it, while TmX is capped by Tmax.
%      The gap is pinned near its floor and the linear phase is slow. AssertTmAGapSatisfiable
%      refuses the run outright when the floor rises above Tmax.
%
%  P3. HETERODIMER, against a threshold solved from phiDimer at one temperature. The reference is
%      the LIMITING primer's pool, since that is the strand whose loss stops the exponential phase,
%      and the two primer balances decouple because each primer's own template duplex is the only
%      coupling term and it is subtracted out before the quadratic is formed. The threshold is
%      therefore independent of template to the same order as stage 11, and reduces to
%      -RT*ln(phiDimer/f0) in the dilute limit, within 0.01 kcal/mol of the exact root.
%
%  P4. 3'-ANCHORED CROSS-DIMER, tested in both orientations because either primer's 3' terminus can
%      extend. Changes shape under proofreading exactly as stage 12 does, and for the same reason.
%
%  P5. AMPLICON 3'-OH COLLAPSE, LOG ONLY.  The single strand that accumulates ends in a free 3'-OH,
%      and so does the limiting-primer strand. Either can fold back on its own body and be extended.
%      Each accepted pair is screened on both strands by ScaffoldCollapseCheck: a stem plus loop
%      scored by nearest neighbours, a persistence test on the fold's melting temperature against
%      AnnealingTemperature - HpTmMargin, and the extension efficiency of the folded end. The result
%      is a tier per strand (REJECT, FLAG, TRANSIENT, NONE) with the probability of an extendable
%      collapse, written to Pair_Log (cp* for the plus strand, cm* for the minus strand) and to
%      Candidate_Primers. It never rejects a pair and never enters Score.
%      The screen also scores the bimolecular state, in which the 3'-OH pairs with the same stem on a
%      second copy of the strand, at the strand's ceiling concentration: the excess primer for the
%      plus strand and the limiting primer for the minus strand. It reports a Taq channel and a
%      proofreading channel separately, with an optional blend, and a tier over all of them.
%
% ____________________________________________________________________________________________________________
% THE SPECIFICITY TARGET
%
% TWO tolerances, because there are two competitions and they do not share a denominator.
%
% phi governs MISPRIMING, stages 14 and 15. It is the tolerated molar ratio of unwanted to wanted
% complex AT ONE SITE, so the denominator is the primer-template duplex. A misprime competes for
% template; template therefore belongs inside the threshold, and a primer that holds its own site
% tightly earns a proportionally more permissive off-target allowance. The scan records of order
% 122 anchored registers per primer across both strands on a 7 kb circle, so the total off-target
% fraction is about 122*phi rather than phi.
%
% phiDimer governs DIMERS, stages 11 and P3. It is the tolerated fraction of the PRIMER POOL held
% as dimer, so the denominator is the primer. A dimer competes for primer and is indifferent to how
% much template is present. It is set looser than phi because this gate covers sequestration alone.
% Extensibility, the concern that actually ruins a reaction, is carried by dGEndMin at stages 12
% and P4. Declare a proofreading exonuclease and phiDimer is tightened to phi, because dGEndMin can
% no longer carry extensibility on an anchored register alone.
%
% ALL FOUR STILL SOLVE AT THE SAME TEMPERATURE, the reaction's own annealing step, and each derives
% the primer-template free energy its threshold rests on rather than accepting one from the caller.
% Evaluating each primer's structures at an offset below ITS OWN melting temperature is reasonable
% for a symmetric reaction and wrong for this one: LATE-PCR holds the two melting temperatures
% dTrfmin to dTrfmax apart by design, so the limiting primer would be screened several degrees
% hotter than the step it has to survive and its thresholds would not be commensurable with the
% excess primer's.
%
% The two contiguous-run gates sit outside both tolerances on purpose. They are combinatorial, they
% carry no concentration and no temperature, and a run of identities is a property of the two
% sequences and of nothing else. No enzyme and no buffer moves either of them.
%
% ____________________________________________________________________________________________________________
% RANKING
%
% Score, low is better, every term in a stated currency:
%       |dTm - centre of the offset window|                      degrees
%     + WtTmA    * ln(slope lost at this TmA - TmX gap), capped      e-folds of yield
%     - WtMargin * worst dimer headroom, capped                   per kcal/mol
%     - WtMisp   * worst ANCHORED mispriming headroom, capped     per kcal/mol
%     + WtAccess * (forward + reverse opening cost)               per kcal/mol
%     + Wt3End   * grippiest 3' end of the pair, capped           per kcal/mol
%     - WtProf   * worst internal-stability profile, floored      per kcal/mol
%     + WtFlank  * flank/(fFlex + rFlex), amplicon overhang         fraction of the window
%     + WtHairpin* worst folded fraction at the annealing step    dimensionless, 0 to 1
%     - Wt3Base  * worst 3'-terminal base, -log10(geometric mean epsilon over its three mispairs)
%
% The sign on the last term is a MINUS and that is not a slip. The statistic is -log10 of the
% geometric mean extension efficiency over the three mispairs the base can form, so a LARGER value
% means sharper discrimination against a mismatched terminus. Score is low-is-better, so better
% discrimination has to pull the score down.
%
% The contiguous-run lengths do NOT enter Score. They are integers counting identities, they have
% no currency the other terms are measured in, and the thermodynamic consequence of a long run is
% already priced by the mispriming headroom at stage 14. They gate and they are logged.
%
% Dimer free energies are clamped at zero before the headroom is measured. A positive dimer dG means
% no duplex forms, and without the clamp the margin would keep growing with how bad the least-bad
% non-dimer register happens to be, so the ranking would reward one non-structure over another.
%
% ____________________________________________________________________________________________________________
% WHAT IS NOT MODELLED
%
% Polymerase capture rate. Every occupancy here is an equilibrium duty cycle, and converting it to a
% per-cycle priming probability needs k_cap times the annealing time, which is not in this suite.
%
% Primer degradation by a proofreading enzyme across a long run.
%
% GAPPED ALIGNMENTS IN THE RUN GATES. Stages 6 and 7 count ungapped perfect runs only. A primer
% whose homology with an off-target site is interrupted by a single-base insertion carries two
% shorter runs by this measure and one long bulged duplex in reality. Stage 14 sees the bulged
% duplex, because MisprimerMismatch bridges internal loops; the run gates do not.
%
% DNA TANDEM MISMATCHES, the 2x2 internal loops that RNA folding calls int22. These are absent
% because they have not been measured. SantaLucia and Hicks state that DNA parameters exist only
% for tandem G.T and that other 2x2 sequences remain undetermined. RNA has a complete periodic
% table of them and it is NOT substituted here, because the two polymers diverge qualitatively even
% on simpler motifs: B-form DNA shows no preference between 3' and 5' dangling ends and carries two
% destabilising values, while A-form RNA strongly prefers 3' and has none.
%
% Sequence dependence of a BULGED base. The review reports A bulges are appreciably more stable
% than C bulges with G and T between, attributing it to unpublished work, and gives no parameters.
%
% HELIX INITIATION FOR A UNIMOLECULAR STRUCTURE. A hairpin stem is scored with dGMismatch, which
% returns a duplex free energy and so charges an initiation at each end. The review computes a
% hairpin as salt-corrected stacks plus the loop term with no initiation at all. Removing both ends
% overshoots Primer3 as far as keeping them undershoots it, so the three models genuinely disagree
% and the question is left open rather than tuned to match one of them. The residual is about
% 0.8 kcal/mol on a typical stem.
%
% ____________________________________________________________________________________________________________
% OUTPUTS
%
%   tab      accepted pairs, ranked best first, so row 1 is what aPCR_Claude offers as its dialog
%            default.
%   results  Forward_Log, Reverse_Log and Pair_Log. EVERY candidate ever built appears exactly once
%            in the matching log, and RejectStage names the first filter it failed, so an empty
%            shortlist is diagnosable rather than mysterious. Parameters carries the settings the
%            run actually used, assembled by PrimerSearchParameters.m.
%
% ____________________________________________________________________________________________________________
% REFERENCES
%
%   SantaLucia J Jr. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor
%       thermodynamics. Proc Natl Acad Sci USA 1998;95:1460-5.  doi:10.1073/pnas.95.4.1460
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
%   Bommarito S, Peyret N, SantaLucia J Jr. Thermodynamic parameters for DNA sequences with
%       dangling ends. Nucleic Acids Res 2000;28:1929-34.  doi:10.1093/nar/28.9.1929
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
%   Wetmur JG. DNA probes: applications of the principles of nucleic acid hybridization.
%       Crit Rev Biochem Mol Biol 1991;26:227-59.  doi:10.3109/10409239109114069
% ____________________________________________________________________________________________________________
%% ---- SETTINGS ----------------------------------------------------------------------------------
% Every constant, flag and threshold this function runs on, and the reason for each value, lives in
% PrimerSearchParameters.m. The call returns them AFTER the polymerase has rewritten the four it
% conditions and after the fold salt has been derived from the buffer, so the Parameters sheet
% describes the configuration that actually executed rather than the one that was typed.
[cfg, Parameters] = PrimerSearchParameters(input);

% Unpacked into local names in one block, so the body below reads in the names a user of the suite
% knows and so this list is a manifest of what the engine consumes. A setting the parameter file no
% longer supplies fails here, by name, before anything runs.
fFlex = cfg.fFlex; rFlex = cfg.rFlex;
fLmin = cfg.fLmin; fLmax = cfg.fLmax; rLmin = cfg.rLmin; rLmax = cfg.rLmax;
Tmin = cfg.Tmin; Tmax = cfg.Tmax;
dTrfmin = cfg.dTrfmin; dTrfmax = cfg.dTrfmax; dTsub = cfg.dTsub;
TmAmaxGap = cfg.TmAmaxGap; TmAoptGap = cfg.TmAoptGap; UseTmAEquiv = cfg.UseTmAEquiv;
TmAPierceGap = cfg.TmAPierceGap; TmAPierceSlope = cfg.TmAPierceSlope;
UseGCWindow = cfg.UseGCWindow; fGCmin = cfg.fGCmin; fGCmax = cfg.fGCmax;
MonoLim = cfg.MonoLim; DiRepLim = cfg.DiRepLim;
Force5primeGC = cfg.Force5primeGC; ForwardForce3primeAT = cfg.ForwardForce3primeAT;
ForwardAvoid3primeT = cfg.ForwardAvoid3primeT; ReverseForce3primeGC = cfg.ReverseForce3primeGC;
ForwardGCSkew = cfg.ForwardGCSkew; ReverseGCSkew = cfg.ReverseGCSkew;
EndWin = cfg.EndWin; UseEndStabilityGate = cfg.UseEndStabilityGate; dG3min = cfg.dG3min;
dGEndMin = cfg.dGEndMin; dGhpmin = cfg.dGhpmin; HpTmMargin = cfg.HpTmMargin;
CollapseLog = cfg.CollapseLog; CollapseWindow = cfg.CollapseWindow; CollapseMinRun = cfg.CollapseMinRun;
CollapseFlagEps = cfg.CollapseFlagEps; CollapseLBase = cfg.CollapseLBase; CollapseLPers = cfg.CollapseLPers;
CollapseBimolecular = cfg.CollapseBimolecular; CollapseExoBlend = cfg.CollapseExoBlend;
CollapseExoEndTol = cfg.CollapseExoEndTol; DimerTemplateBeyond = cfg.DimerTemplateBeyond;
phi = cfg.phi; phiDimer = cfg.phiDimer;
EndKmer = cfg.EndKmer; EndKmerMax = cfg.EndKmerMax;
Use3PrimeRunGate = cfg.Use3PrimeRunGate; MaxRun3 = cfg.MaxRun3;
UseInternalRunGate = cfg.UseInternalRunGate; MaxRunInt = cfg.MaxRunInt;
RunMinReport = cfg.RunMinReport; RunPoolFloor = cfg.RunPoolFloor;
MispMinRun = cfg.MispMinRun; MispSurveyRun = cfg.MispSurveyRun;
MispEndTol = cfg.MispEndTol; MispMargin = cfg.MispMargin;
MispEndWin = cfg.MispEndWin; Rychlik = cfg.Rychlik;
ProofExo = cfg.ProofExo; ProofExoClass = cfg.ProofExoClass;
UseExtensionEfficiency = cfg.UseExtensionEfficiency;
UseAccessibility = cfg.UseAccessibility; AccW = cfg.AccW; AccL = cfg.AccL;
AccSalt = cfg.AccSalt; AccSaltSrc = cfg.AccSaltSrc; AccMaxOpen = cfg.AccMaxOpen;
AccParamFile = cfg.AccParamFile;
AccCacheOn = cfg.AccCacheOn; AccCacheRoot = cfg.AccCacheRoot; AccCacheWrite = cfg.AccCacheWrite;
WtMargin = cfg.WtMargin; MarginCap = cfg.MarginCap; WtMisp = cfg.WtMisp; MispCap = cfg.MispCap;
WtAccess = cfg.WtAccess; WtTmA = cfg.WtTmA; TmACap = cfg.TmACap;
WtFlank = cfg.WtFlank; WtHairpin = cfg.WtHairpin;
Wt3End = cfg.Wt3End; WtProf = cfg.WtProf; Wt3Base = cfg.Wt3Base;
WtEns = cfg.WtEns; EnsCap = cfg.EnsCap;
Verbose = cfg.Verbose; ProgressEvery = cfg.ProgressEvery;

RunClock = tic;
if Verbose
    fprintf('\n%s\n', repmat('=', 1, 96));
    fprintf('LOCAL PRIMER SEARCH ENGINE   started %s\n', ...
        string(datetime('now', 'Format', 'dd-MM-yyyy HH:mm:ss')));
    fprintf('%s\n', repmat('=', 1, 96));
end

%% ---- DEPENDENCY PREFLIGHT ----------------------------------------------------------------------
% MATLAB dispatches on the file name and ignores the name on the function declaration line, so any
% residual mismatch between the two surfaces here rather than half-way through a scan.
% Existence alone is not enough. Several of these functions carry signatures that have changed, and
% an older copy sitting earlier on the path is found instead and then fails deep inside a scan.
% A two-input copy of realcard found ahead of the three-input one that dGMispTarget calls surfaces
% at the first candidate scored and not at startup, which is an hour of scanning spent to learn
% that a file is duplicated. nargin and nargout report the declared counts, so a stale file is
% caught before anything runs. Any signature hit means two copies sit on the path.
%
% ONE TABLE, NOT TWO. A bare existence list (dep) alongside a separate
% signature list (sigreq), and the signature loop opened with "if isempty(which(nm)), continue".
% Anything present in sigreq but absent from dep therefore passed the preflight BY BEING MISSING,
% which inverts the whole purpose of the block. dGMispTarget and dGMispOffTarget appeared in
% neither list, and MisprimerMismatch calls both before it scores a single register.
%
% Columns: name, EXACT declared inputs, minimum declared outputs. See the test below for why
% inputs cannot be a minimum.
sigreq = {'PrimTemp', 7, 3; 'dGPrim', 6, 1; 'dGMismatch', 10, 1; 'dGmaxMismatch', 5, 1; ...
    'dGLoop', 4, 1; 'NucleoRun', 3, 3; 'DimerScan', 8, 4; 'LooperMismatch', 10, 7; ...
    'TmAmplicon', 5, 2; 'FreeMgSpeciated', 3, 2; ...
    ... % FreeMgSpeciated holds the Mg-dNTP chelation root that NNParams and TmAmplicon both
    ... % take their free magnesium from. A missing copy leaves TmAmplicon's useEquiv branch
    ... % dead and NNParams unable to build its salt handles at all, so it is checked here
    ... % rather than discovered at the first anonymous-function call three files away.
    ... % H8: the measured motif tables and their accessors. A stale copy of any of these changes
    ... % every hairpin energy in the run without erroring, which is exactly what this block exists
    ... % to prevent.
    'NNParams', 1, 1; 'dGDangle', 5, 4; 'dGTerminalMismatch', 6, 4; 'dGBulge', 5, 2; ...
    'BestSubRegister', 10, 1; 'ExtensionEfficiency', 3, 1; ...
    'ExtensionEfficiencyTable', 2, 1; ...
    'DimerMeltTm', 9, 2; ...
    ... % ExtensionEfficiencyTable owns the session's choice of efficiency table. A missing
    ... % copy makes ExtensionEfficiency fall back to the built-in Huang values while the
    ... % Parameters sheet still reports whatever EffTableFile was set to, which is the one
    ... % failure mode that would put a wrong table behind a right-looking record.

    ... % NOTE: BestSubRegister gained STRUCT FIELDS rather than arguments, so nargin and nargout
    ... % cannot see the difference. An older copy on the path passes this check and then errors on
    ... % sub.iTrunRight inside SeedScan. If that happens, the copy is stale, not the caller.
    'HomodimerizerMismatch', 11, 6; 'HeterodimerizerMismatch', 13, 5; ...
    'MisprimerMismatch', 17, 6; 'InternalLoopBridge', 9, 2; ...
    'dGMispTarget', 5, 1; 'dGMispOffTarget', 5, 1; 'dGHomo', 5, 1; 'dGHetero', 7, 1; ...
    'realcard', 3, 1; 'cardano', 4, 1};
missing = {}; stale = {};
if Verbose
    fprintf('[1/8] Dependencies   checking %d functions for presence and signature ... ', ...
        size(sigreq, 1));
end
for s = 1:size(sigreq, 1)
    nm = sigreq{s,1};
    if isempty(which(nm))
        missing{end+1} = nm; %#ok<AGROW>
        continue
    end
    try
        gotin = nargin(nm); gotout = nargout(nm);
    catch
        continue    % p-coded or built in, so there is nothing to compare against
    end
    % INPUTS ARE COMPARED EXACTLY, OUTPUTS AS A MINIMUM.
    % A minimum on inputs has a blind spot, and it is not hypothetical: HeterodimerizerMismatch
    % went from 13 declared inputs to 12 when dGf and dGr were removed. A stale 13-input copy
    % satisfies 13 >= 12, passes the check, and then receives the caller's 12 arguments
    % POSITIONALLY, so f0 lands in dGf, phi lands in f0, and every argument after the second is
    % shifted by two. That runs to completion and returns nonsense. Any signature that SHRINKS is
    % invisible to a minimum test, so the declared input count must match the table exactly.
    % Outputs stay a minimum, because a function declaring more outputs than the caller requests is
    % harmless; MATLAB simply does not compute them.
    % nargin returns a negative number for a varargin function, which is why the sign is guarded.
    if ((gotin >= 0)&&(gotin ~= sigreq{s,2}))||((gotout >= 0)&&(gotout < sigreq{s,3}))
        stale{end+1} = sprintf('%s declares %d in / %d out, this engine expects %d in / >= %d out  [%s]', ...
            nm, gotin, gotout, sigreq{s,2}, sigreq{s,3}, which(nm)); %#ok<AGROW>
    end
end
% Bioinformatics Toolbox, reached from DimerScan, EndStability and MisprimerMismatch.
btb = {'seqcomplement', 'seqrcomplement', 'swalign'};
btbmissing = btb(cellfun(@(f) isempty(which(f)), btb));
if ~isempty(btbmissing)
    missing = [missing, btbmissing];
end
if ~isempty(missing)
    error('LocalPrimerSearchEngine:MissingDependency', ...
        ['These functions are not on the MATLAB path: %s.\nRun "which -all <name>" to check for ' ...
        'shadowing. Bioinformatics Toolbox supplies seqcomplement, seqrcomplement and swalign.'], ...
        strjoin(missing, ', '));
end
if ~isempty(stale)
    error('LocalPrimerSearchEngine:StaleDependency', ...
        ['An out-of-date copy of at least one dependency is first on the MATLAB path:\n   %s\n' ...
        'Run "which -all <name>" to find the duplicate, then delete or reorder it.'], ...
        strjoin(stale, sprintf('\n   ')));
end

% The collapse screen is log-only, so a missing copy switches it off with a warning and the search runs.
if CollapseLog&&(isempty(which('ScaffoldCollapseCheck'))||isempty(which('BimolecularTm')))
    warning('LocalPrimerSearchEngine:NoCollapseCheck', ...
        'ScaffoldCollapseCheck.m or BimolecularTm.m is not on the path; the amplicon collapse log is switched off.');
    CollapseLog = false;
end

if Verbose, fprintf('all present.\n'); end

%% ---- SEQUENCE SET-UP ---------------------------------------------------------------------------
PlasmidSeq = upper(char(input.PlasmidSequence));
if any(~ismember(PlasmidSeq, 'ACGT'))
    bad = find(~ismember(PlasmidSeq, 'ACGT'), 1);
    error('LocalPrimerSearchEngine:NonACGT', ...
        ['PlasmidSequence carries a character outside {A,C,G,T} at index %d (''%s''). Nearest-neighbour ' ...
        'tables, strfind and count are defined only on unambiguous DNA; clean the sequence at load ' ...
        'time in CollectaPCRData_Claude.m.'], bad, PlasmidSeq(bad));
end
PlasLen = strlength(PlasmidSeq);
RevComPlas = seqrcomplement(PlasmidSeq);
% Doubled strings, built once, so the k-mer uniqueness test sees the circle rather than a cut line.
PlasDbl = [PlasmidSeq PlasmidSeq];
RevDbl  = [RevComPlas RevComPlas];

% IS THE AMPLICON GATE REACHABLE AT ALL?
% The amplicon melting temperature is fixed by the target. Past about a kilobase the -500/length
% term is dead and TmA is set by GC and salt alone, so no choice of primer moves it. The excess
% primer melting temperature is capped by Tmax. The P2 gate therefore demands
%       tmf >= TmA - TmAmaxGap
% and if that floor sits above Tmax the gate is empty by arithmetic, for every sequence, and the
% engine will grind through the whole scan to return nothing. A gate wanting tmf >= 65.12 against
% a Tmax of 65 is empty by 0.12 degrees and takes the better part of an hour to say so. Check it up
% front instead.
AssertTmAGapSatisfiable(PlasDbl, input, Tmax, TmAmaxGap, dTrfmin, fLmax);

if Verbose
    nFwdTry = (fFlex + 1)*(fLmax - fLmin + 1);
    nRevTry = (rFlex + 1)*(rLmax - rLmin + 1);
    fprintf('[2/8] Template       %d nt circular, GC %.1f%%\n', PlasLen, ...
        100*nnz((PlasmidSeq == 'G')|(PlasmidSeq == 'C'))/PlasLen);
    fprintf('      Scaffold       %d to %d  (%d nt)\n', input.ScaffoldStart, input.ScaffoldEnd, ...
        input.ScaffoldEnd - input.ScaffoldStart + 1);
    fprintf('      Forward window %d to %d, lengths %d-%d   -> %d candidates\n', ...
        input.ScaffoldStart - fFlex, input.ScaffoldStart, fLmin, fLmax, nFwdTry);
    fprintf('      Reverse window %d to %d, lengths %d-%d   -> %d candidates\n', ...
        input.ScaffoldEnd, input.ScaffoldEnd + rFlex, rLmin, rLmax, nRevTry);
    fprintf('      Worst case     %d pairs if every candidate survives\n', nFwdTry*nRevTry);
    OnOffWords = {'off', 'ON'};       % indexing a cell LITERAL is invalid in MATLAB, so name it
    onoff = @(b) OnOffWords{double(logical(b)) + 1};
    fprintf('[3/8] Active gates\n');
    fprintf('      Tm window          %.1f to %.1f C,  excess-limiting offset %.1f to %.1f C\n', ...
        Tmin, Tmax, dTrfmin, dTrfmax);
    fprintf('      Amplicon Tm gap    gate <= %.1f C, ranked towards %.1f C  (monovalent convention)\n', ...
        TmAmaxGap, TmAoptGap);
    fprintf('      GC window          %s  (%.2f to %.2f)\n', onoff(UseGCWindow), fGCmin, fGCmax);
    fprintf('      GC skew  fwd %s / rev %s        5'' G-or-C clamp %s\n', ...
        onoff(ForwardGCSkew), onoff(ReverseGCSkew), onoff(Force5primeGC));
    fprintf('      3'' base clamps     fwd A-or-T %s, fwd avoid-T %s, rev G-or-C %s\n', ...
        onoff(ForwardForce3primeAT), onoff(ForwardAvoid3primeT), onoff(ReverseForce3primeGC));
    fprintf('      3'' end stability   stage-6 gate %s at %.2f kcal/mol;  Rychlik cap %s at %.1f\n', ...
        onoff(UseEndStabilityGate), dG3min, onoff(Rychlik.On), Rychlik.Max);
    % The Rychlik cap does not move with the enzyme, so the class is printed beside it where a
    % reader is most likely to wonder whether it should have.
    fprintf('      Enzyme             3''->5'' exonuclease %s;  extension weighting %s\n', ...
        ProofExoClass, onoff(UseExtensionEfficiency));
    % The table is a persistent inside ExtensionEfficiency, so printing it here is one of the
    % three places the choice is announced. See ExtensionEfficiencyTable.m.
    if isempty(cfg.EffTableFile)
        effTxt = 'built-in Huang et al. 1992';
    else
        [~, effN, effE] = fileparts(cfg.EffTableFile); effTxt = [effN effE];
    end
    fprintf('      Efficiency table   %s;  3''-base spread %.3f decades, Wt3Base %.3f\n', ...
        effTxt, cfg.Base3Spread, Wt3Base);
    fprintf('      Hairpin            dG %.2f kcal/mol, blocking margin %.1f C\n', dGhpmin, HpTmMargin);
    fprintf('      Amplicon collapse  %s (log only), window %d nt, stem run >= %d nt\n', ...
        onoff(CollapseLog), CollapseWindow, CollapseMinRun);
    if isfinite(CollapseExoBlend), blendTxt = sprintf('%.2f', CollapseExoBlend); else, blendTxt = 'not set'; end
    fprintf('      Collapse channels  bimolecular %s, Exo flap <= %d nt, exonuclease blend %s\n', ...
        onoff(CollapseBimolecular), min(CollapseExoEndTol, CollapseWindow - CollapseMinRun), blendTxt);
    fprintf('      Dimer extension    template beyond the paired 3'' end >= %d nt\n', DimerTemplateBeyond);
    fprintf('      Contiguous runs    3'' terminal %s at >= %d nt, internal %s at >= %d nt\n', ...
        onoff(Use3PrimeRunGate), MaxRun3, onoff(UseInternalRunGate), MaxRunInt);
    fprintf('      Accessibility      %s\n', onoff(UseAccessibility));
end

%% ---- CONTIGUOUS-RUN MAPS -----------------------------------------------------------------------
% Built once for the whole circle, then read in O(1) per candidate by stages 6 and 7. See the
% header for what the maps hold and for why one plus-strand table answers the question for a
% forward primer and for a reverse primer alike.
% The length list carries L and L - 1 for every primer length in the search, because stage 6 reads
% the full window and stage 7 reads the same window one base shorter, with the 3'-terminal base
% dropped.
RunKmax = max(fLmax, rLmax);
RunLens = unique([fLmin:fLmax, rLmin:rLmax, (fLmin-1):(fLmax-1), (rLmin-1):(rLmax-1)]);
RunLens = RunLens(RunLens >= RunMinReport);
if Verbose
    fprintf('[4/8] Run maps       k = %d to %d over %d positions, both strands ... ', ...
        RunMinReport, RunKmax, PlasLen);
end
MapClock = tic;
NonUniq = BuildNonUniqueMaps(PlasmidSeq, RevComPlas, RunMinReport, RunKmax);
RunMap  = BuildRunLengthMaps(NonUniq, RunLens, RunMinReport, RunKmax, PlasLen);
if Verbose, fprintf('done in %.1f s.\n', toc(MapClock)); end

% ARE THE TWO RUN GATES SATISFIABLE ON THIS TEMPLATE AT ALL?
% The gates are a joint property of the primer and the TEMPLATE LENGTH, and the dependence is
% exponential, so a setting that is free on a short construct can admit nobody on a plasmid. The
% pool is measured here, with no other screen in front of it, because a gate that admits nobody
% should cost a second and not the whole scan.
AssertRunGatesSatisfiable(RunMap, NonUniq, input, PlasLen, fFlex, rFlex, fLmin, fLmax, ...
    rLmin, rLmax, Use3PrimeRunGate, MaxRun3, UseInternalRunGate, MaxRunInt, RunPoolFloor, Verbose);

% Lowest annealing temperature the reaction reaches. A touchdown recipe ramps down to a constant,
% and the lowest temperature is the permissive one for a hairpin, so that is the one to gate on.
% Resolved here, ahead of the fold, because the fold and every structural screen below are judged
% at the temperature the reaction anneals at. A per-primer offset tm - dTsub is not that
% temperature: it tracks each candidate's own melting temperature, so a 68.9 C reverse primer was
% screened at 63.9 C while the reaction anneals at 58.0 C, and the template was folded at 53.5 C.
% MispTemp returns this number when it is a real temperature and falls back to tm - dTsub when
% AnnealLow cannot resolve one, so a run with no usable programme behaves as it did before.
TAnneal = AnnealLow(input);
% The collapse screen is judged at the annealing temperature, so it needs a resolved one. It is
% log-only, so without one it switches off and the search carries on.
if CollapseLog&&~isfinite(TAnneal)
    warning('LocalPrimerSearchEngine:NoCollapseTemperature', ...
        'No annealing temperature could be resolved; the amplicon collapse log is switched off.');
    CollapseLog = false;
end

%% ---- FOLD THE TEMPLATE ONCE, IF ASKED ----------------------------------------------------------
% RNAplfold returns the unpaired probability of every window in a single pass, so this is one
% external call before the search and an O(1) lookup per candidate inside it. The plus strand is
% what the reverse primer must invade; the minus strand is what the forward primer must invade
% once it exists.
Rgas = 8.314E-3/4.184;          % kcal/mol/K, the same constant the rest of the suite uses
AccOn = UseAccessibility;
LunpPlus = []; MapPlus = []; LunpMinus = []; MapMinus = [];
if ~AccOn && Verbose
    fprintf('[5/8] Accessibility  switched off; every candidate scored with dGopen = 0.\n');
end
if AccOn && isempty(which('DNAAccessibility'))
    warning('LocalPrimerSearchEngine:NoAccessibility', ...
        'DNAAccessibility.m is not on the path; continuing with dGopen = 0 for every candidate.');
    AccOn = false;
end
if AccOn
    AccTemp = MispTemp(TAnneal, Tmin, dTsub); % the temperature the reaction anneals at
    if Verbose
        fprintf(['[5/8] Folding template with RNAplfold at %.1f C, %.3f M, W %d, L %d.\n' ...
            '      One external pass over both strands of a %d nt circle. This is the slow ' ...
            'step and it prints nothing until it returns.\n'], ...
            AccTemp, AccSalt, AccW, AccL, PlasLen);
        fprintf('      Fold salt %.4g M, from %s.\n', AccSalt, AccSaltSrc);
    end
    FoldClock = tic;
    aopts = struct('W', AccW, 'L', AccL, 'uMax', max(fLmax, rLmax) + 2, 'Temp', AccTemp, ...
        'Salt', AccSalt, 'Circular', true, 'NucLen', 12, 'Header', 'Plus', ...
        'ParamFile', AccParamFile);
    dummy = struct('Label', {'calibration'}, 'Start', {1}, 'Stop', {min(rLmax, PlasLen)});
    % ---- cache lookup ----
    AccKey = AccCacheKey(PlasmidSeq, aopts, AccParamFile);
    AccFile = '';
    if AccCacheOn
        if isfield(input, 'PlasmidFileName'), plasName = input.PlasmidFileName; else, plasName = ''; end
        AccFile = AccCachePath(AccCacheRoot, plasName, AccKey);
        [AccHit, LunpPlus, MapPlus, LunpMinus, MapMinus] = AccCacheLoad(AccFile, AccKey, Verbose);
    else
        AccHit = false;
    end
    if AccHit
        if Verbose
            fprintf('      Fold loaded from cache in %.1f s. No RNAplfold run.\n', toc(FoldClock));
        end
    else
    try
        [~, LunpPlus,  MapPlus]  = DNAAccessibility(PlasmidSeq, dummy, aopts);
        aopts.Header = 'Minus';
        [~, LunpMinus, MapMinus] = DNAAccessibility(RevComPlas, dummy, aopts);
        if Verbose
            fprintf('      Fold returned in %.1f s.\n', toc(FoldClock));
        end
        if AccCacheOn && AccCacheWrite
            AccCacheStore(AccFile, AccKey, LunpPlus, MapPlus, LunpMinus, MapMinus, PlasLen, Verbose);
        end
    catch ME
        warning('LocalPrimerSearchEngine:AccessibilityFailed', ...
            'RNAplfold step failed (%s); continuing with dGopen = 0.', ME.message);
        AccOn = false;
    end
    end
end
% The opening probability is read off the fold, so the free energy it is converted into has to use
% the temperature the fold ran at. Both come from one MispTemp call, so the two cannot drift apart.
AccTK = MispTemp(TAnneal, Tmin, dTsub) + 273.15;
% THE RESOLVED TEMPERATURE GOES ON THE PARAMETERS SHEET.
% input.AnnealingTemperature carries a temperature OR a mode flag OR a touchdown recipe, so the
% sheet renders it through AnnealText and no reader can take a number out of that row. Every
% downstream analysis needing the temperature the run ACTUALLY annealed at was therefore guessing,
% and CalibrateAggregateFlux was guessing 58. This row is the resolved number, written here rather
% than in PrimerSearchParameters because AnnealLow needs the input struct and runs on this side.
Parameters(end+1, :) = {"Thermal programme (resolved)", "TAnneal", ...
    string(sprintf('%.4g', TAnneal)), "C", ...
    "Lowest annealing temperature the programme reaches. Resolved from input.AnnealingTemperature, which may be a temperature, a mode flag or a touchdown recipe"};

%% ---- OUTPUT CONTAINERS -------------------------------------------------------------------------
% All three logs exist from the outset, so an empty search does not leave results.Reverse_Log
% unassigned and break writetable in aPCR_Claude.m.
% Candidate_Primers and Parameters ride in the SAME struct as the logs so the caller writes one
% workbook instead of two files. Parameters is filled here rather than at the end because the
% function returns early on several paths, and a run that found nothing is exactly when the reader
% most wants to see what it was asked to find.
tab = table();
results = struct('Forward_Log', table(), 'Reverse_Log', table(), 'Pair_Log', table(), ...
    'Candidate_Primers', table(), 'Parameters', Parameters);
h = waitbar(0, 'Running...');
waitbarcloser = onCleanup(@() CloseWaitbar(h));
pat = [{'G'},{'C'}];
TotalForwardSearchCount = (fFlex + 1)*(fLmax - fLmin + 1);
TotalReverseSearchCount = (rFlex + 1)*(rLmax - rLmin + 1);
ForwardSearchCounter = 0; ReverseSearchCounter = 0; PairSearchCounter = 0;

%% ================================================================== FORWARD PRIMER SEARCH =======
if Verbose, fprintf('[6/8] Forward scan   %d candidates\n', TotalForwardSearchCount); end
ScanClock = tic;
ForwardLog = table();
ForwardCandidate = table();
for i = (input.ScaffoldStart - fFlex):(input.ScaffoldStart)
    shiftplasmid = circshift(PlasmidSeq, -(i-1));
    for j = fLmin: fLmax
        primj = shiftplasmid(1: j);
        fPos = mod(i-1, PlasLen) + 1;
        record = struct('fPos', fPos, ...
            'flen', j, 'fseq', string(primj), 'fAT3', "", 'fGC5', "", ...
            'fGC', NaN,'dflocgc35', NaN, 'fnumocc', NaN, 'fdirep', NaN, ...
            'fKmer3', "", 'fKmerN', NaN, ...
            'fRun3', NaN, 'fRunInt', NaN, 'fRunIntOff', NaN, 'fRunIntSeq', "", ...
            'dGf3end', NaN, ...
            'dGf3stack', NaN, 'dGf5stack', NaN, 'fProf', NaN, ...
            'pfOpen', NaN, 'dGfOpen', NaN, 'fTm', NaN,...
            'dGf', NaN, 'dGffmin', NaN, 'fAlignf', "" , 'dGff', NaN, 'dGffend', NaN, 'ffTm', NaN, ...
            'fstemalign', "", 'dGfhp', NaN, 'fhpTm', NaN, 'flockTm', NaN, 'ffold', NaN, ...
            'fBase3', "", 'fEps3', NaN, ...
            'dGfTcmin', NaN, 'dGfTc', NaN, 'dGfTmin', NaN, 'dGfT', NaN, ...
            'fTcMargin', NaN, 'fTMargin', NaN, ...
            'dGfTcAnch', NaN, 'dGfTAnch', NaN, 'dGfTcEnd', NaN, 'dGfTEnd', NaN, ...
            'dGfTcAnchEff', NaN, 'dGfTAnchEff', NaN, 'fTcEps', NaN, 'fTEps', NaN, ...
            'fTcEpsPair', "", 'fTEpsPair', "", ...
            'fTcAnchMargin', NaN, 'fTAnchMargin', NaN, ...
            'fTcEnsShift', NaN, 'fTEnsShift', NaN, ...
            'fMispPosC', NaN, 'fMispPosT', NaN, 'fMispEndRnC', NaN, 'fMispEndRnT', NaN, ...
            'fMispAlignC', "", 'fMispAlignT', "", ...
            'fSurvC', NaN, 'fSurvT', NaN, 'fAnchC', NaN, 'fAnchT', NaN, ...
            ... % Which strands the mispriming scan actually reached. Scan T is skipped once
            ... % scan C has rejected, and every count above uses 0 for 'nothing found', so
            ... % without these a skipped scan reads exactly like a clean one.
            'fScanC', false, 'fScanT', false, ...
            ... % Aggregate off-target flux, MEASUREMENT ONLY, nothing gates on these.
            ... % Ens is the ensemble free energy over every anchored register and EnsRun the
            ... % same sum restricted to registers carrying MispMinRun identities. Read each
            ... % against the matching AnchEff to see how far the ensemble sits from the
            ... % single worst site, and against its count to know how many sites made it.
            'dGfTcAnchEns', NaN, 'dGfTcAnchEnsRun', NaN, 'fAnchRunC', NaN, ...
            'dGfTAnchEns',  NaN, 'dGfTAnchEnsRun',  NaN, 'fAnchRunT', NaN, ...
            'RejectStage', "Terminal_Bases", 'Status', "Eliminated");
        record.fAT3 = string(primj(end)); record.fGC5 = string(primj(1));
        record.fBase3 = string(primj(end));
        record.fEps3  = TerminalBaseDiscrimination(primj(end));
        % --- 1. terminal bases ---
        fGC5ok = ~Force5primeGC||((primj(1) == 'G')||(primj(1) == 'C'));
        fAT3ok = (~ForwardForce3primeAT||((primj(end) == 'A')||(primj(end) == 'T')))&& ...
            (~ForwardAvoid3primeT||(primj(end) ~= 'T'));
        if fGC5ok&&fAT3ok
            record.RejectStage = "GC_Content";
            % --- 2. global GC ---
            GC = count(primj, pat); gcf = GC/strlength(primj);
            record.fGC = gcf;
            if ~UseGCWindow||((gcf > fGCmin)&&(gcf <= fGCmax))
                record.RejectStage = "Local_GC_Skew";
                % --- 3. half-primer GC asymmetry ---
                locgc3 = count(primj((j-floor(j/2)+1):j),pat)/floor(j/2);
                locgc5 = count(primj(1:floor(j/2)),pat)/floor(j/2);
                record.dflocgc35 = locgc3 - locgc5;
                if ~ForwardGCSkew||(locgc3 < locgc5)
                    record.RejectStage = "Low_Complexity";
                    % --- 4. low complexity ---
                    [fnumocc, fdirep, check] = NucleoRun(primj, MonoLim, DiRepLim);
                    record.fnumocc = fnumocc; record.fdirep = fdirep;
                    if check == true
                        record.RejectStage = "Three_Prime_Kmer";
                        % --- 5. 3'-terminal k-mer uniqueness ---
                        fKmer = primj(end - EndKmer + 1:end);
                        fKmerN = KmerCount(fKmer, PlasDbl, RevDbl, PlasLen);
                        record.fKmer3 = string(fKmer); record.fKmerN = fKmerN;
                        % --- 6 and 7. contiguous non-cognate runs, both O(1) from the maps ---
                        % The forward primer occupies plus(fPos : fPos+j-1) and reads 5' to 3' with
                        % the plus strand, so its 3' terminus sits at the window's HIGH end.
                        [record.fRun3, record.fRunInt, record.fRunIntOff, record.fRunIntSeq] = ...
                            RunLengths(RunMap, NonUniq, primj, fPos, j, PlasLen, true);
                        % Three combinatorial gates share one branch, and a candidate is charged to
                        % the first of them it fails, which is what the cascade order means. All
                        % three quantities are logged whether or not any of them is enforced.
                        if fKmerN > EndKmerMax
                            fRunStage = "Three_Prime_Kmer";
                        elseif Use3PrimeRunGate&&(record.fRun3 >= MaxRun3)
                            fRunStage = "Three_Prime_Run";
                        elseif UseInternalRunGate&&(record.fRunInt >= MaxRunInt)
                            fRunStage = "Internal_Run";
                        else
                            fRunStage = "";
                        end
                        record.RejectStage = fRunStage;
                        if strlength(fRunStage) == 0
                            record.RejectStage = "Three_Prime_Stability";
                            % --- 8. 3'-terminal stability, published scale ---
                            % Logged and scored always; gated only when UseEndStabilityGate is on.
                            [dGf3end, dGf3stack, dGf5stack, fProf] = EndStability(primj, EndWin);
                            record.dGf3end = dGf3end;
                            record.dGf3stack = dGf3stack;
                            record.dGf5stack = dGf5stack;
                            record.fProf = fProf;
                            if (~UseEndStabilityGate)||(dGf3end >= dG3min)
                                record.RejectStage = "Site_Accessibility";
                                % --- 9. accessibility of the MINUS strand site ---
                                % The forward primer is identical to the plus strand, so it binds the
                                % minus strand. A plus coordinate p maps to minus coordinate
                                % PlasLen - p + 1, and the footprint's high end is at fPos.
                                [pfOpen, dGfOpen] = OpeningCost(AccOn, LunpMinus, MapMinus, ...
                                    PlasLen - fPos + 1, j, Rgas, AccTK);
                                record.pfOpen = pfOpen; record.dGfOpen = dGfOpen;
                                if dGfOpen <= AccMaxOpen
                                    record.RejectStage = "Melting_Temperature";
                                    % --- 10. melting temperature at its own concentration ---
                                    [tmf,~,~] = PrimTemp(primj, input.ForwardPrimerConcentration, input.PlasmidConcentration, ...
                                        input.MonovalentCationConcentration, input.BivalentCationConcentration, input.dNTPsConcentration);
                                    record.fTm = tmf;
                                    if isreal(tmf)&&isfinite(tmf)&&(tmf >= Tmin)&&(tmf <= (Tmax))
                                        record.RejectStage = "Homodimer";
                                        % --- 11. homodimer ---
                                        [dGf, dGffmin, fAlignf, dGff, check, dGffend] = HomodimerizerMismatch(primj, phiDimer, input.ForwardPrimerConcentration, ...
                                            input.PlasmidConcentration, input.MonovalentCationConcentration, ...
                                            input.BivalentCationConcentration, input.dNTPsConcentration, tmf, dTsub, TAnneal, DimerTemplateBeyond);
                                        record.dGf = dGf; record.dGffmin = dGffmin; record.fAlignf = AlignRow(fAlignf); record.dGff = dGff;
                                        record.ffTm = DimerMeltTm(primj, primj, ...
                                            input.ForwardPrimerConcentration, [], ...
                                            input.MonovalentCationConcentration, ...
                                            input.BivalentCationConcentration, ...
                                            input.dNTPsConcentration, MispTemp(TAnneal, tmf, dTsub));
                                        if check == true
                                            record.RejectStage = "Three_Prime_Self_Dimer";
                                            % --- 12. 3'-anchored self-dimer ---
                                            % dGffend is the sixth output of HomodimerizerMismatch: the minimum over
                                            % every ungapped register in which the primer 3' terminus is PAIRED, found
                                            % by DimerScan's exhaustive enumeration.
                                            record.dGffend = dGffend;
                                            % With proofreading the gate reads the WHOLE-MOLECULE
                                            % minimum instead, because the exonuclease trims an
                                            % unpaired 3' tail and makes an internal register
                                            % extendable. Both quantities stay in the log, so the
                                            % run is reconstructable either way.
                                            if ProofExo, ffGate = dGff; else, ffGate = dGffend; end
                                            if ffGate >= dGEndMin
                                                record.RejectStage = "Hairpin";
                                                % --- 13. hairpin, judged on melting temperature ---
                                                [fstemalign, check, fdGhp, fhpTm, flockTm, ffold] = ...
                                                    LooperMismatch(primj, input.MonovalentCationConcentration, ...
                                                    input.BivalentCationConcentration, input.dNTPsConcentration, ...
                                                    MispTemp(TAnneal, tmf, dTsub), dGhpmin, TAnneal, HpTmMargin);
                                                record.fstemalign = AlignRow(fstemalign);
                                                record.dGfhp = fdGhp; record.fhpTm = fhpTm;
                                                record.flockTm = flockTm; record.ffold = ffold;
                                                if check == true
                                                    record.RejectStage = "Mispriming";
                                                    % --- 14. whole-plasmid mispriming ---
                                                    [dGfTcmin, dGfTc, dGfTmin, dGfT, check, fmisp] = MisprimerMismatch(shiftplasmid, primj, phi, ...
                                                        input.ForwardPrimerConcentration, input.PlasmidConcentration, ...
                                                        input.MonovalentCationConcentration, input.BivalentCationConcentration, ...
                                                        input.dNTPsConcentration, tmf, dTsub, MispMinRun, MispEndTol, MispEndWin, ...
                                                        UseExtensionEfficiency, Rychlik, TAnneal, MispSurveyRun);
                                                    record.dGfTcmin = dGfTcmin; record.dGfTc = dGfTc;
                                                    record.dGfTmin = dGfTmin; record.dGfT = dGfT;
                                                    % Unanchored global optimum, kept as a diagnostic only.
                                                    record.fTcMargin = dGfTc - dGfTcmin;
                                                    record.fTMargin  = dGfT  - dGfTmin;
                                                    % 3'-anchored optimum. This is what stage 15 and Score consume.
                                                    record.dGfTcAnch = fmisp.dGpTcAnch; record.dGfTAnch = fmisp.dGpTAnch;
                                                    record.dGfTcAnchEff = fmisp.dGpTcAnchEff;
                                                    record.dGfTAnchEff  = fmisp.dGpTAnchEff;
                                                    record.fTcEps = fmisp.dGpTcEps; record.fTEps = fmisp.dGpTEps;
                                                    record.fTcEpsPair = fmisp.dGpTcEpsPair;
                                                    record.fTEpsPair  = fmisp.dGpTEpsPair;
                                                    record.dGfTcEnd  = fmisp.dGpTcEnd;  record.dGfTEnd  = fmisp.dGpTEnd;
                                                    record.fMispPosC = fmisp.dGpTcPos;  record.fMispPosT = fmisp.dGpTPos;
                                                    record.fMispEndRnC = fmisp.dGpTcEndRn; record.fMispEndRnT = fmisp.dGpTEndRn;
                                                    record.fMispAlignC = fmisp.dGpTcAlign; record.fMispAlignT = fmisp.dGpTAlign;
                                                    record.fScanC = fmisp.ScannedC; record.fScanT = fmisp.ScannedT;
                                                    record.fSurvC = fmisp.nSurvC; record.fSurvT = fmisp.nSurvT;
                                                    record.fAnchC = fmisp.nAnchC; record.fAnchT = fmisp.nAnchT;
                                                    record.dGfTcAnchEns    = fmisp.dGpTcAnchEns;
                                                    record.dGfTcAnchEnsRun = fmisp.dGpTcAnchEnsRun;
                                                    record.fAnchRunC       = fmisp.nAnchRunC;
                                                    record.dGfTAnchEns     = fmisp.dGpTAnchEns;
                                                    record.dGfTAnchEnsRun  = fmisp.dGpTAnchEnsRun;
                                                    record.fAnchRunT       = fmisp.nAnchRunT;
                                                    % The margin is measured on the EXTENSION-WEIGHTED
                                                    % free energy, so the quantity thresholded is
                                                    % theta*epsilon and not theta.
                                                    record.fTcAnchMargin = fmisp.dGpTcAnchEff - dGfTcmin;
                                                    record.fTAnchMargin  = fmisp.dGpTAnchEff  - dGfTmin;
                                                    % ENSEMBLE SHIFT: what every OTHER register
                                                    % adds on top of the worst one. Zero when the
                                                    % burden is a single site, negative when it is
                                                    % many. Ranked, never gated.
                                                    % BOTH TERMS ARE RUN-QUALIFIED, AND THEY HAVE
                                                    % TO BE. dGpTcAnchEff is the GATE WINNER:
                                                    % SeedScan updates it only inside if gated, so
                                                    % it describes a register that meets MispMinRun.
                                                    % Subtracting it from the ALL-anchored ensemble
                                                    % would difference two different populations.
                                                    % At MispSurveyRun = MispMinRun the two
                                                    % ensembles are identical by construction and
                                                    % nothing distinguishes them, which is why the
                                                    % error is invisible at the production setting;
                                                    % below that floor the all-anchored sum carries
                                                    % registers the gate says cannot prime, and
                                                    % charging multiplicity for those contradicts
                                                    % the gate. Reading the Run column makes this
                                                    % term floor-invariant in the same way the
                                                    % margin and the rejection already are.
                                                    record.fTcEnsShift = fmisp.dGpTcAnchEnsRun - fmisp.dGpTcAnchEff;
                                                    record.fTEnsShift  = fmisp.dGpTAnchEnsRun  - fmisp.dGpTAnchEff;
                                                    if check == true
                                                        record.RejectStage = "Mispriming_Margin";
                                                        % --- 15. discrimination margin ---
                                                        AssertMispSatisfiable(MispMargin, ...
                                                            [dGfTcmin, dGfTmin], phi, MispTemp(TAnneal, tmf, dTsub), 'forward');
                                                        fm = [record.fTcAnchMargin, record.fTAnchMargin];
                                                        fm = fm(isfinite(fm));
                                                        if ~isempty(fm)&&(min(fm) >= MispMargin)
                                                            record.RejectStage = "";
                                                            record.Status = "Selected";
                                                            ForwardCandidate = [ForwardCandidate; struct2table(record, "AsArray", true)];
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        ForwardLog = [ForwardLog; struct2table(record, "AsArray", true)];
        ForwardSearchCounter = ForwardSearchCounter + 1;
        waitbar((ForwardSearchCounter/TotalForwardSearchCount)/3, h, sprintf('Forward scan: %d of %d. Progress: %d%%', ...
            ForwardSearchCounter, TotalForwardSearchCount, round(100*(ForwardSearchCounter/TotalForwardSearchCount)/3)));
        if Verbose&&(mod(ForwardSearchCounter, ProgressEvery) == 0)
            fprintf('      %5d of %5d  (%3.0f%%)   surviving: %d\n', ForwardSearchCounter, ...
                TotalForwardSearchCount, 100*ForwardSearchCounter/TotalForwardSearchCount, ...
                height(ForwardCandidate));
        end
    end
end
results(1).Forward_Log = ForwardLog;
if Verbose
    fprintf('      done in %.1f s: %d of %d selected\n', toc(ScanClock), ...
        height(ForwardCandidate), TotalForwardSearchCount);
    ReportRejects(ForwardLog);
end

%% ================================================================== REVERSE PRIMER SEARCH =======
if ~isempty(ForwardCandidate)
    if Verbose, fprintf('[7/8] Reverse scan   %d candidates\n', TotalReverseSearchCount); end
    ScanClock = tic;
    ReverseLog = table();
    ReverseCandidate = table();
    for k = input.ScaffoldEnd: (input.ScaffoldEnd + rFlex)
        shiftplasmid = circshift(PlasmidSeq, -k);
        shiftrcplasmid = circshift(RevComPlas, -(PlasLen - k)); % reverse primer site at index 1 here
        for l = rLmin: rLmax
            priml = seqrcomplement(shiftplasmid(end - l + 1: end));
            rPos = mod(k-1, PlasLen) + 1;
            record = struct('rPos', rPos, ...
                'rlen', l, 'rseq', string(priml), 'rGC3', "", 'rGC5', "", ...
                'rGC', NaN, 'drlocgc35', NaN, 'rnumocc', NaN, 'rdirep', NaN, ...
                'rKmer3', "", 'rKmerN', NaN, ...
                'rRun3', NaN, 'rRunInt', NaN, 'rRunIntOff', NaN, 'rRunIntSeq', "", ...
                'dGr3end', NaN, ...
                'dGr3stack', NaN, 'dGr5stack', NaN, 'rProf', NaN, ...
                'prOpen', NaN, 'dGrOpen', NaN, 'rTm', NaN, ...
                'dGr', NaN, 'dGrrmin', NaN, 'rAlignr', "" , 'dGrr', NaN, 'dGrrend', NaN, 'rrTm', NaN, ...
                'rstemalign', "", 'dGrhp', NaN, 'rhpTm', NaN, 'rlockTm', NaN, 'rfold', NaN, ...
                'rBase3', "", 'rEps3', NaN, ...
                'dGrTcmin', NaN, 'dGrTc', NaN, 'dGrTmin', NaN, 'dGrT', NaN, ...
                'rTcMargin', NaN, 'rTMargin', NaN, ...
                'dGrTcAnch', NaN, 'dGrTAnch', NaN, 'dGrTcEnd', NaN, 'dGrTEnd', NaN, ...
                'dGrTcAnchEff', NaN, 'dGrTAnchEff', NaN, 'rTcEps', NaN, 'rTEps', NaN, ...
                'rTcEpsPair', "", 'rTEpsPair', "", ...
                'rTcAnchMargin', NaN, 'rTAnchMargin', NaN, ...
                'rTcEnsShift', NaN, 'rTEnsShift', NaN, ...
                'rMispPosC', NaN, 'rMispPosT', NaN, 'rMispEndRnC', NaN, 'rMispEndRnT', NaN, ...
                'rMispAlignC', "", 'rMispAlignT', "", ...
                'rSurvC', NaN, 'rSurvT', NaN, 'rAnchC', NaN, 'rAnchT', NaN, ...
                ... % Which strands the mispriming scan actually reached; see the forward side.
                'rScanC', false, 'rScanT', false, ...
                ... % Aggregate off-target flux, MEASUREMENT ONLY, nothing gates on these.
                ... % Ens is the ensemble free energy over every anchored register and EnsRun
                ... % the same sum restricted to registers carrying MispMinRun identities.
                'dGrTcAnchEns', NaN, 'dGrTcAnchEnsRun', NaN, 'rAnchRunC', NaN, ...
                'dGrTAnchEns',  NaN, 'dGrTAnchEnsRun',  NaN, 'rAnchRunT', NaN, ...
                'RejectStage', "Terminal_Bases", 'Status', "Eliminated");
            record.rGC5 = string(priml(1)); record.rGC3 = string(priml(end));
            record.rBase3 = string(priml(end));
            record.rEps3  = TerminalBaseDiscrimination(priml(end));
            % --- 1. terminal bases ---
            rGC5ok = ~Force5primeGC||((priml(1) == 'C')||(priml(1) == 'G'));
            rGC3ok = ~ReverseForce3primeGC||((priml(end) == 'C')||(priml(end) == 'G'));
            if rGC5ok&&rGC3ok
                record.RejectStage = "GC_Content";
                % --- 2. global GC ---
                GC = count(priml, pat); gcr = GC/strlength(priml);
                record.rGC = gcr;
                if ~UseGCWindow||((gcr > fGCmin)&&(gcr <= fGCmax))
                    record.RejectStage = "Local_GC_Skew";
                    % --- 3. half-primer GC asymmetry, logged always, enforced only if asked ---
                    rlocgc3 = count(priml((l-floor(l/2)+1):l),pat)/floor(l/2);
                    rlocgc5 = count(priml(1:floor(l/2)),pat)/floor(l/2);
                    record.drlocgc35 = rlocgc3 - rlocgc5;
                    if ~ReverseGCSkew||(rlocgc3 < rlocgc5)
                        record.RejectStage = "Low_Complexity";
                        % --- 4. low complexity ---
                        [rnumocc, rdirep, check] = NucleoRun(priml, MonoLim, DiRepLim);
                        record.rnumocc = rnumocc; record.rdirep = rdirep;
                        if check == true
                            record.RejectStage = "Three_Prime_Kmer";
                            % --- 5. 3'-terminal k-mer uniqueness ---
                            rKmer = priml(end - EndKmer + 1:end);
                            rKmerN = KmerCount(rKmer, PlasDbl, RevDbl, PlasLen);
                            record.rKmer3 = string(rKmer); record.rKmerN = rKmerN;
                            % --- 6 and 7. contiguous non-cognate runs, both O(1) from the maps ---
                            % The reverse primer is the reverse complement of plus(rLo : rPos), so
                            % its 3' terminus sits at the LOW end of that window and the maps are
                            % read from the opposite end to the forward case.
                            rLo = mod(rPos - l, PlasLen) + 1;
                            [record.rRun3, record.rRunInt, record.rRunIntOff, record.rRunIntSeq] = ...
                                RunLengths(RunMap, NonUniq, priml, rLo, l, PlasLen, false);
                            if rKmerN > EndKmerMax
                                rRunStage = "Three_Prime_Kmer";
                            elseif Use3PrimeRunGate&&(record.rRun3 >= MaxRun3)
                                rRunStage = "Three_Prime_Run";
                            elseif UseInternalRunGate&&(record.rRunInt >= MaxRunInt)
                                rRunStage = "Internal_Run";
                            else
                                rRunStage = "";
                            end
                            record.RejectStage = rRunStage;
                            if strlength(rRunStage) == 0
                                record.RejectStage = "Three_Prime_Stability";
                                % --- 8. 3'-terminal stability, published scale ---
                                % Logged and scored always; gated only when UseEndStabilityGate is on.
                                [dGr3end, dGr3stack, dGr5stack, rProf] = EndStability(priml, EndWin);
                                record.dGr3end = dGr3end;
                                record.dGr3stack = dGr3stack;
                                record.dGr5stack = dGr5stack;
                                record.rProf = rProf;
                                if (~UseEndStabilityGate)||(dGr3end >= dG3min)
                                    record.RejectStage = "Site_Accessibility";
                                    % --- 9. accessibility of the PLUS strand site ---
                                    % The reverse primer binds the plus strand over [rPos-l+1, rPos],
                                    % so the footprint's high end is rPos.
                                    [prOpen, dGrOpen] = OpeningCost(AccOn, LunpPlus, MapPlus, rPos, l, Rgas, AccTK);
                                    record.prOpen = prOpen; record.dGrOpen = dGrOpen;
                                    if dGrOpen <= AccMaxOpen
                                        record.RejectStage = "Melting_Temperature";
                                        % --- 10. melting temperature ---
                                        [tmr,~,~] = PrimTemp(priml, input.ReversePrimerConcentration, input.PlasmidConcentration, ...
                                            input.MonovalentCationConcentration, input.BivalentCationConcentration, input.dNTPsConcentration);
                                        record.rTm = tmr;
                                        if isreal(tmr)&&isfinite(tmr)&&(tmr >= Tmin)&&(tmr <= Tmax)
                                            record.RejectStage = "Homodimer";
                                            % --- 11. homodimer ---
                                            [dGr, dGrrmin, rAlignr, dGrr, check, dGrrend] = HomodimerizerMismatch(priml, phiDimer, input.ReversePrimerConcentration, ...
                                                input.PlasmidConcentration, input.MonovalentCationConcentration, ...
                                                input.BivalentCationConcentration, input.dNTPsConcentration, tmr, dTsub, TAnneal, DimerTemplateBeyond);
                                            record.dGr = dGr; record.dGrrmin = dGrrmin; record.rAlignr = AlignRow(rAlignr); record.dGrr = dGrr;
                                            record.rrTm = DimerMeltTm(priml, priml, ...
                                                input.ReversePrimerConcentration, [], ...
                                                input.MonovalentCationConcentration, ...
                                                input.BivalentCationConcentration, ...
                                                input.dNTPsConcentration, MispTemp(TAnneal, tmr, dTsub));
                                            if check == true
                                                record.RejectStage = "Three_Prime_Self_Dimer";
                                                % --- 12. 3'-anchored self-dimer ---
                                                % Sixth output of HomodimerizerMismatch, from DimerScan.
                                                record.dGrrend = dGrrend;
                                                if ProofExo, rrGate = dGrr; else, rrGate = dGrrend; end
                                                if rrGate >= dGEndMin
                                                    record.RejectStage = "Hairpin";
                                                    % --- 13. hairpin, judged on melting temperature ---
                                                    [rstemalign, check, rdGhp, rhpTm, rlockTm, rfold] = ...
                                                        LooperMismatch(priml, input.MonovalentCationConcentration, ...
                                                        input.BivalentCationConcentration, input.dNTPsConcentration, ...
                                                        MispTemp(TAnneal, tmr, dTsub), dGhpmin, TAnneal, HpTmMargin);
                                                    record.rstemalign = AlignRow(rstemalign);
                                                    record.dGrhp = rdGhp; record.rhpTm = rhpTm;
                                                    record.rlockTm = rlockTm; record.rfold = rfold;
                                                    if check == true
                                                        record.RejectStage = "Mispriming";
                                                        % --- 14. whole-plasmid mispriming ---
                                                        [dGrTcmin, dGrTc, dGrTmin, dGrT, check, rmisp] = MisprimerMismatch(shiftrcplasmid, priml, phi, ...
                                                            input.ReversePrimerConcentration, input.PlasmidConcentration, ...
                                                            input.MonovalentCationConcentration, input.BivalentCationConcentration, ...
                                                            input.dNTPsConcentration, tmr, dTsub, MispMinRun, MispEndTol, MispEndWin, ...
                                                            UseExtensionEfficiency, Rychlik, TAnneal, MispSurveyRun);
                                                        record.dGrTcmin = dGrTcmin; record.dGrTc = dGrTc;
                                                        record.dGrTmin = dGrTmin; record.dGrT = dGrT;
                                                        record.rTcMargin = dGrTc - dGrTcmin;
                                                        record.rTMargin  = dGrT  - dGrTmin;
                                                        record.dGrTcAnch = rmisp.dGpTcAnch; record.dGrTAnch = rmisp.dGpTAnch;
                                                        record.dGrTcAnchEff = rmisp.dGpTcAnchEff;
                                                        record.dGrTAnchEff  = rmisp.dGpTAnchEff;
                                                        record.rTcEps = rmisp.dGpTcEps; record.rTEps = rmisp.dGpTEps;
                                                        record.rTcEpsPair = rmisp.dGpTcEpsPair;
                                                        record.rTEpsPair  = rmisp.dGpTEpsPair;
                                                        record.dGrTcEnd  = rmisp.dGpTcEnd;  record.dGrTEnd  = rmisp.dGpTEnd;
                                                        record.rMispPosC = rmisp.dGpTcPos;  record.rMispPosT = rmisp.dGpTPos;
                                                        record.rMispEndRnC = rmisp.dGpTcEndRn; record.rMispEndRnT = rmisp.dGpTEndRn;
                                                        record.rMispAlignC = rmisp.dGpTcAlign; record.rMispAlignT = rmisp.dGpTAlign;
                                                        record.rScanC = rmisp.ScannedC; record.rScanT = rmisp.ScannedT;
                                                        record.rSurvC = rmisp.nSurvC; record.rSurvT = rmisp.nSurvT;
                                                        record.rAnchC = rmisp.nAnchC; record.rAnchT = rmisp.nAnchT;
                                                        record.dGrTcAnchEns    = rmisp.dGpTcAnchEns;
                                                        record.dGrTcAnchEnsRun = rmisp.dGpTcAnchEnsRun;
                                                        record.rAnchRunC       = rmisp.nAnchRunC;
                                                        record.dGrTAnchEns     = rmisp.dGpTAnchEns;
                                                        record.dGrTAnchEnsRun  = rmisp.dGpTAnchEnsRun;
                                                        record.rAnchRunT       = rmisp.nAnchRunT;
                                                        record.rTcAnchMargin = rmisp.dGpTcAnchEff - dGrTcmin;
                                                        record.rTAnchMargin  = rmisp.dGpTAnchEff  - dGrTmin;
                                                        % Run-qualified on both sides, as forward.
                                                        record.rTcEnsShift = rmisp.dGpTcAnchEnsRun - rmisp.dGpTcAnchEff;
                                                        record.rTEnsShift  = rmisp.dGpTAnchEnsRun  - rmisp.dGpTAnchEff;
                                                        if check == true
                                                            record.RejectStage = "Mispriming_Margin";
                                                            % --- 15. discrimination margin ---
                                                            AssertMispSatisfiable(MispMargin, ...
                                                                [dGrTcmin, dGrTmin], phi, MispTemp(TAnneal, tmr, dTsub), 'reverse');
                                                            rm = [record.rTcAnchMargin, record.rTAnchMargin];
                                                            rm = rm(isfinite(rm));
                                                            if ~isempty(rm)&&(min(rm) >= MispMargin)
                                                                record.RejectStage = "";
                                                                record.Status = "Selected";
                                                                ReverseCandidate = [ReverseCandidate; struct2table(record, "AsArray", true)];
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            ReverseLog = [ReverseLog; struct2table(record, "AsArray", true)];
            ReverseSearchCounter = ReverseSearchCounter + 1;
            waitbar(1/3 + (ReverseSearchCounter/TotalReverseSearchCount)/3, h, sprintf('Reverse scan: %d of %d. Progress: %d%%', ...
                ReverseSearchCounter, TotalReverseSearchCount, round(100*(1/3 + (ReverseSearchCounter/TotalReverseSearchCount)/3))));
            if Verbose&&(mod(ReverseSearchCounter, ProgressEvery) == 0)
                fprintf('      %5d of %5d  (%3.0f%%)   surviving: %d\n', ReverseSearchCounter, ...
                    TotalReverseSearchCount, 100*ReverseSearchCounter/TotalReverseSearchCount, ...
                    height(ReverseCandidate));
            end
        end
    end
    results(1).Reverse_Log = ReverseLog;
    if Verbose
        fprintf('      done in %.1f s: %d of %d selected\n', toc(ScanClock), ...
            height(ReverseCandidate), TotalReverseSearchCount);
        ReportRejects(ReverseLog);
    end

    %% ============================================================== PRIMER-PAIR SEARCH ==========
    if ~isempty(ReverseCandidate)
        PairLog = table();
        % Collapse screen: one scan per distinct terminus position, kept for the whole pairing loop.
        % The strand concentration is the ceiling of the free strand: the plus strand cannot exceed the
        % excess primer and the minus strand cannot exceed the limiting primer.
        CollapseCachePlus = containers.Map('KeyType', 'double', 'ValueType', 'any');
        CollapseCacheMinus = containers.Map('KeyType', 'double', 'ValueType', 'any');
        CollapseOptsPlus = struct('StrandConc', input.ForwardPrimerConcentration, 'LBase', CollapseLBase, ...
            'LPers', CollapseLPers, 'HpTmMargin', HpTmMargin, 'FlagEps', CollapseFlagEps, ...
            'Bimolecular', CollapseBimolecular, 'ExoBlend', CollapseExoBlend, 'ExoEndTol', CollapseExoEndTol);
        CollapseOptsMinus = CollapseOptsPlus;
        CollapseOptsMinus.StrandConc = input.ReversePrimerConcentration;
        TotalPairSearchCount = max(height(ForwardCandidate)*height(ReverseCandidate), 1);
        if Verbose
            fprintf('[8/8] Pairing        %d x %d = %d combinations\n', ...
                height(ForwardCandidate), height(ReverseCandidate), TotalPairSearchCount);
        end
        ScanClock = tic;
        for m = 1:height(ForwardCandidate)
            ForwardData = ForwardCandidate(m,:);
            fprimer = char(ForwardData.fseq); tmf = ForwardData.fTm;   % dGf no longer read here
            for n = 1:height(ReverseCandidate)
                ReverseData = ReverseCandidate(n,:);
                rprimer = char(ReverseData.rseq); tmr = ReverseData.rTm;   % dGr likewise
                record = struct('fPos', ForwardData.fPos, 'flen', ForwardData.flen, 'fseq', ForwardData.fseq, ...
                    'fAT3', ForwardData.fAT3, 'fGC5', ForwardData.fGC5, ...
                    'fGC', ForwardData.fGC, 'dflocgc35', ForwardData.dflocgc35, ...
                    'fnumocc', ForwardData.fnumocc, 'fdirep', ForwardData.fdirep, ...
                    'fKmer3', ForwardData.fKmer3, 'fKmerN', ForwardData.fKmerN, ...
                    'fRun3', ForwardData.fRun3, 'fRunInt', ForwardData.fRunInt, ...
                    'fRunIntOff', ForwardData.fRunIntOff, 'fRunIntSeq', ForwardData.fRunIntSeq, ...
                    'dGf3end', ForwardData.dGf3end, ...
                    'dGf3stack', ForwardData.dGf3stack, 'dGf5stack', ForwardData.dGf5stack, ...
                    'fProf', ForwardData.fProf, ...
                    'pfOpen', ForwardData.pfOpen, 'dGfOpen', ForwardData.dGfOpen, ...
                    'fTm', ForwardData.fTm, ...
                    'dGf', ForwardData.dGf, 'dGffmin', ForwardData.dGffmin, 'fAlignf', ForwardData.fAlignf, ...
                    'dGff', ForwardData.dGff, 'dGffend', ForwardData.dGffend, 'ffTm', ForwardData.ffTm, ...
                    'fstemalign', ForwardData.fstemalign, ...
                    'dGfhp', ForwardData.dGfhp, 'fhpTm', ForwardData.fhpTm, ...
                    'flockTm', ForwardData.flockTm, 'ffold', ForwardData.ffold, ...
                    'fBase3', ForwardData.fBase3, 'fEps3', ForwardData.fEps3, ...
                    'dGfTcmin', ForwardData.dGfTcmin, 'dGfTc', ForwardData.dGfTc, ...
                    'dGfTmin', ForwardData.dGfTmin, 'dGfT', ForwardData.dGfT, ...
                    'fTcMargin', ForwardData.fTcMargin, 'fTMargin', ForwardData.fTMargin, ...
                    'dGfTcAnch', ForwardData.dGfTcAnch, 'dGfTAnch', ForwardData.dGfTAnch, ...
                    'dGfTcEnd', ForwardData.dGfTcEnd, 'dGfTEnd', ForwardData.dGfTEnd, ...
                    'fTcAnchMargin', ForwardData.fTcAnchMargin, 'fTAnchMargin', ForwardData.fTAnchMargin, ...
                    'fTcEnsShift', ForwardData.fTcEnsShift, 'fTEnsShift', ForwardData.fTEnsShift, ...
                    'fMispPosC', ForwardData.fMispPosC, 'fMispPosT', ForwardData.fMispPosT, ...
                    'fMispEndRnC', ForwardData.fMispEndRnC, 'fMispEndRnT', ForwardData.fMispEndRnT, ...
                    'fMispAlignC', ForwardData.fMispAlignC, 'fMispAlignT', ForwardData.fMispAlignT, ...
                    'fScanC', ForwardData.fScanC, 'fScanT', ForwardData.fScanT, ...
                    'fSurvC', ForwardData.fSurvC, 'fSurvT', ForwardData.fSurvT, ...
                    'fAnchC', ForwardData.fAnchC, 'fAnchT', ForwardData.fAnchT, ...
                    'rPos', ReverseData.rPos, 'rlen', ReverseData.rlen, 'rseq', ReverseData.rseq, ...
                    'rGC3', ReverseData.rGC3, 'rGC5', ReverseData.rGC5, ...
                    'rGC', ReverseData.rGC, 'drlocgc35', ReverseData.drlocgc35, ...
                    'rnumocc', ReverseData.rnumocc, 'rdirep', ReverseData.rdirep, ...
                    'rKmer3', ReverseData.rKmer3, 'rKmerN', ReverseData.rKmerN, ...
                    'rRun3', ReverseData.rRun3, 'rRunInt', ReverseData.rRunInt, ...
                    'rRunIntOff', ReverseData.rRunIntOff, 'rRunIntSeq', ReverseData.rRunIntSeq, ...
                    'dGr3end', ReverseData.dGr3end, ...
                    'dGr3stack', ReverseData.dGr3stack, 'dGr5stack', ReverseData.dGr5stack, ...
                    'rProf', ReverseData.rProf, ...
                    'prOpen', ReverseData.prOpen, 'dGrOpen', ReverseData.dGrOpen, ...
                    'rTm', ReverseData.rTm, 'dTmrf', NaN, ...
                    'dGr', ReverseData.dGr, 'dGrrmin', ReverseData.dGrrmin, 'rAlignr', ReverseData.rAlignr, ...
                    'dGrr', ReverseData.dGrr, 'dGrrend', ReverseData.dGrrend, 'rrTm', ReverseData.rrTm, ...
                    'rstemalign', ReverseData.rstemalign, ...
                    'dGrhp', ReverseData.dGrhp, 'rhpTm', ReverseData.rhpTm, ...
                    'rlockTm', ReverseData.rlockTm, 'rfold', ReverseData.rfold, ...
                    'rBase3', ReverseData.rBase3, 'rEps3', ReverseData.rEps3, ...
                    'dGrTcmin', ReverseData.dGrTcmin, 'dGrTc', ReverseData.dGrTc, ...
                    'dGrTmin', ReverseData.dGrTmin, 'dGrT', ReverseData.dGrT, ...
                    'rTcMargin', ReverseData.rTcMargin, 'rTMargin', ReverseData.rTMargin, ...
                    'dGrTcAnch', ReverseData.dGrTcAnch, 'dGrTAnch', ReverseData.dGrTAnch, ...
                    'dGrTcEnd', ReverseData.dGrTcEnd, 'dGrTEnd', ReverseData.dGrTEnd, ...
                    'rTcAnchMargin', ReverseData.rTcAnchMargin, 'rTAnchMargin', ReverseData.rTAnchMargin, ...
                    'rTcEnsShift', ReverseData.rTcEnsShift, 'rTEnsShift', ReverseData.rTEnsShift, ...
                    'EnsShiftWorst', NaN, ...
                    'rMispPosC', ReverseData.rMispPosC, 'rMispPosT', ReverseData.rMispPosT, ...
                    'rMispEndRnC', ReverseData.rMispEndRnC, 'rMispEndRnT', ReverseData.rMispEndRnT, ...
                    'rMispAlignC', ReverseData.rMispAlignC, 'rMispAlignT', ReverseData.rMispAlignT, ...
                    'rScanC', ReverseData.rScanC, 'rScanT', ReverseData.rScanT, ...
                    'rSurvC', ReverseData.rSurvC, 'rSurvT', ReverseData.rSurvT, ...
                    'rAnchC', ReverseData.rAnchC, 'rAnchT', ReverseData.rAnchT, ...
                    'dGfrmin', NaN, 'fAlignr', "", 'dGfr', NaN, 'dGfrend', NaN, 'frTm', NaN, ...
                    'ampLen', NaN, 'flankLen', NaN, ...
                    'TmAP', NaN, 'TmXP', NaN, 'dTmAXP', NaN, ...
                    'TmA', NaN, 'ampGC', NaN, 'dTmAX', NaN, ...
                    ... % Amplicon competition as a ranked cost, never a rejection.
                    ... % TmAPenalty is ln(best tested slope / slope at this gap);
                    ... % TmARelYield is the implied fraction of that best slope;
                    ... % TmAExtrap marks a gap beyond the widest Pierce tested.
                    'TmAPenalty', NaN, 'TmARelYield', NaN, 'TmAExtrap', false, ...
                    ... % Amplicon 3'-OH collapse, log only (stage P5). cp* is the plus strand, whose
                    ... % 3'-OH sits at rPos. cm* is the minus strand, whose 3'-OH sits at fPos.
                    'cpTier', "NA", 'cpDGfold', NaN, 'cpDGfoldEng', NaN, 'cpTm', NaN, 'cpTmEng', NaN, ...
                    'cpLoop', NaN, 'cpEps', NaN, 'cpPExt', NaN, 'cpPExtEng', NaN, 'cpJoverC', NaN, 'cpSites', NaN, ...
                    'cpTierAll', "NA", 'cpPExtTot', NaN, 'cpBiShare', NaN, 'cpPExtTaq', NaN, 'cpPExtExo', NaN, ...
                    'cpPExtBlend', NaN, 'cpTmBi', NaN, ...
                    'cmTier', "NA", 'cmDGfold', NaN, 'cmDGfoldEng', NaN, 'cmTm', NaN, 'cmTmEng', NaN, ...
                    'cmLoop', NaN, 'cmEps', NaN, 'cmPExt', NaN, 'cmPExtEng', NaN, 'cmJoverC', NaN, 'cmSites', NaN, ...
                    'cmTierAll', "NA", 'cmPExtTot', NaN, 'cmBiShare', NaN, 'cmPExtTaq', NaN, 'cmPExtExo', NaN, ...
                    'cmPExtBlend', NaN, 'cmTmBi', NaN, ...
                    'Score', NaN, ...
                    'RejectStage', "Tm_Offset", 'Status', "Eliminated");
                % Amplicon geometry on the circular frame: forward 5' terminus to reverse 5' terminus.
                record.ampLen = mod(ReverseData.rPos - ForwardData.fPos, PlasLen) + 1;
                record.flankLen = record.ampLen - input.ScaffoldLength;
                % LATE-PCR: the limiting primer (reverse here, since ReversePrimerConcentration is the
                % smaller of the two) must melt at or above the excess primer, each Tm taken at its own
                % starting concentration. Sanchez et al. 2004 PNAS 101:1933; Pierce et al. 2005 PNAS
                % 102:8609.
                dTmrf = tmr - tmf;
                record.dTmrf = dTmrf;
                if (dTmrf >= dTrfmin)&&(dTmrf <= dTrfmax)
                    record.RejectStage = "Amplicon_Tm";
                    % --- P2. amplicon melting temperature ---
                    % Once the limiting primer is spent, the accumulating excess-primer strand
                    % competes with the excess primer itself for the same template strand, so the
                    % rate and extent of linear synthesis follow TmA - TmX. The amplicon runs from
                    % the forward 5' terminus for ampLen bases on the circle, read off the doubled
                    % plasmid so a wrap costs nothing.
                    ampSeq = PlasDbl(ForwardData.fPos : ForwardData.fPos + record.ampLen - 1);
                    [tmA, ampInfo] = TmAmplicon(ampSeq, input.MonovalentCationConcentration, ...
                        input.BivalentCationConcentration, input.dNTPsConcentration, UseTmAEquiv);
                    record.TmA = tmA; record.ampGC = ampInfo.GC;
                    record.dTmAX = tmA - tmf;
                    % THE GATE IS EVALUATED IN PIERCE'S CONVENTION, NOT IN THE ENGINE'S.
                    % TmAmaxGap and TmAoptGap are Pierce's numbers, and a number means nothing apart
                    % from the convention it was measured in. Pierce computed BOTH sides with
                    % monovalent cation only: the amplicon from the Wetmur %GC form with [M] the
                    % total monovalent, the primers from Allawi-SantaLucia nearest neighbours with
                    % the SantaLucia 1996 salt term, also monovalent. They say so in as many words,
                    % noting that the formulas underestimate Tm in magnesium-containing buffers and
                    % using them anyway for comparison.
                    % The engine does not work in that convention. PrimTemp applies the full
                    % Owczarzy divalent correction with free magnesium from the dNTP chelation
                    % equilibrium, so tmf carries magnesium. Subtracting a monovalent-only TmA from
                    % a magnesium-corrected TmX gives a difference that is neither convention.
                    % MEASURED on a 40 mM monovalent, 2.50 mM Mg, 1.20 mM dNTP buffer with a
                    % 3971 nt amplicon at 41.65% GC:
                    %     Pierce convention, both sides monovalent   75.05 - 55.59 = 19.46
                    %     UseTmAEquiv = false, mixed                 75.05 - 63.99 = 11.06
                    %     UseTmAEquiv = true,  mixed                 85.12 - 63.99 = 21.13
                    % The false setting understates the gap by 8.40 degrees and would wave through a
                    % target that is actually at the edge of Pierce's gate. The true setting is much
                    % closer, because the von Ahsen equivalent lifts TmA by 10.07 while Owczarzy
                    % lifts TmX by 8.47 and the two nearly cancel, but the residual 1.67 is an
                    % accident of this buffer and not a correction.
                    % So the gate is computed from a second pair of numbers, both monovalent, and
                    % Pierce's constants apply to them verbatim. No new constant is introduced and
                    % no signature changes: passing zero divalent and zero dNTP drives PrimTemp into
                    % its monovalent branch, and useEquiv = false does the same for TmAmplicon.
                    % The engine's monovalent branch is Owczarzy 2004 rather than SantaLucia 1996,
                    % which runs 0.76 degrees cooler across these twelve primers. That is inside the
                    % resolution of a gate whose source measured in whole degrees.
                    % record.TmA keeps the magnesium-corrected value, because that is the amplicon's
                    % real melting temperature and the one to quote in a methods section.
                    tmAP = TmAmplicon(ampSeq, input.MonovalentCationConcentration, 0, 0, false);
                    tmfP = PrimTemp(fprimer, input.ForwardPrimerConcentration, ...
                        input.PlasmidConcentration, input.MonovalentCationConcentration, 0, 0);
                    record.TmAP = tmAP; record.TmXP = tmfP;
                    record.dTmAXP = tmAP - tmfP;
                    % AMPLICON COMPETITION NO LONGER EXCLUDES A PAIR.
                    % What survives here is a validity guard: a pair whose monovalent gap cannot
                    % be computed has no scorable amplicon competition and is still flagged. Every
                    % pair with a real gap goes through, however wide, and the cost of a wide gap
                    % is carried by the yield penalty in the Score instead.
                    %
                    % Pierce's thresholds were measured on their own amplicons and their
                    % transferability to a multi-kilobase product is unestablished. A threshold
                    % that cannot be checked should not delete candidates; it should order them.
                    if isfinite(record.dTmAXP)
                    record.RejectStage = "Heterodimer";
                    [dGfrmin, fAlignr, dGfr, check, dGfrendF] = HeterodimerizerMismatch(fprimer, rprimer, ...
                        input.ForwardPrimerConcentration, input.ReversePrimerConcentration, phiDimer, ...
                        input.PlasmidConcentration, input.MonovalentCationConcentration, ...
                        input.BivalentCationConcentration, input.dNTPsConcentration, tmf, dTsub, TAnneal, DimerTemplateBeyond);
                    record.dGfrmin = dGfrmin; record.fAlignr = AlignRow(fAlignr); record.dGfr = dGfr;
                    record.frTm = DimerMeltTm(fprimer, rprimer, ...
                        input.ForwardPrimerConcentration, input.ReversePrimerConcentration, ...
                        input.MonovalentCationConcentration, ...
                        input.BivalentCationConcentration, ...
                        input.dNTPsConcentration, MispTemp(TAnneal, tmf, dTsub));
                    if check == true
                        record.RejectStage = "Three_Prime_Cross_Dimer";
                        % Both orientations, since either 3' terminus is extendable. dGfrendF is the
                        % fifth output of HeterodimerizerMismatch and anchors the FORWARD terminus;
                        % the reverse-anchored value comes from DimerScan directly rather than from a
                        % second HeterodimerizerMismatch call, because swapping the primers there
                        % would also swap their concentration roles inside dGHetero and recompute a
                        % threshold that means nothing.
                        %
                        % THE TEMPERATURE MUST BE THE ONE HeterodimerizerMismatch RESOLVED, and it
                        % resolves it with exactly this helper: the programmed annealing temperature
                        % when one is given, and tmf - dTsub only as a fallback. Writing tmf - dTsub
                        % here instead scores the two orientations of the SAME cross-dimer at
                        % different temperatures whenever an annealing temperature is supplied,
                        % which is the usual case. The reverse-anchored orientation then comes out
                        % several degrees hot, understates its own stability and slips past
                        % dGEndMin, leaving a one-sided blind spot that nothing downstream reports.
                        [~, dGfrendR] = DimerScan(rprimer, fprimer, input.MonovalentCationConcentration, ...
                            input.BivalentCationConcentration, input.dNTPsConcentration, ...
                            MispTemp(TAnneal, tmf, dTsub), [], DimerTemplateBeyond);
                        dGfrend = min(dGfrendF, dGfrendR);
                        record.dGfrend = dGfrend;
                        if ProofExo, frGate = dGfr; else, frGate = dGfrend; end
                        if frGate >= dGEndMin
                            record.RejectStage = "";
                            record.Status = "Selected";
                            % ---- ranking ----
                            % A positive dimer dG means no duplex forms at all, so clamp at 0 before
                            % measuring headroom. Without the clamp the margin keeps growing with how
                            % bad the least-bad non-dimer register happens to be, and the ranking
                            % rewards one non-structure over another. On a 3940 nt target that is
                            % exactly what separates four 21-mer forwards (dGff = +7.3, margin 9.8)
                            % from three 18/19-mers (dGff = +1.5, margin 4.0), when neither group
                            % has a forward homodimer at all.
                            % Rejection is unaffected: HomodimerizerMismatch and
                            % HeterodimerizerMismatch already guard with (dG < 0) && (dG < threshold).
                            ffMargin = min(record.dGff, 0) - record.dGffmin;
                            rrMargin = min(record.dGrr, 0) - record.dGrrmin;
                            frMargin = min(dGfr, 0)       - dGfrmin;
                            margin = [ffMargin, rrMargin, frMargin];
                            margin = margin(isfinite(margin));
                            if isempty(margin)
                                worstmargin = 0;
                            else
                                worstmargin = min(min(margin), MarginCap);
                            end
                            % Worst mispriming headroom across both primers and both strands.
                            % Anchored registers only, since an unextendable flapped duplex has no
                            % business moving a candidate up or down the ranking.
                            mispall = [record.fTcAnchMargin, record.fTAnchMargin, ...
                                record.rTcAnchMargin, record.rTAnchMargin];
                            mispall = mispall(isfinite(mispall));
                            if isempty(mispall)
                                worstmisp = 0;
                            else
                                worstmisp = min(min(mispall), MispCap);
                            end
                            % OFF-TARGET MULTIPLICITY. The worst ensemble shift across both
                            % primers and both strands, as a POSITIVE penalty: the shift is at
                            % most zero, so its negation is the burden. A candidate whose
                            % off-target weight sits on one register scores zero here and is not
                            % moved at all, which is the behaviour this term is for.
                            ensall = [record.fTcEnsShift, record.fTEnsShift, ...
                                record.rTcEnsShift, record.rTEnsShift];
                            ensall = ensall(isfinite(ensall));
                            if isempty(ensall)
                                ensworst = 0;
                            else
                                ensworst = min(max(-min(ensall), 0), EnsCap);
                            end
                            record.EnsShiftWorst = ensworst;
                            % Total opening cost the pair has to pay against template structure.
                            opencost = 0;
                            if isfinite(record.dGfOpen), opencost = opencost + record.dGfOpen; end
                            if isfinite(record.dGrOpen), opencost = opencost + record.dGrOpen; end
                            % AMPLICON OVERHANG, RANKED WEAKLY.
                            % flank is the bases the amplicon carries beyond the scaffold, which is
                            % non-negative because both search windows only ever grow outward. It
                            % is normalised by the whole admissible overhang, so the term is
                            % bounded by WtFlank and that number reads directly against the other
                            % ceilings. The extra nucleotides are unpaired tails rather than a
                            % thermodynamic liability, so this states a preference for a tighter
                            % amplicon without letting geometry overrule the duplex terms. A hard
                            % limit on the overhang belongs in a gate, not in this weight.
                            flexSpan = max(fFlex + rFlex, 1);
                            flankTerm = max(record.flankLen, 0)/flexSpan;
                            % 3'-END STABILITY AND INTERNAL PROFILE, now ranked rather than gated.
                            % Both terms are PENALTIES for a grippy 3' end: dG3end is negative and
                            % more negative is worse, so its magnitude enters with a plus sign;
                            % prof is positive when the 3' end is the looser one, so it enters with
                            % a minus. Worst case over the pair, since either primer can misprime.
                            % Capped at the most stable pentamer that exists on each scale, so one
                            % extreme primer cannot swamp the dimer and mispriming terms.
                            end3 = max([-record.dGf3end, -record.dGr3end, 0]);
                            end3 = min(end3, 6.86);
                            profworst = min([record.fProf, record.rProf, 0]);
                            profworst = max(profworst, -6.0);
                            % AMPLICON COMPETITION, PRICED IN YIELD RATHER THAN IN DEGREES.
                            % The penalty is the natural log of how much early linear-phase slope
                            % this gap gives up against Pierce's best tested gap, interpolated
                            % through their four measured points. Zero at and below TmAoptGap,
                            % because their curve is still climbing there and nothing is known
                            % about narrower gaps, so none is rewarded. Monotone above it. The
                            % quantity is dimensionless and WtTmA reads as Score points per e-fold
                            % of lost slope.
                            [tmaTerm, tmaYield, tmaExtrap] = TmAYieldPenalty(record.dTmAXP, ...
                                TmAPierceGap, TmAPierceSlope, TmAoptGap);
                            % THE PENALTY IS CAPPED BEFORE IT IS WEIGHTED.
                            % Past Pierce's widest tested gap the penalty extends a trend instead
                            % of reporting a measurement, and it grows linearly in degrees with
                            % nothing to stop it. TmACap bounds the term at WtTmA*TmACap, which is
                            % set equal to the dimer ceiling, while leaving the ordering intact
                            % across every gap the search actually produces. record.TmAPenalty
                            % keeps the UNCAPPED value, because that is the physical quantity and
                            % the one to read when deciding whether the cap is placed correctly.
                            record.TmAPenalty = tmaTerm;
                            tmaTerm = min(tmaTerm, TmACap);
                            record.TmARelYield = tmaYield;
                            record.TmAExtrap = tmaExtrap;
                            % Primer availability. Worst folded fraction of the pair at the
                            % annealing step, 0 to 1, so the weight is the whole cost of losing
                            % every primer molecule to its own hairpin.
                            foldworst = max([record.ffold, record.rfold, 0]);
                            % 3'-terminal base. Larger is better discrimination, hence the minus.
                            base3worst = min([record.fEps3, record.rEps3]);
                            record.Score = abs(dTmrf - 0.5*(dTrfmin + dTrfmax)) ...
                                + WtTmA*tmaTerm ...
                                - WtMargin*worstmargin ...
                                - WtMisp*worstmisp ...
                                + WtEns*ensworst ...
                                + WtAccess*opencost ...
                                + Wt3End*end3 ...
                                - WtProf*profworst ...
                                + WtFlank*flankTerm ...
                                + WtHairpin*foldworst ...
                                - Wt3Base*base3worst;
                            % Stage P5, log only. Nothing below reads these fields except the output.
                            if CollapseLog
                                cres = CollapseForPair(ForwardData.fPos, ReverseData.rPos, PlasDbl, PlasLen, ...
                                    input, fFlex, rFlex, CollapseWindow, CollapseMinRun, TAnneal, ...
                                    UseExtensionEfficiency, CollapseOptsPlus, CollapseOptsMinus, CollapseCachePlus, CollapseCacheMinus);
                                record = CollapseToRecord(record, 'cp', cres.plus);
                                record = CollapseToRecord(record, 'cm', cres.minus);
                            end
                            tabrow = table(ForwardData.fPos, ForwardData.fseq, tmf, ForwardData.flen, ForwardData.fGC, ...
                                ReverseData.rPos, ReverseData.rseq, tmr, ReverseData.rlen, ReverseData.rGC, ...
                                dTmrf, record.ampLen, record.flankLen, record.TmA, record.dTmAX, ...
                                record.TmAP, record.TmXP, record.dTmAXP, ...
                                ffMargin, rrMargin, frMargin, ...
                                ForwardData.fKmer3, ForwardData.fKmerN, ReverseData.rKmer3, ReverseData.rKmerN, ...
                                ForwardData.fRun3, ForwardData.fRunInt, ...
                                ReverseData.rRun3, ReverseData.rRunInt, ...
                                record.fTcAnchMargin, record.fTAnchMargin, ...
                                record.rTcAnchMargin, record.rTAnchMargin, ...
                                record.fTcMargin, record.fTMargin, record.rTcMargin, record.rTMargin, ...
                                record.dGfTcEnd, record.dGfTEnd, record.dGrTcEnd, record.dGrTEnd, ...
                                record.dGfOpen, record.dGrOpen, record.Score, ...
                                record.cpTier, record.cmTier, record.cpPExt, record.cmPExt, ...
                                record.cpTm, record.cmTm, ...
                                record.cpTierAll, record.cmTierAll, record.cpPExtTot, record.cmPExtTot, ...
                                'VariableNames', {'F_Start', 'F_Primer', 'F_Tm', 'F_Length', 'F_GC', ...
                                'R_Start', 'R_Primer', 'R_Tm', 'R_Length', 'R_GC', ...
                                'dTm', 'Amplicon_Length', 'Flank_Length', 'Amplicon_Tm', 'dTm_A_X', ...
                                'Amplicon_Tm_Pierce', 'Excess_Tm_Pierce', 'dTm_A_X_Pierce', ...
                                'dG_ff_Margin', 'dG_rr_Margin', 'dG_fr_Margin', ...
                                'F_Kmer3', 'F_Kmer_N', 'R_Kmer3', 'R_Kmer_N', ...
                                'F_Run3', 'F_RunInt', 'R_Run3', 'R_RunInt', ...
                                'F_Misp_Cog_Anch_Margin', 'F_Misp_Oth_Anch_Margin', ...
                                'R_Misp_Cog_Anch_Margin', 'R_Misp_Oth_Anch_Margin', ...
                                'F_Misp_Cog_Glob_Margin', 'F_Misp_Oth_Glob_Margin', ...
                                'R_Misp_Cog_Glob_Margin', 'R_Misp_Oth_Glob_Margin', ...
                                'F_Misp_Cog_dG3end', 'F_Misp_Oth_dG3end', ...
                                'R_Misp_Cog_dG3end', 'R_Misp_Oth_dG3end', ...
                                'F_dGopen', 'R_dGopen', 'Score', ...
                                'Collapse_Plus_Tier', 'Collapse_Minus_Tier', ...
                                'Collapse_Plus_pExt', 'Collapse_Minus_pExt', ...
                                'Collapse_Plus_Tm', 'Collapse_Minus_Tm', ...
                                'Collapse_Plus_TierAll', 'Collapse_Minus_TierAll', ...
                                'Collapse_Plus_pExtTot', 'Collapse_Minus_pExtTot'});
                            tab = [tab; tabrow];
                        end
                    end
                    end   % closes the amplicon gap validity guard
                end
                % Appended exactly once. The original wrote every accepted pair into PairLog twice.
                PairLog = [PairLog; struct2table(record, "AsArray", true)];
                PairSearchCounter = PairSearchCounter + 1;
                waitbar(2/3 + (PairSearchCounter/TotalPairSearchCount)/3, h, sprintf('Pairing: %d of %d. Progress: %d%%', ...
                    PairSearchCounter, TotalPairSearchCount, round(100*(2/3 + (PairSearchCounter/TotalPairSearchCount)/3))));
                if Verbose&&(mod(PairSearchCounter, ProgressEvery) == 0)
                    fprintf('      %5d of %5d  (%3.0f%%)   accepted: %d\n', PairSearchCounter, ...
                        TotalPairSearchCount, 100*PairSearchCounter/TotalPairSearchCount, nnz(strcmp(string(PairLog.Status), "Selected")));
                end
            end
        end
        if Verbose
            fprintf('      done in %.1f s\n', toc(ScanClock));
            ReportRejects(PairLog);
            ReportDimerMelts(PairLog, TAnneal);
            ReportCollapse(PairLog, CollapseLog);
            fprintf('%s\n', repmat('-', 1, 96));
            fprintf('Search complete in %.1f s total.\n', toc(RunClock));
        end
        results(1).Pair_Log = PairLog;
    else
        fprintf('No suitable reverse candidate found! Inspect the RejectStage column of results.Reverse_Log.\n');
        return
    end
else
    fprintf('No suitable forward candidate found! Inspect the RejectStage column of results.Forward_Log.\n');
    return
end
if isempty(tab)
    fprintf(['No suitable primer-pair found! Inspect the RejectStage column of results.Pair_Log. The ' ...
        'filters run in this order: terminal bases, GC content, local GC skew, low complexity, ' ...
        '3''-terminal k-mer, 3''-terminal contiguous run, internal contiguous run, ' ...
        '3''-terminal stability, site accessibility, Tm window, homodimer, 3''-anchored ' ...
        'self-dimer, hairpin, mispriming, mispriming margin.\n']);
    return
else
    % Ranked best first, so row 1 is the recommendation aPCR_Claude offers as its dialog default.
    tab = sortrows(tab, 'Score', 'ascend');
    results.Candidate_Primers = tab;
end
end

%% ____________________________________________________________________________________________________________
function key = AccCacheKey(seqPlus, aopts, paramFile)
%% ACCCACHEKEY  Everything that can change the fold, in one struct
% The fold is a pure function of these. Anything NOT here is either irrelevant to lunp (the dummy
% site list, which only drives the discarded siteTable) or is a hole, and the one hole is named.
key = struct();
key.CacheFormat = 1;                      % bump when the stored layout changes, which invalidates
key.SeqLen      = numel(seqPlus);
key.Seq         = upper(char(seqPlus));   % stored whole; the hash is for the filename, this is the test
key.W           = aopts.W;
key.L           = aopts.L;
key.uMax        = aopts.uMax;
key.Temp        = aopts.Temp;
key.Salt        = aopts.Salt;
key.Circular    = logical(aopts.Circular);
key.NucLen      = aopts.NucLen;
% The energy parameter file. When the caller leaves it empty DNAAccessibility resolves it privately
% and the engine cannot see which file was used, so the cache would not notice it being swapped.
% Recording 'auto' makes that visible rather than silent.
if isempty(paramFile)
    key.ParamFile = 'auto';
    key.ParamBytes = -1;
    key.ParamDate  = '';
else
    key.ParamFile = char(paramFile);
    d = dir(key.ParamFile);
    if isempty(d)
        key.ParamBytes = -1; key.ParamDate = '';
    else
        key.ParamBytes = d(1).bytes;
        key.ParamDate  = char(datetime(d(1).datenum, 'ConvertFrom', 'datenum', ...
            'Format', 'yyyy-MM-dd HH:mm:ss'));
    end
end
end

function p = AccCachePath(root, plasFileName, key)
%% ACCCACHEPATH  Where this exact parameter set lives on disk
% <root>\<Plasmid_Name>\<Plasmid_Name>_n<len>_W<..>_L<..>_u<..>_T<..>C_S<..>M_circ<0|1>_nuc<..>_<hash>.mat
% The parameters sit in the NAME so the folder can be read without opening anything, and the eight
% hex digits make two different sequences of the same length under the same parameters land on
% different files. The name is navigation. Correctness is the key stored inside the file.
if isempty(root)
    if exist('F:\', 'dir'), root = 'F:\aPCR_Computations\DNA_Accessibility';
    else,                   root = 'D:\aPCR_Computations\DNA_Accessibility';
    end
end
plasName = erase(char(plasFileName), '.txt');
if isempty(plasName), plasName = 'Unnamed'; end
plasName = regexprep(plasName, '[^\w\-]', '_');
stem = sprintf('%s_n%d_W%d_L%d_u%d_T%.1fC_S%.4fM_circ%d_nuc%d_%s', ...
    plasName, key.SeqLen, key.W, key.L, key.uMax, key.Temp, key.Salt, ...
    double(key.Circular), key.NucLen, AccKeyHash(key));
p = fullfile(root, plasName, [stem '.mat']);
end

function h = AccKeyHash(key)
%% ACCKEYHASH  Eight hex digits of FNV-1a 32 over the key, for the filename only
% Pure MATLAB on purpose: no Java, no toolbox, identical on every machine that runs the suite.
% THE ARITHMETIC IS DONE IN DOUBLES, NOT uint32 OR uint64, because MATLAB integer types SATURATE
% instead of wrapping and a hash needs wrapping. A double holds integers exactly to 2^53, so the
% multiply is split into 16-bit halves to keep every partial product under 2^40:
%       h*p mod 2^32  =  ( (hiH*p mod 2^16)*2^16 + loH*p )  mod 2^32
% with hiH, loH < 2^16 and p = 16777619 < 2^24, so each product is below 2^40 and exact.
% Verified against the canonical FNV-1a 32 vectors: '' -> 811c9dc5, 'a' -> e40c292c and
% 'foobar' -> bf9cf968, and cross-checked against a reference implementation.
% A collision would put two parameter sets on one filename; the loader still compares the full key
% stored inside the file and calls it a miss, so a collision costs a refold and never a wrong
% answer.
txt = sprintf('%d|%d|%d|%d|%d|%.6f|%.8f|%d|%d|%s|%d|%s|', ...
    key.CacheFormat, key.SeqLen, key.W, key.L, key.uMax, key.Temp, key.Salt, ...
    double(key.Circular), key.NucLen, key.ParamFile, key.ParamBytes, key.ParamDate);
h = sprintf('%08s', Fnv1a32([txt key.Seq]));
end

function hex = Fnv1a32(str)
%% FNV1A32  FNV-1a 32-bit over the bytes of a char row, returned as 8 lowercase hex digits
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

function ok = AccKeyMatch(a, b)
%% ACCKEYMATCH  Exact comparison of two keys, sequence included
ok = false;
f = {'CacheFormat', 'SeqLen', 'W', 'L', 'uMax', 'Circular', 'NucLen', 'ParamBytes'};
for i = 1:numel(f)
    if ~isfield(b, f{i}) || ~isequal(double(a.(f{i})), double(b.(f{i}))), return; end
end
% Temperature and salt are doubles written through a text filename, so compare on a tolerance that
% is far below anything that changes a fold and far above double round-trip noise.
if ~isfield(b, 'Temp')||~isfield(b, 'Salt'), return; end
if abs(a.Temp - b.Temp) > 1E-6 || abs(a.Salt - b.Salt) > 1E-9, return; end
for g = {'ParamFile', 'ParamDate'}
    if ~isfield(b, g{1}) || ~strcmp(char(a.(g{1})), char(b.(g{1}))), return; end
end
if ~isfield(b, 'Seq') || ~strcmp(a.Seq, b.Seq), return; end
ok = true;
end

function [hit, LunpPlus, MapPlus, LunpMinus, MapMinus] = AccCacheLoad(cacheFile, key, Verbose)
%% ACCCACHELOAD  Read a cached fold, or report a miss and say what IS cached
% Every failure path is a miss, never an error. A cache that cannot be read must cost a refold and
% nothing else.
hit = false; LunpPlus = []; MapPlus = []; LunpMinus = []; MapMinus = [];
if isempty(cacheFile), return; end
if ~isfile(cacheFile)
    if Verbose
        fprintf('      Cache MISS. Looking for %s\n', cacheFile);
        d = dir(fullfile(fileparts(cacheFile), '*.mat'));
        if isempty(d)
            fprintf('      Nothing cached for this plasmid yet.\n');
        else
            fprintf('      %d other parameter set(s) already cached here:\n', numel(d));
            for i = 1:min(numel(d), 8), fprintf('        %s\n', d(i).name); end
            if numel(d) > 8, fprintf('        ... and %d more\n', numel(d) - 8); end
        end
    end
    return
end
try
    C = load(cacheFile);
catch ME
    warning('LocalPrimerSearchEngine:AccCacheUnreadable', ...
        'Accessibility cache %s could not be read (%s); refolding.', cacheFile, ME.message);
    return
end
need = {'Key', 'LunpPlus', 'LunpMinus', 'MapN', 'MapOff'};
for i = 1:numel(need)
    if ~isfield(C, need{i})
        warning('LocalPrimerSearchEngine:AccCacheIncomplete', ...
            'Accessibility cache %s is missing %s; refolding.', cacheFile, need{i});
        return
    end
end
if ~AccKeyMatch(key, C.Key)
    % The filename hash matched but the stored key did not. Either a hash collision or a file
    % written by an older format. Both mean refold.
    warning('LocalPrimerSearchEngine:AccCacheKeyMismatch', ...
        ['Accessibility cache %s carries a different parameter set from the one requested, ' ...
        'despite the matching filename. Refolding and overwriting.'], cacheFile);
    return
end
% Rebuild the position map rather than deserialising a function handle. Its form is fixed by
% DNAAccessibility and the two scalars it closes over were recorded at save time, verified there
% against the handle the fold actually returned.
MapPlus  = @(p) mod(p - 1, C.MapN) + 1 + C.MapOff;
MapMinus = MapPlus;
LunpPlus = C.LunpPlus; LunpMinus = C.LunpMinus;
if size(LunpPlus, 2) < key.uMax + 1 || size(LunpMinus, 2) < key.uMax + 1
    warning('LocalPrimerSearchEngine:AccCacheTooNarrow', ...
        ['Cached lunp has %d columns but uMax = %d needs %d; refolding.'], ...
        size(LunpPlus, 2), key.uMax, key.uMax + 1);
    LunpPlus = []; LunpMinus = []; MapPlus = []; MapMinus = [];
    return
end
hit = true;
if Verbose
    fprintf('      Cache HIT.  %s\n', cacheFile);
    if strcmp(key.ParamFile, 'auto')
        fprintf(['      NOTE AccParamFile is empty, so this cache cannot tell whether ' ...
            'dna_mathews2004.par changed since it was written.\n']);
    end
end
end

function AccCacheStore(cacheFile, key, LunpPlus, MapPlus, LunpMinus, MapMinus, N, Verbose)
%% ACCCACHESTORE  Write the fold, after proving the rebuilt map reproduces the real one
% The map is stored as the two scalars it closes over, not as a serialised handle, because an
% anonymous handle is the one thing in this payload whose round trip depends on the MATLAB version.
% Before writing, the scalars are checked against BOTH handles the fold returned, at positions that
% exercise the wrap. A disagreement means the construction in DNAAccessibility changed and this
% function has gone out of sync with it, so nothing is written and the run carries on uncached.
if isempty(cacheFile), return; end
off = MapPlus(1) - 1;
probe = unique([1, 2, N - 1, N, N + 1, 2*N, round(N/2)]);
probe = probe(probe >= 1);
rebuilt = mod(probe - 1, N) + 1 + off;
if ~isequal(rebuilt, arrayfun(MapPlus, probe)) || ~isequal(rebuilt, arrayfun(MapMinus, probe))
    warning('LocalPrimerSearchEngine:AccCacheMapDrift', ...
        ['The position map returned by DNAAccessibility no longer matches ' ...
        'mod(p-1,N)+1+off, which is the form AccCacheStore knows how to record. Nothing was ' ...
        'cached. Update AccCacheStore and AccCacheLoad together before relying on the cache ' ...
        'again.']);
    return
end
d = fileparts(cacheFile);
if ~isfolder(d)
    [ok, msg] = mkdir(d);
    if ~ok
        warning('LocalPrimerSearchEngine:AccCacheNoDir', ...
            'Could not create %s (%s); continuing without caching.', d, msg);
        return
    end
end
Key = key; MapN = N; MapOff = off; %#ok<NASGU>
WrittenOn = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')); %#ok<NASGU>
try
    save(cacheFile, 'Key', 'LunpPlus', 'LunpMinus', 'MapN', 'MapOff', 'WrittenOn', '-v7.3');
    if Verbose
        f = dir(cacheFile);
        fprintf('      Fold cached to %s  (%.1f MB)\n', cacheFile, f(1).bytes/1048576);
    end
catch ME
    warning('LocalPrimerSearchEngine:AccCacheWriteFailed', ...
        'Could not write %s (%s); continuing without caching.', cacheFile, ME.message);
end
end

function AssertTmAGapSatisfiable(PlasDbl, input, Tmax, TmAmaxGap, dTrfmin, fLmax)
%% ASSERTTMAGAPSATISFIABLE  Refuse to start when P1 and P2 cannot both be met
% Two gates squeeze the EXCESS primer from opposite sides, and on a long amplicon they can close
% completely. The engine used to discover this after scanning the whole target and returning an
% empty table, which looks like a bug in the scan and is not.
%
%   P2 pushes the excess primer UP.   TmA - TmX <= TmAmaxGap  needs  TmX >= TmA - TmAmaxGap.
%                                     TmA is a property of the target. Past a kilobase the
%                                     -500/length term is dead, so TmA is set by GC and salt alone
%                                     and no primer choice moves it.
%   P1 pushes the excess primer DOWN. The limiting primer must melt dTrfmin above the excess one
%                                     and is itself capped by Tmax, so TmX <= Tmax - dTrfmin.
%
% Feasible only when  TmA - TmAmaxGap  <=  Tmax - dTrfmin, expressed in ONE set of units.
%
% THE TWO SIDES ARE NOT IN THE SAME UNITS, which is the part that makes this worth a function.
% TmA - TmAmaxGap lives in Pierce's monovalent convention, since that is where TmAmaxGap was
% measured. Tmax is applied to PrimTemp output, which carries magnesium through Owczarzy. The
% bridge is delta, the amount magnesium lifts a primer melting temperature in this buffer,
% measured here on a real candidate from the forward window rather than assumed. On a 40 mM
% monovalent buffer with 2.50 mM magnesium against 1.20 mM dNTP, delta is 8.47 degrees across
% twelve primers with a spread of 0.5, so one representative candidate is enough.
%
% MEASURED INSTANCE on that buffer. TmA 75.05 monovalent, TmAmaxGap 20, delta 8.47, so the excess
% primer has to reach 63.52. A Tmax of 65 less dTrfmin 4 caps it at 61.00. Empty by 2.52 degrees,
% for every sequence on the target, and without this check discoverable only after a full scan
% returns nothing.
% THE SCAFFOLD LENGTH IS READ ON THE CIRCLE, THE SAME WAY THE SEARCH READS AN AMPLICON.
% A scaffold may cross the origin, which CollectaPCRData_Claude.m supports and branches on. The
% linear difference goes NEGATIVE there, the slice below comes back empty, and TmAmplicon stops
% the run with a message about a zero-length amplicon that points at the target rather than at
% the arithmetic. mod(...)+1 is the form the pair loop uses, so the two agree by construction.
PlasLenLocal = numel(PlasDbl)/2;
ampLen = mod(input.ScaffoldEnd - input.ScaffoldStart, PlasLenLocal) + 1;
ampSeq = PlasDbl(input.ScaffoldStart : input.ScaffoldStart + ampLen - 1);
tmA    = TmAmplicon(ampSeq, input.MonovalentCationConcentration, 0, 0, false);

% delta, on a real candidate: the longest forward starting at the scaffold 5' end. Passing zero
% divalent and zero dNTP drives PrimTemp into its monovalent branch, so no signature changes.
probe = PlasDbl(input.ScaffoldStart : input.ScaffoldStart + fLmax - 1);
tmMg  = PrimTemp(probe, input.ForwardPrimerConcentration, input.PlasmidConcentration, ...
    input.MonovalentCationConcentration, input.BivalentCationConcentration, input.dNTPsConcentration);
tmMo  = PrimTemp(probe, input.ForwardPrimerConcentration, input.PlasmidConcentration, ...
    input.MonovalentCationConcentration, 0, 0);
delta = tmMg - tmMo;

floorX = tmA - TmAmaxGap + delta;     % lowest excess-primer Tm P2 will accept, engine units
ceilX  = Tmax - dTrfmin;              % highest excess-primer Tm P1 leaves room above
if ~isfinite(floorX)||~isfinite(ceilX)||(floorX <= ceilX), return; end
warning('LocalPrimerSearchEngine:TmAGapUnsatisfiable', ...
    ['ADVISORY. Pierce''s amplicon criterion and the LATE-PCR offset cannot both be met on this ' ...
    'target, by ANY primer. This no longer stops the run: amplicon competition is ranked rather ' ...
    'than gated, so every pair still appears, ordered by the linear-phase yield it gives up. ' ...
    'Read it as a statement about the target.\n' ...
    'Amplicon %d nt, TmA = %.2f C monovalent, which is the convention TmAmaxGap is measured in.\n' ...
    'Magnesium lifts a primer Tm here by %.2f C, so in the units Tmax is applied in:\n' ...
    '    P2 needs the excess primer at or above  %.2f C   (TmA - TmAmaxGap + delta)\n' ...
    '    P1 leaves room only up to               %.2f C   (Tmax - dTrfmin)\n' ...
    'The window is empty by %.2f C.\n' ...
    'TmA cannot be lowered by choosing different primers: past a kilobase it is set by GC and ' ...
    'salt alone. Three real ways out.\n' ...
    '  1. Raise Tmax to at least %.2f, then check that a limiting primer that stable exists ' ...
    'inside the length and GC windows. This is the only route that keeps both of Pierce''s ' ...
    'criteria intact.\n' ...
    '  2. Raise fLmax so the EXCESS primer can be longer and hotter. This attacks the gap ' ...
    'directly and is what Pierce recommend, that excess primers with TmX closer to TmA give ' ...
    'higher linear rates. It costs nothing in specificity while the extra length stays unique.\n' ...
    '  3. Raise TmAmaxGap to at least %.2f, which is agreeing to a slower linear phase. Pierce ' ...
    'measured early slopes of 11.3, 7.6, 4.2 and 2.3 at gaps of 13, 17, 20 and 23, so past 23 ' ...
    'you are extrapolating beyond anything they tested.\n' ...
    'A scaffold-length amplicon is intrinsically hard for LATE-PCR for exactly this reason. The ' ...
    'gates are reporting that, not malfunctioning.'], ...
    ampLen, tmA, delta, floorX, ceilX, floorX - ceilX, floorX + dTrfmin, ...
    tmA + delta - (Tmax - dTrfmin));
end

function AssertMispSatisfiable(MispMargin, thr, phi, tann, side)
%% ASSERTMISPSATISFIABLE  Stop a run whose mispriming margin no primer can meet
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
%
% MisprimerMismatch clamps its anchored accumulator at zero. A register whose free energy is
% positive does not form, so it is reported as "no stabilising anchored register" rather than as a
% positive number. Stage 15 then computes margin = dGanch - threshold, and since dGanch <= 0 the
% margin cannot exceed -threshold whatever the primer happens to be.
%
% So when MispMargin reaches -threshold the stage rejects every candidate by arithmetic alone,
% including a primer carrying no extendable off-target register anywhere on either strand. The run
% ends with an empty table and the message "relax the primer search criteria", which is the one
% thing that will not help.
%
% MEASURED INSTANCE. At phi = 1E-4 the thresholds come out between -2.994 and -3.027 kcal/mol, so
% the ceiling on the margin is 2.994 to 3.027, against a MispMargin of 3.0. A candidate with NO
% stabilising anchored register on the cognate strand scores the maximum possible margin of 2.996
% and is eliminated for being 0.004 kcal/mol short of a bound it cannot reach.
%
% The two constants are not independent. Demanding m kcal/mol of headroom over the threshold is
% algebraically identical to lowering phi by exp(-m/RT), because the threshold enters as
% -RT*log(phi*K/(1+K*f)). Expressing the specificity target twice, once as a ratio and once as an
% energy, hides the product from whoever sets either one.
%
% INPUTS
%   MispMargin  the constant under test, kcal/mol
%   thr         the thresholds in play for this candidate, [cognate-strand, other-strand]
%   phi         tolerated ratio of off-target duplex to cognate duplex at one site
%   tann        annealing temperature for this candidate, degrees C
%   side        'forward' or 'reverse', for the message
% ____________________________________________________________________________________________________________
thr = thr(isfinite(thr));
if isempty(thr), return; end
ceiling = -max(thr);
if MispMargin < ceiling, return; end
RT = 1.987E-3*(273.15 + tann);
error('LocalPrimerSearchEngine:MispMarginUnsatisfiable', ...
    ['MispMargin = %.3f kcal/mol cannot be met by ANY %s primer at phi = %.3g.\n' ...
    'The margin is bounded above by -threshold = %.3f, because the anchored free energy is ' ...
    'clamped at zero when no stabilising register exists.\n' ...
    'MispMargin and phi set the same quantity. Demanding %.3f kcal/mol of headroom multiplies ' ...
    'phi by exp(-%.3f/RT) = %.3g at %.1f C, so the pair you have set is equivalent to a single ' ...
    'phi of %.3g.\n' ...
    'Set MispMargin = 0 and choose phi, which is the ratio of off-target duplex to cognate duplex ' ...
    'at one site. Roughly %d anchored registers are scored per primer across both strands, so the ' ...
    'total off-target fraction is about that many times phi.'], ...
    MispMargin, side, phi, ceiling, MispMargin, MispMargin, exp(-MispMargin/RT), tann, ...
    phi*exp(-MispMargin/RT), 122);
end

%% ____________________________________________________________________________________________________________
function tAnn = AnnealLow(input)
% Lowest annealing temperature the programmed cycle actually reaches, degrees C.
%
% A hairpin is most folded at the lowest temperature the primer sees while annealing, so that is the
% temperature the hairpin gate must be judged against. A plain run states one temperature. A
% touchdown recipe [ramp-cycles, ramp, start, constant] descends from start by ramp per cycle for
% ramp-cycles, then holds at constant, so the lowest reached is the smaller of the end of the ramp
% and the constant. Returns NaN when nothing usable is present, which sends LooperMismatch back to
% its free-energy fallback.
tAnn = NaN;
if ~isfield(input, 'AnnealingTemperature'), return; end
v = input.AnnealingTemperature;
if ischar(v)||isstring(v)
    v = str2num(char(v)); %#ok<ST2NM>  a bracketed recipe arrives as text from inputdlg
end
if isempty(v)||any(~isfinite(v)), return; end
if isscalar(v)
    % 0, 1 AND 2 ARE MODE FLAGS, NOT TEMPERATURES.
    % CollectaPCRData_Claude.m prompts for the annealing step as
    % "[0 = Programmed, 1 = Variable, 2: Touchdown]" and stores the answer raw. Mode 2 is then
    % replaced by its four-element recipe, so a scalar 2 never reaches here, but 0 and 1 do. Read
    % literally they set the whole structural screen to 0 or 1 degree: a perfect 21-mer duplex
    % scores -36.75 kcal/mol at 0 C against -13.72 at 51 C, and the hairpin gate rejects most of
    % the pool. There is no temperature to recover for those two modes, since the programme lives
    % outside this struct, so the run stops here rather than proceeding on a number that is not a
    % temperature.
    if any(v == [0 1 2])
        modes = {'Programmed', 'Variable', 'Touchdown'};
        error('LocalPrimerSearchEngine:AnnealingModeFlag', ...
            ['AnnealingTemperature is %g, which is the %s MODE FLAG rather than a temperature. ' ...
            'The structure and mispriming stages need the lowest temperature the reaction ' ...
            'actually reaches. Supply that temperature as a number, or supply a touchdown ' ...
            'recipe as [cycles ramp start constant].'], v, modes{v + 1});
    end
    tAnn = v;
elseif numel(v) >= 4
    cyc = v(1); ramp = v(2); start = v(3); const = v(4);
    % THE RAMP IS SIGNED AND A TOUCHDOWN DESCENDS, so a four-element recipe must carry a ramp that
    % is negative or zero. A positive value walks the ramp UPWARD from start, and the mistake is
    % usually masked: min(..., const) returns const whenever the hold sits below the end of the
    % ramp, which is the ordinary recipe shape, so the right answer comes out for the wrong reason.
    % It surfaces only when the ramp is meant to descend PAST the hold, and there it is worth
    % several degrees at the exact temperature the hairpin gate is judged against. The input is
    % rejected rather than coerced with abs(), which would make an ascending programme impossible
    % to express and would hide the error instead of reporting it.
    if ramp > 0
        error('LocalPrimerSearchEngine:TouchdownRampSign', ...
            ['AnnealingTemperature was given as [cycles ramp start constant] with ramp = %+g. ' ...
            'A touchdown descends, so the ramp is negative: write %+g to step down %g degrees ' ...
            'per cycle from %g.'], ramp, -abs(ramp), abs(ramp), start);
    end
    tAnn = min(start + ramp*max(cyc - 1, 0), const);
else
    tAnn = min(v);
end
end

%% ____________________________________________________________________________________________________________
function d = TerminalBaseDiscrimination(base)
% How well a 3'-terminal base refuses to be extended from its own mispairs.
%
% Huang et al. measured the relative rate at which Taq extends each terminal pair. A primer base can
% form three mispairs, and the geometric mean of epsilon over those three is the single number that
% says how much catalytic discrimination that base buys. The logarithm is taken because the values
% span five decades and the quantity is a rate ratio, so a decade is the natural unit. A LARGER
% return is better discrimination: G 5.18, A 4.96, C 3.84, T 3.84.
%
% The table itself is not duplicated here. ExtensionEfficiency owns it, and is called with the
% terminal pair alone so no sequence context enters.
d = NaN;
base = upper(char(base));
if ~ismember(base, 'ACGT'), return; end
wc = struct('A', 'T', 'C', 'G', 'G', 'C', 'T', 'A');
others = setdiff('ACGT', wc.(base));
acc = 0; n = 0;
for q = 1:numel(others)
    try
        e = ExtensionEfficiency(base, others(q));
    catch
        return    % a failure here must not take a candidate down with it
    end
    if ~isfinite(e)||(e <= 0), continue; end
    acc = acc + log10(e); n = n + 1;
end
if n == 0, return; end
d = -acc/n;
end

%% ____________________________________________________________________________________________________________
function n = KmerCount(kmer, plasdbl, revdbl, plaslen)
% Occurrences of a k-mer across BOTH strands of the circular plasmid. The doubled strings let a site
% straddling the origin be counted; only matches starting inside the first copy are tallied, so no
% site is counted twice.
a = strfind(plasdbl, kmer); b = strfind(revdbl, kmer);
n = sum(a <= plaslen) + sum(b <= plaslen);
end

%% ____________________________________________________________________________________________________________
function [p, dGopen] = OpeningCost(accOn, lunp, posmap, endPos, u, Rgas, TK)
% Free energy cost of opening a primer footprint against template secondary structure. Column u+1 of
% the RNAplfold table holds the probability that the stretch of length u ENDING at that row is
% unpaired, so a footprint is one lookup. Returns 0 when accessibility is switched off or
% unavailable, which reproduces the behaviour of the search before this stage existed.
if ~accOn
    p = NaN; dGopen = 0; return
end
p = lunp(posmap(endPos), u + 1);
if isnan(p)||(p <= 0)
    dGopen = Inf;
else
    dGopen = -Rgas*TK*log(p);
end
end

%% ____________________________________________________________________________________________________________
function [dG3end, dG3stack, dG5stack, prof] = EndStability(prim, endwin)
% Terminal-window stabilities at 37 C, 1 M NaCl. Passing monv = 1 M drives both Owczarzy salt terms
% to zero, which is the reference state the published figures are quoted at.
%
% dG3end   3'-terminal window WITH the terminal and initiation rows, which is the PUBLISHED scale.
%          Primer3's Max 3' Stability is this quantity with the sign flipped: our -6.856 for GCGCG
%          is their 6.86, our -0.877 for TATAT is their 0.86. Verified against both of their worked
%          values. Compare dG3min against this one.
% dG3stack same window, stacking only, which is what this function returned before G9. Kept because
%          the per-locus terminal window of stage 12 is stacking-only and the two must be readable
%          side by side, and because older logs carry it.
% dG5stack the 5'-terminal window, stacking only.
% prof     dG3stack - dG5stack, positive when the 3' end is the LOOSER of the two, which is
%          Rychlik's internal-stability recommendation. Initiation cancels between two stacking-only
%          windows, so this difference needs no scale correction.
%
% DO NOT COMPARE dG3min AGAINST dG3stack. Primer3 includes the helix initiation and the
% stacking-only sum does not, so the two scales sit about 2 kcal/mol apart and a threshold applied
% to the wrong one means something other than what its source measured.
prim = char(prim);
n = strlength(prim);
w3 = prim(max(1, n - endwin + 1):end);
w5 = prim(1:min(n, endwin));
dG3end   = dGMismatch(w3, char(seqcomplement(w3)), 1, 0, 0, 37, false, true,  true);
dG3stack = dGMismatch(w3, char(seqcomplement(w3)), 1, 0, 0, 37, false, false, false);
dG5stack = dGMismatch(w5, char(seqcomplement(w5)), 1, 0, 0, 37, false, false, false);
prof = dG3stack - dG5stack;
end

%% ____________________________________________________________________________________________________________
function alignrow = AlignRow(alignment)
% Middle row of a swalign alignment as a string scalar, "" when the alignment is empty. Keeping every
% textual log field a string scalar is what makes vertical concatenation of records legal: the
% original mixed 'Eliminated' (1x10 char) with 'Selected' (1x8 char) inside one table variable, which
% throws a dimension error on the first accepted candidate.
if isempty(alignment)||(size(alignment, 1) < 2)
    alignrow = "";
else
    alignrow = string(alignment(2, :));
end
end

%% ____________________________________________________________________________________________________________
%% COLLAPSEFORPAIR  Amplicon 3'-OH collapse of both strands for one primer pair, log only
% The plus strand runs 5' to 3' from the earliest admissible forward start to rPos, where its 3'-OH
% sits. The minus strand runs from the latest admissible reverse end down to fPos, read as the
% reverse complement, and ends in its 3'-OH at fPos. A scan depends on the terminus alone, so it is
% made once per rPos (plus) or fPos (minus) and kept in the caches, which are handle objects.
% The evaluation then asks only about folds whose stem lies inside the amplicon, past its 5' end.
function cres = CollapseForPair(fPos, rPos, PlasDbl, PlasLen, input, fFlex, rFlex, W, minRun, tann, useEps, optsPlus, optsMinus, cachePlus, cacheMinus)
monv = input.MonovalentCationConcentration;
biv = input.BivalentCationConcentration;
dNTP = input.dNTPsConcentration;
fEarliest = mod(input.ScaffoldStart - fFlex - 1, PlasLen) + 1;
rLatest = mod(input.ScaffoldEnd + rFlex - 1, PlasLen) + 1;
if ~isKey(cachePlus, rPos)
    len = mod(rPos - fEarliest, PlasLen) + 1;
    scaff = PlasDbl(fEarliest:fEarliest + len - 1);
    cachePlus(rPos) = ScaffoldCollapseCheck('scan', scaff, W, minRun, 1, monv, biv, dNTP, tann, useEps, optsPlus);
end
if ~isKey(cacheMinus, fPos)
    len = mod(rLatest - fPos, PlasLen) + 1;
    scaff = char(seqrcomplement(PlasDbl(fPos:fPos + len - 1)));
    cacheMinus(fPos) = ScaffoldCollapseCheck('scan', scaff, W, minRun, 1, monv, biv, dNTP, tann, useEps, optsMinus);
end
% 1-based index of the amplicon's 5' end on each strand
minR0Plus = mod(fPos - fEarliest, PlasLen) + 1;
minR0Minus = mod(rLatest - rPos, PlasLen) + 1;
cres = struct('plus', ScaffoldCollapseCheck('eval', cachePlus(rPos), minR0Plus, 0), ...
    'minus', ScaffoldCollapseCheck('eval', cacheMinus(fPos), minR0Minus, 0));
end

%% ____________________________________________________________________________________________________________
%% COLLAPSETORECORD  Copy one strand's collapse evaluation into the pair record under a prefix
function record = CollapseToRecord(record, pre, ev)
record.([pre 'Tier']) = string(ev.tier);
record.([pre 'DGfold']) = ev.dGfoldMin;
record.([pre 'DGfoldEng']) = ev.dGfoldEng;
record.([pre 'Tm']) = ev.tmFold;
record.([pre 'TmEng']) = ev.tmEng;
record.([pre 'Loop']) = ev.loopLen;
record.([pre 'Eps']) = ev.eps;
record.([pre 'PExt']) = ev.pExt;
record.([pre 'PExtEng']) = ev.pExtEng;
record.([pre 'JoverC']) = ev.jOverC;
record.([pre 'Sites']) = ev.nSites;
record.([pre 'TierAll']) = string(ev.tierAll);
record.([pre 'PExtTot']) = ev.pExtTot;
record.([pre 'BiShare']) = ev.biShare;
record.([pre 'PExtTaq']) = ev.pExtTaq;
record.([pre 'PExtExo']) = ev.pExtExo;
record.([pre 'PExtBlend']) = ev.pExtBlend;
record.([pre 'TmBi']) = ev.tmBi;
% An infinite free energy stands for no fold found, and the log shows it as empty
if ~isfinite(record.([pre 'DGfold'])), record.([pre 'DGfold']) = NaN; end
if ~isfinite(record.([pre 'DGfoldEng'])), record.([pre 'DGfoldEng']) = NaN; end
end

%% ____________________________________________________________________________________________________________
%% REPORTCOLLAPSE  Tier tally of the accepted pairs, per strand. Log only
function ReportCollapse(logtab, isOn)
if ~isOn||isempty(logtab)||~ismember('cpTier', logtab.Properties.VariableNames), return; end
keep = strcmp(cellstr(string(logtab.Status)), 'Selected');
if ~any(keep), return; end
tiers = {'REJECT', 'FLAG', 'TRANSIENT', 'NONE'};
prefixes = {'cp', 'cm'};
names = {'plus', 'minus'};
fprintf('      Amplicon 3''-OH collapse of %d accepted pairs (log only, nothing rejected)\n', nnz(keep));
for q = 1:2
    col = cellstr(string(logtab.([prefixes{q} 'Tier'])(keep)));
    fprintf('        %-5s strand, foldback in the declared channel:', names{q});
    for t = 1:numel(tiers)
        fprintf('  %s %d', tiers{t}, nnz(strcmp(col, tiers{t})));
    end
    fprintf('\n');
    if ismember([prefixes{q} 'TierAll'], logtab.Properties.VariableNames)
        col = cellstr(string(logtab.([prefixes{q} 'TierAll'])(keep)));
        fprintf('        %-5s strand, both enzyme channels and both topologies:', names{q});
        for t = 1:numel(tiers)
            fprintf('  %s %d', tiers{t}, nnz(strcmp(col, tiers{t})));
        end
        fprintf('\n');
    end
end
fprintf(['      REJECT means a fold with its 3'' end paired persists at the annealing temperature less\n' ...
         '      HpTmMargin. Read cpTm, cpPExt and cpJoverC in Pair_Log before turning any of it into a gate.\n']);
end

%% ____________________________________________________________________________________________________________
%% REPORTREJECTS  Print the stage at which candidates were eliminated, commonest first
% The cascade is strictly ordered, so a candidate is charged to the FIRST gate it failed and the
% counts below say where the search is actually spending its rejections. A stage holding most of
% the population is either doing the real work or standing in for a gate further down; which of
% the two it is can be settled by switching the stage off and re-running, since every gate that
% can be disabled logs its quantity whether or not it is enforced.
function ReportDimerMelts(logtab, tAnneal)
% Dimer melting temperatures for the pairs that survived, each at its own working concentration.
% A dimer is bimolecular, so these are not comparable across concentrations unless each is quoted
% at the concentration its primer actually runs at, which is what DimerMeltTm does. Reported, not
% gated: admission is decided on free energy against phiDimer, which already carries the same
% concentration term.
if isempty(logtab)||~ismember('Status', logtab.Properties.VariableNames), return; end
if ~all(ismember({'ffTm', 'rrTm', 'frTm'}, logtab.Properties.VariableNames)), return; end
keep = strcmp(cellstr(string(logtab.Status)), 'Selected');
if ~any(keep), return; end
sel = logtab(keep, :);
[~, ord] = sort(sel.Score);
sel = sel(ord, :);
shown = min(height(sel), 20);
fprintf('      Dimer melting temperatures, each at its own working concentration\n');
fprintf('      %-24s %-28s %8s %8s %8s %9s\n', ...
    'forward', 'reverse', 'self f', 'self r', 'cross', 'headroom');
for q = 1:shown
    v = [sel.ffTm(q), sel.rrTm(q), sel.frTm(q)];
    hot = v(isfinite(v));
    if isempty(hot), head = NaN; else, head = tAnneal - max(hot); end
    fprintf('      %-24s %-28s %8s %8s %8s %9s\n', ...
        char(string(sel.fseq(q))), char(string(sel.rseq(q))), ...
        TmText(v(1)), TmText(v(2)), TmText(v(3)), TmText(head));
end
if height(sel) > shown
    fprintf('      ... %d further pairs in results.Pair_Log\n', height(sel) - shown);
end
fprintf(['      Headroom is the annealing temperature less the warmest of the three. It is an\n' ...
         '      extrapolation of dG(T) well outside the fitted range when it is large, so read it\n' ...
         '      as distance from the annealing step rather than as a measurable temperature.\n']);
end

%% ____________________________________________________________________________________________________________
function s = TmText(v)
if ~isfinite(v), s = 'none'; else, s = sprintf('%.1f', v); end
end

%% ____________________________________________________________________________________________________________
function ReportRejects(logtab)
if isempty(logtab)||~ismember('RejectStage', logtab.Properties.VariableNames)
    return
end
% cellstr throughout. The log columns may arrive as a string array, a cellstr or a categorical
% depending on how the table was assembled, and cellstr with brace indexing is the one spelling
% that prints correctly under all three.
stage = cellstr(string(logtab.RejectStage));
if ismember('Status', logtab.Properties.VariableNames)
    stage = stage(strcmp(cellstr(string(logtab.Status)), 'Eliminated'));
end
if isempty(stage)
    return
end
[u, ~, g] = unique(stage);
cnt = accumarray(g(:), 1);
[cnt, ord] = sort(cnt, 'descend'); u = u(ord);
tot = sum(cnt);
fprintf('      eliminated %d:', tot);
for q = 1:numel(u)
    fprintf('  %s %d (%.0f%%)', u{q}, cnt(q), 100*cnt(q)/tot);
end
fprintf('\n');
end

%% ____________________________________________________________________________________________________________
function [pen, relYield, extrapolated] = TmAYieldPenalty(gap, gapTab, slopeTab, zeroGap)
%% TMAYIELDPENALTY  Linear-phase yield given up at this amplicon-to-excess-primer gap
%
% WHAT IT RETURNS. pen is ln(slope at the reference gap / slope at this gap), zero when a pair
% gives up nothing and growing as the gap widens. relYield = exp(-pen) is the same thing as a
% fraction, which is the number to quote. extrapolated is true when the gap lies beyond the widest
% value in the table, and a penalty carrying that flag extends a trend rather than reporting a
% measurement.
%
% WHY NOT A DISTANCE FROM A TARGET. The quantity being traded is a rate, and rates compare as
% ratios. A penalty linear in degrees says two degrees cost the same wherever they sit, which
% Pierce's own points contradict: the slope falls by 0.099 per degree between 13 and 17 and by
% 0.201 per degree between 20 and 23, twice as fast. Working in log slope puts the cost where the
% measurement puts it, and makes the weight on this term a statement about yield.
%
% BELOW THE REFERENCE GAP THE PENALTY IS ZERO AND NOT NEGATIVE. The reference is the narrowest gap
% anyone tested and the curve is still climbing there, so a narrower gap is unmeasured rather than
% better. Rewarding it would invent data; a symmetric distance from a target invents the opposite.
%
% NOTHING HERE EXCLUDES A PAIR. A gap of any width returns a finite penalty, so the widest gap on
% a long amplicon is ranked last and still ranked. Whether Pierce's relationship holds at
% multi-kilobase length is unestablished, and an unestablished threshold should order candidates
% rather than delete them.
%
% REFERENCE
%   Pierce KE, Sanchez JA, Rice JE, Wangh LJ. Linear-After-The-Exponential (LATE)-PCR: primer
%       design criteria for high yields of specific single-stranded DNA and improved real-time
%       detection. Proc Natl Acad Sci USA 2005;102:8609-14. doi:10.1073/pnas.0501946102
pen = 0; relYield = 1; extrapolated = false;
if ~isfinite(gap), return; end
g = gapTab(:).'; s = slopeTab(:).';
ok = isfinite(g) & isfinite(s) & (s > 0);
g = g(ok); s = s(ok);
if numel(g) < 2, return; end
[g, ord] = sort(g); s = s(ord);
ls = log(s);
lref = interp1(g, ls, min(max(zeroGap, g(1)), g(end)), 'linear');
if gap <= zeroGap
    pen = 0;
elseif gap <= g(end)
    pen = lref - interp1(g, ls, gap, 'linear');
else
    tail = (ls(end) - ls(end-1))/(g(end) - g(end-1));
    pen = lref - (ls(end) + tail*(gap - g(end)));
    extrapolated = true;
end
pen = max(pen, 0);
relYield = exp(-pen);
end

%% ____________________________________________________________________________________________________________
%% MISPTEMP  The temperature the mispriming screen actually used
% AssertMispSatisfiable asks whether the requested margin is reachable given thresholds that
% MisprimerMismatch computed at ITS evaluation temperature. If the assertion is evaluated at a
% different temperature from the thresholds it is testing, it can declare a satisfiable margin
% unsatisfiable or the reverse. This mirrors the fallback inside the screens exactly.
function t = MispTemp(tAnneal, tm, dTsub)
if ~isempty(tAnneal)&&isfinite(tAnneal)
    t = tAnneal;
else
    t = tm - dTsub;
end
end

%% ____________________________________________________________________________________________________________
function CloseWaitbar(h)
% Invoked through onCleanup, so the waitbar closes on an early return or an error as well as on success.
if ~isempty(h)&&isvalid(h)
    close(h);
end
end
%% ____________________________________________________________________________________________________________
function NonUniq = BuildNonUniqueMaps(seqPlus, seqMinus, kmin, kmax)
%% BUILDNONUNIQUEMAPS  For each k and each plus-strand position, is that k-mer repeated anywhere
%
% NonUniq{k}(p) is true when the k-mer beginning at plus-strand position p occurs MORE THAN ONCE
% across the two strands of the circular template, counting its own occurrence. False therefore
% means the primer sharing that k-mer has exactly one home on the duplex.
%
% WHY BOTH STRANDS ARE COUNTED INTO ONE NUMBER. A primer pairs its template through Watson-Crick
% complementarity, so a contiguous match of the primer against the minus strand IS an occurrence of
% that word on the plus strand, and a contiguous match against the plus strand is an occurrence of
% its reverse complement. Counting the two strands together therefore gives a total that is
% invariant under reverse complementation, and one table indexed by plus-strand coordinate answers
% the question for a forward primer and for a reverse primer alike.
%
% A PALINDROMIC WORD IS COUNTED TWICE BY CONSTRUCTION, and that is physically right rather than an
% artefact. If w equals its own reverse complement then w occurs on the minus strand wherever it
% occurs on the plus strand, so a primer carrying it can pair the SENSE strand at its own location
% as well as the antisense strand it was designed for. That is a real second binding opportunity
% and the count should see it.
%
% THE COUNT IS EXACT ON A CIRCLE. Every start position on each strand is enumerated from a doubled
% string, so a word straddling the origin is counted once and only once: N starts per strand, never
% 2N, and never N - k + 1.
%
% MONOTONICITY, which the callers rely on. If a (k+1)-mer occurs twice then both of its k-mer
% halves occur twice, so NonUniq{k} is true wherever NonUniq{k+1} is. A longest-run search can
% therefore assign k in increasing order and let the larger value overwrite the smaller.
%
% COST. One unique-by-rows over a 2N by k char array per value of k, which is a fraction of a
% second in total on a kilobase template and is paid once for the whole run. Every candidate then
% reads its answer out of an array.
N = numel(seqPlus);
dblP = [seqPlus seqPlus];
dblM = [seqMinus seqMinus];
NonUniq = cell(kmax, 1);
for k = kmin:kmax
    idx = (1:N)' + (0:k-1);              % N by k: every k-mer start on the circle
    W   = [dblP(idx); dblM(idx)];        % 2N by k: the same, on both strands
    [~, ~, g] = unique(W, 'rows');
    cnt = accumarray(g, 1);
    NonUniq{k} = cnt(g(1:N)) >= 2;       % plus-strand rows only, but counted against both strands
end
end

%% ____________________________________________________________________________________________________________
function R = BuildRunLengthMaps(NonUniq, lengths, kmin, kmax, N)
%% BUILDRUNLENGTHMAPS  Longest repeated contiguous word inside, and at the ends of, every window
%
% Three maps, all indexed by the plus-strand coordinate of the window's LOW end.
%   R.Max{L}(p)    longest k, between kmin and L, for which some k-mer inside plus(p : p+L-1) is
%                  repeated elsewhere on the duplex.
%   R.HiEnd{L}(p)  the same restricted to the k-mer ENDING at the window's high end, which is a
%                  forward primer's 3' terminus.
%   R.LoEnd(p)     the same restricted to the k-mer STARTING at the window's low end, which is a
%                  reverse primer's 3' terminus. It does not depend on L; a caller holding a primer
%                  of length L takes min(R.LoEnd(p), L), and monotonicity makes that exactly the
%                  longest repeated terminal word that fits inside the primer.
%
% A value of kmin - 1 means no repeated word of kmin or more, which is reported rather than computed
% because every primer shares short words with a kilobase of template and the number carries no
% information.
%
% The sliding OR is done on a cumulative sum rather than by looping over window positions, so the
% whole map for one k and one L costs two array reads and a subtraction per position.
R = struct();
R.Max   = cell(max(lengths), 1);
R.HiEnd = cell(max(lengths), 1);
R.LoEnd = (kmin - 1)*ones(N, 1);
for k = kmin:kmax
    R.LoEnd(NonUniq{k}) = k;
end
for L = lengths(:)'
    mx = (kmin - 1)*ones(N, 1);
    hi = (kmin - 1)*ones(N, 1);
    for k = kmin:min(L, kmax)
        mx(SlideOr(NonUniq{k}, L - k + 1, N)) = k;
        hi(circshift(NonUniq{k}, -(L - k)))   = k;
    end
    R.Max{L} = mx;
    R.HiEnd{L} = hi;
end
R.kmin = kmin; R.kmax = kmax; R.N = N;
end

%% ____________________________________________________________________________________________________________
function h = SlideOr(a, w, N)
%% SLIDEOR  h(p) is true when any of a(p), a(p+1) ... a(p+w-1) is true, circularly
% A prefix sum over the doubled vector turns a window OR into one subtraction per position, so the
% cost does not grow with the window width.
d  = double([a(:); a(:)]);
cs = [0; cumsum(d)];
p  = (1:N)';
h  = (cs(p + w) - cs(p)) > 0;
end

%% ____________________________________________________________________________________________________________
function [run3, runInt, intOff, intSeq] = RunLengths(RunMap, NonUniq, prim, base, L, N, isFwd)
%% RUNLENGTHS  The two contiguous non-cognate run lengths of one candidate, read from the maps
%
% base is the plus-strand coordinate of the LOW end of the primer's footprint and L its length, so
% the footprint is plus(base : base+L-1) read circularly. A forward primer reads 5' to 3' with the
% plus strand and carries its 3' terminus at the HIGH end of that interval; a reverse primer is the
% reverse complement of it and carries its 3' terminus at the LOW end. Every lookup below turns on
% that single difference.
%
%   run3    longest repeated word CARRYING the primer 3' terminus. Gated at MaxRun3. The low-end
%           map does not depend on the window length, so the reverse value is capped at L, which
%           monotonicity makes exactly the longest terminal repeated word that fits in the primer.
%   runInt  longest repeated word lying wholly inside the primer's first L-1 bases, so clear of the
%           3' terminus by at least one base. Gated at MaxRunInt. It is the same window one base
%           shorter, with the base dropped from whichever end the terminus sits at.
%   intOff  1-based offset of that internal run along the PRIMER, counted 5' to 3', and intSeq the
%           word itself. Both are for the log.
%
% WHY THE OFFSET IS WELL DEFINED. Monotonicity again: if a k-mer occurs twice then so does every
% shorter word inside it, so the first position at which a word of length runInt is repeated is the
% start of a genuine run of exactly that length and not of a longer one.
if isFwd
    run3   = RunMap.HiEnd{L}(base);
    intLow = base;                       % plus(base : base+L-2), terminus dropped from the high end
else
    run3   = min(RunMap.LoEnd(base), L);
    intLow = mod(base, N) + 1;           % plus(base+1 : base+L-1), terminus dropped from the low end
end
Lint = L - 1;
intOff = NaN; intSeq = "";
if Lint < RunMap.kmin
    runInt = RunMap.kmin - 1;
    return
end
runInt = RunMap.Max{Lint}(intLow);
% The offset and the word are for the log alone, so a caller that wants only the two lengths, such
% as the satisfiability count, pays nothing for them and need not supply the primer.
if (nargout < 3)||(runInt < RunMap.kmin)
    return
end
k   = runInt;
idx = mod((intLow:(intLow + Lint - k)) - 1, N) + 1;
hit = find(NonUniq{k}(idx), 1);
if isempty(hit)
    return
end
% The window scanned is plus(intLow : intLow+Lint-1) and the run starts at its hit-th position.
% Mapping that plus coordinate back onto the primer differs by orientation: a forward primer
% numbers with the plus strand, a reverse primer against it.
if isFwd
    intOff = hit;
else
    intOff = L - k + 1 - hit;
end
intSeq = string(prim(intOff:intOff + k - 1));
end

%% ____________________________________________________________________________________________________________
function AssertRunGatesSatisfiable(RunMap, NonUniq, input, N, fFlex, rFlex, fLmin, fLmax, ...
    rLmin, rLmax, use3, maxRun3, useInt, maxRunInt, poolFloor, Verbose)
%% ASSERTRUNGATESSATISFIABLE  Refuse a run whose contiguous-run gates admit nobody
%
% The gates are a joint property of the primer and the TEMPLATE LENGTH, and the dependence is
% exponential. The expected number of occurrences of one specific k-mer elsewhere on a circle of N
% bases, counting both strands, is 2N/4^k, so the same constant is free on a 500 nt construct and
% punishing on a 7 kb plasmid, and which of the two it is cannot be read off the constant.
%
% The test here is the two gates ALONE, with no composition and no thermodynamic screen in front of
% them, so it is an upper bound on what the scan can produce and it costs milliseconds. It is
% measured over the actual search windows rather than over the whole circle, because the local
% search looks at a few dozen positions and the circle-wide statistic would not describe them.
%
% Every forward survivor can pair with every reverse survivor here, the amplicon length being free
% within the two flexes, so the raw pair count is the product and a side that admits nobody makes
% the whole search empty by arithmetic. poolFloor is the number of survivors demanded on EACH side.
fPos = mod((input.ScaffoldStart - fFlex : input.ScaffoldStart) - 1, N) + 1;
rPos = mod((input.ScaffoldEnd : input.ScaffoldEnd + rFlex) - 1, N) + 1;
probe = maxRun3:min(maxRun3 + 5, RunMap.kmax);
if isempty(probe), probe = maxRun3; end
Fp = zeros(size(probe)); Rp = zeros(size(probe));
for q = 1:numel(probe)
    t3 = probe(q);
    ti = t3 + (maxRunInt - maxRun3);     % the two thresholds are probed keeping their separation
    Fp(q) = CountPool(RunMap, NonUniq, fPos, fLmin:fLmax, N, true,  use3, t3, useInt, ti);
    Rp(q) = CountPool(RunMap, NonUniq, rPos, rLmin:rLmax, N, false, use3, t3, useInt, ti);
end
if Verbose
    fprintf(['      run gates  3'' terminal >= %d and internal >= %d rejected: %d forward and ' ...
        '%d reverse candidates survive\n'], maxRun3, maxRunInt, Fp(1), Rp(1));
end
if (Fp(1) >= poolFloor)&&(Rp(1) >= poolFloor), return; end
txt = '';
for q = 1:numel(probe)
    txt = [txt sprintf('\n      MaxRun3 %2d, MaxRunInt %2d :  %5d forward,  %5d reverse', ...
        probe(q), probe(q) + (maxRunInt - maxRun3), Fp(q), Rp(q))]; %#ok<AGROW>
end
good = probe((Fp >= poolFloor)&(Rp >= poolFloor));
if isempty(good)
    advice = ['None of the values tested leaves a pool on both sides of this template. The gates ' ...
        'scale as 2N/4^k, so the other lever is the template: a shorter construct makes the same ' ...
        'thresholds easy. Widening fFlex and rFlex also helps, since it offers more positions.'];
else
    advice = sprintf(['MaxRun3 = %d with MaxRunInt = %d is the smallest pair that leaves a pool ' ...
        'on both sides here. Raising either by one nucleotide multiplies the admitted pool by ' ...
        'roughly four, because that is what each extra base does to the chance of a chance ' ...
        'match.'], good(1), good(1) + (maxRunInt - maxRun3));
end
error('LocalPrimerSearchEngine:RunGatesUnsatisfiable', ...
    ['The contiguous-run gates at MaxRun3 = %d and MaxRunInt = %d admit %d forward and %d ' ...
    'reverse candidates on this %d nt template, against a floor of %d on each side.\n' ...
    'Measured on this template over the actual search windows, gates alone, no other screen:%s\n' ...
    '%s\n' ...
    'A random %d-mer is expected to occur %.3f times elsewhere on a circle this size, and a ' ...
    '20-mer carries %d of them, so carrying none is a rare event rather than a design choice.\n' ...
    'Nothing here is a malfunction. Set the thresholds deliberately, or shorten the template, ' ...
    'and run again.'], ...
    maxRun3, maxRunInt, Fp(1), Rp(1), N, poolFloor, txt, advice, ...
    maxRun3, 2*N/4^maxRun3, 21 - maxRun3);
end

%% ____________________________________________________________________________________________________________
function n = CountPool(RunMap, NonUniq, positions, lens, N, isFwd, use3, t3, useInt, ti)
%% COUNTPOOL  Candidates in one search window that both run gates admit
% positions are plus-strand coordinates of the primer 5' terminus as the scan generates them, which
% is the window's LOW end for a forward primer and its HIGH end for a reverse one.
n = 0;
for L = lens(:)'
    for q = 1:numel(positions)
        if isFwd
            base = positions(q);
        else
            base = mod(positions(q) - L, N) + 1;
        end
        [r3, ri] = RunLengths(RunMap, NonUniq, '', base, L, N, isFwd);
        if (use3&&(r3 >= t3))||(useInt&&(ri >= ti)), continue; end
        n = n + 1;
    end
end
end
