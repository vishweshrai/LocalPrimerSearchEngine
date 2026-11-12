%% PrimerHybridizationMap: Where a primer hybridizes on a circular DNA, drawn as a two-strand annulus.
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 07-11-2026 00:24 IST
%
% The script asks for a circular DNA sequence, a primer and the reaction conditions. It scores every
% gapless register of the primer against both strands of the circle with dGMismatch.m, converts the
% free energies into equilibrium occupancies at the stated primer and template concentrations, and
% draws the result as two concentric annuli.
%
% ____________________________________________________________________________________________________________
% FIGURE
%
%   Outer annulus   (+) strand, the sequence exactly as supplied. Position 1 sits at 12 o'clock and
%                   the strand runs 5'->3' clockwise. Colour map: autumn (red = most stable,
%                   yellow = weakest displayed duplex).
%   Inner annulus   (-) strand. Nucleotide i of this ring is the base paired to (+) nucleotide i, so
%                   the strand runs 5'->3' anticlockwise. Colour map: winter (blue = most stable,
%                   green = weakest displayed duplex).
%   Label dots      Every label in the outer ring carries a dot whose SHAPE names its class, and
%                   the shape is the same one that class uses elsewhere on the figure:
%                   a TRIANGLE for a displayed region, as on the rings; a CIRCLE for an extension
%                   threat, as in the gap; a square, diamond or pentagram for a seed. The dot
%                   takes the colour of the quantity its label quotes.
%   Label text      Three lines: strand and the TEMPLATE range, then the PRIMER positions that
%                   pair there after a vertical bar; the free energy of that duplex; its melting
%                   temperature at the free primer. A threat label quotes the 3'-REACHABLE duplex,
%                   so it will not agree with a region label at the same locus, which quotes the
%                   most stable one. The melting line is dropped when the fit degenerates.
%   Seed labels     Labels whose dot is a SQUARE, DIAMOND or PENTAGRAM mark exact-complement
%                   seeds: a square for a run reaching the primer 3' terminus, a diamond for an
%                   internal run, a pentagram for a site carrying both. These sit in the same
%                   label ring as the region labels, and a locus that is already coloured,
%                   already labelled and already circled still gets one, because the seed is a
%                   separate statement about that locus. See SEEDS below.
%   Markers         Triangles outside the (+) ring and inside the (-) ring mark displayed regions.
%                   Black-edged circles in the gap between the rings mark EXTENSION THREATS: sites
%                   whose duplex reaches the primer 3' terminus to within opts.EndTol and can
%                   therefore be extended. These are scored on that 3'-reachable duplex rather
%                   than on the most stable duplex at the site, weighted by the rate at which the
%                   polymerase extends from that terminus, and marked down to opts.ThreatFactor
%                   times the display threshold, so a circle can sit on grey nucleotides.
%                   opts.EndTol IS AN ENZYME PROPERTY. Answer N to the proofreading question and it
%                   stays at its Taq value of 1. Answer Y and it widens to the whole primer, since
%                   the exonuclease trims any 3' overhang back to the outermost paired base, so
%                   far more sites become threats.
%   Centre          Name and length of the circle, primer and its length, annealing temperature,
%                   monovalent and bivalent cation, total dNTP, primer and template concentrations,
%                   the strongest binding region with its free energy and melting temperature,
%                   the display cutoff, and the extension-efficiency table in force.
%
% Both colour maps share one DeltaG axis, so a given shade means the same free energy on either
% ring. Where displayed duplexes overlap on one strand, the ring shows their union, and every
% nucleotide in the union takes its colour from the most stable duplex that covers it. Small
% triangles mark each displayed region, since a 20 nt site on a 7 kb circle spans about one degree.
%
% ____________________________________________________________________________________________________________
% STRAND AND FRAME CONVENTION
%
% A primer binds the (+) strand when it is complementary to it, which means it reads identical to the
% (-) strand. For register r the primer lies antiparallel over (+) coordinates r .. r+n-1. dGMismatch
% takes the primer 5'->3' as seq and the base physically opposing each primer position as seq2:
%   binding the (+) strand   primer base k faces (+) base r+n-k
%   binding the (-) strand   primer base k faces the complement of (+) base r+k-1
% This is the same frame SeedScan in MisprimerMismatch.m uses. With a perfect match, seq2 is the
% complement of seq and dGMismatch reduces to dGPrim.
%
% ____________________________________________________________________________________________________________
% SCORING
%
%  1. Survey. A register is scored only when it holds at least opts.SurveyRun contiguous Watson-Crick
%     pairs. This is deliberately more permissive than the engine's SURVEYRUN, which is minRun and
%     gates sub-registers rather than whole ones; a map should show what binds, not decide what to
%     reject. The survey is vectorised over all registers of both strands, so it is exhaustive over
%     registers rather than seeded, and a register whose longest run is under 3 can hold at most one
%     nearest-neighbour stack and is never stable.
%  2. Best sub-duplex. The Watson-Crick pairs of a register fall into maximal blocks. Every contiguous
%     union of blocks that holds at least one Watson-Crick stack is scored with dGMismatch, which
%     applies the SantaLucia & Hicks (2004) nearest-neighbour stacks, the Table 2 internal-mismatch
%     increments, the Eq. 12 internal-loop length term through dGLoop, and the Owczarzy 2004/2008
%     salt corrections with Mg2+ chelation by dNTP. The register reports its most stable union, and
%     the ring colours only the nucleotides that union actually spans. Unpaired primer overhangs carry
%     no energy, because dangling-end terms are outside dGMismatch.
%  3. Scores are cached by (primer segment, opposing bases). Pure Watson-Crick segments depend only
%     on the primer, so at most n(n+1)/2 of them are ever computed.
%
% ____________________________________________________________________________________________________________
% EXACT-COMPLEMENT SEEDS
%
% Two criteria that are COMBINATORIAL rather than thermodynamic, and they are kept apart from
% everything above for that reason. A register is coloured because its duplex is stable enough to
% be occupied, and circled because that occupancy times an extension rate clears a bound. A seed
% asks a different question, which no single free energy answers cleanly: is there a stretch of
% PERFECT Watson-Crick pairing long enough to NUCLEATE a duplex? Nucleation is what gets a duplex
% started, and an uninterrupted run of exact pairs is the object that does it.
%
%   3'-terminal run >= opts.Seed3Len (7)      counted back from primer position n, so the run
%       reaches the terminus. This is the dangerous one: it nucleates AND leaves a paired 3'-OH,
%       so the polymerase can extend from it at once with no trimming and no mispaired terminus.
%   Internal run >= opts.SeedIntLen (9)       the longest run lying wholly inside primer positions
%       1..n-1, which is what "internal" means here: the primer with its 3' terminal base removed.
%       Such a site cannot prime as it stands, since the terminus is outside the run. It holds
%       primer at the locus, and a proofreading exonuclease turns it into a primable site by
%       trimming the overhang back to the outermost pair. Nine rather than seven, because an
%       internal run has no immediate catalytic consequence and should not be flagged as readily.
%
% A LONG 3' RUN SATISFIES BOTH, and that is the definition working rather than failing: a run of
% ten reaching the terminus leaves nine inside positions 1..n-1. Such a site is labelled once,
% with a pentagram, and both lengths are printed.
%
% THESE FLAGS SUPPRESS NOTHING AND ARE SUPPRESSED BY NOTHING. A site can be a coloured region, an
% extension threat and a seed at the same time, and it then carries a region label, a threat
% circle and a seed label together. The one exclusion is the cognate site, which satisfies both
% criteria by construction; flagging it would put a seed marker on the intended site in every run
% and say nothing. A SECOND perfect match is not the cognate under the rule above, so it is
% flagged here exactly as it is listed as a threat.
%
% NO FREE ENERGY IS QUOTED FOR A RUN ON ITS OWN. dGMismatch charges a helix initiation at each end
% of whatever it is handed, so a 7-mer scored alone comes back positive for almost every sequence
% and would read as reassuring next to a flag that says otherwise. The meaningful numbers at a
% seed locus are the dG and dGanch of the register carrying it, and both are already reported.
%
% ____________________________________________________________________________________________________________
% MELTING TEMPERATURE OF A SITE
%
% Two are reported per site, one for each free energy already there: Tm_site on the free-edge
% optimum that colours the ring, and Tm_anchored on the 3'-reachable duplex the threat circles are
% scored on. Both are quoted at the FREE primer concentration.
%
% A primer over a template site is bimolecular, so it has no melting temperature until a
% concentration is named. Occupancy is theta = K f/(1 + K f), so theta = 1/2 is K f = 1 and the
% condition is dG(Tm) = R Tm ln(f). That is the same statement as the display cutoff below
% evaluated at theta_min = 0.5, so a site's Tm is the temperature at which it would sit exactly on
% this script's own cutoff with the threshold set to half. Nothing new is asserted; the existing
% machinery is read at a different point.
%
% dH and dS come from a second dGMismatch evaluation of the SAME duplex at Ta + opts.dTm, since
% dG(T) is linear in T for a fixed structure. The relation itself lives in BimolecularTm, which
% DimerMeltTm also calls, so there is one implementation rather than two that can drift. A site
% with no anchored duplex carries dGanch = Inf and its Tm_anchored comes back NaN.
%
% ____________________________________________________________________________________________________________
% EQUILIBRIUM AND DISPLAY CUTOFF
%
% Each scored register is a site with K = exp(-DeltaG/RT). The free primer f follows from
%       [P]0 = f + sum_sites [T]strand * K f/(1 + K f)
% summed over the strands physically present: the (+) strand alone for a circular ssDNA template,
% both strands for dsDNA. Registers are treated as independent sites; with the primer in large excess
% over template, f stays within a fraction of a percent of [P]0 and the approximation does not move
% the cutoff. Occupancy is theta = K f/(1 + K f), and a site is displayed when theta >= theta_min,
% which is the same statement as
%       DeltaG <= DeltaG_cut = -RT ln[ theta_min / ((1 - theta_min) f) ]
% For an ssDNA template the (-) ring therefore shows the occupancy a (-) site would reach once a (-)
% strand has been made, at the same free primer.
%
% THE MAP SHOWS BINDING, NOT PRIMING. A site can be occupied and still be catalytically dead, because
% Taq carries no 3'->5' proofreading exonuclease and will not extend a primer whose 3' terminus lies
% outside the duplex (Huang et al. 1992 NAR 20:4567). Whether the terminus is paired is recorded per
% site as ThreePrimePaired and exported in the Sites sheet; it does not gate the colour. Read the
% rings as a hybridization map and the column as the priming filter.
% THAT PREMISE IS ENZYME-SPECIFIC. A 3'->5' exonuclease excises mismatched nucleotides and removes
% unpaired 3' overhanging nucleotides to create blunt ends, so under a proofreading polymerase the
% overhang stops being a filter: it is trimmed back to the outermost paired base and the primer
% extends from there. Declaring proofreading therefore widens opts.EndTol to the whole primer as
% well as dropping every extension efficiency to 1, and the threat list grows accordingly. Expect a
% primer that looks clean under Taq to show threats under Pfu, Q5, Phusion or KOD. That is the
% model telling the truth about the enzyme, not a change of standard.
%
% Outside the model: bulged registers (use dGmaxMismatch / InternalLoopBridge), primer hairpins and
% self-dimers (LooperMismatch, DimerScan), template secondary structure (DNAAccessibility,
% aPCRAccessibility), re-annealing of the complementary strand in dsDNA, and kinetics.
%
% ____________________________________________________________________________________________________________
% REQUIREMENTS
%
% MATLAB R2019b or newer and no toolboxes. dGMismatch.m and dGLoop.m from the aPCR suite must be on
% the path; the script offers a folder picker when they are missing. exportgraphics (R2020a) is used
% when present, with print as the fallback.
%
% ____________________________________________________________________________________________________________
% REFERENCES
%
%  SantaLucia J Jr, Hicks D (2004) Annu Rev Biophys Biomol Struct 33:415-440.
%       doi:10.1146/annurev.biophys.32.110601.141800
%  SantaLucia J Jr (1998) Proc Natl Acad Sci USA 95:1460-1465. doi:10.1073/pnas.95.4.1460
%  Owczarzy R et al. (2004) Biochemistry 43:3537-3554. doi:10.1021/bi034621r
%  Owczarzy R et al. (2008) Biochemistry 47:5336-5353. doi:10.1021/bi702363u
%  Allawi HT, SantaLucia J Jr (1997) Biochemistry 36:10581-10594. doi:10.1021/bi962590c
%  M13mp18 reference sequence, GenBank X02513.1. https://www.ncbi.nlm.nih.gov/nuccore/X02513.1
%
%

%% 0. Options
clear; close all; clc;
opts = struct();
opts.SurveyRun      = 3;                  % contiguous Watson-Crick pairs before a WHOLE REGISTER is
                                          % scored at all. This is NOT the engine's SURVEYRUN, which
                                          % is minRun (5) and gates SUB-registers inside
                                          % BestSubRegister.m rather than whole ones. A display tool
                                          % should be the more permissive of the two, since its job
                                          % is to show what binds rather than to decide what to
                                          % reject. Dropping a register whose longest run is under 3
                                          % is safe regardless: such a register holds at most one
                                          % nearest-neighbour stack and is never stable.
opts.EndTol         = 1;                  % bases of the primer 3' end allowed OUTSIDE the duplex and
                                          % still counted as extendable. Matches MispEndTol in
                                          % LocalPrimerSearchEngine.m. 0 demands a paired 3' terminus.
opts.UseExtEff      = true;               % weight every anchored site by the rate at which the
                                          % polymerase extends from its 3' REGION, so the threat
                                          % score is theta*epsilon and matches what stage 13 of
                                          % LocalPrimerSearchEngine.m thresholds. Not the terminus
                                          % alone: on a table carrying distance columns, such as
                                          % the shipped extension_efficiency_rejali2018.csv, the
                                          % last three primer positions are each priced at their
                                          % own distance from the 3' end and the smallest of those
                                          % efficiencies is the one applied. A matched terminus one
                                          % base ahead of a mismatch is therefore no longer free.
                                          % false scores threats on occupancy alone.
opts.ThreatFactor   = 0.1;                % extension threats are listed and marked down to this
                                          % multiple of the display threshold, so a site that can
                                          % prime is reported even when it is too weak to colour.
opts.dTm            = 5;                  % degrees C between the two points used to recover dH and
                                          % dS for a site. dG(T) is linear in T for a fixed
                                          % structure, so two evaluations fix the pair exactly and
                                          % the melting temperature follows from BimolecularTm. The
                                          % second point is taken only for the two duplexes that are
                                          % reported, the free-edge optimum and the 3'-anchored one,
                                          % rather than for every sub-duplex the survey considers.
opts.ThreatMax      = 15;                 % most threats listed on the console and marked on the map
opts.ThreatLabels   = 5;                  % most threats given a LABEL on the ring, taken in priming
                                          % order. Every threat up to ThreatMax still gets its
                                          % circle in the gap between the rings; this bounds only
                                          % the labelled subset, because the label ring is shared
                                          % with the region and seed labels and SpreadAngles can
                                          % only separate so many before they smear
opts.ExactSubDuplex = true;               % true  places sub-duplex edges at EVERY Watson-Crick pair,
                                          %       which is what BestSubRegister.m in the engine does.
                                          % false places them only at block boundaries: about three
                                          %       times faster, and it misses the optimum in roughly
                                          %       4% of stable registers by up to 0.14 kcal/mol at
                                          %       60 C (measured, 142 registers below -6 kcal/mol).
opts.RingPlus       = [0.875 1.000];      % (+) annulus, inner and outer radius, data units
opts.RingMinus      = [0.705 0.830];      % (-) annulus, inner and outer radius, data units. Both
                                          % rings are 0.125 thick and the 0.045 gap between them
                                          % carries the extension-threat markers. Thicker than the
                                          % 0.10 they used to be: on a 7 kb circle a 20 nt site is
                                          % about one degree of arc, so the only thing that makes
                                          % its colour legible is radial depth.
opts.TextRadius     = 0.655;              % centre text shrinks until every line fits inside this
                                          % circle. Held just inside the (-) ring, so thickening
                                          % the rings has not cost the centre block any room.
opts.FigureSizeCm   = [30 23];
opts.CentreFontSize = 11.5;               % starting size of the centre block, not its final one
opts.CentreFontMax  = 16;                 % how far LayoutCentreText may GROW it to fill
                                          % TextRadius. The old loop only shrank, so the panel was
                                          % as large as the constant above said rather than as
                                          % large as the circle allowed, and a short caption left
                                          % the middle of the figure half empty. Geometry usually
                                          % binds before this does.
opts.TitleScale     = 1.55;               % the template name, as a multiple of CentreFontSize
opts.LabelFontSize  = 11.5;
opts.TickFontSize   = 11;
opts.LegendFontSize = 11;                 % colour-bar ticks and the two swatch captions. Was 9.
opts.BarFontSize    = 12;                 % colour-bar titles and axis label. Was 10.
opts.LabelMinSepDeg = 8;                  % FLOOR on the angular spacing between labels, degrees.
                                          % The working value is computed from how tall a label
                                          % actually is, by LabelSepDeg below, and this only
                                          % raises it. A fixed number cannot be right for both a
                                          % two-line and a three-line label.
opts.LabelLineFactor = 1.30;              % line pitch of a multi-line text object, as a multiple
                                          % of the font size, plus a little clearance. MATLAB does
                                          % not expose the pitch it uses, so this is the one
                                          % number in the separation that is estimated. Raise it
                                          % if labels still touch, lower it if they sit too far
                                          % apart
opts.UnboundColour  = [0.90 0.90 0.90];
opts.ColourLevels   = 256;
opts.ColourFloor    = [];                 % kcal/mol; [] spans the colours from the strongest shown site
opts.MaxMarkers     = 25;                 % above this many regions on a strand, only labelled ones get a triangle
% EXACT-COMPLEMENT SEEDS. Two COMBINATORIAL criteria, deliberately not thermodynamic.
% Everything else on this map is an energy: a register is coloured because its duplex is stable
% enough to be occupied, and marked as a threat because that occupancy times an extension rate
% clears a bound. A seed asks a different question, which no free energy answers cleanly at a
% single temperature: is there a stretch of PERFECT Watson-Crick pairing long enough to nucleate?
% Nucleation is what gets a duplex started, and a run of exact pairs is the object that does it.
% The two runs are kept apart because they fail the reaction in different ways.
%   Seed3Len    a run ENDING at the primer 3' terminus. This is the dangerous one: it both
%               nucleates and leaves a paired 3'-OH, so the polymerase can extend from it
%               immediately. Seven is the usual working figure for a 3' anchor.
%   SeedIntLen  the longest run lying wholly inside primer positions 1..n-1, the primer with its
%               3' terminal base removed, which is what "internal" means here. Such a site cannot
%               prime on its own, since the terminus is outside the run, but it holds primer at
%               the locus and is the structure a proofreading exonuclease converts into a primable
%               one by trimming the overhang. Nine is required rather than seven because an
%               internal run has no catalytic consequence and should not be flagged as readily.
% NO SEED FREE ENERGY IS QUOTED, AND THAT IS DELIBERATE. dGMismatch charges a helix initiation at
% each end of whatever it is handed, so a 7-mer scored alone comes back positive for almost every
% sequence and would read as "harmless" next to a flag that says otherwise. The physically
% meaningful numbers at a seed locus are the register's own dG and dGanch, and both are already
% reported. The seed is reported as what it is: a count of contiguous pairs.
opts.Seed3Len       = 7;                  % contiguous WC pairs ending at the primer 3' terminus
opts.SeedIntLen     = 9;                  % longest contiguous WC run inside primer 1..n-1
opts.SeedMax        = 15;                 % most seed sites labelled on the map and listed
R_GAS = 8.314e-3/4.184;                   % kcal/(mol K). Matched to dGMismatch.m and the threshold
                                          % functions, NOT to the 1.98720e-3 of dGLoop.m and
                                          % PrimTemp.m, which differ in the sixth figure. The free
                                          % energies this script converts come from dGMismatch, and
                                          % RT now also multiplies log(epsilon), which reaches -13.8
                                          % at the detection limit, so the two constants are made
                                          % identical rather than merely close.

dflt = struct('Name', 'M13mp18', 'Strands', '2', 'Primer', 'GACGCTTACTGGAACGGTAAA', ...
    'Ta', '60', 'Monv', '40', 'Biv', '2.5', 'dNTP', '800', 'PrimerConc', '1000', ...
    'TemplateValue', '1', 'TemplateUnit', 'ng/uL', 'ThetaPct', '0.01', 'NumLabels', '5', ...
    'Proofread', 'none');

%% 1. Dependencies
need = {'dGMismatch', 'dGLoop', 'ExtensionEfficiency', 'ExtensionEfficiencyTable', ...
    'BimolecularTm'};
% ExtensionEfficiencyTable is listed because ExtensionEfficiency reaches for it whenever the
% caller names neither a table nor a file, which is how this script calls it. Without the file
% on the path the failure would surface deep inside the threat loop rather than here.
missing = need(cellfun(@(f) isempty(which(f)), need));
if ~isempty(missing)
    apcrDir = uigetdir(pwd, 'Select the aPCR folder holding the suite''s .m files');
    if isequal(apcrDir, 0)
        error('PrimerHybridizationMap:MissingDependency', 'Not on the MATLAB path: %s.', ...
            strjoin(missing, ', '));
    end
    addpath(apcrDir);
    missing = need(cellfun(@(f) isempty(which(f)), need));
    if ~isempty(missing)
        error('PrimerHybridizationMap:MissingDependency', 'Still not on the MATLAB path: %s.', ...
            strjoin(missing, ', '));
    end
end

%% 2. Circular DNA
src = questdlg('How should the circular DNA sequence be supplied?', 'Circular DNA', ...
    'Load a file', 'Paste the sequence', 'Cancel', 'Load a file');
switch src
    case 'Load a file'
        [fn, fp] = uigetfile({'*.txt;*.fa;*.fasta;*.fna;*.seq;*.gb;*.gbk', ...
            'Sequence files (*.txt, *.fasta, *.gb)'; '*.*', 'All files (*.*)'}, ...
            'Select the circular DNA sequence', 'D:\OneDrive - Indian Institute of Science\Codes\Developed_Scripts\aPCR\Plasmids\');
        if isequal(fn, 0)
            fprintf('No sequence file selected.\n');
            return
        end
        [~, fileStem] = fileparts(fn);
        [plusSeq, dflt.Name] = ParseSequenceText(fileread(fullfile(fp, fn)), fileStem);
    case 'Paste the sequence'
        pasted = inputdlg({'Name of the circular DNA:', ...
            'Sequence, 5''->3'' (plain text, FASTA or GenBank):'}, 'Paste the circular DNA', ...
            [1 90; 14 90], {dflt.Name, ''}, struct('Resize', 'on'));
        if isempty(pasted)
            return
        end
        [plusSeq, parsedName] = ParseSequenceText(strjoin(cellstr(pasted{2}), newline), ...
            strtrim(pasted{1}));
        if ~isempty(strtrim(pasted{1}))
            dflt.Name = strtrim(pasted{1});
        else
            dflt.Name = parsedName;
        end
    otherwise
        return
end
bad = ~ismember(plusSeq, 'ACGT');
if any(bad)
    firstBad = find(bad, 1);
    error('PrimerHybridizationMap:NonACGT', ...
        ['The circular DNA carries %d character(s) outside {A,C,G,T}; the first sits at position %d ' ...
        'and is ''%s''. The nearest-neighbour tables are defined on unambiguous DNA only.'], ...
        nnz(bad), firstBad, plusSeq(firstBad));
end
L = numel(plusSeq);
if L < 20
    error('PrimerHybridizationMap:TooShort', 'The circular DNA is %d nt long; at least 20 nt are needed.', L);
end

%% 3. Primer and reaction conditions
prompt = {'Name of the circular DNA:', ...
    'Template strands: 1 = circular ssDNA, 2 = circular dsDNA:', ...
    'Primer sequence, 5''->3'':', ...
    'Annealing temperature (deg C):', ...
    'Monovalent cation concentration (mM):', ...
    'Bivalent cation concentration (mM):', ...
    'Total dNTP concentration (uM):', ...
    'Primer concentration (nM):', ...
    'Template concentration, value:', ...
    'Template concentration, unit (nM or ng/uL):', ...
    'Display threshold, equilibrium occupancy (%):', ...
    'Regions to label on each strand:', ...
    'Polymerase 3''->5'' proofreading exonuclease (none / weak / strong):'};
definput = {dflt.Name, dflt.Strands, dflt.Primer, dflt.Ta, dflt.Monv, dflt.Biv, dflt.dNTP, ...
    dflt.PrimerConc, dflt.TemplateValue, dflt.TemplateUnit, dflt.ThetaPct, dflt.NumLabels, ...
    dflt.Proofread};
rx = inputdlg(prompt, sprintf('Reaction conditions for %s (%d nt)', dflt.Name, L), ...
    repmat([1 70], numel(prompt), 1), definput, struct('Resize', 'on'));
if isempty(rx)
    return
end
dnaName = strtrim(rx{1});
if isempty(dnaName), dnaName = dflt.Name; end
nStrands = ReadNumber(rx{2}, 'Template strands', 1, 2);
if ~ismember(nStrands, [1 2])
    error('PrimerHybridizationMap:Strands', 'Template strands must be 1 or 2.');
end
primer = strrep(upper(regexprep(rx{3}, '[^A-Za-z]', '')), 'U', 'T');
n = numel(primer);
if any(~ismember(primer, 'ACGT'))
    error('PrimerHybridizationMap:PrimerNonACGT', 'The primer carries characters outside {A,C,G,T}.');
end
if (n < 6)||(n >= L)
    error('PrimerHybridizationMap:PrimerLength', 'The primer must be 6 to %d nt long; it is %d nt.', L - 1, n);
end
Ta       = ReadNumber(rx{4}, 'Annealing Temperature', 0, 100);
monv_mM  = ReadNumber(rx{5}, 'Monovalent Cation', eps, 2000);
biv_mM   = ReadNumber(rx{6}, 'Bivalent Cation', 0, 1000);
dNTP_uM  = ReadNumber(rx{7}, 'Total dNTP', 0, 1e5);
p0_nM    = ReadNumber(rx{8}, 'Primer Concentration', eps, 1e7);
tValue   = ReadNumber(rx{9}, 'Template Concentration', eps, 1e7);
tUnit    = lower(regexprep(rx{10}, '\s', ''));
thetaPct = ReadNumber(rx{11}, 'Display Threshold', eps, 100 - eps);
nLabels  = round(ReadNumber(rx{12}, 'Regions to Label', 0, 50));

% A 3' -> 5' PROOFREADING EXONUCLEASE CHANGES WHAT A MISPAIRED 3' REGION MEANS. Neither source of
% extension efficiencies was measured on a proofreading enzyme: Huang et al. 1992 used Taq and
% Rejali et al. 2018, whose table is the shipped default, used KlenTaq 1, the N-terminal deletion
% mutant. Taq carries no 3' -> 5' exonuclease and the deletion removes the 5' -> 3' one, so neither
% enzyme can excise a mispaired nucleotide and both simply fail to extend past it. An enzyme that
% proofreads excises it and carries on, so every epsilon collapses towards 1 and the threat score
% has to fall back to occupancy alone. New England Biolabs describe the mechanism: the domain
% excises mismatched nucleotides in the 3' to 5' direction and removes unpaired 3' overhanging
% nucleotides to create blunt ends, which is the geometry a mispaired primer terminus presents.
% https://www.neb.com/en/products/pcr-qpcr-and-amplification-technologies/dna-polymerase-proofreading
% Answer Y for Pfu, Vent (exo+), Phusion, Q5 and KOD; N for Taq, Stoffel Fragment and Deep Vent
% (exo-). This is the same question CollectaPCRData_Claude.m asks, and it must be answered the same
% way in both places or the map and the search engine will rank sites differently.
% THE CLASS, NOT A FLAG. LocalPrimerSearchEngine.m carries none, weak and strong through
% ProofExoClass and derives its own flag as ~strcmp(class, 'none'), so weak and strong condition
% every gate identically there. They are kept distinct because strong carries a warning about
% single-stranded YIELD that neither program models. This script follows that division exactly, or
% the map and the engine would describe different enzymes. The older Y and N spellings are accepted
% so an operator who types them still gets a defined answer.
proofClass = lower(strtrim(rx{13}));
switch proofClass
    case {'y', 'yes', '1', 'true'},  proofClass = 'weak';
    case {'n', 'no', '0', 'false'},  proofClass = 'none';
    case {'none', 'weak', 'strong'}
    otherwise
        error('PrimerHybridizationMap:BadProofExoClass', ...
            ['The proofreading answer must be ''none'', ''weak'' or ''strong''. Received ''%s''. ' ...
            'This is the same vocabulary PrimerSearchParameters.m uses.'], proofClass);
end
proofread = ~strcmp(proofClass, 'none');
if proofread
    opts.UseExtEff = false;
    % THE END TOLERANCE IS ALSO AN ENZYME PROPERTY, and this is the half that was missing.
    % opts.EndTol says how many primer 3' bases may hang OUTSIDE the duplex and still count as
    % extendable. The value 1 is a Taq number. It is small because Taq cannot remove an overhang,
    % so a register that leaves the terminus dangling is catalytically dead however tightly it
    % binds. The same NEB description that justifies dropping epsilon says the proofreading domain
    % removes unpaired 3' overhanging nucleotides to create blunt ends, which is precisely the
    % geometry EndTol is measuring. Under such an enzyme the overhang is trimmed back to the
    % outermost Watson-Crick pair and the primer extends from there, so NO overhang length
    % disqualifies a register and the tolerance is the whole primer.
    % This is not the same as ranking on the free-edge optimum. The exonuclease stops at the first
    % paired base; it does not keep chewing into the helix. ScanStrand already builds exactly that
    % object when the tolerance admits the register: right edge pinned at the outermost pair, left
    % edge free. Widening EndTol therefore gets the right duplex without a second code path.
    % Trimming does not weaken what remains, because the overhanging bases were contributing no
    % stacking to begin with. What stands between a stable unanchored register and a product is
    % whether trimming outruns dissociation, and that is kinetics. NOT MODELLED. Admitting every
    % register is the conservative reading, and it matches the move LocalPrimerSearchEngine.m makes
    % at H13 when it puts dGEndMin on the whole-molecule dimer minimum.
    opts.EndTol = n - 1;
    warning('PrimerHybridizationMap:ProofreadingPolymerase', ...
        ['A proofreading polymerase was declared, so TWO things changed. Every extension ' ...
        'efficiency is set to 1, and the end tolerance was widened from its Taq value to %d, the ' ...
        'whole primer, because the exonuclease removes unpaired 3'' overhanging nucleotides and ' ...
        'primes from the outermost paired base. Threats are therefore ranked on occupancy of the ' ...
        'trimmed duplex and no overhang length disqualifies a site. This is the conservative ' ...
        'reading for an enzyme that excises the mispair; whether trimming outruns dissociation is ' ...
        'kinetics and is not modelled. Such an enzyme also nibbles free single-stranded primers ' ...
        'over a long run, shortening them and lowering their Tm; nothing here models that ' ...
        'either.'], opts.EndTol);
    if strcmp(proofClass, 'strong')
        warning('PrimerHybridizationMap:StrongProofreadingYield', ...
            ['The class declared was STRONG. Veneziano et al. report that Phusion "was unable to ' ...
            'synthesize fragment larger than 1 kb ssDNA", which they attribute to "enhanced ' ...
            'exonuclease activity", and that "Taq polymerases gave the highest ssDNA ' ...
            'production", "possibly due to lower exonuclease activity". Both are the authors'' ' ...
            'own hypotheses rather than measurements. This warning is about YIELD, which nothing ' ...
            'here models; the map is drawn exactly as for a weak exonuclease. Veneziano R et al. ' ...
            'Sci Rep 2018;8:6548, doi:10.1038/s41598-018-24677-5.']);
    end
end

% WHICH EXTENSION-EFFICIENCY TABLE IS IN FORCE. ExtensionEfficiency falls back to the session's
% choice when the caller names neither a matrix nor a file, which is how this script calls it. That
% choice is set by PrimerSearchParameters during an engine run and survives in
% ExtensionEfficiencyTable, so a map drawn after such a run can be scored on a table this script
% never mentions. Two differences between the tables are large enough to reorder a threat list.
% G.G alone differs by a factor near 880 between the built-in Huang bounds and the Rejali
% measurement. The second matters more: the Huang table carries NO distance columns, so on it a
% mismatch one or two bases behind the terminus costs nothing at all, while on the Rejali table it
% costs between 0.70 and 5.07 kcal/mol at 58 C. A map drawn on the built-in table after a search
% run on the shipped one would therefore keep threats the search displaced, and keep them at full
% strength.
%
% WHICH IS WHY THE SESSION CHOICE IS NO LONGER THE LAST WORD. This script is normally run on its
% own, in a session where no engine run has set that persistent, and until 07-11-2026 that meant
% it dropped past the choice to the built-in table every time. ExtensionEfficiency now resolves
% extension_efficiency_rejali2018.csv from its own folder when nothing else names a table, so a
% map drawn cold prices sites the way a search would. The built-in is reached only by asking for
% it through ExtensionEfficiencyTable('set', 'builtin'), or if that CSV has gone missing, and
% either way it says which of the two happened. The source is resolved once here on a matched pair
% and reported with everything else, so the figure cannot claim a provenance it does not have.
if opts.UseExtEff
    [~, epsInfo] = ExtensionEfficiency('AAA', 'TTT');
    effSource = epsInfo.Source;
    effCite   = epsInfo.Citation;
else
    effSource = 'not applied';
    effCite   = 'not applied';
end

%% 4. Derived quantities
monv = monv_mM*1E-3; biv = biv_mM*1E-3; dNTP = dNTP_uM*1E-6; p0 = p0_nM*1E-9;

% Molar mass of the circle. Residue masses are those of CollectaPCRData_Claude.m. A covalently closed
% circle has no 5' or 3' terminus, so the +79.0 g/mol terminal term used there for linear DNA is left
% out; at 7249 nt the two conventions differ by 0.004%.
resMass = [313.2 289.2 329.2 304.2];                              % A C G T, g/mol
baseCnt = [sum(plusSeq == 'A'), sum(plusSeq == 'C'), sum(plusSeq == 'G'), sum(plusSeq == 'T')];
MWplus  = baseCnt*resMass';
MWminus = fliplr(baseCnt)*resMass';                               % (-) composition is the complement
if nStrands == 1
    MW = MWplus;  topoTxt = 'ssDNA';
else
    MW = MWplus + MWminus;  topoTxt = 'dsDNA';
end
switch tUnit
    case 'nm'
        cT = tValue*1E-9;
    case {'ng/ul', ['ng/' char(181) 'l'], 'ngul'}
        cT = tValue/(MW*1E3);                                     % 1 ng/uL = 1E-3 g/L
    otherwise
        error('PrimerHybridizationMap:TemplateUnit', 'Template unit must be nM or ng/uL, not ''%s''.', rx{10});
end
cT_nM = cT*1E9;
ngul  = cT*MW*1E3;

% Free Mg2+ after chelation by dNTP, Ka = 3E4 M^-1, the same expression dGPrim.m and dGMismatch.m
% evaluate internally. Reported in the centre text only; the salt correction recomputes it.
KaMg = 3E4;
freeBiv = (-(KaMg*dNTP - KaMg*biv + 1) + sqrt((KaMg*dNTP - KaMg*biv + 1)^2 + 4*KaMg*biv))/(2*KaMg);

%% 5. Register scan on both strands
rxn = struct('monv', monv, 'biv', biv, 'dNTP', dNTP, 'Ta', Ta, 'R', R_GAS);
cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
ws = warning('off', 'dGMismatch:TandemMismatch');
hw = waitbar(0, 'Scanning the (+) strand...', 'Name', 'Primer hybridization map');
try
    seedMin = min(opts.Seed3Len, opts.SeedIntLen);
    sitesP = ScanStrand(plusSeq, primer, '+', rxn, opts.SurveyRun, opts.ExactSubDuplex, opts.EndTol, opts.UseExtEff, cache, hw, 0.0, 0.5, seedMin, opts.dTm);
    sitesM = ScanStrand(plusSeq, primer, '-', rxn, opts.SurveyRun, opts.ExactSubDuplex, opts.EndTol, opts.UseExtEff, cache, hw, 0.5, 0.5, seedMin, opts.dTm);
catch scanErr
    warning(ws);
    if ishghandle(hw), delete(hw); end
    rethrow(scanErr);
end
warning(ws);
if ishghandle(hw), delete(hw); end
S = CatSites(sitesP, sitesM);
nSites = numel(S.dG);

%% 6. Equilibrium occupancy and display cutoff
Tk = Ta + 273.15; RT = R_GAS*Tk;
K = exp(-S.dG/RT);
present = (S.Strand == '+')|(nStrands == 2);
Kp = K(present);
massBalance = @(u) u + sum((cT/p0)*(Kp*p0*u)./(1 + Kp*p0*u)) - 1;   % u = free primer / total primer
uFree = fzero(massBalance, [0 1], optimset('TolX', 1E-14));
fFree = uFree*p0;
% Overlapping registers are counted as independent sites, which is only defensible while almost no
% primer is bound. The header asserts that; this measures it rather than assuming it.
boundFrac = 1 - uFree;
if boundFrac > 0.01
    warning('PrimerHybridizationMap:PrimerDepletion', ...
        ['%.2f%% of the primer is bound at equilibrium. Registers are treated as independent ' ...
        'sites, but overlapping registers on one template molecule are mutually exclusive, so ' ...
        'the free primer is underestimated and every occupancy with it. Raise the primer ' ...
        'concentration or lower the template concentration to keep this below 1%%.'], 100*boundFrac);
end
% MELTING TEMPERATURES. A primer over a template site is bimolecular, so the site has no
% melting temperature until a concentration is named. Occupancy is theta = K*f/(1 + K*f), so
% theta = 1/2 is K*f = 1, and the melting condition is dG(Tm) = R*Tm*ln(f) with f the FREE primer.
% This is not a new assertion: dGcut below is -RT*log(thetaMin/((1 - thetaMin)*fFree)), which at
% thetaMin = 0.5 is exactly RT*log(fFree). A site's Tm is therefore the temperature at which it
% would sit on this script's own display cutoff with the threshold set to half.
% Free rather than total primer, because f is what the sites actually see. The two differ by the
% bound fraction reported above, which the depletion warning keeps under one per cent.
% Two temperatures per site, so BimolecularTm converts the whole list in one call. Sites with no
% anchored duplex carry dGanch = Inf and come back as NaN.
S.Tm     = BimolecularTm(S.dG,     S.dG2,     Ta, Ta + opts.dTm, fFree, R_GAS);
S.TmAnch = BimolecularTm(S.dGanch, S.dGanch2, Ta, Ta + opts.dTm, fFree, R_GAS);

theta = 1./(1 + 1./(K*fFree));
thetaMin = thetaPct/100;
dGcut = -RT*log(thetaMin/((1 - thetaMin)*fFree));
shown = theta >= thetaMin;

% Extension threats. theta above is the occupancy of the most stable duplex at a site, which may
% leave the primer 3' terminus outside the helix and so cannot prime. thetaAnch is the occupancy of
% the 3'-ANCHORED duplex at the same site, which can. A perfect full-length register is the primer's
% own cognate site and is reported separately rather than as a threat.
Kanch     = exp(-S.dGanch/RT);                 % dGanch = Inf where no anchored duplex exists
thetaAnch = 1./(1 + 1./(Kanch*fFree));
thetaAnch(~isfinite(thetaAnch)) = 0;

% THE THREAT SCORE IS theta*epsilon, NOT theta. Priming requires the primer to be there AND the
% polymerase to extend from it, and those two factors multiply. primeRate is the rate of productive
% priming at a site relative to a saturated, perfectly matched one, so the cognate sits near 1 and
% everything else is read against that. This is the same quantity stage 13 of
% LocalPrimerSearchEngine.m now thresholds, so the two programs rank sites the same way.
primeRate = thetaAnch.*S.Eps;
threatMin = thetaMin*opts.ThreatFactor;
% A PERFECT MATCH IS NOT AUTOMATICALLY THE COGNATE SITE.
% The cognate site used to be identified as "every primer position paired", which is a property of
% the SEQUENCE and not of the position. If a primer has two full-length perfect matches on the
% circle, the second one is the single most dangerous off-target event that can exist: it binds as
% tightly as the intended site and primes just as well. Filing it under COGNATE excluded it from
% the threat list and from the markers.
% Note that "most stable perfect match" cannot break the tie either, because a full-length perfect
% duplex is scored from the primer alone and every such site returns exactly the same free energy.
% So when more than one exists the honest response is to say the primer is not unique, keep ONE as
% the cognate and treat the rest as threats.
kPerfect = find(S.Perfect);
isCognate = false(size(S.Perfect));
if ~isempty(kPerfect)
    isCognate(kPerfect(1)) = true;
    if numel(kPerfect) > 1
        loc = arrayfun(@(k) sprintf('(%c) %d-%d', S.Strand(k), S.StartPos(k), S.EndPos(k)), ...
            kPerfect, 'UniformOutput', false);
        warning('PrimerHybridizationMap:PrimerNotUnique', ...
            ['This primer matches the template PERFECTLY at %d places: %s. Every full-length ' ...
            'perfect duplex has the same free energy, so none of them can be singled out as the ' ...
            'intended site on thermodynamics. The first is reported as cognate and the rest are ' ...
            'listed as extension threats. A primer that is not unique cannot give a clean ' ...
            'asymmetric reaction whatever the rest of the map looks like.'], ...
            numel(kPerfect), strjoin(loc', ', '));
    end
end
isThreat  = S.Anchored&~isCognate&(primeRate >= threatMin);
[~, ordThreat] = sort(primeRate, 'descend');
ordThreat = ordThreat(isThreat(ordThreat));
ordThreat = ordThreat(1:min(opts.ThreatMax, numel(ordThreat)));

% EXACT-COMPLEMENT SEED SITES.
% Purely combinatorial, so these are flagged independently of every threshold above and a site can
% carry a seed flag together with a displayed-region colour, a threat circle, or both. Nothing is
% suppressed because a site is already marked some other way: the seed is a separate statement
% about the same locus and it gets its own label.
% THE COGNATE SITE IS EXCLUDED, for the same reason it is excluded from the threat list. A perfect
% full-length match satisfies both criteria by construction, so including it would put a seed
% marker on the intended site in every run and say nothing. A SECOND perfect match is not the
% cognate under the rule above, so it is flagged here exactly as it is flagged as a threat.
isSeed3   = (S.End3Run >= opts.Seed3Len)&~isCognate;
isSeedInt = (S.IntRun  >= opts.SeedIntLen)&~isCognate;
isSeedAny = isSeed3|isSeedInt;
% Ranked 3'-terminal first, since only those carry a free 3'-OH, then by how long the run is.
seedRank = zeros(nSites, 1);
seedRank(isSeedInt) = S.IntRun(isSeedInt);
seedRank(isSeed3)   = 1000 + S.End3Run(isSeed3);
[~, ordSeed] = sort(seedRank, 'descend');
ordSeed = ordSeed(isSeedAny(ordSeed));
nSeedAll = numel(ordSeed);
ordSeed  = ordSeed(1:min(opts.SeedMax, nSeedAll));

if nSites > 0
    [~, iBest] = min(S.dG);
else
    iBest = [];
end

%% 7. Per-nucleotide union maps
bestP = inf(1, L); bestM = inf(1, L);
% The melting temperature travels with the free energy rather than being recomputed from the union.
% A region is the union of overlapping duplexes and has no melting temperature of its own, so each
% nucleotide carries the Tm of whichever single duplex gave it its dG. The region then quotes the Tm
% of the duplex that set the region's minimum, which is the site the label's dG already describes.
bestPTm = nan(1, L); bestMTm = nan(1, L);
% The primer register travels the same way and for the same reason. A region spans whatever its
% member duplexes cover, and those need not use the same stretch of the primer, so the register
% quoted on a label is the one belonging to the duplex that set the region's free energy.
bestPpf = nan(1, L); bestPpt = nan(1, L);
bestMpf = nan(1, L); bestMpt = nan(1, L);
for k = find(shown)'
    pos = S.Positions{k};
    if S.Strand(k) == '+'
        hit = S.dG(k) < bestP(pos);
        bestP(pos(hit))   = S.dG(k);
        bestPTm(pos(hit)) = S.Tm(k);
        bestPpf(pos(hit)) = S.PrimerFrom(k);
        bestPpt(pos(hit)) = S.PrimerTo(k);
    else
        hit = S.dG(k) < bestM(pos);
        bestM(pos(hit))   = S.dG(k);
        bestMTm(pos(hit)) = S.Tm(k);
        bestMpf(pos(hit)) = S.PrimerFrom(k);
        bestMpt(pos(hit)) = S.PrimerTo(k);
    end
end
regP = RegionTable(bestP, bestPTm, bestPpf, bestPpt, '+');
regM = RegionTable(bestM, bestMTm, bestMpf, bestMpt, '-');

%% 8. Figure
minusSym = ['(' char(8722) ')'];             % U+2212, a real minus rather than a hyphen

mapP = autumn(opts.ColourLevels);            % red (strongest) -> yellow (weakest shown)
mapM = winter(opts.ColourLevels);            % blue (strongest) -> green (weakest shown)
gHi = dGcut;
if ~isempty(opts.ColourFloor)
    gLo = opts.ColourFloor;                      % sites below the floor take the extreme colour
elseif any(shown)
    gLo = floor(2*min(S.dG(shown)))/2;
else
    gLo = gHi - 5;
end
if gHi - gLo < 1, gLo = gHi - 1; end
cIdx = @(g) min(max(1 + round((g - gLo)/(gHi - gLo)*(opts.ColourLevels - 1)), 1), opts.ColourLevels);

fig = figure('Color', 'w', 'Units', 'centimeters', 'Name', 'Primer hybridization map', ...
    'NumberTitle', 'off', 'WindowStyle', 'normal');
% PLACE THE WINDOW SO THE WHOLE OF IT IS ON SCREEN.                  
% This used to read the default figure's Position and substitute only the SIZE:
%       pos0 = get(fig, 'Position');  set(fig, 'Position', [pos0(1:2) opts.FigureSizeCm]);
% A figure Position is [left bottom width height] measured from the BOTTOM-left of the display, so
% keeping the default bottom and raising the height from about 11 cm to 23 cm grows the window
% UPWARDS. On any display where that bottom edge plus 23 cm passes the top of the screen the title
% bar and the first centimetres of the canvas end up above it, which is what put the top of the map
% off screen. The window also has to clear its own chrome, the title bar plus the menu and tool
% bars, which the Position rectangle does NOT include.
FitFigureOnScreen(fig, opts.FigureSizeCm);
% LAYOUT. The limits are set to what the artwork actually occupies rather than to a round number.
% The rings reach 1.045 with their markers, the position ticks sit at 1.085 and a region label runs
% outward from 1.235. Dropping the repeated "kcal mol^-1" from every label (it is on the colour-bar
% axis already) is what makes the horizontal budget close: the labels used to reach past 1.8 and
% forced limits of +/-1.70, which left the ring filling barely half the canvas and stranded the
% colour bars about five centimetres clear of it.
% THE FIGURE IS NOT NECESSARILY THE SIZE WE ASKED FOR.
% MATLAB clamps an on-screen figure to the display, and it does not clamp both dimensions equally.
% A request for 30 x 23 cm came back as 26.8 x 18.2 cm on one machine, turning an aspect of 1.30
% into 1.48. Every normalised position below used to be a constant tuned for the requested aspect,
% so on that machine the drawing letterboxed inside its own axes and the rings fell to 39 per cent
% of the width. The size is therefore READ BACK and the layout computed from what we actually got.
set(fig, 'Units', 'centimeters');
figPos = get(fig, 'Position');
figW = figPos(3); figH = figPos(4);
if (abs(figW - opts.FigureSizeCm(1)) > 0.5)||(abs(figH - opts.FigureSizeCm(2)) > 0.5)
    warning('PrimerHybridizationMap:FigureClamped', ...
        ['The window is %.1f x %.1f cm where %.1f x %.1f cm was requested, because the display ' ...
        'cannot hold the larger figure. The layout has been recomputed for the size actually ' ...
        'given, and the saved file is written at the requested size regardless.'], ...
        figW, figH, opts.FigureSizeCm(1), opts.FigureSizeCm(2));
end

% Bands reserved at the bottom for the caption and at the right for the colour bars, as a fraction
% of the figure. Both are expressed in CENTIMETRES first and converted, so they do not stretch or
% collapse when the figure is resized.
% No caption band. The footnote that used to sit here restated the strand convention, the colour
% meaning and the marker meaning, all of which the rings, the colour bars and the legend already
% carry. Removing it gives its whole height back to the drawing.
capCm = 0.40;                                % a bare bottom margin, not a text band
barCm = 5.4;                                 % colour bars, their labels and the two swatches
capFrac = min(0.10, capCm/figH);
barFrac = min(0.42, barCm/figW);

% The drawing occupies xRad horizontally and yRad vertically. xRad is set by how far a region label
% reaches outward from 1.235; yRad by the rings, their ticks and the top and bottom labels.
xRad = 1.70; yRad = 1.34;
availW = (1 - barFrac)*figW;  availH = (1 - capFrac)*figH;
% Choose the axes box to have EXACTLY the aspect the data needs, so axis equal has no slack to
% letterbox. Whichever of width or height binds, the other is trimmed to match.
boxAsp = (xRad/yRad);
if availW/availH > boxAsp
    useH = availH;  useW = useH*boxAsp;      % height binds
else
    useW = availW;  useH = useW/boxAsp;      % width binds
end
axL = (availW - useW)/2/figW;
axB = capFrac + (availH - useH)/2/figH;
ax = axes('Parent', fig, 'Position', [axL axB useW/figW useH/figH]);
hold(ax, 'on'); axis(ax, 'equal'); axis(ax, 'off');
xlim(ax, [-xRad xRad]); ylim(ax, [-yRad yRad]);

% Rings, bound nucleotides and outlines
DrawAnnulus(ax, opts.RingPlus(1), opts.RingPlus(2), opts.UnboundColour);
DrawAnnulus(ax, opts.RingMinus(1), opts.RingMinus(2), opts.UnboundColour);
iP = find(isfinite(bestP)); iM = find(isfinite(bestM));
DrawNucleotides(ax, iP, mapP(cIdx(bestP(iP)), :), opts.RingPlus(1), opts.RingPlus(2), L);
DrawNucleotides(ax, iM, mapM(cIdx(bestM(iM)), :), opts.RingMinus(1), opts.RingMinus(2), L);
aa = linspace(0, 2*pi, 1441);
for r = [opts.RingPlus opts.RingMinus]
    plot(ax, r*cos(aa), r*sin(aa), 'Color', [0.55 0.55 0.55], 'LineWidth', 0.4);
end

% Region markers: inward triangles outside the (+) ring, outward triangles inside the (-) ring. A strand
% with more than opts.MaxMarkers regions marks only its labelled regions, since its slivers are then
% dense enough to find without help.
[~, oP] = sort(regP.dG); [~, oM] = sort(regM.dG);
mkP = 1:numel(regP.dG); mkM = 1:numel(regM.dG);
if numel(mkP) > opts.MaxMarkers, mkP = oP(1:min(nLabels, numel(oP)))'; end
if numel(mkM) > opts.MaxMarkers, mkM = oM(1:min(nLabels, numel(oM)))'; end
for k = mkP
    DrawMarker(ax, regP.xMid(k), L, opts.RingPlus(2) + 0.005, opts.RingPlus(2) + 0.045, ...
        mapP(cIdx(regP.dG(k)), :));
end
for k = mkM
    DrawMarker(ax, regM.xMid(k), L, opts.RingMinus(1) - 0.005, opts.RingMinus(1) - 0.045, ...
        mapM(cIdx(regM.dG(k)), :));
end
% Extension-threat markers, in the gap between the two rings. A threat is deliberately drawn even
% when its nucleotides are grey: the anchored duplex it names can be weaker than the display cutoff
% and still be the only structure on the plasmid that Taq can extend from.
rThreat = [opts.RingPlus(1) - 0.007, opts.RingMinus(2) + 0.007];
for k = ordThreat'
    if S.Strand(k) == '+'
        rr = rThreat(1); cc = mapP(cIdx(S.dGanchEff(k)), :);
    else
        rr = rThreat(2); cc = mapM(cIdx(S.dGanchEff(k)), :);
    end
    aT = XToAngle(S.AnchMid(k), L);
    plot(ax, rr*cos(aT), rr*sin(aT), 'o', 'MarkerSize', 4.2, 'MarkerFaceColor', cc, ...
        'MarkerEdgeColor', [0 0 0], 'LineWidth', 0.7);
end
% (allMarked and markTxt went with the caption they described.)

% Coordinate ticks on the (+) numbering
tickStep = NiceStep(L/12);
tickPos = unique([1, tickStep:tickStep:L]);
tickPos = tickPos((tickPos == 1)|(L - tickPos >= 0.2*tickStep));
for p = tickPos
    a = XToAngle(p - 0.5, L);
    plot(ax, [1.00 1.028]*cos(a), [1.00 1.028]*sin(a), 'Color', [0.45 0.45 0.45], 'LineWidth', 1);
    text(ax, 1.085*cos(a), 1.085*sin(a), sprintf('%d', p), 'FontSize', opts.TickFontSize, ...
        'Color', [0.40 0.40 0.40], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

% Strand direction. Each label sits at 12 o'clock, where clockwise reads left to right, unless a
% displayed region occupies that arc; it then moves to 6 o'clock, where clockwise reads right to left.
DrawDirectionLabel(ax, bestP, mean(opts.RingPlus), '(+)', true, L, opts.TickFontSize);
DrawDirectionLabel(ax, bestM, mean(opts.RingMinus), minusSym, false, L, opts.TickFontSize);

% Labels for the most stable regions of each strand
% Every entry carries its own marker, so the region, threat and seed labels below can share one
% ring and one call to SpreadAngles while still reading as different classes. ONE SHAPE MEANS ONE
% THING ACROSS THE WHOLE FIGURE:
%     triangle                     a displayed region, matching the triangles on the rings
%     circle                       an extension threat, matching the circles in the gap
%     square, diamond, pentagram   an exact-complement seed
% A region label used to carry a circle, which collided with the threat circles once threats were
% labelled as well, and with the legend, where a circle is captioned as the off-target site.
lab = struct('txt', {}, 'cw', {}, 'rgb', {}, 'mk', {}, 'msz', {}, 'edge', {}, 'lw', {});
for k = oP(1:min(nLabels, numel(oP)))'
    lab(end + 1) = struct('txt', {SiteLabel('(+) ', regP.s(k), regP.e(k), ...
        regP.pf(k), regP.pt(k), regP.dG(k), regP.Tm(k))}, 'cw', 360*regP.xMid(k)/L, ...
        'rgb', mapP(cIdx(regP.dG(k)), :), 'mk', '^', 'msz', 5.0, ...
        'edge', [0.25 0.25 0.25], 'lw', 1); %#ok<SAGROW>
end
for k = oM(1:min(nLabels, numel(oM)))'
    lab(end + 1) = struct('txt', {SiteLabel([minusSym ' '], regM.s(k), regM.e(k), ...
        regM.pf(k), regM.pt(k), regM.dG(k), regM.Tm(k))}, 'cw', 360*regM.xMid(k)/L, ...
        'rgb', mapM(cIdx(regM.dG(k)), :), 'mk', '^', 'msz', 5.0, ...
        'edge', [0.25 0.25 0.25], 'lw', 1); %#ok<SAGROW>
end
% EXTENSION-THREAT LABELS. The circles in the gap say WHERE a threat is; they cannot say how much
% of the primer is engaged or which part, and that is the thing that decides whether a site can
% prime at all. These labels carry both ranges for the 3'-REACHABLE duplex, which is a different
% and usually smaller duplex than the one colouring the ring, so the numbers here will not match a
% region label at the same locus. That is the point of listing them separately.
for k = ordThreat(1:min(opts.ThreatLabels, numel(ordThreat)))'
    if S.Strand(k) == '+'
        preT = '(+) '; rgbT = mapP(cIdx(S.dGanchEff(k)), :);
    else
        preT = [minusSym ' ']; rgbT = mapM(cIdx(S.dGanchEff(k)), :);
    end
    lab(end + 1) = struct('txt', {SiteLabel(preT, S.AnchStart(k), S.AnchEnd(k), ...
        S.AnchFrom(k), S.AnchTo(k), S.dGanch(k), S.TmAnch(k))}, ...
        'cw', 360*S.AnchMid(k)/L, 'rgb', rgbT, 'mk', 'o', 'msz', 4.2, ...
        'edge', [0 0 0], 'lw', 1.1); %#ok<SAGROW>
end

% EXACT-COMPLEMENT SEED LABELS.
% These go in the same ring as the region labels rather than in a ring of their own, so that
% SpreadAngles separates a seed label from the region label of the same locus instead of letting
% the two print on top of each other. A site that is already coloured, already labelled as a
% region and already circled as a threat still gets this label: it is a separate statement about
% the locus and nothing above suppresses it.
% NO SEPARATE MARKER IS DRAWN IN THE GAP BETWEEN THE RINGS, and the reason is geometric rather
% than a preference. That gap is 0.045 data units wide, which at this scale is about nine points,
% and the two extension-threat lanes already occupy it; a third class of marker placed there would
% overlap the circles it exists to be told apart from. The leader line carries the position
% instead, since its inner end sits at the seed's own angular midpoint, and the dot at its outer
% end carries the class:
%     square      a 3'-terminal run only
%     diamond     an internal run only
%     pentagram   both, which is what any 3' run of Seed3Len + 1 or longer produces
% The dot takes the colour of the register's whole-site dG, the same quantity the ring is painted
% with, so a seed on a grey arc shows as the palest colour on the scale and is still legible
% against the black edge.
primeCh = char(8242);
for k = ordSeed'
    if S.Strand(k) == '+'
        preS = '(+) '; rgbS = mapP(cIdx(S.dG(k)), :);
    else
        preS = [minusSym ' ']; rgbS = mapM(cIdx(S.dG(k)), :);
    end
    if isSeed3(k)&&isSeedInt(k)
        mkS = 'p'; mszS = 8.5; cwS = 360*S.Seed3Mid(k)/L;
        t1S = [preS RangeTxt(S.Seed3Start(k), S.Seed3End(k))];
        t2S = sprintf('3%s run %d, int %d', primeCh, S.End3Run(k), S.IntRun(k));
    elseif isSeed3(k)
        mkS = 's'; mszS = 5.4; cwS = 360*S.Seed3Mid(k)/L;
        t1S = [preS RangeTxt(S.Seed3Start(k), S.Seed3End(k))];
        t2S = sprintf('3%s run %d nt', primeCh, S.End3Run(k));
    else
        mkS = 'd'; mszS = 5.8; cwS = 360*S.IntMid(k)/L;
        t1S = [preS RangeTxt(S.IntStart(k), S.IntEnd(k))];
        t2S = sprintf('int run %d nt', S.IntRun(k));
    end
    lab(end + 1) = struct('txt', {{t1S, t2S}}, 'cw', cwS, 'rgb', rgbS, ...
        'mk', mkS, 'msz', mszS, 'edge', [0 0 0], 'lw', 1.1); %#ok<SAGROW>
end
% The ring holds 360/LabelMinSepDeg labels before SpreadAngles has to push one past its neighbour
% and the order stops being readable. Say so rather than silently drawing a smear.
if numel(lab) > floor(360/max(opts.LabelMinSepDeg, 1E-6))
    warning('PrimerHybridizationMap:LabelCrowding', ...
        ['%d labels are being placed around a ring that holds %d at %g degrees apart. Lower the ' ...
        'number of region labels or opts.SeedMax, or reduce opts.LabelMinSepDeg.'], ...
        numel(lab), floor(360/opts.LabelMinSepDeg), opts.LabelMinSepDeg);
end
% HOW FAR APART LABELS HAVE TO SIT IS A CONSEQUENCE, NOT A CONSTANT. Two labels collide when the
% arc between them is shorter than a label is tall, so the separation they need follows from the
% line count, the font size and the scale of the axes. opts.LabelMinSepDeg was chosen when a label
% was two lines. It became three when the primer register and the melting temperature were added,
% and labels on the western flank began to overprint. The option is now the FLOOR; the working
% value is computed, and opts.LabelMinSepDeg only raises it.
if isempty(lab)
    labLines = 1;                            % nothing to space; keep the arithmetic defined
else
    labLines = max(cellfun(@numel, {lab.txt}));
end
minSep = max(opts.LabelMinSepDeg, ...
    LabelSepDeg(ax, labLines, opts.LabelFontSize, 1.235, opts.LabelLineFactor));

hLab = gobjects(0); hLead = gobjects(0); hMark = gobjects(0); kept = [];
if ~isempty(lab)
    cwTrue = [lab.cw];
    cwLab = SpreadAngles(cwTrue, minSep);
    for k = 1:numel(lab)
        % text() refuses a cell holding an empty character vector and takes the figure down with
        % it, after every expensive step has already completed. One filter here means no future
        % label builder can cost a whole run over a blank line.
        txt = lab(k).txt;
        txt = txt(~cellfun(@isempty, txt));
        if isempty(txt), continue; end
        a0 = pi/2 - cwTrue(k)*pi/180; a1 = pi/2 - cwLab(k)*pi/180;
        hLead(end + 1) = plot(ax, ...
            [1.05*cos(a0), 1.12*cos(a0), 1.18*cos(a1)], ...
            [1.05*sin(a0), 1.12*sin(a0), 1.18*sin(a1)], ...
            'Color', [0.30 0.30 0.30], 'LineWidth', 1); %#ok<SAGROW>
        hMark(end + 1) = plot(ax, 1.18*cos(a1), 1.18*sin(a1), lab(k).mk, ...
            'MarkerSize', lab(k).msz, 'MarkerFaceColor', lab(k).rgb, ...
            'MarkerEdgeColor', lab(k).edge, 'LineWidth', lab(k).lw); %#ok<SAGROW>
        [ha, va] = LabelAlignment(a1);
        hLab(end + 1) = text(ax, 1.235*cos(a1), 1.235*sin(a1), txt, ...
            'FontSize', opts.LabelFontSize, 'HorizontalAlignment', ha, ...
            'VerticalAlignment', va, 'Interpreter', 'tex'); %#ok<SAGROW>
        kept(end + 1) = k; %#ok<SAGROW>
    end
end

% THE LABEL RING IS SETTLED BY MEASUREMENT, NOT BY A MODEL OF IT. Two things have to come right
% and they fight each other: no label may run past the axis limit, and no two labels may overprint.
%
% WIDTH. xRad is chosen before any text exists, because the axes box is sized from it, so it can
% only ever be a guess. Adding the primer register to every label broke that guess and clipped the
% westmost one. The drawn labels are therefore asked for their extent and the limits open if any
% of them overruns. axis equal holds the data aspect at 1:1, so the rings stay circular and simply
% occupy less of the box.
%
% OVERLAP. The arc needed to clear a label depends on WHERE ON THE RING it sits. Two labels due
% west stack fully vertically, so the arc between them is the clearance. At a polar angle a the
% same arc buys only arc*abs(cos a) of vertical separation, so labels off the east-west line need
% more of it: the pair at 6088 and 5917 sits near 147 and 156 degrees, where abs(cos a) is about
% 0.87, and an arc computed as though it were 1.0 left them touching. LabelSepDeg models the
% east-west case and is a good starting point, not an answer.
%
% Rather than carry the geometry any further, the drawn boxes are asked whether they actually
% overlap and the separation is opened until they do not. Opening it moves labels, which changes
% both the widths and the overlaps, so the two tests share one loop and it runs until both are
% satisfied. The ring is finite, so the loop also stops when the labels have filled it, which is
% the condition the crowding warning above already describes.
if ~isempty(hLab)
    for pass = 1:10
        drawnow;
        ext = LabelExtents(hLab);
        reach = max([abs(ext(:, 1)); abs(ext(:, 1) + ext(:, 3))]);
        if reach > xRad
            xRad = reach*1.02;                   % a two per cent hair so nothing touches the edge
            xlim(ax, [-xRad xRad]);
            drawnow;
            ext = LabelExtents(hLab);
        end
        if ~AnyOverlap(ext), break; end
        if minSep*numel(hLab) >= 345, break; end  % the ring is full; nothing more can be gained
        % Start from whichever is larger: a fifteen per cent opening, or the east-west requirement
        % recomputed at the scale the limits now have.
        minSep = max(minSep*1.15, ...
            LabelSepDeg(ax, labLines, opts.LabelFontSize, 1.235, opts.LabelLineFactor));
        cwLab = SpreadAngles(cwTrue, minSep);
        for q = 1:numel(kept)
            k = kept(q);
            a0 = pi/2 - cwTrue(k)*pi/180; a1 = pi/2 - cwLab(k)*pi/180;
            set(hLead(q), 'XData', [1.05*cos(a0), 1.12*cos(a0), 1.18*cos(a1)], ...
                'YData', [1.05*sin(a0), 1.12*sin(a0), 1.18*sin(a1)]);
            set(hMark(q), 'XData', 1.18*cos(a1), 'YData', 1.18*sin(a1));
            [ha, va] = LabelAlignment(a1);
            set(hLab(q), 'Position', [1.235*cos(a1), 1.235*sin(a1), 0], ...
                'HorizontalAlignment', ha, 'VerticalAlignment', va);
        end
    end
    drawnow;
    if AnyOverlap(LabelExtents(hLab))
        warning('PrimerHybridizationMap:LabelsStillOverlap', ...
            ['Labels still overlap after opening their spacing to %.1f degrees, which is as far ' ...
            'as %d of them will go on one ring. Lower the number of region labels, lower ' ...
            'opts.ThreatLabels or opts.SeedMax, or reduce opts.LabelFontSize.'], ...
            minSep, numel(hLab));
    end
end

% Centre text. Each line is its own text object so that wide lines can sit near the equator of the
% circle and narrow ones near its poles; the font shrinks until every line fits inside TextRadius.
drawnow;
centreLines = {['\bf' TexEscape(dnaName) '\rm'], sprintf('%d nt, Circular %s', L, topoTxt), '', ...
    sprintf('\\bfPrimer\\rm (%d nt)', n)};
chunkLen = ceil(n/ceil(n/32));                  % primers above 32 nt split into near-equal lines
for c0 = 1:chunkLen:n
    chunk = primer(c0:min(c0 + chunkLen - 1, n));
    if c0 == 1, chunk = ['5' char(8242) '-' chunk]; end
    if c0 + chunkLen - 1 >= n, chunk = [chunk '-3' char(8242)]; end
    centreLines{end + 1} = ['\fontname{Courier}' '\bf' chunk '\rm' '\fontname{Helvetica}'];
end
centreLines = [centreLines, {'', ...
    ['Annealing Temperature: ' sprintf('%.1f', Ta) ' ^{\circ}C'], ...
    ['Monovalent Cation: ' Num(monv_mM, 4) ' mM'], ...
    ['Bivalent Cation: ' Num(biv_mM, 4) ' mM (' Num(freeBiv*1E3, 3) ' mM free)'], ...
    ['Total dNTP: ' Num(dNTP_uM, 4) ' \muM'], ...
    ['Primer: ' Num(p0_nM, 4) ' nM'], ...
    ['Template: ' Num(cT_nM, 3) ' nM (' Num(ngul, 3) ' ng \muL^{-1})'], ...
    '', '\bfStrongest Binding\rm'}];
if ~isempty(iBest)
    if S.Strand(iBest) == '+', sb = '(+)'; else, sb = minusSym; end
    centreLines = [centreLines, { ...
        sprintf('%s Strand, %s, %d/%d Bases Paired', sb, RangeTxt(S.StartPos(iBest), S.EndPos(iBest)), ...
        sum(S.Pairing{iBest} == '|'), S.Len(iBest)), ...
        ['\DeltaG_{Hyb} = ' FmtSigned(S.dG(iBest), 2) ' kcal mol^{-1}'], ...
        ['T_{m} = ' Num(S.Tm(iBest), 4) ' ^{\circ}C at ' Num(fFree*1E9, 5) ' nM free primer']}];
else
    centreLines = [centreLines, {'No Register Passed the Survey'}];
end
centreLines = [centreLines, {'', ['Shown: \theta \geq ' Num(thetaPct, 3) '%  (\DeltaG_{Hyb} \leq ' ...
    FmtSigned(dGcut, 2) ' kcal mol^{-1})'], ...
    ['\epsilon Data: ' TexSafe(effCite)]}];
% The template name is the one thing a reader looks for first, so it is set larger than the body.
% The factor rides through the shrink loop, so it stays proportional however much the block has to
% give up to fit inside TextRadius.
centreScale = ones(numel(centreLines), 1);
centreScale(1) = opts.TitleScale;
LayoutCentreText(ax, centreLines, opts.CentreFontSize, opts.TextRadius, 'Helvetica', centreScale, ...
    opts.CentreFontMax);

% Colour bars on a shared DeltaG axis
% The bars sit in the band reserved for them, positioned from the axes' real right edge rather
% than from a constant, so they stay beside the rings at any figure size.
barL   = axL + useW/figW + 0.018;
barW   = 0.030;
barGap = 0.034;
barH   = 0.40*useH/figH;
barB   = axB + 0.34*useH/figH;
% THE SWATCH STACK HANGS BELOW THE BARS AND swDrop IS ITS REACH. Every swatch and caption below is
% positioned as barB minus a fixed offset, and swDrop is the largest of those offsets, the foot of
% the lowest caption. If the stack would reach into the caption band the whole thing is lifted
% until it clears, which is cheaper than letting the two overlap.
%
% KEEP swDrop IN STEP WITH THE OFFSETS BELOW. It was once written as its own number, 0.336, and a
% fourth swatch was added without raising it. The lift then cleared only the old stack, the last
% caption landed at a negative position, and annotation rejects any position outside 0 to 1, so
% the figure died at the final step. One constant, used by the lift and by the last caption.
% WHICH SWATCHES THE LEGEND CARRIES. Only the classes actually drawn on this map. A legend that
% names a shape the figure does not contain sends a reader hunting for something that is not
% there, which is what a square captioned "Exact-Complement Seed" did on a primer with no seeds.
% The stack is therefore built here, before anything is positioned, and its height follows from
% how many entries survive.
swOff0  = 0.070;                             % first swatch, below barB
swPitch = 0.098;                             % slot to slot
swGap2  = 0.062; swCapH2 = 0.048; swCapW2 = 0.200;   % two-line caption
swGap1  = 0.050; swCapH1 = 0.030; swCapW1 = 0.170;   % one-line caption
swatch = struct('kind', {}, 'face', {}, 'edge', {}, 'lw', {}, 'cap', {});
j = 0;
% The grey of an unbound nucleotide is always on the rings, so this one is unconditional.
j = j + 1;
swatch(j).kind = 'wide'; swatch(j).face = opts.UnboundColour;
swatch(j).edge = [0.55 0.55 0.55]; swatch(j).lw = 0.4;
swatch(j).cap  = {'No Displayed', 'Duplex'};
if ~isempty(regP.s)||~isempty(regM.s)
    j = j + 1;
    swatch(j).kind = 'tri'; swatch(j).face = [0.65 0.65 0.65];
    swatch(j).edge = [0.25 0.25 0.25]; swatch(j).lw = 0.7;
    swatch(j).cap  = {'Displayed Region'};
end
if ~isempty(ordThreat)
    j = j + 1;
    swatch(j).kind = 'circle'; swatch(j).face = [0.65 0.65 0.65];
    swatch(j).edge = [0 0 0]; swatch(j).lw = 0.7;
    swatch(j).cap  = {'Extendable', 'Off-Target Site'};
end
if ~isempty(ordSeed)
    j = j + 1;
    swatch(j).kind = 'square'; swatch(j).face = [0.65 0.65 0.65];
    swatch(j).edge = [0 0 0]; swatch(j).lw = 0.7;
    swatch(j).cap  = {'Exact-Complement', 'Seed'};
end
% swDrop is the foot of the LAST caption, and the lift below uses it. It is derived from the same
% constants the drawing loop uses, so the two cannot drift: that drift is what drove a caption to
% a negative position once already, and annotation rejects anything outside 0 to 1.
if numel(swatch(end).cap) > 1, swLastGap = swGap2; else, swLastGap = swGap1; end
swDrop = swOff0 + (numel(swatch) - 1)*swPitch + swLastGap;

if barB - swDrop < capFrac + 0.01
    barB = capFrac + swDrop + 0.010;
end
if barB - swDrop < 0
    % Unreachable while the lift above works, and worth stating anyway: annotation's own complaint
    % names only a position, not the reason, and it arrives after the whole figure has been drawn.
    error('PrimerHybridizationMap:SwatchStackOverruns', ...
        ['The legend swatch stack reaches %.3f below the colour bars, which puts its last ' ...
        'caption at %.3f of the figure. annotation accepts 0 to 1 only. Raise swDrop to the ' ...
        'largest offset used in the swatch block, or shorten the stack.'], swDrop, barB - swDrop);
end
cbP = axes('Parent', fig, 'Position', [barL barB barW barH]);
image(cbP, 1, linspace(gLo, gHi, opts.ColourLevels), reshape(mapP, [opts.ColourLevels 1 3]));
set(cbP, 'YDir', 'normal', 'XTick', [], 'YTickLabel', {}, 'TickDir', 'out', ...
    'FontSize', opts.LegendFontSize, 'Box', 'on');
ylim(cbP, [gLo gHi]);
title(cbP, '(+)', 'FontWeight', 'normal', 'FontSize', opts.BarFontSize);
cbM = axes('Parent', fig, 'Position', [barL + barGap, barB, barW, barH]);
image(cbM, 1, linspace(gLo, gHi, opts.ColourLevels), reshape(mapM, [opts.ColourLevels 1 3]));
set(cbM, 'YDir', 'normal', 'XTick', [], 'YAxisLocation', 'right', 'TickDir', 'out', ...
    'FontSize', opts.LegendFontSize, 'Box', 'on');
ylim(cbM, [gLo gHi]);
title(cbM, minusSym, 'FontWeight', 'normal', 'FontSize', opts.BarFontSize);
ylabel(cbM, '\DeltaG_{hyb} (kcal mol^{-1})', 'FontSize', opts.BarFontSize, 'Interpreter', 'tex', 'FontWeight', 'bold');
% Swatches tucked immediately under the bars, again measured from the bars rather than fixed.
% One loop over the spec built above, so a class that is absent from the map leaves no gap and the
% stack simply ends sooner.
%
% THE TRIANGLE IS A PLOTTED MARKER, NOT AN ANNOTATION. annotation accepts line, arrow,
% doublearrow, textarrow, textbox, ellipse and rectangle, and nothing else, so there is no
% triangular counterpart to the ellipse and rectangle swatches. A one-marker axes draws the real
% marker instead, which also makes the swatch identical to the thing it stands for.
%
% The square swatch is a rectangle of the same box as the ellipse, which is a square on screen
% because the box compensates the figure aspect. It stands for all three seed shapes; which shape
% means what is on the console and in the header, since three captions here would cost more height
% than the distinction is worth.
swC = barL + (barGap + barW)/2;              % centre of the pair of bars
swW = barGap + barW;
for k = 1:numel(swatch)
    yk = barB - (swOff0 + (k - 1)*swPitch);
    switch swatch(k).kind
        case 'wide'
            annotation(fig, 'rectangle', [barL, yk, swW, 0.022], ...
                'FaceColor', swatch(k).face, 'Color', swatch(k).edge, 'LineWidth', swatch(k).lw);
        case 'circle'
            annotation(fig, 'ellipse', [swC - 0.0055, yk, 0.0110, 0.0145], ...
                'FaceColor', swatch(k).face, 'Color', swatch(k).edge, 'LineWidth', swatch(k).lw);
        case 'square'
            annotation(fig, 'rectangle', [swC - 0.0055, yk, 0.0110, 0.0145], ...
                'FaceColor', swatch(k).face, 'Color', swatch(k).edge, 'LineWidth', swatch(k).lw);
        case 'tri'
            axTri = axes('Parent', fig, 'Position', [swC - 0.0090, yk, 0.0180, 0.0180], ...
                'Visible', 'off', 'XLim', [-1 1], 'YLim', [-1 1]);
            hold(axTri, 'on');
            plot(axTri, 0, 0, '^', 'MarkerSize', 7, 'MarkerFaceColor', swatch(k).face, ...
                'MarkerEdgeColor', swatch(k).edge, 'LineWidth', swatch(k).lw);
    end
    if numel(swatch(k).cap) > 1
        capGap = swGap2; capH = swCapH2; capW = swCapW2;
    else
        capGap = swGap1; capH = swCapH1; capW = swCapW1;
    end
    annotation(fig, 'textbox', ...
        [swC - capW/2, barB - (swOff0 + (k - 1)*swPitch + capGap), capW, capH], ...
        'String', swatch(k).cap, 'EdgeColor', 'none', 'FontSize', opts.LegendFontSize, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

%% 9. Command Window summary
fprintf('\nPrimer 5''-%s-3'' (%d nt) on %s (%d nt, circular %s)\n', primer, n, dnaName, L, topoTxt);
fprintf('Registers scored: %d on (+), %d on (-); %d distinct duplexes evaluated with dGMismatch.\n', ...
    numel(sitesP.dG), numel(sitesM.dG), cache.Count);
fprintf('Free primer at equilibrium: %.5g nM of %.5g nM (%.3g%% bound).\n', fFree*1E9, p0_nM, 100*boundFrac);
fprintf('Display cutoff: occupancy >= %.3g%%, i.e. dG <= %.2f kcal/mol at %.1f C.\n', thetaPct, dGcut, Ta);
fprintf('Melting temperatures are quoted at the free primer, %.5g nM, where occupancy is one half.\n', fFree*1E9);
fprintf('Extension efficiency table: %s.  Proofreading class: %s.\n', effSource, proofClass);
fprintf('Displayed regions: %d on (+), %d on (-).\n', numel(regP.dG), numel(regM.dG));
fprintf('Extension threats (3''-reachable, endTol %d%s, cognate excluded): %d at theta*eps >= %.3g%%', ...
    opts.EndTol, Tern(proofread, ' = whole primer, proofreading', ''), nnz(isThreat), 100*threatMin);
if nnz(isThreat) > numel(ordThreat), fprintf(', %d listed', numel(ordThreat)); end
fprintf('.\n');
for strandSym = '+-'
    j = find((S.Strand == strandSym)&S.Anchored&~isCognate);
    if isempty(j)
        fprintf('  (%c) no extendable off-target register at all.\n', strandSym);
    else
        [gb, jb] = min(S.dGanch(j)); kb = j(jb);
        fprintf(['  (%c) most stable anchored off-target : dG %7.2f  Tm %7.2f  theta %9.3g  %d-%d, ' ...
            'primer %d-%d, eps from %s, eps %.3g -> theta*eps %.3g%s\n'], strandSym, gb, ...
            S.TmAnch(kb), thetaAnch(kb), S.AnchStart(kb), S.AnchEnd(kb), S.AnchFrom(kb), S.AnchTo(kb), ...
            S.EpsPair{kb}, S.Eps(kb), primeRate(kb), ...
            Tern(S.TerminalPaired(kb), '', ' (3'' base mispaired)'));
        [~, jr] = max(primeRate(j)); kr = j(jr);
        if kr ~= kb
            fprintf(['  (%c) fastest priming off-target     : dG %7.2f  Tm %7.2f  theta %9.3g  %d-%d, ' ...
                'primer %d-%d, eps from %s, eps %.3g -> theta*eps %.3g\n'], strandSym, ...
                S.dGanch(kr), S.TmAnch(kr), thetaAnch(kr), S.AnchStart(kr), S.AnchEnd(kr), S.AnchFrom(kr), ...
                S.AnchTo(kr), S.EpsPair{kr}, S.Eps(kr), primeRate(kr));
        end
    end
end
fprintf('\n');
fprintf('MOST STABLE DISPLAYED SITES  (whole duplex; may leave the primer 3'' end unpaired)\n');
iShown = find(shown);
[~, ord] = sort(S.dG(iShown));
iShown = iShown(ord);
if isempty(iShown), fprintf('  none above the display threshold.\n'); end
for k = iShown(1:min(10, numel(iShown)))'
    fprintf('  (%c) %5d-%-5d  dG %7.2f kcal/mol  Tm %7.2f C  theta %.4g   3''-end paired: %d\n', S.Strand(k), ...
        S.StartPos(k), S.EndPos(k), S.dG(k), S.Tm(k), theta(k), S.ThreePrimePaired(k));
    fprintf('        5''-%s-3''  primer\n           %s\n        3''-%s-5''  template\n', ...
        S.Segment{k}, S.Pairing{k}, S.Opposing{k});
end

fprintf(['\nEXTENSION THREATS  (3''-reachable duplex, the sites the polymerase can extend from)\n' ...
    'Ranked by theta*epsilon, the relative rate of productive priming. dG is the reachable duplex,\n' ...
    'not the whole-site dG above; dGeff adds the kcal/mol that the terminal mispair costs in rate.\n']);
if proofread
    fprintf(['  Proofreading was declared, so endTol is the whole primer: the exonuclease trims any\n' ...
        '  3'' overhang back to the outermost paired base and primes from there. Every epsilon is 1,\n' ...
        '  so this is an occupancy ranking over every register with a pair.\n']);
elseif ~opts.UseExtEff
    fprintf('  opts.UseExtEff is FALSE, so every epsilon is 1 and this is an occupancy ranking.\n');
end
if isempty(ordThreat)
    fprintf('  none at theta*eps >= %.3g%%.\n', 100*threatMin);
end
for k = ordThreat'
    fprintf(['  (%c) %5d-%-5d  dG %7.2f  Tm %7.2f  dGeff %7.2f  theta %9.3g  eps %8.3g  theta*eps %9.3g\n' ...
        '        primer %d-%d of %d, eps from %s%s%s\n'], ...
        S.Strand(k), S.AnchStart(k), S.AnchEnd(k), S.dGanch(k), S.TmAnch(k), S.dGanchEff(k), thetaAnch(k), ...
        S.Eps(k), primeRate(k), S.AnchFrom(k), S.AnchTo(k), n, S.EpsPair{k}, ...
        Tern(S.TerminalPaired(k), ' (3'' base paired)', ' (3'' base MISPAIRED)'), ...
        Tern(theta(k) >= thetaMin, '', ', below the display cutoff'));
    fprintf('        5''-%s-3''  primer\n           %s\n        3''-%s-5''  template\n', ...
        S.AnchSegment{k}, S.AnchPairing{k}, S.AnchOpposing{k});
end
fprintf(['\nEXACT-COMPLEMENT SEEDS  (contiguous perfect runs; combinatorial, not thermodynamic)\n' ...
    'A 3'' run of >= %d nt reaches the primer terminus, so it nucleates AND leaves a paired 3''-OH\n' ...
    'the polymerase can extend from at once. An internal run of >= %d nt lies wholly inside primer\n' ...
    'positions 1-%d, so it cannot prime as it stands; it holds primer at the locus, and a\n' ...
    'proofreading exonuclease converts it into a primable site by trimming the 3'' overhang.\n' ...
    'No free energy is quoted for a run on its own: dGMismatch charges a helix initiation at each\n' ...
    'end, so a 7-mer scored alone comes back positive for almost any sequence and would read as\n' ...
    'reassuring. The dG and dGanch of the register carrying the run are the meaningful numbers.\n'], ...
    opts.Seed3Len, opts.SeedIntLen, n - 1);
fprintf('Seed sites: %d (%d with a 3'' run, %d with an internal run)', ...
    nnz(isSeedAny), nnz(isSeed3), nnz(isSeedInt));
if nSeedAll > numel(ordSeed), fprintf(', %d labelled', numel(ordSeed)); end
if any(isCognate)
    fprintf('. The cognate site satisfies both by construction and is excluded.\n');
else
    fprintf('. There is no cognate site on this template.\n');
end
if isempty(ordSeed)
    fprintf('  none.\n');
end
for k = ordSeed'
    if isSeed3(k)&&isSeedInt(k)
        cls = 'both      ';
    elseif isSeed3(k)
        cls = '3'' run    ';
    else
        cls = 'internal  ';
    end
    fprintf('  (%c) %s', S.Strand(k), cls);
    if isSeed3(k)
        fprintf('%5d-%-5d  %2d nt, primer %d-%d', S.Seed3Start(k), S.Seed3End(k), ...
            S.End3Run(k), n - S.End3Run(k) + 1, n);
    else
        fprintf('%5d-%-5d  %2d nt, primer %d-%d', S.IntStart(k), S.IntEnd(k), ...
            S.IntRun(k), S.IntFrom(k), S.IntTo(k));
    end
    fprintf('  site dG %7.2f  Tm %7.2f', S.dG(k), S.Tm(k));
    if S.Anchored(k)
        fprintf('  dGanch %7.2f  theta*eps %9.3g', S.dGanch(k), primeRate(k));
    else
        fprintf('  not 3''-reachable');
    end
    fprintf('%s%s\n', Tern(shown(k), ', displayed', ''), Tern(isThreat(k), ', threat', ''));
    if isSeed3(k)&&isSeedInt(k)
        fprintf('        internal run %d nt at primer %d-%d, %d-%d\n', S.IntRun(k), ...
            S.IntFrom(k), S.IntTo(k), S.IntStart(k), S.IntEnd(k));
    end
    fprintf('        5''-%s-3''  primer\n           %s\n        3''-%s-5''  template\n', ...
        S.Segment{k}, S.Pairing{k}, S.Opposing{k});
end

kc = find(isCognate);
if ~isempty(kc)
    fprintf('\nCOGNATE SITE(S)  (every primer position paired)\n');
    for k = kc'
        fprintf('  (%c) %5d-%-5d  dG %7.2f kcal/mol  Tm %7.2f C  theta %.4g\n', S.Strand(k), ...
            S.StartPos(k), S.EndPos(k), S.dG(k), S.Tm(k), theta(k));
    end
end
fprintf('\n');

%% 10. Export
choice = questdlg('Save the map (SVG and PNG) and the site table (XLSX)?', 'Export', 'Save', 'Skip', 'Save');
if strcmp(choice, 'Save')
    defName = regexprep([dnaName '_' primer '_HybMap'], '[^\w\-]', '_');
    [sf, sp] = uiputfile({'*.svg', 'Scalable Vector Graphics (*.svg)'}, 'Save the hybridization map', ...
        [defName '.svg']);
    if ~isequal(sf, 0)
        [~, baseName] = fileparts(sf);
        outStem = fullfile(sp, baseName);
        % exportgraphics arrived in R2020a and accepts .svg only in recent releases, so SVG falls back to
        % print whenever exportgraphics is absent or refuses the extension. The fallback carries NO
        % renderer flag: -vector is the R2021b replacement for -painters and does not exist in the
        % releases where exportgraphics is missing, so passing it threw from inside this catch block
        % and killed the export after the file name had already been chosen. -dsvg is vector output
        % by construction, so no flag is needed.
        % THE SAVED FILE MUST NOT DEPEND ON THE SCREEN.
        % exportgraphics and print both render at the figure's own size, and the window may have
        % been clamped to the display. The paper is therefore set to the figure's ASPECT, scaled up
        % to at least the requested width, so the file is large and proportioned exactly like what
        % is on screen. Scaling the aspect rather than forcing the requested box is deliberate:
        % forcing it would reintroduce the letterboxing the layout above exists to avoid.
        outW = max(figW, opts.FigureSizeCm(1));
        outH = outW*figH/figW;
        set(fig, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual', ...
            'PaperPosition', [0 0 outW outH], 'PaperSize', [outW outH]);
        hasEG = ~isempty(which('exportgraphics'));
        try
            if ~hasEG, error('PrimerHybridizationMap:NoExportgraphics', 'Using print.'); end
            exportgraphics(fig, [outStem '.svg'], 'ContentType', 'vector');
        catch
            print(fig, [outStem '.svg'], '-dsvg');
        end
        % print is used for the raster too, because it honours PaperPosition while exportgraphics
        % renders at the on-screen size and would hand back the clamped figure.
        print(fig, [outStem '.png'], '-dpng', '-r1200');
        strandCol = arrayfun(@(c) c, S.Strand, 'UniformOutput', false);
        Tsites = table(strandCol, S.Register, S.PrimerFrom, S.PrimerTo, S.StartPos, S.EndPos, S.Len, ...
            S.dG, S.Tm, theta, shown, S.ThreePrimePaired, S.Segment, S.Pairing, S.Opposing, ...
            S.Anchored, S.dGanch, S.TmAnch, thetaAnch, S.Eps, S.EpsPair, S.dGanchEff, primeRate, ...
            S.AnchStart, S.AnchEnd, S.AnchFrom, S.AnchTo, ...
            S.TerminalPaired, isThreat, S.Perfect, S.AnchSegment, S.AnchPairing, S.AnchOpposing, ...
            S.End3Run, S.IntRun, isSeed3, isSeedInt, S.Seed3Start, S.Seed3End, ...
            S.IntFrom, S.IntTo, S.IntStart, S.IntEnd, ...
            'VariableNames', {'Strand', 'Register', 'PrimerFrom', 'PrimerTo', 'Start', 'End', 'Length', ...
            'dG_kcal_per_mol', 'Tm_site_C', 'Occupancy', 'Displayed', 'PrimerThreePrimePaired', ...
            'PrimerSegment_5to3', ...
            'Pairing', 'Template_3to5', ...
            'Anchored', 'dGanch_kcal_per_mol', 'Tm_anchored_C', 'OccupancyAnchored', 'Epsilon', ...
            'EpsilonPair', ...
            'dGanchEff_kcal_per_mol', 'PrimingRate', 'AnchStart', 'AnchEnd', ...
            'AnchPrimerFrom', 'AnchPrimerTo', 'ThreePrimeBasePaired', 'ExtensionThreat', 'Cognate', ...
            'AnchSegment_5to3', 'AnchPairing', 'AnchTemplate_3to5', ...
            'ThreePrimeRun', 'InternalRun', 'Seed3Prime', 'SeedInternal', ...
            'Seed3Start', 'Seed3End', 'IntPrimerFrom', 'IntPrimerTo', 'IntStart', 'IntEnd'});
        Tsites = sortrows(Tsites, 'dG_kcal_per_mol');
        Tthreat = sortrows(Tsites(Tsites.ExtensionThreat, :), 'PrimingRate', 'descend');
        Tseed = sortrows(Tsites(Tsites.Seed3Prime|Tsites.SeedInternal, :), ...
            {'Seed3Prime', 'ThreePrimeRun', 'InternalRun'}, {'descend', 'descend', 'descend'});
        Treg = table([repmat({'+'}, numel(regP.dG), 1); repmat({'-'}, numel(regM.dG), 1)], ...
            [regP.s; regM.s], [regP.e; regM.e], [regP.len; regM.len], [regP.dG; regM.dG], ...
            'VariableNames', {'Strand', 'Start', 'End', 'Length', 'dGmin_kcal_per_mol'});
        Tcond = table({'Circular DNA'; 'Length (nt)'; 'Template strands'; 'Primer (5''->3'')'; ...
            'Primer length (nt)'; 'Annealing temperature (C)'; 'Monovalent cation (mM)'; ...
            'Bivalent cation (mM)'; 'Free bivalent cation (mM)'; 'Total dNTP (uM)'; 'Primer (nM)'; ...
            'Template (nM)'; 'Template (ng/uL)'; 'Template molar mass (g/mol)'; 'Free primer (nM)'; ...
            'Occupancy threshold (%)'; 'dG cutoff (kcal/mol)'; 'Survey run (bp)'; ...
            'End tolerance (nt)'; 'Threat threshold (%)'; 'Extension threats'; ...
            'Extension efficiency applied'; 'Extension efficiency table'; ...
            'Proofreading polymerase'; ...
            '3'' seed length (nt)'; 'Internal seed length (nt)'; 'Seed sites (3'')'; ...
            'Seed sites (internal)'; 'Seed sites (any)'}, ...
            cellfun(@ValueText, {dnaName; L; topoTxt; primer; n; Ta; monv_mM; biv_mM; freeBiv*1E3; ...
            dNTP_uM; p0_nM; cT_nM; ngul; MW; fFree*1E9; thetaPct; dGcut; opts.SurveyRun; ...
            opts.EndTol; 100*threatMin; nnz(isThreat); Tern(opts.UseExtEff, 'yes', 'no'); ...
            effSource; proofClass; opts.Seed3Len; opts.SeedIntLen; nnz(isSeed3); ...
            nnz(isSeedInt); nnz(isSeedAny)}, ...
            'UniformOutput', false), ...
            'VariableNames', {'Parameter', 'Value'});
        try
            writetable(Tsites, [outStem '_Sites.xlsx'], 'Sheet', 'Sites');
            writetable(Tthreat, [outStem '_Sites.xlsx'], 'Sheet', 'ExtensionThreats');
            writetable(Tseed, [outStem '_Sites.xlsx'], 'Sheet', 'Seeds');
            writetable(Treg, [outStem '_Sites.xlsx'], 'Sheet', 'Regions');
            writetable(Tcond, [outStem '_Sites.xlsx'], 'Sheet', 'Conditions');
        catch xlsErr
            fprintf('XLSX export failed (%s); writing CSV files instead.\n', xlsErr.message);
            writetable(Tsites, [outStem '_Sites.csv']);
            writetable(Tthreat, [outStem '_ExtensionThreats.csv']);
            writetable(Tseed, [outStem '_Seeds.csv']);
            writetable(Treg, [outStem '_Regions.csv']);
            writetable(Tcond, [outStem '_Conditions.csv']);
        end
        fprintf('Saved %s.svg, %s.png and the site table.\n', outStem, outStem);
    end
end

%% ____________________________________________________________________________________________________________
%% Local functions

function FitFigureOnScreen(fig, wantCm)
%% FITFIGUREONSCREEN  Size the figure to the request, then put all of it on the display 
% Three things have to be true at once and the old one-liner guaranteed none of them.
%   The window must not be TALLER than the usable screen. MATLAB clamps a figure that is too big,
%   but it clamps the SIZE and leaves the origin alone, so a clamped window can still sit too high.
%   The window must clear its own chrome. A figure Position is the CANVAS, not the frame: the title
%   bar, the menu bar and the tool bar all sit above the rectangle and none of them are in it.
%   Several centimetres of the request therefore have to be given back before placing it.
%   It should land on the monitor the figure is already on, not jump to the primary one, which is
%   why MonitorPositions is consulted rather than ScreenSize alone.
CHROME_CM = 3.0;    % title bar plus menu bar plus tool bar, the generous end of what these cost
TASKBAR_CM = 1.5;   % a docked taskbar or dock along the bottom edge
MARGIN_CM = 0.5;    % breathing room at the sides

oldRootUnits = get(groot, 'Units');
set(groot, 'Units', 'centimeters');
cleanupRoot = onCleanup(@() set(groot, 'Units', oldRootUnits));
try
    mon = get(groot, 'MonitorPositions');
catch
    mon = [];
end
if isempty(mon)||(size(mon, 2) ~= 4)
    mon = get(groot, 'ScreenSize');
end

set(fig, 'Units', 'centimeters');
here = get(fig, 'Position');
cx = here(1) + here(3)/2; cy = here(2) + here(4)/2;
pick = 1;
for k = 1:size(mon, 1)
    if (cx >= mon(k,1))&&(cx <= mon(k,1) + mon(k,3))&& ...
       (cy >= mon(k,2))&&(cy <= mon(k,2) + mon(k,4))
        pick = k; break
    end
end
sL = mon(pick,1); sB = mon(pick,2); sW = mon(pick,3); sH = mon(pick,4);

usableW = sW - 2*MARGIN_CM;
usableH = sH - CHROME_CM - TASKBAR_CM;
w = min(wantCm(1), max(usableW, 1));
h = min(wantCm(2), max(usableH, 1));

% Centre horizontally. Vertically, sit the canvas above the taskbar and leave the chrome room at
% the top, so the frame as a whole is inside the monitor rather than only the canvas.
left   = sL + (sW - w)/2;
bottom = sB + TASKBAR_CM + max(usableH - h, 0)/2;
set(fig, 'Position', [left bottom w h]);
% No warning here on purpose. The caller already reads the size back and reports a shortfall
% through PrimerHybridizationMap:FigureClamped, and one message about one condition is enough.
end

function [seq, name] = ParseSequenceText(txt, fallbackName)
% Plain text, FASTA (first record) or GenBank (ORIGIN block). Every non-letter is stripped, so line
% numbers, spaces and line breaks vanish; letters outside ACGTU are left in place for the caller to
% reject with their position.
name = fallbackName;
ln = regexp(txt, '\r\n|\n|\r', 'split');
lt = strtrim(ln);
if any(startsWith(lt, 'LOCUS')) && any(startsWith(lt, 'ORIGIN'))
    tok = strsplit(lt{find(startsWith(lt, 'LOCUS'), 1)});
    if numel(tok) >= 2, name = tok{2}; end
    if ~isempty(strfind(lower(lt{find(startsWith(lt, 'LOCUS'), 1)}), 'linear')) %#ok<STREMP>
        warning('PrimerHybridizationMap:LinearRecord', ...
            'The GenBank record is annotated linear; the map treats the sequence as circular.');
    end
    o = find(startsWith(lt, 'ORIGIN'), 1);
    e = find(startsWith(lt, '//'), 1, 'last');
    if isempty(e)||(e < o), e = numel(lt) + 1; end
    body = strjoin(ln(o + 1:e - 1), '');
elseif any(startsWith(lt, '>'))
    h = find(startsWith(lt, '>'));
    hdr = strtrim(lt{h(1)}(2:end));
    if ~isempty(hdr)
        tok = strsplit(hdr);
        name = tok{1};
    end
    stop = numel(lt);
    if numel(h) > 1
        stop = h(2) - 1;
        warning('PrimerHybridizationMap:MultiFASTA', 'The FASTA text holds %d records; the first is used.', numel(h));
    end
    body = strjoin(ln(h(1) + 1:stop), '');
else
    body = txt;
end
seq = strrep(upper(regexprep(body, '[^A-Za-z]', '')), 'U', 'T');
end

%% ____________________________________________________________________________________________________________
function S = ScanStrand(plusSeq, primer, strand, rxn, surveyRun, exact, endTol, useEps, cache, hw, wb0, wbSpan, seedMin, dTm)
% Every register of the primer against one strand of the circle, scored as its most stable
% contiguous union of Watson-Crick blocks. Register r covers (+) coordinates r .. r+n-1 (wrapped).
% Positions, StartPos and EndPos are (+) coordinates; StartPos -> EndPos runs clockwise.
L = numel(plusSeq); n = numel(primer);
sNum = EncodeDNA(plusSeq); pNum = EncodeDNA(primer);
compChar = 'TGCA';                                 % complement of A, C, G, T indexed by code 1..4
I = mod((0:L-1)' + (0:n-1), L) + 1;                % I(r,m): (+) coordinate of footprint column m
D = sNum(I);
if strand == '+'
    Mp  = fliplr(D == (5 - fliplr(pNum)));         % Mp(r,k): primer base k pairs with (+) base OPP(r,k)
    OPP = fliplr(I);
else
    Mp  = (D == pNum);                             % Mp(r,k): primer base k pairs with (-) base at OPP(r,k)
    OPP = I;
end

% Survey: longest Watson-Crick run of each register, vectorised over registers
runLen = zeros(L, 1); maxRun = zeros(L, 1);
for k = 1:n
    runLen = (runLen + 1).*Mp(:, k);
    maxRun = max(maxRun, runLen);
end
% EXACT-COMPLEMENT RUNS, vectorised over registers in the same style as the survey above.
% end3(r)   length of the contiguous Watson-Crick run ENDING at primer position n. The walk starts
%           at the terminus and stops at the first unpaired position, so a register whose 3' base
%           is unpaired scores 0 and never qualifies, which is the point of the criterion.
% intRun(r) longest contiguous run lying wholly inside primer positions 1..n-1, with intEnd(r) the
%           position it ends at, so it starts at intEnd - intRun + 1. Ties keep the 5'-most run:
%           arbitrary, but fixed, so the label does not move between runs.
% A LONG 3' RUN SATISFIES BOTH CRITERIA, and that is the definition working rather than failing.
% A run of 10 reaching the terminus leaves 9 inside positions 1..n-1, so the site is flagged on
% both counts and the map says so with one marker rather than two.
end3 = zeros(L, 1); alive = true(L, 1);
for k = n:-1:1
    alive = alive&Mp(:, k);
    end3  = end3 + alive;
end
intRun = zeros(L, 1); intEnd = zeros(L, 1); rl = zeros(L, 1);
for k = 1:(n - 1)
    rl  = (rl + 1).*Mp(:, k);
    hit = rl > intRun;
    intRun(hit) = rl(hit);
    intEnd(hit) = k;
end

% THE SURVEY MUST NOT HIDE A SEED. surveyRun gates which registers are scored at all, and a seed
% criterion shorter than it would be unreachable: the register carrying the seed would be dropped
% before anything looked at it. At the shipped settings this costs nothing, since 3 is already
% below both seed lengths, and it only bites if surveyRun is ever raised past them.
if (nargin < 13)||isempty(seedMin), seedMin = surveyRun; end
cand = find(maxRun >= min(surveyRun, seedMin));
nC = numel(cand);

S = struct('Strand', repmat(strand, nC, 1), 'Register', zeros(nC, 1), 'PrimerFrom', zeros(nC, 1), ...
    'PrimerTo', zeros(nC, 1), 'StartPos', zeros(nC, 1), 'EndPos', zeros(nC, 1), 'Len', zeros(nC, 1), ...
    'dG', zeros(nC, 1), 'dG2', zeros(nC, 1), 'ThreePrimePaired', false(nC, 1), ...
    'dGanch', inf(nC, 1), 'AnchFrom', zeros(nC, 1), 'AnchTo', zeros(nC, 1), ...
    'AnchStart', zeros(nC, 1), 'AnchEnd', zeros(nC, 1), 'AnchLen', zeros(nC, 1), ...
    'AnchMid', zeros(nC, 1), 'Anchored', false(nC, 1), 'TerminalPaired', false(nC, 1), ...
    'Perfect', false(nC, 1), 'dGanch2', inf(nC, 1), 'dGanchEff', inf(nC, 1), 'Eps', ones(nC, 1), ...
    'End3Run', zeros(nC, 1), 'IntRun', zeros(nC, 1), 'IntFrom', zeros(nC, 1), 'IntTo', zeros(nC, 1), ...
    'Seed3Start', zeros(nC, 1), 'Seed3End', zeros(nC, 1), 'Seed3Mid', zeros(nC, 1), ...
    'IntStart', zeros(nC, 1), 'IntEnd', zeros(nC, 1), 'IntMid', zeros(nC, 1));
S.EpsPair = cell(nC, 1);
S.Positions = cell(nC, 1); S.Segment = cell(nC, 1); S.Pairing = cell(nC, 1); S.Opposing = cell(nC, 1);
S.AnchSegment = cell(nC, 1); S.AnchPairing = cell(nC, 1); S.AnchOpposing = cell(nC, 1);
if strand == '+', strandTxt = '(+)'; else, strandTxt = '(-)'; end
upd = max(1, round(nC/100));
m = 0;
for q = 1:nC
    r = cand(q);
    ip  = Mp(r, :);
    idx = OPP(r, :);
    if strand == '+'
        opp = plusSeq(idx);                        % the (+) base physically opposing each primer base
    else
        opp = compChar(sNum(idx));                 % the (-) base physically opposing each primer base
    end
    % Candidate edges. A duplex ends on a base pair, so both edges must be paired. Taking them at
    % every paired position reproduces BestSubRegister.m; taking them at block boundaries only is
    % faster and slightly conservative, since trimming into a block can only remove Watson-Crick
    % stacks, each of which is still stabilising at 60 C.
    if exact
        eL = find(ip); eR = eL;
    else
        d  = diff([0, ip, 0]);
        eL = find(d == 1); eR = find(d == -1) - 1; % maximal Watson-Crick blocks, primer coordinates
    end
    best = inf; bf = 0; bl = 0;
    for f = eL
        for l = eR(eR >= f + 1)
            % at least one Watson-Crick stack inside, else there is no duplex to score
            if ~any(ip(f:l - 1)&ip(f + 1:l)), continue; end
            key = [primer(f:l) '/' opp(f:l)];
            if isKey(cache, key)
                g = cache(key);
            else
                g = dGMismatch(primer(f:l), opp(f:l), rxn.monv, rxn.biv, rxn.dNTP, rxn.Ta);
                cache(key) = g;
            end
            if g < best, best = g; bf = f; bl = l; end
        end
    end

    % 3'-REACHABLE OPTIMUM. Only a duplex the polymerase can start from can be extended (Huang et
    % al. 1992 NAR 20:4567). The right edge is therefore pinned to the outermost Watson-Crick pair,
    % exactly as BestSubRegister.m pins it, and the register qualifies only when that pair lies
    % within endTol of position n. This is a DIFFERENT duplex from the free-edge optimum above and
    % usually a less stable one; the free-edge value answers "how tightly does the primer stick
    % here", this one answers "can it prime here".
    % endTol carries the enzyme. At 1 it is the Taq rule, that the terminus must be in the helix.
    % Widened to the primer length it is the proofreading rule, that any overhang is trimmed off and
    % the duplex starts at the outermost pair. The pinning is what makes the second case right: the
    % exonuclease stops at the first paired base rather than chewing into the helix, so the object
    % built here is the trimmed duplex and not the free-edge optimum.
    ga = inf; af = 0;
    jR = find(ip, 1, 'last');
    if ~isempty(jR)&&(jR >= n - endTol)
        for f = eL(eL <= jR - 1)
            if ~any(ip(f:jR - 1)&ip(f + 1:jR)), continue; end
            key = [primer(f:jR) '/' opp(f:jR)];
            if isKey(cache, key)
                g = cache(key);
            else
                g = dGMismatch(primer(f:jR), opp(f:jR), rxn.monv, rxn.biv, rxn.dNTP, rxn.Ta);
                cache(key) = g;
            end
            if g < ga, ga = g; af = f; end
        end
    end
    if bf > 0
        m = m + 1;
        S.Register(m) = r; S.PrimerFrom(m) = bf; S.PrimerTo(m) = bl; S.Len(m) = bl - bf + 1;
        S.dG(m) = best; S.ThreePrimePaired(m) = (bl == n);
        % SECOND POINT FOR THE MELTING TEMPERATURE. dG(T) is linear in T for a fixed structure,
        % so one further evaluation of THIS duplex fixes dH and dS. Taken only for the duplex
        % that won, not for every union the survey scored, and deliberately uncached: the cache
        % is keyed on sequence alone and holds energies at rxn.Ta.
        S.dG2(m) = dGMismatch(primer(bf:bl), opp(bf:bl), rxn.monv, rxn.biv, rxn.dNTP, rxn.Ta + dTm);
        S.Positions{m} = idx(bf:bl);
        if strand == '+'
            S.StartPos(m) = idx(bl); S.EndPos(m) = idx(bf);
        else
            S.StartPos(m) = idx(bf); S.EndPos(m) = idx(bl);
        end
        pat = repmat('.', 1, bl - bf + 1); pat(ip(bf:bl)) = '|';
        S.Segment{m} = primer(bf:bl); S.Pairing{m} = pat; S.Opposing{m} = opp(bf:bl);
        S.Perfect(m) = (nnz(ip) == n);
        S.TerminalPaired(m) = ip(n);
        % Seed runs and their footprints on the ring. Template coordinates follow the same frame as
        % StartPos/EndPos above: on (+) the primer runs antiparallel, so the 3'-most primer position
        % carries the LOWEST (+) coordinate and the span is written idx(last) -> idx(first).
        S.End3Run(m) = end3(r);
        S.IntRun(m)  = intRun(r);
        if end3(r) > 0
            f3 = n - end3(r) + 1;
            if strand == '+'
                S.Seed3Start(m) = idx(n);  S.Seed3End(m) = idx(f3);
            else
                S.Seed3Start(m) = idx(f3); S.Seed3End(m) = idx(n);
            end
            S.Seed3Mid(m) = mod(S.Seed3Start(m) - 1 + end3(r)/2, L);
        end
        if intRun(r) > 0
            li = intEnd(r); fi = li - intRun(r) + 1;
            S.IntFrom(m) = fi; S.IntTo(m) = li;
            if strand == '+'
                S.IntStart(m) = idx(li); S.IntEnd(m) = idx(fi);
            else
                S.IntStart(m) = idx(fi); S.IntEnd(m) = idx(li);
            end
            S.IntMid(m) = mod(S.IntStart(m) - 1 + intRun(r)/2, L);
        end
        if af > 0
            S.Anchored(m) = true; S.dGanch(m) = ga;
            S.dGanch2(m) = dGMismatch(primer(af:jR), opp(af:jR), rxn.monv, rxn.biv, rxn.dNTP, rxn.Ta + dTm);
            S.AnchFrom(m) = af; S.AnchTo(m) = jR; S.AnchLen(m) = jR - af + 1;
            if strand == '+'
                S.AnchStart(m) = idx(jR); S.AnchEnd(m) = idx(af);
            else
                S.AnchStart(m) = idx(af); S.AnchEnd(m) = idx(jR);
            end
            S.AnchMid(m) = mod(S.AnchStart(m) - 1 + S.AnchLen(m)/2, L);
            patA = repmat('.', 1, jR - af + 1); patA(ip(af:jR)) = '|';
            S.AnchSegment{m} = primer(af:jR); S.AnchPairing{m} = patA; S.AnchOpposing{m} = opp(af:jR);
            % EXTENSION EFFICIENCY. Occupancy of the anchored duplex says how often the primer sits
            % here with its 3' end inside the helix. It does not say whether Taq will add a
            % nucleotide, and for a mispaired 3' region the answer spans four decades. On the
            % shipped table (Rejali NA, Moric E, Wittwer CT, Clin Chem 2018;64(5):801-809) a
            % terminal mismatch costs between 2.74 and 6.56 kcal/mol at 58 C, one base in between
            % 1.79 and 5.07, two bases in between 0.70 and 3.06. epsilon depends on the last three
            % primer positions and the bases opposing them, so it is a property of the REGISTER and
            % does not shift which left edge won above. The penalty is written in the same currency
            % as everything else, dGeff = dG - RT*log(epsilon).
            %
            % THE MAP MUST CHARGE WHAT THE SEARCH CHARGED. These three lines mirror TerminalEps in
            % MisprimerMismatch.m, down to the depth and to the format of the pair string, because
            % a map that prices a site differently from the gate that passed it is worse than no
            % map at all.
            k0 = max(1, n - 2);
            if useEps
                [S.Eps(m), ei] = ExtensionEfficiency(primer(k0:n), opp(k0:n));
                if ei.GoverningPos > 0
                    S.EpsPair{m} = sprintf('%s@%d', ei.GoverningPair, ei.GoverningPos);
                else
                    S.EpsPair{m} = ei.GoverningPair;
                end
            else
                S.Eps(m) = 1;
                S.EpsPair{m} = [primer(n) '.' opp(n)];
            end
            S.dGanchEff(m) = ga - rxn.R*(rxn.Ta + 273.15)*log(S.Eps(m));
        else
            S.AnchSegment{m} = ''; S.AnchPairing{m} = ''; S.AnchOpposing{m} = '';
            S.EpsPair{m} = '';
        end
    end
    if (mod(q, upd) == 0)&&~isempty(hw)&&ishghandle(hw)
        waitbar(wb0 + wbSpan*q/nC, hw, sprintf('Scanning the %s strand: register %d of %d', strandTxt, q, nC));
    end
end
fn = fieldnames(S);
for k = 1:numel(fn)
    S.(fn{k}) = S.(fn{k})(1:m, :);
end
end

%% ____________________________________________________________________________________________________________
function c = EncodeDNA(s)
[~, c] = ismember(s, 'ACGT');
end

%% ____________________________________________________________________________________________________________
function S = CatSites(A, B)
S = struct();
fn = fieldnames(A);
for k = 1:numel(fn)
    S.(fn{k}) = [A.(fn{k}); B.(fn{k})];
end
end

%% ____________________________________________________________________________________________________________
function R = RegionTable(best, bestTm, bestPf, bestPt, strand)
% Connected runs of displayed nucleotides on one circular strand, with the wrap at position L -> 1.
L = numel(best);
mask = isfinite(best);
R = struct('strand', strand, 's', zeros(0, 1), 'e', zeros(0, 1), 'len', zeros(0, 1), ...
    'dG', zeros(0, 1), 'Tm', zeros(0, 1), 'pf', zeros(0, 1), 'pt', zeros(0, 1), ...
    'xMid', zeros(0, 1));
if ~any(mask), return; end
if all(mask)
    s = 1; e = L;
else
    k0 = find(~mask, 1);
    order = [k0:L, 1:k0 - 1];
    d = diff([0, mask(order), 0]);
    s = order(find(d == 1))'; e = order(find(d == -1) - 1)';
end
len = mod(e - s, L) + 1;
dG = zeros(numel(s), 1); Tm = nan(numel(s), 1);
pf = nan(numel(s), 1); pt = nan(numel(s), 1);
for k = 1:numel(s)
    idx = mod((s(k) - 1):(s(k) + len(k) - 2), L) + 1;
    [dG(k), j] = min(best(idx));
    Tm(k) = bestTm(idx(j));
    pf(k) = bestPf(idx(j));
    pt(k) = bestPt(idx(j));
end
R.s = s(:); R.e = e(:); R.len = len(:); R.dG = dG; R.Tm = Tm; R.pf = pf; R.pt = pt;
R.xMid = mod(s(:) - 1 + len(:)/2, L);
end

%% ____________________________________________________________________________________________________________
function a = XToAngle(x, L)
% Nucleotide coordinate x (0 = start of position 1) to polar angle; position 1 at 12 o'clock, clockwise.
a = pi/2 - 2*pi*x/L;
end

%% ____________________________________________________________________________________________________________
function DrawAnnulus(ax, rIn, rOut, rgb)
a = linspace(0, 2*pi, 1441);
patch(ax, 'XData', [rOut*cos(a), rIn*cos(fliplr(a))], 'YData', [rOut*sin(a), rIn*sin(fliplr(a))], ...
    'FaceColor', rgb, 'EdgeColor', 'none');
end

%% ____________________________________________________________________________________________________________
function DrawNucleotides(ax, idx, rgb, rIn, rOut, L)
% One flat-coloured wedge per nucleotide. Each wedge overlaps its neighbours by 0.1 nt, which hides the
% anti-aliasing seam between faces and moves no region edge by more than 0.1 nt.
if isempty(idx), return; end
nsub = max(1, ceil((360/L)/0.5));                  % chords no longer than 0.5 degree
ov = 0.1;
t = (0:nsub)/nsub*(1 + 2*ov) - ov;
ang = XToAngle((idx(:) - 1) + t, L);               % numel(idx) x (nsub+1)
X = [rOut*cos(ang), rIn*cos(fliplr(ang))];
Y = [rOut*sin(ang), rIn*sin(fliplr(ang))];
nv = size(X, 2);
V = [reshape(X.', [], 1), reshape(Y.', [], 1)];
F = reshape(1:numel(X), nv, []).';
patch(ax, 'Faces', F, 'Vertices', V, 'FaceVertexCData', rgb, 'FaceColor', 'flat', 'EdgeColor', 'none');
end

%% ____________________________________________________________________________________________________________
function DrawDirectionLabel(ax, best, r, sym, isPlus, L, fs)
w = L*9/360;                                        % half-width of the label arc, nucleotides
xTop = mod(round(-w):round(w), L) + 1;
xBot = mod(round(L/2 - w):round(L/2 + w), L) + 1;
atTop = ~any(isfinite(best(xTop)))||any(isfinite(best(xBot)));
pr = char(8242);
if atTop
    y = r;
    if isPlus, arrow = ['5' pr '\rightarrow3' pr]; else, arrow = ['3' pr '\leftarrow5' pr]; end
else
    y = -r;
    if isPlus, arrow = ['3' pr '\leftarrow5' pr]; else, arrow = ['5' pr '\rightarrow3' pr]; end
end
if (nargin < 7)||isempty(fs), fs = 10; end
text(ax, 0, y, [sym '  ' arrow], 'FontSize', fs, 'Color', [0.25 0.25 0.25], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'tex');
end

%% ____________________________________________________________________________________________________________
function DrawMarker(ax, xMid, L, rTip, rBase, rgb)
a = XToAngle(xMid, L);
da = 0.016/rBase;
patch(ax, 'XData', [rTip*cos(a), rBase*cos(a - da), rBase*cos(a + da)], ...
    'YData', [rTip*sin(a), rBase*sin(a - da), rBase*sin(a + da)], 'FaceColor', rgb, ...
    'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 0.3);
end

%% ____________________________________________________________________________________________________________
function e = LabelExtents(h)
% Extent of every label, one row of [x y w h] per handle, in data units. get() hands back a bare
% numeric row for a single handle and a cell array for several, so both shapes are normalised here
% rather than at each of the three call sites.
e = get(h, 'Extent');
if ~iscell(e), e = {e}; end
e = vertcat(e{:});
end

%% ____________________________________________________________________________________________________________
function tf = AnyOverlap(e)
% True when any two of the boxes in e intersect. Boxes that merely touch do not count: labels
% pushed to exactly abutting are as far apart as the ring can usefully put them.
tf = false;
n = size(e, 1);
for i = 1:n - 1
    for j = i + 1:n
        if (e(i, 1) < e(j, 1) + e(j, 3))&&(e(j, 1) < e(i, 1) + e(i, 3)) ...
                &&(e(i, 2) < e(j, 2) + e(j, 4))&&(e(j, 2) < e(i, 2) + e(i, 4))
            tf = true;
            return
        end
    end
end
end

%% ____________________________________________________________________________________________________________
function d = LabelSepDeg(ax, nLines, fontPt, rLab, lineFactor)
% Angular separation, in degrees, that clears a label nLines tall at radius rLab.
%
% Two labels overprint when the arc between them is shorter than a label is tall, so the quantity
% that matters is a length in centimetres on the page, not an angle. The label's height is set by
% the font and the line count; the arc per degree is set by the scale of the axes. Both are read
% here rather than assumed, which is why this takes the axes rather than a number: the scale
% changes when the limits open to clear a wide label, and the separation has to follow it.
%
% axis equal holds the data aspect at 1:1, so whichever of width or height binds sets the scale
% for both; min() picks it.
u = get(ax, 'Units'); set(ax, 'Units', 'centimeters');
p = get(ax, 'Position'); set(ax, 'Units', u);
cmPerUnit = min(p(3)/diff(xlim(ax)), p(4)/diff(ylim(ax)));
hCm   = nLines*fontPt*lineFactor/28.3465;        % 1 pt = 1/28.3465 cm
arcCm = rLab*(pi/180)*cmPerUnit;                 % cm of arc per degree at the label radius
if arcCm <= 0, d = 0; else, d = hCm/arcCm; end
end

%% ____________________________________________________________________________________________________________
function cw = SpreadAngles(cw, minSep)
% Clockwise label angles (degrees) pushed apart to at least minSep, preserving their circular order.
% The pass starts just after the widest empty arc, so labels either side of 12 o'clock are separated too.
nL = numel(cw);
if nL < 2, return; end
[s, o] = sort(mod(cw, 360));
gaps = [diff(s), s(1) + 360 - s(end)];
[~, g] = max(gaps);
idx = [g + 1:nL, 1:g];
t = s(idx);
t(t < t(1)) = t(t < t(1)) + 360;
for k = 2:nL
    t(k) = max(t(k), t(k - 1) + minSep);
end
s(idx) = t;
cw(o) = s;
end

%% ____________________________________________________________________________________________________________
function [ha, va] = LabelAlignment(a)
c = cos(a); s = sin(a);
if c > 0.25
    ha = 'left';
elseif c < -0.25
    ha = 'right';
else
    ha = 'center';
end
if strcmp(ha, 'center')
    if s > 0, va = 'bottom'; else, va = 'top'; end
else
    va = 'middle';
end
end

%% ____________________________________________________________________________________________________________
function LayoutCentreText(ax, lines, fs, R, fontName, scale, fsMax)
% Stacks one text object per line about the origin and sizes the block so that the outer corners of
% every line lie inside the circle of radius R. Line pitch is set from the font size (1.3 em, plus
% 0.25 em for lines carrying a sub- or superscript), converted to data units through the plot box,
% so it does not depend on how tightly the renderer reports the text extent. Blank lines take half
% a pitch.
%
% IT FITS IN BOTH DIRECTIONS NOW. This used to shrink only, in 4 per cent steps from fs down to a
% 5 pt floor, which meant a block that already fitted at fs stayed at fs however much room was left
% inside R. The panel was therefore as large as a constant said rather than as large as the circle
% allowed, and every time a line was shortened the middle of the figure got emptier rather than the
% text getting bigger. Phase one comes down until it fits. Phase two goes back up in 2 per cent
% steps and stops at the last size that still fits, or at fsMax.
%
% THE TWO PHASES CANNOT BOTH RUN TO EFFECT. Phase two is entered only when phase one succeeded, so
% a block that bottomed out at the floor is never grown back into the rings it could not clear.
lines = lines(:);
% PER-LINE SIZE. scale multiplies the base size line by line, so a heading can be set larger than
% the body and still take part in the same shrink. Written as a factor rather than as a \fontsize
% directive inside the string, because the loop below rescales fs and a literal size would not
% follow it down.
if (nargin < 6)||isempty(scale)
    scale = ones(numel(lines), 1);
else
    scale = scale(:);
    if numel(scale) ~= numel(lines)
        error('LayoutCentreText:ScaleLength', ...
            'scale has %d entries for %d lines.', numel(scale), numel(lines));
    end
end
if (nargin < 7)||isempty(fsMax), fsMax = 2*fs; end
if fsMax < fs, fsMax = fs; end
blank = cellfun(@isempty, lines);
tall = ~cellfun(@isempty, regexp(lines, '[\^_]\{', 'once'));
oldUnits = get(ax, 'Units');
set(ax, 'Units', 'points'); pp = get(ax, 'Position'); set(ax, 'Units', oldUnits);
upp = max(diff(get(ax, 'XLim'))/pp(3), diff(get(ax, 'YLim'))/pp(4));   % data units per point
hs = [];
for k = 1:numel(lines)
    hs(k, 1) = text(ax, 0, 0, lines{k}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'Interpreter', 'tex', 'FontName', fontName, 'FontSize', fs*scale(k)); %#ok<AGROW>
end

% ---- phase one, come down until the block fits -------------------------------------------------
ok = false; yc = zeros(numel(lines), 1);
for it = 1:60
    [ok, yc, fits, wid] = CentreFits(hs, fs, scale, upp, tall, blank, R);
    if ok, break; end
    % THE FLOOR IS REPORTED, NOT JUST OBEYED. Bottoming out at 5 pt with lines still outside R used
    % to end the loop in silence and hand back a panel nobody could read, which is what a 138
    % character table path did on 07-11-2026. Naming the widest offender turns an unreadable figure
    % into one line of diagnosis.
    if fs <= 5
        off = ~fits&~blank;
        [~, worst] = max(wid(off));
        bad = lines(off);
        warning('PrimerHybridizationMap:CentreTextFloor', ...
            ['The centre panel hit the 5 pt floor with %d of %d lines still reaching outside ' ...
            'TextRadius, so it is set smaller than it should be. The widest is "%s". Shorten ' ...
            'it, or raise opts.TextRadius.'], sum(off), numel(lines), bad{worst});
        break
    end
    fs = 0.96*fs;
    SetCentreSizes(hs, fs, scale);
end

% ---- phase two, go back up until it stops fitting -----------------------------------------------
if ok
    for it = 1:120
        if fs >= fsMax, break; end
        trial = min(fs*1.02, fsMax);
        SetCentreSizes(hs, trial, scale);
        [ok2, yc2] = CentreFits(hs, trial, scale, upp, tall, blank, R);
        if ~ok2
            SetCentreSizes(hs, fs, scale);    % step back to the last size that fitted
            break
        end
        fs = trial; yc = yc2;
    end
end

for k = 1:numel(hs)
    set(hs(k), 'Position', [0 yc(k) 0]);
end
end

%% ____________________________________________________________________________________________________________
function [ok, yc, fits, wid] = CentreFits(hs, fs, scale, upp, tall, blank, R)
%% CENTREFITS  Does the stacked block lie inside the circle of radius R at this font size
% Returns the per-line vertical centres as well, because they are computed here and the caller
% needs the set that belongs to the size it ends up keeping. Blank lines are exempt: they carry no
% ink, so a blank wider than R is not a thing that can happen.
ex  = cell2mat(get(hs, 'Extent'));
wid = ex(:, 3);
hk  = fs*scale.*upp.*(1.3 + 0.25*tall);
hk(blank) = 0.65*fs*upp;
yc  = sum(hk)/2 - [0; cumsum(hk(1:end-1))] - hk/2;
fits = hypot(wid/2, abs(yc) + hk/2) <= R;
ok  = all(fits|blank);
end

%% ____________________________________________________________________________________________________________
function SetCentreSizes(hs, fs, scale)
%% SETCENTRESIZES  One font size per line, kept proportional
% LINE BY LINE, AND IT HAS TO BE. A single set(hs, 'FontSize', fs) flattens every line to the same
% size and throws the per-line scale away on the first step, taking the title down with it.
for k = 1:numel(hs)
    set(hs(k), 'FontSize', fs*scale(k));
end
end

%% ____________________________________________________________________________________________________________
function s = TexSafe(t)
%% TEXSAFE  Make free text safe for the TeX interpreter every centre-panel line runs under
% BACKSLASH FIRST, AND IT BECOMES A FORWARD SLASH. MATLAB's TeX mode has no way to draw a literal
% backslash, and a Windows path carries one before nearly every element: \O, \C and \a are all read
% as control sequences. The interpreter then abandons the string and draws it verbatim, which is
% how the caption came to be rendered as the fourteen letters that spell \epsilon followed by a
% path whose every underscore showed its own escape. Forward slashes keep the string readable and
% keep the interpreter alive.
%
% THE REMAINING FOUR ARE ESCAPED AFTER, NOT BEFORE. Escaping an underscore INTRODUCES a backslash,
% so a backslash pass run second would escape the escapes.
s = strrep(char(t), '\', '/');
for c = {'_', '^', '{', '}'}
    s = strrep(s, c{1}, ['\' c{1}]);
end
end

%% ____________________________________________________________________________________________________________
function s = Tern(c, a, b)
if c, s = a; else, s = b; end
end

%% ____________________________________________________________________________________________________________
function s = ValueText(v)
if ischar(v), s = v; else, s = sprintf('%.6g', v); end
end

%% ____________________________________________________________________________________________________________
function FitText(h, fits)
% Shrinks a text object in 3% steps until fits(Extent) holds or the font reaches 5 pt.
for it = 1:40
    if fits(get(h, 'Extent'))||(get(h, 'FontSize') <= 5), break; end
    set(h, 'FontSize', 0.97*get(h, 'FontSize'));
end
end

%% ____________________________________________________________________________________________________________
function x = ReadNumber(str, label, lo, hi)
x = str2double(str);
if ~isfinite(x)||(x < lo)||(x > hi)
    error('PrimerHybridizationMap:BadInput', '%s must be a number in [%g, %g]; received ''%s''.', ...
        label, lo, hi, str);
end
end

%% ____________________________________________________________________________________________________________
function c = SiteLabel(pre, tStart, tEnd, pFrom, pTo, dG, tm)
% The three lines of a site label, used for displayed regions and for extension threats alike so
% the two read the same way round the ring.
%
%   line 1   strand, the TEMPLATE coordinates, then the PRIMER positions that pair there
%   line 2   free energy of that duplex
%   line 3   its melting temperature at the free primer
%
% The two ranges answer different questions and are both needed: the template range says where on
% the plasmid the primer is sitting, the primer range says how much of the primer is engaged and
% from which end. A seven-base duplex at primer 1-7 leaves the 3'-OH outside the helix and cannot
% prime; the same duplex at primer 13-19 can. The label carries the distinction rather than making
% a reader open the Sites sheet for it.
%
% \DeltaG_{hyb} matches the colour-bar label rather than the centre panel's \DeltaG_{Hyb}, since
% the colour bar is the scale a reader checks the number against.
%
% RETURNS A PLAIN CELL OF LINES. Every caller wraps it once, as struct('txt', {SiteLabel(...)}),
% which is the same idiom the seed labels use with their literal {{t1S, t2S}}. struct() consumes
% one level of cell to mean "one struct element per entry", so returning an already-wrapped cell
% here and wrapping it again at the call site would bury the lines one level too deep and text()
% would reject the result.
t1 = [pre RangeTxt(tStart, tEnd)];
if isfinite(pFrom)&&isfinite(pTo)
    t1 = [t1 ' | ' RangeTxt(pFrom, pTo)];
end
t2 = ['\DeltaG_{hyb} = ' FmtSigned(dG, 2)];
lines = {t1, t2};
if isfinite(tm)
    lines{end + 1} = sprintf('T_{m} = %.1f ^{\\circ}C', tm);
end
% THE LINE IS DROPPED, NOT BLANKED. A degenerate two-point fit gives no temperature, and the
% obvious move is to emit '' for that line. MATLAB's text() rejects an empty character vector
% inside a cell array outright, which aborts the whole figure at the last drawing step, after
% the scan has already run. Such a label is therefore two lines rather than three with a hole.
c = lines;
end

%% ____________________________________________________________________________________________________________
function s = RangeTxt(a, b)
s = sprintf('%d%s%d', a, char(8211), b);
end

%% ____________________________________________________________________________________________________________
function s = FmtSigned(x, d)
s = sprintf(['%.' num2str(d) 'f'], x);
if x < 0, s = [char(8722) s(2:end)]; end
end

%% ____________________________________________________________________________________________________________
function s = Num(x, sig)
s = sprintf(['%.' num2str(sig) 'g'], x);
end

%% ____________________________________________________________________________________________________________
function s = TexEscape(s)
s = regexprep(s, '([\\_^{}])', '\\$1');
end

%% ____________________________________________________________________________________________________________
function st = NiceStep(x)
p = 10^floor(log10(x));
for mlt = [1 2 2.5 5 10]
    st = mlt*p;
    if st >= x, return; end
end
end
