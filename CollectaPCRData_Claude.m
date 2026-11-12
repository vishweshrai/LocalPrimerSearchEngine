function [input, decider, mode] = CollectaPCRData_Claude(masterpath)
%% COLLECTAPCRDATA_CLAUDE  Reaction, template and primer inputs for the aPCR suite
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 05-11-2026 21:10 IST
% Prompts for the plasmid and scaffold, the reaction conditions and the polymerase, converts every
% entry to SI-consistent units, and returns the input struct the engine and the simulator both read.
% Units in: ng/uL, uM, mM, mU/uL, nmol. Units out: g/L and M throughout, U/uL for activity.
input = struct('PlasmidFileName', '', ...
    'PlasmidSequence', '', ...
    'PlasmidStranded', NaN, ...
    'PlasmidConcentration', NaN, ...
    'ScaffoldStart', NaN, ...
    'ScaffoldEnd', NaN, ...
    'ScaffoldLength', NaN, ...
    'ForwardPrimer', '', ...
    'ReversePrimer', '', ...
    'ForwardPrimerConcentration', NaN, ...
    'ReversePrimerConcentration', NaN, ...
    'MonovalentCationConcentration', NaN, ...
    'BivalentCationConcentration', NaN, ...
    'dNTPsConcentration', NaN, ...
    'AnnealingTemperature', NaN,...
    'ExtensionTemperature', NaN, ...
    'NumAmplificationCycles', NaN, ...
    'PolymeraseName', '', ...
    'PolymeraseActivity', NaN, ...
    'ExtensionTime', NaN, ...
    'UnitDNTP', NaN, ...
    'UnitTime', NaN, ...
    'PolymeraseEfficiency', NaN, ...
    'MaxLaps', NaN, ...
    'ProofreadingExo', NaN, ...
    'ProofreadingClass', '', ...
    'FivePrimeExo', NaN);
% masterpath = fileparts(masterpath);
[plasfilename, plasfilelocation] = uigetfile({'*.txt';}, 'Select the Custom Genome Sequence File',[masterpath, '\Plasmids\']);
input.PlasmidFileName = plasfilename; input.PlasmidSequence = char(extractFileText(fullfile(plasfilelocation, plasfilename)));
SeqCheck = upper(char(input.PlasmidSequence)); bad = ~ismember(SeqCheck, 'ACGT');
if any(bad)
    firstBad = find(bad, 1);
    error('LocalPrimerSearchEngine:NonACGT', ...
        ['PlasmidSequence carries %d character(s) outside {A,C,G,T}; the first sits at ' ...
        'index %d and is ''%s''. Nearest-neighbour tables and other peripheral functions are only ' ...
        'defined on unambiguous DNA. Clean the sequence at loading time in collectaPCRdata.m. ' ...
        'Stripping characters here would desynchronise ScaffoldStart and ScaffoldEnd.'], ...
        nnz(bad), firstBad, SeqCheck(firstBad));
else
    input.PlasmidSequence = SeqCheck;
end
promptb = {'Is it double-stranded (2) or single stranded (1)?', 'Enter plasmid concentration (ng/uL):'};
dlgtitle = 'Input Plasmid Details';
fieldsize = [1 50; 1 50];
definput = {'2','1'};
dna = inputdlg(promptb,dlgtitle,fieldsize,definput);
input.PlasmidStranded = str2double(dna{1}); plasdens = str2double(dna{2})*1E-3;
an = count(input.PlasmidSequence,'A'); tn = count(input.PlasmidSequence,'T');
cn = count(input.PlasmidSequence,'C'); gn = count(input.PlasmidSequence,'G');
input.PlasmidConcentration = plasdens/(input.PlasmidStranded*((an*313.2)+(tn*304.2)+(cn*289.2)+ ...
    (gn*329.2)+79.0));
promptc = {'Do you have the primers already designed?: (Y/N)'};
dlgtitle = 'Primer Details';
fieldsize = [1 50];
definput = {'N'};
primdes = inputdlg(promptc,dlgtitle,fieldsize,definput);
if ismember(primdes, {'N', 'n', 'No', 'no'})
    promptd = {'Do you want the primers to be designed for a scaffold of particular sequence?: (Y/N)'};
    dlgtitle = 'Primer Details';
    fieldsize = [1 50];
    definput = {'Y'};
    scafseq = inputdlg(promptd,dlgtitle,fieldsize,definput);
    if ismember(scafseq, {'Y', 'y', 'Yes', 'yes'})
        promptj = {'Do you have the desired scaffold sequence?'};
        dlgtitle = 'Scaffold Details';
        fieldsize = [1 50];
        definput = {'N'};
        scafseq = inputdlg(promptj,dlgtitle,fieldsize,definput);
        if ismember(scafseq, {'Y', 'y', 'Yes', 'yes'})
            [scaffilename, scaffilelocation] = uigetfile({'*.txt';}, 'Select the Scaffold Sequence File',[masterpath, '\Plasmids\']);
            scaffseq = char(extractFileText(fullfile(scaffilelocation, scaffilename)));
            decider = 'Local'; mode = 'Sequence';
            input.ScaffoldLength = strlength(scaffseq);
            [~, Alignment, Start] = swalign(input.PlasmidSequence, scaffseq, ALPHABET = 'NT');
            if contains(Alignment(2,:), ' ') || contains(Alignment(2,:), '.')
                fprintf('Scaffold sequence has mismatches, and cannot be processed further!\n');
                return
            elseif strlength(Alignment(3,:)) == input.ScaffoldLength
                input.ScaffoldStart = Start(1,1); input.ScaffoldEnd = input.ScaffoldStart + input.ScaffoldLength - 1;
            else
                [~, Alignment, Start] = swalign([input.PlasmidSequence input.PlasmidSequence], scaffseq, ALPHABET ...
                    = 'NT');
                if strlength(Alignment(3,:)) == input.ScaffoldLength
                    input.ScaffoldStart = Start(1,1); input.ScaffoldEnd = input.ScaffoldStart + input.ScaffoldLength - ...
                        strlength(input.PlasmidSequence) - 1;
                else
                    fprintf('Scaffold sequence has mismatches, and cannot be processed further!\n');
                    return
                end
            end
        else
            prompte = {'Input the start and end indices of the scaffold?'};
            dlgtitle = 'Scaffold Details';
            fieldsize = [1 50];
            definput = {'[2173 6112]'};
            seq = inputdlg(prompte,dlgtitle,fieldsize,definput);
            input.ScaffoldStart = sscanf(seq{1}, '[%d', 1); input.ScaffoldEnd = sscanf(seq{1}, '[%*d %d]', 1);
            if input.ScaffoldStart > input.ScaffoldEnd
                input.ScaffoldLength = input.ScaffoldEnd + strlength(input.PlasmidSequence) - input.ScaffoldStart + 1;
            else
                input.ScaffoldLength = input.ScaffoldEnd - input.ScaffoldStart + 1;
            end
            decider = 'Local'; mode = 'Positions';
        end
    else
        promptf = {'Input the length of the scaffold?'};
        dlgtitle = 'Scaffold Details';
        fieldsize = [1 50];
        definput = {'3700'};
        seq = inputdlg(promptf,dlgtitle,fieldsize,definput);
        input.ScaffoldLength = str2double(seq{1});
        decider = 'Global'; mode = 'None';
    end
else
    promptg = {'Input the forward primer sequence:', 'Input the reverse primer sequence:'};
    dlgtitle = 'Primer Details';
    fieldsize = [1 50; 1 50];
    definput = {'GGAACGGTAAATTCAGAGACT';'CACATTAATTGCGTTGCGCTCAC'};
    primerseq = inputdlg(promptg,dlgtitle,fieldsize,definput);
    input.ForwardPrimer = primerseq{1}; input.ReversePrimer = primerseq{2};
    decider = 'Yield_Test'; mode = 'None';
    [~, fAlignA, Start] =  swalign(input.ForwardPrimer, input.PlasmidSequence, 'GapOpen', 1000, 'Alphabet', 'NT');
    if strlength(fAlignA(1,:)) == strlength(input.ForwardPrimer)
        input.ScaffoldStart = Start(2);
    elseif contains(fAlignA(2,:), ' ') || contains(fAlignA(2,:), '.')
        fprintf('Forward Primer sequence has mismatches, and cannot be processed further!\n');
        return
    else
        [~, ~, Start] =  swalign(input.ForwardPrimer, [input.PlasmidSequence, input.PlasmidSequence], 'GapOpen', 1000, 'Alphabet', 'NT');
        input.ScaffoldStart = Start(2);
    end
    [~, rAlignS, Start] =  swalign(input.ReversePrimer, seqrcomplement(input.PlasmidSequence), 'GapOpen', 1000, 'Alphabet', 'NT');
    if strlength(rAlignS(1,:)) == strlength(input.ReversePrimer)
        input.ScaffoldEnd = strlength(input.PlasmidSequence) - Start(2) +1 ;
    elseif contains(rAlignS(2,:), ' ') || contains(rAlignS(2,:), '.')
        fprintf('Reverse Primer sequence has mismatches, and cannot be processed further!\n');
        return
    else
        [~, ~, Start] =  swalign(input.ReversePrimer, seqrcomplement([input.PlasmidSequence, input.PlasmidSequence]), 'GapOpen', 1000, 'Alphabet', 'NT');
        input.ScaffoldEnd = strlength(input.PlasmidSequence) - Start(2) +1 ;
    end
    if input.ScaffoldStart > input.ScaffoldEnd
        input.ScaffoldLength = input.ScaffoldEnd + strlength(input.PlasmidSequence) - input.ScaffoldStart + 1;
    else
        input.ScaffoldLength = input.ScaffoldEnd - input.ScaffoldStart + 1;
    end
end
promptf = {'Enter forward primer concentration (uM):', 'Enter forward/reverse primer concentration ratio:',...
    'Enter monovalent cation concentration (mM):',...
    'Enter bivalent cation concentration (mM):', 'Enter dNTPs concentration (mM):',...
    'Enter annealing temperature (oC) [0 = Programmed, 1 = Variable, 2: Touchdown]:',...
    'Enter extension temperature (oC):', 'Enter the number of amplification cycles:'};
dlgtitle = 'Input Reaction Conditions';
fieldsize = [1 50; 1 50; 1 50; 1 50; 1 50; 1 50; 1 50; 1 50];
definput = {'1', '50', '40', '2.5', '1.2', '55', '72', '30'};
react = inputdlg(promptf,dlgtitle,fieldsize,definput);
input.ForwardPrimerConcentration = str2double(react{1})*1E-6;
input.ReversePrimerConcentration = input.ForwardPrimerConcentration/str2double(react{2});
input.MonovalentCationConcentration = str2double(react{3})*1E-3;
input.BivalentCationConcentration = str2double(react{4})*1E-3;
input.dNTPsConcentration = str2double(react{5})*1E-3;
input.ExtensionTemperature = str2double(react{7});
input.AnnealingTemperature = str2double(react{6});
input.NumAmplificationCycles = str2double(react{8});
if input.AnnealingTemperature == 2
    % The slot order below matches the default recipe. The earlier label read
    % [ramp-cycles, start, ramp, constant], which put -0.5 in the start slot and 63 in the ramp
    % slot; aPCRCycleAlgebra reads slot 2 as the per-cycle increment and slot 3 as the starting
    % temperature, and warns if the two arrive the other way round.
    promptf = {'Enter touchdown recipie [ramp-cycles, ramp, start, constant]: '};
    dlgtitle = 'Input Touchdown Recipie';
    fieldsize = [1 100];
    definput = {'[10, -0.5, 63, 58]'};
    touchdown = inputdlg(promptf,dlgtitle,fieldsize,definput);
    input.AnnealingTemperature = touchdown{1};
end
%% Polymerase Details
% The simulator caps synthesis per cycle at what the enzyme can physically turn over:
%       Cap [mol/L] = Activity [U/uL] * 1E6 [uL/L] * UnitDNTP/UnitTime [mol/(U min)] * ExtensionTime
% One unit is the amount incorporating UnitDNTP mol of dNTP into acid-insoluble material in
% UnitTime minutes at the supplier's assay temperature. Confirm both numbers against the datasheet
% of the lot you are running; 10 nmol in 30 min is the common convention but it is not universal.
%
% MaxLaps is how many complete turns the polymerase makes on a circular template before falling
% off. A circle presents no 5' terminus, so one lap always happens. A second lap needs strand
% displacement, and each further lap adds a rung to the run-around ladder. The values below are
% engineering estimates from the published properties of each enzyme, not measured numbers; treat
% them as starting points and override them here.
%
% ____________________________________________________________________________________________________________
% THE TWO NUCLEASE ACTIVITIES, AND WHY THE SEARCH ENGINE NEEDS TO KNOW ABOUT THEM
%
% These are different enzymes living in the same protein and they do unrelated things. Collecting
% them as one "exo" flag loses the distinction that matters.
%
% 3' -> 5' PROOFREADING EXONUCLEASE. This one decides whether a primer sitting on an off-target site
% with a MISPAIRED 3' terminus can be extended at all. Without it the terminus is simply not
% extended: Taq turns over a mispaired 3' end at 1E-3 to below 1E-6 of a matched one, and that
% factor is the whole basis of the extension-efficiency weighting in stage 13 of
% LocalPrimerSearchEngine.m (Huang MM, Arnheim N, Goodman MF, Nucleic Acids Res 1992;20:4567,
% doi:10.1093/nar/20.17.4567). With proofreading, the mispaired nucleotide is excised and synthesis
% proceeds, so that factor collapses towards 1 and the mispriming screen MUST fall back to pure
% occupancy. New England Biolabs state the mechanism plainly: the domain excises mismatched
% nucleotides in the 3' to 5' direction and removes unpaired 3' overhanging nucleotides to create
% blunt ends, which is exactly the geometry a mispaired primer terminus presents.
% https://www.neb.com/en/products/pcr-qpcr-and-amplification-technologies/dna-polymerase-proofreading
%
% A proofreading enzyme therefore makes the reaction LESS specific at the priming step, not more,
% whatever its advantage in copying fidelity. It also nibbles free single-stranded primers over the
% course of a run, which shortens them and lowers their Tm. Nothing in this suite models that.
%
% 5' -> 3' EXONUCLEASE. This one has no bearing on whether a mispaired terminus extends. It decides
% what happens when a polymerase extending from an upstream site reaches a primer already bound
% downstream: with the activity, that primer is degraded; without it, the polymerase either stalls
% or displaces the strand. It is the reason MaxLaps is 1 for the two Taqs and larger for the
% enzymes that lack it, so it is already encoded implicitly below and is now recorded explicitly.
%
% FOUR OF THE FIVE NAMED ENZYMES LACK 3' -> 5' PROOFREADING. Stoffel Fragment is the N-terminally
% truncated Taq that lacks the 5' -> 3' domain; Deep Vent (exo-) is the archaeal enzyme with its
% proofreading domain inactivated and no 5' -> 3' domain to begin with. NEB LongAmp Taq is the
% exception and is the only preset at class weak: its datasheet declares 3' - 5' Exonuclease Yes,
% so it reaches here from the datasheet rather than by inference from the blend composition.
% Choose 6 = Other and answer the two questions yourself for anything else, answering strong for
% Pfu, Vent (exo+), Phusion, Q5, KOD and the other high-fidelity enzymes.
promptp = {['Choose the polymerase:' newline ...
    '1 = TaKaRa Taq' newline ...
    '2 = Quantabio AccuStart II Taq' newline ...
    '3 = Qiagen/ABI Stoffel Fragment' newline ...
    '4 = NEB Deep Vent (exo-)' newline ...
    '5 = NEB LongAmp Taq (M0323)' newline ...
    '6 = Other']};
dlgtitle = 'Polymerase Selection';
fieldsize = [1 60];
definput = {'1'};
polychoice = inputdlg(promptp,dlgtitle,fieldsize,definput);
polyidx = str2double(polychoice{1});
%                1 TaKaRa   2 AccuStart   3 Stoffel   4 DeepVent(exo-)   5 LongAmp   6 Other
PolyName = {'TaKaRa Taq', 'Quantabio AccuStart II Taq', 'Stoffel Fragment', ...
    'NEB Deep Vent (exo-)', 'NEB LongAmp Taq', 'Other'};
PolyAct = {'25', '25', '60', '20', '100', '25'};   % mU/uL, from the recommended unit load per 50 uL
PolyLaps = {'1', '1', '3', '5', '1', '1'};         % see the lap note below
PolyUnitDNTP = {'10', '10', '10', '10', '10', '10'};   % nmol of dNTP per unit
PolyUnitTime = {'30', '30', '30', '30', '30', '30'};   % minutes
PolyProof    = {'none', 'none', 'none', 'none', 'weak', 'none'};   % 3'->5' exonuclease class
Poly5to3     = {'Y', 'Y', 'N', 'N', 'Y', 'Y'};     % 5'->3' exonuclease; Stoffel and Deep Vent lack it
% ---- NEB LongAmp Taq (M0323), entry 5, read off the vial and the Properties & Usage block -------
% ACTIVITY. The vial reads 2500 U/mL, which is 2.5 U/uL, and the standard protocol directs 2 uL
% into a 50 uL reaction. That is 5 U per 50 uL, and the column holds mU/uL IN THE REACTION, so
%       5 U / 50 uL = 0.1 U/uL = 100 mU/uL
% which is four times the TaKaRa Taq entry. The conversion for any new row is the same:
% mU/uL = (units per 50 uL) * 20.
% UNIT DEFINITION. 10 nmol of dNTP into acid-insoluble material in 30 minutes at 75 C, matching the
% convention the other rows already carry. The unit assay runs on 15 nM PRIMED M13 DNA, which is
% this suite's own template family, so the activity transfers without a template correction.
% LAPS ARE AN ESTIMATE HERE, AND MORE SO THAN FOR THE OTHERS. LongAmp is the only entry declaring
% BOTH 5' - 3' Exonuclease Yes AND Strand Displacement +, and the two activities pull opposite ways
% at the point where the polymerase comes round the circle and meets the 5' end of the strand it
% just laid down. The 5' -> 3' domain degrades that end and nick-translates, so the 3' terminus
% advances while the 5' terminus retreats and the product never grows past one circumference. The
% displacement activity instead peels the strand off and keeps synthesising, which is rolling
% circle and gives a product longer than the circle. Which one dominates is set by the blend ratio,
% and LongAmp is Taq-based with a minor Deep Vent component, so most encounters are Taq encounters
% and nick translation should win. Hence 1. A band or smear above the expected amplicon size on a
% circular template settles it the other way; override this field if you see one.
if isnan(polyidx)||(polyidx < 1)||(polyidx > 6)
    polyidx = 6;
end
promptq = {'Polymerase name:', 'Enter polymerase activity (mU/uL):', ...
    'Enter extension time per cycle (min):', ...
    'Enter unit definition, nmol dNTP incorporated per unit:', ...
    'Enter unit definition window (min):', ...
    'Relative activity at the extension temperature (0 to 1):', ...
    'Maximum laps on a circular template:', ...
    ['3'' -> 5'' PROOFREADING exonuclease: none / weak / strong   (N and Y still accepted)' newline ...
    'none   Taq, AccuStart, Stoffel Fragment, Deep Vent (exo-).' newline ...
    'weak   NEB LongAmp Taq, which is a Taq and Deep Vent blend.' newline ...
    'strong Pfu, Vent exo+, Phusion, Q5, KOD.' newline ...
    'Anything but none disables the extension-efficiency weighting in the mispriming screen, ' ...
    'because a mispaired 3'' terminus is excised rather than refused, and the screen falls ' ...
    'back to occupancy alone. weak and strong condition the gates identically; strong adds a ' ...
    'warning about single-stranded yield at this amplicon length.'], ...
    ['Has 5'' -> 3'' exonuclease? (Y/N)' newline ...
    'Taq yes; Stoffel Fragment and the Vent family no. Governs what happens to a ' ...
    'primer bound downstream of an extending polymerase, and so MaxLaps above.']};
dlgtitle = 'Polymerase Throughput';
fieldsize = [1 60; 1 60; 1 60; 1 60; 1 60; 1 60; 1 60; 1 60; 1 60];
definput = {PolyName{polyidx}, PolyAct{polyidx}, '4.5', PolyUnitDNTP{polyidx}, ...
    PolyUnitTime{polyidx}, '1', PolyLaps{polyidx}, PolyProof{polyidx}, Poly5to3{polyidx}};
poly = inputdlg(promptq,dlgtitle,fieldsize,definput);
input.PolymeraseName = poly{1};
input.PolymeraseActivity = str2double(poly{2})*1E-3;
input.ExtensionTime = str2double(poly{3});
input.UnitDNTP = str2double(poly{4})*1E-9;
input.UnitTime = str2double(poly{5});
input.PolymeraseEfficiency = str2double(poly{6});
input.MaxLaps = str2double(poly{7});
input.ProofreadingClass = ProofClass(poly{8}, 'none');
input.ProofreadingExo   = ~strcmp(input.ProofreadingClass, 'none');
input.FivePrimeExo    = YesNo(poly{9}, true);
if input.ProofreadingExo
    warning('CollectaPCRData:ProofreadingPolymerase', ...
        ['%s was declared to carry a %s 3'''' -> 5'''' proofreading exonuclease. The mispriming ' ...
        'screen will run on OCCUPANCY ALONE, because the Huang et al. 1992 extension-efficiency ' ...
        'table describes an enzyme that cannot remove a mispaired 3'''' terminus and this one ' ...
        'can. Expect the screen to reject far more candidates, which is the correct behaviour: ' ...
        'a proofreading enzyme is less discriminating at the priming step. Note also that it ' ...
        'degrades free single-stranded primers over a long run, shortening them and lowering ' ...
        'their Tm, and that nothing in this suite models that.'], ...
        input.PolymeraseName, input.ProofreadingClass);
end
end

%% ____________________________________________________________________________________________________________
function cls = ProofClass(txt, dflt)
% 3' -> 5' exonuclease class: 'none', 'weak' or 'strong', with the older Y and N spellings kept.
%   none    no 3' -> 5' activity. Taq, AccuStart, Stoffel Fragment, Deep Vent (exo-).
%   weak    a Taq and Deep Vent blend such as NEB LongAmp Taq.
%   strong  Pfu, Vent exo+, Phusion, Q5, KOD.
% A BARE Y MAPS TO weak, NOT strong. Y used to mean only that proofreading was declared, and weak
% conditions the gates in exactly the way that flag used to, so an old answer keeps its old
% meaning. The extra yield warning that strong carries is opt-in rather than inherited.
t = lower(strtrim(char(txt)));
switch t
    case {'none', 'n', 'no', '0', 'false'},        cls = 'none';
    case {'weak', 'w', 'y', 'yes', '1', 'true'},   cls = 'weak';
    case {'strong', 's'},                          cls = 'strong';
    otherwise,                                     cls = dflt;
end
end

%% ____________________________________________________________________________________________________________
function tf = YesNo(txt, dflt)
% Y/y/Yes/yes/1/true -> true, N/n/No/no/0/false -> false, anything else -> dflt.
t = lower(strtrim(char(txt)));
if ismember(t, {'y', 'yes', '1', 'true'})
    tf = true;
elseif ismember(t, {'n', 'no', '0', 'false'})
    tf = false;
else
    tf = dflt;
end
end