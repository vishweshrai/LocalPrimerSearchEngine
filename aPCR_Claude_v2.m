%% *Main Script for Running aPCR Design and Simulations*
% Written by Vishwesh Rai, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 04-11-2026 23:06 IST
% ____________________________________________________________________________________________________________
% Solves three types of problems:
% 1. Simulates aPCR if primer sequences are already available (asks for plasmid details, primer sequences, and
% reaction conditions),
% 2. Designs primers for a user defined gene sequence (asks for plasmid details, gene details, and reaction conditions)
% 3. Designs primers for a user defined scaffold length (asks for plasmid details, scaffold length, and reaction conditions)
% Annealing temperature is NOT defined by the user, but suggested by the program. The extension temperature is governed
% by the polymerase the user desires to use. On the same lines, the minimum primer melting temperature is not requested
% from the user.
% ____________________________________________________________________________________________________________
% Revision history is NOT kept here. See aPCR_Claude_v2_<YYYYMMDD>.txt beside this file.
% ____________________________________________________________________________________________________________
%% System Compatibility Check
clear; close all; clear functions; clc;
fprintf('Checking system compatibility...\n');
ToolBox = ver;
NeededToolboxes = {'Bioinformatics Toolbox', 'Text Analytics Toolbox'};
NeededLicences = {'Bioinformatics_Toolbox', 'Text_Analytics_Toolbox'};
for t = 1:numel(NeededToolboxes)
    if license('test', NeededLicences{t})
        fprintf('%s is licensed and available.\n', NeededToolboxes{t});
    end
    if ~any(strcmpi({ToolBox.Name}, NeededToolboxes{t}))
        fprintf(['%s is not installed. Open the Add-On Explorer from the Home tab, install it, ' ...
            'restart MATLAB and run this script again.\n'], NeededToolboxes{t});
        return
    end
end
%% Collect User Input
file = matlab.desktop.editor.getActiveFilename;
[filepath, ~, ~] = fileparts(file);
masterpath = filepath; parentDir = '';
while true
    [upperDir, folderName, ~] = fileparts(masterpath);
    if strcmpi(folderName, 'aPCR')
        parentDir = masterpath;
        break
    end
    if (isempty(upperDir)||strcmp(upperDir, masterpath))
        break
    end
    masterpath = upperDir;
end
[input, decider, mode] = CollectaPCRData_Claude(masterpath);
tab = table(); primerseq = {};
%% Name the Annealing Programme
% One string, used by the run summary and by every exported filename. The old script referred to
% tramp in both places without ever assigning it.
if ischar(input.AnnealingTemperature)||isstring(input.AnnealingTemperature)
    tramp = '_Touchdown';
elseif input.AnnealingTemperature == 0
    tramp = '_Programmed';
elseif input.AnnealingTemperature == 1
    tramp = '_Variable';
else
    tramp = ['_Fixed' num2str(input.AnnealingTemperature)];
end
%% Create Logging Directories
waqt = char(datetime('now','Format','yyyyMMddHHmm'));
if exist('F:\', 'dir')
    resFold = 'F:\aPCR_Computations\Simulation_Results\';
else
    resFold = 'D:\aPCR_Computations\Simulation_Results\';
end
res = [resFold, erase(input.PlasmidFileName,'.txt'), '\'];
if isfolder(res) ~= 1
    mkdir(res);
end
%% Decide Whether to Launch Local or Global Search and Design Primers
if ismember(decider, {'Global', 'Local'})
    if ismember(decider, 'Global')
        [tab, srchlog] = GlobalPrimerSearchEngine(input);
    elseif ismember(decider, {'Local'})
        [tab, srchlog] = LocalPrimerSearchEngine(input);
    end
    fprintf('Primer search results will be saved at %s\n', [res, 'Primer_Search\', decider, '\']);
    if isfolder([res, 'Primer_Search\', decider]) ~= 1
        mkdir([res, 'Primer_Search\', decider]);
    end
    % ONE WORKBOOK PER RUN. The path was being rebuilt five times from the same four pieces, and a
    % run that found a pair left two files that had to be opened side by side to read one result.
    % Everything now goes into a single book: the three logs, the ranked candidates when there are
    % any, and the settings the run used. Candidate_Primers is ABSENT rather than empty when no
    % pair qualified, because an empty sheet reads like a failure of the writer.
    xlsFile = [res, 'Primer_Search\', decider, '\', num2str(input.ScaffoldLength), ...
        '_nt_', waqt, '_Log.xlsx'];
    % BOTH ENGINES RETURN REAL FORWARD AND REVERSE LOGS, so neither is skipped any more. The
    % guard here used to write them for the Local path alone, because
    % GlobalPrimerSearchEngine_Claude evaluated every reverse primer inside the forward loop and
    % returned Reverse_Log empty. GlobalPrimerSearchEngine screens the two families in separate
    % passes and both logs are complete. Without them an empty Pair_Log says only that no pair
    % qualified, and not which stage ate the candidates, which is the one thing worth knowing.
    % A FREE-POSITION SCAN OVERFLOWS A SHEET'S USEFUL SIZE. It logs one row per candidate ever
    % built, of order 36000 forward and 65000 reverse on a 7249 nt circle, and writetable to xlsx
    % is slow enough at that row count to dominate the run. Past LogRowCap the log goes to a CSV
    % beside the workbook instead, which carries the same columns and opens in seconds.
    LogRowCap = 100000;
    for lg = {'Forward_Log', 'Reverse_Log'}
        LogTab = srchlog.(lg{1});
        if isempty(LogTab), continue; end
        if height(LogTab) <= LogRowCap
            writetable(LogTab, xlsFile, 'Sheet', lg{1});
        else
            csvFile = strrep(xlsFile, '.xlsx', ['_', lg{1}, '.csv']);
            writetable(LogTab, csvFile);
            fprintf('%s has %d rows, past the %d-row sheet cap; written to %s\n', ...
                lg{1}, height(LogTab), LogRowCap, csvFile);
        end
    end
    writetable(srchlog.Pair_Log, xlsFile, 'Sheet', 'Pair_Log');
    if ~isempty(tab)
        writetable(tab, xlsFile, 'Sheet', 'Candidate_Primers');
    end
    % GlobalPrimerSearchEngine_Claude does not build a Parameters table, so this is guarded rather
    % than assumed. A missing sheet is better than a crash on the path that still works.
    if isstruct(srchlog)&&isfield(srchlog, 'Parameters')&&~isempty(srchlog.Parameters)
        writetable(srchlog.Parameters, xlsFile, 'Sheet', 'Parameters');
    end
    fprintf('Wrote %s\n', xlsFile);
end
%% Run Simulations with Chosen Primers
if ismember(decider, {'Global', 'Local'})
    if ~isempty(tab)
        promptcont = {'Do you want to simulate aPCR for your chosen primers? [Y/N]'};
        dlgtitle = 'Simulations';
        fieldsize = [1 50];
        definput = {'N'};
        conti = inputdlg(promptcont,dlgtitle,fieldsize,definput); contichoice = conti{1};
        if ismember(contichoice, {'N', 'n', 'No', 'no'})
            fprintf('You chose to NOT continue with simulations!\n');
            return
        elseif ismember(contichoice, {'Y', 'y', 'Yes', 'yes'})
            prompti = {'Input the chosen forward primer sequence:', ...
                'Input the chosen reverse primer sequence:'};
            dlgtitle = 'Confirm Chosen Sequences for Simulation';
            fieldsize = [1 50; 1 50];
            definput = {tab.F_Primer{1}, tab.R_Primer{1}};
            primchoice = inputdlg(prompti,dlgtitle,fieldsize,definput);
            primerseq = {primchoice{1}, primchoice{2}};
            input.ForwardPrimer = primerseq{1}; input.ReversePrimer = primerseq{2};
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
            % SCAFFOLDLENGTH FOLLOWS THE COORDINATES IT DESCRIBES.
            % The two lines above replace ScaffoldStart and ScaffoldEnd with the coordinates of the
            % primers the search chose, which are not in general the scaffold the user entered.
            % Leaving ScaffoldLength at its old value puts the struct in contradiction with itself:
            % the sequence written out below is cut from the new coordinates while the file is
            % named with the old length, and aPCRCycleAlgebra_Claude receives a length that does
            % not match the Start and End beside it. The circular form is the one
            % CollectaPCRData_Claude.m uses on the Yield_Test route, so the two paths into the
            % simulator now mean the same thing by the same arithmetic.
            if input.ScaffoldStart > input.ScaffoldEnd
                input.ScaffoldLength = input.ScaffoldEnd + strlength(input.PlasmidSequence) - input.ScaffoldStart + 1;
            else
                input.ScaffoldLength = input.ScaffoldEnd - input.ScaffoldStart + 1;
            end
        end        
    else
        primerseq = {};
    end
elseif strcmp(decider, 'Yield_Test')
    primerseq{1} = input.ForwardPrimer; primerseq{2} = input.ReversePrimer;
end
%% Equilibrium and Kinetic Simulations
if ~isempty(primerseq)
    if strcmp(input.ForwardPrimer, '')
        input.ForwardPrimer = primerseq{1}; input.ReversePrimer = primerseq{2};
    end
    promptz = {'Choose approximation mode for simulation: Equilibrium or Kinetic (E/K)'};
    dlgtitle = 'Equilibrium or Kinetic Approximation';
    fieldsize = [1 50];
    definput = {'E'};
    simdata = inputdlg(promptz,dlgtitle,fieldsize,definput); simchoice = simdata{1};

    simfig2 = []; simresources = struct(); simother = struct([]);
    if ismember(simchoice, {'E', 'e'})
        sim = '_Equilibrium';
        [simfig, simfig2, simamplicon, simother, simresources, simpool] = aPCRCycleAlgebra_Claude(input);
    elseif ismember(simchoice, {'K', 'k'})
        sim = '_Kinetic';
        if exist('kineticapcr_Claude', 'file') ~= 2
            fprintf('kineticapcr_Claude.m is not on the path. Rebuild it or choose the equilibrium mode.\n');
            return
        end
        [simfig, simfig2, simamplicon, simother, simresources, simpool] = kineticapcr_Claude(input, tramp);
    else
        fprintf('Approximation mode must be E or K.\n');
        return
    end
    if input.ScaffoldStart <= input.ScaffoldEnd
        scaffoldSeq = input.PlasmidSequence(input.ScaffoldStart:input.ScaffoldEnd);
    else
        scaffoldSeq = [input.PlasmidSequence(input.ScaffoldStart:end), input.PlasmidSequence(1:input.ScaffoldEnd)];
    end

    decider = 'Yield_Test';
    fprintf('Simulation results will be saved at %s.\n', [res, decider, '\']);
    ScaffoldIDx = [num2str(input.ScaffoldStart) '-' num2str(input.ScaffoldEnd)];
    if isfolder([res, decider]) ~= 1
        mkdir([res, decider]);
    end
    stem = [res, decider, '\', num2str(input.ScaffoldLength), '_nt_', waqt, sim, tramp];

    %% Write the Scaffold Sequence
    fid = fopen([res, decider, '\', erase(input.PlasmidFileName,'.txt'), ...
        erase(['_' mode '_'], '_None'), num2str(input.ScaffoldLength), '_nt_', 'IDx_', ScaffoldIDx 'scaffoldSeq.txt'], "w+");
    if fid < 0
        fprintf('Could not open the scaffold sequence file for writing.\n');
    else
        fprintf(fid, '%s', scaffoldSeq);
        fclose(fid);
    end

    %% Export the Figures
    exportgraphics(simfig, [stem, 'T.svg'], 'ContentType', 'vector');
    if ~isempty(simfig2)&&isgraphics(simfig2)
        exportgraphics(simfig2, [stem, 'T_Byproducts.svg'], 'ContentType', 'vector');
    end
    %% Write the Run Summary
    % SimParSum prints the parameter block, the outcome and the resource diagnostics, and writes
    % the workbook when it is handed a path. The kinetic branch has no resources struct yet, so it
    % gets the shorter call.
    if exist('SimParSum', 'file') == 2
        book = [stem, 'T_Run_Summary.xlsx'];
        if ismember(simchoice, {'E', 'e'})
            summarytxt = SimParSum(input, tramp, sim, simamplicon, simresources, simother, book);
        else
            summarytxt = SimParSum(input, tramp, sim, simamplicon);
        end
        fid = fopen([stem, 'T_Parameters.txt'], "w+");
        if fid < 0
            fprintf('Could not open the parameter log for writing.\n');
        else
            fprintf(fid, '%s', summarytxt);
            fclose(fid);
        end
    else
        fprintf('SimParSum.m is not on the path; no run summary written.\n');
    end
    fprintf('Mock aPCR run complete!\n');
else
    fprintf('%s search returned no sequences! Relax the primer search criteria, and re-run the script!\n', decider);
end
