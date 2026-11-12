function [maxnumocc, maxdirep, check] = NucleoRun(sequence, monolim, direplim)
%% NUCLEORUN  Low-complexity screen for a primer candidate
% Low-complexity screen for a primer candidate.
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 05-11-2026 18:49 IST
%
% maxnumocc : length, in nucleotides, of the longest homopolymer run.
% maxdirep  : length, in repeat units, of the longest perfect dinucleotide repeat.
% check     : false when either limit is exceeded.
%
% Defaults follow the Primer3 convention for PRIMER_MAX_POLY_X (Untergasser et al. 2012, NAR 40:e115).
% Pass monolim = 2 to reproduce the behaviour of the previous version, which rejected any run of three.
%
% Two defects in the previous version are corrected here:
%   1. The dinucleotide table was built as a 16x2 char matrix and then indexed as dinucleotides(k), a linear
%      index into column-major storage. Every pattern it produced was a single base, so the second pass simply
%      repeated the mononucleotide test and the dinucleotide screen never executed.
%   2. maxnumocc was reassigned on every loop iteration rather than tracking a running maximum, so the value
%      returned described whichever base happened to be examined last.
% Scanning directly is also cheaper than the regexp pass it replaces, since the candidate is under 30 nt.

if nargin < 2 || isempty(monolim)
    monolim = 4;
end
if nargin < 3 || isempty(direplim)
    direplim = 4;
end

% CASE AND ALPHABET ARE FIXED ON ENTRY, BECAUSE EVERY FAILURE HERE IS A SILENT PASS.
% Both scans below compare raw characters. A uniformly lowercase sequence still finds its runs,
% but a MIXED-CASE one breaks a run at each change of case, so 'AAaa' reports a longest run of two
% and clears a gate that 'AAAA' fails. An N or an IUPAC ambiguity code is counted as an ordinary
% base and can form a run of its own, or silently interrupt a real one. Neither condition
% announces itself: the function returns check = true and the candidate enters the pool on a run
% length that was never measured.
sequence = upper(char(sequence));
if isempty(sequence)
    error('NucleoRun:EmptySequence', ...
        'NucleoRun received an empty sequence. A candidate with no bases has no complexity to screen.');
end
% Written out rather than through ismember, which converts a char row to numeric before
% comparing and is both slower and noisier about it.
offAlphabet = ~((sequence == 'A')|(sequence == 'C')|(sequence == 'G')|(sequence == 'T'));
if any(offAlphabet)
    badChars = unique(sequence(offAlphabet));
    error('NucleoRun:NonACGT', ...
        ['The sequence carries %d character(s) outside {A,C,G,T}, namely ''%s''. The homopolymer ' ...
        'and dinucleotide scans treat any character as an ordinary base, so an ambiguity code ' ...
        'would be counted into a run and the candidate would clear this gate on a run that is ' ...
        'not one.'], nnz(offAlphabet), badChars);
end

% numel rather than strlength: the scans index sequence(k) as a char row, and numel is the length
% that indexing actually uses.
seqlen = numel(sequence);
check = true;

% Longest homopolymer run, in nucleotides.
maxnumocc = 1; runlen = 1;
for k = 2:seqlen
    if sequence(k) == sequence(k-1)
        runlen = runlen + 1;
    else
        maxnumocc = max(maxnumocc, runlen);
        runlen = 1;
    end
end
maxnumocc = max(maxnumocc, runlen);

% Longest perfect dinucleotide repeat, in repeat units. Units of the form XX are skipped, since a homopolymer
% is already counted above and would otherwise be reported twice under two different names.
maxdirep = 0;
for k = 1:(seqlen - 1)
    unit = sequence(k:k+1);
    if unit(1) == unit(2)
        continue
    end
    numocc = 1; p = k + 2;
    while (p + 1 <= seqlen)&&strcmp(sequence(p:p+1), unit)
        numocc = numocc + 1; p = p + 2;
    end
    maxdirep = max(maxdirep, numocc);
end

if (maxnumocc > monolim)||(maxdirep > direplim)
    check = false;
end
end
