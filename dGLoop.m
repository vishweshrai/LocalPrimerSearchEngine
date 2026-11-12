function dGloop = dGLoop(looptype, n, tann, N)
%% DGLOOP  Loop-initiation free energy for a DNA bulge, internal loop or hairpin
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 04-11-2026 16:58 IST
% ____________________________________________________________________________________________________________
% WHAT IT RETURNS
%
% DeltaG (kcal/mol, 1 M NaCl) for the initiation of a loop holding n unpaired nucleotides, read
% from Table 4 of SantaLucia and Hicks 2004. looptype selects the column: 'bulge', 'internal' or
% 'hairpin'. The loop size n counts unpaired nucleotides on BOTH strands for an internal loop,
% which is how Table 4 is indexed.
%
% WHERE THE NUMBERS COME FROM
%
% The three columns are read from loop_length_params.csv through NNParams, as N.Loop.Internal,
% N.Loop.Bulge and N.Loop.Hairpin, indexed by N.Loop.N. Nothing is held here. A caller that
% already has the parameter structure passes it as the fourth argument, which costs nothing; a
% caller without one leaves it out and this function fetches it from the NNParams cache.
%
% LENGTHS THE TABLE DOES NOT LIST
%
% Table 4 lists 1 to 10, then 12, 14, 16, 18, 20, 25 and 30. Sizes inside that span but absent
% from it (11, 13, 15, 17, 19, 21 to 24, 26 to 29) are reached by linear interpolation between the
% two bracketing entries, which the review does not state and which is used here as a practical
% fallback. Beyond 30 the Jacobson-Stockmayer entropy extrapolation of Equation 7 applies,
%
%       DeltaG37(n) = DeltaG37(x) + c * R * 310.15 * ln(n/x)
%
% with x the longest tabulated size and c the coefficient carried in N.Loop.JSCoeff. The review
% prefers 2.44, measured from DNA kinetics, over the theoretical 1.75.
%
% SIZES A MOTIF DOES NOT ADMIT
%
% An internal loop of 1 is impossible and an internal loop of 2 is a single mismatch, which takes
% the mismatch nearest-neighbour parameters of Table 2 rather than this function (Table 4, footnote
% f). Hairpin loops shorter than 3 are sterically prohibited. Those cells are empty in the CSV, and
% a request for one raises.
%
% TEMPERATURE
%
% Loop enthalpies are zero by construction (Table 4, footnote a), so the motif is purely entropic:
% DeltaS = -DeltaG37*1000/310.15 and DeltaG(T) = DeltaG37*(T/310.15) with T in Kelvin. Pass tann in
% degrees Celsius to get the value at the annealing temperature; omit it or pass [] for the raw
% 37 degree, 1 M NaCl literature value.
%
% SALT
%
% No salt dependence is published for the loop-length term itself. The published salt corrections
% apply to the terminal mismatches and Watson-Crick stacks that flank a loop, which are handled by
% dGTerminalMismatch and dGMismatch separately.
%
% REFERENCE
%   SantaLucia J Jr, Hicks D. The thermodynamics of DNA structural motifs. Annu Rev Biophys Biomol
%       Struct 2004;33:415-40. doi:10.1146/annurev.biophys.32.110601.141800
% ____________________________________________________________________________________________________________

if (nargin < 4)||isempty(N)
    N = NNParams(struct('Verbose', false));
end
if ~isfield(N, 'Loop')||isempty(N.Loop.N)
    error('dGLoop:NoLoopTable', ...
        ['The parameter structure carries no loop table. loop_length_params.csv must sit beside ' ...
        'NNParams.m; without it there are no Table 4 values to read.']);
end

sizes = N.Loop.N(:).';
switch lower(looptype)
    case 'bulge'
        tab = N.Loop.Bulge(:).';
    case 'internal'
        tab = N.Loop.Internal(:).';
        if n == 2
            error('dGLoop:SingleMismatch', ...
                ['Internal loops of 2 are single mismatches. Use the mismatch nearest-neighbour ' ...
                'parameters (Table 2, dGMismatch.m), not Table 4.']);
        end
    case 'hairpin'
        tab = N.Loop.Hairpin(:).';
    otherwise
        error('dGLoop:BadType', 'looptype must be "bulge", "internal" or "hairpin".');
end
if n < 1
    error('dGLoop:BadLength', 'Loop length must be at least 1.');
end

if any(sizes == n)
    dGloop = tab(sizes == n);
    if isnan(dGloop)
        error('dGLoop:NotAdmitted', ...
            'Loop length %d is not allowed for %s loops (Table 4, SantaLucia and Hicks 2004).', ...
            n, looptype);
    end
elseif n < sizes(end)
    lo = max(sizes(sizes < n)); hi = min(sizes(sizes > n));
    dGlo = tab(sizes == lo); dGhi = tab(sizes == hi);
    dGloop = dGlo + (dGhi - dGlo)*(n - lo)/(hi - lo);
else
    x = sizes(end); dGx = tab(end);
    dGloop = dGx + N.Loop.JSCoeff*N.Rgas*310.15*log(n/x);
end

if (nargin >= 3)&&~isempty(tann)
    % A STRUCT HERE MEANS THE PARAMETER STRUCTURE LANDED IN THE TEMPERATURE SLOT.
    % Calling dGLoop(type, n, N) instead of dGLoop(type, n, [], N) is the easy slip to make once
    % this function takes four arguments, and without this test it surfaces as an arithmetic
    % error from the line below rather than as the argument error it is.
    if ~isnumeric(tann)||~isscalar(tann)
        error('dGLoop:BadTemperature', ...
            ['tann must be a scalar temperature in degrees Celsius, or [] for the 37 degree ' ...
            'value. To pass the parameter structure, use dGLoop(looptype, n, [], N).']);
    end
    dGloop = dGloop*((tann + 273.15)/310.15);
end
end
