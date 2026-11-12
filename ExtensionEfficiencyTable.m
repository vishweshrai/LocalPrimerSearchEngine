function out = ExtensionEfficiencyTable(action, value)
%% EXTENSIONEFFICIENCYTABLE  Owns this session's choice of extension-efficiency table
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 07-11-2026 00:02 IST
% ____________________________________________________________________________________________________________
% WHY THIS FILE EXISTS
%
% ExtensionEfficiency already accepts an override through opts.File, and nothing passed one. The
% three call sites that matter are LocalPrimerSearchEngine.m:1877, MisprimerMismatch.m:597 and
% PrimerHybridizationMap.m:1374, all of which call the two-argument form from inside helper
% functions that have no access to the configuration. Threading an options argument down to them
% would change MisprimerMismatch's signature for the third time this month and would still leave
% PrimerHybridizationMap describing a different table from the one the search ran on.
%
% So the choice is held here instead, set once per run by PrimerSearchParameters.m, and consulted
% by ExtensionEfficiency whenever a caller supplies neither opts.Table nor opts.File. This is the
% same pattern AccSaltEquivalent uses for its NNParams struct, and for the same reason.
%
% THE STATE IS HIDDEN, SO IT IS ANNOUNCED IN THREE PLACES. The Parameters sheet carries
% EffTableFile, the engine prints the active table in its startup block, and every call returns the
% source string in info.Source. A run whose table cannot be identified from its own log would be a
% worse defect than the one this fixes.
%
% ____________________________________________________________________________________________________________
% USAGE
%   ExtensionEfficiencyTable('set', path)   choose a CSV
%   ExtensionEfficiencyTable('set', '')     resolve the shipped table, which is the default
%   ExtensionEfficiencyTable('set', 'builtin')  the built-in Huang table, asked for deliberately
%   ExtensionEfficiencyTable('get')         the current choice, '' when nothing has chosen
%   ExtensionEfficiencyTable('clear')       back to the default, which resolves the shipped table
%
% AN EMPTY CHOICE NO LONGER MEANS THE BUILT-IN TABLE. It means nothing has chosen, and
% ExtensionEfficiency then resolves extension_efficiency_rejali2018.csv from its own folder. The
% old behaviour put a standalone caller in a fresh session on Huang's censored bounds without
% anyone asking for them, because this persistent is only ever set by an engine run. 'builtin' is
% the way to ask for that table on purpose, and it is passed through to ExtensionEfficiency
% untouched rather than being checked for existence here.
%
% INPUTS
%   action   'get' (default), 'set' or 'clear'
%   value    for 'set', the full path to an extension-efficiency CSV, '' for the shipped default,
%            or the word 'builtin' for the built-in Huang table
%
% OUTPUTS
%   out      the choice in force after the call
% ____________________________________________________________________________________________________________

persistent chosen
if isempty(chosen), chosen = ''; end
if (nargin < 1)||isempty(action), action = 'get'; end

switch lower(strtrim(char(action)))
    case 'get'
        out = chosen;
    case 'clear'
        chosen = ''; out = chosen;
    case 'set'
        if (nargin < 2)||isempty(value)
            chosen = '';
        else
            v = strtrim(char(value));
            if strcmpi(v, 'builtin')
                chosen = 'builtin';
            elseif ~isfile(v)
                error('ExtensionEfficiencyTable:NoFile', ...
                    ['No extension-efficiency table at %s. Pass a full path to a CSV with ' ...
                    'columns PrimerBase, OpposingBase and Epsilon, '''' to resolve the shipped ' ...
                    'table, or ''builtin'' for the built-in Huang et al. 1992 values.'], v);
            else
                chosen = v;
            end
        end
        out = chosen;
    otherwise
        error('ExtensionEfficiencyTable:BadAction', ...
            'action must be ''get'', ''set'' or ''clear''. Received ''%s''.', char(action));
end
end
