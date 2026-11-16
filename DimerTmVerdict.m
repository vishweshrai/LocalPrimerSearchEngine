function ok = DimerTmVerdict(rule, tm, dGreg, tmMax, dGok)
%% Dimer Tm Verdict
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 14-11-2026 20:43 IST
%
% Revision history is not kept here. See DimerTmVerdict_<YYYYMMDD>.txt beside this file.
% ____________________________________________________________________________________________________________
% What This Function Does
%
% Decides whether a homodimer or a 3'-anchored self-dimer is acceptable, under the rule chosen by
% DimerRule in PrimerSearchParameters. Primer3 limits both by the melting temperature of the
% dimer (PRIMER_MAX_SELF_ANY_TH and PRIMER_MAX_SELF_END_TH, 47 C by default), so the Tm rule rejects
% a dimer that melts above tmMax. The free-energy rule is the one the engine used before, passed in
% as a verdict because it differs by stage.
%
%   'dG'    the free-energy verdict decides
%   'Tm'    the melting temperature decides. A candidate with no register passes. A candidate whose
%           melting temperature is undefined, which happens when the best register moves between the
%           two probe temperatures, is decided by the free-energy verdict, because a number that
%           cannot be computed must neither pass nor reject on its own
%   'both'  rejected if either rule rejects
% ____________________________________________________________________________________________________________
% Usage
%   ok = DimerTmVerdict(rule, tm, dGreg, tmMax, dGok)
%
%   rule    'dG', 'Tm' or 'both'
%   tm      melting temperature of the register, degrees C, NaN when undefined
%   dGreg   free energy of the register, kcal/mol, 0 when the scan found no register
%   tmMax   highest melting temperature allowed, degrees C
%   dGok    verdict of the free-energy rule, true when it passes
%
% Output
%   ok      true when the candidate passes
% ____________________________________________________________________________________________________________
% References
%   Primer3 manual, PRIMER_MAX_SELF_ANY_TH and PRIMER_MAX_SELF_END_TH. https://primer3.org/manual.html,
%       read 14-11-2026
% ____________________________________________________________________________________________________________

rule = char(rule);
if strcmp(rule, 'dG')
    ok = logical(dGok);
    return
end
if ~isfinite(dGreg)||(dGreg == 0)
    tmOk = true;
elseif isfinite(tm)
    tmOk = (tm <= tmMax);
else
    tmOk = logical(dGok);
end
if strcmp(rule, 'Tm')
    ok = tmOk;
elseif strcmp(rule, 'both')
    ok = tmOk&&logical(dGok);
else
    error('DimerTmVerdict:BadRule', 'rule must be ''dG'', ''Tm'' or ''both''. Received ''%s''.', rule);
end
end
