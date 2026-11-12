function realroots = realcard(sol, f0, fmin)
% REALCARD  Extract real roots between fmin and f0 from Cardano output.
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 05-11-2026 13:30 IST
%
%   realroots = realcard(sol, f0)
%   realroots = realcard(sol, f0, fmin)
%
%   Inputs:
%       sol  : output of cardano(a,b,c,d), size [m n 3]
%       f0   : upper bound for accepted real roots
%       fmin : optional lower bound, default 0. Roots must satisfy fmin < x <= f0.
%
%   Output:
%       realroots : matrix of same size as sol(:,:,1)
%                   containing:
%                      - the largest real root in (fmin, f0], if any,
%                      - NaN otherwise.

% SELECTION RULE, stated because it is a physical assumption and not a mathematical one: when more
% than one root falls in range, the LARGEST is returned. That is correct for these mispriming and
% dimerisation balances only because the primer is in five to six orders of magnitude excess over
% the template, so the free-primer concentration sits just below its total and is also the largest
% root. It is not correct in general. Where the model supplies a real bound, pass it as fmin and the
% choice stops being a heuristic. For the two-site cubic in dGMispTarget each template molecule binds
% at most two primers, so [f] >= [f]0 - 2*[T*]0 exactly, and that bound alone isolates the physical
% root.

% CHANGES
%  1. imag(roots_ij) == 0 was an exact test. Roots that pass through any complex arithmetic pick up
%     an imaginary part of order machine epsilon and were being discarded as complex, which returned
%     NaN and, in dGMispTarget, silently produced a NaN threshold. NaN < anything is false, so the
%     candidate was then accepted with no mispriming check at all. The test is now relative.
%  2. The optional fmin bound was added, see above.
%  3. sol may now be complex-typed throughout (cardano preallocates complex), which is what makes
%     the imaginary-part test meaningful in the first place.

    if nargin < 3 || isempty(fmin)
        fmin = 0;
    end
    if ndims(sol) < 3 || size(sol,3) ~= 3
        error('realcard:BadInput', 'Input must be the output of cardano.m with 3 roots per element.');
    end

    sz = size(sol);
    m = sz(1); n = sz(2);
    realroots = NaN(m,n);

    for i = 1:m
        for j = 1:n
            roots_ij = squeeze(sol(i,j,:));
            % Relative test: a genuinely complex root has an imaginary part comparable to its real
            % part, so any modest relative tolerance separates the two cases cleanly.
            magnitude = max(abs(roots_ij));
            if magnitude == 0 || ~isfinite(magnitude)
                imtol = 0;
            else
                imtol = 1e-8*magnitude;
            end
            % THE BOUNDS ARE READ PER ELEMENT, BECAUSE THE CALLER MAY SUPPLY THEM PER ELEMENT.
            % dGMispTarget forms fmin = max(0, p0 - 2*target) elementwise, so both bounds are
            % arrays whenever p0 or target is. Comparing one element's three roots against the
            % whole array selects the wrong entries and, when the shapes do not conform, raises an
            % indexing error from inside this loop rather than at the call that caused it.
            if isscalar(f0),   f0ij   = f0;   else, f0ij   = f0(i,j);   end
            if isscalar(fmin), fminij = fmin; else, fminij = fmin(i,j); end
            real_roots = real(roots_ij(abs(imag(roots_ij)) <= imtol));
            real_roots = real_roots(real_roots > fminij & real_roots <= f0ij);
            if ~isempty(real_roots)
                realroots(i,j) = max(real_roots);
            end
        end
    end
end
