function dGmisptarget = dGMispTarget(dGptarget, phi, p0, target, Ta)
    %% DGMISPTARGET  Mispriming threshold on a template molecule that also carries the cognate site
    % Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
    % 05-11-2026 13:30 IST
    % Threshold DeltaG for a primer binding a second, non-cognate site on the SAME template molecule
    % that carries its cognate site. Four occupancy states exist (empty, cognate only, misprime only,
    % both), which is why this solves a cubic where dGMispOffTarget solves a quadratic.
    %
    % VERIFIED against Limiting Off-Target Priming.pdf, section on mispriming of both strands by the
    % forward primer. With
    %       phi = ([T*fm] + [fpT*fm]) / [fpT*]
    % the occupancy sum collapses to [T*] = [T*]0/(1 + Kp[f](1+phi)) and the acceptance rule becomes
    %       Km = phi*Kp/(1 + Kp[f])
    % EXACTLY, not approximately, because the doubly-occupied state sits in the numerator of phi.
    % Substituting into [f]0 = [f] + [fpT*] + [T*fm] + 2[fpT*fm] gives
    %       Kp^2*(1+phi)*[f]^3
    %     + (Kp*(2+phi) + Kp^2*[T*]0*(1+2phi) - Kp^2*[f]0*(1+phi))*[f]^2
    %     + (1 + Kp*[T*]0*(1+phi) - Kp*[f]0*(2+phi))*[f]
    %     - [f]0 = 0
    % which is A, B, C, D below, term for term. No algebraic change was needed here.
    %
    % CHANGES
    %  1. The physical lower bound on the free primer is now passed to realcard. Each template
    %     molecule binds at most two primers, so [f] >= [f]0 - 2*[T*]0 exactly. That isolates the
    %     physical root outright instead of relying on realcard's "take the largest" heuristic,
    %     which is only valid while the primer is in large excess.
    %  2. Two Newton steps polish the selected root against the original cubic. Delta for these
    %     coefficients runs at 1e-42 to 1e-55 while a ~ 1e21 and d ~ 1e-6, so the closed-form root
    %     carries five to six significant figures rather than sixteen. Polishing costs nothing and
    %     is applied only to the selected root; the two spurious near-degenerate roots must not be
    %     polished, since f' is near zero there and Newton walks away from them.
    %  3. realcard returning NaN, empty or a non-physical value now raises. Previously that produced
    %     a NaN threshold, and NaN < anything is false, so the candidate passed with no mispriming
    %     check at all.
    R = 8.314E-3/4.184;
    Ta = 273.15 + Ta;
    Kp = exp(-dGptarget./(R*Ta));
    A = Kp.^2.*(1+phi);
    B = Kp.*(2+phi) + Kp.^2.*target.*(1+2*phi) - Kp.^2.*p0.*(1+phi);
    C = 1 + Kp.*target.*(1+phi) - Kp.*p0.*(2+phi);
    D = -p0.*ones(size(A));

    fmin = max(0, p0 - 2.*target); % capacity bound: at most two primers bound per template molecule
    sol = realcard(cardano(A,B,C,D), p0, fmin);

    if isempty(sol)||any(~isfinite(sol(:)))||any(sol(:) <= 0)
        error('dGMispTarget:NoPhysicalRoot', ...
            ['No positive finite free-primer root in (%.4e, %.4e] for the two-site cubic ' ...
            '(Kp = %.3e, p0 = %.3e, target = %.3e). Returning NaN here would silently pass the ' ...
            'candidate, so this is raised instead.'], fmin(1), p0(1), Kp(1), p0(1), target(1));
    end

    % Newton polish on the original cubic, selected root only.
    for k = 1:2
        fval = ((A.*sol + B).*sol + C).*sol + D;
        fder = (3.*A.*sol + 2.*B).*sol + C;
        step = fval./fder;
        step(~isfinite(step)) = 0;
        sol = sol - step;
    end

    % THE ROOT THAT IS RETURNED IS THE ROOT THAT IS CHECKED.
    % The guard above tests what realcard selected. What reaches the logarithm below is the
    % POLISHED value, and two Newton steps from a root already good to five significant figures
    % will not move it out of range in practice. Testing only the pre-polish value still leaves the
    % returned quantity unvalidated, and a sol below -1/Kp would turn the argument of the logarithm
    % negative and hand back a COMPLEX threshold, which every downstream comparison treats as false
    % in exactly the way the original NaN did.
    bad = ~isfinite(sol) | (sol <= 0) | ((1 + Kp.*sol) <= 0);
    if any(bad(:))
        k = find(bad, 1);
        error('dGMispTarget:PolishLeftRange', ...
            ['Two Newton steps moved the free-primer root out of its physical range at element ' ...
            '%d: %s, with Kp = %.3e and p0 = %.3e. The pre-polish root was accepted, so the ' ...
            'cubic itself is sound and the polish is what failed.'], ...
            k, mat2str(sol(k), 6), Kp(min(k, numel(Kp))), p0(min(k, numel(p0))));
    end

    dGmisptarget = -R*Ta.*log(phi.*Kp./(1+Kp.*sol));
end
