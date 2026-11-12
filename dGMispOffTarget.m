function dGmisprimoftarg = dGMispOffTarget(dGptarget, phi, p0, target, Ta)
    %% DGMISPOFFTARGET  Mispriming threshold on a template molecule carrying no cognate site
    % Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
    % 04-11-2026 23:10 IST
    % Threshold DeltaG for a primer binding a site on a molecule that carries NO cognate site for
    % that primer, i.e. mispriming on the opposite strand. Two independent template species, each at
    % concentration target: the correct strand binds with Ka, the off-target strand with the unknown
    % Km being bounded here.
    %
    % DERIVATION, and it matches Limiting Off-Target Priming.pdf line for line
    %   Acceptance rule: [off-target duplex] = phi * [cognate duplex]
    %       target*Km*f/(1+Km*f) = phi*target*Ka*f/(1+Ka*f)
    %   ->  Km = phi*Ka/(1 + Ka*f*(1-phi))                      <- the returned expression
    %   Primer balance: p0 = f + target*Ka*f*(1+phi)/(1+Ka*f)
    %   Clearing (1+Ka*f):
    %       Ka*f^2 + f*[1 + Ka*target*(1+phi) - Ka*p0] - p0 = 0
    %
    % CHANGE: the leading coefficient was Ka.^2. The derivation document writes Ka, B and C in the
    % code already matched that balance exactly, and dimensional analysis settles it on its own,
    % since A*f^2 must carry units of molarity and Ka^2*f^2 does not. The consequence was severe: at
    % 1 uM primer, 1 nM template and 55 C the free primer solved to about 2.5e-14 M instead of
    % 1e-6 M, and the returned threshold to -13.9 instead of -7.1 kcal/mol, which left the
    % off-target mispriming screen effectively disabled.
    %
    % CHANGE: the quadratic root now uses the numerically stable form. With C = -p0 < 0 exactly one
    % root is positive, but (-B + sqrt(B^2-4AC))/(2A) loses precision to cancellation when B > 0.
    R = 8.314E-3/4.184;
    Ta = 273.15 + Ta;
    Ka = exp(-dGptarget./(R*Ta));
    A = Ka;
    B = 1 + Ka.*(1+phi).*target - Ka.*p0;
    C = -p0;
    disc = B.^2 - 4.*A.*C;
    signB = ones(size(B)); signB(B < 0) = -1;
    q = -0.5.*(B + signB.*sqrt(disc));
    fsol = max(q./A, C./q); % A > 0 and C < 0, so one root is positive and one negative
    dGmisprimoftarg = -R.*Ta.*log(phi.*Ka./(1+Ka.*(1-phi).*fsol));
end
