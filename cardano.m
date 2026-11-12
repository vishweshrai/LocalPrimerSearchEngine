function sol = cardano(a,b,c,d)
% CARDANO  Solve cubic equations ax^3 + bx^2 + cx + d = 0 using Cardano's method.
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 04-11-2026 23:10 IST
% Works for scalar, vector, or 2-D matrix inputs of identical size.
% Returns sol(:,:,k) where k = 1,2,3 corresponds to each root.

% CHANGES
%  1. The degeneracy tolerance was the absolute constant 1E-14, assigned to a variable named eps
%     which also shadowed the builtin. For the cubics this suite actually solves, a is of order
%     Kp^2 ~ 1e21 while d = -p0 ~ 1e-6, so Delta lands between 1e-42 and 1e-55 across the whole
%     parameter range and |Delta| <= 1e-14 was true on EVERY call. Cases 1 and 2 were unreachable
%     and the multiple-root branch ran every time. The tolerance is now relative to the two terms
%     whose difference forms Delta, which is the only scale-free way to ask whether they cancel.
%
%  2. Subcase 3B carried the wrong sign on the depressed-to-original shift. For t^3 + p*t + q with
%     Delta = 0 the roots are t = cbrt(q/2) twice and t = -2*cbrt(q/2) once, and
%     x = t - b/(3a) = t + base. The code computed r1 - base and -2*r1 - base. Combined with
%     defect 1 this was the operative error. At dGp = -16 kcal/mol, 1 uM primer, 1 nM template and
%     55 C the free-primer root came out at 3.330e-7 M against a true 9.990e-7 M, which moved the
%     dGMispTarget threshold from -7.06 to -7.77 kcal/mol, i.e. 0.7 kcal/mol too permissive on
%     every sense-strand mispriming test.
%
%  3. Case 3 indexed as a(idx(pz)), where idx is a full-size logical mask and pz is a shorter
%     logical vector over the already-selected elements. idx(pz) logically indexes idx by pz and
%     the result is then used to index a, which is meaningless. It is harmless at 1x1, which is why
%     nothing had broken. Full-size masks are now carried throughout and the three roots are built
%     as separate arrays and concatenated, so no page indexing is needed at all.
%
%  4. The acos argument in Case 1 is clamped to [-1, 1]. It equals +/-1 exactly at Delta = 0, and
%     roundoff pushes it outside, at which point acos returns complex and all three supposedly real
%     roots come back complex, so realcard returns NaN. With the tolerance fixed, Case 1 is now
%     reachable for these cubics and this matters.
%
%  5. sol is preallocated complex. Case 2 writes complex conjugate roots into it, and realcard's
%     imag() test depends on those imaginary parts surviving.
%
%  6. a == 0 is rejected instead of producing Inf and NaN silently.
%
%  7. The original header claimed N-D support, which sol(:,:,k) does not provide. The header now
%     states what the code does.

    % --- Input validation ---
    if ~isequal(size(a), size(b), size(c), size(d))
        error('cardano:SizeMismatch', 'Inputs a, b, c, d must have identical dimensions.');
    end
    if ndims(a) > 2 %#ok<ISMAT>
        error('cardano:NotTwoDimensional', 'Inputs must be scalar, vector or 2-D.');
    end
    if any(a(:) == 0)
        error('cardano:DegenerateLeadingCoefficient', ...
            'a == 0 for at least one element; that equation is not cubic.');
    end

    % --- Precompute parameters ---
    p = (3.*a.*c - b.^2) ./ (3.*a.^2);
    q = (2.*b.^3 - 9.*a.*b.*c + 27.*a.^2.*d) ./ (27.*a.^3);
    term1 = (q./2).^2;
    term2 = (p./3).^3;
    Delta = term1 + term2;
    base = -b./(3.*a);

    % Relative degeneracy tolerance. Delta is the difference of two quantities of size `scale`, and
    % each carries a few units in the last place from the formation of p and q, so anything below a
    % few tens of machine epsilon times `scale` is numerically indistinguishable from zero.
    scale = max(abs(term1), abs(term2));
    tol = 64*eps(class(Delta)).*scale;

    % --- Initialize output arrays ---
    sz = size(a);
    r1 = complex(nan(sz)); r2 = complex(nan(sz)); r3 = complex(nan(sz));

    % ======================
    % Case 1: Delta < -tol (three distinct real roots)
    % ======================
    idx = Delta < -tol;
    if any(idx(:))
        arg = (3.*q(idx)./(2.*p(idx))).*sqrt(-3./p(idx));
        arg = min(1, max(-1, arg)); % see change note 4
        theta = (1/3).*acos(arg);
        r = 2.*sqrt(-p(idx)./3);
        bb = base(idx);
        r1(idx) = bb + r.*cos(theta - 2*pi*(0/3));
        r2(idx) = bb + r.*cos(theta - 2*pi*(1/3));
        r3(idx) = bb + r.*cos(theta - 2*pi*(2/3));
    end

    % ======================
    % Case 2: Delta > tol (one real + two complex conjugates)
    % ======================
    idx = Delta > tol;
    if any(idx(:))
        A = -q(idx)./2 + sqrt(Delta(idx));
        B = -q(idx)./2 - sqrt(Delta(idx));
        w = -1/2 + 1i*sqrt(3)/2;
        bb = base(idx);
        cA = nthroot(A,3); cB = nthroot(B,3);
        r1(idx) = bb + cA + cB;
        r2(idx) = bb + w.*cA + w^2.*cB;
        r3(idx) = bb + w^2.*cA + w.*cB;
    end

    % ======================
    % Case 3: |Delta| <= tol (multiple roots)
    % ======================
    idx = (Delta >= -tol) & (Delta <= tol);

    % --- Subcase 3A: triple root. Both terms vanish only when p = q = 0. ---
    idx3a = idx & (scale <= realmin(class(Delta)));
    if any(idx3a(:))
        r1(idx3a) = base(idx3a);
        r2(idx3a) = base(idx3a);
        r3(idx3a) = base(idx3a);
    end

    % --- Subcase 3B: one single + one double root ---
    idx3b = idx & ~idx3a;
    if any(idx3b(:))
        rr = nthroot(q(idx3b)./2, 3);
        bb = base(idx3b);
        r1(idx3b) = rr + bb;        % double root
        r2(idx3b) = rr + bb;        % double root
        r3(idx3b) = -2.*rr + bb;    % single root
    end

    sol = cat(3, r1, r2, r3);
end
