function out = ScaffoldCollapseCheck(mode, varargin)
%% Amplicon 3'-OH Foldback and Self-Priming Screen
% Written for the aPCR suite, Laboratory of Nanobiology, IISc Bangalore, INDIA
% 12-11-2026 02:16 IST
% ____________________________________________________________________________________________________________
% Physics
%
% In asymmetric (LATE) PCR the accumulating single strand ends, at its 3' terminus, where the
% opposite primer's 5' end sits. That 3'-OH can fold back onto the strand's own body. When the fold
% leaves the terminal base paired and leaves template beyond the stem, polymerase extends the
% 3'-OH and the strand lengthens in an unintended way.
%
% Each fold is a stem plus a hairpin loop (SantaLucia and Hicks 2004). The stem is scored by
% nearest neighbours through dGMismatch, with no bimolecular initiation. The loop takes the Table 4
% length term from dGLoop up to 30 nt and the polymer increment described under Loop Model beyond
% that. A loop of 4 or more nt adds the terminal mismatch of its first and last bases (Equations 9 and 10). A
% triloop takes no terminal mismatch, takes a penalty on an A-T closing pair, and a loop of 3 or 4
% nt adds the sequence bonus where one is tabulated (Equation 8). The hairpin routine of the
% mispriming engine, LooperMismatch, applies the same terms.
%
% No dangling-end term is applied, in line with the decision recorded in dGDangle.
%
% K = exp(-dGfold/RT) is the fold equilibrium constant against the unfolded coil. The extension
% efficiency eps comes from the 3'-terminal three bases (Rejali 2018) and is a kinetic factor on the
% folded end. dGeff = dGfold - RT*ln(eps) ranks the folds and sets the pass or fail decision.
%
% Folds of one 3' tail compete for the same bases and exclude each other. Their partition function
% is Z = 1 + sum(K), the fraction of tails folded is sum(K)/Z, and the probability of an extendable
% collapse is pExt = sum(K*eps)/Z.
%
% A fold counts only when template lies beyond its stem: at least one unpaired base must sit 5' of
% the stem inside the active amplicon.
%
% Scan range. The window start r0 runs to the last position that keeps the window inside the strand,
% and a tail base takes part only if pairing it leaves a loop of at least 3 nt. The scan therefore
% reaches the tight hairpins at the very 3' end, whose loop is shorter than about 2*(W - stem) + 3 nt.
% ____________________________________________________________________________________________________________
% Loop Model
%
% Hairpin loops longer than 30 nt take the worm-like-chain increment of Aksel et al. 2024 and Dunn
% et al. 2015. The mean squared end-to-end distance of a chain of contour length Lc = n*LBase is
%
%       d2 = 2*Lp*Lc*(1 - (Lp/Lc)*(1 - exp(-Lc/Lp)))
%
% with LBase = 0.6 nm per nt and Lp = 0.9 nm (the expression is Rivetti 1998 equation 5, the values
% are those of the Aksel et al. supplement; Rivetti measured 1.3 nm, which moves the fold Tm by about
% 0.1 degrees C). The effective concentration
% of one chain end around the other is j(n) = (1/NA)*(3/(2*pi*d2))^(3/2), in mol/L. The loop
% penalty for n > 30 is
%
%       dG(n) = dG_table(30) + R*T*(ln j(30) - ln j(n))
%
% so the measured table is kept where it applies. The engine form, dG_table(30) + 2.44*R*T37*ln(n/30)
% from dGLoop, is carried beside it in dGfoldEng, tmEng and pExtEng for comparison and drives nothing.
% ____________________________________________________________________________________________________________
% Persistence and Tiers (log-only)
%
% A hairpin is unimolecular, so its melting temperature has no concentration dependence. tmFold comes
% from dGfold at tann and tann + 10 degrees C (two-point van't Hoff). A fold persists when
% tmFold >= tann - HpTmMargin. Tier of one fold:
%   REJECT     persists and eps = 1, the last three 3' bases are paired
%   FLAG       persists and FlagEps <= eps < 1, a 3' mismatch limits extension
%   TRANSIENT  does not persist and eps >= FlagEps, extendable but melted at tann
%   NONE       otherwise
% The tier of a strand is the highest tier among its folds with template beyond the stem. jOverC
% compares the effective concentration of the loop with the strand concentration; above 1 the
% hairpin outcompetes strand-strand binding.
% ____________________________________________________________________________________________________________
% Bimolecular Channel
%
% The 3'-OH can also pair with the same stem on a second copy of the strand. That complex has no loop
% and owes helix initiation, so its free energy dGbi is the stem scored through dGMismatch with the
% initiation and terminal terms applied, which is the duplex BestSubRegister already scores between
% two molecules. With [S] the concentration of partner strands, referred to the 1 M standard state,
%
%       Kb = [S]*exp(-dGbi/RT)          Ku = exp(-dGfold/RT)
%
% Ku carries the loop term. With the loop anchored at its 30 nt table value, a stem of free energy
% dGstem gives Ku = exp(-dGstem/RT)*exp(-dGtab30/RT)*j(n)/j(30), and the same stem on a second strand
% gives Kb = [S]*exp(-dGinit/RT)*exp(-dGstem/RT), where dGinit is the initiation and terminal term that
% dGbi carries over the stem. Apart from the terminal mismatch of the loop,
%
%       Kb/Ku = ([S]/j(n))*exp(-dGinit/RT)*j(30)*exp(dGtab30/RT)
%
% The last two factors come from the table and equal about 88 at 58 C. The bimolecular state therefore
% overtakes foldback once j(n) falls below 88*[S]*exp(-dGinit/RT). At [S] = 1 uM and dGinit between 1.7
% and 2.05 kcal/mol that is a loop of about 1800 to 2500 nt.
%
% The states of one 3' tail exclude each other: Z = 1 + sum(Ku) + sum(Kb). The extendable fraction is
% sum(eps*(Ku + Kb))/Z, and the share carried by the bimolecular states is
% sum(eps*Kb)/sum(eps*(Ku + Kb)). Partner strands are taken at their ceiling concentration and binding
% does not deplete them, so the share is an upper bound.
%
% A bimolecular complex persists when its melting temperature reaches tann - HpTmMargin. The
% temperature comes from BimolecularTm, Tm = dH/(dS + R*ln(ceff)), with dG at tann and tann + 10 C.
% Tails and sites are equimolar, so the half-bound point is K*[S]/2 = 1 and ceff = [S]/2. Polymerase
% extends the 3'-OH along the partner's body, so the template rule is the same as for foldback: at
% least one unpaired base must sit 5' of the stem inside the active amplicon.
%
% [S] is the ceiling of the free strand that carries the stem. For the plus strand it is the excess
% primer concentration. For the minus strand it is the limiting primer concentration, because the
% minus strand cannot exceed the amount of limiting primer extended.
% ____________________________________________________________________________________________________________
% Enzyme Channels
%
% Two channels are evaluated for every strand, because the end an enzyme will extend depends on it.
%   Taq   has no 3' to 5' exonuclease. A site counts when at most endTol bases of the tail hang unpaired
%         beyond the stem, and its efficiency is eps from the Rejali table.
%   Exo   a proofreading exonuclease trims an unpaired 3' flap back to the outermost paired base and
%         excises a mispaired terminus. Efficiency is 1 and a flap of up to ExoEndTol bases counts.
% pExtTaq and pExtExo are the extendable fractions in each channel, over the same partition function.
% A blend in which a fraction phi of the 3' end engagements is served by the proofreading enzyme gives
%
%       pExtBlend = (1 - phi)*pExtTaq + phi*pExtExo
%
% which is exact when every engagement is served by one enzyme molecule. NEB LongAmp Taq is a blend of
% Taq and Deep Vent DNA polymerases whose proportions are not published, and NEB reports a fidelity
% two-fold higher than Taq alone, which does not fix phi. ExoBlend is therefore NaN unless the caller
% supplies a value, and pExtBlend stays empty until then. The strand tier tierAll takes the highest rank
% over both channels and both topologies, so a design that fails under either enzyme is visible.
% The headline fields tier, pExt, pExtEng and eps describe foldback in the declared channel: Taq when
% useEps is true and Exo when it is false.
% ____________________________________________________________________________________________________________
% Usage
%   scanRes = ScaffoldCollapseCheck('scan', maxScaffSeq, tailWin, minRun, endTol, monv, biv, dNTP, tann, useEps, opts)
%   evalRes = ScaffoldCollapseCheck('eval', scanRes, minR0, dGCollapseMax)
%
%   maxScaffSeq  the longest admissible strand, 5' to 3', whose 3' end carries the OH
%   tailWin      length of the 3'-terminal window searched for stems, nt
%   minRun       shortest contiguous Watson-Crick run that counts as a stem
%   endTol       unpaired bases tolerated at the 3' end before the fold is a flap
%   monv, biv, dNTP  monovalent, divalent and dNTP concentrations, M
%   tann         annealing temperature, degrees C
%   useEps       true applies the extension efficiency, false takes eps = 1
%   opts         optional struct, every field optional: StrandConc (M, default 1E-6), LBase (nm, 0.6),
%                LPers (nm, 0.9), HpTmMargin (degrees C, 5.0), FlagEps (1E-3), Bimolecular (true),
%                ExoBlend (phi in [0 1], default NaN), ExoEndTol (nt, default 15, clamped to the window
%                less minRun)
%   minR0        index on maxScaffSeq where the active amplicon begins
%   dGCollapseMax  folds with dGeff at or below this value fail, kcal/mol
%
% Output (eval) dGfoldMin, dGeffMin, dGensEff, pFold, pExt, loopLen, targetOffset, eps, epsPair,
%   align, nSites, pass, plus the log-only fields dGfoldEng, tmFold, tmEng, jOverC, pExtEng, tier,
%   tierRank (0 none, 1 transient, 2 flag, 3 reject), pExtTot (foldback plus bimolecular, declared
%   channel), biShare (bimolecular share of pExtTot), pExtTaq, pExtExo, pExtBlend, tmBi (melting
%   temperature of the strongest bimolecular stem) and tierAll, tierAllRank (worst tier over both
%   channels and both topologies)
% ____________________________________________________________________________________________________________
% References
%   SantaLucia J Jr, Hicks D. The thermodynamics of DNA structural motifs. Annu Rev Biophys Biomol
%       Struct 2004;33:415-40. doi:10.1146/annurev.biophys.32.110601.141800
%   Rejali NA, Moric E, Wittwer CT. The effect of single mismatches on primer extension. Clin Chem
%       2018;64(5):801-809. doi:10.1373/clinchem.2017.282285
%   Aksel T, Navarro EJ, Fong N, Douglas SM. Design principles for accurate folding of DNA origami.
%       Proc Natl Acad Sci USA 2024;121(48):e2406769121. doi:10.1073/pnas.2406769121
%   Dunn KE, Dannenberg F, Ouldridge TE, Kwiatkowska M, Turberfield AJ, Bath J. Guiding the folding
%       pathway of DNA origami. Nature 2015;525:82-86. doi:10.1038/nature14860
%   Rivetti C, Walker C, Bustamante C. Polymer chain statistics and conformational analysis of DNA
%       molecules with bends or sections of different flexibility. J Mol Biol 1998;280:41-59.
%       doi:10.1006/jmbi.1998.1830
%   New England Biolabs. LongAmp Taq DNA Polymerase, M0323. Product page, accessed 12-11-2026:
%       https://www.neb.com/en-us/products/m0323-longamp-taq-dna-polymerase
% ____________________________________________________________________________________________________________

switch lower(char(mode))
    case 'scan'
        out = ScanMaxScaffold(varargin{:});
    case 'eval'
        out = EvalAmpliconSpan(varargin{:});
    otherwise
        error('ScaffoldCollapseCheck:BadMode', 'Unknown mode "%s". Use ''scan'' or ''eval''.', char(mode));
end
end

%% Scan Maximum Scaffold
function scanRes = ScanMaxScaffold(scaffSeq, W, minRun, endTol, monv, biv, dNTP, tann, useEps, opts)
scaffSeq = upper(char(scaffSeq));
Lmax = numel(scaffSeq);
W = min(W, floor((Lmax - 3)/2));
RT = 1.987E-3*(273.15 + tann);
N = NNParams(struct('Verbose', false));
if nargin < 10
    opts = [];
end
opts = LoopOptions(opts);

% The scan runs at the wider of the two flap tolerances. Each site records its own flap, so the Taq
% channel can drop the sites that only the exonuclease would extend.
scanTol = max(endTol, min(opts.ExoEndTol, W - minRun));

emptySite = struct('r0', 0, 'iTarget5', 0, 'iTarget3', 0, ...
    'iTail5', 0, 'iTail3', 0, 'loopLen', 0, ...
    'runLen', 0, 'dGstem', Inf, 'dGloop', Inf, ...
    'dGfold', Inf, 'eps', 1, 'dGeff', Inf, ...
    'epsPair', "", 'align', "", ...
    'dGfoldEng', Inf, 'tmFold', NaN, 'tmEng', NaN, ...
    'jOverC', NaN, 'tierRank', 0, ...
    'flap', 0, 'epsTaq', 1, 'dGbi', Inf, 'tmBi', NaN);

if (W < minRun)||(Lmax < (2*W + 3))
    scanRes = MakeScan(Lmax, W, RT, tann, endTol, useEps, opts, emptySite([]));
    return;
end

% 3'-terminal window, 5' to 3', where tailSeq(W) is the exposed 3'-OH base
tailSeq = scaffSeq(Lmax - W + 1:Lmax);
% Antiparallel frame: tailSeq(1:W) pairs against scaffSeq(r0:r0 + W - 1) read 3' to 5'
rcScaff = char(seqrcomplement(scaffSeq));

% Largest r0 that keeps the window inside the strand. A tail base k pairs only when it leaves a loop
% of at least 3 nt, which holds for k >= (4 + r0 + 2*W - Lmax)/2 and is enforced inside the loop below.
maxR0 = Lmax - W + 1;
if maxR0 < 1
    scanRes = MakeScan(Lmax, W, RT, tann, endTol, useEps, opts, emptySite([]));
    return;
end

% Physical stems keyed by anchor so shifted windows are not counted twice
byAnchor = containers.Map('KeyType', 'double', 'ValueType', 'any');

for r0 = 1:maxR0
    rcIdx = Lmax - (r0 + W - 1) + 1;
    basb = rcScaff(rcIdx:rcIdx + W - 1);
    isp = (tailSeq == basb);

    % Tail bases that would leave a loop under 3 nt cannot pair
    kMin = max(1, ceil((4 + r0 + 2*W - Lmax)/2));
    if kMin > 1
        isp(1:min(kMin - 1, W)) = false;
    end

    % The 3'-OH base must pair, or sit within scanTol of a pair. Tested first because it is cheap
    if ~any(isp(max(1, W - scanTol):W))
        continue;
    end

    % Skip windows with no run long enough to be a stem
    if MaxContigRun(isp) < minRun
        continue;
    end

    sub = BestSubRegister(tailSeq, basb, isp, minRun, monv, biv, dNTP, tann, false, minRun);
    if ~isfinite(sub.dGright)
        continue;
    end

    % The 3'-OH base must pair, or sit within scanTol of a pair
    iL = sub.iRight;
    jR = sub.jRight;
    if jR < (W - scanTol)
        continue;
    end
    flap = W - jR;

    % Stem bounds on scaffSeq: tailSeq(k) sits at Lmax - W + k and opposes r0 + W - k
    iTail5 = Lmax - W + iL;     % 5'-most paired base of the tail, closes the loop
    iTail3 = Lmax - W + jR;     % 3'-most paired base of the tail
    iTarget3 = r0 + W - iL;     % 3'-most paired base of the body arm, closes the loop
    iTarget5 = r0 + W - jR;     % 5'-most paired base of the body arm

    % Template must lie beyond the stem, and the loop must hold at least 3 nt
    loopLen = iTail5 - iTarget3 - 1;
    if (iTarget5 <= 1)||(loopLen < 3)
        continue;
    end

    % Fold geometry in the closing pair: 5' arm base on top, 3' arm base below
    stem5 = tailSeq(iL:jR);
    stem3 = char(seqcomplement(basb(iL:jR)));
    geom = struct('stem5', stem5, 'stem3', stem3, 'closeTop', stem3(1), 'closeBot', tailSeq(iL), ...
        'loopSeq', scaffSeq(iTarget3 + 1:iTail5 - 1), 'loopLen', loopLen, ...
        'monv', monv, 'biv', biv, 'dNTP', dNTP);

    % Primary fold energy (polymer loop beyond 30 nt) and its melting temperature
    [dGfold, dGstem, dGlp] = FoldEnergy(geom, tann, 'polymer', N, opts);
    if ~isfinite(dGfold)||(dGfold > 13.5)
        continue;
    end
    g10 = FoldEnergy(geom, tann + 10, 'polymer', N, opts);
    tmFold = FoldTm(dGfold, g10, tann, 10);

    % Engine loop form, logged beside the primary
    dGfoldEng = FoldEnergy(geom, tann, 'engine', N, opts);
    g10Eng = FoldEnergy(geom, tann + 10, 'engine', N, opts);
    tmEng = FoldTm(dGfoldEng, g10Eng, tann, 10);
    jOverC = exp(LnJEff(loopLen, opts))/opts.StrandConc;

    % Bimolecular state of the same stem on a second strand: initiation and terminal terms, no loop
    dGbi = sub.dGright;
    tmBi = NaN;
    if opts.Bimolecular
        gbi10 = dGMismatch(stem5, stem3, monv, biv, dNTP, tann + 10, false, true, true);
        tmBi = BimolecularTm(dGbi, gbi10, tann, tann + 10, opts.StrandConc/2, N.Rgas);
    end

    % Extension efficiency on the 3'-terminal three bases of the tail, kept for the Taq channel
    k0 = max(1, W - 2);
    pB = tailSeq(k0:W);
    tB = char(seqcomplement(basb(k0:W)));
    [epsTaq, einfo] = ExtensionEfficiency(pB, tB);
    epsVal = 1.0;
    pairTxt = "matched";
    if useEps
        epsVal = epsTaq;
        if einfo.GoverningPos > 0
            pairTxt = string(sprintf('%s@%d', einfo.GoverningPair, einfo.GoverningPos));
        else
            pairTxt = string(einfo.TerminalPair);
        end
    end

    dGeff = dGfold - RT*log(max(epsVal, realmin));
    tierRank = TierRank(epsVal, tmFold, tann, opts);

    bars = repmat(' ', 1, jR - iL + 1);
    bars(isp(iL:jR)) = '|';
    alnStr = string(sprintf('5''-%s-3'' / %s / 3''-%s-5'' (loop = %d nt)', ...
        stem5, bars, stem3, loopLen));

    site = struct('r0', iTarget5, 'iTarget5', iTarget5, 'iTarget3', iTarget3, ...
        'iTail5', iTail5, 'iTail3', iTail3, 'loopLen', loopLen, ...
        'runLen', sub.runRight, 'dGstem', dGstem, 'dGloop', dGlp, ...
        'dGfold', dGfold, 'eps', epsVal, 'dGeff', dGeff, ...
        'epsPair', pairTxt, 'align', alnStr, ...
        'dGfoldEng', dGfoldEng, 'tmFold', tmFold, 'tmEng', tmEng, ...
        'jOverC', jOverC, 'tierRank', tierRank, ...
        'flap', flap, 'epsTaq', epsTaq, 'dGbi', dGbi, 'tmBi', tmBi);

    % Keep the more stable of two windows that describe the same stem
    anchorKey = iTarget5*1E6 + iTail3;
    if ~isKey(byAnchor, anchorKey)
        byAnchor(anchorKey) = site;
    else
        prev = byAnchor(anchorKey);
        if dGeff < prev.dGeff
            byAnchor(anchorKey) = site;
        end
    end
end

vals = values(byAnchor);
if isempty(vals)
    sites = emptySite([]);
else
    sites = [vals{:}];
end
scanRes = MakeScan(Lmax, W, RT, tann, endTol, useEps, opts, sites);
end

%% Make Scan Result
function scanRes = MakeScan(Lmax, W, RT, tann, endTol, useEps, opts, sites)
scanRes = struct('Lmax', Lmax, 'W', W, 'RT', RT, 'tann', tann, 'endTol', endTol, ...
    'useEps', logical(useEps), 'opts', opts, 'sites', sites);
end

%% Evaluate Amplicon Span
function evalRes = EvalAmpliconSpan(scanRes, minR0, dGCollapseMax)
% Keeps the scanned folds whose body arm starts past the amplicon's 5' end by at least one base
evalRes = struct('dGfoldMin', Inf, 'dGeffMin', Inf, 'dGensEff', Inf, ...
    'pFold', 0, 'pExt', 0, 'loopLen', NaN, 'targetOffset', NaN, ...
    'eps', NaN, 'epsPair', "", 'align', "", ...
    'nSites', 0, 'pass', true, ...
    'dGfoldEng', Inf, 'tmFold', NaN, 'tmEng', NaN, 'jOverC', NaN, ...
    'pExtEng', 0, 'tier', "NONE", 'tierRank', 0, ...
    'pExtTot', 0, 'biShare', 0, 'pExtTaq', 0, 'pExtExo', 0, 'pExtBlend', NaN, ...
    'tmBi', NaN, 'tierAll', "NONE", 'tierAllRank', 0);

if isempty(scanRes.sites)
    return;
end

keep = [scanRes.sites.iTarget5] > minR0;
allSites = scanRes.sites(keep);
if isempty(allSites)
    return;
end

RT = scanRes.RT;
opts = scanRes.opts;
tierNames = {'NONE', 'TRANSIENT', 'FLAG', 'REJECT'};

% Headline set: the sites the declared enzyme channel can extend
isTaq = ([allSites.flap] <= scanRes.endTol);
if scanRes.useEps
    sites = allSites(isTaq);
else
    sites = allSites;
end

if ~isempty(sites)
    evalRes.nSites = numel(sites);
    [~, idxBest] = min([sites.dGeff]);
    win = sites(idxBest);

    boltz = sum(exp(-[sites.dGeff]/RT));
    dGens = -RT*log(max(boltz, realmin));

    % Mutually exclusive folds of one tail: Z = 1 + sum(K), pExt = sum(K*eps)/Z
    Kfold = exp(-[sites.dGfold]/RT);
    Z = 1 + sum(Kfold);

    evalRes.dGfoldMin = win.dGfold;
    evalRes.dGeffMin = win.dGeff;
    evalRes.dGensEff = dGens;
    evalRes.pFold = sum(Kfold)/Z;
    evalRes.pExt = sum(Kfold.*[sites.eps])/Z;

    % Log-only fields: engine loop form, persistence, and the strand tier
    KfoldEng = exp(-[sites.dGfoldEng]/RT);
    evalRes.pExtEng = sum(KfoldEng.*[sites.eps])/(1 + sum(KfoldEng));
    evalRes.dGfoldEng = win.dGfoldEng;
    evalRes.tmFold = win.tmFold;
    evalRes.tmEng = win.tmEng;
    evalRes.jOverC = win.jOverC;
    evalRes.tierRank = max([sites.tierRank]);
    evalRes.tier = string(tierNames{evalRes.tierRank + 1});
    evalRes.loopLen = win.loopLen;
    evalRes.targetOffset = win.iTarget5 - minR0 + 1;    % 1-based offset from the amplicon 5' end
    evalRes.eps = win.eps;
    evalRes.epsPair = win.epsPair;
    evalRes.align = win.align;
    evalRes.pass = (win.dGeff > dGCollapseMax);
end

% Channel totals over every physical state of the tail: foldback, bimolecular and unfolded
Ku = exp(-[allSites.dGfold]/RT);
if opts.Bimolecular
    Kb = opts.StrandConc*exp(-[allSites.dGbi]/RT);
else
    Kb = zeros(size(Ku));
end
Zall = 1 + sum(Ku) + sum(Kb);
epsT = [allSites.epsTaq].*isTaq;
evalRes.pExtTaq = sum((Ku + Kb).*epsT)/Zall;
evalRes.pExtExo = sum(Ku + Kb)/Zall;
if isfinite(opts.ExoBlend)
    evalRes.pExtBlend = (1 - opts.ExoBlend)*evalRes.pExtTaq + opts.ExoBlend*evalRes.pExtExo;
end
if scanRes.useEps
    epsDecl = epsT;
else
    epsDecl = ones(size(Ku));
end
totDecl = sum((Ku + Kb).*epsDecl);
evalRes.pExtTot = totDecl/Zall;
if totDecl > 0
    evalRes.biShare = sum(Kb.*epsDecl)/totDecl;
end
if opts.Bimolecular
    [~, ib] = max(Kb.*epsDecl);
    evalRes.tmBi = allSites(ib).tmBi;
end

% Worst tier over foldback and bimolecular states, in the Taq and the Exo channel
rk = 0;
for q = 1:numel(allSites)
    rk = max(rk, TierRank(epsT(q), allSites(q).tmFold, scanRes.tann, opts));
    rk = max(rk, TierRank(1, allSites(q).tmFold, scanRes.tann, opts));
    if opts.Bimolecular
        rk = max(rk, TierRank(epsT(q), allSites(q).tmBi, scanRes.tann, opts));
        rk = max(rk, TierRank(1, allSites(q).tmBi, scanRes.tann, opts));
    end
end
evalRes.tierAllRank = rk;
evalRes.tierAll = string(tierNames{rk + 1});
end

%% Max Contiguous Run
function rmax = MaxContigRun(isp)
rmax = 0; cur = 0;
for k = 1:numel(isp)
    if isp(k)
        cur = cur + 1;
        if cur > rmax, rmax = cur; end
    else
        cur = 0;
    end
end
end

%% Special Loop Bonus
% Sequence bonus of a tabulated triloop or tetraloop. Keyed by the closing pair's top base, the loop
% 5' to 3', then the closing pair's bottom base. Zero when the loop is not tabulated.
function g = SpecialLoop(closeTop, loopSeq, closeBot, tann, N)
g = 0;
if ~isfield(N, 'Loop')||~isfield(N.Loop, 'SpecialKey')||isempty(N.Loop.SpecialKey)
    return
end
key = upper([closeTop loopSeq closeBot]);
idx = find(strcmp(N.Loop.SpecialKey, key), 1);
if isempty(idx), return; end
g = N.Loop.SpecialdH(idx) - (tann + 273.15)*N.Loop.SpecialdS(idx)/1000;
end

%% Loop Options
% Fills any missing option field with its default
function opts = LoopOptions(opts)
def = struct('StrandConc', 1E-6, 'LBase', 0.6, 'LPers', 0.9, 'HpTmMargin', 5.0, 'FlagEps', 1E-3, ...
    'Bimolecular', true, 'ExoBlend', NaN, 'ExoEndTol', 15);
if (nargin < 1)||(isempty(opts))
    opts = def;
    return;
end
fn = fieldnames(def);
for k = 1:numel(fn)
    if (~isfield(opts, fn{k}))||(isempty(opts.(fn{k})))
        opts.(fn{k}) = def.(fn{k});
    end
end
end

%% Effective Concentration of a Loop
% Natural log of j(n) in mol/L for a worm-like chain of n nt
function lnj = LnJEff(n, opts)
Lc = n*opts.LBase;
d2 = 2*opts.LPers*Lc*(1 - (opts.LPers/Lc)*(1 - exp(-Lc/opts.LPers)));
jPerM3 = (3/(2*pi*d2*1E-18))^1.5;
lnj = log(jPerM3/6.02214076E23/1000);
end

%% Hairpin Loop Term
% Table value up to 30 nt. Beyond 30 nt the polymer model anchors at the 30 nt table value.
function g = HairpinLoop(n, tann, N, opts, model)
if (n > 30)&&(strcmp(model, 'polymer'))
    Tk = tann + 273.15;
    g = dGLoop('hairpin', 30, tann, N) + Tk*N.Rgas*(LnJEff(30, opts) - LnJEff(n, opts));
else
    g = dGLoop('hairpin', n, tann, N);
end
end

%% Fold Energy
% Stem plus hairpin loop at temperature T for the geometry in geom
function [dGfold, dGstem, dGlp] = FoldEnergy(geom, T, model, N, opts)
dGstem = dGMismatch(geom.stem5, geom.stem3, geom.monv, geom.biv, geom.dNTP, T, false, false, false);
dGlp = HairpinLoop(geom.loopLen, T, N, opts, model);
if geom.loopLen == 3
    dGlp = dGlp + SpecialLoop(geom.closeTop, geom.loopSeq, geom.closeBot, T, N);
    if any(geom.closeTop == 'AT')
        dGlp = dGlp + N.Loop.ATClosureG;
    end
else
    dGlp = dGlp + dGTerminalMismatch(geom.closeTop, geom.loopSeq(1), geom.closeBot, geom.loopSeq(end), T, N);
    if geom.loopLen == 4
        dGlp = dGlp + SpecialLoop(geom.closeTop, geom.loopSeq, geom.closeBot, T, N);
    end
end
dGfold = dGstem + dGlp;
end

%% Fold Melting Temperature
% Two-point estimate from dG at T and T + dT, degrees C. No concentration term for a unimolecular fold.
function tm = FoldTm(g1, g2, T, dT)
dS = (g1 - g2)/dT;
if abs(dS) <= 1E-7
    tm = -99;
else
    tm = (g1 + (T + 273.15)*dS)/dS - 273.15;
end
end

%% Tier Rank
% 3 reject, 2 flag, 1 transient, 0 none
function rk = TierRank(epsVal, tm, tann, opts)
persists = (tm >= tann - opts.HpTmMargin);
if (persists)&&(epsVal >= 0.999)
    rk = 3;
elseif (persists)&&(epsVal >= opts.FlagEps)
    rk = 2;
elseif (~persists)&&(epsVal >= opts.FlagEps)
    rk = 1;
else
    rk = 0;
end
end
