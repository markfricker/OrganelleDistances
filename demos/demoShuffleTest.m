%DEMOSHUFFLETEST  organelleShuffleTest on a synthetic cell: a random
% ER-like line network plus elongated "mitochondria", placed either
% hugging the ER (attraction) or at random (null), with the CDF /
% shuffled-envelope plot for each.

rootDir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(rootDir, 'src'), fullfile(rootDir, 'utils'));

p   = organelleDistanceParamsDefault();
cal = p.pixelSize;
sz  = [300 300];
rng(11);

% --- cell + nucleus -------------------------------------------------------
[X, Y]   = meshgrid(1:sz(2), 1:sz(1));
cellMask = ((X-150)/140).^2 + ((Y-150)/120).^2 <= 1;
nucleus  = (X-170).^2 + (Y-140).^2 <= 35^2;

% --- ER: random straight tubules clipped to the cytoplasm -----------------
er = false(sz);
for k = 1:25
    p0 = [randi(sz(2)) randi(sz(1))];
    th = rand * pi;
    t  = -80:0.5:80;
    xs = round(p0(1) + t*cos(th));  ys = round(p0(2) + t*sin(th));
    ok = xs >= 1 & xs <= sz(2) & ys >= 1 & ys <= sz(1);
    er(sub2ind(sz, ys(ok), xs(ok))) = true;
end
er = er & cellMask & ~nucleus;

% --- mitochondria: 3x9 px rods --------------------------------------------
rod = @(r, c, vertical) rodPix(r, c, vertical, sz);
nMito = 40;
cyto  = cellMask & ~nucleus;
dEr   = bwdist(er);

% attraction: rods whose centre lies 2-4 px from ER
cand = find(cyto & dEr >= 2 & dEr <= 4);
pixAttract = cell(nMito, 1);
occupied   = false(sz);
k = 0;
while k < nMito                     % rejection-place: inside cytoplasm, no overlaps
    [r, c] = ind2sub(sz, cand(randi(numel(cand))));
    idx = rod(r, c, rand > 0.5);
    if all(cyto(idx)) && ~any(occupied(idx))
        k = k + 1;
        pixAttract{k} = idx;
        occupied(idx) = true;
    end
end

% null: the same rods shuffled by the null model itself
pixNull = shuffleObjectsInWindow(pixAttract, cyto);

opts = struct('excludeMask', nucleus, 'rngSeed', 1, ...
    'restrictTargetToCell', p.shuffleRestrictTarget, 'dihedral', true);

figure('Name', 'organelleShuffleTest demo', 'Color', 'w');
cases = {pixAttract, 'ER-associated'; pixNull, 'Random (null)'};
for iCase = 1:2
    stats = table(cases{iCase,1}, ones(nMito,1), ...
        'VariableNames', {'organellePixelIdxList','cellID'});
    [curves, summary] = organelleShuffleTest({stats}, er, double(cellMask), ...
        1, cal, p.shuffleNSimulations, 'demo', opts);
    S = summary{1};  C = curves{1};

    lab = zeros(sz);
    for k = 1:nMito, lab(cases{iCase,1}{k}) = k; end
    subplot(2, 2, iCase);
    imshow(double(cat(3, er, lab > 0, bwperim(cellMask) | bwperim(nucleus))));
    title(cases{iCase,2});

    subplot(2, 2, 2 + iCase); hold on
    fill([C.r; flipud(C.r)], [C.cdfNullLo; flipud(C.cdfNullHi)], [0.6 0.9 0.6], ...
        'EdgeColor', 'none', 'DisplayName', 'shuffle 95% (pointwise)');
    plot(C.r, C.cdfNullMean, 'r-', 'LineWidth', 1.5, 'DisplayName', 'shuffle mean');
    plot(C.r, C.cdfObs, 'b-', 'LineWidth', 2, 'DisplayName', 'observed');
    xlabel('min organelle-to-ER distance (\mum)'); ylabel('cumulative fraction');
    title(sprintf('SDI = %.3f, p_{closer} = %.3f, p_{further} = %.3f', ...
        S.sdi, S.pCloser, S.pFurther));
    legend('Location', 'southeast'); box on
    fprintf('%-14s  SDI=%.3f  pCloser=%.3f  pFurther=%.3f  meanObs=%.3f  meanNull=%.3f um  fallback=%.3f\n', ...
        cases{iCase,2}, S.sdi, S.pCloser, S.pFurther, S.meanDistObs, S.meanDistNull, S.fallbackFraction);
end


function idx = rodPix(r, c, vertical, sz)
if vertical
    [R, C] = ndgrid(r + (-4:4), c + (-1:1));
else
    [R, C] = ndgrid(r + (-1:1), c + (-4:4));
end
ok  = R >= 1 & R <= sz(1) & C >= 1 & C <= sz(2);
idx = sub2ind(sz, R(ok), C(ok));
end
