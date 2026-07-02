function demoMitoMitoRays(labelIm, opts)
% DEMOMITOMITORAYS
% Demo for spline normals + ray intersection between mitochondria (i.e.
% inter-organelle nearest-neighbour rays), the mito-mito analogue of
% demoMitoErRays.
%
% Inputs
%   labelIm  : label matrix, one integer per mitochondrion (0 = background)
%   opts     : optional struct with fields
%       .numSamples   (default: 200)    number of contour samples per object
%       .rayStride    (default: 3)      plot every Nth ray
%       .maxRange     (default: 100)    max ray length in pixels
%       .chunkSize    (default: 200)    intersection chunk size
%
% This demo:
%   1) extracts every object's perimeter
%   2) for every object, builds the "all other objects" target polyline
%      (otherObjectsPolyline) and casts spline-smoothed outward normal rays
%      against it (normalsToNearestIntersectionSpline)
%   3) plots the label image, contours, nearest-neighbour rays coloured by
%      distance, and highlights each object's single nearest-neighbour ray

    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'numSamples') || isempty(opts.numSamples)
        opts.numSamples = 200;
    end
    if ~isfield(opts, 'rayStride') || isempty(opts.rayStride)
        opts.rayStride = 3;
    end
    if ~isfield(opts, 'maxRange') || isempty(opts.maxRange)
        opts.maxRange = 100;
    end
    if ~isfield(opts, 'chunkSize') || isempty(opts.chunkSize)
        opts.chunkSize = 200;
    end

    imSize = size(labelIm);
    nO = max(labelIm(:));
    if nO < 2
        error('demoMitoMitoRays:tooFewObjects', 'Need at least 2 labelled objects.');
    end

    perimList = cell(nO, 1);
    centroids = zeros(nO, 2);
    for k = 1:nO
        mask = labelIm == k;
        if ~any(mask(:))
            continue
        end
        [rows, cols] = find(mask);
        centroids(k, :) = [mean(cols), mean(rows)];  % [x y]
        [r0, c0] = find(mask, 1);
        B = bwtraceboundary(mask, [r0 c0], 'N', 8);
        perimList{k} = sub2ind(imSize, B(:,1), B(:,2));
    end

    figure('Color', 'w');
    ax = axes();
    hold(ax, 'on');
    imshow(label2rgb(labelIm, 'lines', 'k', 'shuffle'), 'Parent', ax);
    axis(ax, 'image');
    set(ax, 'YDir', 'reverse');

    nnDist = nan(nO, 1);
    nnLine = cell(nO, 1);

    for iO = 1:nO
        if isempty(perimList{iO})
            continue
        end
        [targetPoly, segOwner] = otherObjectsPolyline(perimList, imSize, iO, centroids, opts.maxRange * 3);
        if isempty(targetPoly)
            continue
        end
        srcContour = contourFromPerimeterIdx(perimList{iO}, imSize);

        rayOpts = struct('numSamples', opts.numSamples, 'maxRange', opts.maxRange, ...
            'chunkSize', opts.chunkSize, 'pad', 5, 'segOwner', segOwner);
        [dist, hitPts, sampled] = normalsToNearestIntersectionSpline(srcContour, targetPoly, rayOpts);

        rayIdx = 1:opts.rayStride:size(sampled, 1);
        for i = rayIdx
            if isnan(dist(i))
                continue
            end
            p = sampled(i, :);
            q = hitPts(i, :);
            plot(ax, [p(2), q(2)], [p(1), q(1)], 'c-', 'LineWidth', 0.5);
        end

        [minD, iMin] = min(dist);
        if isfinite(minD)
            nnDist(iO) = minD;
            nnLine{iO} = [sampled(iMin, :); hitPts(iMin, :)];
        end
    end

    for iO = 1:nO
        if isempty(nnLine{iO})
            continue
        end
        L = nnLine{iO};
        plot(ax, [L(1,2) L(2,2)], [L(1,1) L(2,1)], 'r-', 'LineWidth', 2);
        plot(ax, L(1,2), L(1,1), 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
    end

    title(sprintf('Mito-mito nearest-neighbour rays (%d objects, red = NN)', nO));
    hold(ax, 'off');

    fprintf('Nearest-neighbour distances (px):\n');
    for iO = 1:nO
        fprintf('  object %d: %.2f\n', iO, nnDist(iO));
    end
end
