function demoMitoErRays(erMask, mitoMask, opts)
% DEMOMITOERRAYS
% Demo for spline normals + ray intersection from mitochondria to ER.
%
% Inputs
%   erMask   : binary image of ER
%   mitoMask : binary image of mitochondria
%   opts     : optional struct with fields
%       .numSamples   (default: 400)   number of mito contour samples
%       .rayStride    (default: 5)     plot every Nth ray
%       .maxRange     (default: 500)   max ray length in pixels
%       .chunkSize    (default: 200)   intersection chunk size
%
% This demo:
%   1) extracts the largest contour from each mask
%   2) computes spline-smoothed mitochondrial normals
%   3) intersects each normal ray with the ER boundary
%   4) plots the masks, contours, rays, and intersection points

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'numSamples') || isempty(opts.numSamples)
        opts.numSamples = 400;
    end
    if ~isfield(opts, 'rayStride') || isempty(opts.rayStride)
        opts.rayStride = 5;
    end
    if ~isfield(opts, 'maxRange') || isempty(opts.maxRange)
        opts.maxRange = 500;
    end
    if ~isfield(opts, 'chunkSize') || isempty(opts.chunkSize)
        opts.chunkSize = 200;
    end

    erMask = logical(erMask);
    mitoMask = logical(mitoMask);

    % Extract the largest contour from each binary mask
    erContour = getLargestContour(erMask);
    mitoContour = getLargestContour(mitoMask);

    if isempty(erContour)
        error('Could not extract an ER contour.');
    end
    if isempty(mitoContour)
        error('Could not extract a mitochondria contour.');
    end

    % Compute ray intersections using spline-smoothed mito contour
    rayOpts = struct();
    rayOpts.numSamples = opts.numSamples;
    rayOpts.maxRange = opts.maxRange;
    rayOpts.chunkSize = opts.chunkSize;
    rayOpts.pad = 5;

    [distances, intersectionPoints, mitoSampled, normals] = ...
        normalsToNearestIntersectionSpline(mitoContour, erContour, rayOpts);

    % Plot
    figure('Color', 'w');
    ax = axes();
    hold(ax, 'on');

    % Background display: combined masks
    bg = zeros([size(erMask), 3]);
    bg(:,:,1) = 0.15 * double(mitoMask);   % red for mitochondria
    bg(:,:,2) = 0.15 * double(erMask);     % green for ER
    bg(:,:,3) = 0.15 * double(erMask);     % blue for ER

    imshow(bg, 'Parent', ax);
    axis(ax, 'image');
    set(ax, 'YDir', 'reverse');

    % Plot original contours
    plot(erContour(:,2), erContour(:,1), 'c-', 'LineWidth', 1.5);
    plot(mitoContour(:,2), mitoContour(:,1), 'y-', 'LineWidth', 1.5);

    % Plot spline-sampled mito contour
    plot(mitoSampled(:,2), mitoSampled(:,1), 'w.', 'MarkerSize', 4);

    % Plot normals/rays and intersection points
    rayIdx = 1:opts.rayStride:size(mitoSampled,1);
    for i = rayIdx
        if isnan(distances(i))
            continue;
        end

        p = mitoSampled(i, :);
        n = normals(i, :);
        q = intersectionPoints(i, :);

        % ray segment
        plot([p(2), q(2)], [p(1), q(1)], 'r-', 'LineWidth', 1);

        % start point and intersection point
        plot(p(2), p(1), 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 3);
        plot(q(2), q(1), 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 4);

        % optional tiny normal tick at start
        % plot([p(2), p(2)+8*n(2)], [p(1), p(1)+8*n(1)], 'r-');
    end

    title(sprintf('Mitochondria normals to ER intersections (%.0f rays)', numel(rayIdx)));
    legend({'ER contour','Mito contour','Mito samples','Ray','Mito point','Intersection'}, ...
           'TextColor', 'w', 'Location', 'bestoutside');

    hold(ax, 'off');

    % Optional: print a quick summary
    validDistances = distances(isfinite(distances));
    fprintf('Computed %d valid ray intersections.\n', numel(validDistances));
    if ~isempty(validDistances)
        fprintf('Mean distance: %.2f px\n', mean(validDistances));
        fprintf('Median distance: %.2f px\n', median(validDistances));
    end
end

function contour = getLargestContour(mask)
% GETLARGESTCONTOUR  Return the longest boundary from a binary mask
%
% Output is an Nx2 array [row col]

    B = bwboundaries(mask, 8, 'noholes');
    if isempty(B)
        contour = zeros(0,2);
        return;
    end

    lengths = cellfun(@(c) size(c,1), B);
    [~, idx] = max(lengths);
    contour = B{idx};
end