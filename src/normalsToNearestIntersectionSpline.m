function [distances, intersectionPoints, mitoSampled, normals] = normalsToNearestIntersectionSpline(mitoContour, erContour, opts)
% NORMALSTONEARESTINTERSECTIONSPLINE
% Spline-smooth mitochondrial contour, compute outward normals, and find
% the nearest intersection of each normal ray with the ER boundary.
%
% Inputs
%   mitoContour : Mx2 [row col] ordered contour
%   erContour   : Kx2 [row col] ordered closed contour
%   opts        : struct with optional fields
%       .numSamples  (default: size(mitoContour,1))
%       .maxRange    (default: 500)
%       .chunkSize   (default: 200)
%       .pad         (default: 5)
%
% Outputs
%   distances         : Nx1 distance along the normal ray
%   intersectionPoints : Nx2 [row col]
%   mitoSampled       : Nx2 spline-resampled mito contour
%   normals           : Nx2 unit outward normals at mitoSampled

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'numSamples') || isempty(opts.numSamples)
        opts.numSamples = size(mitoContour, 1);
    end
    if ~isfield(opts, 'maxRange') || isempty(opts.maxRange)
        opts.maxRange = 500;
    end
    if ~isfield(opts, 'chunkSize') || isempty(opts.chunkSize)
        opts.chunkSize = 200;
    end
    if ~isfield(opts, 'pad') || isempty(opts.pad)
        opts.pad = 5;
    end

    mitoContour = double(mitoContour);
    erContour = double(erContour);

    % Close contours if needed
    if any(mitoContour(1,:) ~= mitoContour(end,:))
        mitoContour = [mitoContour; mitoContour(1,:)];
    end
    if any(erContour(1,:) ~= erContour(end,:))
        erContour = [erContour; erContour(1,:)];
    end

    % Spline-resample mito contour and compute normals
    [mitoSampled, normals] = computeSplineNormals(mitoContour, opts.numSamples);

    % Flip normals outward using centroid test
    centroid = mean(mitoSampled, 1);
    flipMask = sum((mitoSampled - centroid) .* normals, 2) < 0;
    normals(flipMask, :) = -normals(flipMask, :);

    numPoints = size(mitoSampled, 1);

    % ER segments
    segStart = erContour(1:end-1, :);
    segEnd   = erContour(2:end, :);
    segVec   = segEnd - segStart;

    segMinRow = min(segStart(:,1), segEnd(:,1));
    segMaxRow = max(segStart(:,1), segEnd(:,1));
    segMinCol = min(segStart(:,2), segEnd(:,2));
    segMaxCol = max(segStart(:,2), segEnd(:,2));

    distances = NaN(numPoints, 1);
    intersectionPoints = NaN(numPoints, 2);

    for startIdx = 1:opts.chunkSize:numPoints
        endIdx = min(startIdx + opts.chunkSize - 1, numPoints);
        idx = startIdx:endIdx;

        p = mitoSampled(idx, :);   % m x 2
        n = normals(idx, :);       % m x 2
        m = size(p, 1);

        rayEnd = p + opts.maxRange * n;

        rayMinRow = min(p(:,1), rayEnd(:,1)) - opts.pad;
        rayMaxRow = max(p(:,1), rayEnd(:,1)) + opts.pad;
        rayMinCol = min(p(:,2), rayEnd(:,2)) - opts.pad;
        rayMaxCol = max(p(:,2), rayEnd(:,2)) + opts.pad;

        % Candidate ER segments for this chunk
        candidateMask = ~(segMaxRow.' < rayMinRow | segMinRow.' > rayMaxRow | ...
                          segMaxCol.' < rayMinCol | segMinCol.' > rayMaxCol);

        segmentMask = any(candidateMask, 1);
        if ~any(segmentMask)
            continue;
        end

        s1 = segStart(segmentMask, :);   % k x 2
        v  = segVec(segmentMask, :);     % k x 2
        k  = size(s1, 1);

        % Pairwise ray/segment intersections for this chunk:
        % p + t*n = s1 + u*v
        % denom = cross(n, v)
        denom = n(:,1) * v(:,2).' - n(:,2) * v(:,1).';   % m x k

        rhsRow = s1(:,1).' - p(:,1);  % m x k
        rhsCol = s1(:,2).' - p(:,2);  % m x k

        t = (rhsRow .* v(:,2).' - rhsCol .* v(:,1).') ./ denom;
        u = (rhsRow .* n(:,2)   - rhsCol .* n(:,1))   ./ denom;

        valid = candidateMask(:, segmentMask) & abs(denom) > 1e-12 & ...
                t > 1e-6 & u >= -1e-9 & u <= 1 + 1e-9;

        t(~valid) = inf;

        [tMin, ~] = min(t, [], 2);
        hitMask = isfinite(tMin);

        if any(hitMask)
            distances(idx(hitMask)) = tMin(hitMask);
            intersectionPoints(idx(hitMask), :) = p(hitMask, :) + tMin(hitMask) .* n(hitMask, :);
        end
    end
end