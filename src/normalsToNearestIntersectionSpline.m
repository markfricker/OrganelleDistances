function [distances, intersectionPoints, sourceSampled, normals, hitOwner] = normalsToNearestIntersectionSpline(sourceContour, targetContour, opts)
% NORMALSTONEARESTINTERSECTIONSPLINE
% Spline-smooth a source contour, compute outward normals, and find the
% nearest intersection of each normal ray with a target contour.
%
% The target contour may be a single closed curve (e.g. an ER skeleton
% polyline) or several disjoint closed curves concatenated with a
% [NaN NaN] separator row between each (e.g. the perimeters of every
% other object in a frame) -- NaN-bounded segments are automatically
% excluded from intersection candidates.
%
% Inputs
%   sourceContour : Mx2 [row col] ordered contour
%   targetContour : Kx2 [row col] ordered contour(s); NaN-separate
%                   multiple disjoint closed curves
%   opts          : struct with optional fields
%       .numSamples  (default: size(sourceContour,1))
%       .maxRange    (default: 500)
%       .chunkSize   (default: 200)
%       .pad         (default: 5)
%       .segOwner    ((size(targetContour,1)-1) x 1, optional) owner id
%                     per target segment (e.g. the object label that
%                     contributed each segment of a multi-curve target).
%                     When supplied, hitOwner reports the owner id of the
%                     nearest-hit segment instead of its raw index.
%
% Outputs
%   distances          : Nx1 distance along the normal ray
%   intersectionPoints : Nx2 [row col]
%   sourceSampled       : Nx2 spline-resampled source contour
%   normals            : Nx2 unit outward normals at sourceSampled
%   hitOwner           : Nx1 owner id (or raw segment index if
%                         opts.segOwner not supplied) of the nearest hit;
%                         NaN where no intersection was found.

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'numSamples') || isempty(opts.numSamples)
        opts.numSamples = size(sourceContour, 1);
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

    sourceContour = double(sourceContour);
    targetContour = double(targetContour);

    % Close the source contour if needed (target may be multi-curve /
    % NaN-separated and is used as-is).
    if any(sourceContour(1,:) ~= sourceContour(end,:))
        sourceContour = [sourceContour; sourceContour(1,:)];
    end

    % Spline-resample source contour and compute normals
    [sourceSampled, normals] = computeSplineNormals(sourceContour, opts.numSamples);

    % Flip normals outward using centroid test
    centroid = mean(sourceSampled, 1);
    flipMask = sum((sourceSampled - centroid) .* normals, 2) < 0;
    normals(flipMask, :) = -normals(flipMask, :);

    numPoints = size(sourceSampled, 1);

    % Target segments
    segStart = targetContour(1:end-1, :);
    segEnd   = targetContour(2:end, :);
    segVec   = segEnd - segStart;
    nSeg     = size(segStart, 1);

    haveOwner = isfield(opts, 'segOwner') && ~isempty(opts.segOwner);
    if haveOwner
        segOwner = opts.segOwner(:);
        if numel(segOwner) ~= nSeg
            error('normalsToNearestIntersectionSpline:segOwnerSize', ...
                'opts.segOwner must have one entry per target segment (size(targetContour,1)-1).');
        end
    else
        segOwner = (1:nSeg).';
    end

    segMinRow = min(segStart(:,1), segEnd(:,1));
    segMaxRow = max(segStart(:,1), segEnd(:,1));
    segMinCol = min(segStart(:,2), segEnd(:,2));
    segMaxCol = max(segStart(:,2), segEnd(:,2));

    distances = NaN(numPoints, 1);
    intersectionPoints = NaN(numPoints, 2);
    hitOwner = NaN(numPoints, 1);

    for startIdx = 1:opts.chunkSize:numPoints
        endIdx = min(startIdx + opts.chunkSize - 1, numPoints);
        idx = startIdx:endIdx;

        p = sourceSampled(idx, :);   % m x 2
        n = normals(idx, :);         % m x 2

        rayEnd = p + opts.maxRange * n;

        rayMinRow = min(p(:,1), rayEnd(:,1)) - opts.pad;
        rayMaxRow = max(p(:,1), rayEnd(:,1)) + opts.pad;
        rayMinCol = min(p(:,2), rayEnd(:,2)) - opts.pad;
        rayMaxCol = max(p(:,2), rayEnd(:,2)) + opts.pad;

        % Candidate target segments for this chunk
        candidateMask = ~(segMaxRow.' < rayMinRow | segMinRow.' > rayMaxRow | ...
                          segMaxCol.' < rayMinCol | segMinCol.' > rayMaxCol);

        segmentMask = any(candidateMask, 1);
        if ~any(segmentMask)
            continue;
        end

        s1 = segStart(segmentMask, :);   % k x 2
        v  = segVec(segmentMask, :);     % k x 2
        globalSegIdx = find(segmentMask);

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

        [tMin, hitCol] = min(t, [], 2);
        hitMask = isfinite(tMin);

        if any(hitMask)
            distances(idx(hitMask)) = tMin(hitMask);
            intersectionPoints(idx(hitMask), :) = p(hitMask, :) + tMin(hitMask) .* n(hitMask, :);
            hitSegGlobal = globalSegIdx(hitCol(hitMask));
            hitOwner(idx(hitMask)) = segOwner(hitSegGlobal);
        end
    end
end