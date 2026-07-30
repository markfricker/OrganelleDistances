function [sampledContour, normals] = computeSplineNormals(contour, numSamples, smoothSpanPx)
% COMPUTESPLINENORMALS
% Spline-interpolate a closed contour, resample it at numSamples output
% points, and return a unit tangent-normal pair at each.
%
% contour        : Nx2 [row col], closed (first point = last point)
% numSamples     : number of output points (ray/plot density).
% smoothSpanPx   : (optional) arc-length window, in pixels, used to
%                  estimate each tangent via a centered finite difference.
%                  Decoupled from numSamples: the spline itself always
%                  interpolates exactly through every input contour point
%                  (no curve-smoothing happens there), so at full output
%                  density a naive gradient() between adjacent OUTPUT
%                  points would just be differencing immediately-adjacent
%                  raw pixels -- exactly the per-pixel staircase noise a
%                  "smooth" control is supposed to avoid. This window
%                  fixes the tangent-estimation baseline at a constant
%                  pixel distance regardless of how densely numSamples
%                  samples the same curve.
%                  Default (omitted/[]): totalArcLength/numSamples, i.e.
%                  the legacy behaviour where smoothing tracked output
%                  spacing -- kept only for callers that don't know about
%                  this parameter.
%
% sampledContour : numSamples x 2 [row col]
% normals        : numSamples x 2 unit normals [row col]

    if nargin < 3
        smoothSpanPx = [];
    end

    contour = double(contour);

    % Catch a non-finite contour here, at the point of entry, with a
    % diagnostic naming exactly which point(s) are bad -- griddedInterpolant
    % (called from interp1 below) only reports "sample points must be
    % finite" with no indication of which value or where it came from.
    if any(~isfinite(contour(:)))
        badRows = find(any(~isfinite(contour), 2));
        error('computeSplineNormals:nonFiniteContour', ...
            ['contour has %d non-finite [row col] point(s), first at row %d: [%s]. ' ...
            'This means the caller passed in a contour containing NaN/Inf -- ' ...
            'if built via contourFromPerimeterIdx, check the PerimeterIdxList ' ...
            'entry and imSize used to decode it (a non-finite or out-of-range ' ...
            'linear index decodes to NaN via ind2sub without erroring).'], ...
            numel(badRows), badRows(1), mat2str(contour(badRows(1),:)));
    end

    if any(contour(1,:) ~= contour(end,:))
        contour = [contour; contour(1,:)];
    end

    x = contour(1:end-1, 2);   % col
    y = contour(1:end-1, 1);   % row

    % Arc-length parameter, in pixels, before normalising to [0,1]
    ds = sqrt(diff(x).^2 + diff(y).^2);
    t = [0; cumsum(ds)];
    if t(end) == 0
        error('Contour has zero length.');
    end
    totalLenPx = t(end);

    if isempty(smoothSpanPx)
        smoothSpanPx = totalLenPx / numSamples;   % legacy-equivalent default
    end
    halfWindowNorm = (smoothSpanPx / 2) / totalLenPx;

    t = t / totalLenPx;

    % Closed, uniformly spaced output sample locations
    tq = linspace(0, 1, numSamples + 1).';
    tq(end) = [];

    % Cyclic padding to reduce endpoint artifacts with interp1(...,'spline').
    % Must cover at least the smoothing half-window (plus a small margin),
    % not just a fixed handful of original segments -- otherwise a wide
    % smoothing window queries interp1 outside the padded domain.
    avgSegNorm = 1 / (numel(t) - 1);
    padNeeded  = ceil((halfWindowNorm + 2 * avgSegNorm) / avgSegNorm);
    pad        = min(max(padNeeded, 3), numel(t) - 1);

    tPad = [t(end-pad:end-1) - 1; t; t(2:pad+1) + 1];
    xPad = [x(end-pad:end-1); x; x(2:pad+1)];
    yPad = [y(end-pad:end-1); y; y(2:pad+1)];

    % Defense in depth -- the contour itself was already checked above, but
    % this catches anything introduced by the arc-length computation itself
    % (e.g. a duplicate-point run producing a degenerate ds/t) before it
    % reaches griddedInterpolant's opaque error.
    if any(~isfinite(tPad)) || any(~isfinite(xPad)) || any(~isfinite(yPad))
        error('computeSplineNormals:nonFiniteInterpInput', ...
            ['Non-finite value(s) in the padded spline-interpolation input. ' ...
            'numel(t)=%d, t range=[%s, %s], any(~isfinite(t))=%d, ' ...
            'any(~isfinite(x))=%d, any(~isfinite(y))=%d.'], ...
            numel(t), num2str(min(t)), num2str(max(t)), ...
            any(~isfinite(t)), any(~isfinite(x)), any(~isfinite(y)));
    end

    xs = interp1(tPad, xPad, tq, 'spline');
    ys = interp1(tPad, yPad, tq, 'spline');

    % Tangent via a centered finite difference over the fixed smoothing
    % window, independent of how densely tq itself is sampled.
    tqPlus  = tq + halfWindowNorm;
    tqMinus = tq - halfWindowNorm;
    xPlus  = interp1(tPad, xPad, tqPlus,  'spline');
    xMinus = interp1(tPad, xPad, tqMinus, 'spline');
    yPlus  = interp1(tPad, yPad, tqPlus,  'spline');
    yMinus = interp1(tPad, yPad, tqMinus, 'spline');

    dx = xPlus - xMinus;
    dy = yPlus - yMinus;

    tangent = [dy, dx];  % [row col] order
    tangent = tangent ./ (sqrt(sum(tangent.^2, 2)) + eps);

    % Rotate tangent by +90 degrees to get normals
    normals = [-tangent(:,2), tangent(:,1)];
    normals = normals ./ (sqrt(sum(normals.^2, 2)) + eps);

    sampledContour = [ys, xs];
end
