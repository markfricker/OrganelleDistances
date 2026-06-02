function [sampledContour, normals] = computeSplineNormals(contour, numSamples)
% COMPUTESPLINENORMALS
% Resample a closed contour using a spline-smoothed periodic approximation,
% then return unit tangent-normal pairs.
%
% contour: Nx2 [row col], closed (first point = last point)
% sampledContour: numSamples x 2 [row col]
% normals: numSamples x 2 unit normals [row col]

    contour = double(contour);

    if any(contour(1,:) ~= contour(end,:))
        contour = [contour; contour(1,:)];
    end

    x = contour(1:end-1, 2);   % col
    y = contour(1:end-1, 1);   % row

    % Arc-length parameter
    ds = sqrt(diff(x).^2 + diff(y).^2);
    t = [0; cumsum(ds)];
    if t(end) == 0
        error('Contour has zero length.');
    end
    t = t / t(end);

    % Closed, uniformly spaced sample locations
    tq = linspace(0, 1, numSamples + 1).';
    tq(end) = [];

    % Cyclic padding to reduce endpoint artifacts with interp1(...,'spline')
    pad = min(3, numel(t) - 1);

    tPad = [t(end-pad:end-1) - 1; t; t(2:pad+1) + 1];
    xPad = [x(end-pad:end-1); x; x(2:pad+1)];
    yPad = [y(end-pad:end-1); y; y(2:pad+1)];

    xs = interp1(tPad, xPad, tq, 'spline');
    ys = interp1(tPad, yPad, tq, 'spline');

    % Tangent from spline-smoothed contour
    dx = gradient(xs);
    dy = gradient(ys);

    tangent = [dy, dx];  % [row col] order
    tangent = tangent ./ (sqrt(sum(tangent.^2, 2)) + eps);

    % Rotate tangent by +90 degrees to get normals
    normals = [-tangent(:,2), tangent(:,1)];
    normals = normals ./ (sqrt(sum(normals.^2, 2)) + eps);

    sampledContour = [ys, xs];
end