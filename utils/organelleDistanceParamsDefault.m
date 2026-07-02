function p = organelleDistanceParamsDefault()
%ORGANELLEDISTANCEPARAMSDEFAULT  Factory default parameters for the
% spline ray-cast organelle distance library (organelleD2ErCompute,
% organelleD2OrganelleCompute).
%
%   p = organelleDistanceParamsDefault()
%
%   Ray-cast core (normalsToNearestIntersectionSpline / computeSplineNormals)
%   ----------------------------------------------------------------------
%   span         Perimeter sample-density control (px). Higher span ->
%                fewer, smoother spline samples around each object's
%                perimeter (numSamples = max(8, round(nPerim/span))).
%                span<=1 keeps one sample per raw traced perimeter pixel.
%   chunkSize    Rays processed per vectorised batch (perf tuning only;
%                does not affect results).
%   pad          px padding added to each ray's coarse bounding-box
%                candidate-segment filter.
%
%   organelleD2ErCompute (mito -> ER)
%   ----------------------------------------------------------------------
%   d2ErRadius     Outward ray search length toward the ER (px).
%
%   organelleD2OrganelleCompute (mito -> mito)
%   ----------------------------------------------------------------------
%   nnRadius        Outward ray search length toward other objects (px).
%                   Also used, with a margin, as a coarse centroid-distance
%                   prefilter when building each object's target polyline.
%   neighborRadius  µm; radius used for the neighbour-count / percent-
%                   perimeter-near clustering readout.
%
%   Acquisition
%   ----------------------------------------------------------------------
%   pixelSize    µm per pixel (calibration).

p.span        = 3;      % px
p.chunkSize   = 200;
p.pad         = 5;      % px

p.d2ErRadius  = 50;     % px search length toward the ER

p.nnRadius        = 50; % px search length toward other objects
p.neighborRadius  = 1;  % um; "crowded" cutoff for neighbour-count/pct-near

p.pixelSize   = 0.09;   % um/px
end
