function contour = contourFromPerimeterIdx(perimIdx, imSize)
%CONTOURFROMPERIMETERIDX  Ordered perimeter linear indices -> [row col] contour.
%
%   contour = contourFromPerimeterIdx(perimIdx, imSize)
%
% INPUTS
%   perimIdx – Nx1 ordered linear pixel indices (e.g. one entry of a
%              PerimeterIdxList cell column, traced by bwtraceboundary).
%   imSize   – [nY nX] size of the frame the indices are linear into.
%
% OUTPUT
%   contour  – Nx2 [row col] ordered contour, NOT closed (caller/
%              downstream functions close it if needed).

[r, c] = ind2sub(imSize, perimIdx(:));
contour = [double(r), double(c)];
end
