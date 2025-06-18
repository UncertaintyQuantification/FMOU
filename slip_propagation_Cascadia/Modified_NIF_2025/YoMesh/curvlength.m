function totlen = curvlength(pts)

% CURVLENGTH Calculate the length of a curve.
%   CURVLENGTH(PTS) calculates the length of a curve defined
%   by number of points PTS. PTS should either be a nx2 (2D)
%   or a nx3 (3D) matrix.
%
%   (Algorithm)
%   The distance between each neighbouring points are summed.
%
%   (Example)
%   t = linspace(0,2*pi,100);
%   x = cos(t); y = sin(t);
%   pts = [x;y]';
%   l = curvlength(pts)
%   figure; plot(pts(:,1),pts(:,2)); axis equal;
%
%   See also DISTANCE.
%

totlen = 0;
for k = 1:size(pts,1)-1
    l = distancept(pts(k,:),pts(k+1,:));
    totlen = totlen + l;
end
