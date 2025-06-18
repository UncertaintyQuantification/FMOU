function vec = getnormal(p1,p2,p3)

% getnormal Get normal vector of a plane that contains three points
%   in 3D space
%
%   vec = getnormal(PT1,PT2,PT3) calculates the unit normal
%   vector of a plane that passes three points PT1, PT2,
%   PT3. Points can be either a 1 x 3 or a 3 x 1 vector.
%
%   This function is essentially identical to planenormvec.
%   The only difference is that it outputs the unit vector
%   (norm(vec) = 1).
%
%   (Example)
%   p1 = [3,4,5];
%   p2 = [8,-4,0];
%   p3 = [0,0,1];
%   vec = getnormal(p1,p2,p3);
%   figure;
%   foo = [p1; p2; p3; p1];
%   plot3(foo(:,1),foo(:,2),foo(:,3));
%   hold on;
%   pm = mean(foo(1:3,:));
%   plot3(pm(1),pm(2),pm(3),'*r');
%   quiver3(pm(1),pm(2),pm(3),vec(1),vec(2),vec(3),3);
%
%   See also planenormvec.
%
%   16 Jul 2010, Yo Fukushima, DPRI, Kyoto Univ.
%

%% ChangeLogs
%  16 Jul 2010: first creation by modifying planenormvec.
%

% transpose if necessary
if size(p1,2) == 1
    p1 = p1';
end
if size(p2,2) == 1
    p2 = p2';
end
if size(p3,2) == 1
    p3 = p3';
end

p = cross(p1-p3,p2-p3);
if p(1)==0 & p(2)==0 & p(3)==0
    vec = [0 0 0];
    error([mfilename ': Three points are on a line.']);
    return
end

vec = p./norm(p);

return
