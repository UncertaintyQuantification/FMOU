function l = distancept(x,y)

% DISTANCEPT Calculate the distance.
%   DISTANCEPT(X,Y) calculates the distance between two
%   points X and Y.
%
%   When X is a 1 x 2 (2D) or a 1 x 3 (3D) vector and
%   Y is a n x 2 matrix (for 2D) or n x 3 matrix (for 3D),
%   where n is the number of points, it returns the
%   distances between X and all the points in Y.
%
%   When X is a matrix of the same size as Y, then it
%   returns distances between each corresponding points
%   (i.e., each row of the matrix) in X and Y.
%
%   (Example)
%   x = [5 5 5];
%   foo = linspace(0,10,10)';
%   y = [foo,foo,foo];
%   l = distancept(x,y);
%   figure;plot(foo,l,'.');
%
%   See also TOTALDISTANCE.
%
%   13 Jul 2005, Yo Fukushima

%% calculate distance %%
if size(x,2) == 2
    if size(x,1) == 1
        l = sqrt((x(1)-y(:,1)).^2+(x(2)-y(:,2)).^2);
    else
        l = sqrt((x(:,1)-y(:,1)).^2+(x(:,2)-y(:,2)).^2);
    end
elseif size(x,2) == 3
    if size(x,1) == 1
        l = sqrt((x(1)-y(:,1)).^2+(x(2)-y(:,2)).^2+(x(3)-y(:,3)).^2);
    else
        l = sqrt((x(:,1)-y(:,1)).^2+(x(:,2)-y(:,2)).^2+(x(:,3)-y(:,3)).^2);
    end
else
    error('Number of dimensions should be 2 or 3.');
end
