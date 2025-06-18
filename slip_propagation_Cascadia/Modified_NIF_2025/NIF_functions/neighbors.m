%
% Find adjacent element for given element in the triangular mesh list
%
% function [k,edge] = neighbors(i,tri)
% input:  i -- center element
%       tri -- triangular mesh connectivity list, N x 3
% output: k -- adjacent element number, 1 x m, m <= 3, m - total # of neighbors
%      edge -- edge list that has neighbor in the order [1 2 3]
%
% written by Z. Liu                  May 18, 2005
%
function [k,edge] = neighbors(i,tri)
if i > size(tri,1)
   error([' Element ', num2str(i), ' exceeds matrix dimensions']);
elseif size(tri,2) < 3
   error([' Connectivity list has to have at least 3 vertices!']);
elseif size(tri,1)== 3 & size(tri,2) > 3
   error([' vertex list should be n x 3 matrice ! ']);
elseif size(tri,2) > 3
   error([' Current version only allow triangular mesh !']);
end
[v1,v2,v3]= deal(tri(i,1),tri(i,2),tri(i,3));

% find element adjacent to side 3, which is across vertice 3
m1 = find((tri(:,1)== v1 & tri(:,2)== v2) | (tri(:,2) == v1 & tri(:,3) == v2) | ...
   (tri(:,3)== v1 & tri(:,1) == v2));
m2 = find((tri(:,2)== v1 & tri(:,1)== v2) | (tri(:,3) == v1 & tri(:,2) == v2) | ...
   (tri(:,1)== v1 & tri(:,3) == v2));
n3 = union(m1,m2);
n3 = setdiff(n3,i);
ed3 = 3*ones(1,length(n3)); % row vector
ed3 = unique(ed3);

% find element adjacemnt to side 2, across vertice 2
m1 = find((tri(:,1)== v3 & tri(:,2)== v1) | (tri(:,2) == v3 & tri(:,3) == v1) | ...
   (tri(:,3)== v3 & tri(:,1) == v1));
m2 = find((tri(:,2)== v3 & tri(:,1)== v1) | (tri(:,3) == v3 & tri(:,2) == v1) | ...
   (tri(:,1)== v3 & tri(:,3) == v1));
n2 = union(m1,m2);
n2 = setdiff(n2,i);
ed2 = 2*ones(1,length(n2));
ed2 = unique(ed2);

% find element adjacemnt to side 1, across vertice 1
m1 = find((tri(:,1)== v2 & tri(:,2)== v3) | (tri(:,2) == v2 & tri(:,3) == v3) | ...
   (tri(:,3)== v2 & tri(:,1) == v3));
m2 = find((tri(:,2)== v2 & tri(:,1)== v3) | (tri(:,3) == v2 & tri(:,2) == v3) | ...
   (tri(:,1)== v2 & tri(:,3) == v3));
n1= union(m1,m2);
n1 = setdiff(n1,i);
ed1 = 1*ones(1,length(n1));;
ed1 = unique(ed1);

% k, ed - row vector, in the edge order [1, 2, 3] across vertices [1 2 3]
k = [n1 n2 n3];
edge = [ed1 ed2 ed3];
