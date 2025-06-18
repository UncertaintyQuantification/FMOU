%
% function [Lap, Lap_inv] = tri_laplacian(vertex, tri, allow)
% input: vertex    - n x 3, n: vertice number
%        tri       - m x 3, m: element number
%        allow       = 1, allow slip at the boundary, ~=1, does not allow
% output: Lap      - m x m
%         Lap_inv  - m x m
%
% Original version, zliu,  May 23, 2005 
% Minor cleanup,    zliu,  March 16, 2006 
%
function [Lap, Lap_inv] = tri_laplacian(vertex, tri, allow)
nvertex = size(vertex,1);
nele = size(tri,1);
fprintf('calculate laplacian matrix for %d fault patches ...\n', nele);

% calculate center coordinates of elements
a = vertex(tri(:,1:3),1); x = reshape(a, nele, 3);
a = vertex(tri(:,1:3),2); y = reshape(a, nele, 3);
a = vertex(tri(:,1:3),3); z = reshape(a, nele, 3);
cxyz = [sum(x,2)/3.0, sum(y,2)/3.0, sum(z,2)/3.0];
%
Lap = zeros(nele);
for i=1:nele
   [k,edge] = neighbors(i,tri); % neighbors and corresponding edge list of i, 1 x N vector
   
   % calculate distance between centers
   diff = repmat(cxyz(i,:),size(k,2),size(k,1)) - cxyz(k,:);
   h = sqrt(sum(diff.^2,2)); 
   
   % Laplacian term
   Li = sum(h);
   Lap(i,i) = -2.0/Li*sum(1.0./h);
   for j = 1: size(k,2)
      Lap(i,k(j)) = 2/Li/h(j);
   end
   
   % adjust for elements at the boundary, if requiring slip tailing to zero close to boundary
   % use image elements right outside boundary elements
   if allow ~= 1
      nb = size(k,1);
      if nb == 1 | nb == 2
         % find center coordinates at each edge, edge order [3, 1, 2]
         r0 = [x(i,:); y(i,:); z(i,:)];  % (vert1, vert2, vert3), along column (x;y;z)
         r1 = [r0(:,2) r0(:,3) r0(:,1)]; % (vert2, vert3, vert1)
         redg = (r0 + r1)/2.0;
         redg = redg'; % edge center in order 3, 1, 2 (across from vertex 3, 1, 2)
         redg = [redg(2,:);redg(3,:);redg(1,:)]; % into order 1,2,3
         diffe = repmat(cxyz(i,:),3, 1) - redg;
         dist = sqrt(sum(diffe.^2,2)); % same order as redg
         % 
         [c, ii] = setdiff([1 2 3], edge);
         % Have image element as "line element"
         h1(ii) = dist(c);
         h1(edge) = h;
         Li = sum(h1);
         Lap(i,i) = -2.0/Li*sum(1.0./h1);
         for j=1: size(edge,2)
            Lap(i,k(j)) = 2/Li/h(j); 
         end
      end
   end  
end
Lap = sparse(Lap);
Lap_inv = inv(Lap);
