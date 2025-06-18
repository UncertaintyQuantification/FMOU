function [nd,ind1,ind2] = intersectnd(nd1,nd2,tol)

% INTERSECTND Find nodes shared by two meshes.
%   [nd,ind1,ind2] = intersectnd(nd1,nd2) returns the nodes
%   common to both nd1 and nd2 (both nx2 or nx3 array).
%
%   ind1 and ind2 are the indice of the common nodes in
%   nd1 and nd2.
%
%   [c,ind1,ind2] = intersectnd(...,tol) uses tolerance
%   number to find the common nodes. Default is 1.
%
%   See also INTERSECT.
%
%   2 Jul 2005, Yo Fukushima

if ~exist('tol')
    tol = 1;
end

ind1 = []; ind2 = []; nd = [];
if size(nd1,2) == 3
    for k = 1:size(nd1,1)
        ind = find(nd1(k,1)-tol < nd2(:,1) & nd2(:,1) < nd1(k,1)+tol &...
            nd1(k,2)-tol < nd2(:,2) & nd2(:,2) < nd1(k,2)+tol &...
            nd1(k,3)-tol < nd2(:,3) & nd2(:,3) < nd1(k,3)+tol);
        if ~isempty(ind)
            ind1 = [ind1;k];
            ind2 = [ind2;ind];
            nd = [nd; nd2(ind,:)];
        end
    end
elseif size(nd1,2) == 2
    for k = 1:size(nd1,1)
        ind = find(nd1(k,1)-tol < nd2(:,1) & nd2(:,1) < nd1(k,1)+tol &...
            nd1(k,2)-tol < nd2(:,2) & nd2(:,2) < nd1(k,2)+tol);
        if ~isempty(ind)
            ind1 = [ind1;k];
            ind2 = [ind2;ind];
            nd = [nd; nd2(ind,:)];
        end
    end
end

