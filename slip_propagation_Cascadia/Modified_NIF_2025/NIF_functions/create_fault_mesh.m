function create_fault_mesh(contourfile, intv, minLat, maxLat, minLong, maxLong, minDepth, maxDepth, savefile)
%
% This script creates a mesh of a specified fault, from a contour
% file containing (lat,lon,depth).
% 
% inputs: 
%
% contourfile: file containing lat, long, depth entries defining a fault.  See subduction_contours.txt for an example.
%
% intv: aproximate dimension of meshed triangles (meshing interval), in km
%
% minLat, maxLat, minLong, maxLong, minDepth, maxDepth: 
% These variables define a volume of interest for creating the mesh, and ignore contour points outside this volume
% minDepth should be the shallowest point, in km, and should be a negative value
% maxDepth should be the deepest point, in km, and should be a negative value
%
% savefile: string name for a .mat file including the resulting mesh variables
%     mesh variables in output file:
%     el: 3 by n array defining each triangle using indices that refer to nodes
%     nd: nodes of the triangular mesh, in local coordinates
%     nd_ll: nodes of the triangular mesh, in lat-long coordinates
%     origin: used to define the local coordinates
%
% originally created by Yo Fukushima, 9 Jul 2010 - please acknowledge
% citation for first usage of this code: 
% Fukushima, Y., V. Cayol, and P. Durand (2005), Finding realistic dike models from interferometric synthetic aperture radar data: 
% The February 2000 eruption at Piton de la Fournaise, J. Geophys. Res., 110, B03206, doi:10.1029/2004JB003268.

%% ChangeLogs
% 8 Jul 2010 (YF): initial creation
% 9 Jul 2010 (YF): work around with griddata problem of not interpolating
%                  on edges
% 28 Jul 2010 (YF) changed to eliminate ends of mesh for inversion
% 20 Jun 2015 (NB) changed to be a function with inputs, cleaned up

%% load contour data
contr = load(contourfile);
contr = contr(~isnan(contr(:,1)),:); % remove nan
%correct for positive depth
%contr(:, 3) = -contr(:, 3);

%rotate for testing :)
%contr(:,  1)=contr(:,  1)+contr(:,  2)/3-16;

%% eliminate points outside area of interest
Ind = find(contr(:,1) > minLat & contr(:,1) < maxLat);
contr = contr(Ind,:);
Ind = find(contr(:,2) > minLong &  contr(:,2) < maxLong);
contr = contr(Ind,:);
Ind = find(contr(:,3) <minDepth  &  contr(:,3) > maxDepth);
contr = contr(Ind,:);


%% convert to local km
% PS change - Use local origin fom grab_data...
origin = [mean(contr(:,1)); mean(contr(:,2))];

contrkm = llh2local([contr(:,2)'; contr(:,1)'],[origin(2) ;origin(1)])';
contrkm = [contrkm, contr(:,3)];

figure;scatter(contrkm(:,1),contrkm(:,2),10,contrkm(:,3),'filled'); colorbar; 
xlabel('Easting'); ylabel('Northing'); daspect([1,1,1]);


%% identify the points of the shallowest and deepest contours
% nodet: the shallowest points
% nodeb: the deepest points
ind = (contrkm(:,3) == contrkm(1,3));% hardwired since first element is shallowest
nodet0 = contrkm(ind,:);
ind = (contrkm(:,3) == contrkm(end,3));
nodeb0 = contrkm(ind,:);
% here I remove edge points to make sure that griddata works properly
%nodet0(end,:) = [];
%nodet0(1,:) = [];
%nodeb0(end,:) = [];
%nodet0(1,:) = [];


%% set the intervals of shallowest and deepest points as specified 
 N = round(curvlength(nodet0)./intv);
 nodet = curvspace(nodet0,N);
 N = round(curvlength(nodeb0)./intv);
 nodeb = curvspace(nodeb0,N);


%% mesh generation
%IMPORTANT: meshfrac2 assumes a linear progression of depth from the
%shallow limit to the deep limit.  This is BAD.
[nd0,el] = meshfrac2(nodet,nodeb,intv);


% nd0  are coordinates of nodes on the mesh
% el   are the elements, the values specify the x,y,z indices of the nodes

%% shift the height of the nodes based on the contour info
% the mesh previously made is a straight-line interpolation from the top
% nodes to bottom nodes
[xi,yi,zi] = griddata(contrkm(:,1),contrkm(:,2),contrkm(:,3),...
    nd0(:,1),nd0(:,2),'linear');
ind = find(nd0(:,3) ~= nd0(1,3) & nd0(:,3) ~= nd0(end,3));
nd = nd0;
nd(ind,:) = [xi(ind),yi(ind),zi(ind)];

figure;
trisurf(el,nd(:,1),nd(:,2),nd(:,3)); colorbar; view(2); daspect([1,1,1]);
xlabel('Easting'); ylabel('Northing'); daspect([1,1,1]);


%% optionally... convert nodes to lon/lat
%nd_ll = local2llh([nd0(:,1)'; nd0(:,2)'],origin)';
% PS fix to change in origin
nd_ll = local2llh([nd(:,1)'; nd(:,2)'],[origin(2) ;origin(1)])';
nd_ll = [nd_ll, nd(:,3)];

figure;
trisurf(el,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3)); colorbar; view(2); 
xlabel('Longitude'); ylabel('Latitude');
%since I shifted to using origin from grabdata
%daspect([1/cosd(abs(origin(2))),1,110]);
daspect([1/cosd(abs(origin(1))),1,110]);

%% save
save(savefile, 'el', 'nd', 'nd_ll', 'origin')




