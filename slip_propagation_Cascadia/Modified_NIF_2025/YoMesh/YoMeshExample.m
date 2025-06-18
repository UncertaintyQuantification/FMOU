%   (Example 1)
nodet = [1000,1000,-1000;1100,1050,-950;...
    1180,1120,-900;1250,1220,-870];
nodeb = [777,635,-2000;963,789,-1960;1080,888,-1900;...
    1186,955,-1870;1338,1076,-1850;...
    1476,1154,-1800;1628,1250,-1770];
intv = 100;
[nd,el] = meshfrac(nodet,nodeb,intv);
figure;trisurf(el,nd(:,1),nd(:,2),nd(:,3));axis equal;
xlabel('x');ylabel('y');zlabel('z');

%   (Example 2: introducing vertical curvature)
nodet = [1000,1000,-1000;1100,1050,-950;...
    1180,1120,-900;1250,1220,-870];
nodeb = [777,635,-2000;963,789,-1960;1080,888,-1900;...
    1186,955,-1870;1338,1076,-1850;...
    1476,1154,-1800;1628,1250,-1770];
intv = 100;
curvaxispts = [0,0,0;1000,1500,0];
vertcurv = 30;
[nd,el] = meshfrac(nodet,nodeb,intv,curvaxispts,vertcurv);
figure;trisurf(el,nd(:,1),nd(:,2),nd(:,3));axis equal;
xlabel('x');ylabel('y');zlabel('z');


%  Load Hawaii map and geometries
[origin] = [-155.273256184167 19.339057464198 ]; % Site of MANE
load hawaii_line
line_xy = llh2local([coast(:,1)'; coast(:,2)'; zeros(1,length(coast))],origin)*1000;
plot(line_xy(1,:),line_xy(2,:),'k'); hold on
Kilauea_geometries; % load geometries of Kilauea tectonic features

% Hawaii map -- Hilinas
nodet = [hilina zeros(size(hilina,1),1)];
x = [-10:1:20]'*1000; y = x*sind(30)-10000;
nodeb = [x y -6000*ones(length(x),1)];
intv = 1000;
[nd,el] = meshfrac(nodet,nodeb,intv);
h = trisurf(el,nd(:,1),nd(:,2),nd(:,3)); axis equal; hold on
% set(h,'EdgeColor','none'); 
plot(line_xy(1,:), line_xy(2,:), 'k')
xlabel('x');ylabel('y');zlabel('z');

%  Hilinas with curve
curvaxispts = [x(1) x(end) 0; y(1) y(end) 0];
vertcurv = 60;
[nd,el] = meshfrac(nodet,nodeb,intv,curvaxispts,vertcurv);
h = trisurf(el,nd(:,1),nd(:,2),nd(:,3)); axis equal; hold on
% set(h,'EdgeColor','none'); 
plot(line_xy(1,:), line_xy(2,:), 'k')
xlabel('x');ylabel('y');zlabel('z');

% East Rift 
junk = ones(size(east_rift,1),1);
nodet = [east_rift 0*junk];
nodeb = [east_rift+[0*junk -4000*junk] -9000*junk];
intv = 1000;
curvaxispts = [east_rift 0*junk];
vertcurv = 20;
[nd,el] = meshfrac(nodet,nodeb,intv,curvaxispts,vertcurv);
h = trisurf(el,nd(:,1),nd(:,2),nd(:,3)); axis equal; hold on
% set(h,'EdgeColor','none'); 
plot(line_xy(1,:), line_xy(2,:), 'k')
xlabel('x');ylabel('y');zlabel('z');

% Decollement
junk = ones(size([SW_rift; east_rift],1),1);
nodet = [[flipud(SW_rift); east_rift] 0*junk ];
x = [-20:1:50]'*1000; y = x*sind(30)-20000;
nodeb = [x y -9000*ones(length(x),1)];
intv = 1000;
[nd,el] = meshfrac(nodet,nodeb,intv);
% figure;h = trisurf(el,nd(:,1),nd(:,2),0*nd(:,3)-9000); axis equal; hold on
h = trisurf(el,nd(:,1),nd(:,2),nd(:,3)-9000); axis equal; hold on
% set(h,'EdgeColor','none'); 
plot(line_xy(1,:), line_xy(2,:), 'k')
xlabel('x');ylabel('y');zlabel('z');
