% This is used to run the Modified NIF for real data
% Super slow. I already saved the results in the folder
% directly extract the results.

% try to get the NIF package set up and run the tutorial
addpath('./objects');
addpath('./YoMesh');
addpath('./geodesy');
addpath('./tridisloc3d');
addpath('./rot_subs');
addpath('./pdco-master/code');
addpath('./NIF_functions');

load GreyNegMapbig

% Optional codes to create input files:
% To create the setup file:
% 1) First run create_fault_mesh.m. 
create_fault_mesh('example/Cascadia_contours2012.txt', 20, 43.5, 49.5, -130, -110, -10, -70, 'example/myfault.mat');

% 2) Now run GF_setup.m to create the elastic Green's functions. This file 
% assumes a homogenous elastic half space and a triangular mesh input as output by create_fault_mesh.m

load example/GPS_info.mat;

gparms; % this will return a structure called gp
disp(gp);

NIFmain(gp); 


disp(gp.svfn);
load Cascadia_setup4_f1b1W0_.mat; % Users can directly call the file to load all results

% load example/coastline.dat;
% slip_hist(10, 70, 8, 5, [-125 -122], [44.5 49.4], gp.svfn, coastline);
% plotlist=['P374'; 'P376'; 'NWBG'; 'P409'; 'P417'; 'P430'; 'ELSR'; 'P436'];
% GPS_fit_plot(plotlist, gp.svfn)

% disp(size(slip_b));
% disp(size(rate_b));


%  csvwrite('saved_data/est_frame_large_mesh.csv',Trans);
%  csvwrite('saved_data/Green_fun_large_mesh.csv',Kern_GPS);
%  csvwrite('saved_data/est_slip_large_mesh.csv',slip_b);
%  csvwrite('saved_data/est_slip_rate_large_mesh.csv',rate_b);
%  csvwrite('saved_data/nd_ll_large_mesh.csv',nd_ll);
%  csvwrite('saved_data/el_large_mesh.csv',el);

