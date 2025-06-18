function GF_setup(inputmesh, slip_dir, nu, lat_G, lon_G, h_G, ID_G, savefile)

%% NIF Green's Function setup script
%% modified by Noel Bartlow for inclusion in NIF codes (June 2015) from:
%% Script to analyze Loma Prieta post seismic data
%% Including spatial variation in slip
%% Paul Segall July-August 1997
%%
%% When you use this code with triangular Green's functions (tridisloc3d) you must cite:
% Thomas, A. L. (1993), Poly3D: A three-dimensional, polygonal element, 
% displacement discontinuity boundary element computer program with 
% applications to fractures, faults, and cavities in the Earth's crust, 
% M.S. thesis, Dep. of Geol. and Environ. Sci., Stanford Univ., Stanford, Calif.
%
% inputs: 
%
%   inputmesh: this is a .mat file containing triangular mesh information, as
% output by create_fault_mesh.m
%
%   slip_dir: a vector containg the direction of fault slip on each triangular
% patch, projected on to a horizontal plane, in degrees clockwise from
% north. Essentially, the local plate convergence direction.  I suggest
% getting this input from the MORVEL calculator at http://geoscience.wisc.edu/~chuck/MORVEL/motionframe_mrvl.html
%
%   nu: Poisson's ratio, something near 0.25
% 
%   lat_G, lon_G, h_G: the latitude, longitude, and height coordinates of every station for which data will be input into the NIF
%
%   ID_G: a list of 4 character station IDs.  This is not used for Green's
%   function generation, but needs to be included in the saved setup file,
%   and must match the 4 character IDs used in the NIF data input
%   structure.
%     
%   savefile: name of output .mat file

%% 


%-(1). 
    % Read in plate boundary mesh 
    % choice of origin read from mesh file
    load(inputmesh)
    
for i=1:length(el)
    trilat(i)=mean(nd_ll(el(i, :),2));
    trilon(i)=mean(nd_ll(el(i, :),1));
end
  

%%

    
    %solve for strike and dip of each element of mesh - they vary slightly
    for i = 1:size(el,1)
    norm_vec(i, :) = -getnormal(nd(el(i,1),:),nd(el(i,2),:),nd(el(i,3),:));
    [dip(i),strike(i)] = normal2dipstrike(norm_vec(i, :));
    end
  %  strike=strike-180; %relative to north, this is correct
    alpha = slip_dir-strike; 
    ss_fact = -cosd(alpha);
    ds_fact = sind(alpha);
    %%
    
    
	llh_G = [lat_G, lon_G, h_G]';

                                                                                         
	xy_obs_GPS = llh2localxy(llh_G(:,:),origin);

	Nsites = size(xy_obs_GPS,1);


%-(2). set constants and hyper parameters -------------------------

	deg2rad = pi/180;

        %shear modulus in Pa, completely does not matter for displacement or strain, just for stress 
        MU=3e10;

	Nframe = 3;


%-(?). Generating Green's functions ---------------------------
%tridisloc3d is a traingular dislocation code
%originates from the MS thesis of Andy Thomas, Stanford University, 1993
%PLEASE CITE

	xyz_obs_GPS = [xy_obs_GPS zeros(size(xy_obs_GPS,1),1)];

	for i = 1:length(el)
	    % strike slip
		[Gss(:,:,i), D, S] = tridisloc3d(xyz_obs_GPS', ...
				nd', el(i,:)', [1 0 0]', MU, nu);
	    % dip slip
		[Gds(:,:,i), D, S] = tridisloc3d(xyz_obs_GPS', ...
				nd', el(i,:)', [0 1 0]', MU, nu);
    end
    

% to get dipslip (ds) and strike slip (ss) components do this:
	Gds_e = squeeze(Gds(1,:,:));
	Gds_n = squeeze(Gds(2,:,:));
	Gds_u = squeeze(Gds(3,:,:));

	Gss_e = squeeze(Gss(1,:,:));
	Gss_n = squeeze(Gss(2,:,:));
	Gss_u = squeeze(Gss(3,:,:));


    %%
    
    %1 for strike slip component, 2 for dip slip component
	G1 = []; G2 = [];
	for is = 1:Nsites
		G1 = [G1; Gss_e(is,:); Gss_n(is,:); Gss_u(is,:)];
		G2 = [G2; Gds_e(is,:); Gds_n(is,:); Gds_u(is,:)];
    end
    

	Nbasis1 = size(G1,2);


	Kern_GPS = G1.*repmat(ss_fact, Nsites*3, 1) + G2.*repmat(ds_fact, Nsites*3, 1);

     
%- Laplacian

	[Lap, Lap_inv] = tri_laplacian(nd, el, 0);

	disp('Kernel: done')


%-(7). Saving setup environment ------------------------------------

	Ncomps  = size(Kern_GPS,1);

	save(savefile, 'el', 'ID_G', 'Kern_GPS', 'G1', 'G2', 'Lap', 'lat_G', 'lon_G', 'Nbasis1', 'Ncomps', 'nd_ll', 'origin')
    save('example/debugfile')



