gp = struct();
% scratch directory?
%gp.scratch_dir = '/Users/nbartlow/NIF/scratch/';
%gp.scratch_dir = '/Users/liuxubo/Documents/MATLAB/Modified_NIF_2025/scratch/';
gp.scratch_dir = './scratch/';
%Main data structure?
gp.data_mat  = './example/ETS2011_data2'; 
%station secular velocities (Ev, Nv, Uv) and annual (*a*) and seminannual (*s*)
%sine (*1) and cosine (*2) terms
gp.secular_mat='./example/ETS2011_secular';
%setup file (green's functions, mesh)
gp.setup_mat   ='./example/ETS2011_setup4.mat';


%%%%%gp.setup_mat   ='./example/ETS2011setup_coarse.mat';
gp.save_dir    = './';
%this is the beginning of the file name.  It will end with various flags
%and the name of the setup file.  Must be unique.

gp.savefile   = strcat('./Cascadia_setup4_');


% Filter parameters
%Smoothing (psuedo observations) applie to just slip ('slip') or
%just slip-rate ('rate'), or both ('both').  'Both' is also the default
%option if any other input is given.
gp.smooth_mode = 'both';

%Weight up/down certain stations
% This is useful if you have reason to insist the model fit a certain
% station, for example one placed by itself in a key area.  Low values
% indicate more weight put on the selected data.
%station weighting factor, 1 for no weighting, lower for upweighting
gp.w=1;
%list of GPS stations to weight
gp.staweight=[];

%Frame correction?  0 for none, 3 for 3 component frame (common mode) correction.  No rotation implemented. 
gp.Nframe = 3;
% gp.Nframe = 0;

% Hyperparameters
%These should be selected by optimizing for maximum likelihood estimation (MLE).  
%This code is set up to do a grid search, see MLE section below.
%Sigma - weighting of data noise.  From MLE.
gp.sigma_hat   = 1;
%GPS random walk amplitude, m/sqrt(year)
gp.tau_GPS_hat     = 1e-3/sqrt(1); %1 mm/sqrt(year)
%gp.tau_GPS_hat     = 1e-10; %1 mm/sqrt(year)
%alpha, temporal smoothing.  Higher values = less smoothing.
gp.alpha_hat   = 60; 
% gamma, spatial smoothing.  Higher values = less smoothing.
gp.gamma_hat   = 0.00015;

%MLE SECTION
gp.mle=0; %1 for mle (maximum likelihood) grid search, 0 for normal inversions
%gp.mle=1; %1 for mle (maximum likelihood) grid search, 0 for normal inversions
%range for grid search.  I recomend log spacing with factors of 2.
%gp.alpha=[25 50 100];%[100, 200, 400]; 
gp.alpha=logspace(1,2, 5);
%gp.gamma=[0.005, 0.01, 0.02, 0.04];%[0.00625, 0.0125]; 
gp.gamma=[0.0002, 0.0003, 0.0005, 0.0006, 0.001];

% Initial variances
%Should be very low for slip and rate.  I suggest 1e-7 and 1e-12
gp.var0_slip  = 1.0e-4; % units m
gp.var0_rate  = 1.0e-15; % units m/year (double check)
%initial coodinates variance, in mm
%should be large to accomindate error in specified postiton
gp.var0_coord = 100000; 


%NON NEGATIVITY SECTION
% Parameters governing projection behavior:
% These flags turn on/off non-negativity constraints on slip rate.  This
% constraint is accomplised by projecting the negip-rate elements 
% of xhat (the predicted state) to 0.  ative slThis should be done in the forwards
% and backwards steps, however using only forward may be OK for testing.
% Projection will make the code run a little bit slower.
% Project x_hat in forward_nnls?
gp.proj_use_fwd = 1;
% Project x_hat in backward_nnls?
gp.proj_use_bwd = 1;
% W constrols how the non-negativity projection is carried out
% W=I projects all negative components to 0, and doesn't change the others.
% Wid = 0 uses (Sigma)^-1 weighting in the projection, so all components
% change, not just the negative one.
% Wid=1 (W=I) will run faster, but Wid=0 is usually a better data fit.
% Set W = I?
% gp.proj_Wid = 1; % This is from the author
 gp.proj_Wid = 0; % This will apply pdco

%FROM HERE DOWN, DO NOT CHANGE
%full save file name
gp.svfn = sprintf('%s%sf%db%dW%d_%s.mat',...
		       gp.save_dir,gp.savefile,gp.proj_use_fwd,gp.proj_use_bwd,...
		       gp.proj_Wid);

%minimum norm instead of Laplacian smoothing?  Generally this is a bad
%idea, keep it at 0.
gp.minnorm=0;

%Strong constrains mean non-negativity is applied to slip.  You can't do
%that because the filter is unable to correct if too much slip comes in the
%model due to data noise.  Soft constraints are non-negativity applied to
%slip-rate only.  Since slip is only the integral of slip rate in the
%predicted state, but not in the final data-updated state, there can still
%be some small negative slip.  However the slip rate is always positive.
% Weak constraints only?
%Strong constraints don't work.  This must be 1.
gp.proj_weak = 1;

%include psuedo-obsertations (for smoothing) in max likelihood
%This should normally be 0.
gp.like_inc_pseudo = 0; % 0 to exclude them as usual
