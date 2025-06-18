function	[ts3]=mergegpsts(ts1,ts2);
%
% merges two gps time series,  according to the merge function in the 
% ident toolbox, there is no smart way to do this other than loadding all of the
% two objects data into the arrays and then calling gpsts to create a new object.
%
%

n1=length(ts1.d); n2=length(ts2.d);
d=[ts1.d; ts2.d;];

% The data covariance is why this is tricky.  not enough memory to
% just concatenate the arrays, so we would be basically redoing all the
% loops in readGPSTimeSeries3.m to make a cell array and then do cell2blkdiag
% Maybe there's no faster way to do this because of the memory limits for big arrays?
dcov=[ts1.dcov, zeros(n1,n2);
      zeros(n2,n1), ts2.dcov]);

sites=[ts1.sites, ts2.sites];
epochs=[ts1.epochs,ts2.epochs];
epochindex=[ts1.epochindex,ts2.epochindex];
apcoords=[ts1.apcoords,ts2.apcoords];

masterlist siteindex

ts3=gpsts(d,dcov,cellstr(masterlist),epochs,siteindex,cat(1,epochindex{:}),apcoords);

