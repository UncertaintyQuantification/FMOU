function [KM]=keepstations(numstats,keepstats)
%keepstations    [KM]=keepstations(numstats,keepstats)
%
%Makes a matrix that keeps the stations whose indices are
%listed in keepstats, and discards the rest. numstats is
%the original total number of stations.
%
%Apply KM as follows:
%
%d_new = KM * d;
%
%dcov_new = KM * dcov * KM';

i=length(keepstats)*3;
keepstats=keepstats(:)';
J=(keepstats-1)*3+1;
KM=sparse(1:i,[J;J+1;J+2],1,i,numstats*3);
