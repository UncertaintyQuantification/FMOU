function [T,L,Lcov,B,Bcov]=baseline(ts,sta1,sta2)
%BASELINE  Returns baseline length and covariance given a time 
%          series object and two station names.
%    [T,L,Lcov]=baseline(ts,sta1,sta2) calculates baseline 
%    lengths between two stations.
%
%    This function uses a linearization that assumes that the 
%    baseline length only changes in response to movement in 
%    the direction between the two stations.
%
% v1.0, 04/16/2002, PFC.
%
% Note: this is not the same function as the baseline.m for
% use with the time series object.  JRM 4/21/03

%Find times when both stations appear

    T=intersect(ts(sta1,:).t,ts(sta2,:).t);

%Get data and covariance for the pair

    d=ts({sta1,sta2},{T}).d;
    dcov=ts({sta1,sta2},{T}).dcov;

%Form differencing operator

    n=length(d);
    I=(1:n/2)';
    J=repmat(0:3:(n/2-3),3,1);
    K=I+J(:);
    D=sparse(I,K,1,n/2,n)+sparse(I,K+3,-1,n/2,n);

%Form length operator

    coords=ts({sta1,sta2},:).apcoords;
    u=coords(:,1)-coords(:,2);
    u=u/norm(u);
    I=1:n/6;
    J=1:3:n/2;

M=sparse(I,J,u(1),n/6,n/2)+sparse(I,J+1,u(2),n/6,n/2)+sparse(I,J+2,u(3),n/6,n/2);

%Apply operators

    B=D*d;
    L=M*B;
    Bcov=D*dcov*D';
    Lcov=M*Bcov*M';
