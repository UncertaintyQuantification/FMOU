function [RS]=RemoveStation(numstat,stat)
%RemoveStation   [RS]=RemoveStation(numstat,stat)
%
%Creates a matrix that gets rid of a station's three displacement
%components in a GPS data vector, such that d2=RS*d and dcov2=RS*dcov*RS'.

m=(numstat-1)*3;
n=numstat*3;
RS=eye(m,n);
i=(stat-1)*3+1;
RS(:,i+3:n)=RS(:,i:n-3);
RS(:,i:i+2)=zeros(m,3);