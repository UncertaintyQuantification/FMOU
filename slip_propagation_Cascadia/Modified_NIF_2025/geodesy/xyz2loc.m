function [xyz,t,T]=xyz2loc(stations,i)
%XYZ2LOC     [xyz,t,T]=xyz2loc(stations,i)
%            Transforms XYZ ECEF receiver coordinates to a local ENU system.
%
%This function takes receiver coordinates and best fits a plane through them to form
%a local horizontal.  The positive x-axis is aligned with local east. The origin of the 
%local system is the ith station (or, if 'i' is not specified, the closest station to the
%mean of the station coordinates
%
%INPUTS:    'stations' is a 3xn matrix of XYZ receiver coordinates (m).
%           'i' is the reference station.
%
%OUTPUTS:   'xyz' is a 3xn matrix of ENU receiver coordinates (m).
%           't' is the translation vector.
%           'T' is the transformation matrix.
%
%Note: xyz = T*(stations-t*ones(1,n))
%
%Peter Cervelli, June 14, 1997.

%Define constants

	n=size(stations,2);

%Best fit plane to stations	

	g=[stations(1:2,:)',ones(n,1)];
	P=g\stations(3,:)';

%Form orthonormal basis oriented at best fitting plane

	v(1,:)=[1,0,P(1)];
	v(2,:)=[0,1,P(2)];
	v(1,:)=v(1,:)/norm(v(1,:));
	v(2,:)=v(2,:)-(v(1,:)*v(2,:)')*v(1,:);
	v(2,:)=v(2,:)/norm(v(2,:));
	v(3,:)=cross(v(1,:),v(2,:));	
	v(3,:)=v(3,:)/norm(v(3,:));	
	
%Transform stations to new basis, calculate distance from plane
	
	xyz=v*stations;

%Pick best station for origin (if necessary)

	if nargin<2	
		m=mean(xyz,2);
		dm=xyz-m*ones(1,n);
		[null,i]=min(sqrt(sum(dm.^2)));
	end

%Translate origin

	t=xyz(:,i);
	xyz=xyz-t*ones(1,n);
	
%Form rotation matrix that will rotate positive x-axis to local east

	%Use latitudes and longitudes to calculate azimuth to farthest station from reference		

		[null,j]=max(sqrt(sum((xyz(1:2,:)-xyz(1:2,i)*ones(1,n)).^2)));

		% This was using the old version of xyz2llh.m:
		%[p0,l0]=xyz2llh(stations(1,i),stations(2,i),stations(3,i));		% returns [lat, lon] in RADIANS
		%[p1,l1]=xyz2llh(stations(1,j),stations(2,j),stations(3,j));	

	% **** Make compatible with new version of xyz2llh.m ****
                [llh_0]=xyz2llh([stations(1,i);stations(2,i);stations(3,i)]);             % returns [lon;lat;h] in DEGREES 
                [llh_1]=xyz2llh([stations(1,j);stations(2,j);stations(3,j)]);     

                p0 = llh_0(2,:) * pi / 180;
		l0 = llh_0(1,:) * pi / 180;

		p1 = llh_1(2,:) * pi / 180;
		l1 = llh_1(1,:) * pi / 180;
	% *******************************************************


		az=atan2(cos(p1).*sin(l1-l0),(cos(p0).*sin(p1)-sin(p0).*cos(p1).*cos(l1-l0)));

	%Compare to current local azimuth
		
		az2=atan2(xyz(1,j),xyz(2,j));
		e=az2-az;
	
	%Create rotation matrix and rotate

		R=[cos(e),-sin(e),0;sin(e),cos(e),0;0,0,1];
		xyz=R*xyz;

%Form transformation matrix

	T=R*v;

%Form translation vector

	t=inv(v)*t;	
