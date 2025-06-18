function [e,ecov]=xyz2enu(d,dcov,origin)
%XYZ2ENU   Transforms from global cartestian to local cartesian.
%   [E,ECOV]=xyz2enu(D,DCOV,ORIGIN) transforms data vector D and
%   data covariance DCOV from a global cartesian (XYZ) coordinate
%   system to a local coordinate system aligned with the geographic
%   directions at the position ORIGIN.  D should be either 3nx1 or
%   3xn (n = number of individual vectors).  DCOV should be 3nx3n.
%   ORIGIN should be a vector of length 2 or 3.  If length 2, ORIGIN
%   is taken as a longitude, latitude pair (degrees); if length 3,
%   ORIGIN is taken as an XYZ station position. E is matrix (vector)
%   of transformed coordinates the same size as input D.  ECOV is a
%   matrix containing the transformed covariance.
%
%   E=xyz2enu(D,ORIGIN) behaves as above without a data covariance
%   matrix.

%-------------------------------------------------------------------------------
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%
%   Aug 24, 2001  Peter Cervelli        Standardized code
%   Nov 04, 2000  Peter Cervelli		Original Code
%
%-------------------------------------------------------------------------------

%Parse input arguments

    if nargin==2
        origin=dcov;
    end

    if length(origin)>2
        origin=xyz2llh(origin(:));
    end
   
    [i,j]=size(d);
    d=d(:);  

%Make transformation matrix

	Tm=xyz2enum(origin,length(d)/3);

%Transform

    e=reshape(Tm*d,i,j);
	if nargout==2 & nargin > 2
	    ecov=Tm*dcov*Tm';
	end

