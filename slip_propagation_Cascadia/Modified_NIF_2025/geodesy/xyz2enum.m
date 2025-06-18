function T=xyz2enum(origin,n)
%XYZ2ENUM   Returns a global to local transformation matrix.
%   T=xyz2enum(ORIGIN) Returns a transformation matrix that
%   tranforms coordinates in a global ECEF cartesian system
%   into to a local coordinate system aligned with the geographic
%   directions at the position ORIGIN.  ORIGIN should contain a
%   longitude and a latitude pair (degrees). T is 3x3.
%
%   T=xyz2enum(ORIGIN,N) returns repeated block diagonally N times.

%-------------------------------------------------------------------------------
%   Record of revisions:
%
%   Date           Programmer            Description of Change
%   ====           ==========            =====================
%
%   Aug 24, 2001   Peter Cervelli        Standardized code
%   July 27, 2001  Peter Cervelli        Eliminated unnecessary memory allocation.
%                                        One element of the 3x3 transformation matrix
%                                        is zero.  Originally, this element still was allocated
%                                        as non-zero in the sparse block diagonal matrix, T.
%                                        This resulted in unnecessary memory usage.
%   April, 2000    Peter Cervelli		 Original Code
%
%-------------------------------------------------------------------------------

%Check input arguments

    if nargin < 2
        n=1;
    end
	
	if prod(size(origin))<2
        error('Input origin must have 2 elements, longitude and latitude (degrees).');
	end
	
	if prod(size(n))~=1
        error('Input n must be a scalar.');
	end

%Convert to radians and evaluate trigonometric functions

    origin=origin*pi/180;
    s=sin(origin);
    c=cos(origin);

%Make the (sparse) transformation matrix

    N=3*n;
    i=[1:3:N 2:3:N 3:3:N 1:3:N 2:3:N 3:3:N 2:3:N 3:3:N]; 
    j=[1:3:N 1:3:N 1:3:N 2:3:N 2:3:N 2:3:N 3:3:N 3:3:N];
    T=repmat([-s(1) -s(2)*c(1) c(2)*c(1) c(1) -s(2)*s(1) c(2)*s(1) c(2) s(2)],n,1);
    T=sparse(i,j,T(:));