function [XYZ] = llh2xyz(llh,datum)
%LLH2XYZ  Calculates global cartesisan coordinates from longitude, latitude, and height.
%   XYZ=llh2xyz(llh,DATUM) calculates global cartestian coordinates
%   given the 3xn (n = number of coordinate triples) matrix LLH that contains
%   longitude (deg), latitude (deg), and height (m) on the ellipsoid
%   specified by DATUM. DATUM can either be a vector the first two elements
%   of which give da and df, or it can be a string containing the name of a
%   datum that is resolved by the function DATUMS function.
%
%   Note that longitude occupies the first row of LLH.
%
%   See DATUMS for more information on datum parameters.

%-------------------------------------------------------------------------------
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%
%   Aug 24, 2001  Peter Cervelli        Standardized code
%   Unknown       Peter Cervelli		Original Code
%
%-------------------------------------------------------------------------------

%Check input arguments

    if nargin==1
        da=0;
        df=0;
    else
        if strcmp(class(datum),'char')
            datum=datums(datum);
            if any(isnan(datum))
                error('Could not resolve datum name.')
            end
        end
        da=datum(1);
        df=datum(2);
    end
	
	if size(llh,1)~=3
        error('Input llh MUST be 3xn.')
	end           

%Ellipsoid parameters

	a = 6378137 - da;
	f = 1/298.257223563 - df;
	b = (1-f)*a;

%Convert degrees to radians

	phi = llh(2,:)*pi/180;
	lam = llh(1,:)*pi/180;
   
%Convert llh to xyz
   
	N = a^2./sqrt( a^2*cos(phi).^2 + b^2*sin(phi).^2);
	XYZ(1,:) = (N+llh(3,:)).*cos(phi).*cos(lam);
	XYZ(2,:) = (N+llh(3,:)).*cos(phi).*sin(lam);
	XYZ(3,:) = ( b^2*N/a^2 + llh(3,:) ).*sin(phi);