function llhB=shiftdatum(llhA,datumA,datumB)
%SHIFTDATUM   Shifts coordinates from one datum to another.
%   LLH_B=shiftdatum(LLH_A,DATUM_A,DATUM_B) shifts the coordinates,
%   specified as longitude (deg), latitude (deg) and height (m),
%   from DATUM_A to DATUM_B.  LLH_A should be a 3xn matrix (n is the
%   number of coorinate triples).  DATUM_A and DATUM_B can either be
%   row vectors containing datum parameters or strings containing
%   the names of datums that are resolved by the function DATUMS.
%
%   The datum shift is performed as follows: LLH_A is converted to
%   XYZ_A using the ellipsoid of DATUM_A.  XYZ_B is calculated by adding
%   the difference between the geocenter of DATUM_A and DATUM_B to XYZ_A.
%   Then XYZ_B is converted to LLH_B using the ellipsoid of DATUM_B.
%
%   See DATUMS for more information on datum parameters.

%-------------------------------------------------------------------------------
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%
%   Aug 20, 2001  Peter Cervelli        Standardized code
%   Unknown       Peter Cervelli		Original Code
%
%-------------------------------------------------------------------------------

%Resolve datum names if necessary

    if strcmp(class(datumA),'char')
        datumA=datums(datumA);
    end

    if strcmp(class(datumB),'char')
        datumB=datums(datumB);
    end

    if any(isnan(datumA)) | any (isnan(datumB))
        error('Could not resolve datum name.')
    end

%Convert from LLH to XYZ in ellipsoid of datumA

    xyzA=llh2xyz(llhA,datumA);

%Shift geocenter

    delta=datumA(3:5)-datumB(3:5);
    xyzB=xyzA+repmat(delta(:),1,size(xyzA,2));

%Convert from XYZ to LLH on ellipsoid of datumB

    llhB=xyz2llh(xyzB,datumB);

