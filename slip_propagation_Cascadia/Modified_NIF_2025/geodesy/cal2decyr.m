function decyr=cal2decyr(cal,pivot)
%CAL2DECYR   Returns decimal years given a calendar date
%   DECYR=cal2decyr(CAL,PIVOT) where CAL is 
%   [year(2 or 4 digit), month, day, hour, minute, second] 
%   If CAL has less than 6 columns, the missing columns are 
%   assumed to be zero. PIVOT, if specified, is the earliest 
%   2-digit year to be considered as being in the 1900s. 
%   Default is 80.

%-------------------------------------------------------------------------------
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%
%   Oct 09, 2001  Andy Hooper           Add optional pivot input
%   Aug 20, 2001  Peter Cervelli        Standardized code
%   Unknown       Peter Cervelli		Original Code
%
%-------------------------------------------------------------------------------

if nargin<2
    pivot=80;
end

%Convert to four digit years and pad with zeros as necessary

    cal(:,1)=yy2yyyy(cal(:,1),pivot);
    cal(:,size(cal,2)+1:6)=0;

%Calculate decimal year

    decyr=cal(:,1)+(datenum(cal(:,1),cal(:,2),cal(:,3),cal(:,4),cal(:,5),cal(:,6))-datenum(cal(:,1),1,1))./(datenum(cal(:,1)+1,1,0)-datenum(cal(:,1),1,0));
