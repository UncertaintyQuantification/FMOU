function doy=decyr2doy(decyr)
%DECYR2DOY    doy=decyr2doy(decyr)
%
%Converts decimal years to days of year.
%

% May 2, 2005       JRM     added comment about +1
% Sep 14, 2007      EKD     added .242199 to days in year
%     Removed ... 9/21

%daysinyear=isleapyear(decyr)+365.242199;
daysinyear=isleapyear(decyr)+365;

% Need to have +1 in following line because first day of the year is
% actually a fraction of a day until midnight.  For instance, noon on the
% first day would be doy 0.5.
doy=str2num(sprintf('%3.6f\n',(decyr-floor(decyr)).*daysinyear))+1;
