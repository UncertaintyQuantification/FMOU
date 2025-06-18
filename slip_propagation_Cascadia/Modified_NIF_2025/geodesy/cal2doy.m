function doy=cal2doy(cal)
%CAL2DOY   doy=cal2doy(cal)  
%
%Returns day of year given a calendar date.
%Input 'cal' should be: [year month day [hr min sec]].
%(i.e., the dates are stored in the rows of 'cal'.)


if size(cal,2) == 6
    doy=datenum(cal(:,1),cal(:,2),cal(:,3),cal(:,4),cal(:,5),cal(:,6))-datenum(cal(:,1),1,0);
else 
    doy=datenum(cal(:,1),cal(:,2),cal(:,3))-datenum(cal(:,1),1,0);
end   
