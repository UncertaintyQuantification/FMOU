function cal=sec2cal(sec)
%SEC2CAL  Returns calendar date given seconds before 2000/01/01 12:00.
%     cal=sec2cal(sec)

%This version is slightly slower but more numerically stable than
%the original: cal=datevec(sec/86400+730486.5)

%Make sure sec is a column vector

    sec=sec(:);
    
%Calculate calendar day

    fsec=floor(sec);
    epsilon=sec-fsec;
    
    fsec=fsec+63114033600;
    rs=mod(fsec,86400);
    rh=mod(rs,3600);
    rm=mod(rh,60);
    
    cal=datevec((fsec-rs)/86400);
    cal(:,4)=(rs-rh)/3600;
    cal(:,5)=(rh-rm)/60;
    cal(:,6)=rm+round(epsilon*1e8)/1e8;
