function R = getmat_enu2sdo(normalvec)

% getmat_enu2sdo Get the rotation matrix that does enu2sdo
%
%  R = getmat_enu2sdo(normalvec)
%
%    normalvec: normal vector of a plane, 1x3 or 3x1 array
%       NOTE: you cannot have a horizontal plane (normalvec=[0,0,1])
%       because strike direction is ambiguous. Treat this case manually and
%       separately.
%
%    R: rotation matrix, 3x3 array
%
%  Strike, dip and normal directions are defined by the normal vector
%  of a plane (normalvec).
%  See help normal2dipstrike for sign (etc) conventions of dip and strike.
%
%  If you want to get the rotation matrix that does sdo2enu, use simply
%  inv(R).
%
%  See also sdo2enu, enu2sdo, axisrot, dipstrike2normal, normal2dipstrike.
%
%  16 Jul 2010, Yo Fukushima, DPRI, Kyoto Univ.
%

%% ChangeLogs
%  16 Jul 2010: first creation


%% Algorithm
% I use axisrot.m. I first rotate around the strike direction to get the
% normal direction to coincide with the up direction. Then, I rotate around
% the up direction so that the strike and dip directions coincide with east
% and north components.


%% Settings
a = normalvec(1);
b = normalvec(2);
c = normalvec(3);

if normalvec(1) == 0 & normalvec(2) == 0 % for horizontal plane
    error('getmat_enu2sdo does not accept a vertical normal vector.');
end


%% Rotation around dip direction
[dip,strike] = normal2dipstrike(normalvec);
rotvec = [-b,a,0];
R1 = axisrot(rotvec,dip);


%% Rotation around strike direction
rotvec = [0,0,1];
if b == 0
    if a > 0
        ang = -90;
    else
        ang = 90;
    end
elseif a > 0 & b > 0 % EN quadrant
    ang = -atand(a/b)+180;
elseif a > 0 & b < 0 % SE quadrant
    ang = -atand(a/b);
elseif a < 0 & b < 0 % SW quadrant
    ang = -atand(a/b);
elseif a < 0 & b > 0 % NW quadrant
    ang = -atand(a/b)+180;
else
    error('something wrong in the code or your input...');
end

R2 = axisrot(rotvec,ang);


%% Output
R = R2*R1;

return

