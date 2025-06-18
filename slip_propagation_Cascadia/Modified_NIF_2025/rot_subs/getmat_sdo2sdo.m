function R = getmat_sdo2sdo(normalvec1,normalvec2)

% getmat_sdo2sdo Get the rotation matrix that does sdo2sdo
%
%  R = getmat_sdo2sdo(normalvec1,normalvec2)
%
%    normalvec1: normal vector of one plane, 1x3 or 3x1 array
%       NOTE: you cannot have a horizontal plane (normalvec=[0,0,1])
%       because strike direction is ambiguous. Treat this case manually and
%       separately.
%
%    normalvec2: normal vector of the second plane, 1x3 or 3x1 array
%
%    R: rotation matrix, 3x3 array
%
%  Strike, dip and normal directions are defined by the normal vector
%  of a plane (normalvec).
%  See help normal2dipstrike for sign (etc) conventions of dip and strike.
%
%  See also sdo2enu, enu2sdo, axisrot, dipstrike2normal, normal2dipstrike.
%
%  16 Jul 2010, Yo Fukushima, DPRI, Kyoto Univ.
%

%% ChangeLogs
%  16 Jul 2010: first creation


%% Get rotation matrix for converting sdo1 to enu
R1 = inv(getmat_enu2sdo(normalvec1));

%% Get rotation matrix for converting enu to sdo2
R2 = getmat_enu2sdo(normalvec2);

%% Output
R = R2*R1;

return

