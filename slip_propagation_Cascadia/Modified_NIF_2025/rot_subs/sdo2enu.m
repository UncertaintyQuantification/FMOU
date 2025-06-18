function enu = sdo2enu(sdo,normalvec)

% sdo2enu Get the [E,N,U] components of a vector from [strike,dip,normal] components
%
%  enu = sdo2enu(sdo,normalvec)
%
%    sdo: [strike,dip,opening] components of a vector. m x 3 array where
%       each row corresponds to a vector. The m vectors are simply converted
%       independently.
%       Striking direction is positive, dipping direction is positive, the
%       same direction as normalvec is positive for opening.
%    normalvec: normal vector of a plane, 1x3 or 3x1 array
%       NOTE: you cannot have a horizontal plane (normalvec=[0,0,1])
%       because strike direction is ambiguous. Treat this case manually and
%       separately.
%
%    enu: [east,north,up] components of a vector. m x 3 array.
%
%  Strike, dip and normal directions are defined by the normal vector
%  of a plane (normalvec).
%  See help normal2dipstrike for sign (etc) conventions of dip and strike.
%
%  Example:
%    normalvec = [1 -1 0.1];
%    [dip,strike] = normal2dipstrike(normalvec)
%    enu = [1 0 0; 0 -1 1; 0 0 -1; 0.2 0.5 2]
%    sdo = enu2sdo(enu,normalvec)
%    enu2 = sdo2enu(sdo,normalvec)
%
%  See also enu2sdo, getmat_enu2sdo, dipstrike2normal, normal2dipstrike.
%
%  16 Jul 2010, Yo Fukushima, DPRI, Kyoto Univ.
%

%% ChangeLogs
%  16 Jul 2010: first creation


%% Get rotation matrix
R = inv(getmat_enu2sdo(normalvec));


%% Output
enu = R*sdo';
enu = enu';

