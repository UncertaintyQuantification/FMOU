function normal = dipstrike2normal(dip,strike)

% dipstrike2normal Get normal vector of a plane from dip and strike
%
%   normal = dipstrike2normal(dip,strike)
%
%     dip: dip angle in degrees, 0<=dip<=180
%     strike: strike angle in degrees, no restriction on the range.
%             The given value is wrapped between 0 and 2*pi.
%     normal: unit normal vector of the plane, 1 x 3 array
%
%   See help normal2dipstrike for the sign (etc) convention.
%     
%   Example:
%     normal = [-1 0.1 2]
%     normal = normal./norm(normal)
%     [dip,strike] = normal2dipstrike(normal)
%     normal = dipstrike2normal(dip,strike)
%
%   See also normal2dipstrike, dipstrike2norm.
%
%   16 Jul 2010, Yo Fukushima, DPRI, Kyoto Univ.
%

%% ChangeLogs
%  16 Jul 2010: first creation


%% Check
if ~(0<=dip<=180)
    error('dip out of bounds.');
end


%%
a0 = cosd(strike);
b0 = -sind(strike);
c0 = tand(90-dip);

normal = [a0,b0,c0];
normal = normal./norm(normal);

return




