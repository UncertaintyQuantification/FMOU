function [dip,strike] = normal2dipstrike(normal)

% NORMAL2DIPSTRIKE Convert from a normal vector of a plane to
%    dip and strike angles (in degrees).
%
%    [dip,strike] = normal2dipstrike(normal);
%
%    normal: 1x3 or 3x1 vector, a normal vector of a plane (norm of the
%            vector can be anything)
%    dip:    dip angle of the plane in degrees, 0<=dip<=180
%    strike: strike angle of the plane in degrees, 0<=strike<360
%
%    Sign (etc) conventions:
%    dip = 0,180 and dip = 90 correspond to horizontal and vertical planes,
%    respectively.
%    strike = 0 corresponds to north. strike is measured clockwise.
%    When the normal vector is heading upward (3rd component of normal is
%    positive), then 0<dip<=90, and strike is defined in such a way that 
%    the hanging wall block is on the right as viewed by an observer 
%    looking along strike (Aki and Richards convention).
%    When the normal vector is heading downward (3rd component of normal is
%    negative), then 90<dip<180, and strike would have 180 deg difference
%    compared to the Aki and Richards convention. This convention allows
%    in strike varying smoothly for dip around 90 deg.
%    When the 1st and 2nd components of normal are simultaneously zero
%    (normal vector is vertical), then strike will be nan.
%
%    Example:
%       normal = [1,1,1];
%       [dip,strike] = normal2dipstrike(normal)
%
%    See also dipstrike2normal, norm2dipstrike.
%
%    16 Jul 2010, Yo Fukushima, DPRI, Kyoto Univ.
%

%% ChangeLogs
%  16 Jul 2010: first creation by modifying norm2dipstrike


%%
a = normal(1);
b = normal(2);
c = normal(3);

if a==0 & b==0 % case the plane is z = const
    
    strike = NaN;
    dip = 0;
    return
    
else
    
    %% dip
    if c >= 0
        dip = atand(sqrt(a.^2+b.^2)./c);
    else
        dip = atand(sqrt(a.^2+b.^2)./c) + 180;        
    end
    
    
    %% strike
    if a == 0
        if b > 0
            strike = 270;
        else
            strike = 90;
        end
        return
    elseif a > 0 & b > 0 % EN quadrant
        strike = -atand(b/a) + 360;
    elseif a > 0 & b <= 0 % SE quadrant
        strike = -atand(b/a);
    elseif a < 0 % SW & NW quadrants
        strike = -atand(b/a) + 180;
    else
        error('something wrong in the code or your input...');
    end
    
end

return

