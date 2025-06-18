% TRIDISLOC3D    [U D S] = tridisloc3d(XYZ, ND, EL, COMP, MU, NU);
%
%  Outputs displacement, and optionally displacement derivatives and stress, at
%  observation coordinates resulting from slip or opening on triangular
%  fault patches.
% 
%   Inputs:
%   XYZ - Observation coordinates, must be 3xn. n = number of observation
%     coordinates. coordinates are [E; N; U].
%   ND - Locations of triangular fault patch nodes 3xNnodes. [E; N; U].
%   EL - Indices into nd defining which three node points form each
%     triangular fault patch 3xNpatches. [nd1; nd2; nd3].
%   COMP - A vector indicating the components of slip allowed on each patch.
%     Dimensions are the samee as for EL. Each column contains [SS; DS; OP],
%     which is the same ordering convention as in disloc3d.
%   MU - Shear modulus. Set to an arbitrary scalar if S is not requested.
%   NU - Poisson's ratio (usually about 0.25).
% 
%  The coordinate system is as follows: 
%    east = positive X, and north = positive Y.
%  Depths should be given negative for both XY and ND (in contrast to the
%  mixed convention in disloc3d).
%
%   Output: 
%     U - 3xn. Each column contains the displacement [east; north; up].
%   (optional: greater computation time)
%     D - 9xn. Each column contains the 9 displacement spatial derivatives
%       Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz.
%     S - 6xn. Each column contains the 6 independent stress components
%       Sxx, Sxy, Sxz, Syy, Syz, and Szz.
% 
%   The units of the displacements will be the units of the slip (opening); 
%   the units of the observation coordinates should be consistent with 
%   the units of the model.
 