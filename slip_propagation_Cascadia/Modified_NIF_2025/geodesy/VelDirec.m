function [dir2]=VelDirec(pole_lat, pole_lon, lat, lon)
%convert Euler poles to velocity directions
%approximates a spherical earth
%output is [East; North] components of direction

%lon, lat, height
pole_llh=[pole_lon; pole_lat; 0];
obs_llh=[lon; lat; 0];

pole_XYZ=llh2xyz(pole_llh);
obs_XYZ=llh2xyz(obs_llh);

dir=cross(pole_XYZ, obs_XYZ-pole_XYZ);
dir=dir/norm(dir);

%unit north at observation
north=[0, 0, 1];
unitE=cross(north, obs_XYZ);
unitE=unitE/norm(unitE);
unitN=cross(obs_XYZ/norm(obs_XYZ), unitE);

dirN=dot(dir, unitN);
dirE=dot(dir, unitE);

dir2=[dirE; dirN];

end