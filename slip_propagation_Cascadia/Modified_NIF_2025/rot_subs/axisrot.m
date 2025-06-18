function R = axisrot(rotvec,ang)

% AXISROT3D Create a matrix for rotation around a vector.
%
%   A = AXISROT3D(ROTVEC,ANG) creates a matrix for rotating ANG
%   degrees around direction specified by ROTVEC.
%
%   (Example)
%   rotvec = [1 1 1];
%   p = [1 0 0]';
%   figure;
%   for ang = 0:10:3600;
%     A = axisrot(rotvec,ang);
%     q = A*p;
%     plot3(q(1),q(2),q(3),'.');
%     axis([-1 2 -1 2 -1 2]);axis square;
%     xlabel('x');ylabel('y');zlabel('z');
%     view([-1 -1 1]);
%     drawnow;
%     pause(0.1);
%   end
%
%   See also ROLLPITCHYAW.
%

% normalize %
if norm(rotvec) ~= 1
   rotvec = rotvec./norm(rotvec);
end

a = rotvec(1);
b = rotvec(2);
c = rotvec(3);
t = ang*pi/180;
R = ones(3,3);

R(1,1)=a^2*(1-cos(t))+cos(t);
R(2,2)=b^2*(1-cos(t))+cos(t);
R(3,3)=c^2*(1-cos(t))+cos(t);
R(1,2)=a*b*(1-cos(t))+c*sin(t);
R(2,1)=a*b*(1-cos(t))-c*sin(t);
R(1,3)=a*c*(1-cos(t))-b*sin(t);
R(3,1)=a*c*(1-cos(t))+b*sin(t);
R(2,3)=b*c*(1-cos(t))+a*sin(t);
R(3,2)=b*c*(1-cos(t))-a*sin(t);
