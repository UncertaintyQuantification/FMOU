function vec = planenormvec(p1,p2,p3)

% PLANENORMVEC Calculate the normal vector of a plane in 3D.
%
%   VEC = PLANENORMVEC(PT1,PT2,PT3) calculates the normal
%   vector of a plane that passes three points PT1, PT2,
%   PT3. When the plane is expressed as ax + by + cz = 1,
%   VEC(1) = a, VEC(2) = b, VEC(3) = c.
%   Points can be either a 1 x 3 or a 3 x 1 vector.
%
%   (Example)
%   p1 = [3,4,5];
%   p2 = [8,-4,0];
%   p3 = [0,0,1];
%   vec = planenormvec(p1,p2,p3);
%   x = -10:0.5:10;
%   y = -10:0.5:10;
%   [xi,yi] = ndgrid(x,y);
%   z = (1-vec(1)*xi-vec(2)*yi)/vec(3);
%   figure;
%   mesh(x,y,z');hidden off;hold on;
%   xlabel('x');ylabel('y');zlabel('z');
%   plot3(p1(1),p1(2),p1(3),'*b',p2(1),p2(2),p2(3),'*r',...
%      p3(1),p3(2),p3(3),'*k');
%   

% Yo Fukushima
%
% 26 Aug 2005, completely rewritten to use cross product
% 10 Mar 2005, first edition

% transpose if necessary
if size(p1,2) == 1
    p1 = p1';
end
if size(p2,2) == 1
    p2 = p2';
end
if size(p3,2) == 1
    p3 = p3';
end

p = cross(p1-p3,p2-p3);
if p(1)==0 & p(2)==0 & p(3)==0
    vec = [0 0 0];
    disp([mfilename ': Three points are on a line.']);
    return
end
f = p(1)*p1(1)+p(2)*p1(2)+p(3)*p1(3);
if f == 0
    vec = p./norm(p);
else
    vec = p./f;
end

return


%% OLDER VERSION %%
% % transpose if necessary
% if size(p1,2) == 1
%     p1 = p1';
% end
% if size(p2,2) == 1
%     p2 = p2';
% end
% if size(p3,2) == 1
%     p3 = p3';
% end
%
% n = origin(p1,p2,p3);
% if origin(p1,p2,p3) ~= 0 % when it passes [0,0,0]
%     if n == 1
%         vec = calcul(p2,p3);
%     elseif n == 2
%         vec = calcul(p1,p3);
%     else
%         vec = calcul(p1,p2);
%     end
%     return
% end
% 
% %% when it does not include [0 0 0]
% 
% % main part
% A = [p1;p2;p3];
% if cond(A) ~= inf % when not singular
%     d = [1;1;1];
%     vec = inv(A)*d;
% else
%     error('Three points on a line: badly conditioned.');
%     return
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%
% 
% function val = origin(p1,p2,p3)
% 
% % if p1=[0 0 0], it returns 1.
% % if p2=[0 0 0], it returns 2.
% % if p3=[0 0 0], it returns 3.
% % else it returns 0.
% 
% if (p1(1)==0 & p1(2)==0 & p1(3)==0)
%     val = 1;
% elseif (p2(1)==0 & p2(2)==0 & p2(3)==0)
%     val = 2;
% elseif (p3(1)==0 & p3(2)==0 & p3(3)==0)
%     val = 3;
% else
%     val = 0;
% end
% 
% return
% 
% %%%%%%%%%%%%%%%%%%%%%%%
% 
% function vec = calcul(pt1,pt2)
% 
% if pt1(1)==0 & pt1(2)==0
%     vec = [-pt2(2)/pt2(1),1,0];
% elseif pt2(1)==0 & pt2(2)==0
%     vec = [-pt1(2)/pt1(1),1,0];
% else
%     B = [pt1(1),pt1(2);pt2(1),pt2(2)];
%     if cond(B) ~= inf
%         d = [-pt1(3);-pt2(3)];
%         foo = inv(B)*d;
%         vec = [foo(1),foo(2),1];
%     else
%         error('Three points on a line: badly conditioned.');
%     end
% end
% return
% 
