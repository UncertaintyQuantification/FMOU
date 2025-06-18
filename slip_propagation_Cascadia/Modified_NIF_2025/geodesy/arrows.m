function [X,Y] = arrows(xy,d)
%--------------Arrows------------------%
% From DM_Quiver.m (Peter's code)
% function [X,Y] = arrows(xy,d)
% Inputs:
%  xy - station locations (ENU) (2 or 3 x n)
%  d - data vector (3n x 1)
%
%  Outputs:
%   X,Y - plot(X,Y) produces all requested arrows on a map


%---Arrow shafts---%
   ns = size(xy,2);
   
   U=d(1:3:end,end)'*1e2;
   V=d(2:3:end,end)'*1e2;
   X=[xy(1,:);xy(1,:)+U;repmat(NaN,1,ns)];
   Y=[xy(2,:);xy(2,:)+V;repmat(NaN,1,ns)];
%---Arrow heads---%
   alpha=0.2; beta=0.33; Up=U; Vp=V;
   L=sqrt(sum([X(1,:)-X(2,:);Y(1,:)-Y(2,:)].^2));
   I=find(L>3);
   Up(I)=Up(I)./(L(I)/3);
   Vp(I)=Vp(I)./(L(I)/3);
   X=[X;X(2,:)-alpha*(Up+beta*(Vp+eps));X(2,:);X(2,:)-alpha*(Up-beta*(Vp+eps));repmat(NaN,1,ns)];
   Y=[Y;Y(2,:)-alpha*(Vp-beta*(Up+eps));Y(2,:);Y(2,:)-alpha*(Vp+beta*(Up+eps));repmat(NaN,1,ns)];
