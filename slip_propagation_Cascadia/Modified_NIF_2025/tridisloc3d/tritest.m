function varargout = tritest(fn,varargin)
% Test tridisloc3d against disloc3d. To run as is, type dbg('Run'). To make
% more extensive comparisons, modify the code in 'cmp'.

%   AMB 01/2011
  [varargout{1:nargout}] = feval(fn,varargin{:});
    
% 1. localsoft and my version produce slightly different results (~1e-8
%    relative error, for example). localsoft and zipped mex file agree. so i
%    bet it's just a machine- or compiler-dependent compilation issue.
% 2. Third component of xy is in fact used, contrary to localsoft's
%    disloctest_mex.m doc.
  
function Run
  [t d Xo Yo] = cmp(1,1);
  cmpplots(t,d,Xo,Yo,0,4);

function [t d Xo Yo] = cmp(ewre,figno)
% Compare tridisloc3d and disloc3d.

  if (nargin < 2) figno = 1; end
  % element-wise relative error in plots?
  if (nargin < 1) ewre = 1; end  
  
  % Dislocation geometry
  dipa = 10;
  strikea = 30;
  deptha = -1000;
  cd = cosd(dipa);
  sd = sind(dipa);
  cs = cosd(strikea);
  ss = sind(strikea);
  xe = 1.5;
  ye = 1;
  ze = xe/cd*sd;
  % Slip. Strike seems quite good; dip seems worst.
  slip = [1 1 1];
  
  % Triangles to cover the square [0,1]x[0,1]
  x = [0 xe xe 0];
  y = [0 0 ye ye];
  xy = [x; y];
  Rs = [cs -ss; ss cs];
  xy = Rs*xy;
  x = xy(1,:);
  y = xy(2,:);  
  z = deptha + [-ze 0 0 -ze];
  tri = [1 2 4; 2 3 4]';
  nt = size(tri,2);

  % Obs
  no = 41;
  if (true)
    % Horizontal plane cutting halfway between bottom and top of dislocation
    Xo = linspace(-1,2,no);
    Yo = linspace(-1,2,no);
    [X Y] = meshgrid(Xo,Yo);
    xo = X(:);
    yo = Y(:);
    zo = mean(z)*ones(no^2,1);
    obs = [xo yo zo]';
  else
    % Horizontal plane just below and around the dislocation elements.
    Xo = linspace(-1,xe+1,no);
    Yo = linspace(-1,ye+1,no);
    [X Y] = meshgrid(Xo,Yo);
    alpha = X/xe;
    Z = -0.1 + (deptha - ze)*(1 - alpha) + deptha*alpha;
    Z(Z > 0) = 0;
    xy = [X(:) Y(:)]';
    xy = Rs*xy;
    obs = [xy' Z(:)]';
  end
  
  % disloc3d m vector
  length = ye;
  width = xe/cd;
  depth = -deptha;
  dip = -dipa;
  strike = -strikea;
  en = Rs*[xe ye/2]';
  east = en(1);
  north = en(2);
  ss = -slip(1);
  ds = slip(2);
  op = slip(3);
  m = [length width depth dip strike east north ss ds op]';
  
  % tridisloc3d inputs
  xyz = [x; y; z];
  comp = repmat(slip',1,nt);
  
  % Call
  mu = 1;
  nu = 0.25;
  tic;
  [d.U d.D d.S] = disloc3d(m,obs,mu,nu);
  fprintf(1,'d et = %f\n',toc);
  tic;
  [t.U t.D t.S] = tridisloc3d(obs,xyz,tri,comp,mu,nu);
  fprintf(1,'t et = %f\n',toc);
  tic;
  if(true)
    tic;  
    U = tridisloc3d(obs,xyz,tri,comp,mu,nu);
    fprintf(1,'t U et = %f\n',toc);
    tic;
    U = disloctest_mex(obs,xyz,tri,comp([2 1 3],:),nu);
    fprintf(1,'orig U et = %f\n',toc);
  end
  
  cmpplots(t,d,Xo,Yo,ewre,figno);
  
  A = {obs xyz tri comp mu nu};
  MatIo('WriteMats',A,'forvg.mats');
  
function cmpplots(t,d,xo,yo,ewre,figno)
  % Compare
  figno = figno - 1;
  figure(figno + 1); clf; cmpplot(d,t,'U',xo,yo,ewre);
  figure(figno + 2); clf; cmpplot(d,t,'D',xo,yo,ewre);
  figure(figno + 3); clf; cmpplot(d,t,'S',xo,yo,ewre);
  
function cmpplot(d,t,fld,xo,yo,ewre)
  n = size(d.(fld),1);
  switch(n)
   case 3, sm = 2; sn = 2;
   case 6, sm = 2; sn = 3;
   case 9, sm = 3; sn = 3;
  end
  if (ewre)
    fprintf(1,'%s median re:\n',fld);
    for (i = 1:n)
      sp(sm,sn,i);
      re = Re(d.(fld)(i,:),t.(fld)(i,:));
      imagesc(xo,yo,reshape(re,length(yo),length(xo)));
      axis equal; axis xy; axis tight; caxis([-16 0]); cb;
      fprintf(1,'%6.2f\n',median(re));
    end
  else
    for (i = 1:n)
      sp(2*sm,sn,i);
      imagesc(xo,yo,reshape(d.(fld)(i,:),length(xo),length(yo)));
      cb; axis equal; axis xy; axis tight;
      ca = caxis;
      sp(2*sm,sn,sm*sn+i);
      imagesc(xo,yo,reshape(t.(fld)(i,:),length(xo),length(yo)));
      cb; axis equal; axis xy; axis tight;
      caxis(ca);
    end
  end
  
function re = Re(a,b)
  % True r.e., but problematic for very small numbers
  %den = abs(b);
  % Error relative to max magnitude, but that changes with the problem.
  %den = max(abs(b(:)));
  % R.e. guadred by a quantity. Can be hard to interpret.
  %den = max(abs(b),sqrt(eps)*max(b(:)));
  % Our problem has max magnitudes of O(1), so let's go with that.
  den = 1;
  re = log10(abs((b - a)./den));
  
function ae = Ae(a,b)
  ae = abs(b - a);
