function varargout = qp_Wcov(R,x_hat,lb,ub,solver)
% x = qp_Wcov(R,x_hat,lb,ub)
%   Solve
%     min 1/2 x' H^{-1} x - x_hat' H^{-1} x
%      st lb <= x <= ub,
% where R = chol(H).

  % Handle call to SnoptFnBc
  if (nargin == 1)
    % x = R
    [varargout{1:nargout}] = SnoptFnBc(R);
    return;
  end
  
  if (nargin < 5) solver = 'pdco'; end
  
  x_hat = x_hat(:);
  lb = lb(:);
  nx = length(x_hat);
  ox = ones(nx,1);
  zx = zeros(nx,1);
  
  switch (solver)
   case 'pdco'
    pdObj = [-(R \ (x_hat' / R)'); zx];
    % v0.7: Previous formulation used inv. What was I thinking?
    pdMat = [eye(nx) -full(R)'];
    b = zeros(nx,1);
    bl = [lb; -inf*ox];
    bu = [ub;  inf*ox];
    d1 = [1e-4*ox; ox];
    d2 = 1e-12; % v0.7: This needs to be much tighter for the new pdMat.
    o = pdcoSet();
    o.OptTol = eps;
    o.FeaTol = o.OptTol;
    o.MaxIter = 60;
    o.Print = 0;
    o.Method = 1;
    
    x0 = [x_hat; (x_hat' / R)'];
    y0 = zx;
    z0 = [zx; zx];
    xsize = max(abs(x0));
    zsize = xsize;
    for (i = 1:10) % Allow up to 10 calls to pdco to solve this problem
      [x y0 z0 inform] = pdco(pdObj,pdMat,b,bl,bu,d1,d2,o,x0,y0,z0,xsize,zsize);
      if (inform ~= 1) break; end
      x0 = x;
      xsize = max(abs(x));
      zsize = max(abs(z0));
      if (i == 10 && inform == 1)
	warning('qp_Wcov: Did not converge.\n');
      end
    end
    fprintf(1,'norm(x - R''*y) = %e\n',norm(x(1:nx) - R'*x(nx+1:end)));
    x = x(1:nx);
    
   case 'snopt'
    addpath /work/ambrad/icme/Research/snopt7.2-9/snopt7/matlab/

    global sn;
    sn.nx = nx;
    y = (x_hat' / R)';
    sn.f = -(R \ y);
    x = x_hat;

    A = [zeros(1,2*nx); speye(nx) -R'];
    [iAfun jAvar] = find(A);
    A = full(A(sub2ind([nx+1 2*nx],iAfun,jAvar)));
    
    lb = [lb; -inf*ox];
    ub = [ub;  inf*ox];
    
    Flo = [-inf; zx];
    Fhi = [ inf; zx];
    
    iGfun = ones(2*nx,1);
    jGvar = (1:2*nx)';
    
    snprintfile 'qp_Wcov_snopt.txt';
    snset('Major feasibility tolerance 1e-12');
    snset('Minor feasibility tolerance 1e-12');
    [x F inform xmul Fmul] = snopt(...
	[x; y],lb,ub,Flo,Fhi,...
	func2str(@(x)qp_Wcov(x)),...
	A,iAfun,jAvar,iGfun,jGvar);
    fprintf(1,'norm(x - R''*y) = %e\n',norm(x(1:nx) - R'*x(nx+1:end)));
    x = x(1:nx);
    snprintfile off;
    clear global sn;
  end
  
  varargout{1} = x;
  
function [F G] = SnoptFnBc(x)
  global sn;
  y = x(sn.nx+1:end);
  x = x(1:sn.nx);
  F = [0.5*dot(y,y) + dot(sn.f,x);
       zeros(sn.nx,1)];
  G = [sn.f; y];
