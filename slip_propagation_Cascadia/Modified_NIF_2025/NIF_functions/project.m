function x = project(opts,x_hat,Sigma,lb,ub)
% x = project(s,x_hat,Sigma,lb,ub)
%   Project x_hat onto the bound constraints lb and ub. Solve
%     min 1/2 x' Sigma^{-1} x - x_hat' Sigma^{-1} x
%      st lb <= x < ub.
% Wid = 1 uses Sigma = I; the provided Sigma can be anything. Wid = 0
% uses the provided Sigma.

  fprintf(1,'qp start...');
  myt = tic;
  if (opts.Wid)
    % Explicit solution because the problem is separable
    x = min(ub,max(lb,x_hat));
  else
    % Use pdco to solve the nonseparable bound-constrained
    % problem. pdco will sometimes emit the warning:
    %    "Linesearch failed (nf too big)";
    % in this setting, this warning is benign. I'm setting the
    % tolerance to be very strict to be safe, and often it cannot be
    % fully achieved.
    
    R = chol(full(Sigma));
    x = qp_Wcov(R,x_hat,lb,ub,'pdco');
    fprintf(1,'||g||_2 = %e\n',Check(x,Sigma,x_hat,lb,ub));
  end
  fprintf(1,'...et = %f  ||x - x_hat||_2 = %e',...
	  toc(myt), norm(x - x_hat));
  if (~opts.Wid)
    fprintf(1,'  ||x - x(Wid=1)||_2 = %e\n',...
	    norm(x - min(ub,max(lb,x_hat))));
  else
    fprintf(1,'\n');
  end

function ng = Check(x,Sigma,x_hat,lb,ub)
  g = Sigma \ (x - x_hat);
  ng = norm(g(abs(x - lb) > 1e-2 & abs(x - ub) > 1e-2));
  