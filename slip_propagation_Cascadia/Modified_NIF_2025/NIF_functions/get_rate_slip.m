
function [slip_f, rate_f, slip_b, rate_b, ...
			sig_slip_f, sig_rate_f, ...
			sig_slip_b, sig_rate_b, ...
			Cslip_b, sig_Cslip_b] = ...
		get_rate_slip(x_f, covx_f, Kern_GPS, Nepochs)
   
%%%
    %load last_test
	Nbasis = size(Kern_GPS,2);

        slip_f = zeros(Nbasis, Nepochs);
        rate_f = zeros(Nbasis, Nepochs);
        slip_b = zeros(Nbasis, Nepochs);
        rate_b = zeros(Nbasis, Nepochs);
        sig_slip_f = zeros(Nbasis, Nepochs);
        sig_rate_f = zeros(Nbasis, Nepochs);
        sig_slip_b = zeros(Nbasis, Nepochs);
        sig_rate_b = zeros(Nbasis, Nepochs);
       

        for k = 1:Nepochs
            
            k
                slip_f(:,k) = x_f{k,k}(1:2:2*Nbasis-1);
                rate_f(:,k) = x_f{k,k}(2:2:2*Nbasis);
                slip_b(:,k) = x_f{k,Nepochs}(1:2:2*Nbasis-1);
                rate_b(:,k) = x_f{k,Nepochs}(2:2:2*Nbasis);
                sig_slip_f(:,k) = sqrt(diag(covx_f{k,k}(1:2:2*Nbasis-1,1:2:2*Nbasis-1)));
                sig_rate_f(:,k) = sqrt(diag(covx_f{k,k}(2:2:2*Nbasis,2:2:2*Nbasis)));
                sig_slip_b(:,k) = sqrt(diag(covx_f{k,Nepochs}(1:2:2*Nbasis-1,1:2:2*Nbasis-1)));
                sig_rate_b(:,k) = sqrt(diag(covx_f{k,Nepochs}(2:2:2*Nbasis,2:2:2*Nbasis)));

		tmp = Kern_GPS*covx_f{k,Nepochs}(1:2:2*Nbasis-1,1:2:2*Nbasis-1)*Kern_GPS';
		sig_Cslip_b(:,k) = sqrt(diag(tmp));
        end
       

        %calculate displacement from only negative part
        %slip_b(slip_b>0)=0;
        
        Cslip_b = Kern_GPS*slip_b;


