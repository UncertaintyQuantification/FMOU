
	function [TT, Trans] ...
		= get_frame(x_f, Nepochs, Nsites, Nbasis, Ncomps, Nframe);


        subRtr = [];  Trans = [];

        if Nframe ~= 0

                for i = 1:Nsites 
                        subRtr = [subRtr; eye(3)];
                end

                %for i = 1:nT
                %        subRtr = [subRtr; zeros(2,3)];
                %end

                for k = 1:Nepochs
                        p = 2*Nbasis + Ncomps;
                    %% XYZ components for network 
                        TT(:,k) = x_f{k,Nepochs}(p+1:p+3)';
                    %% ENU components for network 
                        Trans = [Trans, subRtr*TT(:,k)];
                end

        end

