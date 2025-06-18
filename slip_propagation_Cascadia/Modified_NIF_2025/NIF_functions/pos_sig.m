        function [pos_adj, pos_adj_f, pos_sigma] ...
		= pos_sig(x_f, covx_f, CascadiaNorth, ...
                        delta_k, Nbasis, Nsitecomps)

%%
%%        function [pos_adj, pos_sigma] = pos_sig(t, x, ...
%%              Nbasis, Ncomps, files)
%%
%%  Reconstruct slip-slip and uncertainty at each epoch from
%%  state vector and covariance
%%


%%  Determine some constants

        Nepochs = size(CascadiaNorth,2);
        statedim = length(x_f(1,1));

        pos_adj = zeros(Nsitecomps, Nepochs);

        kk_last = size(x_f);
        Mepochs = kk_last(1);

        k_last = delta_k*(kk_last-1) + 1;

    for k = 1:delta_k:Nepochs

        kk = (k-1)/delta_k + 1;

        [k kk]

        p1 = 2*Nbasis;

        pos_adj(:,k) = x_f{kk,Mepochs}(p1+1:p1+Nsitecomps)';
        pos_adj_f(:,k) = x_f{kk,kk}(p1+1:p1+Nsitecomps)';

        tmpvar = covx_f{kk,Mepochs}(p1+1:p1+Nsitecomps,p1+1:p1+Nsitecomps);
        pos_sigma(:,k) = sqrt(diag(tmpvar));

    end

