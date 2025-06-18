	function [d_p, H_p]  = smooth_slip_rate2(statedim, ...
		Lap, Nbasis, smooth_mode, Npseudo, delta_t)

% 		 [d_p, H_p]  = smooth_rate(x, R1, Ntype, Nbasis1, Npseudo)
%
%  function to implement spatial smoothing as pseduo-observations
%  Note there is some question as to whether we should be smoothing
%  velocity or slip, or whether it in fact matters.
%  For now smooth slip
%
%  Version to omit the first basis function for each fault
% 	P.Segall 12/1/1997
%	corrected 3/1/1999, PSegall



%% pseudo-observation for the smoothing roughness

	d_p = zeros(Npseudo,1);

%% submatrix

	sub = zeros(Nbasis, statedim);	% smooth slip velocity
for i = 1:Nbasis
      if strcmp(smooth_mode,'slip') 
          sub(i,2*i-1) = (1);
      else
          if strcmp(smooth_mode,'rate')
               sub(i,2*i) = (1);  
          else %both case, also default
              %slip
               sub(i,2*i-1) = (1);
               %rate, additional factor to make it similar in magnitude to
               %slip
               sub(i,2*i) = delta_t;
          end
      end
end



%% Form H_p

	H_p = Lap*sub;