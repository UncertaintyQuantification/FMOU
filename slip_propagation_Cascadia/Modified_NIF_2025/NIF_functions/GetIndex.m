      function [sta_index] = GetIndex(sta_list, obs_sta, flg)
%     function [sta_index] = GetIndex(sta_list, obs_sta, flg)
%
% searches through the vector of station names 'sta_list' 
% to find the one that matches the given station name 'obs_sta'.
% returns the  index number of the matching component of the 
% sta_list array.  


	[Nsites, four] = size(sta_list);
	sta_index = 0;

	tmp = 0;
	for k=1:Nsites
		if obs_sta == sta_list(k,:)
			tmp = tmp + 1;
			if tmp >= 2
				sta_index = [sta_index; k];
			else
				sta_index = k;
			end
		end
	end

%if nargin == 2
%	if    sta_index == 0
%		disp(['station  ',obs_sta,'  not found'])
%	end
%end



