function noise_matrix = get_noise_matrix(twix)

% This function prepares the k-space raw data into an ordered
% multi-dimensional matrix.
%final data dimensions: following Siemens PRS coordinate
%line-col-par-cha-sli-ave-phs-eco-rep-set-seg-ida-idb-idc-idd-ide

twix_obj = twix;
twix_obj_size = numel(twix_obj);
if 1%read noise data
    if iscell(twix_obj)
        noise_data = double(twix_obj{1}.noise());
    else
        noise_data = double(twix_obj.noise());
    end
    %noise_data is has a size of 512-cha-128-1-1-2
    noise_data = mean(noise_data,6);
    noise_data = permute(noise_data,[1,3,2]);%512*128*16
    noise_data = reshape(noise_data, [size(noise_data,1)*size(noise_data,2), size(noise_data,3)]);%65536*16
    noise_mean = mean(noise_data, 1);%1*16
    temp =  noise_data - ones(size(noise_data,1),1)* noise_mean; %65536*16
    % noise_matrix = (temp' * temp); % 16*65536 * 65536*16 here it's E[n1* x n2]
    a = temp.';
    % noise_norm = var(a, [], 2)*(size(a,2)-1);
    %scale = (mean(1./noise_norm));
    scale = size(a,2)-1;
    noise_matrix = (1/scale)*(a * a'); % 16*65536 * 65536*16 here it's E[n1 x n2*], correct one on 08/31/2014
end

