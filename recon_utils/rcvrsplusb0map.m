function [Cgre,B0map] = rcvrsplusb0map(twix,twix_gre,noise_matrix,pulseq_gre,N)
% Estimation of coil sensitivity maps and static off-resonance maps     
% Inputs:
%   twix:           twix of imaging data
%   twix_gre:       twix of GRE data
%   noise_matrix:   Noise data
%   N:              Imaging matrix size
%
% Outputs:
%   Cgre:           Coil sensitivity maps
%   B0map:          Static off-resonance maps
%
% Copyright (c) 2026 Paul Dubovan and Xingwang Yong

    kspace = squeeze(twix_gre.image(:, :, :, 1, :, 1, 1, :, 1, 1, 1, 1, 1, 1, 1, 1));  %there are some dimensions that we dont care
    kspace = permute(kspace,[1,3,2,4,5]); %we reorder dimensions
    
    im = ifft2cmulti_norm(kspace); %IFFT to get images (these are actual GRE images)
    kspace = fft2cmulti_norm(fov_reduced(im)); %I removed FOV (oversampling 2X) and transfomred again to k-space.
    kfirst_echo = squeeze(kspace(:,:,:,:,1)); %We only need one echo for coil sensitivty, I typically use the first.
    ksecond_echo = squeeze(kspace(:,:,:,:,2));
    
    %% Interpolate k-space data to match imaging resolution
    if ~pulseq_gre 
        sz_preinterp = twix_gre.hdr.MeasYaps.sKSpace.lBaseResolution;
    else    
        sz_preinterp = twix_gre.image.NLin;
    end

    if exist('N')
        sz_postinterp = N;
    else
        sz_postinterp = twix.hdr.MeasYaps.sKSpace.lBaseResolution;
    end    

    if sz_postinterp > sz_preinterp
        pad_size = (sz_postinterp - sz_preinterp)/2;

        kfirst_echo = padarray(kfirst_echo, [pad_size pad_size], 'both');
        ksecond_echo = padarray(ksecond_echo, [pad_size pad_size], 'both');
    
        im_echo1 = ifft2cmulti_norm(kfirst_echo);
        im_echo2 = ifft2cmulti_norm(ksecond_echo);
    elseif sz_postinterp < sz_preinterp
        truncate_size = (sz_preinterp - sz_postinterp)/2;

        kfirst_echo = kfirst_echo(truncate_size+1:end-truncate_size,truncate_size+1:end-truncate_size,:,:);
        ksecond_echo = ksecond_echo(truncate_size+1:end-truncate_size,truncate_size+1:end-truncate_size,:,:);

        im_echo1 = ifft2cmulti_norm(kfirst_echo);
        im_echo2 = ifft2cmulti_norm(ksecond_echo);
    else
        im_echo1 = ifft2cmulti_norm(kfirst_echo);
        im_echo2 = ifft2cmulti_norm(ksecond_echo);
    end

    im_echo1 = permute(im_echo1,[1 2 4 3]);
    im_echo2 = permute(im_echo2,[1 2 4 3]);

    %%
    % Prewhiten the GRE data before sensitivity map estimation
    kfirst_echo = permute(kfirst_echo, [1 2 4 3]);
    kfirst_echo = my_prewhiten(kfirst_echo,noise_matrix);
    kfirst_echo = permute(kfirst_echo, [1 2 4 3]);

    ksecond_echo = permute(ksecond_echo, [1 2 4 3]);
    ksecond_echo = my_prewhiten(ksecond_echo,noise_matrix);
    ksecond_echo = permute(ksecond_echo, [1 2 4 3]);

    %% Coil estimation happens here. The function coil_espirit_2D loop goes over slices to call a coil senstivity estimator (ESPIRiT)
    eigThresh_2 = 0.8; %This threshold is like a mask. Higher means that coil sensivities take the whole background. Very low, small reduced support.
    Cgre = coil_espirit_2D_loop(permute(kfirst_echo,[1,2,4,3]),eigThresh_2); %there is anothe parameter to play inside this function, eigthres_1. It determines the maximum number
    Cgre = permute(Cgre, [1 2 4 3]);

    % echo times in ms
    TE_1 = twix_gre.hdr.MeasYaps.alTE{1}*10^-3;
    TE_2 = twix_gre.hdr.MeasYaps.alTE{2}*10^-3;
    
    im_echo1_cc = im_echo1.*conj(Cgre);
    im_echo1_cc = sum(im_echo1_cc,4);
    
    im_echo2_cc = im_echo2.*conj(Cgre);
    im_echo2_cc = sum(im_echo2_cc,4);

    Cgre = permute(Cgre, [1 2 4 3]);
    
    im_echo1_sos = sos(im_echo1,4);
    im_echo2_sos = sos(im_echo2,4);
    
    threshold = 0.1*max(im_echo1_sos(:))
    mask_v2 = im_echo1_sos > threshold;
    
    im_echo1_cc = im_echo1_cc.*mask_v2;
    im_echo2_cc = im_echo2_cc.*mask_v2;
    
    phi1 = angle(im_echo1_cc); % Phase of first echo
    phi2 = angle(im_echo2_cc); % Phase of second echo
    
    diffphs = phi2-phi1;
    
    %%
    % Other variation
    %im = angle(sum(conj(im_echo1).*im_echo2,4));
    %im = im.*imMask;

    Nslices = size(im_echo1,3);
    for Nsl = 1:Nslices
        phs_diff(:,:,Nsl) = unwrap_phase(diffphs(:,:,Nsl));
    end
        
    B0map = phs_diff/(TE_2-TE_1); %rad/ms

    % Requires smoothing, sigma value of at least 3 seems to provide enough
    % smoothing to eliminate noisy skull region
    B0map = imgaussfilt(B0map, 3.5); 

end

