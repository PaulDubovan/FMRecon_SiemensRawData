% B0map generation from GRE dicom data using spm toolbox function
% magnitude data (_e1.nii or _e2.nii file)
mag_fn = '';
% phase data (_ph.nii file)
phase_diff_fn = '';
% Supply used echo times [ms]
EchoTimes =  [e1 e2];
deltaTE = diff(EchoTimes);

%%

B0map_spm = calc_B0_map_spm(phase_diff_fn,mag_fn,EchoTimes);

%% higher res
new_size = [Nx Ny];
B0map_spm = imresize(B0map_spm, [new_size]);

%% lower res
kspace = fft2cmulti_norm(B0map_spm); %I removed FOV (oversampling 2X) and transfomred again to k-space.

truncate_size = (new_size-old_size)/2;

kspace = kspace(truncate_size+1:end-truncate_size,truncate_size+1:end-truncate_size,:,:);

B0map_spm = ifft2cmulti_norm(kspace);

