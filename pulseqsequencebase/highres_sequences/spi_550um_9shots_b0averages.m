% this is an experimental spiral sequence

addpath(genpath('/path to vdspiral'))
addpath(genpath('/path to safe_pns_prediction'))
addpath(genpath('path to pulseq v142, modifications for v151 and up in progress'))
addpath(genpath('FMRecon_SiemensRawData/utils'))

clear
close all


%% input
isTestRun = 0; % for test only, set Nslices and number of diffusion directions to 1 to accelerate script

% system
sys_type                  = 'Connectome2'; % prisma, skyra, Connectome2, C2_simulate_prisma, trio, prisma_XA30A
slew_safety_magrin        = 0.95;
grad_safety_magrin        = 0.9;
spiral_slew_safety_margin = 0.7; % to reduce PNS. decrease this would not lengthen readout too much
spiral_grad_safety_margin = 0.5;
lowPNS_slew_safety_margin = 0.3;   % used in rewinder to go back to k-space center.
lowPNS_grad_safety_margin = grad_safety_magrin;
diff_slew_safety_margin   = 0.45; % decrease this to reduce PNS, this would not lengthen TE too much
diff_grad_safety_margin   = 0.97;
spiSystemMargin           = 1e-6; % this is needed to avoid slewRate violation due to rouding errors

% geometry
fov            = 224e-3;
Nx             = 408;%320;
Nslices        = 30;
sliceThickness = 2e-3;

% timing related
Nechoes                  = 1;
Ndummies                 = 1;
TR                       = 6000e-3;
TE                       = 19e-3;
rfRingdownTime           = 100e-6;
rfDeadtime               = 100e-6;
fsSpoilerPhaseDispersion = 3*pi; % make it small, to avoid PNS. Handbook, eq. 10.87, here I use delta_phi much smaller than that used in example (5pi)
newfsSpoilerDuration     = 2500e-6; % not use the shortest possible duration, to avoid PNS
oldfsSpoilerDuration     = 5000e-6; % not use the shortest possible duration, to avoid PNS
trig_dur                 = 20e-6;
delay_after_trig         = 200e-6; % Gradient-free time for probe excitation [Unit: s] https://github.com/SkopeMagneticResonanceTechnologies/Pulseq-Sequences/blob/2576a44edb2f4123774565862472fe76ecd10435/sequences/PulseqBase.m#L102C9-L102C60
isDelayBetweenShots      = 0;
crusher_area             = 1e3; % unit is 1/m


% diff
bvalue0_threshold    = 0; % when b-value smaller than this, crusher would be used instead of diffusion gradient
bvalues              = [bvalue0_threshold 1e3        ];
% Define acquisition order as a cell array
% Each row: {b-value index, bTable row indices or 0 for b=0}
% acquisition_order = {
%     1, 0;           % 1st: b=1, 1 direction (crusher, no bTable)
%     2, 1:5;        % 2nd: b=1000, 5 directions (bTable rows 1-5)
%     1, 0;           % 3rd: b=1, 1 direction (crusher, no bTable)
%     2, 6:10        % 4th: b=1000, 5 directions (bTable rows 6-10)
% };
acquisition_order = {
    1, zeros(1,15);           % 1st: b=1, 1 direction (crusher, no bTable)
};
isIterativeFindGdiff = true;

% spiral
Narms             = 9; % number of spiral interleaves when designing spiral
Nshots            = 9; % controls acceleration factor
N_out_in          = 1;
isSpiralIn        = 0;
isEsnails         = 0;
isUseBrianVDScode = 1; % Sometimes Miki Lustig's spiral design may violate slew rate constraint, hence Brian's code is preferred
isVds             = 0;
if isEsnails
    isUseBrianVDScode = 0;
    isVds             = 1;
end
low_res_img_resolution  = 8e-3; % esnails option, center calibration region
acc_factor_kcenter      = 2; % esnails option
ksp_edge_start_location = 0.25; % esnails option, from here to the kmax we use the conventional sampling density
isFixROtoHalfTE         = 0; % fix readout duration to half TE to satisfy CPMG
golden_angle            = pi*222.4969/180; % https://www.sciencedirect.com/science/article/abs/pii/S1090780715002037, https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.27085
rot_angle               = 2*pi/Nshots;
angle_beteen_outin_pair = golden_angle;%2*pi/N_out_in;%golden_angle;
adcDwell                = 2e-6;
adcDwellSpiralDesign    = 100e-9;2e-6; % only useful if using Brian's spiral design code. It is suggested
isADCduringRewinder     = 0;


% rf and crusher
tRFex     = 3e-3;
tRFref        = 3e-3;
fspS          = 0.5;
fspS2         = 0.05; % should be different from 1st one if not CPMG
fcrusher180XY = 0.75;
fcrusher180XY2 = 0.25;

% options preferred to set to true
isStrictTimingCheck     = 0;
isBridgeZgrad           = 0; % if true, bridge z gradients that are around rf180
isSliceInterleave       = 1;

% options preferred to set to false, i.e. these are not fully tested
isJeffDuynStyleMeasTraj = 0; % Simple Correction Method fork-Space Trajectory Deviations in MRI, https://www.sciencedirect.com/science/article/pii/S1090780798913969?via%3Dihub
isNewFatSat         = 1;
isXYCrusher180      = 0; % x,y crusher around 180 refocus pulse
isOrthogonalCrusher = 0;
isSSgradReversal    = 0; % slice selection gradient reversal, SSGR. may lead to larger pns
if isSSgradReversal
    tRFref     = 2*tRFref;
end


% debug options, rarely used and might be buggy. Should set to false unless
% you know what you are doing
isNoGradInFirstEcho = 0; % set XY grad to 0 during 1st echo
isDifferentCrusher = 0;


if isTestRun
   Nslices = 1;
   msg = 'Test run, set Nslices to 1';
   fprintf('<strong>%s</strong>\n',msg);
    warndlg(msg)
end


%% constants
% MAX_ADCS  = 128; % on ve11c, the upper limit for this is 100, on XA30 128 works
MAX_ADCS    = 100; % on ve11c, the upper limit for this is 100, on XA30 128 works
rf90_phase  = 0; % MZ: we need to maintain these as variables because we will overwrtite phase offsets for multiple slice positions
rf180_phase = pi/2;
if isStrictTimingCheck
    timingCheckTol = 1e-10;
else
    timingCheckTol = 1e10;
end


%% parse some inputs
assert( (spiral_slew_safety_margin+spiSystemMargin)<=1,'Slew safety margin too large');
assert( (spiral_grad_safety_margin+spiSystemMargin)<=1,'Grad safety margin too large');

if isSpiralIn
    N_spi_one_pair = 2;
else
    N_spi_one_pair = 1;
end
% bval and bvec
bTable      = xlsread('bvec_18dir.xlsx');
bTable(:,3) = -bTable(:,3);
% bTable = [0 0 1;0 0 1];warning('Debug, use 0 0 1 diff direction')


% Initialize output arrays
num_bvolumes = sum(cellfun(@(x) length(x), acquisition_order(:, 2))); % Total volumes
bval = zeros(num_bvolumes, 1);
bvec = zeros(num_bvolumes, 3);
current_volume = 1;
% Track volume ranges for output
volume_ranges = cell(size(acquisition_order, 1), 1); % Store [start, end] for each segment
% Generate bval and bvec based on acquisition order
for iter = 1:size(acquisition_order, 1)
    bval_idx = acquisition_order{iter, 1}; % Index into bvalues
    btable_rows = acquisition_order{iter, 2}; % bTable rows or 0 (for b=0)
    num_dirs = length(btable_rows); % Number of directions from bTable rows
    
    % Assign b-value
    bval(current_volume:current_volume+num_dirs-1) = bvalues(bval_idx);
    
    % Assign directions
    if bvalues(bval_idx) <= bvalue0_threshold
        % For b=0 (b=1), use crusher direction [0 0 1]
        v = repmat([0 0 1], num_dirs, 1);
    else
        % For b>0, use specified bTable rows
        v = bTable(btable_rows, :);
    end
    bvec(current_volume:current_volume+num_dirs-1, :) = v;
    
    % Store volume range for this segment
    volume_ranges{iter} = [current_volume, current_volume+num_dirs-1];
    
    current_volume = current_volume + num_dirs;
end
% Validation
assert(length(bval) == num_bvolumes, 'bval size mismatch');
assert(size(bvec, 1) == num_bvolumes, 'bvec size mismatch');
% Optimized output: Print b-value and volume range in acquisition order
fprintf('Generated %d b-volumes in the following order:\n', num_bvolumes);
for iter = 1:size(acquisition_order, 1)
    bval_idx = acquisition_order{iter, 1};
    range = volume_ranges{iter};
    if range(1) == range(2)
        fprintf('  Volume %d: b=%d\n', range(1), bvalues(bval_idx));
    else
        fprintf('  Volumes %d-%d (%d volumes): b=%d\n', range(1), range(2), range(2)-range(1)+1, bvalues(bval_idx));
    end
end
% Summary of total volumes per b-value
fprintf('Summary:\n');
for b = 1:length(bvalues)
    count = sum(abs(bval - bvalues(b)) < eps); % Count volumes for each b-value
    fprintf('  b=%d: %d volumes\n', bvalues(b), count);
end


nFrames = numel(bval) + Ndummies;

% Insert bval=0 and bvec=[0 0 0] for the dummy and refscan, i.e. dummy and
% ref scan are treated as b=0. 
% Note that b value of 0 is not strictly 0, it is a very small b value, 
% i.e. bvalue0_threshold. And bvec is not [0 0 0] either, it is [0 0 1]
nImagingFrame = numel(bval);
bval_all_frames = [bvalue0_threshold*ones(nFrames-nImagingFrame,1);bval];
% bvec = [zeros(nFrames-nImagingFrame,3);bvec];
bvec = [repmat([0 0 1],nFrames-nImagingFrame,1);bvec];
maxbval = max(bval_all_frames);
bFactor_scale = sqrt(bval_all_frames./maxbval);
assert(~ismember([0 0 0], bvec, 'rows'), 'In the current implementation, b value=0 (i.e. bvec=0 0 0) is treated as a small b value, thus its direction cannot be [0 0 0], it should be set to [0 0 1], where the last 1 means z direction diff gradient, which is used as crusher. Please change your btable accordingly.');


% system defaults
if strcmp(sys_type,'prisma') || strcmp(sys_type,'C2_simulate_prisma') || strcmp(sys_type,'prisma_XA30A')
    physical_slew_max = 200;
    physical_grad_max = 80;
    B0=2.89; % 1.5 2.89 3.0
elseif strcmp(sys_type,'Connectome2')
    physical_slew_max = 598.802;
    physical_grad_max = 500;
    B0=2.89;
elseif strcmp(sys_type,'skyra')
    physical_slew_max = 180;
    physical_grad_max = 43;
    B0=2.89;
elseif strcmp(sys_type,'trio')
    physical_slew_max = 170;
    physical_grad_max = 38;
    B0=2.89;
else
    error('Undefined')
end
rfDeadTime =  100e-6;
adcDeadTime = 20e-6;
%     adcRasterTime = 2e-6;
adcRasterTime = 100e-9;
rfRasterTime = 1e-6;
gradRasterTime = 10e-6;
blockDurationRaster = 10e-6;
sys = mr.opts('MaxGrad',physical_grad_max*grad_safety_magrin,'GradUnit','mT/m',...
    'MaxSlew',physical_slew_max*slew_safety_magrin,'SlewUnit','T/m/s',...
    'rfDeadTime', rfDeadTime, ...
    'rfRingdownTime', rfRingdownTime, ...
    'adcDeadTime', adcDeadTime,...
    'adcRasterTime', adcRasterTime,...
    'rfRasterTime', rfRasterTime,...
    'gradRasterTime', gradRasterTime,...
    'blockDurationRaster', blockDurationRaster,...
    'B0',B0);
sys_lowPNS = mr.opts('MaxGrad',physical_grad_max*lowPNS_grad_safety_margin,'GradUnit','mT/m',...
    'MaxSlew',physical_slew_max*lowPNS_slew_safety_margin,'SlewUnit','T/m/s',...
    'rfDeadtime', rfDeadTime, ...
    'rfRingdownTime', rfRingdownTime, ...
    'adcDeadTime', adcDeadTime,...
    'adcRasterTime', adcRasterTime,...
    'rfRasterTime', rfRasterTime,...
    'gradRasterTime', gradRasterTime,...
    'blockDurationRaster', blockDurationRaster,...
    'B0',B0);
sys_diff = mr.opts('MaxGrad',physical_grad_max*diff_grad_safety_margin,'GradUnit','mT/m',...
    'MaxSlew',physical_slew_max*diff_slew_safety_margin,'SlewUnit','T/m/s',...
    'rfDeadtime', rfDeadTime, ...
    'rfRingdownTime', rfRingdownTime, ...
    'adcDeadTime', adcDeadTime,...
    'adcRasterTime', adcRasterTime,...
    'rfRasterTime', rfRasterTime,...
    'gradRasterTime', gradRasterTime,...
    'blockDurationRaster', blockDurationRaster,...
    'B0',B0);
sys_spi = mr.opts('MaxGrad',physical_grad_max*(spiral_grad_safety_margin+spiSystemMargin),'GradUnit','mT/m',...
    'MaxSlew',physical_slew_max*(spiral_slew_safety_margin+spiSystemMargin),'SlewUnit','T/m/s',...
    'rfDeadtime', rfDeadTime, ...
    'rfRingdownTime', rfRingdownTime, ...
    'adcDeadTime', adcDeadTime,...
    'adcRasterTime', adcRasterTime,...
    'rfRasterTime', rfRasterTime,...
    'gradRasterTime', gradRasterTime,...
    'blockDurationRaster', blockDurationRaster,...
    'B0',B0);

seq    = mr.Sequence(sys);      % Create a new sequence object


%% misc
% sanity check
if max(bvalues)>0
    minTE = calc_minTE_spiral(max(bvalues), diff_grad_safety_margin*physical_grad_max, diff_slew_safety_margin*physical_slew_max, sys.gradRasterTime);
    assert(minTE<=TE,'TE too small, minial TE is %.1f ms, current TE is %.1f ms', minTE*1e3, TE*1e3);
    fprintf('<strong>For single echo, minimal TE is %.1f ms, current TE is %.1f ms\n</strong>', minTE*1e3, TE*1e3);
end

% change slice index, for interleaving slices: 1,3,5; 2,4,6
if isSliceInterleave
%     slice_indices_odd = 1:2:Nslices;
%     slice_indices_even = 2:2:Nslices;
%     slice_indices = [slice_indices_odd,slice_indices_even];
    slice_indices = tdr_sliceorder(Nslices,1);
else
    slice_indices = 1:Nslices;
end

% trigger for skope
trig = mr.makeDigitalOutputPulse('ext1','duration', trig_dur); % possible channels: 'osc0','osc1', ext1

%
if isDifferentCrusher
    assert(~isequal(fspS,fspS2), 'Choose to use different cursher, but their area is same.')
end


%% spiral design
%%%%%%%%%%%%%%%%%%%% VDS design %%%%%%%%%%%%%%%%%%%%
if isUseBrianVDScode
    Ninterleaves = Narms;		 % Interleaves
    if isVds
        Fcoeff = [fov*100 -fov*100/2]; 	% FOV decreases linearly from 24 to 12cm.
    else
        Fcoeff = [fov*100];
        
        GmaxFOV = 1/4257 /Fcoeff/adcDwellSpiralDesign;
        GmaxFOV = GmaxFOV * 10;
        if mr.convert(sys_spi.maxGrad,'Hz/m','mT/m') > GmaxFOV
            warning('System Gmax is %.1f mT/m, while Gmax available for spiral design is only %.1f mT/m, this is due to improper dwell time setting when designing spiral, please consider decreasing it\n',mr.convert(sys_spi.maxGrad,'Hz/m','mT/m'),GmaxFOV);
        end
    end
    res_mm = 1e3*fov/Nx; % millimetter
    rmax = 5/res_mm;		% cm^(-1), corresponds to 1mm resolution. % rmax = 1/(2*resolution); (cm^-1)


    if isFixROtoHalfTE
        if Nechoes<=1
            warning('Fixing readout duration to 0.5*TE is only useful for multi echo acquisition. Now this option is ignored.')
        else
            rewinderAndCrusherDur = 4e-3;
            Fcoeff = func_iterative_change_FOV(TE/2-rewinderAndCrusherDur,spiral_slew_safety_margin*physical_slew_max, spiral_grad_safety_margin*physical_grad_max,adcDwellSpiralDesign,Ninterleaves,fov,res_mm,sys.gradRasterTime);
            tmpFcoeff = Fcoeff/(fov*100);
            fprintf('<strong>Normalized FOV at k-space edge is %.2f</strong>\n',tmpFcoeff(1)+tmpFcoeff(2));
        end
    end

    disp('Calculating Gradient');
    [k,g,s,time,r,theta] = vds(spiral_slew_safety_margin*physical_slew_max*100,spiral_grad_safety_margin*physical_grad_max/10,adcDwellSpiralDesign,Ninterleaves,Fcoeff,rmax);

    underSampling_factor = round( sys.gradRasterTime/adcDwellSpiralDesign );
    k = k(1:underSampling_factor:end);
    g = g(1:underSampling_factor:end);
    s = s(1:underSampling_factor:end);

    disp('Plotting Gradient');
    g = [real(g(:)) imag(g(:))];
    % plotgradinfo(g,T);

    g_mT_m = 10*g; % G/cm to mT/m
    s_T_m_s = 1e-2*[real(s(:)),imag(s(:))]; % G/cm/s to T/m/s
    k_std = k*100; % cm-1 to m-1
    k_std = [real(k_std(:)) imag(k_std(:))];
    k_std = [k_std, zeros(size(k_std,1),1)]; % kz=0
    grad_std = mr.convert(g_mT_m,'mT/m');
    slew_std = mr.convert(s_T_m_s,'T/m/s');

else
    res_mm = 1e3*fov/Nx; % millimetter
    if isVds % this is for esnails
        fprintf('<strong>esnails spiral design</strong>\n');
        low_res_img_kmax = 1/low_res_img_resolution/2;
        full_img_kax = 1/(res_mm*1e-3)/2;
        k_ratio_full_sampling = low_res_img_kmax./full_img_kax;
        k_ratio2 = max([k_ratio_full_sampling + 0.01, ksp_edge_start_location]);
        % [k_riv,g_riv,s_riv,time_riv,Ck_riv] = vdSpiralDesign(Narms,0, res_mm,100*[fov,fov,fov/2,fov/2],[0,k_ratio_full_sampling,k_ratio2,1],physical_grad_max*spiral_grad_safety_margin/10,spiral_slew_safety_margin*physical_slew_max/10,1e3*sys.gradRasterTime,[],'pchip');
        sampling_locs = [0,k_ratio_full_sampling,k_ratio2,1];
        sampling_dens = [Narms*fov/acc_factor_kcenter,Narms*fov/acc_factor_kcenter,fov,fov];
        [k_riv,g_riv,s_riv,time_riv,Ck_riv] = vdSpiralDesign(Narms,0, res_mm,100*sampling_dens,sampling_locs,physical_grad_max*spiral_grad_safety_margin/10,spiral_slew_safety_margin*physical_slew_max/10,1e3*sys.gradRasterTime,[],'pchip');
        fprintf('<strong>K-space sampling locations: %s</strong>\n', num2str(sampling_locs,'%.2f '));
        fprintf('<strong>K-space sampling densities: %s</strong>\n', num2str(sampling_dens/fov,'%.2f '));
        seq.setDefinition('k_ratio_full_sampling',k_ratio_full_sampling);

        [~,~,~,time_riv_cds,~] = vdSpiralDesign(Narms,0, res_mm,100*[fov,fov],[0,1],physical_grad_max*spiral_grad_safety_margin/10,spiral_slew_safety_margin*physical_slew_max/10,1e3*sys.gradRasterTime,[],'pchip');
        fprintf('<strong>Esnails time is %.1f ms, cds time is %.1f ms</strong>\n', time_riv,time_riv_cds);
    else
        [k_riv,g_riv,s_riv,time_riv,Ck_riv] = vdSpiralDesign(Narms,0, res_mm,100*[fov,fov],[0,1],physical_grad_max*spiral_grad_safety_margin/10,spiral_slew_safety_margin*physical_slew_max/10,1e3*sys.gradRasterTime,[],'pchip');
    end
    g_mT_m = g_riv(:,1:2)*10; % G/cm to mT/m
    grad_std = mr.convert(g_mT_m,'mT/m');
    s_T_m_s = s_riv(:,1:2)*10;
    slew_std = mr.convert(s_T_m_s,'T/m/s');
    k_std = k_riv*100;

    L = length(s_riv);
    time_array = linspace(0,time_riv,size(g_riv,1));
    figure
    subplot(2,2,1)
    plot(k_riv(:,1), k_riv(:,2));
    title('k-space');% axis([-6 6 -6 6]);
    subplot(2,2,2)
    plot(time_array,g_riv(:,1));
    %axis([0,L,-4.5,4.5]);
    title('gradient waveforms')
    hold on
    plot(time_array,g_riv(:,2), 'r');
    xlabel 'ms'
    legend('gx', 'gy', 'Location', 'NorthEast');
    subplot(2,2,3)
    plot(time_array,(g_riv(:,1).^2 + g_riv(:,2).^2).^0.5, 'r');
    xlabel 'ms'
    %axis([0 L 0 6]);
    title('gradient magnitude')
    subplot(2,2,4), plot(time_array(1:end-1),(s_riv(:,1).^2 + s_riv(:,2).^2).^0.5, 'r');
    xlabel 'ms'
    title('slew-rate magnitude');

end
%%%%%%%%%%%%%%%%%%%% VDS design %%%%%%%%%%%%%%%%%%%%

% gradients and slew rates
% ka = k_std;
% ka = [real(ka); imag(ka)];
ga = grad_std.';
sa = slew_std.';

% [gos, sos]=mr.traj2grad(kopt_smooth);
gos = ga;
sos = sa;

g_for_plot = g_mT_m.';
s_for_plot = s_T_m_s.';
figure;
subplot(211);plot((0:size(g_for_plot,2)-1)*sys.gradRasterTime*1e3,[g_for_plot;abs(g_for_plot(1,:)+1i*g_for_plot(2,:))]');title('gradient with abs(grad) constraint');ylabel 'mT/m';xlabel 'Time (ms)'
subplot(212);plot((0:size(s_for_plot,2)-1)*sys.gradRasterTime*1e3,[s_for_plot;abs(s_for_plot(1,:)+1i*s_for_plot(2,:))]');title('slew rate with abs(slew) constraint');ylabel 'mT/m/ms';xlabel 'Time (ms)'
% subplot(212);plot((0:size(g_for_plot,2)-2)*sys.gradRasterTime*1e3,[s_for_plot;abs(s_for_plot(1,:)+1i*s_for_plot(2,:))]');title('slew rate with abs(slew) constraint');ylabel 'mT/m/ms';xlabel 'Time (ms)';warning('tmp')

spiral_grad_shape=gos;

% extend spiral_grad_shape by repeating the last sample
% this is needed to accomodate for the ADC tuning delay
spiral_grad_shape = [spiral_grad_shape spiral_grad_shape(:,end)];

% readout grad
if isSpiralIn  
    [kx_spi_out,ky_spi_out] = mr.addRamps({k_std(:,1).',k_std(:,2).'},sys_lowPNS);
    gx_spi_out = mr.makeArbitraryGrad('x',mr.traj2grad(kx_spi_out),sys_spi);
    gy_spi_out = mr.makeArbitraryGrad('y',mr.traj2grad(ky_spi_out),sys_spi);
    gx_spi_in = mr.makeArbitraryGrad('x',flip(gx_spi_out.waveform),sys_spi);
    gy_spi_in = mr.makeArbitraryGrad('y',flip(gy_spi_out.waveform),sys_spi);
    
    
    spiral_grad_shape = [gx_spi_out.waveform;gy_spi_out.waveform];    
else
    gx_spi_out = mr.makeArbitraryGrad('x',spiral_grad_shape(1,:),'system',sys_spi,'Delay',sys.adcDeadTime);
    gy_spi_out = mr.makeArbitraryGrad('y',spiral_grad_shape(2,:),'system',sys_spi,'Delay',sys.adcDeadTime);
    
%     gx_spi_out = mr.makeArbitraryGrad('x',mr.traj2grad(k_std(:,1).'),'system',sys_spi,'Delay',sys.adcDeadTime);
%     gy_spi_out = mr.makeArbitraryGrad('y',mr.traj2grad(k_std(:,2).'),'system',sys_spi,'Delay',sys.adcDeadTime);
 
    gx_spi_out.first=0;
    gy_spi_out.first=0;
end
base_delay_between_shots = ceil(  mr.calcDuration(gx_spi_out,gy_spi_out)/Nshots            /sys.blockDurationRaster) * sys.blockDurationRaster;
base_delay_between_shots = 1e-3;

fprintf('<strong>Spiral out gradient waveform duration is %.1f ms</strong>\n',mr.calcDuration(gx_spi_out,gy_spi_out)*1e3);
% gx_spi_in = mr.makeArbitraryGrad('x',-spiral_grad_shape(1,end:-1:1),sys_spi);
% gy_spi_in = mr.makeArbitraryGrad('y',-spiral_grad_shape(2,end:-1:1),sys_spi);
% gx_spi_in = mr.makeArbitraryGrad('x',-spiral_grad_shape(1,end:-1:1),sys);
% gy_spi_in = mr.makeArbitraryGrad('y',-spiral_grad_shape(2,end:-1:1),sys);


% gx_spi_in = mr.makeArbitraryGrad('x',spiral_grad_shape(1,end:-1:1),sys);
% gy_spi_in = mr.makeArbitraryGrad('y',spiral_grad_shape(2,end:-1:1),sys);

% get rewinder, move from (kx,ky) to (-kx,ky)
% mr.addRamps
% k0 = k_std(end-1:end,:);
% k0(:,3) = 0;
% k0 = k0.';
% kend = -k0;
% kend(:,[1 2]) = kend(:,[2 1]);
% [kout, success] = mr.calcRamp(k0,kend,sys);
% 
% gx_rewinder = mr.makeArbitraryGrad('x',mr.traj2grad(kout(1,:)),sys);
% gy_rewinder = mr.makeArbitraryGrad('y',mr.traj2grad(kout(2,:)),sys);
            % % sometimes we cannot do it in a single step, have to go back to
            % % center first and then to -kx,-ky
            % try
            %     [gx_spi_in_prephasing, times, amplitudes] = mr.makeExtendedTrapezoidArea('x',gx_spi_out.last,gx_spi_in.first, -2*gx_spi_out.area, sys_lowPNS);
            %     [gy_spi_in_prephasing, times, amplitudes] = mr.makeExtendedTrapezoidArea('y',gy_spi_out.last,gy_spi_in.first, -2*gy_spi_out.area, sys_lowPNS);
            % 
            %     if mr.calcDuration(gx_spi_in_prephasing) > mr.calcDuration(gy_spi_in_prephasing)
            %         tmp = mr.align('left',gx_spi_in_prephasing,'right',gy_spi_in_prephasing);
            %         gx_spi_in_prephasing = tmp{1};
            %         gy_spi_in_prephasing = tmp{2};
            %     else
            %         tmp = mr.align('left',gy_spi_in_prephasing,'right',gx_spi_in_prephasing);
            %         gx_spi_in_prephasing = tmp{2};
            %         gy_spi_in_prephasing = tmp{1};
            %     end
            % catch
            %     gx_spi_in_prephasing = [];
            % end

            if 0%isempty(gx_spi_in_prephasing) % go back to 0,0 first then to -kx,-ky
                [gx_spi_in_prephasing1, times, amplitudes] = mr.makeExtendedTrapezoidArea('x',gx_spi_out.last,0, -gx_spi_out.area, sys_lowPNS);
                [gy_spi_in_prephasing1, times, amplitudes] = mr.makeExtendedTrapezoidArea('y',gy_spi_out.last,0, -gy_spi_out.area, sys_lowPNS);

                [gx_spi_in_prephasing2, times, amplitudes] = mr.makeExtendedTrapezoidArea('x',0,gx_spi_in.first, -gx_spi_out.area, sys_lowPNS);
                [gy_spi_in_prephasing2, times, amplitudes] = mr.makeExtendedTrapezoidArea('y',0,gy_spi_in.first, -gy_spi_out.area, sys_lowPNS);

                if mr.calcDuration(gx_spi_in_prephasing2) > mr.calcDuration(gy_spi_in_prephasing2)
                    tmp = mr.align('left',gx_spi_in_prephasing2,'right',gy_spi_in_prephasing2);
                    gx_spi_in_prephasing2 = tmp{1};
                    gy_spi_in_prephasing2 = tmp{2};
                else
                    tmp = mr.align('left',gy_spi_in_prephasing2,'right',gx_spi_in_prephasing2);
                    gx_spi_in_prephasing2 = tmp{2};
                    gy_spi_in_prephasing2 = tmp{1};
                end

                gx_spi_in_prephasing = mr.addGradients({gx_spi_in_prephasing1,gx_spi_in_prephasing2},'system',sys_lowPNS);
                gy_spi_in_prephasing = mr.addGradients({gy_spi_in_prephasing1,gy_spi_in_prephasing2},'system',sys_lowPNS);
            end



% TODO, the optimal rewinder duration may be different between out-in and
% in-out
gx_rewinder_dur = ceil(abs(spiral_grad_shape(1,end))/sys.maxSlew/seq.gradRasterTime)*seq.gradRasterTime;
gy_rewinder_dur = ceil(abs(spiral_grad_shape(2,end))/sys.maxSlew/seq.gradRasterTime)*seq.gradRasterTime;
gx_gy_rewinder_max_dur = max(gx_rewinder_dur,gy_rewinder_dur);
gx_rewinder_out = mr.makeExtendedTrapezoid('x','times',[0,gx_gy_rewinder_max_dur],'amplitudes',[spiral_grad_shape(1,end),0]);
gy_rewinder_out = mr.makeExtendedTrapezoid('y','times',[0,gx_gy_rewinder_max_dur],'amplitudes',[spiral_grad_shape(2,end),0]);
gx_rewinder_in  = mr.makeExtendedTrapezoid('x','times',[0,gx_gy_rewinder_max_dur],'amplitudes',[0,spiral_grad_shape(1,end)]);
gy_rewinder_in  = mr.makeExtendedTrapezoid('y','times',[0,gx_gy_rewinder_max_dur],'amplitudes',[0,spiral_grad_shape(2,end)]);
if Nechoes > 1
    try
        gx_to_ksp_center = mr.makeExtendedTrapezoidArea('x',gx_spi_out.waveform(end),0, -gx_spi_out.area, sys_lowPNS);
        gy_to_ksp_center = mr.makeExtendedTrapezoidArea('y',gy_spi_out.waveform(end),0, -gy_spi_out.area, sys_lowPNS);
    catch
        tic
        [gx_to_ksp_center,gy_to_ksp_center,tmps,tmpg,tmpt] = grid_search_rewinder(physical_grad_max*lowPNS_grad_safety_margin, physical_slew_max*lowPNS_slew_safety_margin, gx_spi_out.waveform(end),gy_spi_out.waveform(end), -gx_spi_out.area,-gy_spi_out.area,sys);
        elapsedTime = toc;
        fprintf('<strong>Failed to find the rewinder, now do grid search</strong>\n')
        fprintf('<strong>Uses %.1f seconds to grid search rewinder</strong>\n', elapsedTime);
        fprintf('<strong>Rewinder duration=%.1f ms, grad=%.1f mT/m, s=%.1f T/m/s</strong>\n',tmpt*1e3,tmpg,tmps);
    end

    gx_rewinder_out = gx_to_ksp_center;
    gy_rewinder_out = gy_to_ksp_center;
end
deltak=1/fov;
gz_spoil=mr.makeTrapezoid('z',sys,'Area',deltak*Nx*4);
% gz_spoil=mr.makeTrapezoid('z',sys,'Area',deltak*112*4);warning('For debug, set Nx in gz spoil to 112')
gx_spoil=mr.makeExtendedTrapezoid('x','times',[0 mr.calcDuration(gz_spoil)],'amplitudes',[spiral_grad_shape(1,end),0]); %todo: make a really good spoiler
gy_spoil=mr.makeExtendedTrapezoid('y','times',[0 mr.calcDuration(gz_spoil)],'amplitudes',[spiral_grad_shape(2,end),0]); %todo: make a really good spoiler
% [gx_prephasing, times, amplitudes] = mrcustomized.makeExtendedTrapezoidArea_xingwang('x', 0, gx_in.first, gx_in.area, sys);
% [gy_prephasing, times, amplitudes] = mrcustomized.makeExtendedTrapezoidArea_xingwang('y', 0, gy_in.first, gy_in.area, sys);


%% fat sat and rf90 rf180
warning('OFF', 'mr:restoreShape'); % restore shape is not compatible with spirals and will throw a warning from each plot() or calcKspace() call

% Create fat-sat pulse
% (in Siemens interpreter from January 2019 duration is limited to 8.192 ms, and although product EPI uses 10.24 ms, 8 ms seems to be sufficient)
sat_ppm=-3.45;
sat_freq=sat_ppm*1e-6*B0*sys.gamma;
rf_fs = mr.makeGaussPulse(110*pi/180,'system',sys,'Duration',8e-3,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
if isNewFatSat
    gz_fs = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',fsSpoilerPhaseDispersion/sliceThickness);
else
    gz_fs = mr.makeTrapezoid('z',sys_lowPNS,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm
%     if (gz_fs.riseTime+gz_fs.fallTime+gz_fs.flatTime) < oldfsSpoilerDuration
%         gz_fs = mr.makeTrapezoid('z',sys_lowPNS,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4,'duration',oldfsSpoilerDuration);
%     end
end
gx_fs_pre = mr.makeTrapezoid('x',sys,'Area',fsSpoilerPhaseDispersion/sliceThickness,'duration',newfsSpoilerDuration);
gy_fs_pre = mr.makeTrapezoid('y',sys,'Area',fsSpoilerPhaseDispersion/sliceThickness,'duration',newfsSpoilerDuration);
gz_fs_pre = mr.makeTrapezoid('z',sys,'Area',fsSpoilerPhaseDispersion/sliceThickness,'duration',newfsSpoilerDuration);
gx_fs_post = mr.makeTrapezoid('x',sys,'delay',mr.calcDuration(rf_fs),'Area',-fsSpoilerPhaseDispersion/sliceThickness,'duration',newfsSpoilerDuration);
gy_fs_post = mr.makeTrapezoid('y',sys,'delay',mr.calcDuration(rf_fs),'Area',-fsSpoilerPhaseDispersion/sliceThickness,'duration',newfsSpoilerDuration);
gz_fs_post = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',-fsSpoilerPhaseDispersion/sliceThickness,'duration',newfsSpoilerDuration);



spoiler_amp = 3*8*42.58*10e2;
est_rise = 500e-6;
est_flat = 2500e-6;

gp_r = mr.makeTrapezoid('x','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',sys_lowPNS);
gp_p = mr.makeTrapezoid('y','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',sys_lowPNS);
gp_s = mr.makeTrapezoid('z','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',sys_lowPNS);

gn_r = mr.makeTrapezoid('x','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',sys_lowPNS);
gn_p = mr.makeTrapezoid('y','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',sys_lowPNS);
gn_s = mr.makeTrapezoid('z','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',sys_lowPNS);


if 0
    % % % % Create 90  degree slice selection pulse and gradient
    % % % [rf90, gs_ex, gzReph] = mr.makeSincPulse(pi/2,'system',sys,'Duration',tRFex,...
    % % %     'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,...
    % % %     'PhaseOffset',rf90_phase);
    % % % % Create 180 degree slice refocusing pulse and gradients
    % % % [rf180, gz180] = mr.makeSincPulse(pi,'system',sys,'Duration',tRFref,...
    % % %     'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,...
    % % %     'PhaseOffset',rf180_phase,'use','refocusing');
    % % % if isSSgradReversal
    % % %     GSref = mr.makeTrapezoid('z',sys,'amplitude',-1*gs_ex.amplitude,'FlatTime',tRefwd);
    % % % else
    % % %     GSref = mr.makeTrapezoid('z',sys,'amplitude',gs_ex.amplitude,'FlatTime',tRefwd);
    % % % end
elseif 0
    % % % 
    % % % [rf90, gs_ex,gzReph] = mr.makeSincPulse(pi/2,sys,'Duration',tRFex,...
    % % %     'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',rf90_phase);
    % % % GSex = mr.makeTrapezoid('z',sys,'amplitude',gs_ex.amplitude,'FlatTime',tExwd);
    % % % 
    % % % [rf180, gz180] = mr.makeSincPulse(pi,sys,'Duration',tRFref,...
    % % %     'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',rf180_phase,'use','refocusing');
    % % % % rf180.delay = sys.rfDeadTime;
    % % % % tRefwd = rf180.delay + tRFref + sys.rfRingdownTime;
    % % % GSref = mr.makeTrapezoid('z',sys,'amplitude',GSex.amplitude,'FlatTime',tRefwd);
    % % % % if isSSgradReversal
    % % % %     GSref = mr.makeTrapezoid('z',sys,'amplitude',-1*GSex.amplitude,'FlatTime',tRefwd);
    % % % % else
    % % % %     GSref = mr.makeTrapezoid('z',sys,'amplitude',GSex.amplitude,'FlatTime',tRefwd);
    % % % % end
else
    dG = 250e-6;

    tEx=tRFex;
    tRef=tRFref;

    flipex=90*pi/180;
    [rf90, gs_ex,gzReph] = mr.makeSincPulse(flipex,sys_lowPNS,'Duration',tEx,...
        'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',rf90_phase,'use','excitation');
    % assert(isequal(rf90.delay,sys.rfRingdownTime),'rf delay should be equal to ringdowntime, otherwise rf is not symmetric w.r.t to gradient, which may cause artifacts');
    tExwd = rf90.delay + tEx + sys.rfRingdownTime;
    GSex = mr.makeTrapezoid('z',sys,'amplitude',gs_ex.amplitude,'FlatTime',tExwd,'riseTime',dG);

    flipref=180*pi/180;
    [rf180,gz180] = mr.makeSincPulse(flipref,sys_lowPNS,'Duration',tRef,...
        'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',rf180_phase,'use','refocusing');
    % assert(isequal(rf180.delay,sys.rfRingdownTime),'rf delay should be equal to ringdowntime, otherwise rf is not symmetric w.r.t to gradient, which may cause artifacts'); % dark band 2nd echo
    tRefwd = rf180.delay + tRFref + sys.rfRingdownTime; % rf180.delay already includes rfDeadTime
    GSref = mr.makeTrapezoid('z',sys,'amplitude',GSex.amplitude,'FlatTime',tRefwd,'riseTime',dG); % TODO, should we use a fixed rise time here or not? Perhaps does not matter as long as symmetric around rf180
    GSref = gz180;
end

% spi_simrf; % slice profile etc. directly copied from pulseq demos


%% crushers
% GS7 is left crusher of 180, GS5 is right crusher for 180.
% GS4 is slice selection gradient for 180 (only flat)
% These gradients are copied from writeTSE.m
% AGSex=gs_ex.area/2;
AGSex=GSex.area/2;
GSspr = mr.makeTrapezoid('z',sys_lowPNS,'area',AGSex*(1+fspS));
GSspex = mr.makeTrapezoid('z',sys_lowPNS,'area',AGSex*fspS);

if isBridgeZgrad
    GS4times=[0 GSref.flatTime];
    GS4amp=[GSref.amplitude GSref.amplitude];
    GS4 = mr.makeExtendedTrapezoid('z','system',sys,'times',GS4times,'amplitudes',GS4amp);

    GS5times=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
    GS5amp=[GSref.amplitude GSspr.amplitude GSspr.amplitude 0];
    [GS5times,ia] = unique(GS5times); % may not be unique if GSspr.flatTime==0
    GS5amp = GS5amp(ia);
    GS5 = mr.makeExtendedTrapezoid('z','system',sys,'times',GS5times,'amplitudes',GS5amp);

    GS7times=[0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
    GS7amp=[0 GSspr.amplitude GSspr.amplitude GSref.amplitude];
    [GS7times,ia] = unique(GS7times);
    GS7amp = GS7amp(ia);
    GS7 = mr.makeExtendedTrapezoid('z','system',sys_lowPNS,'times',GS7times,'amplitudes',GS7amp);


    % for later echoes
    GSspr_2 = mr.makeTrapezoid('z',sys,'area',AGSex*(1+fspS2),'Duration',mr.calcDuration(GS7));
    GS5times=[0 GSspr_2.riseTime GSspr_2.riseTime+GSspr_2.flatTime GSspr_2.riseTime+GSspr_2.flatTime+GSspr_2.fallTime];
    GS5amp=[GSref.amplitude GSspr_2.amplitude GSspr_2.amplitude 0];
    [GS5times,ia] = unique(GS5times);
    GS5amp = GS5amp(ia);
    GS5_2 = mr.makeExtendedTrapezoid('z','system',sys,'times',GS5times,'amplitudes',GS5amp);

    GS7times=[0 GSspr_2.riseTime GSspr_2.riseTime+GSspr_2.flatTime GSspr_2.riseTime+GSspr_2.flatTime+GSspr_2.fallTime];
    GS7amp=[0 GSspr_2.amplitude GSspr_2.amplitude GSref.amplitude];
    [GS7times,ia] = unique(GS7times);
    GS7amp = GS7amp(ia);
    GS7_2 = mr.makeExtendedTrapezoid('z','system',sys_lowPNS,'times',GS7times,'amplitudes',GS7amp);

else
    GS7 = GSspr;
    GS4 = GSref;
    GS5 = GSspr;

    GSspr_2 = mr.makeTrapezoid('z',sys,'area',AGSex*(1+fspS2),'Duration',mr.calcDuration(GS7));
    GS7_2 = GSspr_2;
    GS5_2 = GSspr_2;
end


if isXYCrusher180
    GXcrusher180 = mr.makeTrapezoid('x',sys,'area',AGSex*fcrusher180XY,'duration',mr.calcDuration(GS7));
    GYcrusher180 = mr.makeTrapezoid('y',sys,'area',AGSex*fcrusher180XY,'duration',mr.calcDuration(GS7));
    GXcrusher180_2 = mr.makeTrapezoid('x',sys,'area',AGSex*fcrusher180XY2,'duration',mr.calcDuration(GS7));
    GYcrusher180_2 = mr.makeTrapezoid('y',sys,'area',AGSex*fcrusher180XY2,'duration',mr.calcDuration(GS7));
else
    GXcrusher180 = mr.makeTrapezoid('x',sys,'area',0*fcrusher180XY,'duration',mr.calcDuration(GS7));
    GYcrusher180 = mr.makeTrapezoid('y',sys,'area',0*fcrusher180XY,'duration',mr.calcDuration(GS7));
    GXcrusher180_2 = mr.makeTrapezoid('x',sys,'area',0*fcrusher180XY2,'duration',mr.calcDuration(GS7));
    GYcrusher180_2 = mr.makeTrapezoid('y',sys,'area',0*fcrusher180XY2,'duration',mr.calcDuration(GS7));
end

if isOrthogonalCrusher && ~isXYCrusher180
    error('Must set isXYCrusher180=1 when trying to use orthogonal crusher')
end
if isOrthogonalCrusher
    GXcrusher180_2 = mr.scaleGrad(GXcrusher180_2,0);
    GYcrusher180_2 = mr.scaleGrad(GYcrusher180_2,-1);
end

gz180_crusher_1 = mr.makeTrapezoid('z',sys_lowPNS,'area',crusher_area);
gz180_crusher_2 = mr.makeTrapezoid('y',sys_lowPNS,'area',crusher_area);
gz180_crusher_3 = mr.makeTrapezoid('x',sys_lowPNS,'area',crusher_area);



%% ADC
if isADCduringRewinder
    adcTime = sys.gradRasterTime*numel(gx_spi_out.waveform);
else
    adcTime = sys.gradRasterTime*size(spiral_grad_shape,2);
end
adcSamplesDesired = round( adcTime/adcDwell );
adcSamplesPerSegment = ceil( adcSamplesDesired/MAX_ADCS );
adcSamplesPerSegment = ceil(adcSamplesPerSegment/100)*100; % must be integer multiple of 100ns
if adcSamplesDesired/adcSamplesPerSegment > MAX_ADCS
    %warning('Current #ADC in a block is %.1f, exceeding maximal #ADC events allowed, which is %d', adcSamplesDesired/adcSamplesPerSegment,MAX_ADCS);
    error('Current #ADC in a block is %.1f, exceeding maximal #ADC events allowed, which is %d', adcSamplesDesired/adcSamplesPerSegment,MAX_ADCS);
end
% adcSegments=round(adcSamplesDesired/adcSamplesPerSegment);
adcSegments=fix(adcSamplesDesired/adcSamplesPerSegment);
adcSamples=adcSegments*adcSamplesPerSegment;
adcDwell=round(adcTime/adcSamples/100e-9)*100e-9; % on Siemens adcDwell needs to be aligned to 100ns (if my memory serves me right)
if adcDwell < 1e-6
    error('ADC dwell < 1us, it is highly likely you would fail on Siemens scanner!')
end
adcSegmentDuration=adcSamplesPerSegment*adcDwell; % with the 100 samples above and the 100ns alignment we automatically fullfill the segment alignment requirement
% if mod(adcSegmentDuration, sys.gradRasterTime)>eps % this is too harsh
if abs(   adcSegmentDuration/sys.gradRasterTime - round(adcSegmentDuration/sys.gradRasterTime)    ) > 1e-10
    error('ADC segmentation model results in incorrect segment duration');
end
% update segment count
adcSegments=floor(adcTime/adcSegmentDuration);
adcSamples=adcSegments*adcSamplesPerSegment;
adc = mr.makeAdc(adcSamples,sys,'Dwell',adcDwell,'Delay',sys.adcDeadTime);
fprintf('<strong>ADC duration is %.1f ms\n</strong>', adc.duration*1e3);
adc2_delay = adcDwell * (adcSamplesDesired - adcSamples); % for spiral in
adc2 = mr.makeAdc(adcSamples,sys,'Dwell',adcDwell,'Delay',adc2_delay); % TODO, fix adc not aligned to raster

%% Calculate delay times
rf90centerInclDelay  = rf90.delay + mr.calcRfCenter(rf90);
rf180centerInclDelay = rf180.delay + mr.calcRfCenter(rf180);

delayTE1 = TE/2   -   (mr.calcDuration(gs_ex,rf90)-rf90centerInclDelay)   -  mr.calcDuration(gzReph) - rf180centerInclDelay;
delayTE1 = round(delayTE1/sys.gradRasterTime)*sys.gradRasterTime;
assert(  abs(delayTE1-(TE/2   -   (mr.calcDuration(gs_ex,rf90)-rf90centerInclDelay)   -  mr.calcDuration(gzReph) - rf180centerInclDelay)) < timingCheckTol, 'Timing wrong'  )

delayTE2 = TE/2   -  (mr.calcDuration(GSref)-rf180centerInclDelay) - mr.calcDuration(trig) - delay_after_trig;
delayTE2 = round(delayTE2/sys.gradRasterTime)*sys.gradRasterTime;
assert(  abs(delayTE2-(TE/2   -  (mr.calcDuration(GSref)-rf180centerInclDelay)-mr.calcDuration(trig)-delay_after_trig))  < timingCheckTol, 'Timing wrong' )

delayTE1_b0 = TE/2   -   (mr.calcDuration(gs_ex,rf90)-rf90centerInclDelay)   -  mr.calcDuration(gzReph) - rf180centerInclDelay - mr.calcDuration(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
delayTE1_b0 = round(delayTE1_b0/sys.gradRasterTime)*sys.gradRasterTime;

delayTE2_b0  = TE/2   -  (mr.calcDuration(GSref)-rf180centerInclDelay) - mr.calcDuration(trig) - delay_after_trig - mr.calcDuration(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
delayTE2_b0  = round(delayTE2_b0/sys.gradRasterTime)*sys.gradRasterTime;

if isBridgeZgrad
    error('not implemented yet')
    delayTE2_b0  = TE/2   -  mr.calcRfCenter(rf180) - sys.rfRingdownTime  - mr.calcDuration(GS5) - mr.calcDuration(trig) - delay_after_trig;
    delayTE2_b0  = round(delayTE2_b0/sys.gradRasterTime)*sys.gradRasterTime;
    assert(  abs(delayTE2_b0-(TE/2   -  mr.calcRfCenter(rf180) - sys.rfRingdownTime  - mr.calcDuration(GS5) - mr.calcDuration(trig) - delay_after_trig))  < timingCheckTol, 'Timing wrong' )
else
%     delayTE2_b0  = TE/2   -  (mr.calcDuration(GSref)-rf180centerInclDelay) - mr.calcDuration(trig) - delay_after_trig - mr.calcDuration(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
%     delayTE2_b0  = round(delayTE2_b0/sys.gradRasterTime)*sys.gradRasterTime;
end

delay_bewteen_RO_and_next180 = TE/2 - mr.calcDuration(gx_spi_out,gy_spi_out,adc) - mr.calcDuration(gx_rewinder_out,gy_rewinder_out) - mr.calcDuration(GS7)  -  rf180centerInclDelay;
delay_bewteen_RO_and_next180 = round(delay_bewteen_RO_and_next180/sys.gradRasterTime)*sys.gradRasterTime;

if isBridgeZgrad
    delayTE_laterEchoes = TE/2 - mr.calcRfCenter(rf180) - sys.rfRingdownTime  - mr.calcDuration(GS5) - mr.calcDuration(trig) - delay_after_trig;
else
    delayTE_laterEchoes = TE/2 - (mr.calcDuration(GSref)-rf180centerInclDelay)  - mr.calcDuration(GS5) - mr.calcDuration(trig) - delay_after_trig;
end
delayTE_laterEchoes = ceil(delayTE_laterEchoes/sys.gradRasterTime)*sys.gradRasterTime;

%%%%%%%%%%%%%%%%% for non-CPMG %%%%%%%%%%%%%%%%%
% see Haacke Green Bible Fig. 18.20 for a detailed timing
tau2 = TE/2 + mr.calcDuration(gx_spi_out,gy_spi_out) + mr.calcDuration(gx_rewinder_out,gy_rewinder_out) + mr.calcDuration(GS7_2) + 0.5*mr.calcDuration(GS4);
tau1 = TE/2;
delay_after_2nd_rf180 = tau2 - tau1; % 2nd secondary spin echo
delay_after_2nd_rf180 = tau1; % stimulated echo
% delay_after_2nd_rf180 = tau2; % SE4
% delay_after_2nd_rf180 = tau2 + tau1; % SE5
delay_after_2nd_rf180 = ceil(delay_after_2nd_rf180/sys.blockDurationRaster)*sys.blockDurationRaster;
%%%%%%%%%%%%%%%%% for non-CPMG %%%%%%%%%%%%%%%%%


                
                


if isNewFatSat
    fatSatModuleDur = mr.calcDuration(gx_fs_pre,gy_fs_pre,gz_fs_pre) + mr.calcDuration(gx_fs_post,gy_fs_post,gz_fs_post);
    fatSatModuleDur = mr.calcDuration(gp_r,gp_p,gp_s) + mr.calcDuration(rf_fs,gn_r,gn_p,gn_s);
else
    fatSatModuleDur = mr.calcDuration(gz_fs);
end
time_1_slice = fatSatModuleDur +  rf90centerInclDelay + 0.5*TE + (Nechoes-1)*TE + 0.5*TE + mr.calcDuration(gx_spi_out,gy_spi_out,adc)  +  mr.calcDuration(gx_spoil,gy_spoil,gz_spoil);
delayTR = ceil(  (TR- Nslices*time_1_slice)  /  seq.gradRasterTime  )*seq.gradRasterTime;

fprintf('<strong>Minimal TR is %.2f, current TR is %.2f\n</strong>', time_1_slice*Nslices, TR);
delay_after_1slice = TR/Nslices - time_1_slice;
assert(delay_after_1slice>=0, 'At least, you should increase TR to %.2f seconds',time_1_slice*Nslices);
delay_after_1slice = ceil(delay_after_1slice/sys.blockDurationRaster)*sys.blockDurationRaster;
fprintf('<strong>Delay after each shot is %.2f\n</strong>',delay_after_1slice)
assert(delayTE1>=0);
assert(delayTE2>=0);
assert(all(delayTR>=0));
if Nechoes>1
    assert(delay_bewteen_RO_and_next180>=0, 'You should increase TE to at least %.1f ms', 1e3*(TE+2*abs(delay_bewteen_RO_and_next180)));
end

% TODO, fix the expression for rampup, we may not use gmax here
small_delta = min([delayTE1,delayTE2])-ceil(sys_diff.maxGrad/sys_diff.maxSlew/sys_diff.gradRasterTime)*sys_diff.gradRasterTime;
big_delta   = delayTE1 + mr.calcDuration(GSref);
maxAvailBvalue = calc_bval_trap(physical_grad_max*diff_grad_safety_margin,small_delta,big_delta,ceil(sys_diff.maxGrad/sys_diff.maxSlew/sys_diff.gradRasterTime)*sys_diff.gradRasterTime);
assert(max(bvalues)<=maxAvailBvalue,'Required b-value too large, max achievable is %.1f, required is %.1f',maxAvailBvalue,max(bvalues));


%% define blocks
% loop structure
% 
% volume ( dummies, b0, b1k(direction1, direction2), b2k(direction1, direction2)... )
%   spiral arms ( 1st arm, 2nd arm(rotated version of 1st arm), 3rd arm...)
%       |<------------TR----------->|
%       slices (1,3,5,..., 2,4,6,...)
%           echoes (1st echo, 2nd echo, ...)
%               spiral-out and spiral-in pair (pair#1, pair#2). Currently this is only for single echo
% 
% scan time is #volumes * #spiral arms * TR
% 
% at the end, an optional trajectory measurement can be performed


adc_backup = adc;
dummy_adc = mr.makeDelay(mr.calcDuration(adc));
isAlreadyShowDelayAfterOneSlice = false;
isAlreadyShowBvalueDueToCrusher = false; 
isFirstTrigger = true;
for iterFrame = 1:nFrames
    isImagingFrame = true;
    if iterFrame<=Ndummies % dummy       
        isImagingFrame = false;
    end
    
    %     % TODO, fix the expression for rampup, we may not use gmax here
    %     if isIterativeFindGdiff
    %         gdiff_dur = min([delayTE1,delayTE2]);
    %         [gamp,rup] =  findDiffGradAmp(sys_diff.maxGrad,sys_diff.maxSlew,bval(iterFrame),gdiff_dur,big_delta,sys_diff.gradRasterTime);
    %         g = mr.convert(gamp,'mT/m');
    %     else
    %         g = sqrt(bval(iterFrame)*1e6/bFactCalc_trap(1,small_delta,big_delta,ceil(sys_diff.maxGrad/sys_diff.maxSlew/sys_diff.gradRasterTime)*sys_diff.gradRasterTime));
    %     end
    
    g=sqrt(maxbval*1e6/bFactCalc(1,small_delta,big_delta)); % for now it looks too large!
    gr=ceil(g/sys_diff.maxSlew/sys.gradRasterTime)*sys.gradRasterTime;
    gDiff=mr.makeTrapezoid('z','amplitude',g,'riseTime',gr,'flatTime',small_delta-gr,'system',sys_diff);
    assert(mr.calcDuration(gDiff)<=delayTE1,'TE too small');
    assert(mr.calcDuration(gDiff)<=delayTE2,'TE too small'); 
    
    g_x=g.*bvec(iterFrame,1);
    g_y=g.*bvec(iterFrame,2);
    g_z=g.*bvec(iterFrame,3);
    
    if bval_all_frames(iterFrame)<=bvalue0_threshold%((sum(bvec(iterFrame,:))==0)||(sum(bvec(iterFrame,:))==1)) % b=0 or dwi with diffusion gradient on one axis we keep using the older version
        gDiff_x=gDiff; gDiff_x.channel='x';
        gDiff_y=gDiff; gDiff_y.channel='y';
        gDiff_z=gDiff; gDiff_z.channel='z';
    else
        [azimuth,elevation,r] = cart2sph(g_x,g_y,g_z);
        polar= -(pi/2-elevation);
        
        Gr=mr.rotate('z',azimuth,mr.rotate('y',polar,gDiff));
        if size(Gr,2)==3
            gDiff_x=Gr{1,2};
            gDiff_y=Gr{1,3};
            gDiff_z=Gr{1,1};
        else
            if size(Gr,2)==2
                diffusion_blank=find( bvec(iterFrame,:)==0);
                switch diffusion_blank
                    case 2
                        gDiff_x=Gr{1,2};
                        gDiff_z=Gr{1,1};
                        gDiff_y=gDiff; gDiff_y.channel='y'; gDiff_y.amplitude=0; gDiff_y.area=0; gDiff_y.flatArea=0;
                    case 1
                        gDiff_z=Gr{1,1};
                        gDiff_y=Gr{1,2};
                        gDiff_x=gDiff; gDiff_x.amplitude=0; gDiff_x.area=0; gDiff_x.flatArea=0;gDiff_x.channel='x';
                    case 3
                        gDiff_x=Gr{1,2};
                        gDiff_y=Gr{1,1};
                        gDiff_z=gDiff; gDiff_z.amplitude=0; gDiff_z.area=0; gDiff_z.flatArea=0;gDiff_z.channel='z';
                end
            end
        end
    end
    
    
    if isImagingFrame
        number_of_shots = Nshots;
    else
        number_of_shots = 1;
    end
    
    for shot_ind = 1:number_of_shots
        % seq.addBlock(mr.makeLabel('SET','SLC', 0));
        % seq.addBlock(mr.makeLabel('INC','SEG', 1));
        phi = (shot_ind-1) * rot_angle;

        isDummy = ~isImagingFrame;
        if isDummy
            adc = dummy_adc;
        else
            adc = adc_backup;
        end                
        


        % gx_spi_out_rot = mrcustomized.rotate_xingwang('z',phi,sys,gx_spi_out);
        % gy_spi_out_rot = mrcustomized.rotate_xingwang('z',phi,sys,gy_spi_out);
        % gx_spi_out_rot = gx_spi_out_rot{1};
        % gy_spi_out_rot = gy_spi_out_rot{1};




        if 0
            % gx_spi_in = mr.makeArbitraryGrad('x',-spiral_grad_shape(1,end:-1:1),sys);
            % gy_spi_in = mr.makeArbitraryGrad('y',-spiral_grad_shape(2,end:-1:1),sys);
            % gx_spi_in = mrcustomized.rotate_xingwang('z',phi,sys,gx_spi_in);
            % gy_spi_in = mrcustomized.rotate_xingwang('z',phi,sys,gy_spi_in);
            % gx_spi_in = gx_spi_in{1};
            % gy_spi_in = gy_spi_in{1};
            % 
            % [gx_spi_in_prephasing, times, amplitudes] = mr.makeExtendedTrapezoidArea('x', 0, gx_spi_in.first, gx_spi_out_rot.area, sys);
            % [gy_spi_in_prephasing, times, amplitudes] = mr.makeExtendedTrapezoidArea('y', 0, gy_spi_in.first, gy_spi_out_rot.area, sys);
            % [gx_rewinder_out2, times, amplitudes] = mr.makeExtendedTrapezoidArea('x', gx_spi_out_rot.last,gx_spi_in.first, 0, sys);
            % [gy_rewinder_out2, times, amplitudes] = mr.makeExtendedTrapezoidArea('y', gy_spi_out_rot.last,gy_spi_in.first, 0, sys);
        else

        end


        for slc_ind=slice_indices
            % seq.addBlock(mr.makeLabel('INC','SLC', 1));
            % seq.addBlock(mr.makeLabel('SET','SLC', s));
            slice_start_time = sum(seq.blockDurations);

            %%%%%% fat sat %%%%%%
            if isNewFatSat
                seq.addBlock(gp_r,gp_p,gp_s);
                seq.addBlock(rf_fs,gn_r,gn_p,gn_s);
            else
                seq.addBlock(rf_fs,gz_fs);
            end

            rf90.freqOffset   = gs_ex.amplitude*sliceThickness*(slc_ind-1-(Nslices-1)/2);
            rf90.phaseOffset  = rf90_phase-2*pi*rf90.freqOffset*mr.calcRfCenter(rf90); % compensate for the slice-offset induced phase
            rf180.freqOffset  = GSref.amplitude*sliceThickness*(slc_ind-1-(Nslices-1)/2);
            rf180.phaseOffset = rf180_phase-2*pi*rf180.freqOffset*mr.calcRfCenter(rf180); % compensate for the slice-offset induced phase

            %%%%%% rf90 %%%%%%
            seq.addBlock(rf90,gs_ex);
            rf90_blk_start_time = prevBlcokStartTime(seq);
            seq.addBlock(gzReph);

            %%%%%% left diff encoding, rf180, right diff encoding%%%%%%
            if bval_all_frames(iterFrame)<=bvalue0_threshold
                seq.addBlock(mr.makeDelay(delayTE1_b0));
                seq.addBlock(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
                t_start_crusher_left = prevBlcokStartTime(seq);
                seq.addBlock(rf180,GSref);
                seq.addBlock(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
                t_start_crusher_right = prevBlcokStartTime(seq);
                seq.addBlock(mr.makeDelay(delayTE2_b0));
                
                if isAlreadyShowBvalueDueToCrusher
                    % do nothing
                else
                    isAlreadyShowBvalueDueToCrusher = true; % only show once
                    
                    b0_small_delta = gz180_crusher_1.riseTime + gz180_crusher_1.flatTime;
                    b0_big_delta = t_start_crusher_right - t_start_crusher_left;
                    bvalue_due_to_crusher = calc_bval_trap(mr.convert(gz180_crusher_1.amplitude,'Hz/m','mT/m'),b0_small_delta,b0_big_delta,gz180_crusher_1.riseTime);
                    total_bval_crusher = 3*bvalue_due_to_crusher; % x,y,z         
                    fprintf('<strong>b-value due to crusher gradient is: %.3f s/mm^2</strong>\n',total_bval_crusher);
                end                
            else
                seq.addBlock(mr.makeDelay(delayTE1),mr.scaleGrad(gDiff_x,bFactor_scale(iterFrame)),mr.scaleGrad(gDiff_y,bFactor_scale(iterFrame)),mr.scaleGrad(gDiff_z,bFactor_scale(iterFrame)));
                seq.addBlock(rf180,GSref);
                rf180_blk_start_time = prevBlcokStartTime(seq);
                seq.addBlock(mr.makeDelay(delayTE2),mr.scaleGrad(gDiff_x,bFactor_scale(iterFrame)),mr.scaleGrad(gDiff_y,bFactor_scale(iterFrame)),mr.scaleGrad(gDiff_z,bFactor_scale(iterFrame)));
            end
            
            
            
            % Store previous trigger time
            if isFirstTrigger
                isFirstTrigger = false;
            else                
                prev_trigger_time = trigger_start_time;
            end            
            trigger_start_time = sum(seq.blockDurations);            
            % Calculate time difference between triggers
            if exist('prev_trigger_time', 'var')
                trigger_time_diff = trigger_start_time - prev_trigger_time;
                assert( round(trigger_time_diff*1e6 )==round(TR*1e6/Nslices), 'Time between triggers is not equal to TR/Nslices, expected %.2f us, got %.2f us', TR*1e6/Nslices, trigger_time_diff*1e6);                
            end
            
            
            
            seq.addBlock(trig);
            seq.addBlock(mr.makeDelay(delay_after_trig));


            %%%%%% acq %%%%%%
            if isNoGradInFirstEcho
                seq.addBlock(mrcustomized.rotate_xingwang('z',phi,sys,mr.scaleGrad(gx_spi_out,0),mr.scaleGrad(gy_spi_out,0),adc));
                adc_blk_start_time = prevBlcokStartTime(seq);
                if 1 == Nechoes % spoil at the end
                    seq.addBlock(mrcustomized.rotate_xingwang('z',phi,sys,mr.scaleGrad(gx_spoil,0),mr.scaleGrad(gy_spoil,0),mr.scaleGrad(gz_spoil,0)));
                else
                    seq.addBlock(mrcustomized.rotate_xingwang('z',phi,sys,mr.scaleGrad(gx_rewinder_out,0),mr.scaleGrad(gy_rewinder_out,0)));
                end
            else
                
                
                if 1 == Nechoes % spoil at the end
                    % seq.addBlock(mrcustomized.rotate_xingwang('z',phi,sys,gx_spoil,gy_spoil,gz_spoil));
                    for iter_outin = 1:N_out_in
                        % seq.addBlock(mrcustomized.rotate_xingwang('z',phi+(iter_outin-1)*golden_angle,sys,gx_spi_out,gy_spi_out,adc));
                        % seq.addBlock(mrcustomized.rotate_xingwang('z',phi+(iter_outin-1)*golden_angle,sys,gx_spi_in,gy_spi_in,adc))
                        if isDelayBetweenShots && shot_ind>1
                           seq.addBlock( mr.makeDelay(base_delay_between_shots*(shot_ind-1)) ); 
                        end                        
                        seq.addBlock(mrcustomized.rotate_xingwang('z',phi + (iter_outin-1)*angle_beteen_outin_pair,sys,gx_spi_out,gy_spi_out,adc));
%                         adc_blk_start_time = prevBlcokStartTime(seq);
                        if isSpiralIn
                            if isADCduringRewinder
                                seq.addBlock(mrcustomized.rotate_xingwang('z',phi + (iter_outin-1)*angle_beteen_outin_pair,sys,gx_spi_in,gy_spi_in,adc))
                            else
                                tmp = mrcustomized.rotate_xingwang('z',phi + (iter_outin-1)*angle_beteen_outin_pair,sys,gx_spi_in,gy_spi_in);
                                seq.addBlock(   mr.align('left',tmp{1},tmp{2}, ...
                                    'right',adc)   )
                            end
                        else
                            seq.addBlock(mrcustomized.rotate_xingwang('z',phi + (iter_outin-1)*angle_beteen_outin_pair,sys,gx_rewinder_out,gy_rewinder_out));
                        end
                    end
                else
                    seq.addBlock(mrcustomized.rotate_xingwang('z',phi + (iter_outin-1)*angle_beteen_outin_pair,sys,gx_rewinder_out,gy_rewinder_out));
                end
            end
%             assert( abs( adc_blk_start_time - (rf180_blk_start_time+rf180centerInclDelay) - TE/2 )     <   timingCheckTol, 'Timing between 180 and adc is not TE/2'  )
            main_echo_readout_end_time = sum(seq.blockDurations);

            %%%%%% later echoes %%%%%%
            for iter_echoes = 2:Nechoes
                prev_rf180_blk_start_time = rf180_blk_start_time;

                seq.addBlock(mr.makeDelay(delay_bewteen_RO_and_next180));
                if isDifferentCrusher
                    seq.addBlock(GS7_2,GXcrusher180_2,GYcrusher180_2);
                else
                    seq.addBlock(GS7,GXcrusher180_2,GYcrusher180_2);
                end
                seq.addBlock(rf180,GS4);
                rf180_blk_start_time = prevBlcokStartTime(seq);
                if isDifferentCrusher
                    seq.addBlock(GS5_2,GXcrusher180_2,GYcrusher180_2);
                else
                    seq.addBlock(GS5,GXcrusher180_2,GYcrusher180_2);
                end
                assert( abs(   rf180_blk_start_time - prev_rf180_blk_start_time - TE  ) < timingCheckTol, 'Time between 180 are not TE')

                seq.addBlock(mr.makeDelay(delayTE_laterEchoes));
                seq.addBlock(trig);
                seq.addBlock(mr.makeDelay(delay_after_trig));
                seq.addBlock(mrcustomized.rotate_xingwang('z',phi + (iter_outin-1)*angle_beteen_outin_pair,sys,gx_spi_out,gy_spi_out,adc));
                adc_blk_start_time = prevBlcokStartTime(seq);
                assert( abs( adc_blk_start_time - (rf180_blk_start_time+rf180centerInclDelay) - TE/2 )<timingCheckTol, 'Timing between 180 and adc is not TE/2'  )
                if iter_echoes == Nechoes % spoil at the end
                    seq.addBlock(mrcustomized.rotate_xingwang('z',phi + (iter_outin-1)*angle_beteen_outin_pair,sys,gx_spoil,gy_spoil,gz_spoil));
                else
                    seq.addBlock(mrcustomized.rotate_xingwang('z',phi + (iter_outin-1)*angle_beteen_outin_pair,sys,gx_rewinder_out,gy_rewinder_out));
                end
            end





            % seq.addBlock(mrcustomized.rotate_xingwang('z',phi,sys,gx_out,gy_out,adc));
            % seq.addBlock(mrcustomized.rotate_xingwang('z',phi,sys,gx_rewinder,gy_rewinder));

            % seq.addBlock(GS7_2);
            % seq.addBlock(rf180,GS4);
            % seq.addBlock(GS5_2);
            % seq.addBlock(mr.makeDelay(delay_after_2nd_rf180));
            % seq.addBlock(mrcustomized.rotate_xingwang('z',phi,sys,gx_rewinder_reverse,gy_rewinder_reverse));
            % seq.addBlock(mrcustomized.rotate_xingwang('z',phi,sys,gx_in,gy_in,adc));


            % seq.addBlock(mr.makeLabel('INC','SLC', 1));
            % seq.addBlock(mr.makeLabel('INC','SEG', 1));
%             if delay_after_1slice > 0
%                 seq.addBlock(mr.makeDelay(delay_after_1slice));
%             end
            slice_end_time = sum(seq.blockDurations);
            time_1_slice = slice_end_time - slice_start_time;
            delay_after_1slice = TR/Nslices - time_1_slice;
            assert(delay_after_1slice>=0, 'TR too short, at least, you should increase TR to %.2f seconds',time_1_slice*Nslices);
            delay_after_1slice = round(delay_after_1slice/sys.blockDurationRaster)*sys.blockDurationRaster;
            if ~isAlreadyShowDelayAfterOneSlice
                isAlreadyShowDelayAfterOneSlice = true;
                fprintf('<strong>Delay after each slice is %.2f ms\n</strong>',1000*delay_after_1slice)                        
            end
            seq.addBlock(mr.makeDelay(delay_after_1slice));
        end % iter slices
    end % iter shots
end % iter volume


num_imaging_blocks = numel(seq.blockEvents);
imaging_blocks = cell2mat(seq.blockEvents);
imaging_blocks = reshape(imaging_blocks,[],num_imaging_blocks);
num_adcs_imaging = sum(  imaging_blocks(6,:)~=0  ); % these are LARGE ADCs concatenated from small adcs (atom adc)
num_adcs_expected = num_bvolumes * Nshots * Nslices * Nechoes * N_out_in*N_spi_one_pair;
assert(isequal(num_adcs_imaging,num_adcs_expected), '#ADC estimation failed, you may have trouble when reshaping data');


if isJeffDuynStyleMeasTraj
    % measure traj, reference, Simple Correction Method for k-Space Trajectory
    % Deviations in MRI, https://www.sciencedirect.com/science/article/pii/S1090780798913969?via%3Dihub        
    delay_between_shots = delay_after_1slice+fatSatModuleDur;250e-3; % give some time for T1 recovery
    for iterFrame = 1:nFrames
        isImagingFrame = true;
        if iterFrame<=Ndummies % dummy
            isImagingFrame = false;
        end        
        if ~(  isImagingFrame && bval_all_frames(iterFrame)==bvalue0_threshold  ) % only do traj measurement for b=0
            continue;
        end        
        for shot_ind = 1:Nshots
            phi = (shot_ind-1) * rot_angle;  
            
            tmp_grad = mrcustomized.rotate_xingwang('z',phi ,sys, gx_spi_out,gy_spi_out);
            assert( numel (tmp_grad)==2              ,  'Assuming 2 elements in output')
            assert( strcmp(tmp_grad{1}.channel,'x')  ,  'Assuming 1st output is x gradient')
            assert( strcmp(tmp_grad{2}.channel,'y')  ,  'Assuming 2st output is y gradient')
            gx_spi_out_new = tmp_grad{1};
            gy_spi_out_new = tmp_grad{2};
                        
            tmp_grad = mrcustomized.rotate_xingwang('z',phi ,sys, gx_rewinder_out,gy_rewinder_out);
            assert( numel (tmp_grad)==2              ,  'Assuming 2 elements in output')
            assert( strcmp(tmp_grad{1}.channel,'x')  ,  'Assuming 1st output is x gradient')
            assert( strcmp(tmp_grad{2}.channel,'y')  ,  'Assuming 2st output is y gradient')
            gx_rewinder_out_new = tmp_grad{1};
            gy_rewinder_out_new = tmp_grad{2};                        
            
            for slc_ind = slice_indices
                slice_start_time = sum(seq.blockDurations);
                
                rf90.freqOffset   = gs_ex.amplitude*sliceThickness*(slc_ind-1-(Nslices-1)/2);
                rf90.phaseOffset  = rf90_phase-2*pi*rf90.freqOffset*mr.calcRfCenter(rf90); % compensate for the slice-offset induced phase
                rf180.freqOffset  = GSref.amplitude*sliceThickness*(slc_ind-1-(Nslices-1)/2);
                rf180.phaseOffset = rf180_phase-2*pi*rf180.freqOffset*mr.calcRfCenter(rf180); % compensate for the slice-offset induced phase
                
                %%%%%%%%%%%%%%%%%%%%%%%% traj measurement %%%%%%%%%%%%%%%%%%%%%%%%
                % x slice select, x grad
                % y slice select, y grad
                % x slice select, x grad=0
                % y slice select, y grad=0
                channels = {'x', 'y', 'x', 'y'};
                gx_scales = [1,   0,   0,   0];
                gy_scales = [0,   1,   0,   0];
                for i = 1:4
                    current_channel = channels{i};
                    tmp_gx_scale    = gx_scales(i);
                    tmp_gy_scale    = gy_scales(i);
                    
                    gs_ex.channel  = current_channel;
                    gzReph.channel = current_channel;
                    GSref.channel  = current_channel;
                    
                    % fat sat
%                     if isNewFatSat
%                         seq.addBlock(gx_fs_pre,gy_fs_pre,gz_fs_pre);
%                         seq.addBlock(rf_fs,gx_fs_post,gy_fs_post,gz_fs_post);
%                     else
%                         seq.addBlock(rf_fs,gz_fs);
%                     end
                    seq.addBlock(rf90,gs_ex);
                    rf90_blk_start_time = prevBlcokStartTime(seq);
                    seq.addBlock(gzReph);
                    % left diff encoding, rf180, right diff encoding
                    seq.addBlock(mr.makeDelay(delayTE1),mr.scaleGrad(gDiff_x,bFactor_scale(iterFrame)),mr.scaleGrad(gDiff_y,bFactor_scale(iterFrame)),mr.scaleGrad(gDiff_z,bFactor_scale(iterFrame)));
                    seq.addBlock(rf180,GSref);
                    rf180_blk_start_time = prevBlcokStartTime(seq);
                    assert( abs(   (rf180_blk_start_time+rf180centerInclDelay) - (rf90_blk_start_time+rf90centerInclDelay) - TE/2   ) < timingCheckTol, 'Time between RF90 and RF180 is not TE/2' );
                    seq.addBlock(mr.makeDelay(delayTE2),mr.scaleGrad(gDiff_x,bFactor_scale(iterFrame)),mr.scaleGrad(gDiff_y,bFactor_scale(iterFrame)),mr.scaleGrad(gDiff_z,bFactor_scale(iterFrame)));
                    seq.addBlock(trig);
                    seq.addBlock(mr.makeDelay(delay_after_trig));
                    % readout
                    seq.addBlock(mr.scaleGrad(gx_spi_out_new,     tmp_gx_scale), mr.scaleGrad(gy_spi_out_new,     tmp_gy_scale),adc);
                    adc_blk_start_time = prevBlcokStartTime(seq);
                    seq.addBlock(mr.scaleGrad(gx_rewinder_out_new,tmp_gx_scale), mr.scaleGrad(gy_rewinder_out_new,tmp_gy_scale)    );
                    assert( abs( adc_blk_start_time - (rf180_blk_start_time+rf180centerInclDelay) - TE/2 )     <   timingCheckTol, 'Timing between 180 and adc is not TE/2'  )
                    seq.addBlock(mr.makeDelay(delay_between_shots));
                    slice_end_time = sum(seq.blockDurations);
                    time_1_slice = slice_end_time - slice_start_time;
                end               
            end
        end
    end
end
num_imaging_blocks = numel(seq.blockEvents);
imaging_blocks = cell2mat(seq.blockEvents);
imaging_blocks = reshape(imaging_blocks,[],num_imaging_blocks);
num_trig_total = sum(  imaging_blocks(7,:)~=0  );
assert(seq.version_major==1&&seq.version_minor==4&&seq.version_revision==1,'It seems blockEvents #7 == trig is only valid for pulseq version 1.4.x');
assert( numel(seq.extensionStringIDs)==1 && strcmp( seq.extensionStringIDs{1},'TRIGGERS' ) , 'Currently, the method to count triggers assume there is only extension')

%%
% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end
scanTime = sum(seq.blockDurations);
predictedScanTime = TR*Nshots*num_bvolumes;
% assert( abs(scanTime-predictedScanTime)<1, 'Scan time calculation is wrong' )
fprintf('<strong>Sequence duration is %02d:%02d\n</strong>', floor(scanTime/60),ceil(rem(scanTime, 60)));
% fprintf('<strong>Single shot duration is %02d:%02d\n</strong>', floor(scanTime/Nshots/60),ceil(rem(scanTime/Nshots, 60)));


nsegments_in_dummy  = 1;
n_triggers_in_dummy = Ndummies*Nslices*nsegments_in_dummy;
n_useless_triggers  = n_triggers_in_dummy;
skope_msg = [];
t = sprintf('<strong>---skope: number of useless triggers=%d</strong>\n',n_useless_triggers);
skope_msg = [skope_msg, t]; 
t = sprintf('<strong>---skope: number of triggers=%d</strong>\n',Nslices*Nshots*nImagingFrame);
skope_msg = [skope_msg, t]; 
t = sprintf('<strong>---skope: time betweem trigger and end of main readout=%.2f ms</strong>\n',1e3*(main_echo_readout_end_time-trigger_start_time));
skope_msg = [skope_msg, t];
t = sprintf('<strong>---skope: number of triggers (including trig from Jeff Duyn method)=%d</strong>\n',num_trig_total);
skope_msg = [skope_msg, t];
fprintf(skope_msg);


if isTestRun
   fprintf('<strong>Test run, set Nslices, num_directions_per_b to 1</strong>\n');
end

%%
if 0
    rep = seq.testReport;
    fprintf([rep{:}]); % as for January 2019 TR calculation fails for fat-sat
end

%%
seq.setDefinition('FOV', [fov fov sliceThickness*Nslices]);
seq.setDefinition('Name', 'spider');
seq.setDefinition('MaxAdcSegmentLength', adcSamplesPerSegment); % this is important for making the sequence run automatically on siemens scanners without further parameter tweaking
seq.setDefinition('NumShots',Nshots);
seq.setDefinition('NumArms',Narms);
seq.setDefinition('NumSlices',Nslices);
seq.setDefinition('MatrixSize',Nx);
seq.setDefinition('bValues',bvalues);
% seq.setDefinition('TR',TR); % may error if Siemens interpreter assumes a fixed TR
% seq.setDefinition('TE',TE);
seq.setDefinition('RepetitionTime',TR);
seq.setDefinition('EchoTime',TE);
seq.setDefinition('Resolution',fov/Nx);
seq.setDefinition('bTable',bTable);
seq.setDefinition('numSamplesPerArm',adc.numSamples);
seq.setDefinition('num_bvolumes',num_bvolumes);
seq.setDefinition('Nechoes',Nechoes);
seq.setDefinition('slice_indices',slice_indices);
seq.setDefinition('sys_type',sys_type);
seq.setDefinition('grad_safety_magrin',grad_safety_magrin);
seq.setDefinition('slew_safety_magrin',slew_safety_magrin);
seq.setDefinition('spiral_slew_safety_margin',spiral_slew_safety_margin);
seq.setDefinition('spiral_grad_safety_margin',spiral_grad_safety_margin);
seq.setDefinition('lowPNS_slew_safety_margin',lowPNS_slew_safety_margin);
seq.setDefinition('lowPNS_grad_safety_margin',lowPNS_grad_safety_margin);
seq.setDefinition('diff_grad_safety_margin',diff_grad_safety_margin);
seq.setDefinition('diff_slew_safety_margin',diff_slew_safety_margin);
seq.setDefinition('bvec',bvec);
seq.setDefinition('bval',bval);
seq.setDefinition('low_res_img_resolution',low_res_img_resolution);
seq.setDefinition('num_adcs_imaging',num_adcs_imaging);
seq.setDefinition('N_out_in',N_out_in);
seq.setDefinition('N_spi_one_pair',N_spi_one_pair);
seq.setDefinition('bvalue0_threshold',bvalue0_threshold);
seq.setDefinition('adcDwell',adc.dwell);
seq.setDefinition('base_delay_between_shots',base_delay_between_shots);
seq.setDefinition('isDelayBetweenShots',isDelayBetweenShots);
% seq.write(result_seqname);   % Output sequence for scanner

% the sequence is ready, so let's see what we got
% seq.plot('timeRange',[0,2*TR]);             % Plot sequence waveforms


%% forbbiden frequency calculation
% warning('Turning off warning for grad spectrum check')
% warning off
if isTestRun
    gradSpectrum(seq,sys,sys_type);
end
% warning on


%% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
ktraj_adc = ktraj_adc(:,1:adc.numSamples*num_adcs_imaging);
t_adc = t_adc(:,1:adc.numSamples*num_adcs_imaging);

% plot k-spaces
%figure; plot(t_ktraj, ktraj'); title('k-space components as functions of time'); % plot the entire k-space trajectory
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); title('2D k-space');

ktraj_adc_splited = reshape(ktraj_adc, 3, adc.numSamples,N_spi_one_pair,N_out_in,Nechoes,Nslices,Nshots,num_bvolumes);
t_adc_splited = reshape(t_adc, 1, adc.numSamples,N_spi_one_pair,N_out_in,Nechoes,Nslices,Nshots,num_bvolumes);
figure;
plot(ktraj(1,:),ktraj(2,:),'b');
hold on;
ind_echo = Nechoes;
ind_slice = Nslices;
ind_volume = num_bvolumes;
for ind_shot = 1:Nshots
    tmpx = ktraj_adc_splited(1,:,:,:,ind_echo,ind_slice,ind_shot,ind_volume);
tmpy = ktraj_adc_splited(2,:,:,:,ind_echo,ind_slice,ind_shot,ind_volume);
    plot(tmpx(:),tmpy(:),'r.')
end
xlim(1.2*Nx/2*1/fov*[-1 1])
ylim(1.2*Nx/2*1/fov*[-1 1])

kx = ktraj_adc_splited(1,:,:,:,1,1,1,1,1,1,1,1);
ky = ktraj_adc_splited(2,:,:,:,1,1,1,1,1,1,1,1);
t  = t_adc_splited(1,:,:,:,1,1,1,1,1,1,1,1,1);
kr = sqrt(kx.*kx+ky.*ky);
figure;
plot(kr(:))
figure;
plot3(t(:)-t(1),kx(:),ky(:))
xlabel 'Time (s)'
ylabel 'kx (1/m)'
zlabel 'ky (1/m)'
view([-17.3000,12.3829]);
title('All pairs in 1st shot');

% in 1st shot, show all out-in pairs
figure;
for ind_outin = 1:N_out_in   
    tmpx = ktraj_adc_splited(1,:,:,ind_outin,ind_echo,ind_slice,1,ind_volume);
    tmpy = ktraj_adc_splited(2,:,:,ind_outin,ind_echo,ind_slice,1,ind_volume);
    plot(tmpx(:),tmpy(:),'.')
    hold on;
end
title('All pairs in 1st shot');
legend_labels = arrayfun(@(x) sprintf('out-in pair %d', x), 1:N_out_in, 'UniformOutput', false);
legend(legend_labels);

% in 1st out-in pair, show all shots
figure;
for ind_shot = 1:Nshots
    tmpx = ktraj_adc_splited(1,:,:,1,ind_echo,ind_slice,ind_shot,ind_volume);
    tmpy = ktraj_adc_splited(2,:,:,1,ind_echo,ind_slice,ind_shot,ind_volume);
    plot(tmpx(:),tmpy(:),'.')
    hold on;
end
legend_labels = arrayfun(@(x) sprintf('shot %d', x), 1:Nshots, 'UniformOutput', false);
legend(legend_labels);
title('All shots in 1st out-in pair');


% in 1st out-in pair, show 1st shot
figure;
subplot(131)
for ind_shot = 1%:Nshots
    tmpx = ktraj_adc_splited(1,:,:,1,ind_echo,ind_slice,ind_shot,ind_volume);
    tmpy = ktraj_adc_splited(2,:,:,1,ind_echo,ind_slice,ind_shot,ind_volume);
    plot(tmpx(:),tmpy(:),'.')
    hold on;
end
legend_labels = arrayfun(@(x) sprintf('shot %d', x), 1:Nshots, 'UniformOutput', false);
legend(legend_labels);
title('1st shot in 1st out-in pair');
axis square;
subplot(132)
for ind_shot = 1%:Nshots
    tmpx = ktraj_adc_splited(1,:,1,1,ind_echo,ind_slice,ind_shot,ind_volume);
    tmpy = ktraj_adc_splited(2,:,1,1,ind_echo,ind_slice,ind_shot,ind_volume);
    plot(tmpx(:),tmpy(:),'.')
    hold on;
end
legend_labels = arrayfun(@(x) sprintf('shot %d', x), 1:Nshots, 'UniformOutput', false);
legend(legend_labels);
title('1st spi-out in 1st shot');
axis square;
subplot(133)
if isSpiralIn
    for ind_shot = 1%:Nshots
        tmpx = ktraj_adc_splited(1,:,2,1,ind_echo,ind_slice,ind_shot,ind_volume);
        tmpy = ktraj_adc_splited(2,:,2,1,ind_echo,ind_slice,ind_shot,ind_volume);
        plot(tmpx(:),tmpy(:),'.')
        hold on;
    end
    legend_labels = arrayfun(@(x) sprintf('shot %d', x), 1:Nshots, 'UniformOutput', false);
    legend(legend_labels);
    title('1st spi-in in 1st shot');
    axis square;
end

try
    % kx_1shot = ktraj_adc_splited(1,:,1,1,1,1);
    % ky_1shot = ktraj_adc_splited(2,:,1,1,1,1);
    % karea = voronoidens(kx_1shot+1i*ky_1shot);
    % kdens = 1./karea;
    % kdist = abs(kx_1shot+1i*ky_1shot);
    % 
    % index = find((kx_1shot(2:end).*kx_1shot(1:end-1))<0);
    % ky_PE = ky_1shot(index);
    % kx_PE = kx_1shot(index);
    % dist_PE = diff(kdist(index));
catch

end
% seq.install('siemens');
% return


%% calculate PNS and check whether we are still within limits
if strcmp(sys_type,'prisma')
    if ispc
        [pns,tpns]=seq.calcPNS('pathto.ascfile');
    else
        [pns,tpns]=seq.calcPNS('pathto.ascfile');
    end
elseif strcmp(sys_type,'prisma_XA30A')
    if ispc
        [pns,tpns]=seq.calcPNS('pathto.ascfile');
    else
        [pns,tpns]=seq.calcPNS('pathto.ascfile');
    end
elseif strcmp(sys_type,'skyra')
    if ispc
        [pns,tpns]=seq.calcPNS('pathto.ascfile');
    else
        [pns,tpns]=seq.calcPNS('pathto.ascfile');
    end
elseif strcmp(sys_type,'Connectome2') || strcmp(sys_type,'C2_simulate_prisma')
    if ispc
        [pns,tpns]=seq.calcPNS('pathto.ascfile');
    else
        [pns,tpns]=seq.calcPNS('pathto.ascfile');
    end
else
    if ~isGEscanner
        error('Undefined system')
    end
end
if max(tpns)>0.85
    warning('PNS=%.2f too high, the sequence may not run on the scanner',max(tpns))
end

%% helper functions
function b=bFactCalc(g, delta, DELTA)
% see DAVY SINNAEVE Concepts in Magnetic Resonance Part A, Vol. 40A(2) 39–65 (2012) DOI 10.1002/cmr.a
% b = gamma^2  g^2 delta^2 sigma^2 (DELTA + 2 (kappa - lambda) delta)
% in pulseq we don't need gamma as our gradinets are Hz/m
% however, we do need 2pi as diffusion equations are all based on phase
% for rect gradients: sigma=1 lambda=1/2 kappa=1/3
% for trapezoid gradients: TODO
sigma=1;
%lambda=1/2;
%kappa=1/3;
kappa_minus_lambda=1/3-1/2;
b= (2*pi * g * delta * sigma)^2 * (DELTA + 2*kappa_minus_lambda*delta);
end

function b = bFactCalc_trap(g, delta, DELTA, rup)
% handbook, page 287, example 9.2
b = (2*pi*g)^2 * (   delta^2*(DELTA-delta/3)  +  rup^3/30-delta*rup^2/6  );
end

function [amp,rup,small_delta ] =  findDiffGradAmp(gmax,slew,bvalue,gdiff_dur,big_delta,gradRasterTime)
% given timing, find the grad amp for a given b-value
%
% Input
%   gmax - unit is mT/m
%   slew - unit is mT/m/s
%   vbalue - unit is s/mm^2
%   for all the other time inputs, unit is second
%
% Output
%   amp - unit is mT/m

% Xingwang Yong, 20231022

% Sanity check
rup = ceil(gmax/slew/gradRasterTime)*gradRasterTime;
small_delta = gdiff_dur - rup;
small_delta = ceil(small_delta/gradRasterTime)*gradRasterTime;
maxbval = calc_bval_trap(gmax,small_delta,big_delta,rup);
if maxbval < bvalue
    error('Required b-value too large.')
end

% Bisection search
gupper = gmax;
glower = 0;
max_iter = 100;
iter = 0;
amp = 0;
while iter<max_iter
    iter = iter + 1;
    gmiddle = (gupper+glower)/2;

    rup = ceil(gmiddle/slew/gradRasterTime)*gradRasterTime;
    small_delta = gdiff_dur - rup;
    small_delta = ceil(small_delta/gradRasterTime)*gradRasterTime;
    bval = calc_bval_trap(gmiddle,small_delta,big_delta,rup);

    if abs(bval-bvalue) < 0.1
        amp = gmiddle;
        break;
    end

    if bval > bvalue
        gupper = gmiddle;
    else
        glower = gmiddle;
    end
end
if amp==0 % should never happen
    error('Something is wrong')
end
end

function bvalue = calc_bval_trap(grad_max,small_delta,big_delta,rup)
% grad unit is mT/m, time unit is seconds
% handbook, page 287, example 9.2
gamma = 2*pi*42.57e3;           % rad/mT
bvalue = (gamma.*grad_max/1000).^2*( small_delta^2*(big_delta-small_delta/3) + rup^3/30-small_delta*rup^2/6  );
end

function [gx_to_ksp_center,gy_to_ksp_center,tmps,tmpg,t] = grid_search_rewinder(gmax, smax, gx_start,gy_start, gx_area,gy_area,sys)
% grid search to find shortest time trapzoid rewinder, ghetto but works
%
% Input
%   grad unit is mT/m
%   slew unit is mT/m/s
%   units of all other parameters follow pulseq convention
%
% Output
%   ...
%

% Xingwang Yong, 20231103
% TODO, write doc

% sarr = 1:1:smax;
% garr = gmax:-1:1;
sarr = linspace(1,smax,20);
garr = linspace(1,gmax,20);
num_s = numel(sarr);
num_g = numel(garr);
result = cell(num_s,num_g,2);
sys_arr = cell(num_s,num_g);
for iter_s = 1:num_s
    sys.maxSlew = mr.convert(sarr(iter_s),'T/m/s');
    for iter_g = 1:num_g
        sys.maxGrad = mr.convert(garr(iter_g),'mT/m');
        sys_arr{iter_s,iter_g} = sys;
    end
end


time = 1e5*ones(num_s,num_g);
for iter_s = 1:num_s
    for iter_g = 1:num_g
        try
            gx_to_ksp_center = mr.makeExtendedTrapezoidArea('x',gx_start,0, gx_area, sys_arr{iter_s,iter_g});
            gy_to_ksp_center = mr.makeExtendedTrapezoidArea('y',gy_start,0, gy_area, sys_arr{iter_s,iter_g});
            result{iter_s,iter_g} = {gx_to_ksp_center,gy_to_ksp_center};
            time(iter_s,iter_g) = mr.calcDuration(gx_to_ksp_center,gy_to_ksp_center);
        catch
            result{iter_s,iter_g} = 0;
        end
    end
end


[t,I] = min(time(:));
[I_row, I_col] = ind2sub(size(time),I);
tmp_result = result{I_row,I_col};
gx_to_ksp_center = tmp_result{1};
gy_to_ksp_center = tmp_result{2};
tmps = sarr(I_row);
tmpg = garr(I_col);
end