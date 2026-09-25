% Dixon epi blip rewound acquisition (debra)
% debra with se-epi


% based on https://github.com/pulseq/pulseq/blob/master/matlab/demoSeq/writeEpiSpinEchoRS.m

addpath(genpath('/path to vdspiral'))
addpath(genpath('/path to safe_pns_prediction'))
addpath(genpath('path to pulseq v142, modifications for v151 and up in progress'))
addpath(genpath('FMRecon_SiemensRawData/utils'))

clear
close all

%% user input
isTestRun = 0; % for test only, set Nslices to 1 to accelerate script

% geometry
fov        = 224e-3;
Nx         = 408;
Nx_org     = Nx;
Ny         = Nx;   % Define FOV and resolution
thickness  = 2e-3; % slice thinckness
Nslices    = 30;
Nx_low_org = 64;
Nx_low     = Nx_low_org;
Ny_low_org = Nx_low_org;
Ny_low     = Ny_low_org;

% timing
TR       = 6;
TE       = 38e-3;
TE_low   = TE;
R        = 1;
R_low    = R;
RSegment = 6;
partFourierFactor = 6/8;  % partial Fourier factor: 1: full sampling 0: start with ky=0
partFourierFactor_fe        =  1;     % frequency-encoding partial Fourier factor
esp               = 0.75e-3; % the duration of gx
isOverwriteESP    = 0;
isAsymmEcho       = 0;
if isAsymmEcho
    isOverwriteESP = 1;
end
dTE_esp_ratio     = 2;
deltaTE           = esp*dTE_esp_ratio;
crusher_d         = 0.95e-3;
crusher_area             = 1e3; % unit is 1/m

% debra
kqShiftPeriod        = 3;
isDebra              = 0; % if set 0, then this sequence is epi
isProlongedDebraBlip = 0; % for temp test only, should always set to 0, see Erpeng Dai's 2022 IEEE Transactions on Medical Imaging. This blip duration is equal to echo spacing, instead of blip_dur
debra_nechoes        = 2;
isEchoTimeShift      = 1; % for multishot
if ~isDebra
    debra_nechoes=1;
end
isConventionalSegmentSpacing = 1; % in conventional multi segment epi, the ky distance between 2 lines is Ry*Rsegment. If this is set to 0, then distance is Ry
isKyShiftBetweenDiffFrames = 0;

% fat saturation
isFatSat         = 1;
isSSgradReversal = 0; % slice selection gradient reversal, SSGR. 

% dummy, ref, acs
isUseGreACS  = 0;
isRefscan    = 0; % for ghost correction
nDummyFrames = 1; % to get into steady state
is2ndEcho    = 0;

% debug options, preferred to be set to false
isWaterSat = 0; % for debug the fat not shift problem
if isWaterSat
    Nslices = 1;
    msg = 'water saturation, set #slices=1';
    warning(msg);
    warndlg(msg);
end


% diff
bvalue0_threshold    = 0; % when b-value smaller than this, crusher would be used instead of diffusion gradient
bvalues              = [bvalue0_threshold   1e3   ];
% Define acquisition order as a cell array
% Each row: {b-value index, bTable row indices or 0 for b=0}
% acquisition_order = {
%     1, 0;           % 1st: b=1, 1 direction (crusher, no bTable)
%     2, 1:5;        % 2nd: b=1000, 5 directions (bTable rows 1-5)
%     1, 0;           % 3rd: b=1, 1 direction (crusher, no bTable)
%     2, 6:10        % 4th: b=1000, 5 directions (bTable rows 6-10)
% };
acquisition_order = {
    1, 0;           % 1st: b=1, 1 direction (crusher, no bTable)
    2, 1:9;        % 2nd: b=1000, 5 directions (bTable rows 1-5)
    1, 0;
    2, 10:18;
    1, 0;
};


% rf
tRFex       = 3e-3;
if isSSgradReversal
    % tRFex = tRFex*1.5;
end
tRFref      = 3e-3;
spoilFactor = 1.5;             % spoiling gradient around the pi-pulse

% misc
pe_enable = 1;     % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
ro_os     = 1;     % oversampling factor (in contrast to the product sequence we don't really need it)
ro_os_low           =   1.00;           % oversampling factor (in contrast to the product sequence we don't really need it) 2.60 (if pe_pf: 2.70) (if pe_pf/TE60: 1.60)
pislquant = 10; % If it is not GE, this will be reset to 0 regardless of the setting here. For GE, number of shots/ADC events used for receive gain calibration.
trig_dur                 = 20e-6;
delay_after_trig         = 200e-6; % Gradient-free time for probe excitation [Unit: s] https://github.com/SkopeMagneticResonanceTechnologies/Pulseq-Sequences/blob/2576a44edb2f4123774565862472fe76ecd10435/sequences/PulseqBase.m#L102C9-L102C60
fatChemShift = 3.5e-6; % 3.5 ppm


%% system
sys_type                  = 'Connectome2'; % prisma, skyra, Connectome2, C2_simulate_prisma, trio, prisma_XA30A, premier
slew_safety_magrin        = 0.7;
grad_safety_magrin        = 0.9;
lowPNS_slew_safety_margin = 0.3;
lowPNS_grad_safety_margin = grad_safety_magrin;
diff_slew_safety_margin   = 0.45; % decrease this to reduce PNS, this would not lengthen TE too much
diff_grad_safety_margin   = 0.97;

if strcmp(sys_type,'prisma') || strcmp(sys_type,'C2_simulate_prisma') || strcmp(sys_type,'prisma_XA30A')
    physical_slew_max = 200;
    physical_grad_max = 80;
    B0=2.89; % 1.5 2.89 3.0
elseif strcmp(sys_type,'premier')
    physical_slew_max = 200;
    physical_grad_max = 70;%80;
    B0=3;
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

isGEscanner = strcmp(sys_type,'premier');
if ~isGEscanner
    pislquant = 0;
end
if isGEscanner
    % RF/gradient delay (sec).
    % Conservative choice that should work across all GE scanners.
    psd_rf_wait = 200e-6;  % section 5.4 in PulseqOnGE_v1.0.pdf
    
    rfDeadTime =  100e-6;
    rfRingdownTime = 60e-6 + psd_rf_wait;
    adcDeadTime = 20e-6;
    adcRasterTime = 2e-6;
    rfRasterTime = 2e-6;
    gradRasterTime = 4e-6;
    blockDurationRaster = 4e-6;
else % is siemens
    rfDeadTime =  100e-6;
    rfRingdownTime = 100e-6;
    adcDeadTime = 20e-6;
    %     adcRasterTime = 2e-6;
    adcRasterTime = 100e-9;
    rfRasterTime = 1e-6;
    gradRasterTime = 10e-6;
    blockDurationRaster = 10e-6;
end
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
lims = sys;

seq=mr.Sequence(sys);         % Create a new sequence object



fatOffresFreq = sys.gamma*sys.B0*fatChemShift; % Hz
% TE_ = 1/fatOffresFreq*[1 2]; % fat and water in phase for both echoes
shortest_TE_outofphase = 0.5*1/fatOffresFreq;
TE_outofphase = (1:2:500)*shortest_TE_outofphase;
[~,ind_TE_outofphase_nearCurrentTE] = min(abs(TE_outofphase-TE));
candidate_TEs = TE_outofphase(ind_TE_outofphase_nearCurrentTE-5:ind_TE_outofphase_nearCurrentTE+5);
fprintf('TEs that will lead to out-of-phase for 1st echo (unit is millisecond): \n%s\n',num2str(1e3*candidate_TEs) );
dTE_pi_difference = 0.5*1/fatOffresFreq;
if isOverwriteESP
    msg = 'Overwring esp and deltaTE';
    warndlg(msg);
    warning(msg);
    if dTE_esp_ratio==1
        deltaTE  = 1/3/  (fatChemShift*sys.B0*sys.gamma); % 1/3=2pi/3, 2=4pi/3
    elseif dTE_esp_ratio==2        
        deltaTE  = 2/3/  (fatChemShift*sys.B0*sys.gamma); % 1/3=2pi/3, 2=4pi/3
    else
        error('undefined')
    end
    esp = deltaTE/dTE_esp_ratio;    
    esp = round(esp/sys.gradRasterTime)*sys.gradRasterTime;
end
fprintf('<strong>ESP is %.2f ms\n</strong>\n',esp*1e3);




%% parse some inputs
crusher_d = round(crusher_d/sys.gradRasterTime)*sys.gradRasterTime;

if isAsymmEcho
    delay_time_pi_2_phase = 0.25/(fatChemShift*sys.gamma*sys.B0);
    delay_time_pi_2_phase = round(  delay_time_pi_2_phase /  sys.blockDurationRaster ) * sys.blockDurationRaster;
else
    delay_time_pi_2_phase = 0;
end

dkyScale_debra = 1;
if debra_nechoes==2 && dTE_esp_ratio==1
    error('This is not feasiable');
elseif debra_nechoes==2 && dTE_esp_ratio==2
    debra_shift_base = 3;
elseif debra_nechoes==2 && dTE_esp_ratio==3
    dkyScale_debra = 2;
    debra_shift_base = 4;
elseif debra_nechoes==3 && dTE_esp_ratio==1
    dkyScale_debra = 2;
    debra_shift_base = 1;
elseif debra_nechoes==3 && dTE_esp_ratio==2
    debra_shift_base = 5;
elseif debra_nechoes==3 && dTE_esp_ratio==3
    dkyScale_debra = 2;
    debra_shift_base = 7;
elseif debra_nechoes==1 % simple epi
    debra_shift_base = 1;
end


if isConventionalSegmentSpacing
    ky_dist_sameEcho_sameSeg = debra_nechoes*R*RSegment;
else
    ky_dist_sameEcho_sameSeg = debra_nechoes*R;
end


if isTestRun
    Nslices = 1;
%     num_directions_per_b = 1;
    msg = 'Test run, set Nslices to 1';
    warndlg(msg)
end


% bval and bvec
bTable      = xlsread('bvec_18dir.xlsx');
% tmpStruct = load('diffusion_6_directions.mat');
% bTable     = tmpStruct.diffusion_directions;
% bTable     = bTable.';
% bTable(1,:)=[]; % delete [0 0 0]
bTable(:,3) = -bTable(:,3);
% bTable = [0 0 1];warning('Debug, use 0 0 1 diff direction')
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



nFrames = numel(bval) + nDummyFrames;
if isRefscan
    nFrames = nFrames + 1;
end

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


%%
% Create fat-sat pulse
% B0=2.89; % 1.5 2.89 3.0
sat_ppm=-3.45;
sat_freq=sat_ppm*1e-6*B0*lims.gamma;
rf_fs = mr.makeGaussPulse(110*pi/180,'system',lims,'Duration',8e-3,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
rf_fs.phaseOffset=-2*pi*rf_fs.freqOffset*mr.calcRfCenter(rf_fs); % compensate for the frequency-offset induced phase
gz_fs = mr.makeTrapezoid('z',sys_lowPNS,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm



spoiler_amp = 3*8*42.58*10e2;
est_rise = 500e-6;
est_flat = 2500e-6;

gp_r = mr.makeTrapezoid('x','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',sys_lowPNS);
gp_p = mr.makeTrapezoid('y','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',sys_lowPNS);
gp_s = mr.makeTrapezoid('z','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',sys_lowPNS);

gn_r = mr.makeTrapezoid('x','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',sys_lowPNS);
gn_p = mr.makeTrapezoid('y','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',sys_lowPNS);
gn_s = mr.makeTrapezoid('z','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',sys_lowPNS);



% Create 90 degree slice selection pulse and gradient
[rf, gz, gzReph] = mr.makeSincPulse(pi/2,'system',sys_lowPNS,'Duration',tRFex,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4);

% Create 90 degree slice refocusing pulse and gradients
[rf180, gz180] = mr.makeSincPulse(pi,'system',sys_lowPNS,'Duration',tRFref,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4,'PhaseOffset',pi/2,'use','refocusing');
if isSSgradReversal
    gz180 = mr.scaleGrad(gz180, -1);
end

gz180_crusher_1 = mr.makeTrapezoid('z',sys_lowPNS,'area',crusher_area);
gz180_crusher_2 = mr.makeTrapezoid('y',sys_lowPNS,'area',crusher_area);
gz180_crusher_3 = mr.makeTrapezoid('x',sys_lowPNS,'area',crusher_area);


if 0
%     [rf90_sinc, gz] = mr.makeSincPulse(pi/2,'system',sys,'Duration',3e-3,'use','excitation',...
%     'PhaseOffset',pi/2,'apodization',0.4,'timeBwProduct',4,'SliceThickness',thickness_mm*1e-3);
rf90_sinc = rf;
thickness_mm = thickness*1e3;

[bw,f0,M_xy_sta,F1]=mr.calcRfBandwidth(rf90_sinc );
[M_z,M_xy,F2]=mr.simRf(rf90_sinc);

%
figure; plot(F1,abs(M_xy_sta),F2,abs(M_xy),F2,M_z);
axis([f0-2*bw, f0+2*bw, -0.1, 1.2]);
legend({'M_x_ySTA','M_x_ySIM','M_zSIM'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('STA vs. simulation, flip angle 90°');

figure; plot(F2,atan2(abs(M_xy),M_z)/pi*180);
axis([f0-2*bw, f0+2*bw, -5, 100]);
xlabel('frequency offset / Hz');
ylabel('flip ange [°]');
legend({'SINC'});
grid on;
title('Achieved flip angle for the nominal 90° flip');

figure; plot(F2,real(M_xy),F2,imag(M_xy));
axis([f0-2*bw, f0+2*bw, -1.2, 1.2]);
legend({'M_xSIM','M_ySIM'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('Real and imag. parts of transverse magnetisation, 90° flip');

%
sl_th=mr.aux.findFlank(F2(end:-1:1)/gz.amplitude,M_xy(end:-1:1),0.5)-mr.aux.findFlank(F2/gz.amplitude,M_xy,0.5);
figure; plot(F2/gz.amplitude*1000,abs(M_xy),'LineWidth',1.5); title('simulated slice profile, 90° flip, SINC'); xlabel('through-slice pos, mm');
hold on; yline(0.5,'-.'); xline([-0.5 0.5]*thickness_mm,'--'); legend({'slice profile','half-amplitude line','desired thickness'});
fprintf('actual slice thickness : %.3f mm\n',sl_th*1e3);





% [rf180_sinc, gz] = mr.makeSincPulse(pi,'system',sys,'Duration',4e-3,'use','refocusing',...
%     'apodization',0.3,'timeBwProduct',6,'SliceThickness',thickness_mm*1e-3);
rf180_sinc = rf180;
gz = gz180;

[bw,f0,M_xy_sta,F1]=mr.calcRfBandwidth(rf180_sinc);
[M_z,M_xy,F2,ref_eff]=mr.simRf(rf180_sinc);


figure; plot(F2,atan2(abs(M_xy),M_z)/pi*180);
axis([f0-2*bw, f0+2*bw, -5, 190]);
xlabel('frequency offset / Hz');
ylabel('flip ange [°]');
legend({'SINC'});
grid on;
title('Achieved flip angle for the nominal 180° flip');

figure; plot(F2,abs(ref_eff)); 
axis([f0-2*bw, f0+2*bw, -0.1, 1.1]);
xlabel('frequency offset / Hz');
ylabel('efficiency');
legend({'SINC'});
title('refocusing efficiency'); 

sl_th=mr.aux.findFlank(F2(end:-1:1)/gz.amplitude,ref_eff(end:-1:1),0.5)-mr.aux.findFlank(F2/gz.amplitude,ref_eff,0.5);
figure; plot(F2/gz.amplitude*1000,abs(ref_eff),'LineWidth',1.5); title('simulated slice profile, 180° flip, SINC'); xlabel('through-slice pos, mm');
hold on; yline(0.5,'-.'); xline([-0.5 0.5]*thickness_mm,'--'); legend({'slice profile','half-amplitude line','desired thickness'});
xlim([-10 10]); ylim([0 1.05]);
fprintf('actual slice thickness : %.3f mm\n',sl_th*1e3);

figure; plot(F2,angle(ref_eff)); 
axis([f0-2*bw, f0+2*bw, -3.2, 3.2]);
xlabel('frequency offset / Hz');
ylabel('efficiency');
legend({'SINC'});
title('refocusing efficiency phase (~2x RF phase)'); 

end





% [~, gzr_t, gzr_a]=mr.makeExtendedTrapezoidArea('z',gz180.amplitude,0,-gzReph.area+0.5*gz180.amplitude*gz180.fallTime,lims);
% gz180n=mr.makeExtendedTrapezoid('z','system',lims,'times',[0 gz180.riseTime gz180.riseTime+gz180.flatTime+gzr_t]+gz180.delay, 'amplitudes', [0 gz180.amplitude gzr_a]);
% [~, gzr1_t, gzr1_a]=mr.makeExtendedTrapezoidArea('z',0,gz180.amplitude,spoilFactor*gz.area,sys_lowPNS);
% [~, gzr2_t, gzr2_a]=mr.makeExtendedTrapezoidArea('z',gz180.amplitude,0,-gzReph.area+spoilFactor*gz.area,sys_lowPNS);
% if gz180.delay>(gzr1_t(4)-gz180.riseTime)
%     gz180.delay=gz180.delay-(gzr1_t(4)-gz180.riseTime);
% else
%     rf180.delay=rf180.delay+(gzr1_t(4)-gz180.riseTime)-gz180.delay;
%     gz180.delay=0;
% end
% gz180n=mr.makeExtendedTrapezoid('z','system',sys_lowPNS,'times',[gzr1_t gzr1_t(4)+gz180.flatTime+gzr2_t]+gz180.delay, 'amplitudes', [gzr1_a gzr2_a]);
% gzCrusher = mr.makeTrapezoid('z',sys_lowPNS,'area',spoilFactor*gz.area);

% define the output trigger to play out with every slice excitatuion
trig=mr.makeDigitalOutputPulse('osc0','duration', 100e-6); % possible channels: 'osc0','osc1','ext1'
% trigger for skope
trig = mr.makeDigitalOutputPulse('ext1','duration', trig_dur); % possible channels: 'osc0','osc1', ext1

% Define other gradients and ADC events
deltak  = 1/fov;
if isConventionalSegmentSpacing
    deltaky = RSegment*R*deltak*dkyScale_debra;      % Rsegement*R
    deltaky_low         =   RSegment*R_low*deltak;
else
    deltaky =          R*deltak*dkyScale_debra;
    deltaky_low         =   R_low*deltak;
end
kWidth  = Nx*partFourierFactor_fe*deltak;
kWidth_org          =   Nx_org*deltak;
kWidth_low          =   Nx_low*deltak;
kWidth_low_org      =   Nx_low_org*deltak;

% Phase blip in shortest possible time
% blip_dur = ceil(2*sqrt(deltaky/lims.maxSlew)/10e-6/2)*10e-6*2; % we round-up the duration to 2x the gradient raster time
if isProlongedDebraBlip
    % TODO, if isConventionalSegmentSpacing==0, here something may need to
    % be changed
    blip_dur            =   ceil(2*sqrt(RSegment*R*1*deltak/sys.maxSlew)/sys.gradRasterTime/2)*sys.gradRasterTime*2; % round-up the duration to 2x the gradient raster time
else
    if isConventionalSegmentSpacing
        blip_area_debra = RSegment*R*debra_shift_base*deltak;
    else
        blip_area_debra = R*debra_shift_base*deltak;
    end
    largest_blip = max(abs(blip_area_debra),abs(deltaky));
    blip_dur     = ceil( 2*sqrt( largest_blip/sys.maxSlew )/sys.gradRasterTime/2 ) *sys.gradRasterTime*2; % round-up the duration to 2x the gradient raster time
end
% the split code below fails if this really makes a trpezoid instead of a triangle...
gy = mr.makeTrapezoid('y',lims,'Area',-deltaky,'Duration',blip_dur); % we use negative blips to save one k-space line on our way towards the k-space center
%gy = mr.makeTrapezoid('y',lims,'amplitude',deltak/blip_dur*2,'riseTime',blip_dur/2, 'flatTime', 0);

readoutTime=esp - blip_dur;
% fprintf('<strong>Readout bandwidth per pixel = 1/readout_duration is %.1f\n</strong>',1/readoutTime); % not accurate

% readout gradient is a truncated trapezoid with dead times at the beginnig
% and at the end each equal to a half of blip_dur
% the area between the blips should be defined by kWidth
% we do a two-step calculation: we first increase the area assuming maximum
% slewrate and then scale down the amlitude to fix the area
extra_area=blip_dur/2*blip_dur/2*lims.maxSlew; % check unit!;
gx = mr.makeTrapezoid('x',lims,'Area',kWidth+extra_area,'duration',readoutTime+blip_dur);
% flat_time/total_time = (1-ab_ratio)
ramp_ratio = 1-gx.flatTime/(gx.flatTime+gx.riseTime+gx.fallTime);
fprintf('<strong>Ramp ratio = (1-flat_time/total_time)*100%% = %.1f%%\n</strong>',...
    100*ramp_ratio);
actual_area=gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
gx.amplitude=gx.amplitude/actual_area*kWidth;
gx.area = gx.amplitude*(gx.flatTime + gx.riseTime/2 + gx.fallTime/2);
gx.flatArea = gx.amplitude*gx.flatTime;
fprintf('<strong>Readout Gx amplitue = %.2f mT/m\n</strong>',mr.convert(gx.amplitude,'Hz/m','mT/m'));
fprintf('<strong>Readout Gx slew     = %.2f T/m/s\n</strong>',1e-3*mr.convert(gx.amplitude,'Hz/m','mT/m')/gx.riseTime);


extra_area_org      =   blip_dur/2*blip_dur/2*lims.maxSlew; % check unit!;
%         gx_org              =   mr.makeTrapezoid('x',lims,'Area',kWidth_org+extra_area_org,'duration',readoutTime+blip_dur);
gx_org              =   mr.makeTrapezoid('x',lims,'Area',kWidth_org+extra_area_org);
actual_area_org     =   gx_org.area-gx_org.amplitude/gx_org.riseTime*blip_dur/2*blip_dur/2/2-gx_org.amplitude/gx_org.fallTime*blip_dur/2*blip_dur/2/2;
gx_org.amplitude    =   gx_org.amplitude/actual_area_org*kWidth_org;
gx_org.area         =   gx_org.amplitude*(gx_org.flatTime + gx_org.riseTime/2 + gx_org.fallTime/2);
gx_org.flatArea     =   gx_org.amplitude*gx_org.flatTime;



% calculate ADC
% we use ramp sampling, so we have to calculate the dwell time and the
% number of samples, which are will be qite different from Nx and
% readoutTime/Nx, respectively.
adcDwellNyquist=deltak/gx.amplitude/ro_os;
% round-down dwell time to 100 ns
% adcDwell=floor(adcDwellNyquist*1e7)*1e-7;
adcDwell=floor(adcDwellNyquist*(1/sys.adcRasterTime))*sys.adcRasterTime; % GE raster is not 1e-7
fprintf('<strong>adcDwell=%.2f us</strong>\n',adcDwell*1e6);
fprintf('<strong>Bandwidth per pixel = 1/adc_dwell/Nx = %.1f Hz/pixel\n</strong>',1/adcDwell/Nx);
if ~isGEscanner
    minimalAdcDwell = 1.3;
    tolerance = 1e-10; 
    assert(  adcDwell*1e6-minimalAdcDwell  >=-tolerance  ,'ADC dwell too small. Currently, on Siemens scanner, the minimal adcDwell is %.2f us, while the current adc dwell is %.2f us',minimalAdcDwell,adcDwell*1e6)
end
adcSamples=floor(readoutTime/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
% MZ: no idea, whether ceil,round or floor is better for the adcSamples...
adc = mr.makeAdc(adcSamples,sys,'Dwell',adcDwell,'Delay',blip_dur/2);
% realign the ADC with respect to the gradient
time_to_center=adc.dwell*((adcSamples-1)/2+0.5); % I've been told that Siemens samples in the center of the dwell period
adc.delay=round((gx.riseTime+gx.flatTime/2-time_to_center)*1e6)*1e-6; % we adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us
adc.delay=round((gx.riseTime+gx.flatTime/2-time_to_center)*(1/sys.rfRasterTime))*sys.rfRasterTime; % Xingwang: in the mr.checkTiming, this is aligned to raster of rf
% adc.delay=round((gx.riseTime+gx.flatTime/2-time_to_center)*(1/sys.gradRasterTime))*sys.gradRasterTime; % we adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us
% this rounding actually makes the sampling points on odd and even readouts
% to appear misalligned. However, on the real hardware this misalignment is
% much stronger anyways due to the grdient delays

% FOV positioning requires alignment to grad. raster... -> TODO


% 2nd echo
gz180_c1            =   mr.makeTrapezoid('z',lims,'Amplitude',0*19.05*42.58*10e2,'Duration',crusher_d); % 2nd ref % v10
gz180_c2            =   mr.makeTrapezoid('y',lims,'Amplitude',-1*19.05*42.58*10e2,'Duration',crusher_d); % 2nd ref % v10
gz180_c3            =   mr.makeTrapezoid('x',lims,'Amplitude',1*19.05*42.58*10e2,'Duration',crusher_d); % 2nd ref % v10

blip_dur_low = ceil(2*sqrt(deltaky_low/sys.maxSlew)/sys.gradRasterTime/2)*sys.gradRasterTime*2; % v3
gy_low       = mr.makeTrapezoid('y',sys,'Area',-deltaky_low,'Duration',blip_dur_low); % v3

readoutTime_low = R_low/R*(readoutTime+blip_dur)-blip_dur_low;

extra_area_low      =   blip_dur_low/2*blip_dur_low/2*lims.maxSlew; % check unit!;
gx_low              =   mr.makeTrapezoid('x',lims,'Area',kWidth_low+extra_area_low,'duration',readoutTime_low+blip_dur_low);
actual_area_low     =   gx_low.area-gx_low.amplitude/gx_low.riseTime*blip_dur_low/2*blip_dur_low/2/2-gx_low.amplitude/gx_low.fallTime*blip_dur_low/2*blip_dur_low/2/2;
gx_low.amplitude    =   gx_low.amplitude/actual_area_low*kWidth_low;
gx_low.area         =   gx_low.amplitude*(gx_low.flatTime + gx_low.riseTime/2 + gx_low.fallTime/2);
gx_low.flatArea     =   gx_low.amplitude*gx_low.flatTime;

extra_area_low_org      =   blip_dur_low/2*blip_dur_low/2*lims.maxSlew; % check unit!;
gx_low_org              =   mr.makeTrapezoid('x',lims,'Area',kWidth_low_org+extra_area_low_org);
actual_area_low_org     =   gx_low_org.area-gx_low_org.amplitude/gx_low_org.riseTime*blip_dur_low/2*blip_dur_low/2/2-gx_low_org.amplitude/gx_low_org.fallTime*blip_dur_low/2*blip_dur_low/2/2;
gx_low_org.amplitude    =   gx_low_org.amplitude/actual_area_low_org*kWidth_low_org;
gx_low_org.area         =   gx_low_org.amplitude*(gx_low_org.flatTime + gx_low_org.riseTime/2 + gx_low_org.fallTime/2);
gx_low_org.flatArea     =   gx_low_org.amplitude*gx_low_org.flatTime;

adcDwellNyquist_low =   deltak/gx_low.amplitude/ro_os_low;
% round-down dwell time to 100 ns
adcDwell_low        =   floor(adcDwellNyquist_low*(1/sys.adcRasterTime))*sys.adcRasterTime;
adcSamples_low      =   floor(readoutTime_low/adcDwell_low/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
% MZ: no idea, whether ceil,round or floor is better for the adcSamples...
adc_low             =   mr.makeAdc(adcSamples_low,'Dwell',adcDwell_low,'Delay',blip_dur_low/2);
% realign the ADC with respect to the gradient
time_to_center_low  =   adc_low.dwell*((adcSamples_low-1)/2+0.5); % Siemens samples in the center of the dwell period
% adc_low.delay       =   round((gx_low.riseTime+gx_low.flatTime/2-time_to_center_low)*1e6)*1e-6; % we adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us
adc_low.delay       =   round((gx_low.riseTime+gx_low.flatTime/2-time_to_center_low)*(1/sys.rfRasterTime))*sys.rfRasterTime;  % Xingwang: in the mr.checkTiming, this is aligned to raster of rf
% adc_low.id          =   2; % v4 (for separate adc events)

Ny_low_pre                          =   round(Ny_low/2);  % PE steps prior to ky=0, excluding the central line
if isConventionalSegmentSpacing
    Ny_low_pre                          =   round(Ny_low_pre/RSegment/R_low);
else
    Ny_low_pre                          =   round(Ny_low_pre/1/R_low);
end
Ny_low_post                         =   round((partFourierFactor-1/2)*Ny_low+1); % PE lines after the k-space center including the central line
if isConventionalSegmentSpacing
    Ny_low_post                         =   round(Ny_low_post/RSegment/R_low);
else
    Ny_low_post                         =   round(Ny_low_post/1/R_low);
end
Ny_low_meas                         =   Ny_low_pre+Ny_low_post;




%% assembly blocks, #0 gre acs
segmentID = 1;
if isUseGreACS
    [seq,segmentID] = addacs(Nx, fov, thickness, sys, sys_lowPNS, adc, seq,segmentID,pislquant,Nslices);
end

% num_imaging_blocks = numel(seq.blockEvents);
% imaging_blocks = cell2mat(seq.blockEvents);
% imaging_blocks = reshape(imaging_blocks,[],num_imaging_blocks);
% num_adcs = sum(  imaging_blocks(6,:)~=0  ); % these are LARGE ADCs concatenated from small adcs (atom adc)


%% water sat pulses
cestPrepParms.tp      = 0.1;
cestPrepParms.td      = 0.005;
cestPrepParms.shape   = 'gauss'; % gauss, block, sinc
cestPrepParms.B1cwpe  = 2; % uT
cestPrepParms.isSpoil = true;%false;
cestPrepParms.npulses = 10;
satPulse = makeSaturationPulseFromCWPE(cestPrepParms.shape, cestPrepParms.B1cwpe, cestPrepParms.tp,cestPrepParms.td,sys);
% offsets_Hz = [192500	2000	0	32	-32	64	-64	96	-96	128	-128	192	-192	256	-256	256	-256	320	-320	320	-320	384	-384	384	-384	416	-416	416	-416	448	-448	448	-448	448	-448	448	-448	448	-448	448	-448	480	-480	480	-480	512	-512	512	-512	576	-576	640	-640	768	-768	1280	2560	3840	5120	6400	7680	8960	10240];
offsets_Hz = [	0];%	32	-32	64	-64	96	-96	128	-128	192	-192	256	-256	256	-256	320	-320	320	-320	384	-384	384	-384	416	-416	416	-416	448	-448	448	-448	448	-448	448	-448	448	-448	448	-448	480	-480	480	-480	512	-512	512	-512	576	-576	640	-640	768	-768	1280	2560	3840	5120	6400	7680	8960	10240];
seq.setDefinition('CESTfreqOffsets',offsets_Hz)

cest_prep_dur = cestPrepParms.npulses * ( cestPrepParms.tp + cestPrepParms.td );
if cestPrepParms.isSpoil
    % spoilers
    spoilRiseTime = 1e-3;
    spoilDuration = 4500e-6+ spoilRiseTime; % [s]
    % create pulseq gradient object
    [gxSpoil, gySpoil, gzSpoil] = makeSpoilerGradients(sys, spoilDuration, spoilRiseTime);
    spoil_dur = mr.calcDuration(gxSpoil, gySpoil, gzSpoil);
    
    cest_prep_dur = cest_prep_dur + spoil_dur;
end



%%
adc.id = seq.registerAdcEvent(adc);
adc_low.id = seq.registerAdcEvent(adc_low);



%% seq blocks
segid_delayTR = segmentID + 1;
segidbase = segmentID + 2;
slice_indices = tdr_sliceorder(Nslices,1);
adc_backup = adc;
adc_low_backup = adc_low;
adc_dummy  = mr.makeDelay(mr.calcDuration(adc));
adc_low_dummy = mr.makeDelay(mr.calcDuration(adc_low));
gx_backup  = gx;
all_kyShiftFactors_frame = [];
% all_kyShiftFactors_frame = zeros(nImagingFrame,1);
tic
isAlreadyExecuted = false;
isAlreadyShowDelayAfterOneSlice = false;
isAlreadyShowBvalueDueToCrusher = false;
isFirstTrigger = true;
for iterFrame = 1:nFrames
    bScale = bFactor_scale(iterFrame);
    
    isImagingFrame = true;    
    if iterFrame<=nFrames-nImagingFrame
        isImagingFrame = false;
    end
    
    if isKyShiftBetweenDiffFrames && isImagingFrame        
        kyShiftFactor_frame = R*mod(iterFrame-(nFrames-nImagingFrame),kqShiftPeriod);
        kyShiftFactor_frame = 1*mod(iterFrame-(nFrames-nImagingFrame),kqShiftPeriod); % FIXME, if only 1 seg, set this to the segment shift, to have unique and uniform shift across diffusion frame        
    else
        kyShiftFactor_frame = 0;
    end
    all_kyShiftFactors_frame(end+1) = kyShiftFactor_frame;
    
    
    
    
    % phase encoding and partial Fourier
    Ny_pre  = round((partFourierFactor-1/2)*Ny-1);  % PE steps prior to ky=0, excluding the central line
    if isConventionalSegmentSpacing
        Ny_pre  = round(Ny_pre/RSegment/R);
        Ny_post = round(Ny/2+1); % PE lines after the k-space center including the central line
        Ny_post = round(Ny_post/RSegment/R);
    else
        Ny_pre  = round(Ny_pre/   1     /R);
        Ny_post = round(Ny/2+1); % PE lines after the k-space center including the central line
        Ny_post = round(Ny_post/  1     /R);
    end
    Ny_meas = Ny_pre+Ny_post;
    
    if isImagingFrame
        number_of_shots = RSegment;
    else % dummy or ref, a single shot
        number_of_shots = 1;
    end
    
    for iterSeg = 1:number_of_shots
        gx = gx_backup;
        
        % Pre-phasing gradients
        gxPre = mr.makeTrapezoid('x',lims,'Area',-gx_org.area/2);
        if 1
            base_shift_between_seg = round(ky_dist_sameEcho_sameSeg/RSegment); % unique and uniform
            assert(base_shift_between_seg>=1,'shift too small, no shift between segments!')            
        else
            base_shift_between_seg = 1; % beneficial for the uniqueness of the ky lines            
        end
        gyPre  = mr.makeTrapezoid('y',sys,'Area',Ny_pre*deltaky/dkyScale_debra-(iterSeg-1)*base_shift_between_seg*deltak-kyShiftFactor_frame*deltak); 
        %         gyPre  = mr.makeTrapezoid('y',sys,'Area',(Ny_pre*deltaky/dkyScale_debra-(Nmulti-1)*deltak));warning('tmp'); % beneficial for the uniqueness of the ky lines
        %         gyPre  = mr.makeTrapezoid('y',sys,'Area',(Ny_pre*deltaky-(Nmulti-1)*round(debra_nechoes*RSegment*R/3)*deltak)); % beneficial for the uniformity of the ky lines
        if iterFrame==1 && iterSeg==1
            gyPreFirst = gyPre;
        end
        
        
        [gxPre,gyPre]=mr.align('right',gxPre,'left',gyPre);
        % relax the PE prepahser to reduce stimulation
        gyPre = mr.makeTrapezoid('y',lims,'Area',gyPre.area,'Duration',mr.calcDuration(gxPre,gyPre));
        gyPre.amplitude=gyPre.amplitude*pe_enable;
        
        % split the blip into two halves and produnce a combined synthetic gradient
        gy_parts = mr.splitGradientAt(gy, blip_dur/2, lims);
        [gy_blipup, gy_blipdown, ~]=mr.align('right',gy_parts(1),'left',gy_parts(2),gx);
        gy_blipdownup=mr.addGradients({gy_blipdown, gy_blipup}, lims);
        
        % pe_enable support
        gy_blipup.waveform=gy_blipup.waveform*pe_enable;
        gy_blipdown.waveform=gy_blipdown.waveform*pe_enable;
        gy_blipdownup.waveform=gy_blipdownup.waveform*pe_enable;
        
        % Calculate delay times
        durationToCenter = (Ny_pre+0.5)*mr.calcDuration(gx);
        rfCenterInclDelay=rf.delay + mr.calcRfCenter(rf);
        rf180centerInclDelay=rf180.delay + mr.calcRfCenter(rf180);
        
        delayTE1=ceil((TE/2 - mr.calcDuration(rf,gz) + rfCenterInclDelay - mr.calcDuration(gzReph) -  rf180centerInclDelay)/lims.gradRasterTime)*lims.gradRasterTime;
        % we do not need to include ETS in delayTE2 calculation, since delayTE2 is defined in 1st shot, where there is not ETS
        delayTE2=ceil((TE/2 - mr.calcDuration(rf180,gz180) + rf180centerInclDelay - mr.calcDuration(trig) - delay_after_trig- durationToCenter)/lims.gradRasterTime)*lims.gradRasterTime;
        assert(delayTE1>=0);
        assert(delayTE2>=0);
        
        delayTE1_b0 = delayTE1 - mr.calcDuration(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
        delayTE2_b0 = delayTE2 - mr.calcDuration(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
        assert(delayTE1_b0>=0);
        assert(delayTE2_b0>=0);        
        
        [gxPre,gyPre]=mr.align('right',gxPre,'left',gyPre);       
        
        sign_gyPre = sign(gyPre.area);
        if isDebra && ~isequal(   sign(gyPre.area),sign(gy.area)  ) ...
                && ~(debra_nechoes==3 && dTE_esp_ratio==1) % this is a special case where gy and gypre should blip differently
            gy_blipdownup = mr.scaleGrad(gy_blipdownup,-1); % gyPre and gy should blip towards the same direction
            gy_blipdown   = mr.scaleGrad(gy_blipdown,-1);
            gy_blipup     = mr.scaleGrad(gy_blipup,-1);
        end
        
        if isDebra
            if isProlongedDebraBlip
                dur = esp;
            else
                dur = blip_dur;
            end
            gyBlip_debra = mr.makeTrapezoid('y',sys,'Area',-sign_gyPre*debra_shift_base*deltaky/dkyScale_debra,'Duration',dur);
            if debra_nechoes==3&&dTE_esp_ratio==1
                gyBlip_debra = mr.scaleGrad(gyBlip_debra,-1);
            end
        else
            gyBlip_debra = gy;
        end
        
        gyBlip_debra_parts = mr.splitGradientAt(gyBlip_debra, blip_dur/2, sys);
        [gyBlip_debra_up,gyBlip_debra_down,~]=mr.align('right',gyBlip_debra_parts(1),'left',gyBlip_debra_parts(2),gx);
        % now for inner echos create a special gy gradient, that will ramp down to 0, stay at 0 for a while and ramp up again
        gyBlip_debra_down_up=mr.addGradients({gyBlip_debra_down, gyBlip_debra_up}, sys);
        
        gyBlip_down_debraUp =mr.addGradients({gy_blipdown, gyBlip_debra_up}, sys);
        gyBlip_debraDown_up =mr.addGradients({gy_blipup, gyBlip_debra_down}, sys);
        
        small_delta=delayTE2-ceil(sys_diff.maxGrad/sys_diff.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
        big_delta=delayTE1+mr.calcDuration(rf180,gz180);
        
        %g=sqrt(bval(iterFrame)*1e6/bFactCalc(1,small_delta,big_delta)); % for now it looks too large!
        g=sqrt(maxbval*1e6/bFactCalc(1,small_delta,big_delta)); % for now it looks too large!
        gr=ceil(g/sys_diff.maxSlew/lims.gradRasterTime)*lims.gradRasterTime;
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
        
        % Calculate the echo time shift for multishot EPI (QL)
        actual_esp          = gx.riseTime + gx.flatTime + gx.fallTime;
        TEShift             = actual_esp/RSegment;
        TEShift             = round(TEShift/sys.gradRasterTime)*sys.gradRasterTime;%round(TEShift,5); % from Berkin: roundn didn't work for the latest matlab, changed to round (sign -/+) % v2
        TEShift_before_echo = (iterSeg-1)*TEShift;
%         if TEShift_before_echo == 0
%             TEShift_before_echo = 0.00001; % apply the minimum duration for the no delay case
%         end
        TEShift_after_echo  = (RSegment-(iterSeg-1))*TEShift;
%         dETS_before         = mr.makeDelay(TEShift_before_echo);
%         dETS_after          = mr.makeDelay(TEShift_after_echo);



        
        
        
        % dummy shots before turning on ADC, to reach steady state
        isDummyShot = iterFrame <= nDummyFrames;
        if isDummyShot
            adc = adc_dummy;
            adc_low = adc_low_dummy;
        else
            adc = adc_backup;
            adc_low = adc_low_backup;
        end
        
        % First frame is EPI calibration/reference scan (blips off)
        isRefShot = (iterFrame == nDummyFrames+1 ) & isRefscan;
        
        blipsOn             = ~isDummyShot & ~isRefShot;
        blipsOn = blipsOn + (blipsOn== 0)*eps;        % non-zero scaling so that the trapezoid shape is preserved in the .seq file, for GE
        gyPre               = mr.scaleGrad(gyPre,blipsOn);
        gy_blipup           = mr.scaleGrad(gy_blipup,blipsOn);
        gy_blipdown         = mr.scaleGrad(gy_blipdown,blipsOn);
        gy_blipdownup       = mr.scaleGrad(gy_blipdownup,blipsOn);
        gyBlip_down_debraUp = mr.scaleGrad(gyBlip_down_debraUp,blipsOn);
        gyBlip_debraDown_up = mr.scaleGrad(gyBlip_debraDown_up,blipsOn);
        gyBlip_debra_down   = mr.scaleGrad(gyBlip_debra_down,blipsOn);
        
                        
        for iterSlice=slice_indices
            gx = gx_backup;
            slice_start_time = sum(seq.blockDurations);
            
            rf.freqOffset=gz.amplitude*thickness*(iterSlice-1-(Nslices-1)/2);
            rf.phaseOffset=-2*pi*rf.freqOffset*mr.calcRfCenter(rf); % compensate for the slice-offset induced phase
            rf180.freqOffset=gz180.amplitude*thickness*(iterSlice-1-(Nslices-1)/2);
            rf180.phaseOffset=pi/2-2*pi*rf180.freqOffset*mr.calcRfCenter(rf180); % compensate for the slice-offset induced phase
            
            segID_due_to_diff_encoding = 1;%~isb0 * iterFrame;
            segID_due_to_shots = iterSeg; % porbably we will not need this since we can pass pge2.validate. but for safety we just use it, since we have a bunch of available TRID
            segmentID = segidbase + isDummyShot + 2*isRefShot + segID_due_to_diff_encoding + segID_due_to_shots;
            lblTRID = mr.makeLabel('SET', 'TRID', segmentID);
            
            if isWaterSat % for debug, supress water
                if isGEscanner
                    error('not ready for GE yet, needs to change label, segment id.')
                end
                
                % CEST prep
                satPulse.freqOffset = offsets_Hz; % set frequency offset of the pulse
                for np = 1:cestPrepParms.npulses
                    seq.addBlock(satPulse) % add sat pulse
                    if cestPrepParms.td > 0
                        seq.addBlock(mr.makeDelay(cestPrepParms.td)); % add delay
                    end
                end
                if cestPrepParms.isSpoil % spoiling before readout
                    seq.addBlock(gxSpoil,gySpoil,gzSpoil);
                end
            end % of water saturation
            
            
            
            
            if isFatSat
%                 seq.addBlock(rf_fs,gz_fs,lblTRID);
                seq.addBlock(gp_r,gp_p,gp_s,lblTRID);
                seq.addBlock(rf_fs,gn_r,gn_p,gn_s);
            end
            
            if isFatSat
                seq.addBlock(rf,gz);
            else
                seq.addBlock(rf,gz,lblTRID);
            end
            seq.addBlock(gzReph);
            
            
            if bval_all_frames(iterFrame)<=bvalue0_threshold
                seq.addBlock(mr.makeDelay(delayTE1_b0));
                seq.addBlock(gz180_crusher_1,gz180_crusher_2,gz180_crusher_3);
                t_start_crusher_left = prevBlcokStartTime(seq);
                seq.addBlock(rf180,gz180);
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
                seq.addBlock(mr.makeDelay(delayTE1),mr.scaleGrad(gDiff_x,bScale),mr.scaleGrad(gDiff_y,bScale),mr.scaleGrad(gDiff_z,bScale));
                seq.addBlock(rf180,gz180);
                seq.addBlock(mr.makeDelay(delayTE2),mr.scaleGrad(gDiff_x,bScale),mr.scaleGrad(gDiff_y,bScale),mr.scaleGrad(gDiff_z,bScale));
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
            
            if isEchoTimeShift && TEShift_before_echo>0
                seq.addBlock(mr.makeDelay(TEShift_before_echo)); % echotimeshift for multishot
            end
            
            gyPreScale = gyPre.area / gyPreFirst.area;
            % seq.addBlock(gxPre,gyPre);
            seq.addBlock(gxPre,mr.scaleGrad(gyPreFirst,gyPreScale)); % for GE
            if isAsymmEcho
                seq.addBlock(   mr.makeDelay(delay_time_pi_2_phase)   ); % delay the readout, water-fat pi/2 difference between water and fat, good for IDEAL algorithm
            end
            for i=1:Ny_meas
                if i==1
                    seq.addBlock(gx,gy_blipup,adc); % Read the first line of k-space with a single half-blip at the end
                elseif i==Ny_meas
                    if mod(i,debra_nechoes)==1
                        seq.addBlock(gx,gyBlip_debra_down,adc);
                    else
                        seq.addBlock(gx,gy_blipdown,adc); % Read the last line of k-space with a single half-blip at the beginning
                    end
                else
                    if mod(i,debra_nechoes)==0
                        seq.addBlock(gx,gyBlip_down_debraUp,adc);
                    elseif mod(i,debra_nechoes)==1
                        seq.addBlock(gx,gyBlip_debraDown_up,adc);
                    elseif mod(i,debra_nechoes)==2
                        seq.addBlock(gx,gy_blipdownup,adc);
                    end
                end
                gx = mr.scaleGrad(gx,-1);   % Reverse polarity of read gradient
            end
            
            main_echo_readout_end_time = sum(seq.blockDurations);
            
            
            if isEchoTimeShift && TEShift_after_echo>0
                if isGEscanner
                    error('FIXME, if delay==0, GE segment violation?');
                end
                seq.addBlock(mr.makeDelay(TEShift_after_echo))
            end 
            
            
if is2ndEcho
    
            isDebug2ndEcho = 0;
            
            % 2nd echo
            if isDebug2ndEcho
                [ktraj_adc1, t_adc1, ktraj1, t_ktraj1, t_excitation1, t_refocusing1] = seq.calculateKspacePP();
                ind1 = size(ktraj_adc1,2);
                ind2 = size(ktraj1,2);
            end
            
            gy_parts_low                        =   mr.splitGradientAt(gy_low, blip_dur_low/2, lims);
            [gy_blipup_low, gy_blipdown_low,~]  =   mr.align('right',gy_parts_low(1),'left',gy_parts_low(2),gx_low);
            gy_blipdownup_low                   =   mr.addGradients({gy_blipdown_low, gy_blipup_low}, lims);
            % pe_enable support
            gy_blipdown_low.waveform            =   gy_blipdown_low.waveform*pe_enable;
            gy_blipup_low.waveform              =   gy_blipup_low.waveform*pe_enable;
            gy_blipdownup_low.waveform          =   gy_blipdownup_low.waveform*pe_enable;            
            
            gxPre_low                           =   mr.makeTrapezoid('x',lims,'Area',-gx_low_org.area/2);
            gyPre_low                           =   mr.makeTrapezoid('y',lims,'Area',(Ny_low_pre*deltaky_low-(iterSeg-1)*R_low*deltak));
            if iterFrame==1 && iterSeg==1 && iterSlice==slice_indices(1)
                gyPre_low_first = gyPre_low;
            end
            [gxPre_low,gyPre_low]               =   mr.align('right',gxPre_low,'left',gyPre_low);
            gyPre_low                           =   mr.makeTrapezoid('y',lims,'Area',gyPre_low.area,'Duration',mr.calcDuration(gxPre_low,gyPre_low));
            gyPre_low.amplitude                 =   gyPre_low.amplitude*pe_enable;
            
            % v3
            gyPre_low_post_area                 =   -1*(Ny_low_post-0.5)*deltaky_low; % this seems to be correct
            
            if(mod(iterSeg,2) ~= 0)
                gyPre_low_post_area             =   (Ny_low_post-1)*deltaky_low; % this seems to be correct
            end
            
            gxPre_low_post_area_nav             =   -1*(0.5*gx_low.flatArea + gx_low.amplitude*gx_low.fallTime/2); % partFourierFactor should be larger than 0.5 % v9
            gxPre_low_post_area                 =   -1*(abs(gxPre_low_post_area_nav*2) - (0.5*gx_low_org.flatArea + gx_low_org.amplitude*gx_low_org.fallTime/2)); % partFourierFactor should be larger than 0.5 % v9
            gxPre_low_post_area_nav             =   -1*(abs(gxPre_low_post_area_nav*2) - (0.5*gx_low_org.flatArea + gx_low_org.amplitude*gx_low_org.fallTime/2)); % partFourierFactor should be larger than 0.5 % v9
            
            if(mod(Ny_low_meas,2) == 0)
                gxPre_low_post_area             =   0.5*gx_low_org.flatArea + gx_low_org.amplitude*gx_low_org.riseTime/2; % v9
            end
            
            gxPre_low_post                      =   mr.makeTrapezoid('x',lims,'Area', gxPre_low_post_area, 'Duration', 0.0012); % Fixme: gxPre_post_area is still not centered % v9
            gyPre_low_post                      =   mr.makeTrapezoid('y',lims,'Area', gyPre_low_post_area, 'Duration', 0.0012); % QL, XW
            % gzPre_low_post                      =   mr.makeTrapezoid('z',lims,'Area', -kz*(caipi_factor-1)/2, 'Duration', 0.0012); % v10
            
            % if(mod(Ny_low_meas,2) ~= 0)
            % gzPre_low_post.amplitude        =   -gzPre_low_post.amplitude;
            % end
            
            %             [gxPre_low_post,gyPre_low_post]     =   mr.align('right',gxPre_low_post,'left',gyPre_low_post); % v10_uc
            % [gxPre_low_post,gyPre_low_post,gzPre_low_post] =   mr.align('right',gxPre_low_post,'left',gyPre_low_post,gzPre_low_post); % v10
            [gxPre_low_post,gyPre_low_post] =   mr.align('right',gxPre_low_post,'left',gyPre_low_post);
            % relax the PE prepahser to reduce stimulation
            %             gyPre_low_post                      =   mr.makeTrapezoid('y',lims,'Area',gyPre_low_post.area,'Duration',mr.calcDuration(gxPre_low_post,gyPre_low_post)); % v10_uc
            gyPre_low_post                      =   mr.makeTrapezoid('y',lims,'Area',gyPre_low_post.area,'Duration',mr.calcDuration(gxPre_low_post,gyPre_low_post)); % v10
            
            gyPre_low_post.amplitude            =   gyPre_low_post.amplitude*pe_enable;
            % gzPre_low_post.amplitude            =   gzPre_low_post.amplitude*pe_enable; % v10
                        
            % post gradient to go back to k-space center            
            if mod(Ny_meas,2)==0
                gxPre_post = mr.scaleGrad(gxPre,-1);
            else
                gxPre_post = gxPre;
            end
            gxPre_post = mr.makeTrapezoid('x',sys_lowPNS,'Area',gxPre_post.area);
            
            
            gyPre_post_area = (Ny_post-0.5)*deltaky; % this seems to be correct
            %         gyPre_post      = mr.makeTrapezoid('y',sys,'Area', gyPre_post_area); % (QL, XW)
            if (~isDummyShot && ~isRefShot) && ~isAlreadyExecuted
                isAlreadyExecuted = true; % run once only
                [~, ~, ktraj1,] = seq.calculateKspacePP();
                gyAreaFirstImagingShot = ktraj1(2,end);
                gyPreFirst_imaging = gyPre;
            end
            if ~isDummyShot && ~isRefShot
                gyShiftAreaRelativeToFirstGy = gyPre.area - gyPreFirst_imaging.area;
                gyPre_post_area = -gyAreaFirstImagingShot-gyShiftAreaRelativeToFirstGy;
                gyPre_post_scale = gyPre_post_area./gxPre_post.area;
            else
                gyPre_post_scale = eps;
            end
            gyPre_post      = mr.scaleGrad(gxPre_post,gyPre_post_scale);
            gyPre_post.channel = 'y';
% % %             gyPre_post  = mr.makeTrapezoid('y',sys_lowPNS,'Area',gyPre_post.area); % not compatiable for GE, we can first re-make gxPre_post, and scale corresponding gy

            % Calculate delay times
            durationToCenter        = (Ny_pre + 0.5) * mr.calcDuration(gx);
            durationToCenter2       = (Ny_post + 0.5) * mr.calcDuration(gx); % v4
            durationToCenter_low    = (Ny_low_pre + 0.5) * mr.calcDuration(gx_low); % v4
            
            delayTE3                = ceil((TE_low/2 - mr.calcDuration(rf180,gz180) + rf180centerInclDelay - durationToCenter2)/lims.gradRasterTime)*lims.gradRasterTime; % v4 % v10_uc
%             delayTE3                = ceil((TE_low/2 -  rf180centerInclDelay - durationToCenter2)/lims.gradRasterTime)*lims.gradRasterTime; % v4 % v10_uc
            delayTE3                = delayTE3 - mr.calcDuration(gxPre_post,gyPre_post) - mr.calcDuration(gz180_c1);
            assert(delayTE3>=0);
            
            delayTE4                = ceil((TE_low/2 - mr.calcDuration(rf180,gz180) + rf180centerInclDelay - durationToCenter_low)/lims.gradRasterTime)*lims.gradRasterTime; % v10_uc
            delayTE4                = delayTE4 - mr.calcDuration(gxPre_low,gyPre_low); % v4 % v5_2
            delayTE4                = delayTE4 - mr.calcDuration(gz180_c1); % v4
            assert(delayTE4>=0);
            
            gyPre_low         = mr.scaleGrad(gyPre_low,blipsOn);
            gy_blipup_low     = mr.scaleGrad(gy_blipup_low,blipsOn);
            gy_blipdown_low   = mr.scaleGrad(gy_blipdown_low,blipsOn);
            gy_blipdownup_low = mr.scaleGrad(gy_blipdownup_low,blipsOn);
            
            
            seq.addBlock(gxPre_post,gyPre_post); % go back to k-space center
            if isDebug2ndEcho
                [~, ~, ktraj1] = seq.calculateKspacePP();
                assert(  abs(  ktraj1(1,end)  )   < 1e-7)
                assert(  abs(  ktraj1(2,end)  )   < 1e-7)
            end
            seq.addBlock(mr.makeDelay(delayTE3));
            seq.addBlock(gz180_c1,gz180_c2,gz180_c3);
            seq.addBlock(rf180,gz180); % (QL) it used to be gz180n
            seq.addBlock(gz180_c1,gz180_c2,gz180_c3);
            
            seq.addBlock(mr.makeDelay(delayTE4));
            
            % seq.addBlock(mr.align('left',gyPre_slt,'right',gxPre_low));
            gyPreScale_low = gyPre_low.area / gyPre_low_first.area; % scaling for GE
            seq.addBlock(  mr.align('left',mr.scaleGrad(gyPre_low_first,gyPreScale_low),'right',gxPre_low)  ); % scaling for GE

            for i = 1:Ny_low_meas
                gx_low_temp = mr.scaleGrad(gx_low,(-1)^(i-1));
                if i==1
                    seq.addBlock(gx_low_temp,gy_blipup_low,adc_low); % Read the first line of k-space with a single half-blip at the end
                elseif i==Ny_low_meas
                    seq.addBlock(gx_low_temp,gy_blipdown_low,adc_low); % Read the last line of k-space with a single half-blip at the beginning
                else
                    seq.addBlock(gx_low_temp,gy_blipdownup_low,adc_low); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                end
            end
            %             seq.addBlock(gxPre_low_post,gyPre_post_slt); % go back to k-space center % v3 % v10%FIXME            
            
            if isDebug2ndEcho
                [ktraj_adc1, t_adc1, ktraj1, t_ktraj1, t_excitation1, t_refocusing1] = seq.calculateKspacePP();
                
                % plot k-spaces
                figure; plot(ktraj1(1,ind2+1:end),ktraj1(2,ind2+1:end),'b'); % a 2D plot
                axis('equal'); % enforce aspect ratio for the correct trajectory display
                hold;plot(ktraj_adc1(1,ind1+1:end),ktraj_adc1(2,ind1+1:end),'r.'); % plot the sampling points
                yline(Ny_low_org/2*deltak,'k--');
                yline(-Ny_low_org/2*deltak,'k--');
            end
end % of is2ndEcho            
            



%             time_1_slice = logical(isWaterSat)*cest_prep_dur ...
%                 + logical(isFatSat)*mr.calcDuration(rf_fs,gz_fs) ...
%                 + rfCenterInclDelay...
%                 + TE...
%                 + mr.calcDuration(trig) + delay_after_trig...
%                 + logical(isEchoTimeShift)*(TEShift_before_echo+TEShift_after_echo)...
%                 + logical(isAsymmEcho)*delay_time_pi_2_phase;
%             if is2ndEcho
%                 time_1_slice = time_1_slice...
%                 + TE_low + Ny_low_post*mr.calcDuration(gx_low); % FIXME, still not very accurate
%             else
%                 time_1_slice = time_1_slice...
%                     +Ny_post*mr.calcDuration(gx);
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
            
        end % of slice loop
    end % of segment loop    
end % of frame loop, i.e. diffusion volume loop
tSeqGeneration = toc;
fprintf('Sequence generation took %.1f seconds\n', tSeqGeneration);
fprintf('<strong>Readout duration is Ny_meas*ESP = %.2f ms</strong>\n',1000*Ny_meas*esp);

nsegments_in_ref    = 1;
nsegments_in_dummy  = 1;
n_triggers_in_ref   = isRefscan*   Nslices*nsegments_in_ref;
n_triggers_in_dummy = nDummyFrames*Nslices*nsegments_in_dummy;
n_useless_triggers  = n_triggers_in_ref + n_triggers_in_dummy;
skope_msg = [];
t = sprintf('<strong>---skope: number of useless triggers=%d</strong>\n',n_useless_triggers);
skope_msg = [skope_msg, t]; 
t = sprintf('<strong>---skope: number of triggers=%d</strong>\n',Nslices*RSegment*nImagingFrame);
skope_msg = [skope_msg, t]; 
t = sprintf('<strong>---skope: time betweem trigger and end of main readout=%.2f ms</strong>\n',1e3*(main_echo_readout_end_time-trigger_start_time));
skope_msg = [skope_msg, t];
fprintf(skope_msg);
scanTime = sum(seq.blockDurations);
fprintf('<strong>Sequence duration is %02d:%02d\n</strong>', floor(scanTime/60),ceil(rem(scanTime, 60)));

%% estimate adc size
num_imaging_blocks = numel(seq.blockEvents);
imaging_blocks = cell2mat(seq.blockEvents);
imaging_blocks = reshape(imaging_blocks,[],num_imaging_blocks);
num_adcs = sum(  imaging_blocks(6,:)~=0  ); % these are LARGE ADCs concatenated from small adcs (atom adc)
num_adcs_expected = 32 + nFrames*Ny_meas;
% assert(isequal(num_adcs,num_adcs_expected), '#ADC estimation failed, you may have trouble when reshaping data');


%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end




%% prepare the sequence output for the scanner
seq.setDefinition('FOV', [fov fov Nslices*thickness]);
seq.setDefinition('Name', 'sedebra');
seq.setDefinition('Nx', Nx);
seq.setDefinition('Ny', Ny);
seq.setDefinition('Ry', R);
seq.setDefinition('Ry_low', R_low);
seq.setDefinition('deltaTE', deltaTE);
seq.setDefinition('Ny_meas', Ny_meas);
seq.setDefinition('partFourierFactor', partFourierFactor);
seq.setDefinition('EchoSpacing', esp);
seq.setDefinition('pislquant', pislquant);
seq.setDefinition('bvec',bvec);
seq.setDefinition('bval',bval);
seq.setDefinition('nImagingFrame',nImagingFrame);
seq.setDefinition('Nslices',Nslices);
seq.setDefinition('numSamples',adc.numSamples);
seq.setDefinition('RSegment',RSegment);
seq.setDefinition('isGEscanner',isGEscanner);
seq.setDefinition('isRefscan',isRefscan);
seq.setDefinition('debra_nechoes',debra_nechoes);
seq.setDefinition('isProlongedDebraBlip',isProlongedDebraBlip);
seq.setDefinition('isConventionalSegmentSpacing',isConventionalSegmentSpacing);
seq.setDefinition('isUseGreACS',isUseGreACS);
seq.setDefinition('dTE_esp_ratio',dTE_esp_ratio);
seq.setDefinition('kqShiftPeriod',kqShiftPeriod);
seq.setDefinition('all_kyShiftFactors_frame',all_kyShiftFactors_frame);
seq.setDefinition('bvalue0_threshold',bvalue0_threshold);
seq.setDefinition('isFatSat',isFatSat);
seq.setDefinition('isSSgradReversal',isSSgradReversal);
seq.setDefinition('Nx_low_org',Nx_low_org);
seq.setDefinition('Ny_low_org',Ny_low_org);
seq.setDefinition('numSamplesLow',adc_low.numSamples);
seq.setDefinition('Ny_low_meas',Ny_low_meas);
seq.setDefinition('Nx_low_org',Nx_low_org);
seq.setDefinition('Ny_low_org',Ny_low_org);
seq.setDefinition('dkyScale_debra',dkyScale_debra);
seq.setDefinition('isDebra',isDebra);
seq.setDefinition('base_shift_between_seg',base_shift_between_seg);
seq.setDefinition('adcDwell',adc.dwell);
seq.setDefinition('TEShift',TEShift);
seq.setDefinition('is2ndEcho',is2ndEcho);
seq.setDefinition('EchoTime',TE);
seq.setDefinition('ro_os',ro_os);

%% do some visualizations
block_lbl = seq.evalLabels('evolution','block');


if isTestRun
% seq.plot();             % Plot sequence waveforms

% trajectory calculation
[ktraj_adc1, t_adc1, ktraj1, t_ktraj1, t_excitation1, t_refocusing1] = seq.calculateKspacePP();
kmax = 1/fov*Nx/2;


% plot k-spaces
% figure; plot(t_ktraj1, ktraj1'); % plot the entire k-space trajectory
% hold on; plot(t_adc1,ktraj_adc1(1,:),'.'); % and sampling points on the kx-axis
figure; plot(ktraj1(1,:),ktraj1(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc1(1,:),ktraj_adc1(2,:),'r.'); % plot the sampling points
xline(kmax,'k--');
xline(-kmax,'k--');
yline(kmax,'k--');
yline(-kmax,'k--');


% if isUseGreACS
%     start_ind = 1 + adc_backup.numSamples*(32*Nslices+pislquant);
% else
%     start_ind = 1;
% end
% ktraj_adc_epi = ktraj_adc1(:,start_ind:end);
% 
% nsegments_in_ref = 1;
% number_of_elements_in_ref = Nslices*nsegments_in_ref*Ny_meas*adc_backup.numSamples;
% if is2ndEcho
%     number_of_elements_in_ref = number_of_elements_in_ref + Nslices*nsegments_in_ref*Ny_low_meas*adc_low_backup.numSamples;
% end
% ktraj_adc_epi_refonly = ktraj_adc_epi(:,1:number_of_elements_in_ref);

% ktraj_adc_epi_high_low = ktraj_adc_epi(:,1+number_of_elements_in_ref:end);
% ktraj_adc_epi_high_low = reshape(ktraj_adc_epi_high_low,3,[],Nslices,RSegment,nImagingFrame);
% ktraj_adc_epi_high = ktraj_adc_epi_high_low(:,1:Ny_meas*adc.numSamples,:,:,:);
% ktraj_adc_epi_high = reshape(ktraj_adc_epi_high,3,adc.numSamples,Ny_meas,Nslices,RSegment,nImagingFrame);
% 
% if is2ndEcho
%     ktraj_adc_epi_low = ktraj_adc_epi_high_low(:,1+Ny_meas*adc.numSamples:end,:,:,:);
%     ktraj_adc_epi_low = reshape(ktraj_adc_epi_low,3,adc_low_backup.numSamples,Ny_low_meas,Nslices,RSegment,nImagingFrame);
% end
% 
% % plot ky-t traj
% % frameInd = 2:nImagingFrame; % nFrames is not correct, since dummy does not have adc
% sliceInd = 1;
% ky = squeeze(   ktraj_adc_epi_high(2,end,:,sliceInd,:,end-nImagingFrame+1:end)   );
% ky = round(  ky/deltak  );
% figure;
% plot(ky(:),'-o');
% hold on;
% yline(0,'m--');
% ylim_values = ylim(); % Get y-axis limits
% y_top = ylim_values(2); % Set text at the top of the plot
% if RSegment > 1
%     for i = 1:RSegment
%         if i < RSegment
%             xline(0.5 + i * Ny_meas, 'k--'); % Draw segment boundaries
%         end
%         % Add text at the top center of each segment
%         text((i-0.5) * Ny_meas, y_top, sprintf('Segment %d', i), ...
%             'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
%             'FontSize', 12, 'Color', 'k', 'FontWeight', 'bold');
%     end
% end
% xlabel 'index'
% ylabel 'ky'
% title 'ky vs. time'
% yline(Nx/2,'b--')
% yline(-Nx/2,'b--')
% 
% 
% 
% diff_frame_ind = 1;
% fprintf('<strong>In diffusion frame #%d, measured ky lines in %d shots=%d, unique ky lines=%d</strong>\n',...
%     diff_frame_ind,size(ky,2),numel(ky(:,:,diff_frame_ind)),numel(unique(ky(:,:,diff_frame_ind))))
% 
% % show ky coverage
% figure; hold on;
% % Define colors for each echo
% colors = lines (debra_nechoes);
% for iter_echo = 1:debra_nechoes
%     ky_echo = ky(iter_echo:debra_nechoes:end,:,diff_frame_ind);
%     ky_echo = ky_echo(:);
%     scatter(iter_echo*ones(size(ky_echo)),ky_echo,  20, colors(iter_echo, :), 'filled', 'MarkerFaceAlpha',1);
%     
%     fprintf('<strong>In diffusion frame #%d, echo #%d, measured ky lines in %d shots=%d, unique ky lines=%d</strong>\n',...
%         diff_frame_ind,iter_echo,size(ky,2),numel(ky_echo),numel(unique(ky_echo)))    
% end
% xlabel('echo index');
% ylabel 'ky'
% legend(arrayfun(@(x) sprintf('Echo %d', x), 1:debra_nechoes, 'UniformOutput', false));
% grid on;
% hold off;
% title 'ky sampling pattern (sum over segments)'



end


%%
seq.write('epise_rs.seq');

% seq.install('siemens');

% seq.sound(); % simulate the seq's tone

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits
if isTestRun
    rep = seq.testReport;
    fprintf([rep{:}]); % as for January 2019 TR calculation fails for fat-sat
end

if ~isTestRun
    return
end

%%
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
if ~isGEscanner && max(tpns)>0.85
    warning('PNS=%.2f too high, the sequence may not run on the scanner',max(tpns))
end


% seq.install('siemens');
% seq.sound(); % simulate the seq's tone


% check forbidden frequencies
if ~isGEscanner && isTestRun
    tic
    fprintf('Checking frequencies... ');
    gradSpectrum(seq,sys,sys_type);
    toc
end



%%
function [seq,segmentID] = addacs(Nx, fov, thickness, sys, sys_lowPNS, adc, seq,segmentID,pislquant,Nslices)
rfSpoilingInc=117;              % RF spoiling increment

bandwidth = 250;        % Hz/pixel
flip_angle = 10;        % degrees


% when Tread is defined as 1/bandwidth, there is "an assertion failed"
% message when writing the sequence at the end (for bandwidth=300)
% maybe Tread becomes 3.33... ms, which causes some problem?
Tread_acs = 1 / bandwidth;  % readout duration in sec
% tmpDwell = adc.dwell;%ceil(5e-6/sys.adcRasterTime)*sys.adcRasterTime; % GE needs to be multiple of 2us
dwell_time = compute_dwell_time(adc.numSamples, sys.adcRasterTime, sys.blockDurationRaster);
%  Tread_acs = adc.numSamples * adc.dwell;
Tread_acs = adc.numSamples * dwell_time;
fprintf('ACS bandwidth is: %.2f Hz/pixel\n',1/Tread_acs);


Tpre = ceil(1e-3/sys.gradRasterTime)*sys.gradRasterTime;            % duration of Gz/Gy gradient blips
Ndummy = 50;            % dummy shots to enter steady state


%--------------------------------------------------------------------------
% calculate timing: ACS
%--------------------------------------------------------------------------

TE_acs = ceil(6e-3/sys.blockDurationRaster)*sys.blockDurationRaster;
TR_acs = ceil(12e-3/sys.blockDurationRaster)*sys.blockDurationRaster;

Ny_acs = 32;
Nz_acs = 32;

slicethickness = thickness;
[rf, gzRF] = mr.makeSincPulse(flip_angle*pi/180,'Duration',3e-3,...
    'SliceThickness',slicethickness,'apodization',0.42,'timeBwProduct',4,'system',sys);

rf_acs = rf;
gzRF_acs = gzRF;

% Define other gradients and ADC events
deltak = 1 ./ fov;

gx_acs = mr.makeTrapezoid('x', sys, 'FlatArea', Nx*deltak(1), 'FlatTime', Tread_acs);    % readout gradient

gxPre_acs = mr.makeTrapezoid('x', sys, 'Area', -gx_acs.area/2,'Duration',Tpre);        % Gx prewinder

gxSpoil_acs = mr.makeTrapezoid('x', sys, 'Area', gx_acs.area,'Duration',Tpre);         % Gx spoiler

gzReph_acs = mr.makeTrapezoid('z','Area',-gzRF.area/2,'Duration',mr.calcDuration(gxPre_acs),'system',sys_lowPNS);

% ACS adc object: match the number of samples in the epi readout
% adc = mr.makeAdc(numSamples, sys, 'Duration', Tread, 'Delay', blipDuration/2);  % + dwell/2? TODO
adc_acs = mr.makeAdc(adc.numSamples,sys, 'Duration', gx_acs.flatTime, 'Delay', gx_acs.riseTime);            % adc event


areaY = ((0:Ny_acs-1)-Ny_acs/2)*deltak(1);


delayTE_acs = ceil( (TE_acs - mr.calcDuration(rf) + mr.calcRfCenter(rf) + rf.delay - mr.calcDuration(gxPre_acs)  ...
    - mr.calcDuration(gx_acs)/2)/seq.gradRasterTime)*seq.gradRasterTime;

if delayTE_acs < 0
    assert('error: acs delay TE is negative')
else
    disp(['acs delay TE: ', num2str(1e3*delayTE_acs), ' ms'])
end


delayTR_acs = ceil((TR_acs - mr.calcDuration(rf) - mr.calcDuration(gxPre_acs) ...
    - mr.calcDuration(gx_acs) - mr.calcDuration(gxSpoil_acs) - delayTE_acs)/seq.gradRasterTime)*seq.gradRasterTime;

if delayTR_acs < 0
    assert('error: acs delay TR is negative')
else
    disp(['acs delay TR: ', num2str(delayTR_acs*1e3), ' ms'])
end


%--------------------------------------------------------------------------
% ACS data: Loop over phase encodes and define sequence blocks
%--------------------------------------------------------------------------
gyPre_acs = mr.makeTrapezoid('y',sys, 'Area', max(abs(areaY)), 'Duration', Tpre);
% gyReph_acs = mr.makeTrapezoid('y',sys, 'Area', -areaY(floor(ny/2)), 'Duration', Tpre);
peScales = areaY/gyPre_acs.area;


rf_phase = 0;
rf_inc = 0;

dummy_adc_acs = mr.makeDelay(mr.calcDuration(adc_acs));
adc_acs_backup = adc_acs;
segidbase = segmentID + 1;
segidmax = segidbase;
slice_indices = tdr_sliceorder(Nslices,1);

% step1, receive gain quantification, for GE only, and only for 1 slice
if pislquant > 0 % for GE only
    for s = slice_indices(ceil(end/2))
        rf_acs.freqOffset=gzRF_acs.amplitude*thickness*(s-1-(Nslices-1)/2);
        % for iY = (-Ndummy-pislquant+1):Ny_acs
        for iY = (-Ndummy-pislquant+1):0 % for receiver gain calib, do not need imaging data
            
            isDummyTR = iY <= -pislquant;
            isReceiveGainCalibrationTR = iY < 1 & iY > -pislquant;
            
            if isDummyTR
                adc_acs = dummy_adc_acs;
            else
                adc_acs = adc_acs_backup;
            end
            
            pesc = (iY>0) * peScales(max(iY,1));  % phase-encode gradient scaling
            pesc = pesc + (pesc == 0)*eps;        % non-zero scaling so that the trapezoid shape is preserved in the .seq file
            
            
            % excitation pulse
            rf_acs.phaseOffset = rf_phase/180*pi;
            adc_acs.phaseOffset = rf_phase/180*pi;
            
            segmentID = segidbase + isDummyTR + 2*isReceiveGainCalibrationTR;
            segidmax = max(segmentID,segidmax);
            seq.addBlock(rf_acs, gzRF_acs, mr.makeLabel('SET', 'TRID', segmentID));
            
            rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase+rf_inc, 360.0);
            
            % Encoding
            seq.addBlock(gxPre_acs, mr.scaleGrad(gyPre_acs,pesc),gzReph_acs);        % Gz, Gy blips, Gx pre-winder
            seq.addBlock(mr.makeDelay(delayTE_acs));    % delay until readout
            
            seq.addBlock(gx_acs, adc_acs);                       % readout
            
            seq.addBlock(mr.scaleGrad(gyPre_acs,-pesc), gxSpoil_acs);% -Gz, -Gy blips, Gx spoiler
            seq.addBlock(mr.makeDelay(delayTR_acs))     % wait until end of TR
        end
    end
    segidbase = segidmax + 1;
    segidmax = segidbase;
end

% step 2, real imaging
pislquant = 0; % we do not need it anymore
for s = slice_indices
    rf_acs.freqOffset=gzRF_acs.amplitude*thickness*(s-1-(Nslices-1)/2);
    for iY = (-Ndummy-pislquant+1):Ny_acs
        
        isDummyTR = iY <= -pislquant;
        isReceiveGainCalibrationTR = iY < 1 & iY > -pislquant;
        
        if isDummyTR
            adc_acs = dummy_adc_acs;
        else
            adc_acs = adc_acs_backup;
        end
        
        pesc = (iY>0) * peScales(max(iY,1));  % phase-encode gradient scaling
        pesc = pesc + (pesc == 0)*eps;        % non-zero scaling so that the trapezoid shape is preserved in the .seq file
        
        
        % excitation pulse
        rf_acs.phaseOffset = rf_phase/180*pi;
        adc_acs.phaseOffset = rf_phase/180*pi;
        
        segmentID = segidbase + isDummyTR + 2*isReceiveGainCalibrationTR;
        segidmax = max(segmentID,segidmax);
        seq.addBlock(rf_acs, gzRF_acs, mr.makeLabel('SET', 'TRID', segmentID));
        
        rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase = mod(rf_phase+rf_inc, 360.0);
        
        % Encoding
        seq.addBlock(gxPre_acs, mr.scaleGrad(gyPre_acs,pesc),gzReph_acs);        % Gz, Gy blips, Gx pre-winder
        seq.addBlock(mr.makeDelay(delayTE_acs));    % delay until readout
        
        seq.addBlock(gx_acs, adc_acs);                       % readout
        
        seq.addBlock(mr.scaleGrad(gyPre_acs,-pesc), gxSpoil_acs);% -Gz, -Gy blips, Gx spoiler
        seq.addBlock(mr.makeDelay(delayTR_acs))     % wait until end of TR
    end
end


% disp(['ACS scan time: ', num2str((Ndummy + Ny_acs * Nz_acs) * TR_acs), ' sec'])
fprintf('ACS scan time: %.1f s\n',sum(seq.blockDurations));

segmentID = segidmax + 1;
gz_dummy = mr.makeTrapezoid('z','system',sys_lowPNS,'amplitude',eps,'duration',2.5);
% for GE, each segment must have at least 2 blocks
seq.addBlock(gz_dummy,mr.makeLabel('SET', 'TRID', segmentID));
seq.addBlock(gz_dummy);
end


function dwell_time = compute_dwell_time(numSamples, adcraster, blkRaster)
% Computes a valid dwell_time based on the following conditions:
% 1. dwell_time must be an integer multiple of adcraster.
% 2. numSamples * dwell_time must be an integer multiple of blkRaster.
% 3. Bandwidth per point (BWPP = 1 / (numSamples * dwell_time)) must be within BWPP limits.

% Set BWPP limits
BWPP_lower_limit = 200;  % Hz/pixel
BWPP_upper_limit = 1500; % Hz/pixel

% Compute valid range for dwell_time
min_dwell_time = 1 / (numSamples * BWPP_upper_limit);
max_dwell_time = 1 / (numSamples * BWPP_lower_limit);

% Start with the largest possible k that satisfies the BWPP limits
max_k = floor(max_dwell_time / adcraster);
min_k = ceil(min_dwell_time / adcraster);

for k = max_k:-1:min_k
    candidate = k * adcraster; % Ensure dwell_time is a multiple of adcraster
    BWPP = 1 / (numSamples * candidate);
    
    ratio = (numSamples * candidate) / blkRaster;
    tolerance = 1e-2;
    isDivisible = abs(ratio - round(ratio)) < tolerance;
    if isDivisible && BWPP >= BWPP_lower_limit && BWPP <= BWPP_upper_limit
        dwell_time = candidate;
        return;
    end
end

% If no valid dwell_time is found, throw an error
error('No valid dwell_time found within the given constraints.');
end

%%
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