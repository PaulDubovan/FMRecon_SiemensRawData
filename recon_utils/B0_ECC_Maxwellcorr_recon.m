function [data,RXTime] = B0_ECC_Maxwellcorr_recon(data,twix,simfolder,start_index,dt,ispulseq,isEPI,param)    
% Estimation of zeroth order eddy current phase and gradient waveforms using
% Siemens IDEA simulation files, as well as zeroth order Maxwell terms
% based on the nominal gradient waveforms
% Zeroth order phase modulation of the raw data
%
% Inputs:
%   data:           Raw data array
%   twix:           twix data
%   simfolder:      Path to IDEA simulation data
%   start_index:    Starting index 
%   dt:             Imaging dwell time
%   ispulseq:       Pulseq flag 
%   isEPI:          EPI flag 
%   param:          Supplied parameters  
%
% Outputs:
%   data:           Phase-modulated raw data array
%   RXTime:         Sampling times
%
% Copyright (c) 2026 Paul Dubovan and Xingwang Yong

    sys_type = param.sys_type;
    Nfre = param.Nfre;
    Nphs = param.Nphs;
    Nshots = param.Nshots;
    Nsli = param.Nsli;
    Nechoes = param.Nechoes;
    Nvol = param.Nvol; 
    Ncha = param.Ncha;
    MaxwellCompfromSim = param.MaxwellCompfromSim;
    if isfield(param,'isRefscan')
        isRefscan = param.isRefscan;
    end    
    if isEPI & ispulseq
        Npershot = param.Npershot;
    end
    %% Extract scanner eddy current parameters (amplitude & time constants)
    CoeffsX.Amp = cell2mat(twix.hdr.MeasYaps.sGRADSPEC.asGPAData{1,2}.sB0CompensationX.aflAmplitude);
    CoeffsX.Tau = cell2mat(twix.hdr.MeasYaps.sGRADSPEC.asGPAData{1,2}.sB0CompensationX.aflTimeConstant);
    
    CoeffsY.Amp = cell2mat(twix.hdr.MeasYaps.sGRADSPEC.asGPAData{1,2}.sB0CompensationY.aflAmplitude);
    CoeffsY.Tau = cell2mat(twix.hdr.MeasYaps.sGRADSPEC.asGPAData{1,2}.sB0CompensationY.aflTimeConstant);
    
    CoeffsZ.Amp = cell2mat(twix.hdr.MeasYaps.sGRADSPEC.asGPAData{1,2}.sB0CompensationZ.aflAmplitude);
    CoeffsZ.Tau = cell2mat(twix.hdr.MeasYaps.sGRADSPEC.asGPAData{1,2}.sB0CompensationZ.aflTimeConstant);
     
    %% Calculate the eddy current phase based on xml file
    doInterpToRX = true;
    [RXSamples, RXTime, TXTime, TrigTime, ECCPhase] = SeqSim_mex(simfolder,'eddyPhase',doInterpToRX,CoeffsX,CoeffsY,CoeffsZ);
    %% Calculate zeroth-order Maxwell term
    % Get data from dsv files
    doInterpToRX = true;
    [~, ~, ~, ~, gx, gy, gz] = SeqSim_mex(simfolder,'gradient',doInterpToRX);

    %%
    % Identify all valid ADC timepoints using the scan counter (there should only
    % be jumps for the ep2d and ep3d sequences)
    if twix.image.scancounter(1) == 1
        if exist('isRefscan') && isRefscan && ispulseq
            idx = RXSamples-RXSamples(1)+1;            
        else
            idx = RXSamples(twix.image.scancounter)-RXSamples(1)+1;
        end    
    else
        idx = RXSamples(twix.image.scancounter-1)+1;
    end
    if ~exist('isRefscan') || ~isRefscan || ~ispulseq
        idx = idx(start_index:end);
    end

    if 0 % old slow way
        sz = length(idx);
        for i = 1:sz
            idx_final_old(i,:) = idx(i):idx(i)+(Nfre-1);
        end
    else
        offset = 0:(Nfre-1);
        idx_final_new = idx(:) + int32(offset);
        idx_final = idx_final_new;
    end
    idx_final = idx_final';
    idx_final = idx_final(:);
    
    gx = gx(idx_final);
    gy = gy(idx_final);
    gz = gz(idx_final);
    
    ECCPhase = ECCPhase(idx_final);
    ECCPhase = ECCPhase';

    data = data.*exp(-1i*ECCPhase);
        
    %% Calculation of phase contributions from zeroth order concomitant fields: for asymmetric gradient coils
    %% Scanner SPECIFIC Maxwell compensation
    if strcmp(sys_type,'Connectome2') && MaxwellCompfromSim
        coilParams.alph = 0.5;
        coilParams.z0x = -0.1;
        coilParams.z0y = -0.09;
        coilParams.x0 = 0;
        coilParams.y0 = 0;
  
        gx = gx/1000;
        gy = gy/1000;
        gz = gz/1000;
        B0 = twix.hdr.Dicom.flMagneticFieldStrength; %T
        gam = twix.hdr.Meas.alLarmorConstant(1); %Hz/T
        dt = dt*1e-6;
        %Maxwell_zero = gx.^2 * PsiX^2 /(2*B0) + gy.^2 * PsiY^2 /(2*B0) ; 
        maxwell_zero = computeMaxwellTerms(0,coilParams,gx,gy,gz,B0);
        maxwell_zero = maxwell_zero.zeroth';
    
        if ~isEPI && ispulseq
            maxwell_zero = squeeze(reshape(maxwell_zero, [Nfre*Nphs, Nechoes, Nsli, Nshots, Nvol]));
        elseif isEPI && ispulseq
            maxwell_zero = squeeze(reshape(maxwell_zero, [Nfre*Npershot, Nechoes, Nsli, Nshots, Nvol]));
        else    
            maxwell_zero = squeeze(reshape(maxwell_zero, [Nfre*Nphs, Nechoes, Nsli, Nvol]));
        end  
    
        phs_maxwell = cumsum((2*pi*gam*dt)*maxwell_zero,1);
    
        if ~isEPI && ispulseq
            phs_maxwell = reshape(phs_maxwell, [Nfre*Nphs*Nechoes*Nsli*Nshots*Nvol, 1]);
        else  
            phs_maxwell = reshape(phs_maxwell, [Nfre*Nphs*Nechoes*Nsli*Nvol, 1]);
        end

        data = data.*exp(1i*phs_maxwell);  
    end         

    %%
    RXTime = 1e6*RXTime(idx_final);
end
