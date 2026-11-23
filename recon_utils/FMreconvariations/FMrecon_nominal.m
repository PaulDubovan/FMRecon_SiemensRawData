function [im,im_denoised] = FMrecon_nominal(filename_img,filename_gre,session_filepath,simfolder,scanId,opt,seq_file_path,spm_B0mapPath)
% Handling of raw Siemens twix data and Skope data to be used for image
% reconstruction using an interative expanded encoding model
%
% Inputs:
%   filename_img:       Path to raw imaging data
%   filename_gre:       Path to GRE data
%   session_filepath:   Path to Skope NYOX .h5 file
%   simfolder:          Path to IDEA sequence simulaton files (
%   scanId:             scan ID of Skope acquisition
%   opt:                Preset reconstruction options 
%   seq_file_path:      Path to corresponding Pulseq .seq file 
%   spm_B0mapPath:      Path to B0 map generated using spm toolbox  
%
%
% Outputs:
%   im:                 Array of reconstructed images [x y z Nvol]
%   im_denoised:        Denoised images using tMPPCA 
%
% Copyright (c) 2026 Paul Dubovan and Xingwang Yong

% Default inputs that can be initialized before calling the recon function
    input_vol = opt.input_vol;
    input_sli = opt.input_sli; 

    if ~isfield(opt, 'sys_type') || isempty(opt.sys_type)
        sys_type = 'Connectome2';
    else
        sys_type = opt.sys_type;
    end    

    if ~isfield(opt, 'ispulseq') || isempty(opt.ispulseq)
        ispulseq = 0;
    else
        ispulseq = opt.ispulseq;
    end  

    if ~isfield(opt, 'doparallelcomp') || isempty(opt.doparallelcomp)
        doparallelcomp = 0;
    else
        doparallelcomp = opt.doparallelcomp;
    end  
    
    if doparallelcomp
        assert(isscalar(input_vol),'Expect scalar');
        assert(isscalar(input_sli),'Expect scalar');
    end    
    
    if ~isfield(opt, 'spiraloutin_fact') || isempty(opt.spiraloutin_fact)
        spiraloutin_fact = 1;
    else
        spiraloutin_fact = opt.spiraloutin_fact;
    end 

    if ~isfield(opt, 'useGPU') || isempty(opt.useGPU)
        useGPU = 1;
    else
        useGPU = opt.useGPU;
    end    

    if ~isfield(opt, 'doConventionalRecon') || isempty(opt.doConventionalRecon)
        doConventionalRecon = 1;
    else
        doConventionalRecon = opt.doConventionalRecon;
    end   

    if ~isfield(opt, 'useSingle') || isempty(opt.useSingle)
        useSingle = 0;
    else
        useSingle = opt.useSingle;
    end    

    if ~isfield(opt, 'ishighres') || isempty(opt.ishighres)
        ishighres = 0;
    else
        ishighres = opt.ishighres;
    end 

    if ~isfield(opt, 'docoilcompress') || isempty(opt.docoilcompress)
        docoilcompress = 1;
    else
        docoilcompress = opt.docoilcompress;
    end 

    if ~isfield(opt, 'num_coil') || isempty(opt.num_coil)
        if docoilcompress
            num_coil = 30;
        end
    else
        num_coil = opt.num_coil;
    end   

    if ~isfield(opt, 'useskopefit') || isempty(opt.useskopefit)
        useskopefit = 0;
    else
        useskopefit = opt.useskopefit;
    end 

    if ~isfield(opt, 'dospiraldelay') || isempty(opt.dospiraldelay)
        dospiraldelay = 0;
    else
        dospiraldelay = opt.dospiraldelay;
    end 

    if ~isfield(opt,'retro_undersample_factor') || isempty(opt.retro_undersample_factor)
        retro_undersample_factor = 1;
    else
        retro_undersample_factor = opt.retro_undersample_factor;
    end

    if ~isfield(opt,'doPhsConstraint') || isempty(opt.doPhsConstraint)
        doPhsConstraint = 0;
    else
        doPhsConstraint = opt.doPhsConstraint;
    end

    if ~isfield(opt,'dotrajextrap') || isempty(opt.dotrajextrap)
        dotrajextrap = 0;
    else
        dotrajextrap = opt.dotrajextrap;
    end

    if ~isfield(opt,'dotensordenoise') || isempty(opt.dotensordenoise)
        dotensordenoise = 0;
    else
        dotensordenoise = opt.dotensordenoise;

        if dotensordenoise
            if ~(length(input_sli) > 1) && ~(length(input_vol) > 1)
                dotensordenoise = 0;
                warning('More slices and/or volumes recommended for denoising. Turning off denoising');
            end
        end    
    end

    if dotensordenoise
        if ~isfield(opt,'denoise_patchsize') || isempty(opt.denoise_patchsize)                    
            if length(input_sli) > 5
                denoise_patchsize = [5 5 5];
            elseif length(input_sli) > 3
                denoise_patchsize = [5 5 3];
            else
                denoise_patchsize = [5 5];
            end    
        else
            denoise_patchsize = opt.denoise_patchsize;
        end
    end    

    if ~isfield(opt,'doGNLcorr') || isempty(opt.doGNLcorr)
        doGNLcorr = 0;
    else
        doGNLcorr = opt.doGNLcorr;
    end

    if ~isfield(opt,'encoding_approach') || isempty(opt.encoding_approach)
        encoding_approach = 1;
    else
        encoding_approach = opt.encoding_approach;
    end

    if doGNLcorr
        % Interpolated approach not supported for not constant z
        % coordinates (i.e. not purely axial slice)
        encoding_approach = 0;
    end    

    if ~isfield(opt,'orient') || isempty(opt.orient)
        orient = 1;
    else
        orient = opt.orient;
    end   

    if ~isfield(opt,'dob0reg') || isempty(opt.dob0reg)
        dob0reg = 1;
    else
        dob0reg = opt.dob0reg;
    end

    if ~isfield(opt,'b0_indices') || isempty(opt.b0_indices)
        b0_indices = 0;
    else
        b0_indices = opt.b0_indices;
    end

    % change this based on necessary # iterations of S-loraks
    if ~isfield(opt,'max_iter') || isempty(opt.max_iter)
        max_iter = 10;
    else
        max_iter = opt.max_iter;
    end

    if ~isfield(opt,'b0_rank') || isempty(opt.b0_rank)
        b0_rank = 30;
    else
        b0_rank = opt.b0_rank;
    end

    if ~isfield(opt,'b0_lambda') || isempty(opt.b0_lambda)
        b0_lambda = 0.01;
    else
        b0_lambda = opt.b0_lambda;
    end

    if ~isfield(opt,'dwi_rank') || isempty(opt.dwi_rank)
        dwi_rank = 40;
    else
        dwi_rank = opt.dwi_rank;
    end

    if ~isfield(opt,'dwi_lambda') || isempty(opt.dwi_lambda)
        dwi_lambda = 0.005;
    else
        dwi_lambda = opt.dwi_lambda;
    end

    if ~isfield(opt,'doL2reg') || isempty(opt.doL2reg)
        doL2reg = 0;
    else
        doL2reg = opt.doL2reg;
    end

    if ~isfield(opt,'L2gamma') || isempty(opt.L2gamma)
        L2gamma = 0.1;
    else
        L2gamma = opt.L2gamma;
    end

    if ~isfield(opt,'doCSreg') || isempty(opt.doCSreg)
        doCSreg = 0;
    else
        doCSreg = opt.doCSreg;
    end

    if ~isfield(opt,'lamVal_CS') || isempty(opt.lamVal_CS)
        lamVal_CS = [];
    else
        lamVal_CS = opt.lamVal_CS;
    end

    if ~isfield(opt,'MaxwellCompfromSim') || isempty(opt.MaxwellCompfromSim)
        MaxwellCompfromSim = 0;
    else
        MaxwellCompfromSim = opt.MaxwellCompfromSim;
    end

    if ~isfield(opt,'doSMSwithkz') || isempty(opt.doSMSwithkz)
        doSMSwithkz = 1;
    else
        doSMSwithkz = opt.doSMSwithkz;
    end

    if ~isfield(opt,'docircshift') || isempty(opt.docircshift)
        docircshift = 1;
    else
        docircshift = opt.docircshift;
    end

    if ~isfield(opt,'doGhostCor') || isempty(opt.doGhostCor)
        doGhostCor = 0;
    else
        doGhostCor = opt.doGhostCor;
    end

    if ~isfield(opt,'doGibbsCor') || isempty(opt.doGibbsCor)
        doGibbsCor = 1;
    else
        doGibbsCor = opt.doGibbsCor;
    end

    if ~isfield(opt,'dophaseshift') || isempty(opt.dophaseshift)
        opt.dophaseshift = 0;
    end

    if ~isfield(opt,'pulseq_gre') || isempty(opt.pulseq_gre)
        pulseq_gre = 0;
    else
        pulseq_gre = opt.pulseq_gre;
    end   

    if ~isfield(opt,'doaxisrotate') || isempty(opt.doaxisrotate)
        doaxisrotate = 0;
    else
        doaxisrotate = opt.doaxisrotate;
    end  

    if exist('spm_B0mapPath', 'var') && exist(spm_B0mapPath, 'file')
        if isfield(opt, 'use_spm_B0map') && ~isempty(opt.use_spm_B0map)  
            if opt.use_spm_B0map 
                use_spm_B0map = 1;
            else
                use_spm_B0map = 0;
            end
        else
            use_spm_B0map = 0;
        end    
    else
        use_spm_B0map = 0;
    end    

    [folderPath, fn, ~] = fileparts(filename_img);
    
    % Define paths to Cgre.mat and B0map.mat
    twixPath = fullfile(folderPath, ['twix_____',fn,'.mat']);
    twix_grePath = fullfile(folderPath, 'twix_gre.mat');
    
    % Check if both files exist. If they don't, create them
    if exist(twixPath, 'file') && exist(twix_grePath, 'file')
        % Load the files
        load(twixPath, 'twix');    
        load(twix_grePath, 'twix_gre');   
        disp('Loading twix');
    else
        twix = mapVBVD(filename_img); 
        %twix = twix{2};
        twix_gre = mapVBVD(filename_gre); 
        %twix_gre = twix_gre{2};
        % Save the variables
        save(twixPath, 'twix');
        save(twix_grePath, 'twix_gre');
        disp('Saving twix');
    end  

    % Obtain noise matrix from noise-only acquisition
    noise_matrix = get_noise_matrix(twix{1});

    twix = twix{2};
    twix_gre = twix_gre{2};
    
    % If pulseq sequence, load necessary seq information
    if ispulseq
        seq = mr.Sequence(mr.opts('B0',2.89));
        seq.read(seq_file_path,'detectRFuse');

        % Trajectory ADC sampling times 
        [~, t_adc] = seq.calculateKspacePP();
        t_adc = 1e6*t_adc;

        twix_sig = twix.hdr.Dicom.tSequenceVariant;
        seq_sig = seq.signatureValue;
        
        % Check if twix and pulseq file version match
        if ~isequal(twix_sig,seq_sig) && ~contains(twix_sig,seq_sig) % ideally these two should be equal, but somehow sometimes the former has some extra trailing chars
            error('Signature mismatch')
        end
    end   

    % 3D acquisition identifier
    if twix.hdr.MeasYaps.sKSpace.ucDimension > 2 && ~ispulseq
        is3D = 1;
    else 
        is3D = 0;
    end   

    % Determine if EPI sequence
    if any(twix.image.IsReflected) == 1
        isEPI = 1;
    elseif ispulseq && ~isempty(seq.getDefinition('EchoSpacing'))
        isEPI = 1;
    else
        isEPI = 0;
    end

    % In-plane and slice acceleration
    R_fact = twix.hdr.MeasYaps.sPat.lAccelFactPE;
    SMS_fact = twix.hdr.MeasYaps.sSliceAcceleration.lMultiBandFactor;
    FOVshift = twix.hdr.MeasYaps.sSliceAcceleration.lFOVShiftFactor;
    if ispulseq
        if ~isempty(seq.getDefinition('SMS_fact'))
            SMS_fact = seq.getDefinition('SMS_fact');
        end    
        if ~isempty(seq.getDefinition('caipi_factor'))
            FOVshift = seq.getDefinition('caipi_factor');
        end  
    end    

   % Load acquisition parameters, using seq or twix file
    if ispulseq
        if ~isEPI
            Nfre = twix.image.NCol;
            Nphs = twix.image.NSet;
            Nshots = seq.getDefinition('NumShots');
            Nsli = seq.getDefinition('NumSlices');
            N = seq.getDefinition('MatrixSize');
            if ~isempty(seq.getDefinition('numSamplesPerShot'))
                Nenc = seq.getDefinition('numSamplesPerShot');
            elseif ~isempty(seq.getDefinition('numSamplesPerArm'))  
                Nenc = seq.getDefinition('numSamplesPerArm');
            end    
            Nechoes = seq.getDefinition('Nechoes');
            Nvol = seq.getDefinition('num_bvolumes');
            Ncha = twix.image.NCha;
            TEshift = 1e6*seq.getDefinition('base_delay_between_shots');
        else
            Nshots = seq.getDefinition('RSegment');
            Nsli = seq.getDefinition('Nslices');
            N = seq.getDefinition('Nx');
            Nechoes = 1;
            Nvol = length(seq.getDefinition('bval'));
            Nfre = twix.image.NCol;
            Npershot = seq.getDefinition('Ny_meas');
            Nphs = Npershot*Nshots;
            Ncha = twix.image.NCha;
            Nenc = Nfre*Nphs;
            TEshift = 1e6*seq.getDefinition('TEShift');
            isRefscan = seq.getDefinition('isRefscan');
            nImagingFrame = seq.getDefinition('nImagingFrame');
        end    
        if Nenc ~= Nfre*Nphs
            error('samples mismatch');
        end    
    else
        Nfre = twix.image.NCol;
        Nphs = twix.image.NLin/R_fact;
        if isEPI && (R_fact > 1)
            if mod(Nphs,2) == 0
                Nphs = Nphs - 1;  
            elseif mod(Nphs,1) ~= 0
                Nphs = ceil(Nphs);
            end
        end  
        if is3D
            Nsli = twix.image.NPar;
        else
            Nsli = twix.image.NSli/SMS_fact;
        end    
        Nvol = twix.image.NRep;
        Ncha = twix.image.NCha;
        Nenc = Nfre*Nphs;  

        if is3D
            Nshots = twix.hdr.MeasYaps.sFastImaging.lSegments;
        else
            Nshots = twix.hdr.MeasYaps.sFastImaging.lShots;
        end    
        Nechoes = 1;
    end   

    if ~docoilcompress
        num_coil = Ncha;
    end    

    param.sys_type = sys_type;
    param.Nfre = Nfre;
    param.Nphs = Nphs;
    param.Nshots = Nshots;
    param.Nsli = Nsli;
    param.Nechoes = Nechoes;
    param.Nvol = Nvol; 
    param.Ncha = Ncha;
    param.MaxwellCompfromSim = MaxwellCompfromSim;
    if isEPI & ispulseq
        param.Npershot = Npershot;
    end    
    
    % Assuming filename_gre is the full path to your GRE file
    [folderPath, ~, ~] = fileparts(filename_gre);
    
    % Define paths to Cgre.mat and B0map.mat
    CgrePath = fullfile(folderPath, 'Cgre.mat');
    B0mapPath = fullfile(folderPath, 'B0map.mat');
    
    % Check if both files exist
    if exist(CgrePath, 'file') && exist(B0mapPath, 'file')
        % Load the files
        load(CgrePath, 'Cgre');    
        load(B0mapPath, 'B0map');   
        disp('Loading sensitivity profiles and field maps');
    else
        if ispulseq
            [Cgre,B0map] = rcvrsplusb0map(twix,twix_gre,noise_matrix,ispulseq,pulseq_gre,Nsli*SMS_fact,N);
        else
            [Cgre,B0map] = rcvrsplusb0map(twix,twix_gre,noise_matrix,ispulseq,pulseq_gre,Nsli*SMS_fact);
        end    
        % Save the variables
        save(CgrePath, 'Cgre','-v7.3');
        save(B0mapPath, 'B0map');
        disp('Saving sensitivity profiles and field maps');
    end  

    if use_spm_B0map
        tmpStruct = load(spm_B0mapPath);
        B0map = tmpStruct.B0map_spm;
        B0map = B0map(:,:,tdr_sliceorder(Nsli*SMS_fact,1));
        disp('Using spm field map')

        switch orient
            case 2
                B0map = flip(B0map,2);
                B0map = imrotate(B0map,-90);
            case 3
                B0map = permute(B0map, [2 1 3]);
        end
    end    

    if use_spm_B0map
        B0map = 2*pi*B0map; % Hz to rad/s
    else
        B0map = 1000*B0map; % rad/ms to rad/s
    end  

    % Reversing slices, since it appears pulseq acquires the imaging data
    % in reverse order
    if ispulseq
        Cgre = Cgre(:,:,:,end:-1:1);
        B0map = B0map(:,:,end:-1:1);
    end

    % Load data in the order it was acquired in
    data = twix.image.unsorted();
    % Send channels to the end
    data = permute(data, [1 3 2]);

    % Flip even EPI lines for twix
    if isEPI && ~ispulseq
        idx = find(twix.image.IsReflected);
        data(:,idx,:) = flip(data(:,idx,:),1);         
    end

    % Load dwell time (us)
    if ispulseq
        adc_array = seq.adcLibrary.data(1).array;
        adc_array = adc_array(1, 2);
        dt = 1e6*adc_array;
    else
        dt = twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1}*1e-3; 
    end

    start_index = 1;
    num_acq = prod([Nphs,Nsli,Nshots,Nvol,Nechoes]);
    if isEPI
        num_acq = prod([Nphs,Nsli,Nvol,Nechoes]); %Nphs includes the multiple segments
    end
    if num_acq < twix.image.NAcq
        if ~exist('isRefscan')
            isRefscan = 1;
        end
        param.isRefscan = isRefscan;        
        % Remove reference scans
        start_index = twix.image.NAcq - num_acq + 1;
        data_ref = data(:,1:start_index-1,:); 
        data = data(:,start_index:end,:);
    end 
    
    % Combine data dimensions except for channels
    if ~isEPI && ispulseq
        data = reshape(data, [Nfre*Nphs*Nechoes*Nsli*Nshots*Nvol, Ncha]);
    else  
        data = reshape(data, [Nfre*Nphs*Nechoes*Nsli*Nvol, Ncha]);
    end    

    % Prewhiten the data
    data = my_prewhiten(data,noise_matrix);

    % Undo zeroth order eddy current and Maxwell compensation 
    % Determine sampling times needed for skope interpolation
    % [data,RXTime] = B0_ECC_Maxwellcorr_recon(data,twix,simfolder,start_index,dt,ispulseq,isEPI,param);
    % %load('data.mat');
    % %load('RXTime.mat');
    % dt_RXTime = RXTime(2)-RXTime(1);
    % if round(dt,1) ~= round(dt_RXTime,1)
    %     error('dwelltime mismatch')
    % elseif round(dt,4) ~= round(dt_RXTime,4)
    %     dt = dt_RXTime;
    % end    

    if ~isEPI && ispulseq
        data = squeeze(reshape(data, [Nfre*Nphs, Nechoes, Nsli, Nshots, Nvol, Ncha]));
    elseif isEPI && ispulseq
        data = squeeze(reshape(data, [Nfre*Npershot, Nechoes, Nsli, Nshots, Nvol, Ncha]));
    else    
        data = squeeze(reshape(data, [Nfre*Nphs, Nechoes, Nsli, Nvol, Ncha]));
    end    

    if Nechoes > 1
        data = squeeze(data(:,1,:,:,:,:));
    end   
    if spiraloutin_fact > 1
        Nenc = Nenc/spiraloutin_fact;
        data = data(1:Nenc,:,:,:,:);
    end
    
    % Ghost correction for EPI 
    if exist('isRefscan') && isRefscan && doGhostCor && ispulseq
        data_ref = reshape(data_ref, [Nfre Npershot, Nsli, 1, Ncha]);
        data = reshape(data, [Nfre Npershot, Nsli, Nshots, Nvol, Ncha]);
        AccY = 1;
    
        % Flip even EPI lines for twix
        idx = 2:2:Npershot;
        data(:,idx,:,:,:,:) = flip(data(:,idx,:,:,:,:),1);
        data_ref(:,idx,:,:,:) = flip(data_ref(:,idx,:,:,:),1); 

        for iter_diff_frame = 1:nImagingFrame
            sqz_dim4 = @(A) reshape(A, [size(A,1), size(A,2), size(A,3), size(A,5:numel(size(A)))]);
            % PAT_nav = zeros([Nx, Ny_meas,  1, num_virtual_chan, 2]);
            % PAT_nav(:,1:2:end,1,:,1) = sqz_dim4( data_ref_regrid(:,1:2:end,iter_slice,1,:) );
            % PAT_nav(:,2:2:end,1,:,2) = sqz_dim4( data_ref_regrid(:,2:2:end,iter_slice,1,:) );
            PAT_nav = zeros([Nfre, 3, Nsli, Ncha, 2]);
            PAT_nav(:,1:2:end,:,:,1) = sqz_dim4( data_ref(:,1:2:3,:,1,:) );
            PAT_nav(:,2:2:end,:,:,2) = sqz_dim4( data_ref(:,2:2:3,:,1,:) );                    
            PAT_nav = permute(PAT_nav, [1,2,4,6,7,8,9,5,10,3]);
            for iter_seg = 1:Nshots
                sqz_dim4_dim5 = @(A) reshape(A, [size(A,1), size(A,2), size(A,3), size(A,6:numel(size(A)))]);
                PAT_ref = zeros([Nfre, Npershot,  Nsli, Ncha, 2]);
                PAT_ref(:,1:2:end,:,:,1) = sqz_dim4_dim5(  data(:,1:2:end,:,iter_seg,iter_diff_frame,:)  );
                PAT_ref(:,2:2:end,:,:,2) = sqz_dim4_dim5(  data(:,2:2:end,:,iter_seg,iter_diff_frame,:)  );
                PAT_ref = permute(PAT_ref, [1,2,4,6,7,8,9,5,10,3]);

                % ghost correction function without regridding:
                tmp = ghost_correct_pat_ref_v1_STD_bb_v1(PAT_ref, PAT_nav, AccY);
                tmp = reshape(tmp, size(tmp,1),size(tmp,2),Ncha,size(tmp,10));
                tmp = permute(tmp,[1 2 4 3]);
                navCor_kspace(:,:,:,iter_seg,iter_diff_frame,:) = tmp;
            end
        end
        navCor_kspace(:,idx,:,:,:,:) = flip(navCor_kspace(:,idx,:,:,:,:),1);         
        data = squeeze(reshape(navCor_kspace, [Nfre*Npershot, Nechoes, Nsli, Nshots, Nvol, Ncha]));
    elseif exist('isRefscan') && isRefscan && doGhostCor && ~ispulseq
        data_ref = twix.phasecor.unsorted();
        data_ref = permute(data_ref, [1 3 2]);
        % The reference data lines should already be reflected to be same orientation 
        % idx_ref = find(twix.phasecor.IsReflected);
        % data_ref(:,idx_ref,:) = flip(data_ref(:,idx_ref,:),1);
        Nphs_ref = size(data_ref,2)/Nsli;
        data_ref = reshape(data_ref, [Nfre, Nphs_ref, Nsli, Ncha]);
        data = reshape(data, [Nfre Nphs, Nsli, Nvol, Ncha]);
        AccY = 1;
    
        % Flip even EPI lines for twix
        idx = 2:2:Nphs;
        data(:,idx,:,:,:) = flip(data(:,idx,:,:,:),1);

        for iter_diff_frame = 1:Nvol
            % PAT_nav = zeros([Nx, Ny_meas,  1, num_virtual_chan, 2]);
            % PAT_nav(:,1:2:end,1,:,1) = sqz_dim4( data_ref_regrid(:,1:2:end,iter_slice,1,:) );
            % PAT_nav(:,2:2:end,1,:,2) = sqz_dim4( data_ref_regrid(:,2:2:end,iter_slice,1,:) );
            PAT_nav = zeros([Nfre, 3, Nsli, Ncha, 2]);
            PAT_nav(:,1:2:end,:,:,1) = data_ref(:,1:2:3,:,:);
            PAT_nav(:,2:2:end,:,:,2) = data_ref(:,2:2:3,:,:);                    
            PAT_nav = permute(PAT_nav, [1,2,4,6,7,8,9,5,10,3]);
            PAT_ref = zeros([Nfre, Nphs,  Nsli, Ncha, 2]);
            PAT_ref(:,1:2:end,:,:,1) = data(:,1:2:end,:,iter_diff_frame,:);
            PAT_ref(:,2:2:end,:,:,2) = data(:,2:2:end,:,iter_diff_frame,:);
            PAT_ref = permute(PAT_ref, [1,2,4,6,7,8,9,5,10,3]);

            % ghost correction function without regridding:
            tmp = ghost_correct_pat_ref_v1_STD_bb_v1(PAT_ref, PAT_nav, AccY);
            tmp = reshape(tmp, size(tmp,1),size(tmp,2),Ncha,size(tmp,10));
            tmp = permute(tmp,[1 2 4 3]);
            navCor_kspace(:,:,:,iter_diff_frame,:) = tmp;
        end
        navCor_kspace(:,idx,:,:,:) = flip(navCor_kspace(:,idx,:,:,:),1);         
        data = squeeze(reshape(navCor_kspace, [Nfre*Nphs, Nsli, Nvol, Ncha]));
    end

    % The ep3d sequence acquires subsequent shots for each partition,
    % whereas multi-shot in pulseq (traditional 2D) acquires all slices
    % followed by changing the shot
    if Nshots > 1 && ispulseq && ~isEPI
        data = permute(data, [1 3 2 4 5]);
        data = reshape(data, [Nenc*Nshots Nsli Nvol Ncha]);
    elseif Nshots > 1 && ispulseq && isEPI
        data = permute(data, [1 3 2 4 5]);
        data = reshape(data, [Nenc Nsli Nvol Ncha]);
    end   
    
    % Use only sampling times from first shot
    % if ~isEPI && ispulseq
    %     RXTime = RXTime(1:Nenc);
    % else
    %     RXTime = RXTime(1:(Nenc/Nshots));
    % end    
    % datatime = RXTime - RXTime(1) + dt;  
    
    if ~isEPI && ispulseq
        datatime = (1:Nenc)*dt;
    else
        datatime = (1:(Nenc/Nshots))*dt;
    end  
    datatime = datatime';

    % Verify sampling times 
    if isEPI && ispulseq
        t_adc = t_adc(1:(Nenc/Nshots))';
        datatime_pulseq_epi = t_adc - t_adc(1) + dt;
        tolerance = 1e-2;
        if any(abs(datatime - datatime_pulseq_epi) > tolerance)
            warning('Pulseq and simulation sampling times are not the same');
        end  
    end

    datatime_interp = datatime;
    datatime = 1e-6*datatime;

    % Define imaging coordinates
    if ispulseq
        rec_FOVx = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
        rec_FOVy = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV;
        rec_FOVz = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;
        
        rec_Nx = N;
        rec_Ny = N;
        rec_Nz = Nsli;
        
        deltax = rec_FOVx / rec_Nx;
        deltay = rec_FOVy / rec_Ny;
        deltaz = rec_FOVz / rec_Nz;

        x = (-rec_FOVx/2) + (0:rec_Nx-1) * deltax;
        y = (-rec_FOVy/2) + (0:rec_Ny-1) * deltay;
        z = ((0:rec_Nz-1) - (rec_Nz-1)/2) * (rec_FOVz/rec_Nz);
    elseif ~is3D
        rec_FOVx = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
        rec_FOVy = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV;
        
        rec_Nx = twix.hdr.MeasYaps.sKSpace.lBaseResolution;
        rec_Ny = twix.hdr.MeasYaps.sKSpace.lPhaseEncodingLines;
        rec_Nz = twix.hdr.MeasYaps.sSliceArray.lSize;
        
        deltax = rec_FOVx / rec_Nx;
        deltay = rec_FOVy / rec_Ny;
        deltaz = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;
    
        rec_FOVz = rec_Nz*deltaz;

        x = (-rec_FOVx/2) + (0:rec_Nx-1) * deltax;
        y = (-rec_FOVy/2) + (0:rec_Ny-1) * deltay;
        z = ((0:rec_Nz-1) - (rec_Nz-1)/2) * (rec_FOVz/rec_Nz);
    else
        rec_FOVx = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
        rec_FOVy = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV;
        rec_FOVz = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;
        
        rec_Nx = twix.hdr.MeasYaps.sKSpace.lBaseResolution;
        rec_Ny = twix.hdr.MeasYaps.sKSpace.lPhaseEncodingLines;
        rec_Nz = twix.hdr.MeasYaps.sKSpace.lPartitions;
        
        deltax = rec_FOVx / rec_Nx;
        deltay = rec_FOVy / rec_Ny;
        deltaz = rec_FOVz / rec_Nz;
    
        x = (-rec_FOVx/2) + (0:rec_Nx-1) * deltax;
        y = (-rec_FOVy/2) + (0:rec_Ny-1) * deltay;
        z = ((0:rec_Nz-1) - (rec_Nz-1)/2) * (rec_FOVz/rec_Nz);
    end  

    % Proper orientations to match GRE scan
    if ~ispulseq
        y = -flip(y);
    end    

    if all(diff(twix.image.slicePos(3,:)) == 0)
        z_shift = twix.image.slicePos(3,1);
        z = z + z_shift;
    elseif ~isempty(twix.hdr.Dicom.lScanRegionPositionIso)
        z_shift = twix.hdr.Dicom.lScanRegionPositionIso;
        z = z + z_shift;     
    end

    if ispulseq
        if SMS_fact == 1
            % Flipping z since pulseq imaging data acquired head to foot
            z = flip(z);  
            %Then order z to match order in which slices were acquired
            z = z(:,tdr_sliceorder(Nsli,1));
        else
            rec_Nz = rec_Nz*SMS_fact;
            z = ((0:rec_Nz-1) - (rec_Nz-1)/2) * (rec_FOVz/rec_Nz);

            z = flip(z);

            % Shift NOT needed
            % if all(diff(twix.image.slicePos(3,:)) == 0)
            %     z_shift = twix.image.slicePos(3,1);
            %     z = z + z_shift;
            % elseif ~isempty(twix.hdr.Dicom.lScanRegionPositionIso)
            %     z_shift = twix.hdr.Dicom.lScanRegionPositionIso;
            %     z = z + z_shift;     
            % end

            z_old = z;

            z = (reshape(z, [Nsli SMS_fact]))'; 
            z =  z(:,tdr_sliceorder(Nsli,1));
            if mod(Nsli,2) ~= 0
                half = ceil(Nsli/2);
                z = [z(:,half+1:end), z(:,1:half)];
            end    

            z0 = z;

            tot_sli = Nsli*SMS_fact;
            Cgre = Cgre(:,:,:,tdr_sliceorder(tot_sli,2));
            B0map = B0map(:,:,tdr_sliceorder(tot_sli,2));

            [~, SMS_indices] = ismember(z, z_old);
            SMS_indices = SMS_indices';
        end    
    elseif SMS_fact == 1 && ~is3D
        % order z to match order in which slices were acquired
        z = z(:,tdr_sliceorder(Nsli,1));
    elseif SMS_fact > 1 && ~is3D
        % Establish z-coodinate pairs corresponding to multiband slices
        z_old = z;
        clear z;
        sli_index = start_index;
        for i = 1:Nsli
            switch orient
                case 1
                    z0(1,i) = twix.image.slicePos(3,sli_index);
                case 2
                    z0(1,i) = twix.image.slicePos(1,sli_index);
                case 3
                    z0(1,i) = twix.image.slicePos(2,sli_index);
            end
            z(1,i) = z0(1,i);
            if ~isempty(twix.hdr.Dicom.lScanRegionPositionIso)
                z_shift_4SMS = twix.hdr.Dicom.lScanRegionPositionIso;
                z(1,i) = z0(1,i) - z_shift_4SMS;
            end               
            for sli_group = 2:SMS_fact
                z(sli_group,i) = z(1,i) + (sli_group-1)*Nsli*deltaz;  
            end
            sli_index = sli_index + Nphs; 
        end
        tot_sli = Nsli*SMS_fact;
        Cgre = Cgre(:,:,:,tdr_sliceorder(tot_sli,2));
        B0map = B0map(:,:,tdr_sliceorder(tot_sli,2));

        [~, SMS_indices] = ismember(z, z_old);
        SMS_indices = SMS_indices';
    else
        % Flipping z for 3D acquisition
        z = flip(z); 
        % Ascending order
        Cgre = Cgre(:,:,:,tdr_sliceorder(Nsli,2));
        B0map = B0map(:,:,tdr_sliceorder(Nsli,2));
    end    

    if docoilcompress
        % Coil compression
    
        data = permute(data, [4 1 2 3]);
        Cgre = permute(Cgre, [3 1 2 4]);
    
        sz_Nenc = size(data,2);
        sz_Nx = size(Cgre,2);
        sz_Ny = size(Cgre,3);
    
        [comp_mtx, sv] = findCoilCompressMat(data(:,:,:,1),num_coil,0);
    
        for i = 1:Nsli
            data_tmp = squeeze(data(:,:,i,:));
            data_tmp = data_tmp(:,:);
            data_comp(:,:,i) = comp_mtx(:,:,i)*data_tmp;
        end

        if SMS_fact > 1
            for i = 1:Nsli
                Cgre_tmp = Cgre(:,:,:,SMS_indices(i,:));
                for j = 1:SMS_fact
                    Cgre_tmp2 = Cgre_tmp(:,:,:,j);
                    Cgre_tmp2 = Cgre_tmp2(:,:);
                    Cgre_comp(:,:,SMS_indices(i,j)) = comp_mtx(:,:,i)*Cgre_tmp2;
                end    
            end            
        else    
            for i = 1:Nsli
                Cgre_tmp = Cgre(:,:,:,i);
                Cgre_tmp = Cgre_tmp(:,:);
                Cgre_comp(:,:,i) = comp_mtx(:,:,i)*Cgre_tmp;
            end
        end    
    
        data_comp = reshape(data_comp, [num_coil sz_Nenc Nvol Nsli]);
        Cgre_comp = reshape(Cgre_comp, [num_coil sz_Nx sz_Ny Nsli*SMS_fact]);
    
        data = permute(data_comp, [2 4 3 1]);
        Cgre = permute(Cgre_comp, [2 3 1 4]);
    end 

    % Alternative SMS encoding. Not as reliable as the default 
    if SMS_fact > 1 && ~doSMSwithkz
        N_groups = tot_sli/SMS_fact;
        if docircshift 
            for k = 1:SMS_fact
                % slice indices belonging to this SMS "band"
                start_SMS = (k-1)*N_groups + 1;
                end_idx   = k*N_groups;
    
                % phase-encode shift in voxels
                shift_vox = (k-1)*rec_Ny/(R_fact*FOVshift);
                shift_vox = floor(shift_vox);
    
                % apply circular shift along PE (2nd) dimension
                Cgre(:,:,:,start_SMS:end_idx) = circshift(Cgre(:,:,:,start_SMS:end_idx),shift_vox,2);
                B0map(:,:,start_SMS:end_idx) = circshift(B0map(:,:,start_SMS:end_idx),shift_vox,2);                
            end
        else
            for k = 1:SMS_fact
                start_SMS = (k-1)*N_groups + 1;
                end_SMS = k*N_groups;
    
                shift_vox = (k-1)*rec_Ny/(R_fact*FOVshift);
    
                % integer PE indices for FFT along PE
                ky = 0:rec_Ny-1;                     
                ky = reshape(ky,1,rec_Ny,1,1);    
    
                % phase ramp along PE (dimension 2)
                phase_ramp = exp(1i*2*pi*ky*shift_vox/rec_Ny);
    
                % FFT along PE only (dimension 2)
                Cgre_tmp = fft(Cgre(:,:,:,start_SMS:end_SMS),[],2);
                B0map_tmp = fft(B0map(:,:,start_SMS:end_SMS),[],2);
    
                % multiply by phase ramp
                Cgre_tmp = Cgre_tmp.*phase_ramp;
                B0map_tmp = B0map_tmp.*phase_ramp;
    
                % IFFT along PE only
                Cgre(:,:,:,start_SMS:end_SMS) = ifft(Cgre_tmp,[],2);
                B0map(:,:,start_SMS:end_SMS) = ifft(B0map_tmp,[],2);
            end
        end
        Cgre = permute(Cgre, [1 2 4 3]);
    elseif SMS_fact > 1 && doSMSwithkz
        Cgre = permute(Cgre, [1 2 4 3]);
    end 

    % Identifying slice and volume with the highest signal, to be used for
    % skope delay synchronization
    data_sos = sos(abs(data));
    N_length = length(data_sos(:,1,1));  
    p = 0.10;  % 10 percent of raw data points
    numPoints = round(p * N_length); % Number of points to take
    if isEPI
        % Based around the central echo
        if is3D
            [~, maxIdx] = max(data_sos(:,Nsli/2,1));  % Index of max signal
        else    
            [~, maxIdx] = max(data_sos(:,1,1));  % Index of max signal
        end
        % Compute start and end indices
        startIdx = max(1, maxIdx - floor(numPoints / 2) + 1);
        endIdx = min(N_length, startIdx + numPoints - 1);

        data_10percent = data_sos(startIdx:endIdx,:,:);
    else
        data_10percent = data_sos(1:numPoints,:,:);
    end 
    if ~is3D
        sums = squeeze(sum(abs(data_10percent),1));
        [~, max_idx] = max(sums(:));
        [nsli_max, nvol_max] = ind2sub([Nsli, Nvol], max_idx);
    else
        sums = squeeze(sum(abs(data_10percent),1));
        sums = sum(sums,1);
        [~, nvol_max] = max(sums);
        nsli_max = Nsli/2; %Assume the center partition is approximately the max
    end    

    % Transform imaging grid based on imaging orientation
    coord.x = x;
    coord.y = y;
    coord.z = z;
    switch orient
        case 1     % transverse
            x = coord.x;
            y = coord.y;
            z = coord.z;

            % Create imaging voxel grid
            [X_grid,Y_grid] = meshgrid(single(x),single(y));
            X_grid = X_grid';
            Y_grid = Y_grid';

            % Complete the spatial grid
            grid.x = X_grid/1000;
            grid.x = repmat(grid.x,1,1,SMS_fact);
            grid.y = Y_grid/1000;
            grid.y = repmat(grid.y,1,1,SMS_fact);
            for i = 1:SMS_fact
                grid.z(:,:,i) = repmat(z(i,nsli_max),size(X_grid))/1000;
            end
        case 2     % sagittal
            x = coord.z;
            y = coord.y;
            z = coord.x;

            % Create imaging voxel grid
            [Z_grid,Y_grid] = meshgrid(single(z),single(y));
            Z_grid = Z_grid';
            Z_grid = flip(Z_grid,1);
            Y_grid = Y_grid';

            % Complete the spatial grid
            grid.z = Z_grid/1000;
            grid.z = repmat(grid.z,1,1,SMS_fact);
            grid.y = Y_grid/1000;
            grid.y = repmat(grid.y,1,1,SMS_fact);
            for i = 1:SMS_fact
                grid.x(:,:,i) = repmat(x(i,nsli_max),size(Z_grid))/1000;
            end
        case 3     % coronal
            x = coord.x;
            y = coord.z;
            z = coord.y;

            % Create imaging voxel grid
            [X_grid,Z_grid] = meshgrid(single(x),single(z));
            X_grid = -X_grid;
            Z_grid = -flip(Z_grid);

            % Complete the spatial grid
            grid.x = X_grid/1000;
            grid.x = repmat(grid.x,1,1,SMS_fact);
            grid.z = Z_grid/1000;
            grid.z = repmat(grid.z,1,1,SMS_fact);
            for i = 1:SMS_fact
                grid.y(:,:,i) = repmat(y(i,nsli_max),size(X_grid))/1000;
            end
        otherwise
            error('Oblique orientation not implemented')
    end

    if doaxisrotate
        slice_info = twix.hdr.MeasYaps.sSliceArray.asSlice{1};

        % Extract Siemens slice normal, replacing empty fields with zero
        fields = {'dSag', 'dCor', 'dTra'};
        
        n1 = zeros(3,1);
        
        for k = 1:3
            if isfield(slice_info.sNormal, fields{k}) && ...
                    ~isempty(slice_info.sNormal.(fields{k}))
        
                n1(k) = slice_info.sNormal.(fields{k});
            end
        end
        
        % Normalize
        if norm(n1) == 0
            error('Invalid slice normal: all components are empty or zero.');
        end

        % Original axial slice normal
        n0 = [0; 0; 1];
        
        % Calculate rotation axis and angle
        a = cross(n0, n1);
        a_norm = norm(a);
        
        theta = atan2(a_norm, dot(n0,n1));
        
        if a_norm > 1e-10
        
            a = a / a_norm;
        
            % Skew-symmetric cross-product matrix
            K = [  0,   -a(3),  a(2);
                  a(3),   0,   -a(1);
                 -a(2),  a(1),   0  ];
        
            % Rodrigues rotation matrix
            R = eye(3)*cos(theta) + ...
                (1-cos(theta))*(a*a') + ...
                sin(theta)*K;
        
        else
            R = eye(3);
        end
        
        % Apply rotation to existing spatial grid
        X0 = grid.x;
        Y0 = grid.y;
        Z0 = grid.z;
        
        grid.x = R(1,1)*X0 + R(1,2)*Y0 + R(1,3)*Z0;
        grid.y = R(2,1)*X0 + R(2,2)*Y0 + R(2,3)*Z0;
        grid.z = R(3,1)*X0 + R(3,2)*Y0 + R(3,3)*Z0;
    end    

    if doGNLcorr
        % Switch into LAI coordinate system which appears to be the
        % appropriate system for the GNL code
        switch orient
            case 1
                if ~ispulseq
                    grid.z = -grid.z;
                end
                if all(diff(twix.image.slicePos(1,:)) == 0)
                    x_shift = twix.image.slicePos(1,1)/1000;
                    grid.x = grid.x + x_shift;
                elseif ~isempty(twix.hdr.Dicom.lScanRegionPositionIso)
                    x_shift = twix.hdr.Dicom.lScanRegionPositionIso/1000;
                    grid.x = grid.x + x_shift;     
                end
                if all(diff(twix.image.slicePos(2,:)) == 0)
                    y_shift = twix.image.slicePos(2,1)/1000;
                    grid.y = grid.y + y_shift;
                elseif ~isempty(twix.hdr.Dicom.lScanRegionPositionIso)
                    y_shift = twix.hdr.Dicom.lScanRegionPositionIso/1000;
                    grid.y = grid.y + y_shift;     
                end
            case 2
                if ~ispulseq
                    grid.x = -grid.x;
                end
                if all(diff(twix.image.slicePos(2,:)) == 0)
                    y_shift = twix.image.slicePos(2,1)/1000;
                    grid.y = grid.y + y_shift;
                elseif ~isempty(twix.hdr.Dicom.lScanRegionPositionIso)
                    y_shift = twix.hdr.Dicom.lScanRegionPositionIso/1000;
                    grid.y = grid.y + y_shift;     
                end
                if all(diff(twix.image.slicePos(3,:)) == 0)
                    z_shift = twix.image.slicePos(3,1)/1000;
                    grid.z = grid.z + z_shift;
                elseif ~isempty(twix.hdr.Dicom.lScanRegionPositionIso)
                    z_shift = twix.hdr.Dicom.lScanRegionPositionIso/1000;
                    grid.z = grid.z + z_shift;     
                end
            case 3
                if ~ispulseq
                    grid.y = -grid.y;
                end
                if all(diff(twix.image.slicePos(1,:)) == 0)
                    x_shift = twix.image.slicePos(1,1)/1000;
                    grid.x = grid.x + x_shift;
                elseif ~isempty(twix.hdr.Dicom.lScanRegionPositionIso)
                    x_shift = twix.hdr.Dicom.lScanRegionPositionIso/1000;
                    grid.x = grid.x + x_shift;     
                end
                if all(diff(twix.image.slicePos(3,:)) == 0)
                    z_shift = twix.image.slicePos(3,1)/1000;
                    grid.z = grid.z + z_shift;
                elseif ~isempty(twix.hdr.Dicom.lScanRegionPositionIso)
                    z_shift = twix.hdr.Dicom.lScanRegionPositionIso/1000;
                    grid.z = grid.z + z_shift;     
                end
        end 

        grid.x = grid.x*1000;
        grid.y = grid.y*1000;
        grid.z = grid.z*1000;

        d = 0;
        % Supply Siemens gradient coefficient file
        gradfile = []; 
        method = 'direct';
        POLARITY = 1;

        % Not publicly available GNL unwarping code
        [grid.x, grid.y, grid.z, VJacDetLAI, Gx, Gy, Gz] = mris_gradient_nonlin__unwarp(grid.x, grid.y, grid.z, d, gradfile, method, POLARITY);
        grid.x = grid.x/1000;
        grid.y = grid.y/1000;
        grid.z = grid.z/1000;

        % Revert back to original coordinates
        switch orient
            case 1
                if ~ispulseq
                    grid.z = -grid.z;
                end
                if exist('x_shift')
                    grid.x = grid.x - x_shift;
                end
                if exist('y_shift')
                    grid.y = grid.y - y_shift;
                end  
            case 2
                if ~ispulseq
                    grid.x = -grid.x;
                end
                if exist('y_shift')
                    grid.y = grid.y - y_shift;
                end
                if exist('z_shift')
                    grid.z = grid.z - z_shift;
                end                  
            case 3
                if ~ispulseq
                    grid.y = -grid.y;
                end
                if exist('x_shift')
                    grid.x = grid.x - x_shift;
                end
                if exist('z_shift')
                    grid.z = grid.z - z_shift;
                end                  
        end     
    end

    % Determine maximum number of iterations based on the acquisition 
    if ~isfield(opt, 'maxNit') || isempty(opt.maxNit)
        if isEPI
            if R_fact == 1
                maxNit = 5;
            elseif R_fact == 2
                maxNit = 7;
            elseif R_fact > 2
                maxNit = 10;
            end    
        else
            maxNit = 20;
        end
    else 
        maxNit = opt.maxNit;
    end

    % Extract field dynamics: Now uses Skope NYOX data reader, download
    % from https://skope.swiss/software/ with valid license
    % session = Session(session_filepath);
    % scan = session.get_scan(scanId);
    % 
    % % Use precomputed skope profiles
    % if useskopefit
    %     % Skope z orientation is flipped relative to the imaging coordinate convention 
    %     scan.field_model.positions(:,3) = -scan.field_model.positions(:,3);
    %     if opt.fitOrder == 1
    %         scan.field_model.selected_spatial_basis(5:16) = 0;
    %     elseif opt.fitOrder == 2
    %         scan.field_model.selected_spatial_basis(10:16) = 0;
    %     end    
    % 
    %     ktraj = scan.trajectory;
    %     ktraj = ktraj - ktraj(1,:,:);
    %     kspha = ktraj(:,1:end-4,:);
    %     kconc = ktraj(:,end-3:end,:);   
    % else
    %    % Calculate using harmonicsFromRaw  
    %     probe_mag = scan.magnitude;
    %     probe_phase = scan.phase;
    % 
    %     % New dwell time for NYOX systems
    %     %dt_skope = 1e-6;
    %     dt_skope = 3.584e-6;
    %     probe_positions = scan.field_model.positions;
    %     fieldOffsets = scan.field_model.field_offsets;
    %     gammaProbes = session.system_config.gamma_probes;
        if ispulseq
            gammaMRI = seq.sys.gamma;
        else
            gammaMRI = twix.hdr.Config.Gamma_MHzT(1);
        end    
        B0 = 2.8937;
    %     B0 = session.system_config.B0;
    % 
    %     probe_phase = skope_demodulate_phase(probe_phase,gammaProbes,fieldOffsets,dt_skope);
    % 
        if strcmp(sys_type,'Connectome2')
            opt.coilParams.alph = 0.5;
            opt.coilParams.z0x = -0.1;
            opt.coilParams.z0y = -0.09;
            opt.coilParams.x0 = 0;
            opt.coilParams.y0 = 0;
        elseif strcmp(sys_type,'CimaX')
            opt.coilParams.alph = 0.5;
            opt.coilParams.z0x = 0;
            opt.coilParams.z0y = 0;
            opt.coilParams.x0 = 0;
            opt.coilParams.y0 = 0;
        elseif strcmp(sys_type,'TerraX')
            opt.coilParams.alph = 0.5;
            opt.coilParams.z0x = 0;
            opt.coilParams.z0y = 0;
            opt.coilParams.x0 = 0;
            opt.coilParams.y0 = 0;
        else
            error('Undefined scanner')
        end 
    % 
    %     % Perform fit 
    %     %[kspha, kconc] = harmonicsFromRaw_NYOX(probe_positions,probe_mag,probe_phase,dt_skope,fieldOffsets,gammaProbes,gammaMRI,B0,opt,isEPI);
    % end

    kspha = 2*pi*seq.calculateKspacePP();
    kspha(2:4,:) = kspha;
    kspha([1,4],:) = 0;
    if ~isEPI
        kspha = reshape(kspha, [4 Nenc Nsli*Nshots*Nvol]);
    else 
        kspha = reshape(kspha, [4 Nenc/Nshots Nsli*Nshots*Nvol]);
    end    
    for k = 1:size(kspha,3)
        kspha_vol = kspha(:,:,k);
        kconc(:,:,k) = compMaxwellPhase_nominal(kspha_vol(2:4,:).',dt*1e-6,gammaMRI,opt.coilParams,B0);
    end    
    kspha = permute(kspha, [2 1 3]);
    %kconc = zeros(size(kspha));

    % Polynomial trajectory extrapolation if probes experience low signal
    % during end of readout
    if dotrajextrap
        if size(kspha,2) == 4
            terms = 4;
        else
            terms = 9; % Don't extrapolate third order terms if present
        end  
        
        % How many points to use in the data fitting
        if ~isfield(opt,'good_end') || isempty(opt.good_end)
            good_end = size(kspha,1) - 0.2*size(kspha,1); % Use 80% of measurement window as default
        else
            good_end = opt.good_end;
        end

        dynamics = size(kspha,3);
        for vol = 1:dynamics
            for j = 1:terms
                poly_order = 2;
    
                traj = kspha(:,j,vol);
                t = (0:length(traj)-1)*dt_skope*1e6; 
                good_idx = 1:good_end;        % indices of good region
    
                % 1. Fit a smooth model to the good part
                p = polyfit(t(good_idx), traj(good_idx), poly_order); 
    
                % 2. Predict entire range
                traj_fit = polyval(p, t);
    
                % 3. Replace bad part
                traj_corrected(:,j,vol) = traj;
                traj_corrected(good_end+1:end,j,vol) = traj_fit(good_end+1:end);
            end   
        end    
        traj_corrected = reshape(traj_corrected, [size(kspha)]);
        % Keeping the original kx and ky
        kspha(:,[1,4:terms],:) = traj_corrected(:,[1,4:terms],:);
    end    
    
    Ncoeff_kspha = size(kspha,2);
    Ncoeff_kconc = size(kconc,2);
    skope_pts = size(kspha,1);

    % Remove reference dynamics
    tot_dynamics = Nsli*Nshots*Nvol;
    if tot_dynamics < size(kspha,3)
        kspha = kspha(:,:,end-tot_dynamics+1:end);
        kconc = kconc(:,:,end-tot_dynamics+1:end);
    end  

    % Zeroth order concomitant phase removal for asymmetric gradient coil
    % if strcmp(sys_type,'Connectome2')
    %     if ~MaxwellCompfromSim
    %         G = computeGradFromPhase(kspha(:,2:4,:),dt_skope,gammaProbes);
    %         if ~is3D
    %             G = reshape(G, [skope_pts 3 Nsli Nshots Nvol]);
    %             maxwell_zero = computeMaxwellTerms(0,opt.coilParams,G(:,1,Nsli,:,1),G(:,2,Nsli,:,1),G(:,3,Nsli,:,1),B0);
    %             maxwell_zero = maxwell_zero.zeroth;
    %             for i = 1:Nshots
    %                 maxwell_zero_sli(:,:,:,i,:) = repmat(maxwell_zero(:,:,:,i), [1 1 Nsli 1 1]);
    %             end  
    %             maxwell_zero = repmat(maxwell_zero_sli, [1 1 1 1 Nvol]);
    %             maxwell_zero = reshape(maxwell_zero, [size(kspha,1), 1, 1, size(kspha,3)]);
    %             phs_maxwell = cumsum((2*pi*gammaProbes*dt_skope)*maxwell_zero,1);
    %             kspha(:,1,:) = kspha(:,1,:) - phs_maxwell;
    %         else
    %             G = reshape(G, [skope_pts 3 Nshots Nsli Nvol]);
    %             maxwell_zero = computeMaxwellTerms(0,opt.coilParams,G(:,1,:,:,1),G(:,2,:,:,1),G(:,3,:,:,1),B0);
    %             maxwell_zero = maxwell_zero.zeroth;
    %             maxwell_zero = repmat(maxwell_zero, [1 1 1 1 Nvol]);
    %             maxwell_zero = reshape(maxwell_zero, [size(kspha,1), 1, 1, size(kspha,3)]);
    %             phs_maxwell = cumsum((2*pi*gammaProbes*dt_skope)*maxwell_zero,1);
    %             kspha(:,1,:) = kspha(:,1,:) - phs_maxwell;
    %         end
    %     end
    % end    

    %Dwell time has changed to 3.584 us with the NYOX systems (refer to new dwell time
    %in the manual)
    %tdwell = 1;
    tdwell = 3.584;

    % Determine delay time if not provided
    if ~isfield(opt, 'delSk') || isempty(opt.delSk)
        if isEPI
            % Perform coarse search using 2n+1 steps
            if ~isfield(opt, 'numCoarseSearch') || isempty(opt.numCoarseSearch)
                numCoarseSearch = 4;
            else
                numCoarseSearch = opt.numCoarseSearch;
            end    

            % Estimate the EPI trajectory start time
            % Find the first official kx = 0 crossing and subtract half the
            % echo spacing 
            echo_spc = Nfre*dt;

            if ispulseq && isEPI
                % Likely the case if frequency encoding partial Fourier is on
                if Nfre < N
                    N_oversamp = N*1.25; % Frequency encoding oversampling factor from example EPI pulseq dataset
                    echo_spc = N_oversamp*dt;
                end
            end
            
            switch orient
                case 1
                    kx = kspha(101:end,2,1);
                case 2
                    kx = kspha(101:end,4,1);
                case 3
                    kx = kspha(101:end,2,1);
            end

            t = 1:length(kx); % or your actual time vector

            % Find both crossing types
            zero_idx = find((kx(1:end-1) > 0 & kx(2:end) < 0) | (kx(1:end-1) < 0 & kx(2:end) > 0), 1);
            
            if ~isempty(zero_idx)
                x1 = kx(zero_idx);
                x2 = kx(zero_idx + 1);
                t1 = t(zero_idx);
                t2 = t(zero_idx + 1);
            
                % Linear interpolation to estimate zero crossing time
                t_zero = t1 - x1 * (t2 - t1) / (x2 - x1);
            else
                t_zero = NaN; % no crossing found
            end

            delSk0 = (t_zero + 100)*tdwell - (echo_spc/2);
        else
            % No course search needed for spiral 
            numCoarseSearch = 0;

            switch orient
                case 1
                    kx = kspha(:,2,1);
                case 2
                    kx = kspha(:,4,1);
                case 3
                    kx = kspha(:,2,1);
            end

            % Smooth the profile
            ktraj_smooth = movmean(kx, 15);  
            thresh = 0.2;
            delSk0 = find(abs(ktraj_smooth) > thresh, 1);
            delSk0 = delSk0*tdwell;            
        end
    end

    if SMS_fact > 1 && ~doSMSwithkz
        poly_order = 1;
        for dynamic = 1:size(kspha,3)
            p = polyfit(1:size(kspha,1), kspha(:,4,dynamic), poly_order);
            kspha(:,4,dynamic) = polyval(p, 1:size(kspha,1));
        end
    end 

    % Reshape field dynamics
    if ~is3D
        kspha = reshape(kspha, [skope_pts Ncoeff_kspha Nsli Nshots Nvol]);
        kconc = reshape(kconc, [skope_pts Ncoeff_kconc Nsli Nshots Nvol]);
    
        data_n = squeeze(data(:,nsli_max,nvol_max,:));
        kspha_n = squeeze(kspha(:,:,nsli_max,:,nvol_max)); 
        kconc_n = squeeze(kconc(:,:,nsli_max,:,nvol_max)); 
    else
        kspha = reshape(kspha, [skope_pts Ncoeff_kspha Nshots Nsli Nvol]);
        kconc = reshape(kconc, [skope_pts Ncoeff_kconc Nshots Nsli Nvol]);

        data_n = squeeze(data(:,:,nvol_max,:));
        data_n = reshape(data_n, [Nenc*Nsli, num_coil]);
        kspha_n = squeeze(kspha(:,:,:,:,nvol_max)); 
        kconc_n = squeeze(kconc(:,:,:,:,nvol_max)); 
    end

    opt.daNFull = size(data_n);
    opt.imNFull = [rec_Nx,rec_Ny,SMS_fact];
    
    % Supplying parameters to findDelAuto
    param.ispulseq = ispulseq;
    param.Ncoeff_kspha = Ncoeff_kspha;
    param.Ncoeff_kconc = Ncoeff_kconc;
    param.skope_pts = skope_pts;
    param.datatime_interp = datatime_interp;
    param.is3D = 0;
    param.rec_Nx = rec_Nx;
    param.rec_Ny = rec_Ny;
    param.SMS_fact = SMS_fact;
    param.isEPI = isEPI;
    param.Nenc = Nenc;
    param.ishighres = ishighres;
    param.doPhsConstraint = doPhsConstraint;
    param.encoding_approach = encoding_approach;
    param.subFact = [];
    param.doSMSwithkz = doSMSwithkz;    
    if ishighres
        if isEPI
            param.subFact = 10;
        else
            param.subFact = 5;
        end
    end
    if SMS_fact > 1
        param.z0 = z0(nsli_max);
    end   

    if is3D
        param.is3D = 1;
        % For multi-shot EPI, each shot contains an additional delay, which
        % is to be considered in the timing of the skope field dynamics
        MemBlock = cell2mat(twix.hdr.MeasYaps.sWipMemBlock.alFree);
        index = 8;
        for i = 1:Nshots
            if i == 1
                param.timeshift(i) = 0;
            else
                param.timeshift(i) = MemBlock(index);
                index = index + 1;
            end
        end     
    elseif ispulseq && Nshots > 1
        if isEPI || dospiraldelay
            % Delay interval determined from the seq file
            MemBlock = (0:Nshots-1)*TEshift;
            for i = 1:Nshots
                param.timeshift(i) = MemBlock(i);
            end 
        end
    end      

    if ~isEPI && ispulseq && ~dospiraldelay
        datatime = repmat(datatime, [Nshots 1]);
    elseif is3D
        if exist('MemBlock')
            base_time = datatime(:);  % Ensure it's a column vector
            Nt = length(base_time);
        
            % Expand datatime to Nshots × Nt
            datatime_all = repmat(base_time, Nshots, 1);
        
            % Create a delay vector of same shape
            delay_all = 1e-6*repelem(MemBlock(:), Nt);  % Nshots × Nt column vector
        
            % Apply the delays
            datatime_all_shifted = datatime_all + delay_all;
            datatime = datatime_all_shifted;
            datatime = repmat(datatime, [Nsli 1]);
        else
            datatime = repmat(datatime, [Nshots*Nsli 1]);
        end    
    else
        if exist('MemBlock')
            base_time = datatime(:);  % Ensure it's a column vector
            Nt = length(base_time);
            
            % Expand datatime to Nshots × Nt
            datatime_all = repmat(base_time, Nshots, 1);
            
            % Create a delay vector of same shape
            delay_all = 1e-6*repelem(MemBlock(:), Nt);  % Nshots × Nt column vector
            
            % Apply the delays
            datatime_all_shifted = datatime_all + delay_all;
            datatime = datatime_all_shifted;
        else
            datatime = repmat(datatime, [Nshots 1]);
        end
    end    

    % Synchronization delay calculation 
    if ~isfield(opt, 'delSk') || isempty(opt.delSk)
        if useGPU == 1 && ~is3D
            data_delAuto = gpuArray(data_n);
        else
            data_delAuto = data_n;
        end 
        if useSingle == 1
            data_delAuto = single(data_delAuto);
        end

        if SMS_fact > 1 && ~is3D 
            [delSk, delSk_perIt] = findDelAuto(delSk0,data_delAuto,tdwell,kspha_n,kconc_n,datatime,grid,B0map(:,:,SMS_indices(nsli_max,:)),Cgre(:,:,SMS_indices(nsli_max,:),:),maxNit,param,3,numCoarseSearch);
        else
            [delSk, delSk_perIt] = findDelAuto(delSk0,data_delAuto,tdwell,kspha_n,kconc_n,datatime,grid,B0map(:,:,nsli_max),Cgre(:,:,:,nsli_max),maxNit,param,3,numCoarseSearch);
        end 
    else
        delSk = opt.delSk;
    end    

    % Interpolate skope data, starting with delSk and interpolating at the
    % sampling times
    % if Nshots == 1 || ~isEPI
    %     [kspha,kconc] = interpTrajTime(kspha,kconc,tdwell,delSk,datatime_interp);
    % elseif is3D
    %     for i = 1:Nshots
    %         [kspha_interp(:,:,i,:,:),kconc_interp(:,:,i,:,:)] = interpTrajTime(kspha(:,:,i,:,:),kconc(:,:,i,:,:),tdwell,delSk+param.timeshift(i),datatime_interp);
    %     end    
    %     kspha = kspha_interp;
    %     kconc = kconc_interp;
    %     clear kspha_interp;
    %     clear kconc_interp;
    % else
    %     for i = 1:Nshots
    %         [kspha_interp(:,:,:,i,:),kconc_interp(:,:,:,i,:)] = interpTrajTime(kspha(:,:,:,i,:),kconc(:,:,:,i,:),tdwell,delSk+param.timeshift(i),datatime_interp);
    %     end    
    %     kspha = kspha_interp;
    %     kconc = kconc_interp;
    %     clear kspha_interp;
    %     clear kconc_interp;
    % end    

    if Nshots > 1 && ~is3D
        kspha = squeeze(permute(kspha, [2 1 4 3 5]));
        kconc = squeeze(permute(kconc, [2 1 4 3 5]));
    else
        kspha = squeeze(permute(kspha, [2 1 3 4 5]));
        kconc = squeeze(permute(kconc, [2 1 3 4 5]));
    end    

    % Replicate sampling matrix if more shots or slices (3D) were acquired
    if ~isEPI && ispulseq
        kspha = reshape(kspha, [Ncoeff_kspha Nenc*Nshots Nsli Nvol]);
        kconc = reshape(kconc, [Ncoeff_kconc Nenc*Nshots Nsli Nvol]);
    elseif is3D
        kspha = reshape(kspha, [Ncoeff_kspha Nenc Nsli Nvol]);
        kconc = reshape(kconc, [Ncoeff_kconc Nenc Nsli Nvol]);
    else
        kspha = reshape(kspha, [Ncoeff_kspha Nenc Nsli Nvol]);
        kconc = reshape(kconc, [Ncoeff_kconc Nenc Nsli Nvol]);
    end

    if retro_undersample_factor > 1
        original_datatime = datatime;
    end    

    clear grid;

    switch orient
        case 1     % transverse
            % Assign appropriate spatial grid coordinates 
            grid.x = X_grid/1000;
            grid.x = repmat(grid.x,1,1,1,SMS_fact);
            grid.y = Y_grid/1000;
            grid.y = repmat(grid.y,1,1,1,SMS_fact);
            for i = 1:SMS_fact
                grid.z(:,:,:,i) = repmat(reshape(z(i,:),[1,1,size(z,2)]),size(X_grid))/1000;
            end
            if doaxisrotate
                grid.x = repmat(grid.x,1,1,Nsli,SMS_fact);
                grid.y = repmat(grid.y,1,1,Nsli,SMS_fact);
            end            
        case 2     % sagittal
            grid.z = Z_grid/1000;
            grid.z = repmat(grid.z,1,1,1,SMS_fact);
            grid.y = Y_grid/1000;
            grid.y = repmat(grid.y,1,1,1,SMS_fact);
            for i = 1:SMS_fact
                grid.x(:,:,:,i) = repmat(reshape(x(i,:),[1,1,size(x,2)]),size(Z_grid))/1000;
            end
            if doaxisrotate
                grid.z = repmat(grid.z,1,1,Nsli,SMS_fact);
                grid.y = repmat(grid.y,1,1,Nsli,SMS_fact);
            end            
        case 3     % coronal
            grid.x = X_grid/1000;
            grid.x = repmat(grid.x,1,1,1,SMS_fact);
            grid.z = Z_grid/1000;
            grid.z = repmat(grid.z,1,1,1,SMS_fact);
            for i = 1:SMS_fact
                grid.y(:,:,:,i) = repmat(reshape(y(i,:),[1,1,size(y,2)]),size(X_grid))/1000;
            end
            if doaxisrotate
                grid.x = repmat(grid.x,1,1,Nsli,SMS_fact);
                grid.z = repmat(grid.z,1,1,Nsli,SMS_fact);
            end            
        otherwise
            error('Oblique orientation not implemented')
    end

    if doaxisrotate
        X0 = grid.x;
        Y0 = grid.y;
        Z0 = grid.z;
        
        grid_rotate.x = R(1,1)*X0 + R(1,2)*Y0 + R(1,3)*Z0;
        grid_rotate.y = R(2,1)*X0 + R(2,2)*Y0 + R(2,3)*Z0;
        grid_rotate.z = R(3,1)*X0 + R(3,2)*Y0 + R(3,3)*Z0;
    end    

    grid_OG = grid;

    if doGNLcorr
        disp('Performing GNL correction')

        switch orient
            case 1
                if ~ispulseq
                    grid.z = -grid.z;
                end
                if all(diff(twix.image.slicePos(1,:)) == 0)
                    x_shift = twix.image.slicePos(1,1)/1000;
                    grid.x = grid.x + x_shift;
                elseif ~isempty(twix.hdr.Dicom.lScanRegionPositionIso)
                    x_shift = twix.hdr.Dicom.lScanRegionPositionIso/1000;
                    grid.x = grid.x + x_shift;     
                end
                if all(diff(twix.image.slicePos(2,:)) == 0)
                    y_shift = twix.image.slicePos(2,1)/1000;
                    grid.y = grid.y + y_shift;
                elseif ~isempty(twix.hdr.Dicom.lScanRegionPositionIso)
                    y_shift = twix.hdr.Dicom.lScanRegionPositionIso/1000;
                    grid.y = grid.y + y_shift;     
                end
            case 2
                if ~ispulseq
                    grid.x = -grid.x;
                end
                if all(diff(twix.image.slicePos(2,:)) == 0)
                    y_shift = twix.image.slicePos(2,1)/1000;
                    grid.y = grid.y + y_shift;
                elseif ~isempty(twix.hdr.Dicom.lScanRegionPositionIso)
                    y_shift = twix.hdr.Dicom.lScanRegionPositionIso/1000;
                    grid.y = grid.y + y_shift;     
                end
                if all(diff(twix.image.slicePos(3,:)) == 0)
                    z_shift = twix.image.slicePos(3,1)/1000;
                    grid.z = grid.z + z_shift;
                elseif ~isempty(twix.hdr.Dicom.lScanRegionPositionIso)
                    z_shift = twix.hdr.Dicom.lScanRegionPositionIso/1000;
                    grid.z = grid.z + z_shift;     
                end
            case 3
                if ~ispulseq
                    grid.y = -grid.y;
                end
                if all(diff(twix.image.slicePos(1,:)) == 0)
                    x_shift = twix.image.slicePos(1,1)/1000;
                    grid.x = grid.x + x_shift;
                elseif ~isempty(twix.hdr.Dicom.lScanRegionPositionIso)
                    x_shift = twix.hdr.Dicom.lScanRegionPositionIso/1000;
                    grid.x = grid.x + x_shift;     
                end
                if all(diff(twix.image.slicePos(3,:)) == 0)
                    z_shift = twix.image.slicePos(3,1)/1000;
                    grid.z = grid.z + z_shift;
                elseif ~isempty(twix.hdr.Dicom.lScanRegionPositionIso)
                    z_shift = twix.hdr.Dicom.lScanRegionPositionIso/1000;
                    grid.z = grid.z + z_shift;     
                end
        end 

        if ~doaxisrotate
            switch orient
                case 1
                    grid_GNL.x = repmat(grid.x,1,1,size(grid.z,3),1)*1000;
                    grid_GNL.y = repmat(grid.y,1,1,size(grid.z,3),1)*1000;
                    grid_GNL.z = grid.z*1000;
                case 2
                    grid_GNL.x = grid.x*1000;
                    grid_GNL.y = repmat(grid.y,1,1,size(grid.x,3),1)*1000;
                    grid_GNL.z = repmat(grid.z,1,1,size(grid.x,3),1)*1000;
                case 3
                    grid_GNL.x = repmat(grid.x,1,1,size(grid.y,3),1)*1000;
                    grid_GNL.y = grid.y*1000;
                    grid_GNL.z = repmat(grid.z,1,1,size(grid.y,3),1)*1000;
                otherwise
                    error('Oblique orientation not implemented')
            end
        else
            grid_GNL.x = grid.x*1000;
            grid_GNL.y = grid.y*1000;
            grid_GNL.z = grid.z*1000;
        end    

        if SMS_fact > 1
            grid_GNL.x = grid_GNL.x(:,:,:);
            grid_GNL.y = grid_GNL.y(:,:,:);
            grid_GNL.z = grid_GNL.z(:,:,:);
        end    

        d = 0;
        % Supply Siemens gradient coefficient file
        gradfile = [];   
        method = 'direct';
        POLARITY = 1;

        % Not publicly available GNL unwarping code
        [grid_GNL.x, grid_GNL.y, grid_GNL.z, VJacDetLAI, Gx, Gy, Gz] = mris_gradient_nonlin__unwarp(grid_GNL.x, grid_GNL.y, grid_GNL.z, d, gradfile, method, POLARITY);
        grid_GNL.x = grid_GNL.x/1000;
        grid_GNL.y = grid_GNL.y/1000;
        grid_GNL.z = grid_GNL.z/1000;

        if SMS_fact > 1
            switch orient
                case 1 
                    grid_GNL.x = reshape(grid_GNL.x, size(grid.z));
                    grid_GNL.y = reshape(grid_GNL.y, size(grid.z));
                    grid_GNL.z = reshape(grid_GNL.z, size(grid.z));
                case 2
                    grid_GNL.x = reshape(grid_GNL.x, size(grid.x));
                    grid_GNL.y = reshape(grid_GNL.y, size(grid.x));
                    grid_GNL.z = reshape(grid_GNL.z, size(grid.x));
                case 3    
                    grid_GNL.x = reshape(grid_GNL.x, size(grid.y));
                    grid_GNL.y = reshape(grid_GNL.y, size(grid.y));
                    grid_GNL.z = reshape(grid_GNL.z, size(grid.y));
            end
        end

        switch orient
            case 1
                if ~ispulseq
                    grid_GNL.z = -grid_GNL.z;
                end
                if exist('x_shift')
                    grid_GNL.x = grid_GNL.x - x_shift;
                end
                if exist('y_shift')
                    grid_GNL.y = grid_GNL.y - y_shift;
                end 
            case 2
                if ~ispulseq
                    grid_GNL.x = -grid_GNL.x;
                end
                if exist('y_shift')
                    grid_GNL.y = grid_GNL.y - y_shift;
                end      
                if exist('z_shift')
                    grid_GNL.z = grid_GNL.z - z_shift;
                end                
            case 3
                if ~ispulseq
                    grid_GNL.y = -grid_GNL.y;
                end
                if exist('x_shift')
                    grid_GNL.x = grid_GNL.x - x_shift;
                end
                if exist('z_shift')
                    grid_GNL.z = grid_GNL.z - z_shift;
                end                 
        end
    end

    if SMS_fact > 1 && doSMSwithkz && ~ispulseq
        for i = 1:Nsli
            phs = kspha(4,:,i,:)*(z0(i)/1000);
            phs = squeeze(phs);
            phs = reshape(phs, size(phs, 1), 1, size(phs,2));
            data(:,i,:,:) = data(:,i,:,:).*exp(1i*phs);
        end
    end    

    % Image reconstruction loop over defined slices and volumes
    for volume = input_vol
        for nslice = input_sli
            if Nshots > 1 && ~is3D
                subFact = [];
                if dob0reg
                    doConventionalRecon = 0;
                    useGPU = 1;
                    subFact = 20;  
                    if ismember(volume, b0_indices) 
                        rank = b0_rank;
                        lambda = b0_lambda; 
                    else
                        rank = dwi_rank;
                        lambda = dwi_lambda; 
                    end
                else
                    if ismember(volume, b0_indices) 
                        doConventionalRecon = 1;
                        useGPU = 0;
                        if isEPI
                            subFact = 10;
                        else
                            subFact = 1;
                        end   
                    else
                        doConventionalRecon = 0;
                        useGPU = 1;
                        subFact = 20; 
                        rank = dwi_rank;
                        lambda = dwi_lambda;
                    end
                end
            end       

            if ~is3D
                data_n = squeeze(data(:,nslice,volume,:));
                kspha_n = squeeze(kspha(:,:,nslice,volume)); 
                kconc_n = squeeze(kconc(:,:,nslice,volume)); 
            else
                data_n = squeeze(data(:,:,volume,:));
                kspha_n = squeeze(kspha(:,:,:,volume)); 
                kconc_n = squeeze(kconc(:,:,:,volume)); 

                data_n = reshape(data_n, [Nenc*Nsli, num_coil]);
                kspha_n = reshape(kspha_n, [Ncoeff_kspha, Nenc*Nsli]);
                kconc_n = reshape(kconc_n, [Ncoeff_kconc, Nenc*Nsli]);
            end    

            if doGNLcorr
                clear grid;
                for i = 1:SMS_fact
                    grid.x(:,:,i) = grid_GNL.x(:,:,nslice,i);
                    grid.y(:,:,i) = grid_GNL.y(:,:,nslice,i);
                    grid.z(:,:,i) = grid_GNL.z(:,:,nslice,i);
                end
            elseif doaxisrotate
                clear grid;
                for i = 1:SMS_fact
                    grid.x(:,:,i) = grid_rotate.x(:,:,nslice,i);
                    grid.y(:,:,i) = grid_rotate.y(:,:,nslice,i);
                    grid.z(:,:,i) = grid_rotate.z(:,:,nslice,i);
                end                
            else
                clear grid;
                switch orient
                    case 1
                        grid.x = squeeze(grid_OG.x);
                        grid.y = squeeze(grid_OG.y);
                        for i = 1:SMS_fact
                            grid.z(:,:,i) = grid_OG.z(:,:,nslice,i);
                        end
                    case 2
                        for i = 1:SMS_fact
                            grid.x(:,:,i) = grid_OG.x(:,:,nslice,i);
                        end
                        grid.y = squeeze(grid_OG.y);
                        grid.z = squeeze(grid_OG.z);
                    case 3
                        grid.x = squeeze(grid_OG.x);
                        for i = 1:SMS_fact
                            grid.y(:,:,i) = grid_OG.y(:,:,nslice,i);
                        end
                        grid.z = squeeze(grid_OG.z);
                    otherwise
                         error('Oblique orientation not implemented')
                end
            end    

            opt.daNFull = size(data_n);
            opt.imNFull = [rec_Nx,rec_Ny,SMS_fact];
            
            % Build operator for receiver coils
            if SMS_fact == 1
                R = rcvrOp(Cgre(:,:,:,nslice),false);
            else
                R = rcvrOp(Cgre(:,:,SMS_indices(nslice,:),:),false);
            end
            
            % retrospective undersampling
            if retro_undersample_factor>1
                if isEPI
                    NsamplesPerShot = Nfre*Npershot;
                else    
                    NsamplesPerShot = Nenc;
                end

                datatime = original_datatime;

                tmp_datatime = reshape(datatime,[],NsamplesPerShot,Nshots);
                datatime = tmp_datatime(:,:,1:retro_undersample_factor:end);
                Nshots_effective = size(datatime,3);
                % reshape back to be compatible with the original code
                datatime = reshape(datatime, NsamplesPerShot*Nshots_effective,1);
                
                kspha_n = reshape(kspha_n, [], NsamplesPerShot,Nshots);
                tmp_kspha_n = kspha_n(:,:,1:retro_undersample_factor:end);
                kspha_n = reshape(tmp_kspha_n, [], NsamplesPerShot*Nshots_effective);
                
                kconc_n = reshape(kconc_n, [], NsamplesPerShot,Nshots);
                tmp_kconc_n = kconc_n(:,:,1:retro_undersample_factor:end);
                kconc_n = reshape(tmp_kconc_n, [], NsamplesPerShot*Nshots_effective);
                
                data_n = reshape(data_n,NsamplesPerShot,Nshots,num_coil);
                data_n = data_n(:,1:retro_undersample_factor:end,:);
                data_n = reshape(data_n,NsamplesPerShot*Nshots_effective,num_coil);
                
                opt.daNFull = size(data_n);
            else
                Nshots_effective = Nshots;
            end

            if doConventionalRecon
                % Build sampling operator 
                if SMS_fact > 1 && ~is3D
                    S = sampHighOrder(B0map(:,:,SMS_indices(nslice,:)),datatime,kspha_n,kconc_n,grid,[],useGPU,useSingle,encoding_approach,opt);
                elseif is3D
                    useSingle = 1;
                    useGPU = 0;
                    opt.subFact = 20;
                    S = sampHighOrder(B0map(:,:,nslice),datatime,kspha_n,kconc_n,grid,[],useGPU,useSingle,encoding_approach,opt);
                elseif ishighres
                    opt.subFact = subFact;
                    S = sampHighOrder(B0map(:,:,nslice),datatime,kspha_n,kconc_n,grid,[],useGPU,useSingle,encoding_approach,opt);    
                else
                    S = sampHighOrder(B0map(:,:,nslice),datatime,kspha_n,kconc_n,grid,[],useGPU,[],encoding_approach,opt);
                end     
    
                opFunc = @(x,transp) mrSampFunc(x,transp,S,R,opt); % mrSampFunc creates the net MRI sampling operation
                %[x,FLAG,RELRES,ITER,RESVEC] = lsqr(opFunc,data_n(:),[],maxNit,[],[]);
                % Iteratively solve using the conjugate gradient method
                if doL2reg
                    [x, resSqAll, mseAll, xnormAll, xdiffAll, stopThresh, status] = cgne(opFunc,data_n(:),[],maxNit,[],L2gamma,1);
                elseif doCSreg
                    maxNit_CS = 800;
                    opt.interpThresh = 0;
                    opt.stopOnResInc = 0;
                    opt.resThresh = 1e-4;

                    ls_temp = cgne(opFunc,data_n(:),[],maxNit,opt); 
                    ls = reshape(ls_temp,opt.imNFull);
                    
                    levelsPerDim = [1, 1];
                    isDec = 0; % Non-decimated to avoid blocky edges
                    W = dwtND(levelsPerDim,opt.imNFull,isDec);
                
                    wim = W*ls;
                
                    lambdaMethod = 0;
                    if isempty(lamVal_CS)
                        lamVal_CS = getLambdas(wim,lambdaMethod);
                    end    
                
                    lambda_CS = waveletObj(lamVal_CS(1),wim); % This replicates the scalar lamVal into all wavelet fields using the template wim_zf
                    lambda_CS = fillHighFreqWithLambdas(lamVal_CS,lambda_CS);
                    lambda_CS.low = 0;
                
                    % Do a few cgne iterations to get a better starting guess
                    x0 = cgne(opFunc,data_n(:),[],2); 
                    x0 = reshape(x0,opt.imNFull); 

                    cs_temp = bfista(opFunc,data_n(:),W,lambda_CS,x0(:),maxNit_CS,opt);
                    x = gather(cs_temp);
                else
                    [x, resSqAll, mseAll, xnormAll, xdiffAll, stopThresh, status] = cgne(opFunc,data_n(:),[],maxNit);
                end 

                % Phase Constraint option for partial Fourier EPI
                if doPhsConstraint
                    x0 = reshape(x, opt.imNFull);
                    clear x;
                    kspace_firstit = fft2cmulti_norm(x0);
                    sz = size(kspace_firstit,1);
    
                    % Keep 40% of k-space data
                    keep_frac = 0.40;               
                    keep_N = round(sz * keep_frac);    % number of samples to keep
                    
                    % Centered truncation (keep symmetric around middle)
                    kspace_firstit_trunc = zeros(size(kspace_firstit));
                    mid = floor(sz/2) + 1;
                    half_keep = floor(keep_N/2);
                    idx_keep = (mid-half_keep):(mid+half_keep);
                    
                    % Avoid indexing outside range
                    idx_keep = max(idx_keep,1);
                    idx_keep = min(idx_keep,sz);
                    
                    kspace_firstit_trunc(idx_keep,idx_keep) = kspace_firstit(idx_keep,idx_keep);
            
                    % Obtain smooth phase map from low-res image
                    im_firstit = ifft2cmulti_norm(kspace_firstit_trunc);
                    phs_firstit = unwrap_phase(angle(im_firstit));
                    phs_firstit = imgaussfilt(phs_firstit, 3.5);
    
                    % Extract real component 
                    im_real = x0.*exp(-1i*phs_firstit);
    
                    % Absorb the measured phase factor into the encoding operator
                    if SMS_fact == 1
                        Cgre_phi = Cgre(:,:,:,nslice).*exp(1i*phs_firstit);
                    else
                        Cgre_phi = Cgre(:,:,SMS_indices(nslice,:),:).*exp(1i*phs_firstit);
                    end
                    R_phi = rcvrOp(Cgre_phi,false);
                    
                    % Only a few more iterations needed, given already
                    % using a starting guess
                    phsconst_Nit = 3;
                    opFunc = @(x,transp) mrSampFunc(x,transp,S,R_phi,opt); % mrSampFunc creates the net MRI sampling operation

                    % Iteratively solve using the conjugate gradient method
                    [x, resSqAll, mseAll, xnormAll, xdiffAll, stopThresh, status] = cgne(opFunc,data_n(:),im_real(:),phsconst_Nit);
                end
            else
                 %% %%%%%%%%%%%%%%%%%%Begin S-LORAKS
                % (1) for each shot, do SENSE separately
                % overwrite data size, this is the size of a single shot now
                sz = size(data_n);
                if isEPI
                    opt.daNFull = [Nfre*Npershot,num_coil];
                else    
                    opt.daNFull = [Nenc,num_coil];
                end    
    
                my_params = [];
                my_params.Nshots=Nshots_effective;
                if SMS_fact > 1
                    my_params.B0map=B0map(:,:,SMS_indices(nslice,:));
                else
                    my_params.B0map=B0map(:,:,nslice);
                end  
                my_params.datatime=datatime;
                my_params.kspha_n=kspha_n;
                my_params.kconc_n=kconc_n;
                my_params.grid=grid;
                my_params.useGPU=useGPU;
                my_params.opt=opt;
                my_params.R=R;
                my_params.encoding_approach = encoding_approach;
                my_params.subFact = subFact;                
                % my_params.Nenc is number of samples in each shot
                if isEPI
                    my_params.Nenc=Nfre*Npershot;
                else
                    my_params.Nenc=Nenc;
                end    
                my_params.isEPI=isEPI;
                my_params.is3D=is3D;
                my_params.useSingle = useSingle;
                my_params.SMS_fact = SMS_fact;
                
                A_for = @(in)my_apply_highOrderSamp(in, 'notransp',my_params);
                A_adj = @(in)my_apply_highOrderSamp(in, 'transp',my_params);
                AHA   = @(in)A_adj(A_for(in));
                AHd   = A_adj(data_n);           
                sens_result = pcg(AHA,AHd);
                sens_result = reshape(sens_result,opt.imNFull(1),opt.imNFull(2),Nshots_effective);
    
    
                % (2) S-LORAKS 
                N1 = rec_Nx;
                N2 = rec_Ny;
                Ns = Nshots_effective;
    
                LORAKS_type = 1; % S-matrix (Support + Phase + parallel imaging)
                Radius = 3;
                tol = 1e-4;
    
                [in1,in2] = meshgrid(-Radius:Radius,-Radius:Radius);
                idx = find(in1.^2+in2.^2<=Radius^2);
                patchSize = numel(idx);            
    
                B = @(x) ft2(x); %
                Bh = @(x)  ift2(x); %
    
                P_M = @(x) LORAKS_operators(x,N1,N2,Ns,Radius,LORAKS_type,[]);
                Ph_M = @(x) LORAKS_operators(x,N1,N2,Ns,Radius,-LORAKS_type,[]);
    
                ZD = @(x) padarray(reshape(x,[N1 N2 Ns]),[2*Radius, 2*Radius], 'post');
                ZD_H = @(x) x(1:N1,1:N2,:,:);
    
                z = sens_result;
                if useSingle
                    z = single(z);
                end
                if useGPU
                    z = gpuArray(z);
                end            
                for iter = 1:max_iter
                    z_cur = z;
                    pz = B(z);
                    MM = P_M(pz);
                    l = sum(MM,1);
                    l = find(abs(l)>0);
                    MM = MM(:,l);
    
                    Um = svd_left(MM);
                    nmm = Um(:,rank+1:end)'; % null space
                    Bhr = 0;
    
                    if 0
                        [~,s] = svd(MM);
                        figure;
                        plot(diag(s)./max(diag(s)))
                    end
    
                    if LORAKS_type == 1 % S
                        nf = size(nmm,1);
                        nmm = reshape(nmm,[nf, patchSize, 2*Ns]);
                        nss_h = reshape(nmm(:,:,1:2:end)+1j*nmm(:,:,2:2:end),[nf, patchSize*Ns]);
                        Nis = filtfilt(nss_h,'C',N1,N2,Ns,Radius);
                        Nis2 = filtfilt(nss_h,'S',N1,N2,Ns,Radius);
    
                        L1 = @(x) ZD_H(ifft2(squeeze(sum(Nis.*repmat(fft2(ZD(B(x))),[1 1 1 Ns]),3))));
                        L2 = @(x) ZD_H(ifft2(squeeze(sum(Nis2.*repmat(conj(fft2(ZD(B(x)))),[1 1 1 Ns]),3))));
                        LhL = @(x) 2*Bh(reshape(L1(x)-L2(x),[N1 N2 Ns 1]));
                    end
    
                    % data fitting
                    % M = @(x) AtA(x) + lambda*LhL(x);
                    reshape_input = @(x) reshape(x,[N1 N2 Ns]);
                    vec = @(x) x(:);
                    M = @(x) vec(  AHA(reshape_input(x)) + vec(lambda*LhL(reshape_input(x)))  ); % for pcg
                    tic
                    [z,~] = pcg(M, AHd(:) + lambda*Bhr);
                    fprintf('LORAKS, pcg() took %.1f seconds\n', toc);
                    z = reshape_input(z);
    
                    t = (norm(z_cur(:)-z(:))/norm(z(:)));
    
                    % display the status
                    if ~rem(iter,1)
                        disp(['iter ' num2str(iter) ', relative change in solution: ' num2str(t)]);
                    end
    
                    % check for convergence
                    if t < tol
                        disp('Convergence tolerance met: change in solution is small');
                        break;
                    end    
                end
                im_S_LORAKS = z;
                % %%%%%%%%%%%%%%%%%%End S-LORAKS
            end

            if doGibbsCor
                Gibbs_alpha = 1;
                if doConventionalRecon
                    Gibbs_filt = filtNd(opt.imNFull(1:2),Gibbs_alpha,2,'kb');
                    x = reshape(ifftnc(Gibbs_filt.*fftnc(reshape(x,opt.imNFull),2),2),[size(x,1) 1]);
                else
                    Gibbs_filt = filtNd(opt.imNFull(1:2),Gibbs_alpha,2,'kb');
                    im_S_LORAKS = ifftnc(Gibbs_filt.*fftnc(im_S_LORAKS,2),2);
                end
            end   

            if ~isEPI
                if doConventionalRecon
                    kspace = fftnc(reshape(x,opt.imNFull),2);
                    R = N/2;
                    
                    [x_coordinate,y_coordinate] = meshgrid( ...
                        (-N/2):(N/2-1), ...
                        (-N/2):(N/2-1));
                    
                    %Make slightly larger than the nominal threshold
                    Rmask = 1.02 * R;
                    mask = sqrt(x_coordinate.^2 + y_coordinate.^2) <= Rmask;
    
                    kspace = kspace.*mask;
                    x = reshape(ifftnc(kspace,2),[size(x,1) 1]);
                else
                    kspace = fftnc(im_S_LORAKS,2);
                    R = N/2;
                    
                    [x_coordinate,y_coordinate] = meshgrid( ...
                        (-N/2):(N/2-1), ...
                        (-N/2):(N/2-1));
                    
                    Rmask = 1.02 * R;
                    mask = sqrt(x_coordinate.^2 + y_coordinate.^2) <= Rmask;
    
                    kspace = kspace.*mask;
                    im_S_LORAKS = ifftnc(kspace,2);
                end     
            end     

            if SMS_fact == 1
                if doparallelcomp
                    if doConventionalRecon
                        im = reshape(x, opt.imNFull);
                    else
                        imgs_S_LORAKS = mean(abs(im_S_LORAKS),3);
                        if useGPU
                            imgs_S_LORAKS = gather(imgs_S_LORAKS);
                        end
                        im = imgs_S_LORAKS;
                    end
                else
                    if doConventionalRecon
                        im(:,:,nslice,volume) = reshape(x, opt.imNFull);
                    else
                        imgs_S_LORAKS = mean(abs(im_S_LORAKS),3);
                        if useGPU
                            imgs_S_LORAKS = gather(imgs_S_LORAKS);
                        end
                        im(:,:,nslice,volume) = imgs_S_LORAKS;
                    end
                end
            elseif doSMSwithkz
                im(:,:,SMS_indices(nslice,:),volume) = reshape(x, opt.imNFull);            
            else    
                im(:,:,SMS_indices(nslice,:),volume) = reshape(x, opt.imNFull);
                if docircshift
                    for k = 1:SMS_fact
                        % slice indices belonging to this SMS "band"
                        start_SMS = (k-1)*N_groups + 1;
                        end_idx   = k*N_groups;
    
                        % phase-encode shift in voxels
                        shift_vox = (k-1)*rec_Ny/(R_fact*FOVshift);
                        shift_vox = floor(shift_vox);
    
                        % apply circular shift along PE (2nd) dimension
                        im(:,:,SMS_indices(nslice,k),volume) = circshift(im(:,:,SMS_indices(nslice,k),volume),-shift_vox,2);
                    end
                else
                    for k = 1:SMS_fact
                        start_SMS = (k-1)*N_groups + 1;
                        end_SMS = k*N_groups;
    
                        shift_vox = -(k-1)*rec_Ny/(R_fact*FOVshift);
    
                        % integer PE indices for FFT along PE
                        ky = 0:rec_Ny-1;                   
                        ky = reshape(ky,1,rec_Ny,1,1);     
    
                        % phase ramp along PE (dimension 2)
                        phase_ramp = exp(1i*2*pi*ky*shift_vox/rec_Ny);
    
                        % FFT along PE only (dimension 2)
                        tmp = fft(im(:,:,SMS_indices(nslice,k),volume),[],2);
    
                        % multiply by phase ramp
                        tmp = tmp.*phase_ramp;
    
                        % IFFT along PE only
                        im(:,:,SMS_indices(nslice,k),volume) = ifft(tmp,[],2);
                    end
                end
            end    
        end
    end

    % MPPCA denoising 
    if dotensordenoise
        disp(['Performing tensor MPPCA denoising. Using patch size: ', num2str(denoise_patchsize)]);
        im_denoised = denoise_recursive_tensor(im, [denoise_patchsize]);
    end    
end

function [out] = ft2( kdata )
    %2DFT Summary of this function goes here
    %   Detailed explanation goes here
    [nx, ny, ~] = size(kdata);
    out = fftshift(fft2(ifftshift(kdata))) / sqrt(nx*ny);
end

function [out] = ift2( data )
    %2DFT Summary of this function goes here
    %   Detailed explanation goes here
    [nx, ny, ~] = size(data);
    out = fftshift(ifft2(ifftshift(data))) * sqrt(nx*ny);
end








