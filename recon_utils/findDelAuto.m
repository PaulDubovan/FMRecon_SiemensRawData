function [delSk, delSk_perIt] = findDelAuto(delSk0,data_in,tdwell,phs_spha,phs_conc,datatime,phs_sph_grid,b0map,Crcvr,maxNit_cgne,param,delJumpFact,numCoarseSearch)
    % Find optimal delay for traj data with respect to MRI raw data
    % Presumes a 2D slice 
    % Units for all times should be consistent (tdwell [time], datatime [time], b0map [rad/time])
    %
    % Inputs:
    %   delSk0:     starting guess for delay
    %   data_in:    MRI raw data
    %   tdwell:     dwell time for trajectory samples. Time for trajectory
    %               is presumed to start at 0.
    %   phs_spha:   trajectory using a higher order encoding model. See
    %               sampHighOrder class in matmri repository
    %               (https://doi.org/10.5281/zenodo.4495477)
    %   phs_conc:   phase from concomitant gradient fields. See sampHighOrder.
    %   datatime:   times for each data point in data_in
    %   phs_sph_grd: x, y, and z positions for each voxel. See sampHighOrder.
    %   b0map:      b0map. See sampHighOrder.
    %   Crcvr:      receiver sensitivity profiles. First three dims should
    %               have same size as b0map, with one extra 4th dim to
    %               count the number of receivers.
    %   maxNit_cgne: maximum number of iterations to use for each
    %                update of the image
    %   numCourseSearch: number of samples to perform a course search over
    %               e.g., for a value of 2, the algorithm finds which of
    %               the delays n*(datatime(2)-datatime(1)) gives the lowest
    %               residual when finding the image, where n spans from
    %               -numCourseSearch to +numCourseSearch. The iterative
    %               optimization then starts from that delay.
    %
    % (c) Corey Baron 2021
    % Modificiations made by Paul Dubovan, 2026

    % Options
    if nargin<13
        numCoarseSearch = 2;
    end
    if nargin<12
        delJumpFact = 3;
    end

    % Assign acquisition parameters
    ispulseq = param.ispulseq;
    Nshots = param.Nshots;
    Ncoeff_kspha = param.Ncoeff_kspha;
    Ncoeff_kconc = param.Ncoeff_kconc;
    skope_pts = param.skope_pts;
    datatime_interp = param.datatime_interp;
    is3D = param.is3D;
    isEPI = param.isEPI;
    ishighres = param.ishighres;
    SMS_fact = param.SMS_fact;
    rec_Nx = param.rec_Nx;
    rec_Ny = param.rec_Ny;
    Nenc = param.Nenc;
    Nsli = param.Nsli;
    subFact = param.subFact;
    doPhsConstraint = param.doPhsConstraint;
    encoding_approach = param.encoding_approach;
    doSMSwithkz = param.doSMSwithkz;
    if Nshots > 1 && isEPI
        timeshift = param.timeshift;
    end
    if SMS_fact > 1
        z0 = param.z0;
    end   

    % Check for data types
    useSingle = 0;
    useGPU = 0;
    if isa(data_in,'gpuArray')
        useGPU = 1;
        if isaUnderlying(data_in,'single')
            useSingle = 1;
        end
    else
        if isa(data_in,'single')
            useSingle = 1;
        end
    end

    if SMS_fact > 1 && doSMSwithkz
        data_ref_SMS = data_in;
    end 
    
    tic1 = tic;
    fprintf('findDelAuto: Init guess %.2f; ', delSk0)
    
    % Basic prep computations
    numsamps = numCoarseSearch; % number of samples to try in initial course search
    opt_cgne.stopOnResInc = 1; % this never happens until noise starts amplifying, so it makes a good stopping criteria here.
    % Create receiver operator
    if useSingle
        R = rcvrOp(single(Crcvr),0); 
    else
        R = rcvrOp(Crcvr,0);
    end
    approach = 3;    % Disable interpolation for iterations, since it could affect the results.
    opt.segmentedNdiv = 4;
    dtsamp = datatime(2)-datatime(1);

    % Find numerical derivative of phs_spha and phs_conc
    % Need to reorganize before taking the derivative if multiple shots
    % were acquired
    kern = [1/8 1/4 0 -1/4 -1/8]';
    if is3D
        phs_spha_dt = permute(phs_spha, [2 1 3 4]);
        phs_conc_dt = permute(phs_conc, [2 1 3 4]);
        phs_spha_dt = phs_spha_dt(:,:);
        phs_conc_dt = phs_conc_dt(:,:);
        dspha_dt = conv2(phs_spha_dt', kern, 'same')/tdwell; % units rad/{1,m,m^2,m^3}/time (spatial units depend on spatial order of term)
        dconc_dt = conv2(phs_conc_dt', kern, 'same')/tdwell;
        dspha_dt = dspha_dt';
        dconc_dt = dconc_dt';
        dspha_dt = reshape(dspha_dt, [Ncoeff_kspha skope_pts Nshots Nsli]);
        dconc_dt = reshape(dconc_dt, [Ncoeff_kconc skope_pts Nshots Nsli]);
        dspha_dt = squeeze(permute(dspha_dt, [2 1 3 4]));
        dconc_dt = squeeze(permute(dconc_dt, [2 1 3 4]));
    elseif Nshots > 1 && ~is3D    
        phs_spha_dt = permute(phs_spha, [2 1 3]);
        phs_conc_dt = permute(phs_conc, [2 1 3]);
        phs_spha_dt = phs_spha_dt(:,:);
        phs_conc_dt = phs_conc_dt(:,:);
        dspha_dt = conv2(phs_spha_dt', kern, 'same')/tdwell; % units rad/{1,m,m^2,m^3}/time (spatial units depend on spatial order of term)
        dconc_dt = conv2(phs_conc_dt', kern, 'same')/tdwell;
        dspha_dt = dspha_dt';
        dconc_dt = dconc_dt';
        dspha_dt = reshape(dspha_dt, [Ncoeff_kspha skope_pts Nshots]);
        dconc_dt = reshape(dconc_dt, [Ncoeff_kconc skope_pts Nshots]);
        dspha_dt = permute(dspha_dt, [2 1 3]);
        dconc_dt = permute(dconc_dt, [2 1 3]);
    else
        dspha_dt = conv2(phs_spha(:,:), kern, 'same')/tdwell; % units rad/{1,m,m^2,m^3}/time (spatial units depend on spatial order of term)
        dconc_dt = conv2(phs_conc(:,:), kern, 'same')/tdwell;
    end

    % Peform course search
    resvec_min = inf;
    if numCoarseSearch > 0
        fprintf('coarse search (N=%d)', 2*numsamps+1);
        for courseTest = -numsamps:numsamps
            % Interpolate to match datatimes
            delSk = delSk0+courseTest*tdwell;
            %delSk = delSk0+courseTest*dtsamp;
            if Nshots == 1 || ~isEPI
                [phs_spha_a,phs_conc_a] = interpTrajTime(phs_spha,phs_conc,tdwell,delSk,datatime_interp);
            else 
                % Interpolate all the shots, considering delay shifts as
                % well if necessary
                for i = 1:Nshots
                    [phs_spha_a(:,:,i,:),phs_conc_a(:,:,i,:)] = interpTrajTime(phs_spha(:,:,i,:),phs_conc(:,:,i,:),tdwell,delSk+timeshift(i),datatime_interp);
                end    
            end

            if SMS_fact > 1 && doSMSwithkz && ~ispulseq
                phs = phs_spha_a(:,4)*(z0/1000);
                data_in = data_ref_SMS.*exp(1i*phs);
            end 

            if Nshots > 1 && ~is3D && ~ishighres
                phs_spha_interp = permute(phs_spha_a, [2 1 3]);
                phs_conc_interp = permute(phs_conc_a, [2 1 3]);
                S = sampHighOrder(b0map, datatime(:), phs_spha_interp(:,:), phs_conc_interp(:,:), phs_sph_grid, [], useGPU, useSingle, approach, opt);
            elseif ishighres
                phs_spha_interp = permute(phs_spha_a, [2 1 3]);
                phs_conc_interp = permute(phs_conc_a, [2 1 3]);
                opt.subFact = subFact;
                S = sampHighOrder(b0map, datatime(:), phs_spha_interp(:,:), phs_conc_interp(:,:), phs_sph_grid, [], useGPU, useSingle, encoding_approach, opt);
            elseif is3D
                phs_spha_interp = permute(phs_spha_a, [2 1 3 4]);
                phs_conc_interp = permute(phs_conc_a, [2 1 3 4]);
                opt.subFact = 20;
                S = sampHighOrder(b0map, datatime(:), phs_spha_interp(:,:), phs_conc_interp(:,:), phs_sph_grid, [], useGPU, useSingle, encoding_approach, opt);

                % Kcoeff_scale0p5 = rescale(phs_spha_interp(2:4,:,:,:),-0.5,0.5);
                % Kx = reshape(Kcoeff_scale0p5(1,:,:),[Nenc,Nsli]);
                % Ky = reshape(Kcoeff_scale0p5(2,:,:),[Nenc,Nsli]);
                % Kz = reshape(Kcoeff_scale0p5(3,:,:),[Nenc,Nsli]);
                % Nop = nufftOp([rec_Nx,rec_Ny,Nsli], [Kx(:),Ky(:),Kz(:)],[],false);
                %Nfun = @(x,transp) mrSampFunc(x,transp,Nop,R,opt); % mrSampFunc creates the net MRI sampling operation
            else
                % Update image
                S = sampHighOrder(b0map, datatime(:), phs_spha_a', phs_conc_a', phs_sph_grid, [], useGPU, useSingle, approach, opt);
            end
            opt.daNFull = size(data_in);
            opt.imNFull = [rec_Nx,rec_Ny,SMS_fact];
            opFunc = @(x,transp) mrSampFunc(x,transp,S,R,opt);
            % Use the same number of iterations for all recons
            opt_cgne_c.stopOnResInc = 0;
            
            [x,flag,relres,iter,resvec] = lsqr(opFunc,data_in(:),[],maxNit_cgne,[],[]);
            
            if doPhsConstraint
                x0 = reshape(x, opt.imNFull);
                clear x;
                kspace_firstit = fft2cmulti_norm(x0);
                sz = size(kspace_firstit,1);
    
                % Keep 40% of data
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
        
                im_firstit = ifft2cmulti_norm(kspace_firstit_trunc);
                phs_firstit = unwrap_phase(angle(im_firstit));
                phs_firstit = imgaussfilt(phs_firstit, 3.5);
    
                im_real = x0.*exp(-1i*phs_firstit);
    
                Crcvr_phi = Crcvr.*exp(1i*phs_firstit);
    
                R_phi = rcvrOp(Crcvr_phi,false);
                
                phsconst_Nit = 3;
                opFunc = @(x,transp) mrSampFunc(x,transp,S,R_phi,opt); % mrSampFunc creates the net MRI sampling operation
                %[x,FLAG,RELRES,ITER,RESVEC] = lsqr(opFunc,data_n(:),[],maxNit,[],[]);
                % Iteratively solve using the conjugate gradient method
                [~,flag,relres,iter,resvec] = lsqr(opFunc,data_in(:),[],phsconst_Nit,[],[],im_real(:));
            end
 
            if resvec(end)<resvec_min
                delSk_min = delSk;
                resvec_min = resvec(end);
            end
            fprintf('.');
        end
        fprintf(' %.2f; ', delSk_min)
        delSk = delSk_min;
    else
        delSk = delSk0;
    end
    
    % Get ready to start iterations
    ntry = 0;
    delChange = inf;
    delThresh = 0.005; % us
    fprintf('iterations...');
    delSk_perIt = zeros(1000,1);
    delSk_perIt(1) = delSk;
    while abs(delChange) > delThresh
        % Interpolate to match datatimes
        if Nshots == 1 || ~isEPI
            [phs_spha_a,phs_conc_a] = interpTrajTime(phs_spha,phs_conc,tdwell,delSk,datatime_interp);
            [dspha_dt_a,dconc_dt_a] = interpTrajTime(dspha_dt,dconc_dt,tdwell,delSk,datatime_interp);
        else
            % Interpolate all the shots, considering delay shifts as
            % well if necessary
            for i = 1:Nshots
                [phs_spha_a(:,:,i,:),phs_conc_a(:,:,i,:)] = interpTrajTime(phs_spha(:,:,i,:),phs_conc(:,:,i,:),tdwell,delSk+timeshift(i),datatime_interp);
                [dspha_dt_a(:,:,i,:),dconc_dt_a(:,:,i,:)] = interpTrajTime(dspha_dt(:,:,i,:),dconc_dt(:,:,i,:),tdwell,delSk+timeshift(i),datatime_interp);
            end    
        end

        if SMS_fact > 1 && doSMSwithkz && ~ispulseq
            phs = phs_spha_a(:,4)*(z0/1000);
            data_in = data_ref_SMS.*exp(1i*phs);
        end 

        if Nshots > 1 && ~is3D && ~ishighres
                phs_spha_interp = permute(phs_spha_a, [2 1 3]);
                phs_conc_interp = permute(phs_conc_a, [2 1 3]);
                dspha_dt_a = permute(dspha_dt_a, [2 1 3]);
                dconc_dt_a = permute(dconc_dt_a, [2 1 3]);
                S = sampHighOrder(b0map, datatime(:), phs_spha_interp(:,:), phs_conc_interp(:,:), phs_sph_grid, [], useGPU, useSingle, approach, opt);
        elseif ishighres
                phs_spha_interp = permute(phs_spha_a, [2 1 3]);
                phs_conc_interp = permute(phs_conc_a, [2 1 3]);
                dspha_dt_a = permute(dspha_dt_a, [2 1 3]);
                dconc_dt_a = permute(dconc_dt_a, [2 1 3]);
                opt.subFact = subFact;
                S = sampHighOrder(b0map, datatime(:), phs_spha_interp(:,:), phs_conc_interp(:,:), phs_sph_grid, [], useGPU, useSingle, encoding_approach, opt);
        elseif is3D 
                phs_spha_interp = permute(phs_spha_a, [2 1 3 4]);
                phs_conc_interp = permute(phs_conc_a, [2 1 3 4]);
                dspha_dt_a = permute(dspha_dt_a, [2 1 3 4]);
                dconc_dt_a = permute(dconc_dt_a, [2 1 3 4]);
                opt.subFact = 20;
                S = sampHighOrder(b0map, datatime(:), phs_spha_interp(:,:), phs_conc_interp(:,:), phs_sph_grid, [], useGPU, useSingle, encoding_approach, opt);

                % Kcoeff_scale0p5 = rescale(phs_spha_interp(2:4,:,:,:),-0.5,0.5);
                % Kx = reshape(Kcoeff_scale0p5(1,:,:),[Nenc,Nsli]);
                % Ky = reshape(Kcoeff_scale0p5(2,:,:),[Nenc,Nsli]);
                % Kz = reshape(Kcoeff_scale0p5(3,:,:),[Nenc,Nsli]);
                % Nop = nufftOp([rec_Nx,rec_Ny,Nsli], [Kx(:),Ky(:),Kz(:)],[],false);
                %Nfun = @(x,transp) mrSampFunc(x,transp,Nop,R,opt); % mrSampFunc creates the net MRI sampling operation
        else
            % Update image
            S = sampHighOrder(b0map, datatime(:), phs_spha_a', phs_conc_a', phs_sph_grid, [], useGPU, useSingle, approach, opt);
        end  
        
        if ntry>0
            % Use the same number of iterations for all recons
            opt_cgne.stopOnResInc = 0;
            maxNit_cgne = length(resvec)-1;
        end

        opFunc = @(x,transp) mrSampFunc(x,transp,S,R);
        [x, resvec] = cgne(opFunc,data_in,[],maxNit_cgne,opt_cgne);

        if doPhsConstraint
            x0 = x;
            clear x;
            kspace_firstit = fft2cmulti_norm(x0);
            sz = size(kspace_firstit,1);
    
            % Keep 40% of data
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
    
            im_firstit = ifft2cmulti_norm(kspace_firstit_trunc);
            phs_firstit = unwrap_phase(angle(im_firstit));
            phs_firstit = imgaussfilt(phs_firstit, 3.5);
    
            im_real = x0.*exp(-1i*phs_firstit);
    
            Crcvr_phi = Crcvr.*exp(1i*phs_firstit);
            
            R_phi = rcvrOp(Crcvr_phi,false);
            
            phsconst_Nit = 3;
            opFunc = @(x,transp) mrSampFunc(x,transp,S,R_phi); % mrSampFunc creates the net MRI sampling operation
            %[x,FLAG,RELRES,ITER,RESVEC] = lsqr(opFunc,data_n(:),[],maxNit,[],[]);
            % Iteratively solve using the conjugate gradient method
            [x, resvec] = cgne(opFunc,data_in,im_real,phsconst_Nit,opt_cgne);
        end
        
        % Update delay
        bhat = data_in - opFunc(x,'notransp');
        clear opFunc
        % Increase matrix segments for 3D acquisitions so they can fit in
        % memory
        if is3D
            opt.segmentedNdiv = 10;
            S = sampHighOrder(b0map, datatime(:), phs_spha_interp(:,:), phs_conc_interp(:,:), phs_sph_grid, [], useGPU, useSingle, 3, opt);
            % Revert to default number of segments
            opt.segmentedNdiv = 4;
        elseif ishighres
            S = sampHighOrder(b0map, datatime(:), phs_spha_interp(:,:), phs_conc_interp(:,:), phs_sph_grid, [], useGPU, useSingle, 3, opt);
        end    
        % Combine dimensions of shots/slices if necessary
        if Nshots > 1 || is3D
            S = S.setPhiDiv(dspha_dt_a(:,:),dconc_dt_a(:,:));
        else
            S = S.setPhiDiv(dspha_dt_a',dconc_dt_a');
        end    
        if doPhsConstraint
            opFunc = @(x,transp) mrSampFunc(x,transp,S,R_phi);  
        else
            opFunc = @(x,transp) mrSampFunc(x,transp,S,R);  
        end    
        bhathat = opFunc(x,'notransp');
        clear opFunc
        clear dspha_dt_a
        clear dconc_dt_a   

        delChange = bhathat(:)\bhat(:);
        delChange = delJumpFact*real(delChange);
        delSk = delSk + delChange;
        fprintf(' %.2f,', delSk)
        if (ntry>0) && (sign(delChangePrev) ~= sign(delChange)) && (delJumpFact>1)
            delJumpFact = delJumpFact/2;
        end
        delChangePrev = delChange;
        ntry = ntry+1;
        delSk_perIt(ntry+1) = delSk;
        if ntry > 50
            warning('did not converge after %d iterations\n', ntry)
            break;
        end
        if useGPU
            delSk = gather(delSk);
        end
    end
    delSk_perIt = delSk_perIt(1:ntry+1);
    fprintf(' done in %d sec\n', round(toc(tic1)));
end
