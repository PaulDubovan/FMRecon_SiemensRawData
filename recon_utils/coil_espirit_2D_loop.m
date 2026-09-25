function [C,Mbig,Wbig]=coil_espirit_2D_loop(kspace,eigThresh_2)
    %Mbig,
    %big] = coil_espirit_2D_loop(kspace,eigThresh_2)



%parforsetup(5)
    
parfor kz=1:(size(kspace,3))
        kz
        k2D = squeeze(kspace(:,:,kz,:));
        [sx,sy,Nc] = size(k2D);
        ncalib = 24; % use 24 calibration lines to compute compression
        ksize = [6,6]; % kernel size

        % Threshold for picking singular vercors of the calibration matrix
        % (relative to largest singlular value.

        eigThresh_1 = 0.1;
         %eigThresh_1 = 0.01;

        % threshold of eigen vector decomposition in image space. 
        %eigThresh_2 = 0.6;
        % crop a calibration area
        calib = crop(k2D,[ncalib,ncalib,Nc]);

        %% Compute ESPIRiT EigenVectors
        % Here we perform calibration in k-space followed by an eigen-decomposition
        % in image space to produce the EigenMaps. 


        % compute Calibration matrix, perform 1st SVD and convert singular vectors
        % into k-space kernels
        tic
        [k,S] = dat2Kernel(calib,ksize);
 
        idx = max(find(S >= S(1)*eigThresh_1))
        [M,W] = kernelEig(k(:,:,:,1:idx),[sx,sy]);

        %%
        % crop sensitivity maps 
        maps = M(:,:,:,end).*repmat(W(:,:,end)>eigThresh_2,[1,1,Nc]);
        C(:,:,:,kz)=maps;
        toc
        fprintf('slice %d of %d\n',kz,size(kspace,3))
        Mlast = M(:,:,:,end);
        Wlast = W(:,:,end);
        Mbig(:,:,:,kz)=Mlast;
        Wbig(:,:,kz)=Wlast;
        
end



end

