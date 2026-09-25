function B0map = calc_B0_map_spm(phase_name,mag_name,TEs)
% Template script for FieldMap toolbox
% _____________________________________________________________________
%
% This script gives an example of how to run the FieldMap toolbox
% FieldMap.m without using the GUI. It can be expanded using standard 
% matlab code to for example create multiple field maps or use a single
% fieldmap to unwarp multiple images.
% 
% As it stands the script uses routines from the FieldMap toolbox to 
% create a single field map which is matched to an EPI and then used
% to unwarp it. A structural image is loaded and matched to the unwarped 
% EPI.
%
% For details about the FieldMap toolbox, see FieldMap.man. For a 
% description of the components of the structure IP, see FieldMap.m.
% For an introduction to the theoretcial and practical principles behind 
% the toolbox, see principles.man.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson and Chloe Hutton
% $Id: FieldMap_ngui.m 1358 2008-04-10 11:20:26Z guillaume $

%----------------------------------------------------------------------  
% Set up default parameters and structures 
%----------------------------------------------------------------------
if nargin<3
    TEs = [4.92,7.38];
end
fnames = {phase_name,mag_name};
try
    fnames = convertContainedStringsToChars(fnames);
catch

end
spm('defaults','FMRI');
% addpath(fullfile(get_syncthing_dir,'third_party_repos/spm12-git/toolbox/FieldMap'));
IP = FieldMap('Initialise'); % Gets default params from pm_defaults
IP.et{1} = TEs(1);
IP.et{2} = TEs(2);

%----------------------------------------------------------------------
% Load measured field map data - phase and magnitude or real and imaginary
%----------------------------------------------------------------------

if IP.uflags.iformat=='PM' 
   % for index=1:4
   for index=1:2
       IP.P{index} = spm_vol(fnames{index});
       if index==1
           % doscl=spm_input('Scale this to radians?',1,'b','Yes|No',[1,0]);
           % Finter=spm_figure('FindWin','Interactive');
           % close(Finter);
           if 1 % doscl==1
               % msgbox([mfilename('fullpath'),'.m: The phase image is scaled and writen to disk'])
               tmp=FieldMap('Scale',IP.P{index}.fname);
               IP.P{index} = spm_vol(tmp.fname);
           end
       end
   end
else
    error('Undefined')
   % for index=1:4
   %    IP.P{index} = FieldMap('LoadFile',index);
   % end
end

%----------------------------------------------------------------------
% Or you may want to load a precalculated Hz phase map instead...
%----------------------------------------------------------------------

% [IP.fm, IP.pP] = FieldMap('LoadFieldMap');

%----------------------------------------------------------------------
% Create field map (in Hz) - this routine calls the unwrapping
%----------------------------------------------------------------------

IP.fm = FieldMap('CreateFieldMap',IP);

B0map = IP.fm.fpm;

% %----------------------------------------------------------------------
% % Write out field map
% % Outputs -> fpm_NAME-OF-FIRST-INPUT-IMAGE.img
% %----------------------------------------------------------------------
% 
% FieldMap('Write',IP.P{1},IP.fm.fpm,'fpm_',64,'Smoothed phase map');
% 
% %----------------------------------------------------------------------
% % Convert Hz to voxels and write voxel displacement map 
% % Outputs -> vdm_NAME-OF-FIRST-INPUT-IMAGE.img
% %----------------------------------------------------------------------
% 
% [IP.vdm, IP.vdmP]=FieldMap('FM2VDM',IP);
% 
% %----------------------------------------------------------------------
% % Select an EPI to unwarp
% %----------------------------------------------------------------------
% 
% IP.epiP = FieldMap('LoadEPI');
% 
% %----------------------------------------------------------------------
% % Match voxel displacement map to image
% % Outputs -> mag_NAME-OF-FIRST-INPUT-IMAGE.img
% %----------------------------------------------------------------------
% 
% IP.vdmP = FieldMap('MatchVDM',IP);
% 
% %----------------------------------------------------------------------
% % Unwarp EPI
% %----------------------------------------------------------------------
% 
% IP.uepiP = FieldMap('UnwarpEPI',IP);
% 
% %----------------------------------------------------------------------
% % Write unwarped EPI 
% % Outputs -> uNAME-OF-EPI.img
% %----------------------------------------------------------------------
% 
% IP.uepiP = FieldMap('Write',IP.epiP,IP.uepiP.dat,'u',IP.epiP.dim(4),'Unwarped image');
% 
% %----------------------------------------------------------------------
% % Load a structural image
% %----------------------------------------------------------------------
% 
% IP.nwarp=FieldMap('LoadStructural');
% 
% %----------------------------------------------------------------------
% % Coregister structural with the unwarped image
% %----------------------------------------------------------------------
% 
% FieldMap('MatchStructural',IP);
% 
% %______________________________________________________________________
end