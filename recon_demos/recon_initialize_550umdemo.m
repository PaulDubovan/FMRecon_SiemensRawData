%% Data initialization and reconstruction settings

% Add dependencies 
addpath(genpath('/path/to/FMRecon_SiemensRawData'))
addpath(genpath('/path/to/FMRecon_SiemensRawData/recon_utils')) % <-- Ensure that
%this is set as the last path so that updated files in recon_utils get
%called, and not the ones in the Dependencies folder/repositories
%%
%Raw imaging and computed Skope data can be downloaded from: https://huggingface.co/datasets/pdubovan3/550um_Connectome2_data4recon/tree/main
%% EPI 18 dir
% Path to imaging and skope raw data
filepath = fullfile('/path/to/datarepository');
dataPath = fullfile(filepath, '/dataload_EPI_550um');
% Path to GRE raw data
grePath = filepath;
% Path to seq file
seq_file_path = fullfile(filepath, '/seqfiles/epi_550um_6shots_18dir.seq');
% Path to B0 map
 spm_B0mapPath = fullfile(filepath, 'B0map_spm.mat');

% Calculated delay value for syncrhonization of skope trajectory and raw
% data
del = 567.03;
% Conjugate gradient iterations (no S-LORAKS)
maxNit = 5;

% opt.fitOrder = 1
% opt.retro_undersample_factor = 2;
%% spiral 18 dir
filepath = fullfile('/path/to/datarepository');
dataPath = fullfile(filepath, '/dataload_spiral_550um');
grePath = ffilepath;
seq_file_path = fullfile(filepath, '/seqfiles/spi_550um_9shots_18dir.seq');
spm_B0mapPath = fullfile(filepath, 'B0map_spm.mat');

del = 84.79;
maxNit = 10;

% opt.fitOrder = 2
% opt.retro_undersample_factor = 3;
%%
opt = [];
opt.fitOrder = 1;
opt.fitDoSteps = false;
opt.fitDoWeights = false;
opt.fitDoResAdj = false;
opt.ispulseq = 1;
opt.delSk = del;
opt.maxNit = maxNit;        
opt.numCoarseSearch = [];
opt.useSingle = 1;
opt.ishighres = 1;
opt.use_spm_B0map = 1;
opt.dospiraldelay = 0;
opt.useGPU = 0;
opt.retro_undersample_factor = 2;
opt.encoding_approach = 1;
opt.dob0reg = 1;
opt.doGibbsCor = 0;
opt.b0_indices = [1,11,21];
opt.input_vol = [1:21];
opt.input_sli = [1:15];
opt.input_vol = [1];
opt.input_sli = [4];

%% EPI
% b0 parameters
opt.b0_rank = 50;
opt.b0_lambda = 0.005;

% dwi parameters
opt.dwi_rank = 50;
opt.dwi_lambda = 0.002;

%% spiral
% b0 parameters
opt.b0_rank = 40;
opt.b0_lambda = 0.008;

% dwi parameters
opt.dwi_rank = 40;
opt.dwi_lambda = 0.005;

%%
im_recon = FMrecon_demo(dataPath,grePath,opt,seq_file_path,spm_B0mapPath);

