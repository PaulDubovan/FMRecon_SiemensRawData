## Data organization pipeline for expanded encoding reconstruction using the matMRI toolbox 
### Requires raw scanner twix data, dual echo GRE for coil sensitivity estimation and field map generation (optional), skope measurements, POET simulation files
### Compatible with Siemens ep2d_diff and ep3d_diff sequences, and provided spiral/EPI diffusion Pulseq sequences (.seq files)

### Repositories to include for Pulseq sequence development:
Pulseq: https://github.com/pulseq/pulseq.git  \
SafePNS: https://github.com/filip-szczepankiewicz/safe_pns_prediction  \
VDspiral: https://web.archive.org/web/20220617183021/http://mrsrl.stanford.edu/~brian/vdspiral/  \
SONIQS: https://github.com/mriphysics/SONIQS

### Repositories to include for reconstruction:
Pulseq: https://github.com/pulseq/pulseq.git  \
MatMRI toolbox: https://gitlab.com/cfmm/matlab/matmri.git  \
spm12: https://github.com/spm/spm12.git \
BUDA_SLORAKS: https://github.com/CongyuLiao/BUDA_SLORAKS    \
ESPIRiT: https://people.eecs.berkeley.edu/~mlustig/Software.html    \
mapVBVD: https://github.com/pehses/mapVBVD  \
Tensor-MP-PCA: https://github.com/Neurophysics-CFIN/Tensor-MP-PCA   \
unwrap_phase: https://github.com/mfkasim1/unwrap_phase   \
FieldMonitoring&CS: https://gitlab.com/cfmm/datasets/ss-spiral-dmri-with-fm

### Where to start  
Start by initializing recon_initialize_550um.m. This gives an example of how to initialize the necessary file paths, recon parameters, and calls the reconstruction script. 

### Reconstruction script variations
FMrecon_acquisitionsystemv1: General script using Skope data compatible with the first-generation acquisition system  \
FMrecon_NYOX: General script using Skope data compatible with the NYOX system   \
FMrecon_demo: Demo script for exampled inputted raw imaging and Skope data  \
FMrecon_nominal: General script for inputted nominal k-space trajectories (Pulseq) 

### Pulseq sequence generation examples included
epi_diffusion_general.m
spiral_diffusion_general.m

Disclaimer
The sequences are not medical device applications. The sequences and any resulting data are to be used for research purposes only and not to diagnose, treat, cure or prevent any disease. The distributors make no representation and assume no responsibility for the misuse of these sequences. The distributors do not take any responsibility or liability for damages, harm, loss of data or similar incidents that are in direct or indirect relation to the usage or presence of any of its devices. MR sequences or particular implementations thereof can be subject to intellectual property of one or several parties.