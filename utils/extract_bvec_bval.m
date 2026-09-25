clear

seq_fn = 'path/to/seqfile';
seq = mr.Sequence();
seq.read(seq_fn);

bval = seq.getDefinition('bval');
nImagingFrame = numel(bval); % how many volumes

% during acquisition, the b=0 scan is implemented as a diffuion scan, but
% with a very small b-value, equal to bvalue0_threshold
bvalue0_threshold = seq.getDefinition('bvalue0_threshold');
ind_b0 = bval==bvalue0_threshold;
bval(ind_b0) = 0;

bvec = seq.getDefinition('bvec');
bvec = reshape(bvec, [], 3);
bvec = bvec(end-nImagingFrame+1:end,:); % only the last nImagingFrame values are meaningful, the first few are directions of the dummy scan
bvec(ind_b0,:) = 0;
