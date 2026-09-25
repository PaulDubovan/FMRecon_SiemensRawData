function r = my_prewhiten(rawdata, noise_matrix)
% MY_PREWHITEN Applies noise pre-whitening to MRI data.
%
%   r = my_prewhiten(rawdata, noise_matrix)
%
%   Inputs:
%       rawdata:      Multi-dimensional array where the LAST dimension 
%                     corresponds to the coils (e.g., [Nx, Ny, Ncoils]).
%       noise_matrix: Noise Covariance Matrix (Psi) of size [Ncoils, Ncoils].
%
%   Outputs:
%       r:            Whitened data with the same dimensions as rawdata.
%                     The noise covariance of 'r' will be the Identity matrix.

    sz = size(rawdata);
    nc = sz(end); % Number of coils in the raw data
    nc_noise = size(noise_matrix, 1);
    
    % Verify dimension consistency
    assert(nc == nc_noise, 'Error: Coil count mismatch between data and noise matrix.');
    
    % 1. Reshape and Permute Data
    % Reshape rawdata into [N_points, N_coils] first, then permute to 
    % [N_coils, N_points]. 
    % We place the coil dimension first to enable efficient matrix operations.
    tmp = permute(reshape(rawdata, [], nc), [2, 1]);
    
    % 2. Cholesky Decomposition
    % Decompose the Noise Covariance Matrix (Psi) such that Psi = L * L'
    % 'lower' returns the lower triangular matrix L.
    L = chol(noise_matrix, 'lower'); 
    
    % 3. Apply Whitening
    % We need to compute r = L^-1 * tmp.
    % In MATLAB, the backslash operator (\) solves the linear system L * r = tmp.
    % Since L is lower triangular, MATLAB automatically uses Forward Substitution,
    % which is numerically stable and faster than explicit inversion (inv).
    whitened_tmp = L \ tmp;
    
    % 4. Restore Dimensions
    % Permute back to [N_points, N_coils]
    r = permute(whitened_tmp, [2 1]);
    
    % Reshape to the original input size (e.g., [Nx, Ny, Ncoils])
    r = reshape(r, sz);

end