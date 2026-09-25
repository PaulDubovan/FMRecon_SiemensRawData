function res = my_apply_highOrderSamp(in, tflag, params )
% apply sense(actually it is extended encoding model) for each shot, i.e. each shot gives an image

% Use persistent variable to cache operators, because it would take a very
% long time to construct the high order sampling operators
persistent cached_operators cached_params

% Copyright (c) 2026 Paul Dubovan and Xingwang Yong

% extract input
Nshots        = params.Nshots;
B0map         = params.B0map;
datatime = params.datatime;
kspha_n       = params.kspha_n;
kconc_n       = params.kconc_n;
grid          = params.grid;
useGPU        = params.useGPU;
opt           = params.opt;
R             = params.R;       % sensitivity map
Nenc          = params.Nenc; % number of samples in each shot
isEPI         = params.isEPI;
is3D          = params.is3D;
useSingle     = params.useSingle;
SMS_fact      = params.SMS_fact;
subFact       = params.subFact;
encoding_approach = params.encoding_approach;

sz = size(R.maps);
numCha = sz(end);

% Check if we can reuse cached operators by comparing params directly
if ~isempty(cached_operators) && ~isempty(cached_params) && isequal(params, cached_params)
    % Use cached operators
    operators = cached_operators;
    % fprintf('Using cached operators\n');
else
    % Calculate new operators and cache them
    fprintf('Calculating new operators and caching them\n');
    
    % reshape
    datatime = reshape(datatime,Nenc,Nshots);
    kspha_n_all_shots = reshape(kspha_n,[],Nenc,Nshots);
    kconc_n_all_shots = reshape(kconc_n,[],Nenc,Nshots);
    
    % construct the operators
    operators = cell(1,Nshots);
    for iter_shot = 1:Nshots
        tic
        kspha_n = kspha_n_all_shots(:,:,iter_shot);
        kconc_n = kconc_n_all_shots(:,:,iter_shot);
        datatime_shot = datatime(:,iter_shot);
        
        % Build sampling operator
        nslice = 1;
        if SMS_fact > 1 && ~is3D
            opt.subFact = subFact;
            S = sampHighOrder(B0map,datatime_shot,kspha_n,kconc_n,grid,[],useGPU,useSingle,encoding_approach,opt);
        elseif is3D
            useSingle = 1;
            useGPU = 0;
            opt.subFact = 20;
            S = sampHighOrder(B0map(:,:,nslice),datatime_shot,kspha_n,kconc_n,grid,[],useGPU,useSingle,encoding_approach,opt);
        elseif opt.ishighres
            opt.subFact = subFact;
            S = sampHighOrder(B0map(:,:,nslice),datatime_shot,kspha_n,kconc_n,grid,[],useGPU,useSingle,encoding_approach,opt);
        else
            S = sampHighOrder(B0map(:,:,nslice),datatime_shot,kspha_n,kconc_n,grid,[],useGPU,[],encoding_approach,opt);
        end
        
        operators{iter_shot} = S;
        
        fprintf('Constructing operator for shot %d, took %.1f seconds\n', iter_shot, toc);
    end
    
    % Cache the operators and params
    cached_operators = operators;
    cached_params = params;
    
end


if strcmp(tflag,'transp')
    % Transposed SENSE operator:
    in = reshape(in,Nenc,Nshots,numCha);
    imgs = zeros(   opt.imNFull(2)*opt.imNFull(2),   Nshots);
    
    for iter_shot = 1:Nshots
        sampFunc = @(x) mrSampFunc(x,'transp',operators{iter_shot},R,opt);
        imgs(:,iter_shot) = sampFunc(in(:,iter_shot,:));
    end
    
    res = imgs(:);
elseif strcmp(tflag,'notransp')
    % Forward SENSE operator:
    in = reshape(in,opt.imNFull(1),opt.imNFull(2),Nshots);
    kspace = zeros(Nenc*numCha,Nshots);
    
    for iter_shot = 1:Nshots
        sampFunc = @(x) mrSampFunc(x,'notransp',operators{iter_shot},R,opt);
        kspace(:,iter_shot) = sampFunc(in(:,:,iter_shot));
    end
    
    kspace = reshape(kspace,[Nenc,numCha,Nshots]);
    kspace = permute(kspace,[1 3 2]);
    res = kspace(:);
else
    error('Syntax error, un-recognized 2nd param: %s', tflag);
end


end % of function