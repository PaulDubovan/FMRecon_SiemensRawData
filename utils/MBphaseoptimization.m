% Inputs:
% rf90.signal, rf180.signal = original single-band RF pulses
% Bandsep = your three band positions
% tb = 4

tb = 4;
nStarts = 100;

rf90in = rf90.signal(:);
rf180in = rf180.signal(:);

% Optimize relative phases; phi1 is fixed at zero
objective = @(x) max([ ...
    makePeakRF(rf90in, Bandsep, tb, [0 x(1) x(2)]), ...
    makePeakRF(rf180in, Bandsep, tb, [0 x(1) x(2)])]);

opts = optimset('Display', 'off', ...
                'MaxIter', 3000, ...
                'MaxFunEvals', 10000, ...
                'TolX', 1e-8, ...
                'TolFun', 1e-8);

bestPeak = Inf;
bestX = [0 0];

for k = 1:nStarts
    if k == 1
        x0 = [0 0];  % Baseline
    else
        x0 = 2*pi*rand(1,2);
    end

    [x, fval] = fminsearch(objective, x0, opts);

    if fval < bestPeak
        bestPeak = fval;
        bestX = x;
    end
end

phs = mod([0 bestX], 2*pi);

fprintf('Optimized phases (rad):\n');
disp(phs);

fprintf('Optimized peak RF amplitude: %.6g\n', bestPeak);

% Baseline: all phases zero
baselinePeak = objective([0 0]);
fprintf('Baseline peak RF amplitude: %.6g\n', baselinePeak);
fprintf('Peak amplitude reduction: %.2f%%\n', ...
    100*(1-bestPeak/baselinePeak));

function peak = makePeakRF(pulse, Bandsep, tb, phs)
    n = numel(pulse);
    t = (-n/2:n/2-1)/n;

    p = zeros(1,n);
    for ii = 1:numel(Bandsep)
        p = p + exp(1i*2*pi*Bandsep(ii)*tb*t ...
                    + 1i*phs(ii));
    end

    peak = max(abs(pulse .* p(:)));
end