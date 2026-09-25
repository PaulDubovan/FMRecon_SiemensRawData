function [TE,big_delta,small_delta,delayTE1,delayTE2] = calc_minTE_spiral(bvalue, grad_max, slew_max, gradRasterTime)
% Minimal TE of spin-echo with monopolar diff grad using spiral-out readout
% 
% Input
%   b-value                  - unit is s/mm^2
%   grad_max                 - unit is mT/m
%   slew_max                 - unit is mT/m/ms
%   gradRasterTime(optional) - unit is us, default is 10us
% 
% Output
%   TE          - unit is second
%   big_delta   - unit is second
%   small_delta - unit is second
% 
% Reference
%   b-value formula from handbook, page 287, example 9.2

% Xingwang Yong, 20231003

% TODO, draw seq timing for better understanding


% Time :    0          TE/2             TE 
% Event:  RF90        RF180       start of spiral-out readout

%% parse input
if ~exist('gradRasterTime','var')
    gradRasterTime = 10e-6;
end


%% check pulseq installation
isHasPulseq = true;
try
    mr.Sequence();
catch
    isHasPulseq = false;
    warning('It is suggested to use pulseq for more accurate timing: https://github.com/pulseq/pulseq')
end


%% get duration of 90 and 180, including corresponding gradient
if ~isHasPulseq % this is inaccurate when Gmax is very high
    grad_dur_180 = 6.3e-3; % whole duration of slice selection gradient
    grad_dur_90  = 4e-3; % from center of RF to the end of slice rephasing gradient
else
     sys = mr.opts('MaxGrad',grad_max,'GradUnit','mT/m',...
    'MaxSlew',slew_max,'SlewUnit','T/m/s','gradRasterTime',gradRasterTime,...
    'rfRingdownTime', 30e-6, 'rfDeadtime', 100e-6, 'adcDeadTime', 10e-6);

    flipAngle = 90;
    tRFex = 3e-3;
    tRFref = 3e-3;
    sliceThickness=3e-3;

    tExwd=tRFex+sys.rfRingdownTime+sys.rfDeadTime;
    tRefwd=tRFref+sys.rfRingdownTime+sys.rfDeadTime;

    [rf90, gs_ex, gzReph] = mr.makeSincPulse(flipAngle*pi/180,'system',sys,'Duration',tRFex,...
        'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4);
    GSex = mr.makeTrapezoid('z',sys,'amplitude',gs_ex.amplitude,'FlatTime',tExwd);
        
    rfCenterInclDelay=rf90.delay + mr.calcRfCenter(rf90);
    grad_dur_90 = (mr.calcDuration(gs_ex)-rfCenterInclDelay)   +  mr.calcDuration(gzReph);

    % Create 180 degree slice refocusing pulse and gradients
    [rf180, gz180] = mr.makeSincPulse(pi,'system',sys,'Duration',tRFref,...
        'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'use','refocusing');
    GSref = mr.makeTrapezoid('z',sys,'amplitude',GSex.amplitude,'FlatTime',tRefwd);
   
    grad_dur_180 =  mr.calcDuration(GSref);
end


%% bi-section
readout_te = 0;
% rup = 1e-3*grad_max/slew_max;
% rup = 1e-3*ceil(grad_max/slew_max/gradRasterTime)*gradRasterTime;
big_delta_formula = @(delta,rup) ( grad_dur_180 + delta + rup + readout_te ); % assuming diff gradients symmetric around 180
roundup_to_gradraster = @(t) ceil(t/gradRasterTime)*gradRasterTime;

flattime_max = 100; % unit is second
flattime_min = 0;
max_iter = 100;
for iter = 1:max_iter
    
    flattime = 0.5*(flattime_min+flattime_max);
    flattime = roundup_to_gradraster(flattime);

    fun = @(grad)  -bvalue + tmp_calc_bval_func(grad,slew_max,flattime,big_delta_formula,gradRasterTime);
    range = [0,100*grad_max]; % ensure we always get a solution
    g = fzero(fun,range);

    if g>grad_max % flat time too small
        flattime_min = flattime;
    else
        flattime_max = flattime;
    end

    if iter > max_iter || ...
            abs(flattime_max-flattime_min)<=gradRasterTime % should be accurate enough
        break;
    end

end

grad = g;
bvalue_real = tmp_calc_bval_func(grad,slew_max,flattime,big_delta_formula,gradRasterTime);
rup = 1e-3*ceil(grad/slew_max/gradRasterTime)*gradRasterTime; % unit is second
small_delta = flattime + rup;
big_delta = big_delta_formula(small_delta,rup);




%%
% % 
% bval = calc_bval(grad_max,small_delta,big_delta);
% rel_err = norm(bval-bvalue)/norm(bvalue);
% if rel_err > 0.01
%     error('Re-calculated bvalue is %.2f, deviating from given bvalue %.2f', bval, bvalue);
% end

half_TE = grad_dur_90 + (big_delta-0.5*grad_dur_180);
% for spiral the left side of 180 is more compact, while for epi the right 
% side is more compact
TE = half_TE*2;

delayTE1=  TE/2   -   grad_dur_90  -  0.5*grad_dur_180;
delayTE1 = ceil(delayTE1/gradRasterTime)*gradRasterTime;
delayTE2=  TE/2   -  0.5*grad_dur_180;
delayTE2 = ceil(delayTE2/gradRasterTime)*gradRasterTime;
assert( (small_delta+rup)<=delayTE1 );
assert( (small_delta+rup)<=delayTE2 );
end




%% helper functions
function bvalue = tmp_calc_bval_func(grad,slew,flattime,big_delta_formula,gradRasterTime)
rup = 1e-3*ceil(grad/slew/gradRasterTime)*gradRasterTime; % unit is second
small_delta = flattime + rup;
big_delta = big_delta_formula(small_delta,rup);
bvalue = calc_bval_trap(grad,small_delta,big_delta,rup);
end

function bvalue = calc_bval_trap(grad,small_delta,big_delta,rup)
% grad unit is mT/m, time unit is seconds
% handbook, page 287, example 9.2
gamma = 2*pi*42.57e3;           % rad/mT
bvalue = (gamma.*grad/1000).^2*( small_delta^2*(big_delta-small_delta/3) + rup^3/30-small_delta*rup^2/6  );
end