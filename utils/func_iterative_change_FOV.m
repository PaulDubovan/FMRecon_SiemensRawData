function Fcoeff = func_iterative_change_FOV(desired_duration,smax, gmax,T,Ninterleaves,fov,res_mm,gradRasterTime)
% change FOV (sampling density) to achieve the desired duration
% 
%   slew unit is T/m/s
%   grad unit is mT/m

% Xingwang Yong
% TDODO, write doc


rmax = 5/res_mm;
underSampling_factor = round( gradRasterTime/T );

Fcoeff = fov*100*[1];
[k,g,s,time,r,theta] = vds(smax*100,gmax/10,T,Ninterleaves,Fcoeff,rmax);
k = k(1:underSampling_factor:end);
g = g(1:underSampling_factor:end);
s = s(1:underSampling_factor:end);
calculated_ro_dur = numel(g) * gradRasterTime;

if desired_duration==calculated_ro_dur
    return;
end

if desired_duration < calculated_ro_dur % decrease FOV to shorten time
    upper = 1;
    lower = 0;
    max_iter = 100;
    iter = 0;
    s_result = 0;
    isFound = false;
    while iter<max_iter
        iter = iter + 1;
        middle = (upper+lower)/2;

        Fcoeff = fov*100*[1, -middle];
        [k,g,s,time,r,theta] = vds(smax*100,gmax/10,T,Ninterleaves,Fcoeff,rmax);
        g = g(1:underSampling_factor:end);
        calculated_ro_dur = numel(g) * gradRasterTime;

        dur_diff = desired_duration - calculated_ro_dur;
        if dur_diff > 0 && dur_diff<=0.1e-3
            s_result = middle;
            isFound = true;
            break;
        end

        if calculated_ro_dur > desired_duration
            lower = middle;
        else
            upper = middle;
        end
    end
    assert(isFound,'Fail to iteratively decrease FOV')


    Fcoeff = fov*100*[1, -s_result];
else
    upper = 100;
    lower = 1;
    max_iter = 100;
    iter = 0;
    s_result = 0;
    isFound = false;

    while iter<max_iter
        iter = iter + 1;
        middle = (upper+lower)/2;

        Fcoeff = fov*100*[middle, 1-middle];
        [k,g,s,time,r,theta] = vds(smax*100,gmax/10,T,Ninterleaves,Fcoeff,rmax);
        g = g(1:underSampling_factor:end);
        calculated_ro_dur = numel(g) * gradRasterTime;

        dur_diff = desired_duration - calculated_ro_dur;
        if dur_diff > 0 && dur_diff<=0.1e-3
            s_result = middle;
            isFound = true;
            break;
        end

        if calculated_ro_dur > desired_duration
            upper = middle;
        else
            lower = middle;
        end
    end
    assert(isFound,'Fail to iteratively increase FOV')


    Fcoeff = fov*100*[s_result, 1-s_result];
end
end