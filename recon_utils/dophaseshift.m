function phsCorr = dophaseshift(phsRaw,isEPI)
% Detect abrupt phase deviations in NMR probe phase time courses and shift the phase to
% align with the preceding phase trend, as defined by the sliding window
% size

% Copyright (c) 2026 Paul Dubovan and Xingwang Yong
    if isEPI
        for nvol = 1:size(phsRaw,3)
            for nprobe = 1:size(phsRaw,2)
                phsCorr_a = phsRaw(:,nprobe,nvol);
                start_idx = 1001;
                nPrev = 5;
                nNext = 5;
                jumpFactor = 3;
                
                for k = start_idx:(length(phsCorr_a)-nNext)
        
                    % Previous trend
                    dPrev = diff(phsCorr_a(k-nPrev:k-1));
                    trend = median(abs(dPrev));
                
                    if trend == 0
                        continue
                    end
                
                    % Current difference
                    dCurr = phsCorr_a(k) - phsCorr_a(k-1);
                
                    % Candidate jump
                    if abs(dCurr) > jumpFactor*trend && abs(dCurr) > 2
                
                        % Look at following differences
                        dNext = diff(phsCorr_a(k:k+nNext));
                        jump_sign_curr = sign(dCurr);
                        jump_sign_next = sign(dNext(1));
        
                
                        % Are the following points back on trend?
                        future_on_trend = all(abs(dNext) < jumpFactor*trend);
                
                        if future_on_trend        
                            % Expected value based on previous trend
                            expected = phsCorr_a(k-1) + median(dPrev);
                
                            % Offset causing the jump
                            offset = phsCorr_a(k) - expected;
                
                            % Shift remainder of trajectory
                            phsCorr_a(k:end) = phsCorr_a(k:end) - offset;
                        elseif jump_sign_curr == jump_sign_next
                            % Expected value based on previous trend
                            expected = phsCorr_a(k-1) + median(dPrev);
                
                            % Offset causing the jump
                            offset = phsCorr_a(k) - expected;
                
                            % Shift remainder of trajectory
                            phsCorr_a(k:end) = phsCorr_a(k:end) - offset;
                        else
                
                            % Single-point corruption
                            phsCorr_a(k) = 0.5*(phsCorr_a(k-1) + phsCorr_a(k+1));
                
                        end
                    end
                end
                phsCorr(:,nprobe,nvol) = phsCorr_a;
            end
        end
    else
        for nvol = 1:size(phsRaw,3)
            for nprobe = 1:size(phsRaw,2)
                phsCorr_a = phsRaw(:,nprobe,nvol);
                start_idx = 1001;
                nPrev = 5;
                nNext = 5;
                thresh = 0.5;
                
                for k = start_idx:(length(phsCorr_a)-nNext)
        
                    % Previous trend
                    dPrev = diff(phsCorr_a(k-nPrev:k-1));
                    trend = median(dPrev);
                
                    if trend == 0
                        continue
                    end
                
                    % Current difference
                    dCurr = phsCorr_a(k) - phsCorr_a(k-1);
                
                    % Candidate jump
                    if sign(dCurr) ~= sign(trend) && abs(dCurr) > thresh
                   
                            % Expected value based on previous trend
                            expected = phsCorr_a(k-1) + median(dPrev);
                
                            % Offset causing the jump
                            offset = phsCorr_a(k) - expected;
                
                            % Shift remainder of trajectory
                            phsCorr_a(k:end) = phsCorr_a(k:end) - offset;
                    end
                end
                phsCorr(:,nprobe,nvol) = phsCorr_a;
            end
        end
    end    
end