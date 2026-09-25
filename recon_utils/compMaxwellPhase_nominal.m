function phiTerms = compMaxwellPhase_nominal(linPhase,dt,gam,coilParams,B0)
% phiConcTerms = compMaxwellPhase(linPhase,dt,gam,coilParams,B0)
% Compute phase from 2nd order concomitant grads based on linear phase terms. 
% All units should be SI 
%
% (c) Corey Baron, 2021

% Find Gx, Gy, Gz using derivative of gaussian filter
G = computeGradFromPhase_nominal(linPhase,dt,gam);

% Find conc grad prefactors (i.e., phase divided by spatial terms)
maxOut = computeMaxwellTerms(2,coilParams,G(:,1,:,:),G(:,2,:,:),G(:,3,:,:),B0);

% final output of phi. This will need to encorporate nonLinSphHarm.
% Probably the best way to incorporate the gradient nonlinearity is to
% include it only in the basisFuncSphHarm and basisFuncConcGrad eqns...
% NB: this currently assumes symmetry in x and y for the gradient
% if coilParams.alph ~= 0.5
%     error('coilParams.alph /= 0.5 not yet accounted for')
% end   
% commented out PD 110625 since it should be able to support /= 0.5 alpha
% values
phiTerms = zeros([size(G,1) 4 size(G,3) size(G,4)], 'like', G);
phiTerms(:,1,:,:) = 2*pi*(gam*dt)*cumsum(maxOut.z2,1);
phiTerms(:,2,:,:) = 2*pi*(gam*dt)*cumsum(maxOut.x2+maxOut.y2,1); % Added the y^2 component as well - PD 110625
phiTerms(:,3,:,:) = 2*pi*(gam*dt)*cumsum(maxOut.zx,1); 
phiTerms(:,4,:,:) = 2*pi*(gam*dt)*cumsum(maxOut.zy,1); 


end

