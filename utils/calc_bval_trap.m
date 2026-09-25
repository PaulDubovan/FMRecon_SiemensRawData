function bvalue = calc_bval_trap(grad_max,small_delta,big_delta,rup)
% mT/m, s
% handbook, page 287, example 9.2
gamma = 2*pi*42.57e3;           % rad/mT
bvalue = (gamma.*grad_max/1000).^2*( small_delta^2*(big_delta-small_delta/3) + rup^3/30-small_delta*rup^2/6  ); 
end