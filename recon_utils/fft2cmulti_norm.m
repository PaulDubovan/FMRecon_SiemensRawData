function [res] = fft2cmulti_norm(x)
fctr = size(x,1)*size(x,2);
res = (1/sqrt(fctr)) *fftshift(fft2(ifftshift(x)));
end

