function [res] = ifft2cmulti_norm(x)
    fctr = size(x,1)*size(x,2);
    res= sqrt(fctr)*fftshift(ifft2(ifftshift(x)));
end

