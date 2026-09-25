function r= fov_reduced(Im)
    sizeim = size(Im);
    im = reshape(Im,sizeim(1),[]);
    im = im(1+end/4:3*end/4,:);
    sizeim1 = size(im,1);
    r = reshape(im,[sizeim1,sizeim(2:end)]);
end
    

