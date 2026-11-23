function win = filtNd(sz,alpha,dim,type)
  % win = filtNd(sz,alpha,dim, type)
  %
  % This filter is assumed to be applied in k-space. Accordingly, for size
  % N, the max value will occur at floor(N/2)+1, corresponding to center of
  % k-space after fftshift(fft()).
  %
  % alpha: alpha for type='kb', gaussian fwhm (wrt full width of kspace) for 'gauss'
  % No filtering along dims of k>dim
  %
  % (c) Corey Baron

  switch type
  case 'kb'
    filtFcn = @(szIn,alphaIn) kaisFunc(szIn,alphaIn);
  case 'gauss'
    filtFcn = @(szIn,alphaIn) gaussFunc(szIn,alphaIn);
  otherwise
    error('Unknown filter type')
  end

  if (dim > 1) && length(alpha)==1
    alpha = ones(dim,1)*alpha;
  end

  win = filtFcn(sz(1),alpha(1));
  win = win(:);
  for n=2:dim
    k2 = permute(reshape(filtFcn(sz(n),alpha(n)), [], 1), [2:n, 1]);
    win = repmat(win, [ones(1,n-1), sz(n)]) .* repmat(k2, [sz(1:n-1), 1]);
  end
  if dim<length(sz)
    win = repmat(win, [ones(1,dim), sz(dim+1:end)]);
  end

  win = win/max(win(:));

end

function win = kaisFunc(szIn,alphaIn)
    if ~mod(szIn,2)
        win = kaiser(szIn+1,alphaIn);
        win = win(1:end-1);
    else
        win = kaiser(szIn,alphaIn);
    end
end

function win = gaussFunc(szIn,alphaIn)
    if ~mod(szIn,2)
        win = normpdf(-szIn/2:szIn/2-1,0,szIn*alphaIn/2.355);
    else
        szh = floor(szIn/2);
        win = normpdf(-szh:szh,0,szIn*alphaIn/2.355);
    end
end

