function y = CDP_Operator(x, masks, dims)
    x = reshape(x, dims);   % image comes in as a vector.  Reshape to rectangle
    L = size(masks,3);              % get number of masks
    % Compute measurements
    copies = repmat(x,[1,1,L]);
    y = fft2(masks.*copies);
    y = y(:);
end