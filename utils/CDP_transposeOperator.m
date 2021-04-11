function x = CDP_transposeOperator(y, masks, dims)
    n1 = dims(1);
    n2 = dims(2);
    L = size(masks,3);              % get number of masks
    y = reshape(y, [n1,n2,L]);   % image comes in as a vector.  Reshape to rectangle
    
    x = n1*n2*ifft2(y).*conj(masks);
    x = sum(x,3);
    x = x(:);
end