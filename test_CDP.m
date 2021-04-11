clc;
clear;
close all;
%% [1] Simulate test data  (Coded diffraction pattern)
%%%% [1.1] read image
ImName =  './figures/2K/0058_color.png';
numMasks = 5; 
ImIn = imread(ImName);
image = double(rgb2gray(ImIn));
dims = size(image);
%%%% [1.2] generate masks 
b1 = [-1;1;-1i;1i];
b2 = [repmat(sqrt(0.5),4,1); sqrt(3)];
masks = b1(randi(4,[dims,numMasks])) .* b2(randi(5,[dims,numMasks])); 
%%%% [1.3] Make linear operators that act on a vectorized image (measurement operator)
A = @(x) CDP_Operator(x, masks, dims);
At = @(y) CDP_transposeOperator(y, masks, dims);
%%%% [1.4] generate intensity-only measurements and add noise
b0 = A(image(:));
b0 = abs(b0);
b0 = add_gaussion_noise(b0,10);
%%%% [1.5] parameters for retrieval
para.xt = image(:);
para.dims = dims;
para.mask = masks;
para.b0 = b0;
para.A = A;
para.At = At;
para.numMasks = numMasks;

%% [2] apply different algorithms for reconstruction
%%%% [2.1] Conventional algorithms  (Phasepack)
[rec_conv,algorithms] = runBenchmark2DImageRecovery(para);
%%%% [2.2] prDeep algorithm 
rec_prDeep = runprDeep(para);
%%%% [2.3]  LPR   algorithm
rec_LPR = LPR(para);
%% [3] Reconstruction quality and running time
num_algorithms  = length(rec_conv);
for i = 1 : num_algorithms
    rec_quality_PSNR(1,i) = psnr(rec_conv{1,i}./max(max(rec_conv{1,i})),image./max(image(:)));
    rec_quality_SSIM(1,i) = ssim(rec_conv{1,i}./max(max(rec_conv{1,i})),image./max(image(:)));
    fprintf('The reconstruction quality of %s is: PSNR--%4.2f, SSIM--%4.2f\n',algorithms{1,i},rec_quality_PSNR(1,i),rec_quality_SSIM(1,i))
end   
rec_quality_PSNR(1,num_algorithms+1) = psnr(rec_prDeep./max(max(rec_prDeep)),image./max(image(:)));
rec_quality_SSIM(1,num_algorithms+1) = ssim(rec_prDeep./max(max(rec_prDeep)),image./max(image(:)));
fprintf('The reconstruction quality of prDeep is: PSNR--%4.2f, SSIM--%4.2f\n',rec_quality_PSNR(1,num_algorithms+1),rec_quality_SSIM(1,num_algorithms+1))

rec_quality_PSNR(1,num_algorithms+2) = psnr(rec_LPR./max(max(rec_LPR)),image./max(image(:)));
rec_quality_SSIM(1,num_algorithms+2) = ssim(rec_LPR./max(max(rec_LPR)),image./max(image(:)));
fprintf('The reconstruction quality of LPR is: PSNR--%4.2f, SSIM--%4.2f\n',rec_quality_PSNR(1,num_algorithms+2),rec_quality_SSIM(1,num_algorithms+2))

%% [4] show results
subplot(4,4,1);imshow(image,[]);title('ground truth');
for j=1:num_algorithms    
    subplot(4,4,j+1);imshow(rec_conv{1,j},[]);title(algorithms{1,j});
end
subplot(4,4,num_algorithms+1);imshow(rec_prDeep,[]);title('prDeep');
subplot(4,4,num_algorithms+2);imshow(rec_LPR,[]);title('LPR');
