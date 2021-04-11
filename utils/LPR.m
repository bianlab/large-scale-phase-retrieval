function rec_LPR = LPR(opt)
fprintf('  LPR begin...');
tic;
%% [1]set param
lambda   =    1; % regularization coefficiency  
load('FFDNet_gray');
opt.net = vl_simplenn_tidy(net);
opt.useGPU = true;
if opt.useGPU
   opt.net = vl_simplenn_move(opt.net, 'gpu') ;
end
opt.ffdnetvnorm = true;      
dims = opt.dims;
% xt = imresize(opt.xt,dims);
% xt = xt./max(xt(:));
y = opt.b0;
sigma = [0.0784];
maxiter = [5];

A  = opt.A;
At = opt.At;
AP_solver = struct('algorithm','gerchbergsaxton');  %based on Phasepack                                   
algorithm = {AP_solver};
% AP solver parameters (for initialization)
opts.maxInnerIters = 10;
opts.algorithm = algorithm{1};
opts.maxIters = 8;
opts.maxTime = 100;
opts.tol = 1e-4;
opts.verbose = 0;
opts.recordTimes = true;
opts.recordMeasurementErrors = false;
opts.recordReconErrors = false; 
opts.recordResiduals = true;
opts.xt = opt.xt(:);

%% [2]LPR 
%%% [2.1] initialization
[sol] = solveGerchbergSaxton(A, At,opt.b0(:), ones(dims(1)*dims(2),1), opts);
v0 = abs(reshape(sol,dims));
v0 = v0./(max(max(v0)));
%%% [2.2] start iteration
v = v0; % initialization
k = 1; % current number of iteration
opts.tol = 1;  % Reduce the iterations of AP solver to save time
opts.maxTime = 10;
for isig = 1:length(maxiter) % extension for a series of noise levels
    nsigma = sigma(isig); 
    opt.sigma = nsigma;
    for iter = 1:maxiter(isig)
        %  AP solver    
        yb = abs(A(v(:)));
        y_new = (y-yb);
        [v_new] = solveGerchbergSaxton(A, At,y_new(:), v(:), opts);
        v_new = abs(v_new);
        v_new = v_new./max(max(v_new));
        v = v(:)+lambda*v_new; 
        % EN solver
        v = reshape(v,dims);
        v = ffdnet_denoise(v,[],opt);                 
        v = v./max(max(v));
        
        k = k+1;
    end
end
rec_LPR = v;
toc;
fprintf('complete----running time:%3.2fs \n',toc);
end