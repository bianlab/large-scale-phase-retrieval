function [denoisedv] = ffdnet_denoise(noisyv,orgv,para)

format compact;
global sigmas; % input noise level or input noise level map
useGPU = false;
ffdnetvnorm = true;
if isfield(para,'useGPU'),      useGPU = para.useGPU;             end
if isfield(para,'ffdnetvnorm'), ffdnetvnorm = para.ffdnetvnorm;   end

if ffdnetvnorm
    maxz = max(noisyv(:));
    minz = min(noisyv(:));
    scale = 0.7;
    shift = (1-scale)/2;
    noisyv = (noisyv-minz)/(maxz-minz);
    noisyv = noisyv*scale+shift;

    sigmas = para.sigma/(maxz-minz)*scale;
else
    % set noise level map
    sigmas = para.sigma; % see "vl_simplenn.m".
end


if isfield(para,'net') && ~isempty(para.net)
    net = para.net;
else
    load(fullfile('models','FFDNet_gray.mat'),'net');
    net = vl_simplenn_tidy(net);
    if useGPU
        net = vl_simplenn_move(net, 'gpu') ;
    end
end



    input = single(noisyv);
    
    % perform denoising  
    if useGPU
        input = gpuArray(input);
        res    = vl_simplenn_ffdnet(net,input,[],[],'conserveMemory',true,'mode','test'); % matconvnet default
        output = res(end).x;
        output = gather(output);
    else
        res    = vl_simplenn_ffdnet(net,single(input),[],[],'conserveMemory',true,'mode','test'); % matconvnet default
        output = res(end).x;
    end
    
    denoisedv = double(output);


if ffdnetvnorm

    denoisedv = (denoisedv-shift)/scale;
    denoisedv = denoisedv*(maxz-minz)+minz;

end

end

