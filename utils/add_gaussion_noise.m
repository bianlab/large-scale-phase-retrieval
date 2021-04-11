function [Y] = add_gaussion_noise(X,SNR)
% 给图片添加指定信噪比的高斯白噪声
% X为输入图像，SNR为信噪比
X = double(X);
signal_power = sum(X(:).^2);
noise_power = signal_power/(10^(SNR/10));
picture_size = size(X);
noise_var = noise_power/prod(picture_size(:));
Y = X+sqrt(noise_var)*randn(picture_size);
Y = abs(Y);
end

