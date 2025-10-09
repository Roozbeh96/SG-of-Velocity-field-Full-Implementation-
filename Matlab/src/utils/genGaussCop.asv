function [rand_u,rand_w] = genGaussCop(rho, n)
% for reproducibility, you can comment it if you want to generate different
% numbers every time.
rng(137); 
% Convert target Pearson correlation for Uniforms to underlying Gaussian
% correlation
r = 2 * sin((pi/6) * rho);
P = toeplitz([1 r]); % Correlation matrix
d = size(P, 1); % Dimension
% Generate sample
RANDS = mvnrnd(zeros(d,1), P, n);
%RANDS have Gaussian distribution. We need to uniform it.
RANDS = normcdf(RANDS);
%Now random numbers are uniformly distributed. You can double check with
%the following code:
% histogram(RANDS(:,1),'Normalization','pdf','NumBins',40)
% histogram(RANDS(:,2),'Normalization','pdf','NumBins',40)
rand_u = RANDS(:,1);
rand_w = RANDS(:,2);
end