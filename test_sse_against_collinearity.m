clear; clc;
n_vars = 100;
n_obs = 10000;
n_sd = 5;
s = 10;

% generate covariance/correlation matrix
eigs = diag([s (1/s)*ones(1,n_vars-1)]);
[p,~] = qr(randn(n_vars));
C = p'*eigs*p;
D = inv(diag(sqrt(diag(C))));
CR = D*C*D;

% generate random data with given covariance/correlation
T = chol(C);
X = randn(n_obs,n_vars)*T;
X = X-mean(X);
XR = corr(X);
mean(mean(abs(XR)))

% generate regression terms
b = rand(n_vars,1);
y = X*b;
res = n_sd*randn(n_obs,1)+10;
yr = y+res;
% yr = yr-mean(yr);

% solve regression
lm = fitlm(X,yr);
X = [ones(size(yr)), X]; % assumes intercept
n_preds = n_vars+1;
bh = pinv(X)*yr;
yh = X*bh;
resh = yh-yr;
sum(resh.^2)

c = ones(size(yr))*mean(yr);
c_res = yr-c;
sum(c_res.^2)
