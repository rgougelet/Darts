clc; clear; close all;
corr_step = .1;
v_corr = (corr_step:corr_step:1);
n_sims = 1000;
n_obs = 10000;
n_var = 10;
b = rand(1,n_var);
res = randn(1,n_obs);
Z = randn(n_var,n_obs);
corr_i = nchoosek(1:n_var,2);
% corrs = -1+2*rand(1,length(corr_i));
corrs = rand(1,length(corr_i));

corrm = nan(n_var);
corrm(tril(true(size(corrm)),-1)) = corrs;
corrm = corrm';
corrm(tril(true(size(corrm)),-1)) = corrs;
corrm(eye(size(corrm))==1) = 1;
[V,D] = eig(corrm);
X = V*sqrt(D)*Z;
X = X';
heatmap(corrm)
corrx = corr(X);
figure; heatmap(corrx)

% figure; heatmap(corr(X'))
% for corr_ii = 1:length(corr_i)
% 	i1 = corr_i(corr_ii,1);
% 	i2 = corr_i(corr_ii,2);
% 	covv(i1,i2) = cov(corr_ii);
% 	covv(i2,i1) = cov(corr_ii);
% end
% covv(eye(size(covv))==1) = var(z);
% [L,p] = chol(covv,'lower');
% tril(true(size(covv)),-1) = triu(true(size(covv)),1)
% clear szs
% close all
% for n_var = 2:100
% 	corr_i = nchoosek(1:n_var,2);
% 	szs(n_var) = length(corr_i);
% end
% figure; plot(1./szs)
% mod(1,1./szs)
% for avg_corr = v_corr
% 	corrs = (1:length(corr_i))./length(corr_i);
% 
% 	corr_mat = avg_corr*rand(1,length(corr_i))

% end

% y = b*x;