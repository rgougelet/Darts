function runmc(s,n_sd, n_obs, n_vars, n_boot, svnm)

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

	% generate regression terms
	b = rand(n_vars,1);
	y = X*b;
	res = n_sd*randn(n_obs,1);
	yr = y+res;
	yr = yr-mean(yr);

	% solve regression
	lm = fitlm(X,yr);
	X = [ones(size(yr)), X]; % assumes intercept
	n_preds = n_vars+1;
	bh = pinv(X)*yr;
	yh = X*bh;
	resh = yh-yr;
	mseh = sum(resh.^2)/(n_obs-n_vars-1);
	[q, ~] = qr([X,yh]);
	rmse = sqrt(mseh);
	r_sqr = corr(yr,yh)^2;
	adj_r_sqr = 1 - ((1-r_sqr)*(n_obs-1)) / (n_obs-n_vars-1);

	ssm = sum((yh-mean(yr)).^2);
	sse = sum(resh.^2);

	% get t stats
	bhlm = lm.Coefficients;
	dfm = n_vars-1; % assumes intercept
	dfe = n_obs-n_vars-1; % assumes intercept
	bhsig = sse/dfe;
	bhvc = bhsig*inv(X'*X);
	bhse = sqrt(diag(bhvc));
	bht = (bh-bh(1))./bhse;
	bhp = tcdf(abs(bht),dfe,'upper')*2; % hacky, probably why there's tiny error

	% get f stats
	fhlm = anova(lm, 'summary');
	dfm = n_preds-1; % assumes intercept
	dfe = n_obs-n_vars-1; % assumes intercept
	msm = ssm/dfm;
	mse = sse/dfe;
	F = msm/mse;
	F_p = fcdf(F,dfm,dfe, 'upper');

	% run bootstrap
	bhs = runboot(X,yr,resh,n_boot);

	% save
	clear q lm
	save(svnm,'-v6')
end

