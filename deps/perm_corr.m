function [r,p] = perm_corr(x,y)
	r = corr(x,y);
	for perm_i = 1:10000
		rs(perm_i) = corr(x,y(randperm(length(y))));
	end
	p = nansum(abs(rs)>=abs(r))/10000;
end

