function bhs = runboot(X,yr,resh,n_boot)
	bhs = nan(n_boot,size(X,2));
	for i = 1:n_boot
		sh = randi(length(resh),length(resh),1);
		y = yr+resh(sh);
		bhs(i,:) = X\y;
	end
end

