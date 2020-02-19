function boot = runboot(X,yr,resh,n_boot)
	boot.bhs = nan(n_boot,size(X,2));
	for i = 1:n_boot
		boot.sh = randi(length(resh),length(resh),1);
		boot.y = yr+resh(boot.sh);
		boot.bhs(i,:) = X\boot.y;
	end
end

