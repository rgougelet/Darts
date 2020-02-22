function bhs = runboot(X,y,res,n_boot)
	bhs = nan(n_boot,size(X,2));
	nanIdx = any(isnan(X),2) | isnan(y) | isnan(res); 
	res = res(~nanIdx); X = X(~nanIdx,:); y = y(~nanIdx);
	X = normalize(X);
	for i = 1:n_boot
		if mod(i,1000) == 0
			disp([num2str(100*(i/n_boot)),'% done']);
		end
		sh = randi(length(res),length(res),1);
		y = y+res(sh);
		bhs(i,:) = X\y ;
	end
end

