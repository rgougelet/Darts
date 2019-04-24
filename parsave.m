function parsave(fname,var,varnames)
	if ischar(varnames)
		eval([varnames,' = var;']) 
		save(fname, varnames,'-v7.3')
	end
	if iscell(varnames)
		for var_ind = 1:length(varnames)
			eval([varnames{var_ind},' = var{',num2str(var_ind),'};'])
			if var_ind == 1
				save(fname, varnames{var_ind},'-v7.3')
			else
				save(fname, varnames{var_ind},'-append')
			end
		end
	end
end