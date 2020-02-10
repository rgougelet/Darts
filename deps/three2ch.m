function chanlocs = three2ch(chanlocs, three)
	x = num2cell(three(:,1));
	[chanlocs.X] = x{:};
	y = num2cell(three(:,2));
	[chanlocs.Y] = y{:};
	z = num2cell(three(:,3));
	[chanlocs.Z] = z{:};
end

