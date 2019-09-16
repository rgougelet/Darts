for k = [0 1 2 3 4 5 6 7 8 9 10]
    disp(['checking node ' int2str(k) '...']);
    system(['ssh computing-0-' int2str(k) ' -n killall amica15c >&/dev/null']);
end

disp(['checking computing...']);
system(['ssh computing -n killall amica15c >&/dev/null']);
