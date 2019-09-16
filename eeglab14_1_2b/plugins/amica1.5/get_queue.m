function queue = get_queue(numcpus)
% Returns a free queue with specified number of CPUs.

% Clemens Brunner <clbrunner@ucsd.edu>, 02/23/2012

queue = [];

queues{32} = {'q3', 'q4', 'q5', 'q6', 'q7', 'q8', 'q9', 'q10'};
queues{64} = {'qa1', 'qa2', 'qa3', 'qa4'};
queues{128} = {'qb1', 'qb2'};
queues{256} = {'qc1'};

[status, output] = system('qstat -g c');

c = textscan(output, '%s%f%d%d%d%d%d%d', 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'Headerlines', 2);

for q = queues{numcpus}
    line = find(strcmp(q, c{1}));
    if c{5}(line) == numcpus
        queue = q{1};
        return
    end
end
     