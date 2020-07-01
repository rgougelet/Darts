function subjects = old_create_subject_structure(pipe_dir,srate,overwrite)
if overwrite
	subjects = struct;
	subjs_to_include = {
		'571'
		'579'
		'580'
		'607'
		'608'
		'616'
		'619'
		'621'
		'627'
		'631'
		};
	subjects.PipeDir = pipe_dir;
	subjects.SampleRate = srate;

	for subj_i = 1:length(subjs_to_include)
		subjects(subj_i).ID = subjs_to_include{subj_i};
		subjects(subj_i).PipeDir = pipe_dir;
		subjects(subj_i).SampleRate = srate;
		subjects(subj_i).Set = ...
			dir([pipe_dir, subjects(subj_i).ID,'*_',num2str(subjects(subj_i).SampleRate),'.set']);
	end
	save([pipe_dir,'subject_structure.mat'],'subjects', '-v7.3');
else
	load([pipe_dir,'subject_structure.mat'],'subjects');
end