function pipeline = create_pipeline_structure(pipe_dir,srate,overwrite)
if overwrite
	pipeline = struct;
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
	pipeline.Name = pipe_dir(1:end-1);
	pipeline.PipeDir = pipe_dir;
	pipeline.SampleRate = srate;
	pipeline.SubjectIDs = subjs_to_include;
	
	for subj_i = 1:length(subjs_to_include)
		pipeline.Subjects(subj_i).ID = subjs_to_include{subj_i};
		pipeline.Subjects(subj_i).Set = ...
			dir([pipe_dir, pipeline.Subjects(subj_i).ID,'*_',num2str(pipeline.SampleRate),'.set']);
	end
	save([pipe_dir,'pipeline_structure.mat'],'pipeline', '-v7.3');
else
	load([pipe_dir,'pipeline_structure.mat'],'pipeline');
end