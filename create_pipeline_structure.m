function pipeline = create_pipeline_structure(data_dir,pipe_name,srate)
pipe_dir = [data_dir,pipe_name,'/'];
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
pipeline.Name = pipe_name;
pipeline.DataDir = data_dir;
pipeline.PipeDir = pipe_dir;
pipeline.SampleRate = srate;
pipeline.SubjectIDs = subjs_to_include;
pipeline.nSubjects = length(subjs_to_include);
mkdir(pipe_dir)
for subj_i = 1:length(subjs_to_include)
	pipeline.Subjects(subj_i).ID = subjs_to_include{subj_i};
	pipeline.Subjects(subj_i).Index = subj_i;
	pipeline.Subjects(subj_i).Set = ...
		dir([pipe_dir, pipeline.Subjects(subj_i).ID,'*_',num2str(pipeline.SampleRate),'.set']);
	pipeline.Subjects(subj_i).SampleRate = pipeline.SampleRate;
end
