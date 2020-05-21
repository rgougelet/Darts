%% uipickfiles skeleton script
filePaths = uipickfiles('FilterSpec','STARTING_DIRECTORY', 'Type',{'*.FILE_EXTENSION', 'FILE_EXTENSION Files description'});

if isempty(filePaths)
    error('No files selected...')
end

fileIndex = 1;

fileString = filePaths{fileIndex};
[fileDir,fileName,fileExt] = fileparts(fileString);
