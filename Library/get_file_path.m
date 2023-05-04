function [path] = get_file_path(folder,file)
%UNTITLED grabs a file path to a specified data file
%   Detailed explanation goes here

path = which(strcat(string(folder),'_',string(file),'.ptu'));

if isempty(path)
    path = which(strcat(strcat(string(folder),'tr'),'_',string(file),'.ptu')); % checks if file is labeled with trace
end

if isempty(path)
    path = 'BAD PATH';
    fprintf('\nBad Path\n');
end