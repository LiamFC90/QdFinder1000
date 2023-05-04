function [path] = get_file_path(folder,file)
%UNTITLED grabs a file path to a specified data file
%   Detailed explanation goes here

badpath = 1; %exit condition
while badpath == 1

    path = which(strcat(string(folder),'_',string(file),'.ptu'));

    if isempty(path)
        path = which(strcat(strcat(string(folder),'tr'),'_',string(file),'.ptu')); % checks if file is labeled with trace(tr)
    end

    if isempty(path) %tell user that they entered a bad path and ask for new path input
        fprintf('\nBad Path. Can"t find file. Rekey inputs. Input 55555 folder to exit\n');
        [folder, file] = get_fofi();
    end

    if ~isempty(path) %set exit cond when path is found
        badpath = 0;
    end

    if folder == 55555
        error('User requested exit:get_file_path:BAD PATH:55555')
    end
end

