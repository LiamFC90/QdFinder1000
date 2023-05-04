function [path] = get_file_path(folder,file)
%get_file_path grabs a file path to a specified data file
%

badpath = 1; %exit condition
while badpath == 1

    path = which(strcat(string(folder),'_',string(file),'.ptu'));

    if isempty(path)
        path = which(strcat(strcat(string(folder),'tr'),'_',string(file),'.ptu')); % checks if file is labeled with trace(tr)
    end

    if isempty(path) %tell user that they entered a bad path and ask for new path input
        fprintf('\nBad Path. Can"t find file. Rekey inputs. leave any input empty to exit\n');
        [folder, file] = get_fofi();
    end

    if ~isempty(path) %set exit cond when path is found
        badpath = 0;
    end

    if isempty(folder)
        error('User requested exit:get_file_path:folderinput:BAD PATH')%outdated. moved error handle to get_fofi
    end
    if isempty(file)
        error('User requested exit:get_file_path:fileinput:BAD PATH')%oudated
    end

end

