function [fold,foldpath] = get_AnikanInput()
%get_AnikanInput gathers file info for data files

%% RAW
% E-File
fprintf('Enter data file information:')
[folderE, fileE] = get_fofi();
pathE = get_file_path(folderE, fileE);

% B-File
fprintf('Enter background file information:')
[folderB, fileB] = get_fofi();
pathB = get_file_path(folderB, fileB);

% DC-File
fprintf('Enter dark count file information:')
[folderD, fileD] = get_fofi();
pathD = get_file_path(folderD, fileD);

% CF-Factors
CF_List = get_cf();

%Build data vars
fold = [folderE, fileE, folderB, fileB, folderD, fileD, CF_List, tc, bint,sigmai, sigmaiy, sigmaj, sigmajy];
foldpath = {pathE, pathB, pathD};
end

