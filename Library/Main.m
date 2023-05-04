%% Begin Main program
   %should remain dominated by user input. Prefer all work done in function
   %workspace. This workspace should hold handles determined by params of
   %data file. These are then passed to function workspaces.







%% User Inputs
% E-File
fprintf('Enter data file information:')
[folderE, fileE] = get_fofi();
pathE = get_file_path(folderE, fileE);
savepath = strcat(get_path(), '\Analysis\',string(folderE),'\',string(fileE));

% B-File
fprintf('Enter background file information:')
[folderB, fileB] = get_fofi();
pathB = get_file_path(folderB, fileB);

% DC-File
fprintf('Enter data file information:')
[folderD, fileD] = get_fofi();
pathD = get_file_path(folderD, fileD);

% CF-Factors
CF_List = get_cf();



%% Analysis

