%% Begin Main program
   %should remain dominated by user input. Prefer all work done in function
   %workspace. This workspace should hold handles determined by params of
   %data file. These are then passed to function workspaces.
   %Run this script to execute the program
%% Begin Setup
% Set output config
configparams();





%% User Inputs
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
pathfold = {pathE, pathB, pathD};

%% Analysis

get_dtimeplot(pathE)
get_B_and_DC_Params_6ch(pathD,pathB,tc,Wh,simulated,folderE,fileE,bint)
[Mcounts] = get_centroidPOS(fold, pathfold, 0);


