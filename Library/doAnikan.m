function [Mcounts,foldcent] = doAnikan()
%DOANIKAN Sets up and requests an analysis job
%   Detailed explanation goes here
fprintf('Analysis\n')
%%put these somewhere better %%
bint = input('Enter BIN time(s)? \n>>>');
tc = input('Enter TIME CUT(ns)?\n>>>');
Wh = .96; %histogram bin width
simulated = 0;
sigmai= 180; sigmaiy = 180;%PSF1 width
sigmaj= 180; sigmajy = 180;%PSF2 width

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
foldpath = {pathE, pathB, pathD};

%% Analysis

get_dtimeplot(pathE)
get_B_and_DC_Params_6ch(pathD,pathB,tc,Wh,simulated,folderE,fileE,bint)
[Mcounts,foldcent] = get_centroidPOS(fold, foldpath, 0);


end

