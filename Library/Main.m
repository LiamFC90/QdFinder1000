clc








%% User Inputs
% E-File
fprintf('Enter data file information:\n')
folderE = input('What is the data folder?\n');
fileE = input('What is the data file?\n');
pathE = get_file_path(folderE, fileE);
savepath = strcat(get_path(), '\Analysis\',string(folderE),'\',string(fileE));
% B-File
fprintf('Enter background file information:\n')
folderB = input('What is the background folder?\n');
fileB = input('What is the background file?\n');
pathB = get_file_path(folderB, fileB);
% DC-File
fprintf('Enter data file information:\n')
folderD = input('What is the dark count folder?\n');
fileD = input('What is the dark count file?\n');
pathD = get_file_path(folderD, fileD);

% CF-Factors
fprintf('\nEnter CF factors\n')
CF1 = input('CF1==> ');
CF2 = input('\nCF2==>  ');
CF3 = input('\nCF3==>  ');
CF4 = input('\nCF4==>  ');
CFlist = [CF1 CF2 CF3 CF4];
clear CF1 CF2 CF3 CF4

