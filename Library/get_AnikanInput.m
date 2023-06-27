function [fold,foldpath] = get_AnikanInput(jobtype)
%get_AnikanInput gathers file info for data files
%gathers input based on job type


%% RAW
if strcmp(jobtype, 'rawptu')
    % E-File
    fprintf('Enter data file information:')
    [folderE, fileE] = get_fofi();
    pathE = get_file_path(folderE, fileE);
    printLine(20)
    printBreak(2)

    % B-File
    fprintf('Enter background file information:')
    [folderB, fileB] = get_fofi();
    pathB = get_file_path(folderB, fileB);
    printLine(20)
    printBreak(2)

    % DC-File
    fprintf('Enter dark count file information:')
    [folderD, fileD] = get_fofi();
    pathD = get_file_path(folderD, fileD);
    printLine(20)
    printBreak(2)

    % CF-Factors
    CF_List = get_cf();
    printLine(20)
    printBreak(2)

    %Build data vars
    fold = [folderE, fileE, folderB, fileB, folderD, fileD, CF_List];
    foldpath = {pathE, pathB, pathD};

elseif strcmp(jobtype,'sim')
    %File information
    XposNM = input('Enter X position in nm\n>>>  ');
    YposNM = input('Enter Y position in nm\n>>>  ');
    TauVal = input('Enter Tau in ns\n>>>  ');
    ratevl = input('Enter emission rate\n>>>  ');
    
    foldpath = strcat(cd,'\Data\55555\X_',string(XposNM),'_Y_',string(YposNM),'\Tau_',string(TauVal),'\rate_',string(ratevl));
    foldpath = [strcat(foldpath,'\E.mat'), strcat(foldpath,'\B.mat'), strcat(foldpath,'\DC.mat')];
    fold = [XposNM, YposNM, TauVal, ratevl];
    
    %C:\Users\liamk\MATLAB\Projects\QDwork\Data\55555\X_5_Y_10\Tau_5\rate_0.005


end

