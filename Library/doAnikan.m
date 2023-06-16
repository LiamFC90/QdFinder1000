function [Mcounts,foldcent] = doAnikan(jobtype)
%DOANIKAN Sets up and requests an analysis job
%   Detailed explanation goes here
printLine(20)
printBreak(2)
fprintf('Analysis: jobtype is %s\n',jobtype)
printBreak
%% Gather basic information


%%%foldbasic = [bint tc Wh sigmai sigmaiy sigmaj sigmajy];

%% start job
if strcmp(jobtype,'rawptu')
    [fold,foldpath] = get_AnikanInput(jobtype);
    get_dtimeplot(cell2mat(foldpath(1)))
    %show user decay before asking for tc selection
    [foldbasic] = get_AnikanBasicParams;
    get_B_and_DC_Params_6ch(cell2mat(foldpath(3)),cell2mat(foldpath(2)),foldbasic(2),foldbasic(3),0,fold(1),fold(2),foldbasic(1))
    [Mcounts,foldcent] = get_centroidPOS(fold, foldpath, foldbasic, 0);

elseif strcmp(jobtype,'sim')
    [fold,foldpath] = get_AnikanInput(jobtype);
    % $work$
end

%%% fold = [folderE, fileE, folderB, fileB, folderD, fileD, CF1, CF2, CF3, CF4, tc, bint,sigmai, sigmaiy, sigmaj, sigmajy];
%%% foldpath = {pathE, pathB, pathD};


% get_dtimeplot(pathE)
% get_B_and_DC_Params_6ch(pathD,pathB,tc,Wh,simulated,folderE,fileE,bint)
% [Mcounts,foldcent] = get_centroidPOS(fold, foldpath, 0);


end

