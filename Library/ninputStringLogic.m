function [out] = ninputStringLogic(varargin)
%NINPUTSTRINGLOGIC takes n input string chars and provides outcome logic
%Entered in order (var1, var2, ... , targetvar)
%LiamK 2023







%% validate input
minArgs = 3; %[char1,char2,userchar]
maxArgs = nargin;
narginchk(minArgs,maxArgs)
clear minArgs maxArgs
%% Load and prepare alphastring
txtalpha = "ABCDEFGHIJKLMONPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
pat = characterListPattern(txtalpha);
alphastring = string(extract(txtalpha,pat));
alphastring = [alphastring(1:26), alphastring(27:end)];
clear txtalpha pat
%% start
varargin = string(varargin)';
targetvar = varargin(end);
varargin = varargin(1:end-1);
realsize = numel(varargin);
%% add alpha cap/lowercase as needed
vararginadd = [];
for alphaselect = [1 2]
    for qq = 1:length(alphastring)
        for pp = 1:length(varargin)
            %fprintf('checking if %s matches %s\n',alphastring(qq,alphaselect), varargin(pp))
            if strcmp(alphastring(qq,alphaselect), varargin(pp))
                %fprintf('Match found\n')
                if alphaselect == 1
                    vararginadd = [vararginadd; alphastring(qq,2)];
                elseif alphaselect == 2
                    vararginadd = [vararginadd; alphastring(qq,1)];
                end
            end
        end
    end
end
varargin = [varargin; vararginadd];
n = 1;
while n <= numel(varargin)    
        if strcmp(targetvar,varargin(n))
            out = n;
            n = numel(varargin)+1;
        else
            n = n + 1;
        end
end
%Adjust output logic if case was found to be switched

if exist("out","var") == 0
    warning('Unknown input')
    out = NaN;
else
    if out > realsize
        out = out - realsize;
    end
end

