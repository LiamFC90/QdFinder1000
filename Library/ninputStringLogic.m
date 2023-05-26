function [out] = ninputStringLogic(varargin)
%NINPUTSTRINGLOGIC takes n input string chars and provides outcome logic
%Entered in order

%% validate input
minArgs = 3; %[char1,char2,userchar]
maxArgs = nargin;
narginchk(minArgs,maxArgs)
clear minArgs maxArgs
%% Load and prepare alphastring
txtalpha = "ABCDEFGHIJKLMONPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
pat = characterListPattern(txtalpha);
alphastring = string(extract(txtalpha,pat));
clear txtalpha pat
%% start
varargin = string(varargin)';
targetvar = varargin(end);
varargin = varargin(1:end-1);
%% add alpha cap/lowercase as needed


n = 1;
while n <= numel(varargin)    
        if strcmp(targetvar,varargin(n))
            out = n;
            n = numel(varargin)+1;
        else
            n = n + 1;
        end
end


end

