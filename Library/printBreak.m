function [] = printBreak(varargin)
%PRINTBREAK Prints a line to CL. if integer input is give, print that many
%new lines
narginchk(0,1)
varargin = cell2mat(varargin);
if isempty(varargin)
    varargin = 1;
end
for pp = 1:varargin
    for n = 1:4
        fprintf('________________________')
    end
    fprintf('\n')
end
end

