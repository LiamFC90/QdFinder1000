function [] = printLine(varargin)
%PRINTLINE prints a newline
narginchk(0,1)
varargin = cell2mat(varargin);

if isempty(varargin)
    varargin = 1;
end

for n = 1:varargin
    fprintf('\n')
end
end

