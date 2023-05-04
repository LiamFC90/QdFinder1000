function [folder,file] = get_fofi()
%GET_FOFI Gets folder and file information from the user
%   The user is prompted to enter information to id the target file. Any
%   blank input will return an error and exit the program. 

folder = input('\nWhat is the folder number?\n>>> ');
file = input('What is the file number?\n>>> ');

if isempty(folder)||isempty(file)
    error('Input empty:get_fofi:BADINPUT')
end

