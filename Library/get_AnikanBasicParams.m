function [foldbasic] = get_AnikanBasicParams()
%GET_ANIKANBASICPARAMS gathers and returns basic params



fprintf('Enter basic params now:\n')
printLine

bint = input('Enter BIN time(s) \n>>> ');
tc = input('Enter TIME CUT(ns) \n>>> ');
Wh = .96; %histogram bin width
sigmai = input('Enter PSF X-width in nm [Default:180] \n>>> ');
sigmaj = sigmai;
%sigmai= 180; sigmaiy = 180;%PSF1 width
sigmaiy = input('Enter PSF Y-width in nm [Default:180] \n>>> ');
sigmajy = sigmaiy;
%sigmaj= 180; sigmajy = 180;%PSF2 width

foldbasic = [bint tc Wh sigmai sigmaiy sigmaj sigmajy];

end