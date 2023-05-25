function [] = doSabine()
%DOSABINE Creates a sim file for the requested params
%   Detailed explanation goes here

printBreak  
T = input('Enter laser pulse time in ns(Default: 400):\n>>> ');
printLine
t = input("Enter sim run time in s(Default: 5):\n>>> ");
printLine
sigmai = input('Enter PSF size in nm(Default: 180):\n>>> ');
sigmaiy = sigmai;
printLine
fprintf('Enter the X and Y position of the simulated emitter in nm:\n')
mui = input('X = >>> ');
muiy = input('Y = >>> ');
printLine
Ri = input('Enter emitter lifetime in ns(Default: 50):\n>>> ');
Ri = 1/Ri;
RiCarry = Ri;
RB = input('Enter background lifetime in ns(Default: 3):\n>>> ');
RB = 1/RB;
printLine
fprintf('The next inputs are scaled for ease of use and entry\n')
Pi = input('Enter rate of emission(Default: 5 )\n>>> ');
Pi = Pi/100;
PiCarry = Pi;
printLine
RD1 = input('Enter rate of dark count emission(Default: 5)\n>>> ');
RD1 = RD1/10000000;
RD2 = RD1;
RD3 = RD1;
RD4 = RD1;
printLine
PB1 = input('Enter rate of background emission(Default: 5)\n>>> ');
PB1 = PB1/1000;
PB2 = PB1;
PB3 = PB1;
PB4 = PB1;
printBreak
fprintf('Starting Sim\nThis will take up to a couple hours\n')
fprintf('Creating file 1/3\n')
fprintf('Efile...\n')
%make E file
doSim1E(T,t,sigmai,sigmaiy, mui, muiy, Ri, RB, Pi, RD1, RD2, RD3, RD4, PB1, PB2, PB3, PB4, PiCarry, RiCarry)
fprintf('Made Efile. Saved to ...DATA... under 55555\n')
%Now, set Pi to 0 and run again to create background file
fprintf('Creating file 2/3\n')
fprintf('Bfile...\n')
Pi = 0;
doSim1E(T,t,sigmai,sigmaiy, mui, muiy, Ri, RB, Pi, RD1, RD2, RD3, RD4, PB1, PB2, PB3, PB4, PiCarry, RiCarry)
fprintf('Made Bfile. Saved to ...DATA... under 55555\n')
%now set PB1 to 0 and do it all again
fprintf('Creating file 3/3\n')
fprintf('DCfile...\n')
PB1 = 0; PB2 = 0; PB3 = 0; PB4 = 0;
doSim1E(T,t,sigmai,sigmaiy, mui, muiy, Ri, RB, Pi, RD1, RD2, RD3, RD4, PB1, PB2, PB3, PB4, PiCarry, RiCarry)
fprintf('Made DCfile. Saved to ...DATA... under 55555\n')
fprintf('Done with simulations\n')
end






