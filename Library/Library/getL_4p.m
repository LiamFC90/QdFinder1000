function [lnL, Pi, Pj] = getL_4p(Mcounts,N,x1guess,y1guess,x2guess,y2guess,...
    sigmai,sigmaiy,pfia,...
    sigmaj,sigmajy,pfja,...
    pfBa,pfBb,...
    PB1,PB2,PB3,PB4,...
    pfD1a,pfD1b,pfD2a,pfD2b,...
    pfD3a,pfD3b,pfD4a,pfD4b,...
    key1a,key2a,key3a,key4a,...
    key1b,key2b,key3b,key4b,...
    keyMulti,keyNone,A1,A2,...
    cf1,cf2,cf3,cf4)
%GETL_4P Initial spatial localization is made after a guess is made.
%   This will provide initial param data to main loop to allow a more
%   focused search. 

global x1
global x2
global x3
global x4
global y1
global y2
global y3
global y4
global w

L1 =0;
m1a = Mcounts(1);
m1b = Mcounts(5);
m2a = Mcounts(2);
m2b = Mcounts(6);
m3a = Mcounts(3);
m3b = Mcounts(7);
m4a = Mcounts(4);
m4b = Mcounts(8);
m2plus = Mcounts(9);

sigma1 = sigmai;
sigma1y = sigmaiy;
sigma2 = sigmaj;
sigma2y = sigmajy;

a1i = 0.5*(erf(((w/2)+x1-x1guess)/sqrt(2*(sigma1.^2)))-erf(((-w/2)+x1-x1guess)/sqrt(2*(sigma1.^2))))*...
    0.5*(erf(((w/2)+y1-y1guess)/sqrt(2*(sigma1y.^2)))-erf(((-w/2)+y1-y1guess)/sqrt(2*(sigma1y.^2))));
a2i = 0.5*(erf(((w/2)+x2-x1guess)/sqrt(2*(sigma1.^2)))-erf(((-w/2)+x2-x1guess)/sqrt(2*(sigma1.^2))))*...
    0.5*(erf(((w/2)+y2-y1guess)/sqrt(2*(sigma1y.^2)))-erf(((-w/2)+y2-y1guess)/sqrt(2*(sigma1y.^2))));
a3i = 0.5*(erf(((w/2)+x3-x1guess)/sqrt(2*(sigma1.^2)))-erf(((-w/2)+x3-x1guess)/sqrt(2*(sigma1.^2))))*...
    0.5*(erf(((w/2)+y3-y1guess)/sqrt(2*(sigma1y.^2)))-erf(((-w/2)+y3-y1guess)/sqrt(2*(sigma1y.^2))));
a4i = 0.5*(erf(((w/2)+x4-x1guess)/sqrt(2*(sigma1.^2)))-erf(((-w/2)+x4-x1guess)/sqrt(2*(sigma1.^2))))*...
    0.5*(erf(((w/2)+y4-y1guess)/sqrt(2*(sigma1y.^2)))-erf(((-w/2)+y4-y1guess)/sqrt(2*(sigma1y.^2))));

a1i = a1i*cf1;
a2i = a2i*cf2;
a3i = a3i*cf3;
a4i = a4i*cf4;

a1j = 0.5*(erf(((w/2)+x1-x2guess)/sqrt(2*(sigma2.^2)))-erf(((-w/2)+x1-x2guess)/sqrt(2*(sigma2.^2))))*...
    0.5*(erf(((w/2)+y1-y2guess)/sqrt(2*(sigma2y.^2)))-erf(((-w/2)+y1-y2guess)/sqrt(2*(sigma2y.^2))));
a2j = 0.5*(erf(((w/2)+x2-x2guess)/sqrt(2*(sigma2.^2)))-erf(((-w/2)+x2-x2guess)/sqrt(2*(sigma2.^2))))*...
    0.5*(erf(((w/2)+y2-y2guess)/sqrt(2*(sigma2y.^2)))-erf(((-w/2)+y2-y2guess)/sqrt(2*(sigma2y.^2))));
a3j = 0.5*(erf(((w/2)+x3-x2guess)/sqrt(2*(sigma2.^2)))-erf(((-w/2)+x3-x2guess)/sqrt(2*(sigma2.^2))))*...
    0.5*(erf(((w/2)+y3-y2guess)/sqrt(2*(sigma2y.^2)))-erf(((-w/2)+y3-y2guess)/sqrt(2*(sigma2y.^2))));
a4j = 0.5*(erf(((w/2)+x4-x2guess)/sqrt(2*(sigma2.^2)))-erf(((-w/2)+x4-x2guess)/sqrt(2*(sigma2.^2))))*...
    0.5*(erf(((w/2)+y4-y2guess)/sqrt(2*(sigma2y.^2)))-erf(((-w/2)+y4-y2guess)/sqrt(2*(sigma2y.^2))));

a1j = a1j*cf1;
a2j = a2j*cf2;
a3j = a3j*cf3;
a4j = a4j*cf4;

pfib = 1-pfia;
pfjb = 1-pfja;
pfi = pfia;
pfj = pfja;

% Pj = (A1+A2-Pi*(a1i+a2i+a3i+a4i))/(a1j+a2j+a3j+a4j);
% pfja = (A1-pfia*Pi*(a1i+a2i+a3i+a4i))/(Pj*(a1j+a2j+a3j+a4j));

Pi = (A1-pfja*(A1+A2))/((a1i+a2i+a3i+a4i)*(pfia-pfja));
Pj = (A1+A2-Pi*(a1i+a2i+a3i+a4i))/(a1j+a2j+a3j+a4j);

test = 1;
if Pj>=1
    lnL = NaN;
    test = 0;
end
if Pj<=0
    lnL = NaN;
    test = 0;
end
if pfj>=1
    lnL = NaN;
    test = 0;
end
if pfj<=0
    lnL = NaN;
    test = 0;
end

if test==1

    %combine emitter event and pixel arrival
    P0 = [Pi;1-Pi];
    Aimat = [a1i a2i a3i a4i 1-a1i-a2i-a3i-a4i];
    P1 = P0*Aimat;
    P1 = reshape(P1',[],1);

    %combine above with emitter arrival time
    pfi_mat = [pfia pfib 1-pfia-pfib];
    P2 = P1*pfi_mat;
    P2 = reshape(P2',[],1);

    %combine with second emitter event
    Pjtotal = [Pj 1-Pj];
    P2a = P2*Pjtotal;
    P2a = reshape(P2a',[],1);

    %combine with second emitter pixel arrival
    Ajmat = [a1j a2j a3j a4j 1-a1j-a2j-a3j-a4j];
    P2b = P2a*Ajmat;
    P2b = reshape(P2b',[],1);

    %combine with second emitter arrival time
    pfj_mat = [pfja pfjb 1-pfja-pfjb];
    P2c = P2b*pfj_mat;
    P2c = reshape(P2c',[],1);

    %combine above with background pixel 1 arrival
    PB1_mat = [PB1 1-PB1];
    P3 = P2c*PB1_mat;
    P3 = reshape(P3',[],1);

    %combine above with background pixel 2 arrival
    PB2_mat = [PB2 1-PB2];
    P4 = P3*PB2_mat;
    P4 = reshape(P4',[],1);

    %combine above with background pixel 3 arrival
    PB3_mat = [PB3 1-PB3];
    P5 = P4*PB3_mat;
    P5 = reshape(P5',[],1);

    %combine above with background pixel 4 arrival
    PB4_mat = [PB4 1-PB4];
    P6 = P5*PB4_mat;
    P6 = reshape(P6',[],1);

    %combine above with background arrival time
    PfB_mat = [pfBa pfBb 1-pfBa-pfBb];
    P7 = P6*PfB_mat;
    P7 = reshape(P7',[],1);

    %detector 1 dark counts
    pfD1 = [pfD1a pfD1b 1-pfD1a-pfD1b];
    P8 = P7*pfD1;
    P8 = reshape(P8',[],1);

    %detector 2 dark counts

    pfD2 = [pfD2a pfD2b 1-pfD2a-pfD2b];
    P9 = P8*pfD2;
    P9 = reshape(P9',[],1);

    %detector 3 dark counts

    pfD3 = [pfD3a pfD3b 1-pfD3a-pfD3b];
    P10 = P9*pfD3;
    P10 = reshape(P10',[],1);

    %detector 4 dark counts

    pfD4 = [pfD4a pfD4b 1-pfD4a-pfD4b];
    P11 = P10*pfD4;
    P11 = reshape(P11',[],1);

    PF = P11';

    % sort for observables
    %one photons detected
    Pm1a = PF*key1a;
    Pm1b = PF*key1b;
    Pm2a = PF*key2a;
    Pm2b = PF*key2b;
    Pm3a = PF*key3a;
    Pm3b = PF*key3b;
    Pm4a = PF*key4a;
    Pm4b = PF*key4b;
    P2plus = PF*keyMulti;
    Pnone = PF*keyNone;

    %% checkpoint
    % Pnone=Pnone
    % Ptest = 1-Pm1a-Pm2a-Pm3a-Pm4a-Pm1b-Pm2b-Pm3b-Pm4b-P2plus

    lnL = m1a*log(Pm1a)+m2a*log(Pm2a)+m3a*log(Pm3a)+m4a*log(Pm4a)+...
        m1b*log(Pm1b)+m2b*log(Pm2b)+m3b*log(Pm3b)+m4b*log(Pm4b)+...
        (N-m1a-m2a-m3a-m4a-m1b-m2b-m3b-m4b)*log(1-Pm1a-Pm2a-Pm3a-Pm4a-Pm1b-Pm2b-Pm3b-Pm4b);

end
end

