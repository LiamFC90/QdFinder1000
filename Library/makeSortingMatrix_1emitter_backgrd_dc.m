fprintf('Creating sorting key: SingleEmitter\n')



%combine emitter event and pixel arrival

% P0 = [Pi;1-Pi];
imat0 = [1;0];
% Aimat = [a1i a2i a3i a4i 1-a1i-a2i-a3i-a4i];
imat1 = [1 0 0 0 0];
imat2 = [0 1 0 0 0];
imat3 = [0 0 1 0 0];
imat4 = [0 0 0 1 0];

% P1 = P0*Aimat;
imat1 = imat0*imat1;
imat2 = imat0*imat2;
imat3 = imat0*imat3;
imat4 = imat0*imat4;

clear imat0

% P1 = reshape(P1',[],1);
imat1 = reshape(imat1',[],1);
imat2 = reshape(imat2',[],1);
imat3 = reshape(imat3',[],1);
imat4 = reshape(imat4',[],1);

%combine above with emitter arrival time

% pfi_mat = [pfia pfib 1-pfia-pfib];
imat1a = [1 0 0];
imat2a = imat1a;
imat3a = imat1a;
imat4a = imat1a;
imat1b = [0 1 0];
imat2b = imat1b;
imat3b = imat1b;
imat4b = imat1b;

% P2 = P1*pfi_mat;
imat1a = imat1*imat1a;
imat2a = imat2*imat2a;
imat3a = imat3*imat3a;
imat4a = imat4*imat4a;
imat1b = imat1*imat1b;
imat2b = imat2*imat2b;
imat3b = imat3*imat3b;
imat4b = imat4*imat4b;

clear imat1 imat2 imat3 imat4

% P2 = reshape(P2',[],1);
imat1a = reshape(imat1a',[],1);
imat2a = reshape(imat2a',[],1);
imat3a = reshape(imat3a',[],1);
imat4a = reshape(imat4a',[],1);
imat1b = reshape(imat1b',[],1);
imat2b = reshape(imat2b',[],1);
imat3b = reshape(imat3b',[],1);
imat4b = reshape(imat4b',[],1);

%background detector 1 contribution
B1 = [1 0];
Bmat = ones(length(imat1a),1);
filler = [1 1];

Bmat1 = Bmat*B1;
Bmat1 = reshape(Bmat1',[],1);

clear B1

imat1a = imat1a*filler;
imat2a = imat2a*filler;
imat3a = imat3a*filler;
imat4a = imat4a*filler;
imat1b = imat1b*filler;
imat2b = imat2b*filler;
imat3b = imat3b*filler;
imat4b = imat4b*filler;
imat1a = reshape(imat1a',[],1);
imat2a = reshape(imat2a',[],1);
imat3a = reshape(imat3a',[],1);
imat4a = reshape(imat4a',[],1);
imat1b = reshape(imat1b',[],1);
imat2b = reshape(imat2b',[],1);
imat3b = reshape(imat3b',[],1);
imat4b = reshape(imat4b',[],1);

%background detector 2 contribution
B2 = [1 0];
Bmat = ones(length(imat1a),1);
filler = [1 1];

Bmat2 = Bmat*B2;
Bmat2 = reshape(Bmat2',[],1);

clear B2

Bmat1 = Bmat1*filler;
Bmat1 = reshape(Bmat1',[],1);

imat1a = imat1a*filler;
imat2a = imat2a*filler;
imat3a = imat3a*filler;
imat4a = imat4a*filler;
imat1b = imat1b*filler;
imat2b = imat2b*filler;
imat3b = imat3b*filler;
imat4b = imat4b*filler;
imat1a = reshape(imat1a',[],1);
imat2a = reshape(imat2a',[],1);
imat3a = reshape(imat3a',[],1);
imat4a = reshape(imat4a',[],1);
imat1b = reshape(imat1b',[],1);
imat2b = reshape(imat2b',[],1);
imat3b = reshape(imat3b',[],1);
imat4b = reshape(imat4b',[],1);

%background detector 3 contribution
B3 = [1 0];
Bmat = ones(length(imat1a),1);
filler = [1 1];

Bmat3 = Bmat*B3;
Bmat3 = reshape(Bmat3',[],1);

clear B3

Bmat1 = Bmat1*filler;
Bmat1 = reshape(Bmat1',[],1);
Bmat2 = Bmat2*filler;
Bmat2 = reshape(Bmat2',[],1);

imat1a = imat1a*filler;
imat2a = imat2a*filler;
imat3a = imat3a*filler;
imat4a = imat4a*filler;
imat1b = imat1b*filler;
imat2b = imat2b*filler;
imat3b = imat3b*filler;
imat4b = imat4b*filler;
imat1a = reshape(imat1a',[],1);
imat2a = reshape(imat2a',[],1);
imat3a = reshape(imat3a',[],1);
imat4a = reshape(imat4a',[],1);
imat1b = reshape(imat1b',[],1);
imat2b = reshape(imat2b',[],1);
imat3b = reshape(imat3b',[],1);
imat4b = reshape(imat4b',[],1);

%background detector 4 contribution
B4 = [1 0];
Bmat = ones(length(imat1a),1);
filler = [1 1];

Bmat4 = Bmat*B4;
Bmat4 = reshape(Bmat4',[],1);

clear B4

Bmat1 = Bmat1*filler;
Bmat1 = reshape(Bmat1',[],1);
Bmat2 = Bmat2*filler;
Bmat2 = reshape(Bmat2',[],1);
Bmat3 = Bmat3*filler;
Bmat3 = reshape(Bmat3',[],1);

imat1a = imat1a*filler;
imat2a = imat2a*filler;
imat3a = imat3a*filler;
imat4a = imat4a*filler;
imat1b = imat1b*filler;
imat2b = imat2b*filler;
imat3b = imat3b*filler;
imat4b = imat4b*filler;
imat1a = reshape(imat1a',[],1);
imat2a = reshape(imat2a',[],1);
imat3a = reshape(imat3a',[],1);
imat4a = reshape(imat4a',[],1);
imat1b = reshape(imat1b',[],1);
imat2b = reshape(imat2b',[],1);
imat3b = reshape(imat3b',[],1);
imat4b = reshape(imat4b',[],1);

%background time arrival
pfBa_mat = [1 0 0];
pfBb_mat = [0 1 0];
filler = [1 1 1];

Bmat1a = Bmat1*pfBa_mat;
Bmat1b = Bmat1*pfBb_mat;
Bmat1a = reshape(Bmat1a',[],1);
Bmat1b = reshape(Bmat1b',[],1);
Bmat2a = Bmat2*pfBa_mat;
Bmat2b = Bmat2*pfBb_mat;
Bmat2a = reshape(Bmat2a',[],1);
Bmat2b = reshape(Bmat2b',[],1);
Bmat3a = Bmat3*pfBa_mat;
Bmat3b = Bmat3*pfBb_mat;
Bmat3a = reshape(Bmat3a',[],1);
Bmat3b = reshape(Bmat3b',[],1);
Bmat4a = Bmat4*pfBa_mat;
Bmat4b = Bmat4*pfBb_mat;
Bmat4a = reshape(Bmat4a',[],1);
Bmat4b = reshape(Bmat4b',[],1);

clear pfBa_mat pfBb_mat

imat1a = imat1a*filler;
imat2a = imat2a*filler;
imat3a = imat3a*filler;
imat4a = imat4a*filler;
imat1b = imat1b*filler;
imat2b = imat2b*filler;
imat3b = imat3b*filler;
imat4b = imat4b*filler;
imat1a = reshape(imat1a',[],1);
imat2a = reshape(imat2a',[],1);
imat3a = reshape(imat3a',[],1);
imat4a = reshape(imat4a',[],1);
imat1b = reshape(imat1b',[],1);
imat2b = reshape(imat2b',[],1);
imat3b = reshape(imat3b',[],1);
imat4b = reshape(imat4b',[],1);

%detector 1 dark counts

% pfD1 = [pfD1a pfD1b 1-pfD1a-pfD1b];
dc = ones(length(imat1a),1);

D1a = [1 0 0];
D1b = [0 1 0];
filler = [1 1 1];

d1a = dc*D1a;
d1b = dc*D1b;
d1a = reshape(d1a',[],1);
d1b = reshape(d1b',[],1);

clear D1a D1b

imat1a = imat1a*filler;
imat2a = imat2a*filler;
imat3a = imat3a*filler;
imat4a = imat4a*filler;
imat1b = imat1b*filler;
imat2b = imat2b*filler;
imat3b = imat3b*filler;
imat4b = imat4b*filler;
imat1a = reshape(imat1a',[],1);
imat2a = reshape(imat2a',[],1);
imat3a = reshape(imat3a',[],1);
imat4a = reshape(imat4a',[],1);
imat1b = reshape(imat1b',[],1);
imat2b = reshape(imat2b',[],1);
imat3b = reshape(imat3b',[],1);
imat4b = reshape(imat4b',[],1);
Bmat1a = Bmat1a*filler;
Bmat2a = Bmat2a*filler;
Bmat3a = Bmat3a*filler;
Bmat4a = Bmat4a*filler;
Bmat1b = Bmat1b*filler;
Bmat2b = Bmat2b*filler;
Bmat3b = Bmat3b*filler;
Bmat4b = Bmat4b*filler;
Bmat1a = reshape(Bmat1a',[],1);
Bmat2a = reshape(Bmat2a',[],1);
Bmat3a = reshape(Bmat3a',[],1);
Bmat4a = reshape(Bmat4a',[],1);
Bmat1b = reshape(Bmat1b',[],1);
Bmat2b = reshape(Bmat2b',[],1);
Bmat3b = reshape(Bmat3b',[],1);
Bmat4b = reshape(Bmat4b',[],1);

%detector 2 dark counts
% pfD2 = [pfD2a pfD2b 1-pfD2a-pfD2b];
dc = ones(length(imat1a),1);

D2a = [1 0 0];
D2b = [0 1 0];
filler = [1 1 1];

d2a = dc*D2a;
d2b = dc*D2b;
d2a = reshape(d2a',[],1);
d2b = reshape(d2b',[],1);

clear D2a D2b

d1a = d1a*filler;
d1b = d1b*filler;
d1a = reshape(d1a',[],1);
d1b = reshape(d1b',[],1);
imat1a = imat1a*filler;
imat2a = imat2a*filler;
imat3a = imat3a*filler;
imat4a = imat4a*filler;
imat1b = imat1b*filler;
imat2b = imat2b*filler;
imat3b = imat3b*filler;
imat4b = imat4b*filler;
imat1a = reshape(imat1a',[],1);
imat2a = reshape(imat2a',[],1);
imat3a = reshape(imat3a',[],1);
imat4a = reshape(imat4a',[],1);
imat1b = reshape(imat1b',[],1);
imat2b = reshape(imat2b',[],1);
imat3b = reshape(imat3b',[],1);
imat4b = reshape(imat4b',[],1);
Bmat1a = Bmat1a*filler;
Bmat2a = Bmat2a*filler;
Bmat3a = Bmat3a*filler;
Bmat4a = Bmat4a*filler;
Bmat1b = Bmat1b*filler;
Bmat2b = Bmat2b*filler;
Bmat3b = Bmat3b*filler;
Bmat4b = Bmat4b*filler;
Bmat1a = reshape(Bmat1a',[],1);
Bmat2a = reshape(Bmat2a',[],1);
Bmat3a = reshape(Bmat3a',[],1);
Bmat4a = reshape(Bmat4a',[],1);
Bmat1b = reshape(Bmat1b',[],1);
Bmat2b = reshape(Bmat2b',[],1);
Bmat3b = reshape(Bmat3b',[],1);
Bmat4b = reshape(Bmat4b',[],1);

%detector 3 dark counts
% pfD3 = [pfD3a pfD2b 1-pfD3a-pfD3b];
dc = ones(length(imat1a),1);

D3a = [1 0 0];
D3b = [0 1 0];
filler = [1 1 1];

d3a = dc*D3a;
d3b = dc*D3b;
d3a = reshape(d3a',[],1);
d3b = reshape(d3b',[],1);

clear D3a D3b

d1a = d1a*filler;
d1b = d1b*filler;
d1a = reshape(d1a',[],1);
d1b = reshape(d1b',[],1);
d2a = d2a*filler;
d2b = d2b*filler;
d2a = reshape(d2a',[],1);
d2b = reshape(d2b',[],1);

imat1a = imat1a*filler;
imat2a = imat2a*filler;
imat3a = imat3a*filler;
imat4a = imat4a*filler;
imat1b = imat1b*filler;
imat2b = imat2b*filler;
imat3b = imat3b*filler;
imat4b = imat4b*filler;
imat1a = reshape(imat1a',[],1);
imat2a = reshape(imat2a',[],1);
imat3a = reshape(imat3a',[],1);
imat4a = reshape(imat4a',[],1);
imat1b = reshape(imat1b',[],1);
imat2b = reshape(imat2b',[],1);
imat3b = reshape(imat3b',[],1);
imat4b = reshape(imat4b',[],1);
Bmat1a = Bmat1a*filler;
Bmat2a = Bmat2a*filler;
Bmat3a = Bmat3a*filler;
Bmat4a = Bmat4a*filler;
Bmat1b = Bmat1b*filler;
Bmat2b = Bmat2b*filler;
Bmat3b = Bmat3b*filler;
Bmat4b = Bmat4b*filler;
Bmat1a = reshape(Bmat1a',[],1);
Bmat2a = reshape(Bmat2a',[],1);
Bmat3a = reshape(Bmat3a',[],1);
Bmat4a = reshape(Bmat4a',[],1);
Bmat1b = reshape(Bmat1b',[],1);
Bmat2b = reshape(Bmat2b',[],1);
Bmat3b = reshape(Bmat3b',[],1);
Bmat4b = reshape(Bmat4b',[],1);

%detector 4 dark counts
% pfD4 = [pfD4a pfD4b 1-pfD4a-pfD4b];
dc = ones(length(imat1a),1);

D4a = [1 0 0];
D4b = [0 1 0];
filler = [1 1 1];

d4a = dc*D4a;
d4b = dc*D4b;
d4a = reshape(d4a',[],1);
d4b = reshape(d4b',[],1);

clear D4a D4b

d1a = d1a*filler;
d1b = d1b*filler;
d1a = reshape(d1a',[],1);
d1b = reshape(d1b',[],1);
d2a = d2a*filler;
d2b = d2b*filler;
d2a = reshape(d2a',[],1);
d2b = reshape(d2b',[],1);
d3a = d3a*filler;
d3b = d3b*filler;
d3a = reshape(d3a',[],1);
d3b = reshape(d3b',[],1);

%final packaging of observable arrays for E and B
imat1a = imat1a*filler;
imat2a = imat2a*filler;
imat3a = imat3a*filler;
imat4a = imat4a*filler;
imat1b = imat1b*filler;
imat2b = imat2b*filler;
imat3b = imat3b*filler;
imat4b = imat4b*filler;
imat1a = reshape(imat1a',[],1);
imat2a = reshape(imat2a',[],1);
imat3a = reshape(imat3a',[],1);
imat4a = reshape(imat4a',[],1);
imat1b = reshape(imat1b',[],1);
imat2b = reshape(imat2b',[],1);
imat3b = reshape(imat3b',[],1);
imat4b = reshape(imat4b',[],1);
Bmat1a = Bmat1a*filler;
Bmat2a = Bmat2a*filler;
Bmat3a = Bmat3a*filler;
Bmat4a = Bmat4a*filler;
Bmat1b = Bmat1b*filler;
Bmat2b = Bmat2b*filler;
Bmat3b = Bmat3b*filler;
Bmat4b = Bmat4b*filler;
Bmat1a = reshape(Bmat1a',[],1);
Bmat2a = reshape(Bmat2a',[],1);
Bmat3a = reshape(Bmat3a',[],1);
Bmat4a = reshape(Bmat4a',[],1);
Bmat1b = reshape(Bmat1b',[],1);
Bmat2b = reshape(Bmat2b',[],1);
Bmat3b = reshape(Bmat3b',[],1);
Bmat4b = reshape(Bmat4b',[],1);

%combine imatxx and Bxx and dxx
leng = length(imat1a);
c1a = ones(leng,1);
c2a = c1a;
c3a = c1a;
c4a = c1a;
c1b = c1a;
c2b = c1a;
c3b = c1a;
c4b = c1a;

%   = struct * Ecounts * Bcounts * Dcounts
c1a = c1a.*imat1a+c1a.*Bmat1a+c1a.*d1a;
c2a = c2a.*imat2a+c2a.*Bmat2a+c2a.*d2a;
c3a = c3a.*imat3a+c3a.*Bmat3a+c3a.*d3a;
c4a = c4a.*imat4a+c4a.*Bmat4a+c4a.*d4a;
c1b = c1b.*imat1b+c1b.*Bmat1b+c1b.*d1b;
c2b = c2b.*imat2b+c2b.*Bmat2b+c2b.*d2b;
c3b = c3b.*imat3b+c3b.*Bmat3b+c3b.*d3b;
c4b = c4b.*imat4b+c4b.*Bmat4b+c4b.*d4b;

% apply effects of lockout
for i = 1:leng
    if c1a(i)>0
        c1b(i)=0;
    end
    if c2a(i)>0
        c2b(i)=0;
    end
    if c3a(i)>0
        c3b(i)=0;
    end
    if c4a(i)>0
        c4b(i)=0;
    end
end

%sort observables

Mob = horzcat(c1a,c2a,c3a,c4a,c1b,c2b,c3b,c4b);

Mob1 = Mob.*0;
for ic = 1:size(Mob,1)
    for jc = 1:size(Mob,2)
        if Mob(ic,jc)>=1
            Mob1(ic,jc)=1;
        end
    end
end

% preallocate vars in memory for speed
MobS = sum(Mob1,2);

key1a=zeros(leng,1);
key2a=zeros(leng,1);
key3a=zeros(leng,1);
key4a=zeros(leng,1);
key1b=zeros(leng,1);
key2b=zeros(leng,1);
key3b=zeros(leng,1);
key4b=zeros(leng,1);
keyMulti=zeros(leng,1);
keyNone=zeros(leng,1);

for i = 1:leng
    if MobS(i)==1
        if Mob1(i,1)==1
            key1a(i) = 1;
        end
        if Mob1(i,2)==1
            key2a(i) = 1;
        end
        if Mob1(i,3)==1
            key3a(i) = 1;
        end
        if Mob1(i,4)==1
            key4a(i) = 1;
        end
        if Mob1(i,5)==1
            key1b(i) = 1;
        end
        if Mob1(i,6)==1
            key2b(i) = 1;
        end
        if Mob1(i,7)==1
            key3b(i) = 1;
        end
        if Mob1(i,8)==1
            key4b(i) = 1;
        end
    end
    if MobS(i)>1
            keyMulti(i) = 1;
    end
    if MobS(i)==0
        keyNone(i) = 1;
    end
end


if ispc
    fprintf('Saving sort key: SingleEmitter to ...\\Library...\n')
save("Library\sortingP_idxGuide_1emitter_backgrd_dc.mat",'key1a', 'key2a', 'key3a', 'key4a', ...
    'key1b', 'key2b', 'key3b', 'key4b', "keyNone","keyMulti")
elseif isunix || ismac
    fprintf('Saving sort key: SingleEmitter to .../Library...\n')
save("Library/sortingP_idxGuide_1emitter_backgrd_dc.mat",'key1a', 'key2a', 'key3a', 'key4a', ...
    'key1b', 'key2b', 'key3b', 'key4b', "keyNone","keyMulti")

end

