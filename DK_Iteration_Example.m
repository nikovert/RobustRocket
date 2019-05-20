% D-K Iteration Example
clear all
clc
G=tf(1,[75 1])*[87.8 -86.4;108.2 -109.6];
% Convert into state-space form
G=minreal(ss(G));
% Determine the performance and uncertainty weights in 2 dimensions
Wp=minreal(ss(tf([0.5 0.05],[1 1e-6])*eye(2)));
Wi=minreal(ss(tf([1 0.2],[0.5 1])*eye(2)));
% Form the generalised plant P by calling sysic-function
systemnames='G Wp Wi'; % defines the models used
inputvar='[udelta(2);w(2);u(2)]'; % defines the sizes and names of the input signals
outputvar='[Wi;Wp; -G-w]'; % defines the names of the output signals
input_to_G='[u+udelta]'; % defines the signals fed to the process model (G)
input_to_Wp='[G+w]'; % defines the signals fed to the performance weighting model
input_to_Wi='[u]'; % defines the signals fed to the input uncertainty weighting model
sysoutname='P'; % defines the resulting model name
cleanupsysic='yes'; % toggles whether the above variables are cleared after calling sysic-function
sysic; % Generates a plant model
P=minreal(ss(P));
%Initialize
omega=logspace(-3,3,61); % spaces 61 divided logarithmically
%over 10^-3 and 10^3
blk=[1 1;1 1;2 2]; % Determine the block structure (2 times 1by1
%block or input uncertainty and 1 times 2by2 for performance
nmeas=2; nu=2; d0=1;
D=append(d0,d0,tf(eye(2)),tf(eye(2))); % Initial scaling
% Start iteration

% STEP 1: Find H-infinity controller with the given scalings
[K,Nsc,gamma,info]=hinfsyn(D*P*inv(D),nmeas,nu,'method','lmi','Tolgam',1e-3);
Nf=frd(lft(P,K),omega);

% STEP 2: Compute mu using upper bound;
[mubnds,Info]=mussv(Nf,blk,'c');
bodemag(mubnds(1,1),omega);
hold on
murp=norm(mubnds(1,1),inf,1e-6);

%STEP 3: Fit resulting D-scales
[dsys1,dsysr]=mussvunwrap(Info);
dsys1=dsys1/dsys1(3,3);
d1=fitfrd(genphase(dsys1(1,1)),4);

% GO TO STEP 1
% Alternatively use automatic software
Delta=[ultidyn('D_1',[1 1]) 0;0 ultidyn('D_2',[1 1])] % Diagonal uncertainty
Punc=lft(Delta,P);
opt=dkitopt('FrequencyVector',omega);
[K,clp,bnd,dkinfo]=dksyn(Punc,nmeas,nu,opt);
for in=1:4
    bodemag(dkinfo{1,in}.MussvBnds(1,1),omega)
end
hold off