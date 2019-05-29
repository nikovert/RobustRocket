%% Generate a linearized model of the system and save it as a plant
clear; clc; close all;
% Create a model of the rocket post main engine cutoff
rocket = CreateRocket("rocketConfiguration.txt");

%% Create Plant model
plant_fig = figure('Name', 'Plant','NumberTitle','off'); clf; grid on; legend;

%x = [height; velocity]; 
x0 = [2000; 160];
u0 = 0.1;

%We put uncertainty on our parameter when creating our linear model
[A,B,K] = simple_linear_model_uss(rocket, x0, u0);

C = [1 0; A(2,:)]; % We can only measure the height and the acceleration 

% We enable using the offset by adding a new state x_hat = [x;1]
%   dX/dt = [A K; 0 0] X + [B; 0] u
%   y = [C 0] X + Du
A_hat = [A K; zeros(1,length(A)+1)];
B_hat = [B;0];
C_hat = [1 0 0; A_hat(2,:)];
D_hat = 0;

%Create a continous state-space model
sys = ss(A_hat,B_hat,C_hat,D_hat,'StateName',{'Height' 'Velocity' 'offset'},'InputName','AirbrakeExtensionPrecentile', 'OutputName', {'Height', 'Acceleration'});

%Generate Latex version of transfer function
l_G = tf_to_latex(sys);

% Check stabilizability https://ch.mathworks.com/help/control/ref/ctrbf.html
[Abar,Bbar,Cbar,T,k] = ctrbf(sys.A.NominalValue,sys.B.NominalValue,sys.C.NominalValue);
size_Auc = length(Abar)-sum(k);
Auc = Abar(size_Auc);
if isstable(Auc)
    disp('Plant is stabilizable');
else
    disp('Plant is not stabilizable');
end

%interesting Plot of the initial condition:
figure(plant_fig);
initial(sys,[x0;1], 50, 'b');

G = minreal(sys.Nominal);

%% Define Uncertainties (Slide 9.13)
Uncertainties_fig = figure('Name', 'Uncertainties','NumberTitle','off'); clf; hold on; grid on; legend;
Performance_fig = figure('Name', 'Performance Weights','NumberTitle','off'); clf; hold on; grid on; legend;

inaccuracies_baro_fig = figure('Name', 'Barometric Uncertainties','NumberTitle','off'); clf; hold on; grid on; legend;

inaccuracies_acc_fig = figure('Name', 'Accelerometer Uncertainties','NumberTitle','off'); clf; hold on; grid on; legend;

%We now define our input uncertainty. At low frequencies, our
%motor can accurately actuate the airbrakes, at high frequencies though our
%slew rate becomes noticable and accurately additionally will naturally go down. 
wc = 30; % At about 30Hz our motor cannot accuratly control the airbrakes anymore and we want Wu to be at 50%
Wu = makeweight(0.8,wc,1.7);
figure(Uncertainties_fig);
bodemag(Wu);

l_wu = tf_to_latex(Wu);

%We will add white gausian noise for our sensors.
%Additionally the sensors will will suffer inaccuracies when sampled at high frequencies (above 20Hz), especially the barometer. 
W_sensorNoise_acc = zpk(0.5);
W_sensorNoise_baro = zpk(0.7);
W_sensorNoise = [W_sensorNoise_baro 0; 0 W_sensorNoise_acc];
figure(Uncertainties_fig);
bodemag(W_sensorNoise_acc, '--',W_sensorNoise_baro,'--');

% https://ch.mathworks.com/help/robust/ref/ucover.html#f10-789594
%Using usample we find the correct uncertainties for our Plant
G1 = tf(sys(1,1));
G2 = tf(sys(2,1));

samples1 = usample(sys(1,1),100);
samples2 = usample(sys(2,1),100);
[~,Info1] = ucover(samples1, G1);
[~,Info2] = ucover(samples2, G2);

relerr1 = (G1-samples1)/G1;
relerr2 = (G2-samples2)/G2;

figure(inaccuracies_baro_fig);
bodemag(relerr1, Info1.W1, '-r');
legend 'relative error to nominal plant', 'upper bound';

figure(inaccuracies_acc_fig);
bodemag(relerr2, Info2.W1, '-r');
legend 'relative error to nominal plant', 'upper bound';

W_inaccuracies_baro = Info1.W1;
W_inaccuracies_acc = Info2.W1;
W_inaccuracies = [W_inaccuracies_baro 0; 0 W_inaccuracies_acc];
figure(Uncertainties_fig);
bodemag(W_inaccuracies_baro);
bodemag(W_inaccuracies_acc);

%Performance weight
% We would like our controller to perform well, up to 30HZ
w_c = 35;
Wp = 0.99 * makeweight(1.001,w_c+20, 0.3) * makeweight(0.001,w_c-15, 1.01);
Wp_h = 0.8 * Wp;
Wp_v = Wp;
figure(Performance_fig);
bodemag(Wp_h, Wp_v);

save('extended_linearisedPlant_workspace');  

%% Design a stabilizing controller for the nominal system
%func  = @(k1, k2) max(real(eig(G.A - G.B * [k1 k2] * G.C)));
%lambda = - ones(201) * inf;
k_opt = ss([-1, -1]); % silly initial guess
stab_controller = false;
if stab_controller
    percentChange = 100 / 201^2;
    percentage =  19.9005;
    fprintf('\n%.2f%%', percentage);
    for x = -60:100
        for y = -100:100
            K_stab = ss([x, y]);
            [A_stab,B_stab,C_stab,D_stab] = linmod('FlightModel_stab');
            G_stab = ss(A_stab,B_stab,C_stab,D_stab);
            if isstable(G_stab)
                lambda(101+x,101+y) = max(eig(G_stab));
                if(norm(k_opt.D) > norm(K_stab.D))
                    k_opt = K_stab;
                end
            end
            percentage = percentage + percentChange;
            if percentage < 10
                fprintf('\b\b\b\b\b%.2f%%', percentage)
            else
                fprintf('\b\b\b\b\b\b%.2f%%', percentage)
            end
        end
    end
    fprintf('\b\b\b\b\b\b%.2f%%', 100)
    if max(max(lambda)) > -inf
        disp('Found stabilizing controller');
    else
        disp('No stabilizing controller found');
    end
end
%[X,Y] = meshgrid(-100:100,-100:100);
%surf(X,Y,lambda_hat)

%Choose a stabilizing controller
% K_stab = ss([-60, 60]);
% G_stab = lft(G,-K_stab);
% isstable(G_stab)

% K_stab = -K_stab;
% [A_stab,B_stab,C_stab,D_stab] = linmod('FlightModel_stab');
% G_stab = ss(A_stab,B_stab,C_stab,D_stab)
% isstable(G_stab)


%Plot nominal Plant (Slide 9.11)

% Create unceratinty model using the stabilizing controller
K_stab = k_opt;
[A_stab,B_stab,C_stab,D_stab] = linmod('FlightModel_stab');
G_stab = ss(A_stab,B_stab,C_stab,D_stab);

figure(plant_fig);
%H(i,j) selects the response from the jth input to the ith output.
G_stab = minreal(G_stab([6:7],[8])); % only concider output without the reference output
initial(G_stab,[x0;1], 50, 'r');

save('extended_linearisedPlant_workspace');  

%% H_inf controller design (Slide 9.14)
tracking_fig = figure('Name', 'Tracking Behaviour','NumberTitle','off'); clf; hold on; grid on; legend;
singularvalue_fig = figure('Name', 'Singular values of the closed loop system','NumberTitle','off'); clf; hold on; grid on; legend;
% Modify the system to allow reference to go to our ouput
% A_new = sys.A;
% B_new = [sys.B zeros(length(sys.B),2)];
% C_new = [sys.C; zeros(2,length(sys.C))];
% D_new = [sys.D zeros(2, length(sys.D));  zeros(length(sys.D),1) eye(2)];
% 
% sys_new = uss(A_new, B_new, C_new, D_new);
% 
% [M,DELTA,BLKSTRUCT,NORMUNC] = lftdata(sys_new);

% G = minreal(sys.NominalValue);
[A_P,B_P,C_P,D_P] = linmod('FlightModel');
P = ss(A_P,B_P,C_P,D_P);
%P_stab = ss(A_stab,B_stab,C_stab,D_stab);

%           _____
%          |Delta|
%      |-->|_____|>--|
%    z |             | v
%      |   _______   |  
%      |--|       |<-|
%  e <----|  sys  |<---- w
%      |--|_______|<-|
%      |             |
%    y |    _____    | u
%      |-->|  K  |---|
%          |_____|
%

Iz = [1:3]'; % Create indices for each block.
    nz = length(Iz);
Iv = [1:3]';
    nv = length(Iv);
    
Ie = [4:5]';
    ne = length(Ie);
Iw = [4:7]';
    nw = length(Iw);

Iy = [6:9]';
    ny = length(Iy);
Iu = [8]';
    nu = length(Iu);

Pnomdesign = P([Ie;Iy],[Iw;Iu]); % select [e;y] <- [w;u]
%Pnomdesign_stab = P_stab([Ie;Iy],[Iw;Iu]);

Pnomdesign = minreal(Pnomdesign); % remove non-minimal states.
%Pnomdesign_stab = minreal(Pnomdesign_stab);

nmeas = 4;
nctrl = 1;

[Knom,Gnom,gamma,info] = hinfsyn(Pnomdesign,nmeas,nctrl,...
                                'METHOD','ric',... % Riccati solution
                                'TOLGAM',0.1); % gamma tolerance
figure(singularvalue_fig);
sigma(Gnom,ss(gamma));
legend("Largest singular values of the closed loop system", "upper bound gamma");
ylim([-20,5]);
                            
if isstable(Gnom)
   disp('close loop h_infinity system is stable')
end
disp(['Gamma: ', num2str(gamma)]);
%[Knom_stab,Gnom_stab,gamma_stab,info_stab] = hinfsyn(Pnomdesign_stab,nmeas,nctrl,'METHOD','ric','TOLGAM',0.1); 

% closed-loop analysis and simulation (Slide 9.15)
[Aclp,Bclp,Cclp,Dclp] = linmod('ControlModel');
Gclp = ss(Aclp,Bclp,Cclp,Dclp);
Gclp_nom = Gclp(Ie,Iw);
Gclp_nom = minreal(Gclp_nom); % remove perturbation weight states.

% Plot nominal command responses (Slide 9.16)
Gclp_nom_tracking =  Gclp_nom([1 2]',[3 4]');
figure(tracking_fig);
step(Gclp_nom_tracking(1,1),[0,50]);
step(Gclp_nom_tracking(2,2),[0,50]); legend 'h_ref->h' 'v_ref->v';
% Plot nominal noise responses (Slide 9.17)
% Plot nominal disturbance responses (Slide 9.18)
% Plot nominal step response (Slide 9.19)
save('extended_linearisedPlant_workspace');   
%% Closed-loop robustness analysis (Slide 9.22)
mu_fig_0 = figure('Name', 'Mu_0','NumberTitle','off'); clf; hold on; grid on; legend;
wc_fig = figure('Name', 'Worst Case','NumberTitle','off'); clf; hold on; grid on; legend;
Grob = lft(P,Knom);

% Compute and plot worst-case gain
figure(wc_fig)
wcsigma(Grob)
axis([1 1000 -20 10]);
grid on;

%Grob_stab = lft(P_stab,Knom_stab);

omega = logspace(-1,2,250); % frequency vector
Grob_f = frd(Grob,omega); % frequency response
%Grob_f_stab  = frd(Grob_stab ,omega);

RS_blk = [1 0; 1 0; 1 0]; % {d1, d2 ,d3} are real
muRS = mussv(Grob_f(Iz,Iv),RS_blk);
%muRS_stab = mussv(Grob_f_stab(Iz,Iv),RS_blk);

muNP = svd(Grob_f(Ie,Iw(3:4))); % nominal performance
%muNP_stab = svd(Grob_f_stab(Ie,Iw));

RP_blk = [RS_blk; 4 2]; % {d1, d2 ,d3, D4} are real, D4 is 4x2 real matrix
[muRP,muinfo0] = mussv(Grob_f,RP_blk); % robust performance
%[muRP_stab,muinfo0_stab] = mussv(Grob_f_stab,RP_blk);

% Plot mu analysis (Slide 9.23)
figure(mu_fig_0);
bodemag(muRS, muNP, muRP);
legend 'muRS' 'muNP' 'mRP'; 
%bodemag(muRS_stab, muNP_stab, muRP_stab);
%legend 'muRS_stab' 'muNP_stab' 'mRP_stab';
%% Start D-K iteration
mu_fig_1 = figure('Name', 'Mu_1','NumberTitle','off'); clf; hold on; grid on; legend;
wc_fig_dk = figure('Name', 'Worst Case after D-K iteration','NumberTitle','off'); clf; hold on; grid on; legend;
dk_iter = figure('Name', 'mu Bounds during D-K iteration','NumberTitle','off'); clf; hold on; grid on; legend;
% Normalized error dynamics
delta1 = ultidyn('delta1',[1 1]);
delta2 = ultidyn('delta2',[1 1]);
delta3 = ultidyn('delta3',[1 1]);

delta = blkdiag(delta1, delta2, delta3);
delta.InputName = {'z'};
delta.OutputName = {'v'};

% Generalized uncertain plant P
systemnames = 'P delta';
inputvar = '[w{4};u]';
outputvar = '[P(4:9)]';
input_to_P = '[delta(1:3);w(1:4);u]';
input_to_delta = '[P(1:3)]';
%cleanupsysic = 'yes';
P_uss = sysic;

[Kmu1,clpmu,bnd, dkinfo] = dksyn(P_uss,nmeas,nctrl);

figure(dk_iter);
bodemag(muRP, dkinfo{1}.MussvBnds)
legend("inital mu bound", "mu bound after first iteration");


Gmu1 = lft(P_uss,Kmu1); % repeat the robustness analysis
figure(wc_fig_dk);
wcsigma(Gmu1,{1,1000})
axis([1 1000 -15 5]); grid on;

Gmu1 = lft(P,Kmu1);
Gmu1_f = frd(Gmu1,omega);

muRS1 = mussv(Gmu1_f(Iz,Iv),RS_blk);
muNP1 = svd(Gmu1_f(Ie,Iw));
[muRP1,muinfo1] = mussv(Gmu1_f,RP_blk);

figure(mu_fig_1);
bodemag(muRS1, muNP1, muRP1);
%%
keyboard
%% Manual DK iteration
dk_raw_fig = figure('Name', 'Raw D-scales','NumberTitle','off'); clf; hold on; grid on; legend;
dk_normalized_fig = figure('Name', 'Normalized D-scales','NumberTitle','off'); clf; hold on; grid on; legend;
D0_1_fig = figure('Name', 'D0_1 fitting','NumberTitle','off'); clf; hold on; grid on; legend;
D0_2_fig = figure('Name', 'D0_2 fitting','NumberTitle','off'); clf; hold on; grid on; legend;
D0_3_fig = figure('Name', 'D0_3 fitting','NumberTitle','off'); clf; hold on; grid on; legend;
D0_4_fig = figure('Name', 'D0_4 fitting','NumberTitle','off'); clf; hold on; grid on; legend;

% nmeas = 4;
% nctrl = 1;
% [k,clp,bnd] = dksyn(sys,nmeas,nctrl);

[Dl0,Dr0] = mussvunwrap(muinfo0); % extract D-scales
D0_perf = Dl0(5,5);
D0_1 = Dl0(1,1)/D0_perf; % normalize w.r.t. performance D-scale
D0_2 = Dl0(2,2)/D0_perf;
D0_3 = Dl0(3,3)/D0_perf;
D0_4 = Dl0(4,4)/D0_perf;

% Plot raw D scalings (Slide 9.25)
figure(dk_raw_fig);
bodemag(Dl0(1,1), Dl0(2,2), Dl0(3,3), Dl0(4,4));
legend Dl0(1,1) Dl0(2,2) Dl0(3,3) Dl0(4,4)

% Plot normalized D scalings (Slide 9.27) 
figure(dk_normalized_fig);
bodemag(D0_1, D0_2, D0_3, D0_4);

% Fitting D-scal (Slide 9.28)
D0_1a = fitfrd(genphase(D0_1),0); % 0th order fit
D0_1b = fitfrd(genphase(D0_1),1); % 1st order fit
D0_1c = fitfrd(genphase(D0_1),2); % 2nd order fit
D0_1d = fitfrd(genphase(D0_1),3); % 3rd order fit
figure(D0_1_fig);
bodemag(D0_1,D0_1a, D0_1b, D0_1c, D0_1d);
% Looking at the plots we choose D0_1d

D0_2a = fitfrd(genphase(D0_2),0); % 0th order fit
D0_2b = fitfrd(genphase(D0_2),1); % 1st order fit
D0_2c = fitfrd(genphase(D0_2),2); % 2nd order fit
D0_2d = fitfrd(genphase(D0_2),3); % 3rd order fit
figure(D0_2_fig);
bodemag(D0_2,D0_2a, D0_2b, D0_2c, D0_2d);
% Looking at the plots we choose D0_2d

D0_3a = fitfrd(genphase(D0_3),0); % 0th order fit
D0_3b = fitfrd(genphase(D0_3),1); % 1st order fit
D0_3c = fitfrd(genphase(D0_3),2); % 2nd order fit
D0_3d = fitfrd(genphase(D0_3),3); % 3rd order fit
figure(D0_3_fig);
bodemag(D0_3,D0_3a, D0_3b, D0_3c, D0_3d);
% Looking at the plots we choose D0_3c
D0_4a = fitfrd(genphase(D0_4),0); % 0th order fit
D0_4b = fitfrd(genphase(D0_4),1); % 1st order fit
D0_4c = fitfrd(genphase(D0_4),2); % 2nd order fit
D0_4d = fitfrd(genphase(D0_4),3); % 3rd order fit
figure(D0_4_fig);
bodemag(D0_4,D0_4a, D0_4b, D0_4c, D0_4d);
% Looking at the plots we choose D0_4a
% Plot fitting D-scal (Slide 9.29)

sysD0 = tf(blkdiag(D0_1d,D0_2d,D0_3c));
sysDi0 = inv(tf(sysD0));

%% Design iteration #1:
Pmu1design =    [sysD0, zeros(nz,ne+nmeas);
                    zeros(ne,nz), eye(ne,ne), zeros(ne,nmeas);
                    zeros(nmeas,nz+ne),eye(nmeas,nmeas)] ...
            * P * ...
                [sysDi0,zeros(nv,nw+nctrl);
                    zeros(nw,nv), eye(nw,nw), zeros(nw,nctrl);
                    zeros(nctrl,nv+nw),eye(nctrl,nctrl)];
                
[Kmu1,Gmu1,gamma1,info1] = hinfsyn(Pmu1design,nmeas,nctrl,...
                                    'METHOD','lmi',... % LMI solution
                                    'TOLGAM',0.1); % gamma tolerance
                    
Gmu1 = lft(P,Kmu1); % repeat the robustness analysis

% Compute and plot worst-case gain
figure(wc_fig_dk);
wcsigma(Gmu1)
axis([1e-1 1000 -20 10])

Gmu1_f = frd(Gmu1,omega);

muRS1 = mussv(Gmu1_f(Iz,Iv),RS_blk);
muNP1 = svd(Gmu1_f(Ie,Iw));
[muRP1,muinfo1] = mussv(Gmu1_f,RP_blk);