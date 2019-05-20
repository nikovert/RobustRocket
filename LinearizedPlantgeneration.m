%% Generate a linearized model of the system and save it as a plant
clear; clc; close all;
addpath("ARIS_files");
global env

%Initialize the rocket
rocket_model = rocket(init_rocket("HEIDI_Final.txt"));

%Add the motor values
motor_init(rocket_model);

%Create environment model
env = environment(1400, 40, 86000, rocket_model);

% Braking params
rocket_model.B_BRAKING_MPC = false;
rocket_model.B_BRAKING = true;
u0 = 0.097353;
rocket_model.hard_set_u(u0);

dt = 0.01;
tend = 10;
[t, state] = accent_calc(rocket_model,tend,dt);

figure(1);
hold on;
plot(state(:,10)/rocket_model.Mass,state(:,3))
xlabel('v(kg * m/s)')
ylabel('Height(m)')
title('State Trajectory');
grid on;

%% Create Plant model
figure(2); clf;
grid on; legend;

height = state(end,3);
velocity = state(end,10)/rocket_model.Mass;
x = [height; velocity]; 
dx = simple_non_linear_model(rocket_model, x);
[A,B,K] = simple_linear_model(rocket_model, x, dx);
C = [1 0; A(2,:)]; % We can only measure the height and the acceleration 
C_inv = inv(C);

% We integrate our by adding a new state x_hat = [x;1]
%   dX/dt = [A K; 0 0] X + [B; 0] u
%   y = [C 0] X + Du
A_hat = [A K; zeros(1,length(A)+1)];
B_hat = [B;0];
C_hat = [1 0 0; A_hat(2,:)];
D_hat = 0;

%Create a continous state-space model
sys = ss(A_hat,B_hat,C_hat,D_hat,'StateName',{'Height' 'Velocity' 'offset'},'InputName','AirbrakeExtensionPrecentile', 'OutputName', {'Height', 'Acceleration'});

%Create the corresponding transfer function
G = tf(sys);

%Plot nominal Plant (Slide 9.11)
bode(G);
grid on;

save('nominal_linear_plant', 'G');

%% Define Uncertainties (Slide 9.13 )
figure(3); clf;
hold on; grid on; legend;

%We now define our input uncertainty. At low frequencies, our
%motor can accurately actuate the airbrakes, at high frequencies though our
%slew rate becomes noticable and accurately additionally will naturally go down. 
w_c = 100; % At about 100Hz our motor cannot accuratly control the airbrakes anymore and we want Wu to be at 50%
Wu = 0.5*tf([1 0],[1 w_c]); 
InputUnc = ultidyn('InputUnc',1);
bode(Wu);

%We will add white gausian noise for our sensors.
%Additionally the sensors will will suffer inaccuracies when sampled at high frequencies (above 20Hz), especially the barometer. 
W_sensorNoise_acc = 0.05 * tf(1);
W_sensorNoise_baro = 0.10 * tf(1);
W_sensorNoise = [W_sensorNoise_baro 0; 0 W_sensorNoise_acc];
SensorNoise = ultidyn('SensorNoise',1);
bode(W_sensorNoise_acc, '--',W_sensorNoise_baro,'--'); 

w_c_baro = 20;
w_c_acc = 100;
W_inaccuracies_acc = 2 * tf([1 1 0],[1 1 w_c_baro]); 
W_inaccuracies_baro = 2 * tf([1 0 0],[1 1 w_c_acc]); 
W_inaccuracies = [W_inaccuracies_baro 0; 0 W_inaccuracies_baro];
SensorUnc = ultidyn('SensorUnc',1);
bode(W_inaccuracies_acc, W_inaccuracies_baro);

%Performance weight
% We would like our controller to perform well, up to 100HZ
w_c = 100;
Wp = 150 * tf(1,[1 w_c]);
Wp_h = Wp;
Wp_v = Wp;
bode(Wp);

save('linearisedPlant_workspace');  
%% H_inf controller design (Slide 9.14)
figure(4); clf;
hold on; grid on; legend;

G = minreal(ss(G));
% est = tf([1],[1 0]);
% % Generalized plant P
% systemnames = 'Wu G est W_inaccuracies_baro W_inaccuracies_acc W_sensorNoise_baro W_sensorNoise_acc Wp_h Wp_v';
% inputvar = '[u_cmd; v_ref; h_ref; v1; v2; v3; v4; v5]';
% outputvar = '[Wu(1); G(1); G(2); Wp_h(1); Wp_v(1)]';
% input_to_Wu = '[u_cmd]';
% input_to_G = '[u_cmd + v1]';
% input_to_W_inaccuracies_baro = '[v2]';
% input_to_W_inaccuracies_acc = '[v3]';
% input_to_W_sensorNoise_baro = '[v4]';
% input_to_W_sensorNoise_acc = '[v5]';
% input_to_est = '[G(2)+W_inaccuracies_acc(1)+W_sensorNoise_acc(1)]';
% input_to_Wp_h = '[h_ref-G(1)-W_inaccuracies_baro(1)-W_sensorNoise_baro(1)]';
% input_to_Wp_v = '[v_ref-est(1)]';
% P = sysic;
% P = ss(P);

[A_P,B_P,C_P,D_P] = linmod('FlightModel');
P = ss(A_P,B_P,C_P,D_P);
 
Iz = [1:3]'; % Create indices for each block.
Iv = [1:3]';

Ie = [4:5]';
Iw = [4:7]';

Iy = [6:9]';
Iu = [8]';

Pnomdesign = P([Ie;Iy],[Iw;Iu]); % select [e;y] <- [w;u]
Pnomdesign = minreal(Pnomdesign); % remove non-minimal states.

% Using 'ultidyn' we can create our delta

nmeas = 4;
nctrl = 1;

[Knom,Gnom,gamma,info] = hinfsyn(Pnomdesign,nmeas,nctrl,...
                                'METHOD','ric',... % Riccati solution
                                'TOLGAM',0.1); % gamma tolerance
       
% closed-loop analysis and simulation (Slide 9.15)
[Aclp,Bclp,Cclp,Dclp] = linmod('ControlModel');
Gclp = ss(Aclp,Bclp,Cclp,Dclp);
Gclp_nom = Gclp(Ie,Iw);
Gclp_nom = minreal(Gclp_nom); % remove perturbation weight states.

% Plot nominal command responses (Slide 9.16)
% Plot nominal noise responses (Slide 9.17)
% Plot nominal disturbance responses (Slide 9.18)
% Plot nominal step response (Slide 9.19)

save('linearisedPlant_workspace');  
%% Closed-loop robustness analysis (Slide 9.22)
figure(5); clf;
hold on; grid on; legend;

Grob = lft(P,Knom);
omega = logspace(-4,2,250); % frequency vector
Grob_f = frd(Grob,omega); % frequency response

RS_blk = [1 0; 1 0; 1 0]; % {d1, d2 ,d3} are real
muRS = mussv(Grob_f(Iz,Iv),RS_blk);

muNP = svd(Grob_f(Ie,Iw)); % nominal performance

RP_blk = [RS_blk; 4 2]; % {d1, d2 ,d3, D4} are real, D4 is 4x2 real matrix
[muRP,muinfo0] = mussv(Grob_f,RP_blk); % robust performance

% Plot mu analysis (Slide 9.23)

%% Start D-K iteration
figure(6); clf;
hold on; grid on; legend;
figure(7); clf;
hold on; grid on; legend;

[Dl0,Dr0] = mussvunwrap(muinfo0); % extract D-scales
D0_perf = Dl0(5,5);
D0_1 = Dl0(1,1)/D0_perf; % normalize w.r.t. performance D-scale
D0_2 = Dl0(2,2)/D0_perf;
D0_3 = Dl0(3,3)/D0_perf;
D0_4 = Dl0(4,4)/D0_perf;

% Plot raw D scalings (Slide 9.25)
figure(6);
bode(Dl0(1,1), Dl0(2,2), Dl0(3,3), Dl0(4,4));
legend Dl0(1,1) Dl0(2,2) Dl0(3,3) Dl0(4,4)

% Plot normalized D scalings (Slide 9.27) 
figure(7);
bode(D0_1, D0_2, D0_3, D0_4);

% Fitting D-scal (Slide 9.28)
D0_1a = fitfrd(genphase(D0_1),0); % 0th order fit
D0_1b = fitfrd(genphase(D0_1),1); % 1st order fit
D0_1c = fitfrd(genphase(D0_1),2); % 2nd order fit
D0_1d = fitfrd(genphase(D0_1),3); % 3rd order fit

% Plot fitting D-scal (Slide 9.29)

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
Gmu1_f = frd(Gmu1,omega);

muRS1 = mussv(Gmu1_f(Iz,Iv),RS_blk);
muNP1 = svd(Gmu1_f(Ie,Iw));
[muRP1,muinfo1] = mussv(Gmu1_f,RP_blk);

