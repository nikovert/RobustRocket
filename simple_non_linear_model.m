function dx = simple_non_linear_model(rocket, x,u)
%returns a simplified non linear model of the system based on a
%predetermined atmosphere model.

Mass = rocket.Mass_dry;
A_ref = rocket.A_ref;

Re = norm(x(2))*rocket.CS_Position/13.164e-6; %kinematic viscosity of air
h = rocket.CS_Length*sin(rocket.CS_Angle); %height of the air brake
delta = 0.37*rocket.CS_Position/(Re^0.25); %height of the boundary layer
CD_theta_corrected = rocket.CS_CD0*(1-0.25*delta/h)*sin(rocket.CS_Angle)^3;

v = x(2)/352.2912; % speed in Mach
M            = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
CS_extension = [0, 0.3364, 0.6527, 1]; %percentage of the extension 
[Xval,Yval]  = meshgrid(M,CS_extension);

% 0%
Cd_values_0 = [0.041, 0.07,  0.106, 0.148, 0.209, 0.258, 0.384];
% 33%
Cd_values_1 = [0.044, 0.076, 0.11,  0.166, 0.2256, 0.296, 0.42];
% 66%
Cd_values_2 = [0.041, 0.076, 0.108, 0.155, 0.21,   0.28,  0.38];
% 100%
Cd_values_3 = [0.048, 0.085, 0.132, 0.19,  0.261,  0.347, 0.46];
% Cd table
Cd_table    = [Cd_values_0;Cd_values_1; Cd_values_2; Cd_values_3];

%Cd at specified velocity and airbrake extension 
Cd_inital = interp2(Xval,Yval,Cd_table,v,u,'linear');

Cd = Cd_inital + u * rocket.CS_Area/rocket.A_ref*CD_theta_corrected*3;



%We use environement(1400, 40, 86000, rocket) as constant
g = 9.7920;     %e.g
rho = 1.1407;    %e.rho

Fa = -g - 0.5/Mass*rho*A_ref*Cd*x(2).^2;
dx = [x(2); Fa];

end