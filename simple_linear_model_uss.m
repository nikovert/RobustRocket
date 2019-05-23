function [A,B,K] = simple_linear_model_uss(r, x, u)

% Calculate dx using our non-linear model
dx = simple_non_linear_model(r, x, u);

%We use environement(1400, 40, 86000, rocket) as constant
Earth_M = 5.9724e+24; % e.Earth_M;
h_g = 1400; % e.h_g;
G = 6.6741e-11; % e.G;
Earth_R = 6378000; % e.Earth_R;
rho_g = 1.2250; % e.rho_g;
Temp_grad = 0.0065; % e.Temp_grad;
Temp_g = 313.1500; % e.Temp_g;
R = 287; %e.R;

%From the rocket
A_ref = r.A_ref;
Mass_dry = r.Mass_dry;
CS_Length = r.CS_Length;
CS_Angle = r.CS_Angle;
CS_Position = r.CS_Position;
CS_Area = r.CS_Area;
CS_CD0 = r.CS_CD0;
Xdot = x(2);

%We replace Cd_mandell_simple(Xdot) with:
%Cd values at different mach numbers and Airbrake extention
v = Xdot/352.2912; % speed in Mach
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
Cd = interp2(Xval,Yval,Cd_table,v,u,'linear');

%Using Mathematica we compute the derivatives
dFdz = 2*Earth_M*G*h_g*(Earth_R+h_g*x(1))^(-3)+(-1/2)*A_ref*Mass_dry^(-1) ...
  *rho_g*x(2)^2*(Temp_g^(-1)*(Temp_g+(-1)*Temp_grad*x(1)))^((-1)+ ...
  Earth_M*G*R^(-1)*Temp_grad^(-1)*(Earth_R+h_g*x(1))^(-2))*((-1)* ...
  Temp_grad*(Temp_g+(-1)*Temp_grad*x(1))^(-1)*((-1)+Earth_M*G*R^( ...
  -1)*Temp_grad^(-1)*(Earth_R+h_g*x(1))^(-2))+(-2)*Earth_M*G*h_g* ...
  R^(-1)*Temp_grad^(-1)*(Earth_R+h_g*x(1))^(-3)*log(Temp_g^(-1)*( ...
  Temp_g+(-1)*Temp_grad*x(1))))*(Cd+3*A_ref^(-1)*CS_Area* ...
  CS_CD0*u*(1+(-0.557171E-2)*CS_Length^(-1)*CS_Position*( ...
  CS_Position*x(2))^(-0.25E0)*csc(CS_Angle))*sin(CS_Angle)^3);

dFdv = (-0.208939E-2)*CS_Area*CS_CD0*CS_Length^(-1)*CS_Position^2* ...
  Mass_dry^(-1)*rho_g*u*x(2)^2*(CS_Position*x(2))^(-0.125E1)*( ...
  Temp_g^(-1)*(Temp_g+(-1)*Temp_grad*x(1)))^((-1)+Earth_M*G*R^(-1) ...
  *Temp_grad^(-1)*(Earth_R+h_g*x(1))^(-2))*sin(CS_Angle)^2+(-1)* ...
  A_ref*Mass_dry^(-1)*rho_g*x(2)*(Temp_g^(-1)*(Temp_g+(-1)* ...
  Temp_grad*x(1)))^((-1)+Earth_M*G*R^(-1)*Temp_grad^(-1)*(Earth_R+ ...
  h_g*x(1))^(-2))*(Cd+3*A_ref^(-1)*CS_Area*CS_CD0*u*(1+( ...
  -0.557171E-2)*CS_Length^(-1)*CS_Position*(CS_Position*x(2))^( ...
  -0.25E0)*csc(CS_Angle))*sin(CS_Angle)^3);

dFdu = (-3/2)*CS_Area*CS_CD0*Mass_dry^(-1)*rho_g*x(2)^2*(Temp_g^(-1)*( ...
  Temp_g+(-1)*Temp_grad*x(1)))^((-1)+Earth_M*G*R^(-1)*Temp_grad^( ...
  -1)*(Earth_R+h_g*x(1))^(-2))*(1+(-0.557171E-2)*CS_Length^(-1)* ...
  CS_Position*(CS_Position*x(2))^(-0.25E0)*csc(CS_Angle))*sin( ...
  CS_Angle)^3;

if imag(dFdz)+imag(dFdv)+imag(dFdu) ~= 0
    disp('warning: imaginary part in A or B matrix')
end
dFdz = ureal('dFdz',dFdz,'percent',40);
dFdv = ureal('dFdv',dFdv,'percent',70);
dFdu = ureal('dFdu',dFdu,'percent',40);

A = [0 1; dFdz dFdv];
B = [0; dFdu];
K = dx-A*x-B*u;

end