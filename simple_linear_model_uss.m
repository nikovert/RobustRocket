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

v = Xdot/352.2912; % speed in Mach
Re = norm(Xdot)*CS_Position/13.164e-6; %kinematic viscosity of air
h = CS_Length*sin(CS_Angle); %height of the air brake
delta = 0.37*CS_Position/(Re^0.25); %height of the boundary layer
CD0_cfd = [1.179337962 1.2395853 1.2360735 1.229218227 1.3254];
M  = [0.3 0.5 0.6 0.8 0.9];
CS_CD0_corrected = interp1(M,CD0_cfd,v,'linear','extrap');
CD_theta_corrected = CS_CD0_corrected*(1-0.25*delta/h)*sin(CS_Angle)^3;
Cd = Cd_mandell(Xdot) +  u * CS_Area/A_ref*CD_theta_corrected*3;
if (not(isreal(Cd)))
    disp('whatt?!');
    keyboard
end

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
dFdz = ureal('dFdz',dFdz,'percent',80);
dFdv = ureal('dFdv',dFdv,'percent',90);
dFdu = ureal('dFdu',dFdu,'percent',80);

A = [0 1; dFdz dFdv];
B = [0; dFdu];
K = dx-A*x-B*u;

end