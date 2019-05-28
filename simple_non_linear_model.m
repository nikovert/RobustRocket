function dx = simple_non_linear_model(rocket, x,u)
%returns a simplified non linear model of the system based on a
%predetermined atmosphere model.

%From the rocket
A_ref = rocket.A_ref;
Mass_dry = rocket.Mass_dry;
CS_Length = rocket.CS_Length;
CS_Angle = rocket.CS_Angle;
CS_Position = rocket.CS_Position;
CS_Area = rocket.CS_Area;
CS_CD0 = rocket.CS_CD0;
Xdot = x(2);

Re = norm(x(2))*rocket.CS_Position/13.164e-6; %kinematic viscosity of air
h = rocket.CS_Length*sin(rocket.CS_Angle); %height of the air brake
delta = 0.37*rocket.CS_Position/(Re^0.25); %height of the boundary layer
CD_theta_corrected = rocket.CS_CD0*(1-0.25*delta/h)*sin(rocket.CS_Angle)^3;

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


%We use environement(1400, 40, 86000, rocket) as constant
g = 9.7920;     %e.g
rho = 1.1407;    %e.rho

Fa = -g - 0.5/Mass_dry*rho*A_ref*Cd*x(2).^2;
Fa = ureal('Fa',Fa,'percent',40);
dx = [x(2); Fa];

end