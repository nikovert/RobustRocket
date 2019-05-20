function [Cd]=Cd_mandell_simple(Xdot) 
    %This is a highly simplified model of the original Cd_mandell function
    %It uses hard coded values, instead of reading the values from the
    %rocket model
    
    %% Set constants
    % Rocket dimentions 
    L = 2.7564;         %roro.Length; 
    L_cone = 0.52;       %roro.Cone_L;
    L_cyl = L - L_cone;
    D_cyl = 0.15;       % roro.D;
    R_cyl = D_cyl/2;

    %% Define fin
    fin.a_ref = 0;
    fin.area = 0;
    fin.a_wet = 0;
    fin.c = 0;
    fin.X_b = 0;    % fin location
    fin.sweepc = 0;
    %%
    %Fin Geometry
    fin.n = 3;              %roro.fin_n;
    fin.sweep = 0.7854;     %roro.fin_sweep;
    fin.h = 0.24;           %roro.fin_h;
    fin.topchord = 0.0800;  % roro.fin_top;
    fin.basechord = 0.24;    %roro.fin_base;
    fin.t = 0.0050;         %roro.fin_t;
    fin.a_ref = fin.n*fin.h*fin.t;
    fin.area = (fin.topchord+fin.basechord)/2*fin.h;
    fin.a_wet = fin.n*2*fin.area;
    fin.c = (fin.topchord+fin.basechord)/2;
    fin.X_b =  L - fin.basechord; % fin location

    %A_ref
    A_ref = pi*R_cyl^2;

    F_ratio = L / D_cyl;
    % Equations from Mandell
    rho =  1.1409; % density 
    mu =  1.9094e-05;  % dynamic viscosity 
    C = 351.7589; %speed of sound dry air 15C sea level
    V = norm(Xdot,2);   %ms-1 Mag of  characteristic velocity at center of pressure location 
    M = V/C;
    if M > 1
        disp('ERROR: M should not be more the one');
    end
    Re = rho * V * L/mu;

    % Approx laminar flow over rocket

    %%
    R_ogive = (L_cone^2+ D_cyl^2/4)/D_cyl;
    %A_wet

    fun = @(x) 2*pi*(sqrt(R_ogive.^2 - power((L_cone - x),2))+R_cyl-R_ogive);
    
    x = linspace(0,L_cone,200000);
    A_wet = cumtrapz(x,fun(x)); %Use Cumulative trapezoidal numerical integration.

    A_wet = A_wet + 2*pi*R_cyl*L_cyl; 



    %% Assigning to equations as discribes in mandell
    l_n = L_cone;

    d_n = D_cyl;
    d_b = D_cyl;
    d_f = D_cyl;
    d_d = D_cyl;
    l_c  = 0; % no btail
    l_b = L_cyl;
    Re_c = 5e5;
    Re_c_fins = 5e6;
    T_f = fin.t;
    l_TR = L;
    n=fin.n;
    % mid chord sweep

    temp.x1 = fin.h*tan(fin.sweep);
    temp.x2 = 0;
    temp.x2 = temp.x1 + fin.topchord - fin.basechord;

    fin.sweepc = atan2((fin.basechord/2 + (temp.x2-fin.topchord/2)),fin.h);

    %clear temp fun fun2 % slows down computation by a lot
    l_m = fin.h/acos(fin.sweepc); % length midchord

    A_fe= (fin.topchord+fin.basechord)/2*fin.h;

    A_fp = A_fe + 0.5*d_f*fin.basechord;

    %% ------Viscous Friction------
    % Viscous friction ROCKET FORBODY Cf
    B = Re_c*(0.074/Re^(0.2) - 1.328/sqrt(Re));

    if (Re < Re_c)
        Cf =  1.328/sqrt(Re);
    else
        Cf=0.074/Re^(0.2)-B/Re;
    end

    %  Viscous friction ROCKET FINS Cf_f
    Re_f  = rho * V * l_m/mu;  %Note the V is at the cop not the finneed to recalculate for better results 

    B_f = Re_c_fins*(0.074/Re_f^(0.2) - 1.328/sqrt(Re_f));

    if (Re_f < Re_c_fins)
        Cf_f =  1.328/sqrt(Re_f);
    else
        Cf_f=0.074/Re_f^(0.2)-B_f/Re_f;
    end
    
    %% -------Drag at zero AoA-------
    % Body drag, Box Eq41
    Cd_fb = (1 + 60/(l_TR/d_b)^3+0.0025*l_TR/d_b)*(2.7*l_n/d_b +4*l_b/d_b + 2*(1-d_d/d_b)*l_c/d_b)*Cf;
    % Base drag, Box Eq42
    Cd_b = 0.029*(d_d/d_b)^3/sqrt(Cd_fb);
    % Fin drag, Box Eq44
    Cd_f = 2*Cf_f *(1+2*T_f/l_m)*4*n*A_fp/(pi*d_f^2);
    % Interference drag, Box Eq44
    Cd_i = 2*Cf_f*(1+2*T_f/l_m)*4*n*(A_fp-A_fe)/(pi*d_f^2);
    %The USAF Stability and Control Datcom does not say that the fin-body 
    %interference coefficient (p. 431-434) should be used in the calculations of the 0� fin drag.
    %Consequence unclear. TODO: solve this.
    
    % Total drag coefficient at zero angle of attack
    Cd0 = Cd_fb + Cd_b + Cd_f + Cd_i;
    
    % Drag of camera shells -> base drag and front drag of faring(0.07)
    A_cam = 3*1.66e-4;
    Cd_cam = (0.12 + 0.13*(M)^2 + 0.07)*A_cam/A_ref;
    
    % Launch pin drag % estimated from Mandell 
    A_pin = 1.7080e-04; %roro.L_pinDia * roro.L_pinH;
    Cd_pin = 2*0.8*A_pin/A_ref; 
    Cd0 = Cd0 + Cd_pin;
    % compressibility correction
    %% -------Additional drag at AoA-------
    % Alpha
    alpha = 0; %roro.alpha;
    % Coefficients delta dn eta from windtunnel experiments, See Box p13
    deltaktab=[4 6 8 10 12 14 16 18 20;0.78 0.86 0.92 0.94 0.96 0.97 0.975 0.98 0.982];
    etatab=[4 6 8 10 12 14 16 18 20 22 24;0.6 0.63 0.66 0.68 0.71 0.725 0.74 0.75 0.758 0.77 0.775];
    % error in paper
    etak=interp1(etatab(1,:),etatab(2,:),F_ratio,'linear','extrap');
    deltak=interp1(deltaktab(1,:),deltaktab(2,:),F_ratio,'linear','extrap');
    if etak>1;
        etak=1;
    end
    if deltak>1;
        deltak=1;
    end
    % Body drag at angle alpha
    Cd_b_alpha = 2*deltak*alpha^2 + 3.6*etak*(1.36*l_TR - 0.55*l_n)*alpha^3/(pi*d_b);
    % Fin body interference coefficients
    Rs = R_cyl/(R_cyl+fin.h);
    Kbf=0.8065*Rs^2+1.1553*Rs;
    Kfb=0.1935*Rs^2+0.8174*Rs+1;
    % Fin drag at angle alpha
    Cd_f_alpha = (1.2*A_fp*4/(pi*d_f^2) +3.12*(Kfb +Kbf-1)*A_fe*4/(pi*d_f^2))*alpha^2;
    %The method proposed by Mandel (p. 415 ) proposes to use the exposed 
    %surface SE. For a basic 4 fin rocket at the right roll angle, 
    %this would simply be two times the fin�s planform area.
    %TODO: A better model is needed here for rockets with multiple fins.

    %% -------Total Drag Coefficient-------
    Cd = Cd0 + Cd_b_alpha + Cd_f_alpha + Cd_cam;
    Cd = Cd/sqrt(1-M^2);
    CnXcp = [0 1.9655 33.2709 0.0949 3.8700 3.1976 4.5443e+03]; %roro.CnXcp;
    Cn= CnXcp(1);
    Cd = (Cd*cos(alpha) -0.5*Cn*sin(2*alpha))/(1-sin(alpha)^2);
 
end
