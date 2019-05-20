function motor_init(roro)
    %% data from thrustcurve.org
    % Reads the .eng thrust file of the motors and calculates the reletive 
    % the reletive properties of the motor. 
    
    % Extracting thrust curve from .eng file
    R1=4;
    C1=0;
    motorname = char(roro.motorname);
    motordata = dlmread(motorname,'',R1,C1);
    roro.motordata = [0, 0; motordata]; % at 0 0 at start of date
    
    % Very strange way of reading mass of motor amd propellent from the
    % Opening File (This has to be done again to get header data)
    fid = fopen(motorname);
    % Skipping to the correct line 
    tline = fgets(fid);
    tline = fgets(fid);
    tline = fgets(fid);
    tline = fgets(fid);
    % Reading mass and motor size from the .eng file
    C = strsplit(tline);
    Motor_diameter = str2double((C(2)))*1e-3; %[m]
    Motor_lenght = str2double((C(3)))*1e-3; %[m]
    Mass_prop = str2double((C(5)));
    Mass_motor = str2double((C(6)));
    roro.propM_tot = Mass_prop;
    roro.Mass_motor = Mass_motor - Mass_prop;
    fclose(fid);
    
    %Initally current and total prop mass are equal
    roro.propM_current = Mass_prop;
    % Calculated impulse from motor data using spline
    roro.Motor_impulse = trapz(roro.motordata(:,1),roro.motordata(:,2));
    
    if(strcmp(motorname,'AeroTech_M2400.eng'))
        %M2400 by Gianni
        prop_density = 1.579e+03;
        prop_OD = 0.0855;
        prop_ID = 0.0285;
        prop_h = 0.4572;
        roro.Xcm_prop = roro.Length - prop_h/2;
    elseif(strcmp(motorname,'AeroTech_M2500.eng'))
        %M2500 by Gianni
        prop_density = 1.51e+03;
        prop_OD = 0.0855;
        prop_ID = 0.0285;
        prop_h = 0.6096;
        roro.Xcm_prop = roro.Length - prop_h/2;
    elseif(strcmp(motorname,'AeroTech_M1419.eng'))
        %M1419 by Bogdan
        prop_density = 1717.980263;
        prop_OD = 0.0855;
        prop_ID = 0.0285;
        prop_h = 0.465;
        roro.Xcm_prop = roro.Length - prop_h/2;
    elseif(strcmp(motorname,'Cesaroni_M1060.eng'))
        %M1060 by Bogdan
        prop_density = 1635.268715;
        prop_OD = 0.0855;
        prop_ID = 0.0285;
        prop_h = 0.434;
        roro.Xcm_prop = roro.Length - prop_h/2;
    elseif(strcmp(motorname,'Cesaroni_M2505.eng'))
        %M2505 by Bogdan
        prop_density = 1406.999288;
        prop_OD = 0.0855;
        prop_ID = 0.0285;
        prop_h = 0.465;
        roro.Xcm_prop = roro.Length - prop_h/2;
    elseif(strcmp(motorname,'AeroTech_M2100.eng'))
        %M2100 by Bogdan
        prop_density = 1598.314333;
        prop_OD = 0.0855;
        prop_ID = 0.0285;
        prop_h = 0.484;
        roro.Xcm_prop = roro.Length - prop_h/2;
    elseif(strcmp(motorname,'AeroTech_M4500.eng'))
        %M4500 by Bogdan
        prop_density = 1389.452987;
        prop_OD = 0.0855;
        prop_ID = 0.0285;
        prop_h = 0.483;
        roro.Xcm_prop = roro.Length - prop_h/2;
    elseif(strcmp(motorname,'Cesaroni_M1790.eng'))
        %M1790 by Bogdan
        prop_density = 1568.212851;
        prop_OD = 0.0855;
        prop_ID = 0.0285;
        prop_h = 0.588;
        roro.Xcm_prop = roro.Length - prop_h/2;
    elseif(strcmp(motorname,'AeroTech_L2200.eng'))
        %L2200G by Bogdan
        prop_density = 1.6188e+03;
        prop_OD = 0.0645;
        prop_ID = 0.0280;
        prop_h = 0.567;
        roro.Xcm_prop = roro.Length - prop_h/2;
        
        
    %% Small engines 
    elseif(strcmp(motorname,'AeroTech_L1390G.eng'))
        %L1390G by Gianni
        prop_density = 1.7062e+3;
        prop_OD = 0.06495;
        prop_ID = 0.02223;
        prop_h = 0.395;
        roro.Xcm_prop = roro.Length - prop_h/2;
    elseif(strcmp(motorname,'AeroTech_I245G.eng'))
        %I245G by Gianni
        prop_density = 1.707e+03;
        prop_OD = 0.0332;
        prop_ID = 0.01113;
        prop_h = 0.0461;
        roro.Xcm_prop = roro.Length - prop_h/2;
    elseif(strcmp(motorname,'AeroTech_I300T.eng'))
        %I300T by Gianni
        prop_density = 1.566e+03;
        prop_OD = 0.0332;
        prop_ID = 0.01113;
        prop_h = 0.0461;
        roro.Xcm_prop = roro.Length - prop_h/2;
    elseif(strcmp(motorname,'AeroTech_K540M.eng'))
        %K540M by Gianni
        prop_density = 3.017e+03;
        prop_OD = 0.0475;
        prop_ID = 0.0158;
        prop_h = 0.0461;
        roro.Xcm_prop = roro.Length - prop_h/2;
    elseif(strcmp(motorname,'AeroTech_L1150.eng'))
        %L1150R by Bogdan
        prop_density = 1.6096e+3;
        prop_OD = 0.0645;
        prop_ID = 0.0235;
        prop_h = 0.395;
        roro.Xcm_prop = roro.Length - prop_h/2;
    elseif(strcmp(motorname,'AeroTech_L1040.eng'))
        %L1040DM by Bogdan
        prop_density = 1.664e+03;
        prop_OD = 0.0645;
        prop_ID = 0.0275;
        prop_h = 0.567;
        roro.Xcm_prop = roro.Length - prop_h/2;
    elseif(strcmp(motorname,'AeroTech_L850.eng'))
        %L1150R by Bogdan
        prop_density = 1.75061e+3;
        prop_OD = 0.0645;
        prop_ID = 0.0235;
        prop_h = 0.395;
        roro.Xcm_prop = roro.Length - prop_h/2;
    else
        printf('Motorname not recognized, please update motor_init.m')
    end
    
    roro.prop_density= prop_density;
    roro.prop_OD = prop_OD;
    roro.prop_h = prop_h;
    
    %Initializing propellent inertias w.r.t. cm
    d = roro.Xcm_prop - roro.Xcm; % Note: The mass of the prop should already be updated in roro to get correct Xcm 
    propIx = 0.5*Mass_prop*(prop_OD^2+prop_ID^2)/4;
    propIy = Mass_prop/12*(3*(prop_ID^2+prop_ID^2)/4 + prop_h^2) + Mass_prop*(d);
    roro.Iprop = [propIx, 0, 0; 0, propIy, 0; 0, 0, propIy]; % propIy = propIz 
    
    
    
%     TODO: implement same thing using splines
%     spl = spline(roro.motordata(:,1),roro.motordata(:,2));
% 
%     % integrate the spline
%     spl.coefs = spl.coefs*[diag(1./[4 3 2 1]'),zeros(4,1)];
%     spl.order = 5;
%     dx = diff(spl.breaks)';
%     C = spl.coefs(:,1);
%     for i = 1:4
%       C = C.*dx + spl.coefs(:,i+1);
%     endzeros(3,3); 
%     spl.coefs(2:end,5) = cumsum(C(1:(end-1));
end
