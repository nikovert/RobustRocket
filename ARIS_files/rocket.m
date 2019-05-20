%% Rocket
% This is the class that contains all relevant information for the rocket.
% It includes simulation parameters as well as controller parameters.

% Defines all the properties of the rocket. The dynamics are continuously
% updated, any varibales need for functions should come from this class
classdef rocket < handle & matlab.mixin.Copyable
    
    % Internal airbrake Characteristics that change during the flight 
    properties (Access = private, Hidden = true)
        CS_u_real = 0           % The current Control input (i.e. position of the airbrakes)
        CS_time_prev = 0        % The time of the last control update
    end
    
    % Airbrake and controller characteristics that change during the flight  and may be viewed during the flight
    properties (SetAccess = private, GetAccess = public)
        CS_u = 0;                % The target Control input (i.e. position of the airbrakes)
        
        ControllersampleTime     % The sample time of the MPC
        Ki                       % The weight of the integrated error in height and velocity [Ki_height, Ki_veloctiy]
        wx                       % The weight of the state-reference error
        wu                       % The weight of the input change
    end
    
    % Controller characteristics and logs that are set during the flight 
    properties (Access = public)
        u_list_Controller = [0;0];     % The control input of every MPC run
        u_list_real = [0;0];           % A list of all the inputs over time  
        w_list = [0; 0; 0];            % A list of the integration error over time [w_height; w_velocity; time]
        Controller_called = 0;         % The number of times we called the MPC;
        
        error_height = [];             % The total height error to the reference trajectory at everytime
        error_velocity = [];           % The total velocity error to the reference trajectory at everytime
        
        u_next = 0;                    % The control input to be set at the next time step (control time, not simulation timestep)
        next_t = 5;                    % The next time step (control time, not simulation timestep)
        w = [0; 0];                    % The integrated error ([height, velocity])
    end
    
    % Simulation properties that change during the flight
    properties (Access = public)
        B_BRAKING = false       % Defines whether we simulate braking or not
        B_BRAKING_TABLE = false % [WARNING, discontinued] Use input from table generated offline
        B_BRAKING_MPC = false;  % Defines if we use the MPC controller (also used for the LQR)
        Motor_efficiency        % This state can be used to simulate stronger or weaker burn, set 1 for nominal behaviour
        Takeoff = 0             % This state determines if the rocket has overcome gravity
        Braking                 % State for Forces.m

        % Motor Characteristics: Updated in each iteration of the accent_calc     
        motorname             % The name of the motor model
        motordata             % The motor data file
        Mass_motor            % The inital mass of the motor
        Motor_impulse         % Total Impulse of motor
        propM_tot             % Total mass of the propelent
        Iprop                 % Inertia matrix of prop wrt Cg
        Xcm_prop              % Center of mass of the prop
        prop_OD
        prop_ID
        prop_h               
        prop_density          % Propelant density

        deltaMass = 0;            % The Mass change at each timestep
        impulseGen                % Impulse generated up to the current time 
        propM_current             % Remaning prop mass
        propM_prev = 0;           % Mass for previous time step to calcualte deltaMass

        departureState            % The departure state, logged during the simulation
        t_Burnout = 0;            % The time at which the Motor burntout, logged during the simulation
        t_Takeoff = 0;            % The takeoff time, logged during the simulation

        % Current State Vector with Initial values, Updated in accent_calc
        time = 0;                 % The Simulation time
        X = [0; 0; 0];            % Position x, y, z   
        Q = [1; 0; 0; 0];         % Angle in quarternions  
        P = [0; 0; 0];            % Linear Momentum  
        L = [0; 0; 0];            % Angular momentum
        Xdot = [0; 0; 0];         % Velocity xdot, ydot, zdot 
        Qdot = [0; 0; 0; 0];      % Angular rates
        Pdot = [0; 0; 0];         % Linear Momentum rates = applied force
        Ldot = [0; 0; 0];         % Angular Momentum rates= applied torques
        alpha = 0;                % Angle of the rocket to the z-Axis
        alpha_angle = []          % A log for the angle over time

        % Save derivatives (speeds are nice)
        state_dot = [];          % The change of state, carefull, not updated in real time
        state     = [];          % The current state, carefull, not updated in real time

        % Sensors: N IMU + M pressure sensor readings saved in an array
        imu_t = [];              % [WARNING, discontinued] Time of IMO sample: synch. with gyro/accel Nx[s]
        imu_a = [];              % [WARNING, discontinued] Accelerometer: Accels    ax, ay, az  Nx3*[m/s^2]
        imu_g = [];              % [WARNING, discontinued] Gyroscope, Angular rates gx, gy, gz  Nx3*[rad/s]
        temp  = [];              % [WARNING, discontinued] Temperature
        press = [];              % [WARNING, discontinued] Pressure sensor: Mx3*[time[s], pressure[Pa/m^2], altitude[m]]
        b_sensors_noise = false; % [WARNING, discontinued] true if sensors contain noise
        sensors_sample_time = 0.01; %[WARNING, discontinued]
        
        % if true, stops the simulation after burnout
        %caution, seems to not work well
        b_stop_simulation_at_burnout = false; % If set, accentCalc function will only simulated up to motor Burnout
    end
     
    %Simulation properties that don't change during the flight
    properties (SetAccess = private, GetAccess = public)
        Length          % Length of the Rocket
        D               % Diameter of the Rocket
        Cone_L          % Length of the nosecone
        fin_n           % Number of fins on the rocket
        fin_h           % Height of the fins
        fin_base        % Base width of the fin
        fin_top         % Top width of the fin
        fin_sweep       % Sweep angle of the fin
        fin_t           % 
        Mass_dry        % Rocket and motor housing no prop
        Ibody_dry       % Rocket and motor housing no prop
        Xcm_dry         % Rocket and motor housing no prop
        Rail = 5.2      % According to newest IREC data
        L_pinDia        % Railpin diameter
        L_pinH          % Railpin height
        A_ref           % Area of attack of the rocket
        CS_Area         % Describes the full usable Area in [m^2] of one control surface
        CS_CD0          % The drag coefficient of the used control surface
        CS_Angle = pi/2 % The angle in [rad] at which the control surface is deployed
                        % measured from the vertical -> 90deg = full extension
        CS_Position = 1.54 %Position of the air brakes in [m]
        CS_Length          %Length of the Control Surface in [m]
        CS_slew = 0        %The time it takes for a control input to be executed
    end
    
    %Properties that shouldn't be used and that will be discontinued
    properties (Access = public)
        B_BRAKING_GUI = false   % [WARNING, GUI discontinued] Use GUI for control 
        CS_Delay                % [WARNING, discontinued] The time it takes to make control plan, i.e. before we start to control (not used)
        
        % Adds gaussian noise if true [WARNING, discontinued]
        noise_variance_h = 0.1; % [WARNING, discontinued]
        noise_variance_v = 0.1; % [WARNING, discontinued]
        b_add_noise = false;    % [WARNING, discontinued]
        
        deltat                    % [WARNING, discontinued] Size of time step calcualted in accent_calc
        
        %%%% CONTROL CPP PARAMETERS %%%%
        
        h_values_upper = [] % [WARNING, Control table discontinued]
        h_values_lower = [] % [WARNING, Control table discontinued]
        v_step              % [WARNING, Control table discontinued]
        speed2idxv_inUpperLowerLim % [WARNING, Control table discontinued]
        img_control = []    % [WARNING, Control table discontinued]
        img_size = [-1,-1]  % [WARNING, Control table discontinued]

        % MonteCarlo params
        b_entered_crit_region = []; % [WARNING, discontinued]
    end
    
    methods (Access = public)
        % Constructor takes data from rocket property file
        function obj = rocket(val) % Val is the table of properties val2 motor data
            if nargin > 0
                prop=table2array(val(1:end-1,2));
                motorname=table2array(val(end,1));
                %if isnumeric(val)
                obj.Length = prop(1);
                obj.Cone_L = prop(2);
                obj.D = prop(3);
                obj.fin_n = prop(4);
                obj.fin_h = prop(5);
                obj.fin_base = prop(6);
                obj.fin_top = prop(7);
                obj.fin_sweep = deg2rad(prop(8));
                obj.fin_t = prop(9);
                obj.Mass_dry = prop(10);
                obj.Ibody_dry = [prop(11), 0 ,0; 0 , prop(12), 0; 0, 0, prop(13)];
                obj.Xcm_dry = prop(14);
                obj.L_pinDia = prop(15);
                obj.L_pinH = prop(16);
                obj.CS_Area = prop(18);
                obj.CS_CD0 = prop(19);
                obj.CS_Angle = deg2rad(prop(20));
                obj.CS_Position = prop(21);
                obj.CS_Length = prop(22);
                obj.CS_Delay = prop(23);
                obj.Motor_efficiency = prop(24);
                
                % Launch Rail Heading
                Ra = deg2rad(90-prop(17));
                %Rax = [0.9511;    0.3090;         0];
                Rax = [1;    0;         0];
                obj.Q = [cos(Ra/2) sin(Ra/2)*Rax(1) sin(Ra/2)*Rax(2) sin(Ra/2)*Rax(3)]';

                %motor
                obj.motorname=motorname;
                obj.A_ref = (pi*obj.D^2/4);
                %else
                %error('Value must be numeric')
                %end
            end
        end
        
        
        function setupMPC(obj,Ki,wx,wu,MPCsampleTime)
            % Set the weights of the MPC
            obj.Ki = Ki;
            obj.wx = wx;
            obj.wu = wu;
            obj.ControllersampleTime = MPCsampleTime;
        end
        
        
        function hard_set_u(obj,u)
            % Sets u, only to be used during initialisation
            obj.CS_u = u;
            obj.CS_u_real = u;
        end
        
        
        function set_slew(obj, slew)
            % Set the slew rate
            obj.CS_slew = slew;
        end
        
        
        function set_u(obj,u)
            % Used to build in slew rate constraints
            obj.CS_u = u;
        end
        
        
        function u = get_u(obj)
            % Used to build in slew rate constraints, called from Cd
            if obj.CS_slew ~= 0
                delt_t = obj.time - obj.CS_time_prev;
                obj.CS_time_prev = obj.time;
                
                max_aloud_change = delt_t/obj.CS_slew;
                
                if max_aloud_change >= abs(obj.CS_u-obj.CS_u_real)
                    obj.CS_u_real = obj.CS_u;
                else
                    obj.CS_u_real = obj.CS_u_real + max_aloud_change * sign(obj.CS_u-obj.CS_u_real);
                end
                if obj.CS_u_real > 1 || obj.CS_u_real < 0
                    disp('ERROR: CS_u out of bounds');
                    obj.CS_u_real = max(min(1,obj.CS_u_real),0);
                end
            else
                obj.CS_u_real = obj.CS_u;
            end
            u = obj.CS_u_real;
        end
        
        
        function Cd = Cd(obj,v)
            % Fuctions to calcualte the drag in axial direction, and update
            % flight parameters
            
            % make sure the airbrakes aren't extended before flight
            if obj.t_Burnout > 0
                current_u = obj.get_u;
            else
                current_u = 0;
            end
            
            if (current_u > 0.000001)
                Re = norm(obj.Xdot)*obj.CS_Position/13.164e-6; %kinematic viscosity of air
                h = obj.CS_Length*sin(obj.CS_Angle); %height of the air brake
                delta = 0.37*obj.CS_Position/(Re^0.25); %height of the boundary layer
                CD0_cfd = [1.179337962 1.2395853 1.2360735 1.229218227 1.3254];
                M  = [0.3 0.5 0.6 0.8 0.9];
                CS_CD0_corrected = interp1(M,CD0_cfd,v,'linear','extrap');
                CD_theta_corrected = CS_CD0_corrected*(1-0.25*delta/h)*sin(obj.CS_Angle)^3;
                Cd = Cd_mandell(obj) +  current_u * obj.CS_Area/obj.A_ref*CD_theta_corrected*3;
                if (not(isreal(Cd)) && obj.Takeoff)
                    disp('whatt?!');
                    keyboard
                end
            else
                Cd = Cd_mandell(obj);
                if (not(isreal(Cd)) && obj.Takeoff)
                    disp('whatt?!');
                    keyboard
                end
            end
            if (isinf(Cd) || Cd > 10)
                Cd =10;
            end
        end
       
        
        function CnXcp = CnXcp(obj)
            % Normal force and Cop location
            [Cn_alpha, Xcp, Cda, zeta, Ssm]=Cn_alphaXcp(obj);
            CnXcp = [Cn_alpha*obj.alpha, Xcp, Cda, zeta, Ssm];
            % Bad
            [Cn_alpha, Xcp, Cda, zeta, Ssm, Ssm_B, Ccm]=Cn_alphaXcp(obj);
            CnXcp = [Cn_alpha*obj.alpha, Xcp, Cda, zeta, Ssm, Ssm_B, Ccm];
        end
        
        
        function T = T(obj)
            % Set the thrust curve
            M = obj.motordata;
            T = interp1(M(:,1),M(:,2),obj.time);
            if ( obj.time > M(end,1))
               T = 0;
            end
        end

         
        function T_attime = T_attime(obj,t)
            % Thrust curve mapped to time
            M = obj.motordata;
            if (t < M(2,1))
               T_attime = M(2,2);
            else
               T_attime = interp1(M(:,1),M(:,2),t);
            end
            if ( t > M(end,1))
               T_attime = 0;
            end
        end

        function Re = Re(obj)
            % Re of rocket
            global env
            Re = env.rho * norm(obj.Xdot) * obj.Length/env.mu;
        end

        function Mass = Mass(obj)
            % Update the mass of rocket
            M = obj.motordata;
            if ( obj.time > M(end,1)) % To assure it goes to zero incase of integartion error
               obj.propM_current = 0;
            end
            Mass= obj.Mass_dry + obj.propM_current;               
        end

        
        function Xcm = Xcm(obj) 
            % Current Center of mass of rocket
            M = obj.motordata;
            if ( obj.time > M(end,1)) % To assure it goes to zero incase of integartion error
               Xcm = obj.Xcm_dry;
               return;
            end
            Xcm= obj.Xcm_dry*obj.Mass_dry + obj.Xcm_prop*(obj.propM_current); 
            Xcm = Xcm/obj.Mass;
        end
        
        
        function PrintMass(obj,description) 
            % Print current mass and Center of mass of the rocket
            disp(description);
            disp(['Rocket Mass: ',num2str(obj.Mass),'              Rocket Center of Mass: ',num2str(obj.Xcm)]);
            disp(['Rocket Mass dry: ',num2str(obj.Mass_dry),'            Rocket Center of Mass dry:',num2str(obj.Xcm_dry)]);
        end

        
        function Trim(obj,trim_mass,trim_Xcm) 
            % Change the mass and the center of mass according to an added trim mass
            old_mass = obj.Mass_dry;
            obj.Mass_dry = old_mass + trim_mass;

            old_Xcm = obj.Xcm_dry;
            obj.Xcm_dry = (old_Xcm*old_mass + trim_Xcm*trim_mass)/(old_mass + trim_mass);
        end
        
        function Ibody = Ibody(obj) 
           % Current Inertia of the rocket
           M = obj.motordata;
           Ibody = obj.Ibody_dry;
           if ( obj.time <= M(end,1))
                Ibody = Ibody + obj.Iprop; 
           end
        end
        
        function turnOnBrakes(obj)
           % imidiatly switches on the brakes
           obj.CS_u = 1;
        end
        
        function setXdot(obj,V)
            % Manually set the velocity in z
            obj.Xdot = [0;0;V];
        end

        function setalpha(obj,a)
            % Manually set angle to the x-axis
            obj.alpha = a;
        end

        function setLaunchRailAngle(obj, angle)
            % Launch Rail Heading
            Ra = deg2rad(90-angle);               
            Rax = [1;    0;         0];
            obj.Q = [cos(Ra/2) sin(Ra/2)*Rax(1) sin(Ra/2)*Rax(2) sin(Ra/2)*Rax(3)]';
        end
    end

end
