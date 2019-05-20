function  [t, state] = accent_calc(roro,tend,dt,getControlValueHandle)
    %Function calculates the assent phase of the rocket
    
    if nargin < 3
        dt = 0.1;
    end
    
    if nargin < 4
        use_control = false;
    else
        use_control = true;
    end
    
    % Multiplies Cd to "fit" flight data...
    % global drag_correction_factor; drag_correction_factor=1.885; % for mestral
    global drag_correction_factor; drag_correction_factor=1.0; % for Tell
    % Noise realizations for MonteCarlo
    t_vec = 0:0.1:tend;
    global noise_vals; noise_vals = 0.00001*randn(13, length(t_vec));
    global noise_vals_h; noise_vals_h = sqrt(roro.noise_variance_h)*randn(1, length(t_vec));
    global noise_vals_v; noise_vals_v = sqrt(roro.noise_variance_v)*randn(1, length(t_vec));
    
    global env;  
    state_0 = [roro.X; roro.Q; roro.P; roro.L];
    tspan = [roro.time:dt:tend];
     
    % Event function to stop at max height
    % Included MaxStep to define the maximum stepsize. That measure
    % prevents ode45 from taking too large steps (like 1.3s which happened)
    options = odeset('Events',@event_function,'MaxStep',0.1);
    
    % Solve flight using ODE45
    [t, state]= ode45(@flight,tspan,state_0,options);

    % --------------------------------------------------------------------------
    %% Equations of motion discribed to be sloved in ode45 
    function state_dot = flight(t,state)
     
        if use_control
            roro.set_u(getControlValueHandle(roro,t));
        end
                
        %TODO: put condition on burn data so it does not excecute after
        %bunout
        
        if (t>0)
            roro.deltat = t - roro.time;
            roro.time = t;
            burn_data(roro); % runs each cycle to update motor stats 
        end
        
        % Add some noise (MonteCarlo)
        if (roro.b_add_noise)
            if (roro.t_Burnout > 0)
                if t>1
%                     noise_value = interp1(t_vec, noise_vals', t);       
%                     state = state + noise_value';
                noise_value_h = interp1(t_vec, noise_vals_h', t);    
%                 noise_value_v = interp1(t_vec, noise_vals_v', t);       
%                 state(3) = state(3) + (dt/6.1*1e-3) .*norm(state(8:10)./roro.Mass) .* noise_value_h';  
                state(3) = state(3) + norm(state(8:10)./roro.Mass) .* noise_value_h';   
%                 Xdot(3) = Xdot(3) + noise_value_v';
                end
            end
        end
        
        X = state(1:3);
        Q = state(4:7);
        P = state(8:10);
        L = state(11:13);
        

        roro.X= state(1:3);
        roro.Q= state(4:7);
        roro.P= state(8:10);
        roro.L= state(11:13);
        
        % Rotation matrix for transforming body coord to ground coord
        Rmatrix= quat2rotm(roro.Q');
        
        % Axis wrt earth coord
        YA = Rmatrix*env.YA0'; 
        PA = Rmatrix*env.PA0'; 
        RA = Rmatrix*env.RA0'; 
        CnXcp = roro.CnXcp;
        Cn= CnXcp(1);
        Xcp= CnXcp(2);
        Cda = CnXcp(3); % Damping coefficient
        zeta = CnXcp(4); % Damping ratio
        Ssm = CnXcp(5); % Static stability margin
        %% ------- X Velocity-------
        Xdot=P./roro.Mass;
        
        %% ------- Q Angular velocity--------- in quarternians 
        invIbody = roro.Ibody\eye(3); %inv(roro.Ibody); inverting matrix
        omega = Rmatrix*invIbody*Rmatrix'*L;
        s = Q(1);
        v =[Q(1); Q(2); Q(3)];
        sdot = 0.5*(dot(omega,v));  % DONE (thomas): fixed sign mistake
        vdot = 0.5*(s*omega + cross(omega,v));
        Qdot = [sdot; vdot];
        
        %% -------Angle of attack------- 
        % Angle between velocity vector of the CoP to the roll axis, given in the ground coord        
        % To Do : windmodel in env, Model gives errors 
        if(norm(X) < roro.Rail)
            W = [0, 0, 0]';
        else
            W = env.W;
        end
        
        Vcm = Xdot  + W;
        Xstab = Xcp- roro.Xcm;

        omega_norm = normalize(omega);
        Xperp =Xstab*sin(acos(dot(RA,omega_norm))); % Prependicular distance between omaga and RA
        
        Vomega = Xperp *cross(RA,omega);
        
        V = Vcm + Vomega; % approxamating the velocity of the cop        
        
        Vmag = norm(V);
        
%         if Vmag > 350
%             keyboard
%         end
        
        Vnorm = normalize(V);
        %If the motor performance fluctuates over the first few centimeters
        %the backward velocity of the rocket messes up the program. Thus
        %alpha is kept zero as long as the launch rail is not travelled
        %halfway
        if ((Vmag == 0) || (X(3)<3))
            alpha = 0;
        else
            alpha = acos(dot(Vnorm,RA));
        end
        roro.alpha = alpha;
        roro.alpha_angle = [roro.alpha_angle,alpha]; 
        %% ------- P Forces = rate of change of Momentums-------

        Fthrust = roro.T*RA;
        
        mg = roro.Mass*env.g;
        Fg = [0, 0, -mg]';
        
        %If we get to this point, then somewhere somebody messed up and
        %probably saved the linear momentum instead of velocity
        if(Vmag > 400)
            keyboard;
        end
        
        % Axial Forces
        if(roro.Takeoff)
            Famag = 0.5*env.rho*Vmag^2*roro.A_ref*roro.Cd(Vmag/env.C)*cos(roro.alpha); 
        else
            Famag = 0;
        end
        
        Fa = -Famag*RA;
        
        % Normal Forces
        Fnmag = 0.5*env.rho*Vmag^2*roro.A_ref*Cn + Famag*sin(roro.alpha);
        
        RA_Vplane = cross(RA,Vnorm);
        Fn = Fnmag*(cross(RA,RA_Vplane));
        
        %Simulate the rocket standing on the launch pad by setting the total
        %force equal to the normal force resulting from the influence of
        %gravity
        if (roro.T< mg && not(roro.Takeoff))
            Ftot = [0, 0, mg]';
        else
            Ftot = Fthrust*roro.Motor_efficiency + Fg + Fa + Fn;
            %Takes care of a bug by which the rocket loses thrust after
            %takeoff
            roro.Takeoff = 1;
        end
        %% ------- L Torque-------
        Trqn = Fnmag*Xstab*(RA_Vplane); 
        
        m=diag([1, 1, 0]);
        invR = Rmatrix';
        Trq_da = -Cda*Rmatrix*m*invR*omega;
        %Tqm=(Cda1*omega)*omegaax2; rotational torque by motor
%        r_f = %TODO roll damping 
%        Trmag = 0.5*env.rho*V^2*roro.A_ref*roro.Cld*r_f;
%        Tr = Trmag*RA;
        if(norm(X) < roro.Rail)
            Trq = [0, 0, 0]';
        else
            Trq = Trqn+Trq_da;
        end
        
        %% -------Update rocket state derivatives-------
        %disp(strcat('Current time: ', num2str(t)));
        
        % Add some noise (MonteCarlo)
        if (roro.b_add_noise)
            if (roro.t_Burnout > 0)
%                 noise_value_h = interp1(t_vec, noise_vals_h', t);    
                noise_value_v = interp1(t_vec, noise_vals_v', t);       
%                 state(3) = state(3) + noise_value_h';   
%                 Xdot(3) = Xdot(3) + (dt/6.1*1e-3) .* norm(Xdot) .* noise_value_v';
                Xdot(3) = Xdot(3) + norm(Xdot) .* noise_value_v';
            end
        end
        
        roro.Xdot= Xdot;
        roro.Qdot= Qdot;
        roro.Pdot= Ftot;
        roro.Ldot= Trq;
            
        state_dot =[Xdot; Qdot; Ftot;Trq];
        roro.state_dot = [roro.state_dot,state_dot]; % save result

        roro.state =     [roro.state_dot,[X;Q;P;L]]; % save result
       
        %% -------Burnout time-------
        if(roro.propM_current<0.0000001 && roro.t_Burnout == 0 )
            roro.t_Burnout = t;
        end

        %% -------Takeoff time-------
        if(roro.t_Takeoff==0&&roro.Takeoff==1)
            roro.t_Takeoff = t;
        end
        
        %% --------Braking time------
        if(t > roro.t_Burnout + roro.CS_Delay && roro.t_Burnout ~= 0)
            %Only start the Brakes after a time dealy
            roro.CS_u = 1;
        end
         
        
        %% --------Generate perfect (no noise) Sensors data------
        % IMU (No noise!)
        while  (length(roro.imu_t)>=1) && (t<roro.imu_t(end))
            roro.imu_t = roro.imu_t(1:end-1);
            roro.imu_a = roro.imu_a(1:end-1,:);
        end
%         if (length(roro.imu_t)<1) || (t>roro.imu_t(end) + roro.sensors_sample_time)
            roro.imu_t = [roro.imu_t; t];
            roro.imu_a = [roro.imu_a; ...
                        Ftot(1)/roro.Mass,Ftot(2)/roro.Mass,Ftot(3)/roro.Mass];
            roro.imu_g = [roro.imu_g; roro.Qdot(1), roro.Qdot(2), roro.Qdot(3)];
%         end
        % Pressure sensor (No noise!)
        
        % Debug
%         disp(Ftot)
%         disp(roro.Mass)
        %fprintf('Current accelerations:\n  ax %f\n  ay %f\n  az %f\n', ...
         %       Ftot(1)/roro.Mass, Ftot(2)/roro.Mass, Ftot(3)/roro.Mass);
        %fprintf('Current angular rates:\n  gx %f\n  gy %f\n  gz %f\n\n', ...
         %       roro.Qdot(1), roro.Qdot(2), roro.Qdot(3));
        
        %% Log Data
       
        %logData(roro.alpha, roro.Cd, Cda, roro.Xcm, roro.Mass, Vmag, Xcp, zeta, Ssm, t);
        
    end
    
    function [value,isterminal,direction] = event_function(t,state)
    %% stops ode integration when the max height is reached 
        if (t > 1 && state(10) <= 0) % Linear momentum in z direction is zero
            value = 0; % when value = 0, an event is triggered
        else
            value =1;
        end
        isterminal = 1; % terminate after the first event
        direction = 0;  % get all the zeros
        
        % Stops at burnout
        if (roro.b_stop_simulation_at_burnout)
            if (t > 1 && roro.t_Burnout>0)
                disp('Burnout! Stopping simulation here.');
                value = 0;
            end
        end
    end
end
