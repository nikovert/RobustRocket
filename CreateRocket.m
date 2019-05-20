function rocket = CreateRocket(string)
    val = readtable(string,'Format','%s%f');
    prop=table2array(val(1:end-1,2));
    motorname=table2array(val(end,1));
    rocket = {};
    rocket.Length = prop(1);
    rocket.Cone_L = prop(2);
    rocket.D = prop(3);
    rocket.fin_n = prop(4);
    rocket.fin_h = prop(5);
    rocket.fin_base = prop(6);
    rocket.fin_top = prop(7);
    rocket.fin_sweep = deg2rad(prop(8));
    rocket.fin_t = prop(9);
    rocket.Mass_dry = prop(10);
    rocket.Ibody_dry = [prop(11), 0 ,0; 0 , prop(12), 0; 0, 0, prop(13)];
    rocket.Xcm_dry = prop(14);
    rocket.L_pinDia = prop(15);
    rocket.L_pinH = prop(16);
    rocket.CS_Area = prop(18);
    rocket.CS_CD0 = prop(19);
    rocket.CS_Angle = deg2rad(prop(20));
    rocket.CS_Position = prop(21);
    rocket.CS_Length = prop(22);
    rocket.CS_Delay = prop(23);

    % Launch Rail Heading
    Ra = deg2rad(90-prop(17));
    Rax = [1;    0;         0];
    rocket.Q = [cos(Ra/2) sin(Ra/2)*Rax(1) sin(Ra/2)*Rax(2) sin(Ra/2)*Rax(3)]';
    rocket.A_ref = (pi*rocket.D^2/4);
end