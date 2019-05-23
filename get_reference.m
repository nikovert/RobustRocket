function xr = get_reference(x0)
    trajectory_struct = coder.load('optimal_trajectory_p');

    % Define reference polynomial
    % h(v)
    sp = trajectory_struct.sp;
    % v(t)
    vp = trajectory_struct.vp;
    % t(v)
    tp = trajectory_struct.tp;
    
    tp_S = trajectory_struct.tp_S;
    vp_S = trajectory_struct.vp_S;
    sp_S = trajectory_struct.sp_S;
    
    tp_mu = trajectory_struct.tp_mu;
    vp_mu = trajectory_struct.vp_mu;
    sp_mu = trajectory_struct.sp_mu;
    
    velocity_error_scaling = 50;
    
    h = @(v) polyval(sp, v, sp_S, sp_mu);
    dh = @(v) polyval(polyder(sp), v, sp_S, sp_mu);
    ddh = @(v) polyval(polyder(polyder(sp)), v, sp_S, sp_mu);
    
    f = @(v) dh(v) .* (x0(1) - h(v)) + velocity_error_scaling^2 * (x0(2) - v);
    df = @(v) ddh(v) .* (x0(1) - h(v)) - (dh(v)).^2 - velocity_error_scaling^2;
    ref_v = newton(f, df, x0(2), 1e-4, 1e-8, 100);
    ref_h = h(ref_v);
    xr = [ref_h; ref_v];
end