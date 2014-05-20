function d_dt_theta = shaftSpeed(t,theta,rpm)
%Returns d_dt_theta(t) = shaft angular velocity [radians/second]
rps = rpm * pi/30; % shaft speed [radians/second]


% this is just a hard coded control of the rotor speed
t_trans = 2;
t_start = 1;
t_ss    = 10;
t_up    = t_start + t_trans;
t_down  = t_up + t_ss;

if     t < t_start
    % shut-down
    d_dt_theta = 0;
    
elseif t > t_start && t <= t_up
    % ramp up
%     d_dt_theta = rps/t_trans * t;
    d_dt_theta = max(+(rps/t_trans)*(t - t_up) + rps, 0);
    
elseif t > t_up && t <= t_down
    % steady state
    d_dt_theta = rps;
    
else 
    % ramp down
    d_dt_theta = max(-(rps/t_trans)*(t - t_down) + rps, 0);
    
end

% fprintf(1, '[shaftSpeed] time  = %g [s]\n', t);
% fprintf(1, '[shaftSpeed] theta = %g [rad]\n', theta);

end % shaftSpeed