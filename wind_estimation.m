%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Aohan (Roger) Mei
% Date: 05/10/2020
% E-mail: rogermei@umich.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d_w_est = wind_estimation(omega, d_w, P, d_w_est, dt, u, v, w)
    %% ----------------------------Input---------------------------------%%
    % X is the state of the drone system
    % omega is the vector of rotors angular rotational velocities
    % d_w is the wind velocity with respect to body frame B2 
    % P is the parameter object 
    % d_w_est is the estimated wind velocity from last loop with respect to
    % body frame B2
    % dt is the time interval
    % F_aero is the aero-dynamic force
    %% -------------------Construct f_0a and Omega_a---------------------%%
    f_0a = zeros(3,1);
    Omega_a = zeros(3,3);
    for idx1 = 1:4
        temp_f_0a = [-P.A*P.K_D*P.R*P.rou*omega(idx1)^2/(P.m*abs(omega(idx1)))*u(idx1);
                     -P.A*P.K_D*P.R*P.rou*omega(idx1)^2/(P.m*abs(omega(idx1)))*v(idx1);
                     -P.A*P.K_z*P.R*P.rou*omega(idx1)^2/(P.m*abs(omega(idx1)))*w(idx1) - P.A*P.R^2*P.rou*P.C_Tstat*omega(idx1)^2/P.m];
        f_0a  = f_0a + temp_f_0a;
        temp_Omega_a = [P.A*P.K_D*P.R*P.rou*omega(idx1)^2/(P.m*abs(omega(idx1))), 0, 0;
                        0, P.A*P.K_D*P.R*P.rou*omega(idx1)^2/(P.m*abs(omega(idx1))), 0;
                        0, 0, P.A*P.K_z*P.R*P.rou*omega(idx1)^2/(P.m*abs(omega(idx1)))];
        Omega_a = Omega_a + temp_Omega_a;
    end
    %% --------------------Construct wind estimator----------------------%%
    % Changing rate of the estimated wind
    d_w_est_dot = wind_est_odefun(d_w_est,d_w, Omega_a, P.epsilon_a, P.gamma_a, P.alpha_a);
    d_w_est = d_w_est_dot * dt + d_w_est;
end
%% ------------Construct the ode function for wind estimator-------------%%
function d_w_est_dot = wind_est_odefun(d_w_est, d_w, Omega_a, epsilon_a, gamma_a, alpha_a)
    noise = epsilon_a*randn(3,1);
    e_a = Omega_a*(d_w - d_w_est) + noise;
    d_w_est_dot = gamma_a*Omega_a'*abs(e_a).^alpha_a.*sign(e_a);
end