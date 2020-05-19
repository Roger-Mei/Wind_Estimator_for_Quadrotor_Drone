%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Aohan (Roger) Mei
% Date: 05/12/2020
% E-mail: rogermei@umich.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
%% -------------------------Specification--------------------------------%%
% In this program we will test our wind estimation algorithm under hovering
% mode. We adopt the LLC of Mellinger et al's paper. We adopt eq(9) of
% Gabriele's paper, which is using translational dynamics to estimate the
% wind velocity.
% In this project, we will use 4 frames:
% World frame 1(W1): normal 3-D coordinate frame for implementing Mellinger LLC
% World frame 2(W2): x-north, y-east, z-down coordinate frame for wind
% estimation.
% Body frame according to world frame 1 (B1): frame after applying rotation matrix 
% on World frame 1.
% Body frame according to world frame 2 (B2): frame after applying rotation
% matrix on World frame 2.
%% ---------------------------Initialization-----------------------------%%
% Rotation matrix from world frame 1 to world frame 2
R_12 = [-1, 0, 0;
        0, 1, 0;
        0, 0, -1];
% Rotation matrix from world frame 2 to world frame 1
R_21 = R_12^(-1);
% Load the parameter object
P = model_parameter();
% Set up time span
tspan = [0 60];
dt = 0.01;
% Initialize system state
% X = [x,y,z,roll,pitch,yaw,x_dot,y_dot,z_dot,p_W2,q_W2,r_W2]'
% Note that the state is based on world frame 2!!!!!
X = [0; 0; -5; 0; 0; 0; 0; 0; 0; 0; 0; 0];
% Initialize X_dot (world frame 2)
X_dot = zeros(12,1);
% Define the target trajectory r_T (world frame 1)
r_T = [0; 0; 5];
% Define the target velocity r_T_dot (world frame 1)
r_T_dot = zeros(3,1);
% Define the target acceleration r_T_2dot (world frame 1)
r_T_2dot = zeros(3,1);
% Define the target angular velocity omega_T_B (world frame 1)
omega_T_B = zeros(3,1);
% Set the initial estimated wind velocity to be zero (world frame 1)
d_w_est_W1 = zeros(3,1);
% Set up the data saver
% The data saver's format would be: time stamp, d_w_est(1), d_w_est(2),
% d_w_est(3), d_w(1), d_w(2), d_w(3)
data_saver = zeros(tspan(2)/dt+1, 7);
for t = tspan(1):dt:tspan(2)
    %% ----------------Construct rotation matrix R_B_W2------------------%%
    % Construct rotation matrix R_B_W2 which rotate world frame 2 to body
    % frame 2. Note that the roll, pitch, yaw angle are defined according to
    % world frame 2.
    roll = X(4);
    pitch = X(5);
    yaw = X(6);
    % The rotation matrix is constructed with the order Z-Y-X
    R_B2_W2 = [cos(yaw)*cos(pitch), -sin(yaw)*cos(roll)+cos(yaw)*sin(pitch)*sin(roll), sin(roll)*sin(yaw)+cos(yaw)*sin(pitch)*cos(roll);
             sin(yaw)*cos(pitch), cos(yaw)*cos(roll)+sin(yaw)*sin(pitch)*sin(roll), -cos(yaw)*sin(roll)+sin(yaw)*sin(pitch)*cos(roll);
             -sin(pitch), cos(pitch)*sin(roll), cos(pitch)*cos(roll)];
    %% ----------------------Set up wind model---------------------------%%
    % Note that the wind model is created in world frame 1.
    u_w = 2*sin(2*t); % x-direction
    v_w = 2*sin(1.5*t); % y-direction
    w_w = 0.5*sin(0.2*t); % z-direction
    d_w_W1 = [u_w; v_w; w_w];
    % Transfer the wind model in world frame 2
    d_w_W2 = R_12*d_w_W1;
    d_w_est_W2 = R_12*d_w_est_W1;
    % Calculate the wind velocity with respect to world frame 2 in body
    % frame
    d_w_B2 = R_B2_W2^(-1)*d_w_W2;
    d_w_est_B2 = R_B2_W2^(-1)*d_w_est_W2;
    %% ----------------------Calculate inputs----------------------------%%
    % Get the input U_B1 in body frame 1.
    % The format of U_B1 is [F_z, M_x, M_y, M_z]
    U_B1 = get_input(X, R_21, P, r_T, r_T_dot, r_T_2dot, omega_T_B);
    %% ---------------Update parameters for wind estimation--------------%%
    c = zeros(4,1);
    s = zeros(4,1);
    c_rotor = zeros(4,1);
    % u_B2, v_B2, w_B2 is the wind velocity in body frame B2
    u_B2 = zeros(4,1);
    v_B2 = zeros(4,1);
    w_B2 = zeros(4,1);
    for i = 1:4
        % Update the angle coefficient.c := cos, s := sin.
        c(i) = cos(pi/2*(i-1) + pi/4);
        s(i) = sin(pi/2*(i-1) + pi/4);
        c_rotor(i) = cos(pi*(i-1));
        % Calculate each rotor's translational velocity in body frame B2
        temp_vec1 = [P.l*c(i); P.l*s(i); P.h];
        % Calculate angular velocity in Body frame B2
        temp_vec2 = R_B2_W2'*X(10:12);
        % Calculate the linear velocity in Body frame B2
        temp_vec3 = R_B2_W2'*X(7:9);
        rotor_velocity = cross(temp_vec2, temp_vec1) + temp_vec3;
        % Calculate the rotor velocity in Body frame B2
        u_B2(i) = rotor_velocity(1);
        v_B2(i) = rotor_velocity(2);
        w_B2(i) = rotor_velocity(3);
    end
    % Construct a matrix M_map which maps the input U to the square of
    % rotation speed of each rotor
    K_f = P.rou*P.A*P.R^2*P.C_Tstat;
    K_m = P.rou*P.A*P.R^3*(P.sigma*P.C_D0/8 + P.lambda_stat*P.sigma*P.a*(P.theta_0/6 - P.lambda_stat/4));
    M_map = [K_f, K_f, K_f, K_f;
             K_f*P.l*s(1), K_f*P.l*s(2), K_f*P.l*s(3), K_f*P.l*s(4);
             K_f*P.l*c(1), K_f*P.l*c(2), K_f*P.l*c(3), K_f*P.l*c(4);
             K_m, -K_m, K_m, -K_m];
    % Each rotor's rotation speed
    temp_omega = M_map^(-1)*U_B1;
    % Format of temp_omega: [omega_1^2; omega_2^2; omega_3^2; omega_4^2;]
    for i = 1:length(temp_omega)
        if temp_omega(i) < 200^2
            temp_omega(i) = 200^2;
        elseif temp_omega >= 500^2
            temp_omega(i) = 500^2;
        else
            temp_omega(i) = temp_omega(i);
        end
    end
    temp_omega = sqrt(temp_omega);
    % This omega is according to body frame 2
    omega = temp_omega.*c_rotor;
    %% -----------------------Estimate wind velocity---------------------%%
    d_w_est_B2 = wind_estimation(omega, d_w_B2, P, d_w_est_B2, dt, u_B2, v_B2, w_B2);
    d_w_est_W2 = R_B2_W2 * d_w_est_B2;
    d_w_est_W1 = R_21*d_w_est_W2;
    %% -----------------------Update Dynamics----------------------------%%
    % Note that we update system dynamics in Body frame B2 and you should
    % change them back to world frame W2
    % Update parameters
    C_T = zeros(4,1);
    lambda = zeros(4,1);
    mu = zeros(4,1);
    C_Rm = zeros(4,1);
    C_Q = zeros(4,1);
    C_H = zeros(4,1);
    u_w_B2 = d_w_B2(1);
    v_w_B2 = d_w_B2(2);
    w_w_B2 = d_w_B2(3);
    for j = 1:4
        C_T(j) = P.C_Tstat + P.K_z*(w_B2(j) - w_w_B2)/(P.R*abs(omega(j)));
        lambda(j) = P.lambda_stat - 4/(P.sigma*P.a)*P.K_z*(w_B2(j) - w_w_B2)/(P.R*abs(omega(j)));
        mu(j) = 1/(abs(omega(j))*P.R)*sqrt((u_B2(j) - u_w_B2)^2 + (v_B2(j) - v_w_B2)^2);
        C_Rm(j) = P.sigma*P.a*mu(j)/8*(lambda(j) - 4/3);
        C_Q(j) = P.sigma*P.C_D0/8*(1+mu(j)^2) + P.sigma*P.a*lambda(j)*(P.theta_0/6 - lambda(j)/4);
        C_H(j) = P.K_D*mu(j);
    end
    % Calculate aerodynamic forces
    F_x = zeros(4,1);
    F_y = zeros(4,1);
    F_z = zeros(4,1);
    L = zeros(4,1);
    M = zeros(4,1);
    N = zeros(4,1);
    temp_L_aero = zeros(4,1);
    temp_M_aero = zeros(4,1);
    temp_N_aero = zeros(4,1);
    for idx = 1:4
        % We set this condition just in case of the singularity issue
        % caused by denomenator being zero.
        if sqrt((u_B2(idx)-u_w)^2 + (v_B2(idx) - v_w)^2) == 0
            F_x(idx) = 0;
            F_y(idx) = 0;
            L(idx) = 0;
            M(idx) = 0;
        else
            F_x(idx) = -P.rou*P.A*P.R^2*(u_B2(idx) - u_w_B2)/sqrt((u_B2(idx)-u_w_B2)^2 + (v_B2(idx) - v_w_B2)^2)*C_H(idx)*omega(idx)^2;
            F_y(idx) = -P.rou*P.A*P.R^2*(v_B2(idx) - v_w_B2)/sqrt((u_B2(idx)-u_w_B2)^2 + (v_B2(idx) - v_w_B2)^2)*C_H(idx)*omega(idx)^2;
            L(idx) = -sign(omega(idx))*P.rou*P.A*P.R^3*(u_B2(idx) - u_w_B2)/sqrt((u_B2(idx)-u_w_B2)^2 + (v_B2(idx) - v_w_B2)^2)*C_Rm(idx)*omega(idx)^2;
            M(idx) = -sign(omega(idx))*P.rou*P.A*P.R^3*(v_B2(idx) - v_w_B2)/sqrt((u_B2(idx)-u_w_B2)^2 + (v_B2(idx) - v_w_B2)^2)*C_Rm(idx)*omega(idx)^2;
        end
        F_z(idx) = -P.rou*P.A*P.R^2*C_T(idx)*omega(idx)^2;
        N(idx) = -sign(omega(idx))*P.rou*P.A*P.R^3*C_Q(idx)*omega(idx)^2;
        
        temp_L_aero(idx) = L(idx) + F_z(idx)*P.l*s(idx) - P.h*F_y(idx);
        temp_M_aero(idx) = M(idx) - F_z(idx)*P.l*c(idx) + P.h*F_x(idx);
        temp_N_aero(idx) = N(idx) + F_y(idx)*P.l*c(idx) - F_x(idx)*P.l*s(idx);
    end
    F_Xaero = sum(F_x);
    F_Yaero = sum(F_y);
    F_Zaero = sum(F_z);
    % F_aero is the vector of external aerodynamic forces in body frame
    F_aero = [F_Xaero; F_Yaero; F_Zaero];
    L_aero = sum(temp_L_aero);
    M_aero = sum(temp_M_aero);
    N_aero = sum(temp_N_aero);
    % tau_aero is the external aerodynamic moments in the body frame
    tau_aero = [L_aero; M_aero; N_aero];
    % x_dot, y_dot, z_dot
    X_dot(1:3) = X(7:9);
    % calculate angular velocity in body frame
    omega_B2 = R_B2_W2^(-1)*X(10:12);
    p = omega_B2(1);
    q = omega_B2(2);
    r = omega_B2(3);
     % roll_dot. pitch_dot, yaw_dot
    X_dot(4:6) = [p + tan(pitch)*(q*sin(roll) + r*cos(roll));
                  q*cos(roll) - r*sin(roll);
                  (q*sin(roll) + r*cos(roll))/cos(pitch)];
    % x_2dot, y_2dot, z_2dot in Body frame B2
    acl_B2 = -cross([p;q;r], R_B2_W2^(-1)*X(7:9)) + F_aero/P.m + R_B2_W2^(-1)*[0;0;P.g];
    % Linear acceleration in World frame 2
    X_dot(7:9) = R_B2_W2*acl_B2;
    % Calculate angular accelearation in Body frame B2
    an_acl_B2 = P.I^(-1)*(-cross([p;q;r], P.I*[p;q;r]) + tau_aero);
    % Calculate angular acceleration in World frame W2
    X_dot(10:12) = R_B2_W2*an_acl_B2;
    % Integrate forward with respect to world frame 2
    X = X + X_dot*dt;
    %% --------------------Construct the data saver----------------------%%
    index = t/dt + 1;
    index = round(index);
    data_saver(index, 1) = t;
    data_saver(index, 2) = d_w_est_W1(1);
    data_saver(index, 3) = d_w_est_W1(2);
    data_saver(index, 4) = d_w_est_W1(3);
    data_saver(index, 5) = d_w_W1(1);
    data_saver(index, 6) = d_w_W1(2);
    data_saver(index, 7) = d_w_W1(3);
end
%% --------------------------Plot out the data --------------------------%%
% Plot the x-axis wind
subplot(6,1,1);
x = data_saver(:,1);
y1_est = data_saver(:,2);
y1_act = data_saver(:,5);
plot(x,y1_est,'r','LineWidth',2.5);
hold on;
plot(x, y1_act, 'b','LineWidth',2.5);
grid on;
title('x-axis wind')
legend('estimated wind','actual wind');
% Plot the y-axis wind
subplot(6,1,2);
x = data_saver(:,1);
y2_est = data_saver(:,3);
y2_act = data_saver(:,6);
plot(x,y2_est,'r','LineWidth',2.5);
hold on;
plot(x, y2_act, 'b','LineWidth',2.5);
grid on;
title('y-axis wind')
% Plot the z-axis wind
subplot(6,1,3);
x = data_saver(:,1);
y3_est = data_saver(:,4);
y3_act = data_saver(:,7);
plot(x,y3_est,'r','LineWidth',2.5);
hold on;
plot(x, y3_act, 'b','LineWidth',2.5);
grid on;
title('z-axis wind')

subplot(6,1,4);
x = data_saver(:,1);
y1_est = data_saver(:,2);
y1_act = data_saver(:,5);
y1_error = y1_act - y1_est;
plot(x,y1_error,'r','LineWidth',2.5);
grid on;
title('x-axis error')

subplot(6,1,5);
x = data_saver(:,1);
y2_est = data_saver(:,3);
y2_act = data_saver(:,6);
y2_error = y2_act - y2_est;
plot(x,y2_error,'r','LineWidth',2.5);
grid on;
title('y-axis error')


subplot(6,1,6);
x = data_saver(:,1);
y3_est = data_saver(:,4);
y3_act = data_saver(:,7);
y3_error = y3_act - y3_est;
plot(x,y3_error,'r','LineWidth',2.5);
title('z-axis error')
