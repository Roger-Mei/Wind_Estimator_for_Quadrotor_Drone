%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Aohan (Roger) Mei
% Date: 05/12/2020
% E-mail: rogermei@umich.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = get_input(X_2, R_21, P, r_T, r_T_dot, r_T_2dot, omega_T_B)
    %% -----------------------------Input--------------------------------%%
    % X_2 is the state of the system in world frame 2
    % X = [x,y,z,roll,pitch,yaw,x_dot,y_dot,z_dot,p,q,r]'
    % R_B_W is the rotation matrix
    % R_21 is the rotation matrix which rotates world frame 1 to world
    % frame 2.
    % P is the parameter object
    % r_T is the desired trajectory
    % r_T_dot is the desired velocity
    % r_T_2dot is the desired acceleration
    U = zeros(4,1);
    %% ----------------Construct states in world frame 1-----------------%%
    X = zeros(12,1);
    X(1:3) = R_21*X_2(1:3);
    X(4:6) = R_21*X_2(4:6);
    X(7:9) = R_21*X_2(7:9);
    X(10:12) = R_21*X_2(10:12);
    %% -------Construct new rotation matrix in world frame 1-------------%%
    roll = X(4);
    pitch = X(5);
    yaw = X(6);
    R_B_W = [cos(yaw)*cos(pitch), -sin(yaw)*cos(roll)+cos(yaw)*sin(pitch)*sin(roll), sin(roll)*sin(yaw)+cos(yaw)*sin(pitch)*cos(roll);
             sin(yaw)*cos(pitch), cos(yaw)*cos(roll)+sin(yaw)*sin(pitch)*sin(roll), -cos(yaw)*sin(roll)+sin(yaw)*sin(pitch)*cos(roll);
             -sin(pitch), cos(pitch)*sin(roll), cos(pitch)*cos(roll)];
    %% --------------------------Calculate u1----------------------------%%
    % Position error
    r = X(1:3);
    e_p = r - r_T;
    % Velocity error
    r_dot = X(7:9);
    e_v = r_dot - r_T_dot;
    F_des = -P.K_p*e_p - P.K_v*e_v + P.m*[0;0;P.g] + P.m*r_T_2dot;
    z_B = R_B_W(:,3);
    thrust_out = F_des' * z_B;
    %% -----------------------Calculate u2,u3,u4-------------------------%%
    % Construct R_des
    z_B_des = F_des/norm(F_des);
    x_c = [cos(X(6)),sin(X(6)),0]';
    y_B_des = cross(z_B_des, x_c)/norm(cross(z_B_des, x_c));
    x_B_des = cross(y_B_des, z_B_des);
    R_des = [x_B_des, y_B_des, z_B_des];
    % Calculate the error on orientation
    e_R = 1/2*unskew(R_des'*R_B_W - R_B_W'*R_des);
    % Calculate angular velocity error
    e_omega = X(10:12) - omega_T_B;
    moment_out = -P.K_R*e_R - P.K_w*e_omega;
    U = [thrust_out; moment_out];
end
%% Helper function for unskew symetric matrix
function result = unskew(M)
    result = [M(3,2); M(1,3); M(2,1)];
end