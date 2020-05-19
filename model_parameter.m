function P = model_parameter()
    % Constant parameters
    P.R = 0.10;
    P.A = pi*P.R^2;
    P.l = 0.185;
    P.h = -0.025;
    P.g = 9.81;
    P.theta_0 = 23.9;
    P.m = 0.472;
    P.Ixx = 3.56*10^(-3);
    P.Iyy = 4.02*10^(-3);
    P.Izz = 7.12*10^(-3);
    P.I = [P.Ixx 0 0;
        0 P.Iyy 0;
        0 0 P.Izz];
    P.rou = 1.25;
    P.sigma = 0.1114;
    P.a = 4.6542;
    P.C_D0 = 2.15;
    P.lambda_stat = 0.1056;
    P.C_Tstat = 0.0223;
    P.K_D = 0.06;
    P.K_z = 0.09;
    P.gamma_a = 40;
    P.alpha_a = 0.9;
    P.gamma_g = 60;
    P.italian_l = 40;
    P.gamma_g_prime = 40;
    P.alpha_prime = 40;
    P.italian_l_prime = 20;
    P.epsilon_a = 0.052; % Noise of accelorameter
    P.epsilon_g = 2.5/180*pi; % Noise of gyroscope 
    % Tunable parameters
    P.K_p = 2*eye(3);
    P.K_v = 0.5*eye(3);
    P.K_R = eye(3);
    P.K_w = 0.03*eye(3);
end