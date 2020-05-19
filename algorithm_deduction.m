clear;
clc;
syms u_j v_j w_j u_w v_w w_w
syms rou A R 
syms C_T_stat K_z sigma a C_D_0 theta_0 lambda mu lambda_stat K_D
C_Tj = C_T_stat + K_z*(w_j - w_w)/(R*abs(w_j));
lambda_j = lambda_stat - 4/(sigma*a)*K_z*(w_j - w_w)/(R*abs(w_j));
mu_j = 1/(abs(w_j)*R)*sqrt((u_j - u_w)^2 + (v_j - v_w)^2);
C_Rmj = sigma*a*mu_j/8*(lambda_j - 4/3);
C_Qj = sigma*C_D_0/8*(1+mu_j^2) + sigma*a*lambda_j*(theta_0/6 - lambda_j/4);
C_H_j = K_D * mu_j;
F_X_j = -rou*A*R^2*(u_j - u_w)/sqrt((u_j-u_w)^2 + (v_j - v_w)^2)*C_H_j*w_j^2
F_Y_j = -rou*A*R^2*(v_j - v_w)/sqrt((u_j-u_w)^2 + (v_j - v_w)^2)*C_H_j*w_j^2
F_Z_j = -rou*A*R^2*C_Tj*w_j^2
L_j = -sign(w_j)*rou*A*R^3*(u_j - u_w)/sqrt((u_j-u_w)^2 + (v_j - v_w)^2)*C_Rmj*w_j^2
M_j = -sign(w_j)*rou*A*R^3*(v_j - v_w)/sqrt((u_j-u_w)^2 + (v_j - v_w)^2)*C_Rmj*w_j^2
N_j = -sign(w_j)*rou*A*R^3*C_Qj*w_j^2
-A*R^3*rou*w_j^2*sign(w_j)*((C_D_0*sigma*(((u_j - u_w)^2 + (v_j - v_w)^2)/(R^2*abs(w_j)^2) + 1))/8 + a*sigma*(lambda_stat - (4*K_z*(w_j - w_w))/(R*a*sigma*abs(w_j)))*(theta_0/6 - lambda_stat/4 + (K_z*(w_j - w_w))/(R*a*sigma*abs(w_j))))