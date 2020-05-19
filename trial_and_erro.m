clear; clc;
pf = [20;20;20];
v0 = [0;0;0];
a0 = [0;0;0];
time_horizon = 10;
sample_time = 0.05;
[T,U,Z] = generate_spline_desired_position(pf,v0,a0,time_horizon,sample_time);
plot3(Z(1,:),Z(2,:),Z(3,:))
grid on;