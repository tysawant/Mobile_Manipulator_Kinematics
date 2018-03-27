%% Part A (Robot Arm Kinematics)
%% Question 1 Defining the variables
close all
clear all
clc
syms t1 t2 t3 b d c e f g
theta = [0 0 t1-(pi/2) t2 t3-pi/2 0 0]
dh = [0 (b+d) 0 0 0 0 g]
a = [-c 0 e f 0 0 0]
alpha = [0 -pi/2 0 0 0 -pi/2 0]
T = cell(7,1);
T_mult = cell(7,1);
origins = [];
M = eye(4,4);

%% Question 2
for i = 1:7
    T{i} = simplify(dhparam2matrix(theta(i), dh(i), a(i), alpha(i)));
    T{i}
    M = M*T{i}
    T_mult{i} = M;
end

%% Question 3
syms x_dot y_dot z_dot wx wy wz theta1_dot theta2_dot theta3_dot
%T_mult{8} = subs(T_mult{8},[t1,t2,t3],[th1(t),th2(t),th3(t)])
cart_vel = [diff(T_mult{7}(1:3,4),t1) diff(T_mult{7}(1:3,4),t2) diff(T_mult{7}(1:3,4),t3)]
x = [0 0 0];
y = [1 1 1];
z = [0 0 0];
ang_vel = [x; y; z]
j_q = [cart_vel; ang_vel]
tip_velocities = [x_dot; y_dot; z_dot; wx; wy; wz] == j_q*[theta1_dot; theta2_dot; theta3_dot]
reduced_jacobian = [j_q(1,:);j_q(3,:);j_q(5,:)]
reduced_tip_velocities = [x_dot; z_dot; wy] == reduced_jacobian*[theta1_dot; theta2_dot; theta3_dot]

%% Question 4

syms tau1 tau2 tau3 Fx_tip Fy_tip Fz_tip Mx_tip My_tip Mz_tip

joint_torques = [tau1; tau2; tau3] == j_q.'*[Fx_tip; Fy_tip; Fz_tip; Mx_tip; My_tip; Mz_tip]

%% Question 5
disp('Position is in mm')
disp('Velocities are mm/sec')
disp('Torques are N/mm')
Q5_a = eval(subs(simplify(T_mult{7}),[b,c,d,e,f,g,t1,t2,t3],[424,300,380,328,323,82.4,pi/9,pi/2,pi/6]))

Q5_b = eval(subs(rhs(tip_velocities),[b,c,d,e,f,g,t1,t2,t3,theta1_dot,theta2_dot,theta3_dot],[424,300,380,328,323,82.4,pi/9,pi/2,pi/6,pi/6,pi/6,pi/6]))

Q5_c = eval(subs(rhs(joint_torques),[b,c,d,e,f,g,t1,t2,t3,Fx_tip, Fy_tip, Fz_tip, Mx_tip, My_tip, Mz_tip],[424,300,380,328,323,82.4,pi/9,pi/2,pi/6,30,0,0,0,0,0]))

%% Question 6
Q6_tip_velocities = [0; -100; 0]
disp('Joint Velocities in degrees/sec:')
Q6_joint_velocities = [theta1_dot; theta2_dot; theta3_dot] == vpa(rad2deg(eval(subs(inv(reduced_jacobian)*Q6_tip_velocities,[b,c,d,e,f,g,t1,t2,t3],[424,300,380,328,323,82.4,pi/9,pi/2,pi/6]))))

%% Part B (Mobile Platform Kinematics)
%% Question 7
% castor_l = distance of castor from robot frame
% dist = offset of castor
syms xb_dot yb_dot tb_dot alpha alpha_dot a r wl wr dist beta beta_dot r_sw w_sw
zeta_w = [xb_dot; yb_dot; tb_dot]
R_rob = [cos(alpha) sin(alpha) 0; -sin(alpha) cos(alpha) 0; 0 0 1];

left_wheel_rolling = eval([sin(pi/2+0) -cos(pi/2+0) -(a/2)*cos(0)]*R_rob*zeta_w - r*wl) == 0
left_wheel_sliding = eval([cos(pi/2+0)  sin(pi/2+0) -(a/2)*sin(0)]*R_rob*zeta_w) == 0

right_wheel_rolling = eval([sin(-pi/2+pi) -cos(-pi/2+pi) -(a/2)*cos(pi)]*R_rob*zeta_w - r*wr) == 0
right_wheel_sliding = eval([cos(-pi/2+pi)  sin(-pi/2+pi) -(a/2)*sin(pi)]*R_rob*zeta_w) == 0

wc = 0.5*(wl+wr)
castor_wheel_rolling = eval([sin(pi+beta) -cos(pi+beta) -c*cos(beta)]*R_rob*zeta_w - r*wc) == 0
castor_wheel_sliding = eval([cos(pi+beta)  sin(pi+beta) c*sin(beta)]*R_rob*zeta_w - r_sw*w_sw) == 0

%% Question 8

J1_f = eval([sin(pi/2+0) -cos(pi/2+0) -(a/2)*cos(0); sin(-pi/2+pi) -cos(-pi/2+pi) -(a/2)*cos(pi); sin(pi+beta) -cos(pi+beta) -c*cos(beta)])
C = eval([cos(pi/2+0)  sin(pi/2+0) -(a/2)*sin(0); cos(-pi/2+pi)  sin(-pi/2+pi) -(a/2)*sin(pi); cos(pi+beta)  sin(pi+beta) c*sin(beta)])
J2 = [r*wl; r*wr; r*wc]
constraint_mat = vpa([J1_f;C])
disp('For the reduced constraint matrix, we take 2 rolling constraints and one sliding constraint as the sliding constraints of both driven wheels is the same, and the castor wheel is not a driven wheel')
J2_phi = [r*wl; r*wr; 0]
reduced_constraint_mat = [J1_f(1:2,:);C(1,:)]
constraint_eq = vpa(reduced_constraint_mat*R_rob*zeta_w) == J2_phi

%% Question 9

velocity_kinematics_base = zeta_w == vpa(simplify(eval(inv(R_rob)*inv(reduced_constraint_mat)*J2_phi)))

%% Question 10
disp('velocities are mm/sec- linear and rad/sec - angular')
Q_10_vel = vpa(eval(subs(velocity_kinematics_base,[a,r,alpha,wl,wr],[507,143,pi/6,pi,2*pi])))

%% Part C (Hybrid System Kinematics)
%% Question 11
syms px py pz
Transf_r_w = [cos(alpha) -sin(alpha) 0 px; sin(alpha) cos(alpha) 0 py; 0 0 1 pz; 0 0 0 1];
Transf_r_t = T_mult{7}
tip_pos_w = Transf_r_w*Transf_r_t

%% Question 12
disp('The tip position is in mm')
Transf_r_w_num = vpa(subs(Transf_r_w,[alpha,px,py,pz], [pi/6, 2500, 1500, 0]))
Transf_r_t = subs(simplify(T_mult{7}),[b,c,d,e,f,g,t1,t2,t3],[424,300,380,328,323,82.4,pi/9,pi/2,pi/6])
tip_pos_w_numeric = eval(Transf_r_w_num*Transf_r_t)

%% Question 13

J_A_w = [diff(tip_pos_w(1:3,4),t1) diff(tip_pos_w(1:3,4),t2) diff(tip_pos_w(1:3,4),t3)]
T3_mult_w = Transf_r_w*T_mult{3}
T4_mult_w = Transf_r_w*T_mult{4};
T5_mult_w = Transf_r_w*T_mult{5};

J_B_w = [T3_mult_w(1:3,3),T4_mult_w(1:3,3),T5_mult_w(1:3,3)];
J_w_r = [J_A_w; J_B_w]

J_mob_r = Transf_r_w(1:3,1:3)*inv(reduced_constraint_mat)*r;

J_mob_w = [J_mob_r(1:2,1:2);0,0;0,0;0,0;J_mob_r(3,1:2)]

J_mob_f = [J_mob_r, J_mob_w]

Q_13_eq = [Fx_tip; Fy_tip; Fz_tip; Mx_tip; My_tip; Mz_tip] == J_mob_f*[theta1_dot; theta2_dot; theta3_dot; wl; wr]

%% Question 14
syms tau_lw tau_rw
disp('Torques are N/mm')
torques = [tau1; tau2; tau3; tau_lw; tau_rw] == J_mob_f.'* [Fx_tip; Fy_tip; Fz_tip; Mx_tip; My_tip; Mz_tip];
torques = vpa(subs(torques,[b,c,d,e,f,g,t1,t2,t3,Fx_tip, Fy_tip, Fz_tip, Mx_tip, My_tip, Mz_tip,r,alpha],[424,300,380,328,323,82.4,pi/9,pi/2,pi/6,30*cos(pi/6),30*sin(pi/6),0,0,0,0,143,pi/6]))
