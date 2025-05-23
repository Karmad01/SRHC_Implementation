clc; clear; close all;
set(0, 'defaultFigureWindowState', 'maximized');
set(0, 'defaultAxesFontSize', 20);
set(0, 'defaultAxesLineWidth', 1.5);
set(0, 'defaultAxesTickDirMode', 'manual');
set(0, 'defaultAxesTickDir', 'out');
set(0, 'defaultAxesTickLength', [0.005;0.005]);
set(groot, 'defaultAxesBox', 'on');
set(groot, 'defaultAxesFontName','Century');
set(groot, 'defaultLegendFontName','Century');
set(groot, 'defaultTextinterpreter','latex');  
set(groot, 'defaultTextfontsize',20);
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultColorbarTickLabelinterpreter','latex'); 
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultLineMarkerSize', 25);
set(groot, 'defaultConstantlineLineWidth', 2);
%% System
A_s = 0.9;
A_o = [1 0 0;0 0 -1;0 1 0];
B_s = 0;
B_o = [1;0;1];
n_s = length(A_s);
n_o = length(A_o);
n = n_s+n_o;
A = [A_s 0 0 0;zeros(3,1) A_o];
B = [B_s;B_o];
C = eye(n);
m = 1;
p = n;
% k = 3;
sigma_x0 = eye(n);
sigma_w = 10.*eye(n);
sigma_v = 10.*eye(p);
Q = 100.*eye(n);
Q_terminal = Q;
R = eye(m);
% N = 5;
% N_c = k;
% phi_max = 1;
U_max = 3.2664;
[K,S,P] = dlqr(A,B,Q,R);
% zeta = 391.3082; % Change this!!!!
% epsilon = 10;

% %% Problem matrices
% A_cal = zeros((N+1)*n,n);
% B_cal = zeros((N+1)*n,N*m);
% D_cal = zeros((N+1)*n,N*n);
% Q_cal = blkdiag(kron(eye(N),Q),Q_terminal);
% R_cal = kron(eye(N),R);
% 
% for ii=1:N+1
%     A_cal((ii-1)*n+1:ii*n,:) = A^(ii-1);
%     for jj=1:N
%         if ii~=1
%             if ii>jj
%                 B_cal((ii-1)*n+1:ii*n,(jj-1)*m+1:jj*m) = A^(ii-jj-1)*B;
%                 D_cal((ii-1)*n+1:ii*n,(jj-1)*n+1:jj*n) = A^(ii-jj-1);
%             end
%         end
%     end
% end
% M_cal = R_cal+B_cal.'*Q_cal*B_cal;
% R_AoBo_Nc = reachability_matrix(A_o,B_o,N_c);
% 
% load("Lambda.mat");
% 
% %% Algorithm 1
% A_cal = A_cal(1:N*n,:);
% B_cal = B_cal(1:N*n,:);
% D_cal = D_cal(1:N*n,:);
% Q_cal = Q_cal(1:N*n,1:N*n);
% Lambda_phi = Lambda_phi(1:N*p);
% Lambda_phiphi = Lambda_phiphi(1:N*p,1:N*p);
% Lambda_e_phi = Lambda_e_phi(:,1:N*p);
% Lambda_w_phi = Lambda_w_phi(1:N*n,1:N*p);

x_norm = zeros(1000,100);
u = zeros(1000,100);
cost = zeros(1000,100);
% z_1 = sdpvar(1,1);
% z_2 = sdpvar(1,1);
% z_3 = sdpvar(1,1);
% eta_t = sdpvar(m*N,1);
% Theta_t = sdpvar(N*m,N*p,'full');

for num_samples = 1:1000
disp(["Starting : ", num2str(num_samples)])
tic
t = 0;
x_hat_t1t = 200.*[1;-1;-1;1];
P_t1t = sigma_x0;
x_t = mvnrnd(x_hat_t1t,sigma_x0,1).';
while t<=101
    % y_t = zeros(N*p,1);
    % y_hat_t = zeros(N*p,1);
    % for ii=0:N_c-1
        y_t = C*x_t+mvnrnd(zeros(length(sigma_v),1),sigma_v,1).';
        x_hat_t1t1 = x_hat_t1t+P_t1t*C.'*((C*P_t1t*C.'+sigma_v)\(y_t-C*x_hat_t1t));
        P_t1t1 = P_t1t-P_t1t*C.'*((C*P_t1t*C.'+sigma_v)\C)*P_t1t;
        x_norm(num_samples,t+1) = norm(x_hat_t1t1);
        % if ii==0
        %     Objective = z_1;
        %     Constraints = [];
        %     x_hat_o = x_hat_t1t1(n_s+1:end);
        %     Lambda_x_phi = Lambda_e_phi+x_hat_t1t1*Lambda_phi.';
        %     c1 = (eta_t+Theta_t*Lambda_phi).'*M_cal*(eta_t+Theta_t*Lambda_phi);
        %     c1 = c1+trace(Theta_t.'*M_cal*Theta_t*(Lambda_phiphi-Lambda_phi*Lambda_phi.'));
        %     c1 = c1+2*x_hat_t1t1.'*A_cal.'*Q_cal*B_cal*eta_t;
        %     c1 = c1+2*trace(Theta_t.'*B_cal.'*Q_cal*(D_cal*Lambda_w_phi+A_cal*Lambda_x_phi));
        %     Constraints = [Constraints,c1<=z_1];
        %     for jj=1:N*m
        %         Constraints = [Constraints,abs(eta_t(jj))+norm(Theta_t(jj,:),1)*phi_max<=U_max];
        %     end
        %     if norm(x_hat_o)>=zeta+epsilon
        %         Constraints = [Constraints,norm(A_o^N_c*x_hat_o+R_AoBo_Nc*eta_t(1:N_c*m))<=z_2];
        %         Constraints = [Constraints,norm(R_AoBo_Nc*Theta_t(1:N_c*m,:),inf)<=z_3];
        %         Constraints = [Constraints,z_2+sqrt(n_o)*phi_max*z_3<=norm(x_hat_o)-zeta-epsilon/2];
        %     end
        %     for k = 1:N
        %         for l = 1:N
        %             if l>k
        %                 Constraints = [Constraints,Theta_t((k-1)*m+1:k*m,(l-1)*p+1:l*p)==0];
        %             end
        %         end
        %     end
        %     options = sdpsettings('solver','sdpt3','verbose',0);
        %     sol = optimize(Constraints,Objective,options);
        %     if sol.problem == 0
        %         disp('Solution Found');
        % 
        %     else
        %         disp(['Problem status: ', sol.info]);
        %         yalmiperror(sol.problem);
        %     end
        % end
        % eta_t_opt = value(eta_t);
        % Theta_t_opt = value(Theta_t);
        % y_hat_t(ii*p+1:(ii+1)*p) = C*x_hat_t1t1;
        % u_t1_vec = eta_t_opt+Theta_t_opt*phi_func(y_t-y_hat_t,phi_max);
        % u_t1 = u_t1_vec(ii*m+1:(ii+1)*m);
        u_t1 = phi_func(-K*x_hat_t1t1,U_max);
        u(num_samples,t+1) = u_t1;
        cost(num_samples,t+1) = x_hat_t1t1.'*Q*x_hat_t1t1+u_t1.'*R*u_t1;
        x_t = A*x_t+B*u_t1+mvnrnd(zeros(length(sigma_w),1),sigma_w,1).';
        x_hat_t1t = A*x_hat_t1t1+B*u_t1;
        P_t1t = A*P_t1t1*A.'+sigma_w;
    %end
    t = t+1;
end
elapsed_time = toc;  % Stop timing
disp("Done")
fprintf('Sample %d completed in %.2f seconds.\n', num_samples, elapsed_time);
end

%% Plotting
avg_state_norm_lqg = mean(x_norm, 1);
time_axis_lqg = 0:(size(x_norm, 2)-1);
std_state_norm_lqg = std(x_norm, 0, 1);
cost1 = cost;
u_lqg = u;

for ii=1:101
    %for jj=1:ii-1
        cost1(:,ii+1) = cost1(:,ii+1)+cost1(:,ii);
    %end
end
avg_cost_lqg = mean(cost1,1);
load("sdpt3_100_traj_km_x_t_x_tt_400_new3.mat","avg_state_norm","cost","u","std_state_norm","time_axis");
cost1 = cost;
for ii=1:101
    %for jj=1:ii-1
        cost1(:,ii+1) = cost1(:,ii+1)+cost1(:,ii);
    %end
end
avg_cost = mean(cost1,1);
time_axis = time_axis-1;
avg_u_lqg = mean(u_lqg,1);
avg_u = mean(u,1);
%%
figure;
hold on;
plot(time_axis,avg_state_norm,'r-','LineWidth',2);
plot(time_axis,std_state_norm,'r--','LineWidth',1);
plot(time_axis_lqg, avg_state_norm_lqg, 'b-', 'LineWidth', 2);
plot(time_axis_lqg, std_state_norm_lqg, 'b--', 'LineWidth', 1);
hold off

xlabel('Time step');
ylabel('State norm');
title('Average and Std Dev of State Norm over Time');
legend('Avg: SRHC', 'Std: SRHC','Avg: Clipped LQG', 'Std: Clipped LQG',Location='southeast');
grid on;
ylim([0,max([avg_state_norm_lqg avg_state_norm])]);
xlim([0,102]);

figure;
hold on
plot(time_axis,avg_cost,'r-','LineWidth',2);
plot(time_axis_lqg, avg_cost_lqg, 'b-', 'LineWidth', 2);
hold off
legend('SRHC','Clipped LQG',Location='southeast')
xlabel('Time step');
ylabel('Cost');
title('Average Cost over Time');
grid on;
ylim([0,max([avg_cost_lqg avg_cost])]);
xlim([0,102]);

figure
hold on
plot(time_axis,avg_u,'r-','LineWidth',2);
plot(time_axis_lqg,avg_u_lqg,'b-','LineWidth',2);
hold off
%yline(U_max,'--');
%yline(-U_max,'--')
legend('SRHC','Clipped LQG')
xlabel('Time step');
ylabel('Input');
title('Average Control Input over Time');
grid on;
ylim([min([avg_u avg_u_lqg]),max([avg_u avg_u_lqg])])
xlim([0,102]);