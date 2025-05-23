clc; clear; close all;

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
k = 3;
sigma_x0 = eye(n);
sigma_w = 10.*eye(n);
sigma_v = 10.*eye(p);
Q = eye(n);
Q_terminal = Q;
R = eye(m);
N = 5;
N_c = k;
phi_max = 1;
U_max = 453;
zeta = 391.3082; % Change this!!!!
epsilon = 10;
alpha = 1000;
beta = 5*(453^2);

%% Problem matrices
A_cal = zeros((N+1)*n,n);
B_cal = zeros((N+1)*n,N*m);
D_cal = zeros((N+1)*n,N*n);
Q_cal = blkdiag(kron(eye(N),Q),Q_terminal);
R_cal = kron(eye(N),R);
sigma_w_blk = kron(eye(N),sigma_w);

for ii=1:N+1
    A_cal((ii-1)*n+1:ii*n,:) = A^(ii-1);
    for jj=1:N
        if ii~=1
            if ii>jj
                B_cal((ii-1)*n+1:ii*n,(jj-1)*m+1:jj*m) = A^(ii-jj-1)*B;
                D_cal((ii-1)*n+1:ii*n,(jj-1)*n+1:jj*n) = A^(ii-jj-1);
            end
        end
    end
end
M_cal = R_cal+B_cal.'*Q_cal*B_cal;
R_AoBo_Nc = reachability_matrix(A_o,B_o,N_c);
S_cal = eye((N+1)*n);
S_tilde_cal = eye(N*m);
L_cal = zeros((N+1)*n,1);

load("Lambda.mat");

%% Algorithm 1
A_cal = A_cal(1:N*n,:);
B_cal = B_cal(1:N*n,:);
D_cal = D_cal(1:N*n,:);
Q_cal = Q_cal(1:N*n,1:N*n);
Lambda_phi = Lambda_phi(1:N*p);
Lambda_phiphi = Lambda_phiphi(1:N*p,1:N*p);
Lambda_e_phi = Lambda_e_phi(:,1:N*p);
Lambda_w_phi = Lambda_w_phi(1:N*n,1:N*p);
S_cal = S_cal(1:N*n,1:N*n);
L_cal = L_cal(1:N*n);

x_norm = zeros(10,51);
u = zeros(10,51);
cost = zeros(10,51);

for num_samples = 1:10
disp(["Starting : ", num2str(num_samples)])
tic
t = 0;
x_hat_t1t = 0.*[1;-1;1;-1];
P_t1t = sigma_x0;
x_t = mvnrnd(x_hat_t1t,sigma_x0,1).';
while t<50
    y_t = zeros(N*p,1);
    y_hat_t = zeros(N*p,1);
    eta_t_opt = [];
    Theta_t_opt = [];
    for ii=0:N_c-1
        y_t(ii*p+1:(ii+1)*p) = C*x_t+mvnrnd(zeros(length(sigma_v),1),sigma_v,1).';
        x_hat_t1t1 = x_hat_t1t+P_t1t*C.'*((C*P_t1t*C.'+sigma_v)\(y_t(ii*p+1:(ii+1)*p)-C*x_hat_t1t));
        P_t1t1 = P_t1t-P_t1t*C.'*((C*P_t1t*C.'+sigma_v)\C)*P_t1t;
        x_norm(num_samples,t+ii+1) = norm(x_hat_t1t1);
        if ii==0
            x_hat_o = x_hat_t1t1(n_s+1:end);
            Lambda_x_phi = Lambda_e_phi+x_hat_t1t1*Lambda_phi.';
            [eta_t_opt,Theta_t_opt,sol] = alg2_optimize(N,N_c,A_cal,B_cal,D_cal, ...
    x_hat_t1t1,Q_cal,M_cal,Lambda_phi,Lambda_phiphi,Lambda_x_phi,Lambda_w_phi, ...
    phi_max,U_max,x_hat_o,A_o,R_AoBo_Nc,zeta,epsilon,n_o,m,p,S_cal, ...
    S_tilde_cal,L_cal,alpha,beta,P_t1t1,sigma_w_blk);
            if sol.problem~=0
                E_xt_xt = x_hat_t1t1*x_hat_t1t1.'+P_t1t1;
                alpha_star = 3*trace(A_cal.'*S_cal*A_cal*E_xt_xt)+3*trace(D_cal.'*S_cal*D_cal*sigma_w_blk);
                alpha_star = alpha_star+L_cal.'*A_cal*x_hat_t1t1+norm(L_cal.'*B_cal,1)*U_max;
                alpha_star = alpha_star+3*N*m*(U_max^2)*max(svd(B_cal.'*S_cal*B_cal));
                beta_star = N*m*(U_max^2)*max(svd(S_tilde_cal));
                alpha_t = alpha;
                beta_t = beta;
                alpha_up = alpha_star;
                alpha_down = alpha_t;
                beta_up = beta_star;
                beta_down = beta_t;
                [eta_t_opt,Theta_t_opt,sol] = alg2_optimize(N,N_c,A_cal,B_cal,D_cal, ...
    x_hat_t1t1,Q_cal,M_cal,Lambda_phi,Lambda_phiphi,Lambda_x_phi,Lambda_w_phi, ...
    phi_max,U_max,x_hat_o,A_o,R_AoBo_Nc,zeta,epsilon,n_o,m,p,S_cal, ...
    S_tilde_cal,L_cal,alpha_up,beta_up,P_t1t1,sigma_w_blk);
                nu = 1;
                while ~((abs(alpha_up-alpha_down)<=1e-1&&abs(beta_up-beta_down)<=1e-1)||nu>8)
                    alpha_t = (alpha_up+alpha_down)/2;
                    beta_t = (beta_up+beta_down)/2;
                    [eta1,Theta1,sol] = alg2_optimize(N,N_c,A_cal,B_cal,D_cal, ...
    x_hat_t1t1,Q_cal,M_cal,Lambda_phi,Lambda_phiphi,Lambda_x_phi,Lambda_w_phi, ...
    phi_max,U_max,x_hat_o,A_o,R_AoBo_Nc,zeta,epsilon,n_o,m,p,S_cal, ...
    S_tilde_cal,L_cal,alpha_t,beta_t,P_t1t1,sigma_w_blk);
                    if sol.problem==0
                        alpha_up = alpha_t;
                        beta_up = beta_t;
                        eta_t_opt = eta1;
                        Theta_t_opt = Theta1;
                    else
                        alpha_down = alpha_t;
                        beta_down = beta_t;
                    end
                    nu = nu+1;
                end
            end
        end
        y_hat_t(ii*p+1:(ii+1)*p) = C*x_hat_t1t1;
        u_t1_vec = eta_t_opt+Theta_t_opt*phi_func(y_t-y_hat_t,phi_max);
        u_t1 = u_t1_vec(ii*m+1:(ii+1)*m);
        u(num_samples,t+ii+1) = u_t1;
        cost(num_samples,t+ii+1) = x_hat_t1t1.'*Q*x_hat_t1t1+u_t1.'*R*u_t1;
        x_t = A*x_t+B*u_t1+mvnrnd(zeros(length(sigma_w),1),sigma_w,1).';
        x_hat_t1t = A*x_hat_t1t1+B*u_t1;
        P_t1t = A*P_t1t1*A.'+sigma_w;
        disp(['t = ',num2str(t+ii),'s Done']);
    end
    t = t+N_c;
end
elapsed_time = toc;  % Stop timing
disp("Done")
fprintf('Sample %d completed in %.2f seconds.\n', num_samples, elapsed_time);
end

%% Plotting
s1 = load("alg2_data_1300.mat","x_norm","u","cost","alpha");
s2 = load("alg2_data_1300_50_samples.mat","x_norm","u","cost");
alpha = s1.alpha;
x_norm_clipped = [s1.x_norm;s2.x_norm];
u_clipped = [s1.u;s2.u];
cost_clipped = [s1.cost;s2.cost];
avg_state_norm = mean(x_norm_clipped, 1);
time_axis = 0:(size(x_norm_clipped, 2)-1);
std_state_norm = std(x_norm_clipped, 0, 1);
cost1 = cost_clipped;

for ii=1:50
    %for jj=1:ii
        cost1(:,ii+1) = cost1(:,ii+1)+cost1(:,ii);
    %end
end
avg_cost = mean(cost1,1);
X_norm = zeros(1,length(0:3:45));
U = zeros(1,length(0:3:45));
for ii=0:3:45
    X_norm((ii/3)+1) = mean(sum(x_norm_clipped(:,ii+1:ii+6).^2,2));
    U((ii/3)+1) = mean(sum(u_clipped(:,ii+1:ii+6).^2,2));
end
avg_u = mean(u_clipped,1);

figure;
hold on;
plot(time_axis, avg_state_norm, 'r-', 'LineWidth', 2);
plot(time_axis, std_state_norm, 'r--', 'LineWidth', 1);

xlabel('Time step');
ylabel('State norm');
title('Average and Std Dev of State Norm over Time');
legend('Average', 'Std',Location='southeast');
grid on;
ylim([0,max(avg_state_norm)]);
xlim([0,51]);

figure;
plot(time_axis, avg_cost, 'r-', 'LineWidth', 2);

xlabel('Time step');
ylabel('Cost');
title('Average Cost over Time');
grid on;
ylim([0,max(avg_cost)]);
xlim([0,51]);

figure
plot(time_axis,avg_u,'r-','LineWidth',2);

xlabel('Time step');
ylabel('Input');
title('Average Control Input over Time');
grid on;
ylim([min(avg_u),max(avg_u)])
xlim([0,51]);

figure
plot(0:3:45,X_norm,'r-',LineWidth=2);
yline(alpha,'--');
xlabel('Time step');
ylabel('$E\left[\|X_t\|_{\mathcal{S}}^2+\mathcal{L}^\mathrm{T}X_t\right]$')
grid on
title('$E\left[\|X_t\|_{\mathcal{S}}^2+\mathcal{L}^\mathrm{T}X_t\right]$ over Time for $\alpha=1300$');
ylim([0,1400]);

figure
plot(0:3:45,U,'r-',LineWidth=2);
xlabel('Time step');
ylabel('$E\left[\|U_t\|_{\tilde{\mathcal{S}}}^2\right]$')
grid on
title('$E\left[\|U_t\|_{\tilde{\mathcal{S}}}^2\right]$ over Time for $\alpha=1300$');
%yline(beta);