function [eta_t_opt,Theta_t_opt,sol] = alg2_optimize(N,N_c,A_cal,B_cal,D_cal, ...
    x_hat_t,Q_cal,M_cal,Lambda_phi,Lambda_phiphi,Lambda_x_phi,Lambda_w_phi, ...
    phi_max,U_max,x_hat_o,A_o,R_AoBo_Nc,zeta,epsilon,n_o,m,p,S_cal, ...
    S_tilde_cal,L_cal,alpha,beta,P_tt,sigma_w)
%ALG2_OPTIMIZE Summary of this function goes here
%   Detailed explanation goes here
z_1 = sdpvar(1,1);
z_2 = sdpvar(1,1);
z_3 = sdpvar(1,1);
eta_t = sdpvar(m*N,1);
Theta_t = sdpvar(N*m,N*p,'full');
Objective = z_1;
Constraints = [];
c1 = (eta_t+Theta_t*Lambda_phi).'*M_cal*(eta_t+Theta_t*Lambda_phi);
c1 = c1+trace(Theta_t.'*M_cal*Theta_t*(Lambda_phiphi-Lambda_phi*Lambda_phi.'));
c1 = c1+2*x_hat_t.'*A_cal.'*Q_cal*B_cal*eta_t;
c1 = c1+2*trace(Theta_t.'*B_cal.'*Q_cal*(D_cal*Lambda_w_phi+A_cal*Lambda_x_phi));
Constraints = [Constraints,c1<=z_1];
for jj=1:N*m
    Constraints = [Constraints,abs(eta_t(jj))+norm(Theta_t(jj,:),1)*phi_max<=U_max];
end
if norm(x_hat_o)>=zeta+epsilon
    Constraints = [Constraints,norm(A_o^N_c*x_hat_o+R_AoBo_Nc*eta_t(1:N_c*m))<=z_2];
    Constraints = [Constraints,norm(R_AoBo_Nc*Theta_t(1:N_c*m,:),inf)<=z_3];
    Constraints = [Constraints,z_2+sqrt(n_o)*phi_max*z_3<=norm(x_hat_o)-zeta-epsilon/2];
end
for kk = 1:N
    for ll = 1:N
        if ll>kk
            Constraints = [Constraints,Theta_t((kk-1)*m+1:kk*m,(ll-1)*p+1:ll*p)==0];
        end
    end
end
E_xt_xt = x_hat_t*x_hat_t.'+P_tt;
c2 = (eta_t+Theta_t*Lambda_phi).'*(B_cal.'*S_cal*B_cal)*(eta_t+Theta_t*Lambda_phi);
c2 = c2+trace(Theta_t.'*B_cal.'*S_cal*B_cal*Theta_t*(Lambda_phiphi-Lambda_phi*Lambda_phi.'));
c2 = c2+2*x_hat_t.'*A_cal.'*S_cal*B_cal*eta_t+2*trace(Theta_t.'*B_cal.'*S_cal*(D_cal*Lambda_w_phi+A_cal*Lambda_x_phi));
c2 = c2+L_cal.'*B_cal*(eta_t+Theta_t*Lambda_phi)+trace(A_cal.'*S_cal*A_cal*E_xt_xt);
c2 = c2+trace(D_cal.'*S_cal*D_cal*sigma_w)+L_cal.'*A_cal*x_hat_t;
Constraints = [Constraints,c2<=alpha];
c3 = (eta_t+Theta_t*Lambda_phi).'*S_tilde_cal*(eta_t+Theta_t*Lambda_phi);
c3 = c3+trace(Theta_t.'*S_tilde_cal*Theta_t*(Lambda_phiphi-Lambda_phi*Lambda_phi.'));
Constraints = [Constraints,c3<=beta];
options = sdpsettings('solver','sdpt3','verbose',0);
sol = optimize(Constraints,Objective,options);
if sol.problem == 0
    %disp('Solution Found');
else
    %disp(['Problem status: ', sol.info]);
    yalmiperror(sol.problem);
end
eta_t_opt = value(eta_t);
Theta_t_opt = value(Theta_t);
end

