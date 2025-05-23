function [lambda_phi, lambda_phi_phi, lambda_e_phi, lambda_w_phi] = monte_carlo_simulation(A,C,Sigma_v,Sigma_w,Sigma_xo, phi_max,N)
[P, ~, ~] = idare(A.', C.', Sigma_w, Sigma_v, [], []);
%P_o = P - P*C.'*inv(C*P*C.'+Sigma_v)*C*P;
[n,p] = size(C);

lambda_phi = zeros(p*N,1);
lambda_w_phi = zeros(n*N, p*N);
lambda_e_phi = zeros(n, p*N);
lambda_phi_phi = zeros(p*N, p*N);

num_samples = 1e5;

disp("Starting...")
for i = 1:num_samples
    x = mvnrnd(zeros(n,1),Sigma_xo).';
    x_hat = x;
    phi_stack = zeros(p * N, 1);
    w_stack = zeros(n*N,1);
    for t = 1:100
        w = mvnrnd(zeros(n,1),Sigma_w).';
        x = A*x + w;
        y = C*x + mvnrnd(zeros(p,1),Sigma_v).';
        x_hat = A*x_hat;
        x_hat = x_hat + P*(C.'/(C*P*C.'+Sigma_v))*(y-C*x_hat);
        y_hat = C*x_hat;
       
        diff_phiy = max(min(y - y_hat, phi_max),-phi_max);

        l = mod(t-1,N) + 1;
        phi_stack((l-1)*p+1 : l*p) = diff_phiy;
        w_stack((l-1)*n+1 : l*n) = w;

        if l == 1 
            diffx = x - x_hat;
        end
    end
    lambda_phi = lambda_phi + phi_stack;
    lambda_w_phi = lambda_w_phi + w_stack * phi_stack.';
    lambda_phi_phi = lambda_phi_phi + phi_stack * phi_stack.';
    lambda_e_phi   = lambda_e_phi + diffx * phi_stack.';
end

disp("Done")
lambda_phi = lambda_phi/num_samples;
lambda_w_phi = lambda_w_phi/num_samples;
lambda_phi_phi = lambda_phi_phi/num_samples;
lambda_e_phi = lambda_e_phi/num_samples;