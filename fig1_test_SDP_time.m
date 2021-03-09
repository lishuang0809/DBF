% Copyright 2021. All Rights Reserved
% Code by Shuang Li
% For Paper, "Digital Beamforming Robust to Time-Varying Carrier Frequency Offset"
% arXiv:2103.04948, 2021.
% by Shuang Li, Payam Nayeri, and Michael B. Wakin
% This code is used to reproduce Figure 1.

%%
clc; clear; close all;

NN = 500:500:5000;
Trial = 1;
for t = 1:Trial
    for n = 1:length(NN)
        fprintf('Trial = %d, N = %d\n',t,NN(n))
        N = NN(n);
        
        %% dense matrix
        A = sqrt(N)*ones(N,N); A = A - diag(diag(A));
        tic;
        cvx_begin quiet
        variable x;
        minimize x
        subject to
        X = x*eye(N) + A;
        X == semidefinite(N,N);
        cvx_end
        time_dense(t,n) = toc;
        
        %% block diagonal matrix
        nb = 10; % number of diagonal blocks 
        Asp = sqrt(N/nb)*ones(N/nb,N/nb); As = zeros(N,N);
        for j = 1:nb
            As((j-1)*N/nb+1:j*N/nb,(j-1)*N/nb+1:j*N/nb) = Asp;
        end
        As = As - diag(diag(As));
        tic;
        cvx_begin quiet
        variable x;
        minimize x
        subject to
        X = x*eye(N) + As;
        X == semidefinite(N,N);
        cvx_end
        time_sparse1(t,n) = toc;
        
        nb = 20; % number of diagonal blocks 
        Asp = sqrt(N/nb)*ones(N/nb,N/nb); As = zeros(N,N);
        for j = 1:nb
            As((j-1)*N/nb+1:j*N/nb,(j-1)*N/nb+1:j*N/nb) = Asp;
        end
        As = As - diag(diag(As));
        tic;
        cvx_begin quiet
        variable x;
        minimize x
        subject to
        X = x*eye(N) + As;
        X == semidefinite(N,N);
        cvx_end
        time_sparse2(t,n) = toc;
        
        nb = 50; % number of diagonal blocks 
        Asp = sqrt(N/nb)*ones(N/nb,N/nb); As = zeros(N,N);
        for j = 1:nb
            As((j-1)*N/nb+1:j*N/nb,(j-1)*N/nb+1:j*N/nb) = Asp;
        end
        As = As - diag(diag(As));
        tic;
        cvx_begin quiet
        variable x;
        minimize x
        subject to
        X = x*eye(N) + As;
        X == semidefinite(N,N);
        cvx_end
        time_sparse3(t,n) = toc;
        
        nb = 100; % number of diagonal blocks 
        Asp = sqrt(N/nb)*ones(N/nb,N/nb); As = zeros(N,N);
        for j = 1:nb
            As((j-1)*N/nb+1:j*N/nb,(j-1)*N/nb+1:j*N/nb) = Asp;
        end
        As = As - diag(diag(As));
        tic;
        cvx_begin quiet
        variable x;
        minimize x
        subject to
        X = x*eye(N) + As;
        X == semidefinite(N,N);
        cvx_end
        time_sparse4(t,n) = toc;
        
    end    
end
time_dense_m = time_dense;
time_sparse1_m = time_sparse1;
time_sparse2_m = time_sparse2;
time_sparse3_m = time_sparse3;
time_sparse4_m = time_sparse4;

%%
figure(1); clf;
lw = 3;
fs = 20;
ms = 15;
h1 = plot(NN,time_dense_m,'b-o','LineWidth',lw,'markersize',ms); hold on;
h2 = plot(NN,time_sparse1_m,'r-*','LineWidth',lw,'markersize',ms); hold on;
h3 = plot(NN,time_sparse2_m,'r-+','LineWidth',lw,'markersize',ms); hold on;
h4 = plot(NN,time_sparse3_m,'r-x','LineWidth',lw,'markersize',ms); hold on;
h5 = plot(NN,time_sparse4_m,'r-d','LineWidth',lw,'markersize',ms); hold on;
h6 = plot(NN,NN.^2.5*2.5e-7,'k:','LineWidth',lw); hold on;
xlabel('Size of PSD constraint matrix ($N_s$)','Interpreter','LaTex','FontSize',fs);
ylabel('Running time','FontSize',fs);
h = legend([h1 h2 h3 h4 h5 h6],'Dense','Sparse-10','Sparse-20','Sparse-50','Sparse-100','$2.5\times 10^{-7}N_s^{2.5}$','Interpreter','LaTex');
set(gca,'FontSize',fs);
set(h,'FontSize',fs,'location','northwest');
set(gcf, 'Color', 'w');
grid on;
axis tight;

