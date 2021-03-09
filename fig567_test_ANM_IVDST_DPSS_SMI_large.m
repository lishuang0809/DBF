% Copyright 2021. All Rights Reserved
% Code by Shuang Li
% For Paper, "Digital Beamforming Robust to Time-Varying Carrier Frequency Offset"
% arXiv:2103.04948, 2021.
% by Shuang Li, Payam Nayeri, and Michael B. Wakin
% This code is used to reproduce Figures 5-7.

%%
clc; clear; close all;
%air = randi(1e4)
air = 6641;  
rng(air);
% Generate data matrix X

% Array setup
N = 4; % number of elements
d = 0.5; % element separation in wavelengths
x = -(N-1)*d/2 : d : (N-1)*d/2; % element positions
% Array steering vector (theta in degrees):
asv = @(theta) exp(1i*2*pi*sind(theta)*x);

% Signal setup
ns = 300; % number of samples
ni = 2; % number of interferers
as = -20; % angle to the desired signal (degrees. 0=broadside)
ai = [-60 ; 20];  % angles to interferers
ia = [1 1]; % interference amplitude (1 -> same as intended transmitter)
is = 1;     % signal amplitude
t = (0:ns-1).'; % time samples

% === linear freq drift ===
% fs = -0.0002*t; %0.0005*rand(ns,1);  %fs = 0.0001*t; % intended signal frequency offset
% fi = t*[-0.0001 0.0002]; %rand(ns,1)*[0.0005 0.001];%fi = t*[0.0001 0.0002];%[1/100 1/200]; % interferer frequency offset
% 
% fs = fs - mean(fs);
% fi = fi - mean(fi);

% === random freq drift ===
fis0 = 0.001+0.0001*randn(ni+1,1);
fs = zeros(ns,1);    fs(1) = fis0(1); %0.01;
fi = zeros(ns,ni);   fi(1,:) = fis0(2:end);    %[0.01 0.01];
sigma = 1e-4;
for i = 2:ns
    fs(i) = fs(i-1) + sigma*randn(1,1);
    fi(i,:) = fi(i-1,:) + sigma*randn(1,2);
end

% === constant freq drift ===

% fs = 0.01*ones(ns,1);
% fi = ones(ns,1)*[0.02 0.03]; 

% === zigzag freq drift ===
% fs = [0.0003*(0:99)';0.0003*99-0.0003*(0:99)';0.0003*(0:99)']; %0.0005*rand(ns,1);  %fs = 0.0001*t; % intended signal frequency offset
% 
% fi = [-0.0002*(0:99)' -0.0004*(0:99)';-0.0002*99+0.0002*(0:99)' -0.0004*99+0.0004*(0:99)';-0.0002*(0:99)' -0.0004*(0:99)']; %0.0005*rand(ns,1);  %fs = 0.0001*t; % intended signal frequency offset
% 
% fs = fs - mean(fs);
% fi = fi - mean(fi);


maxdrift = max(max(abs([fs fi])))

fs0 = 0.2;  %1/5;
fi0 = [0.24 0.3];%[0.25 0.35];  %[1/10 1/20];

figure(1); clf;
lw=3;      % linewidth
Fs=25;     % FontSize
hold on;
h1 = plot(t,fs,'r-','linewidth',lw);
h2 = plot(t,fi(:,1),'b-','linewidth',lw);
h3 = plot(t,fi(:,2),'g-','linewidth',lw);
h = legend([h1 h2 h3],'$f_1$','$f_2$','$f_3$','Interpreter','LaTex','FontSize',Fs,'location','southwest');
xlabel('Samples','FontSize',Fs)
ylabel('Frequency offset','FontSize',Fs)
set(gca,'FontSize',Fs);
set(h,'FontSize',Fs,'numcolumn',3);
set(gcf, 'Color', 'w');
grid on;
box on;


noisAmp = 0;  %no noise  %0.01; % noise amplitude
pp = 500; % plot points

% Signals

s = exp(1j*2*pi*t*fs0)*is;% intended transmitter signal
si = exp(1j*2*pi*t*fi0)*diag(ia); % interferer signals

% Signal matrix generation
pds = exp(1i*2*pi*t.*fs); % intended transmitter phase drift
pdi = exp(1i*2*pi*t.*fi); % interferer phase drift
Nois = sqrt(noisAmp)*(randn(ns,N)+1i*randn(ns,N)); % noise matrix
X = (s.*pds)*asv(as) + (si.*pdi)*asv(ai) + Nois; % signal matrix

%% %%%%%%% IVDST+DPSS+SMI %%%%%%%%%%


% Creat DPSS basis
W = 1.1*max(max(abs([fs fi]))); 
L = 2 
Snw = dpss(ns,ns*W,L);

Err_pds = norm(pds-Snw*Snw'*pds)/norm(pds)   


% Given X, estimate dual solution with the fast alg
tic;
maxIter = 200; 
delta = 4;        % stepsize for GD
ratio = 0.99;     % decreasing ratio for stepsize

% Initialization
t0 = 1;
P0 = randn(ns,L,N)+1j*randn(ns,L,N); 
H0 = zeros(ns,ns,N);
for n = 1:N
    H0(:,:,n) = P0(:,:,n)*P0(:,:,n)';
end
Lstar_X = zeros(ns,L,N);
for n = 1:N
    Lstar_X(:,:,n) = Snw.*(X(:,n)*ones(1,L));
end

tt = zeros(1,maxIter+1); 
P_all = zeros(ns,L,N,maxIter);  P_all(:,:,:,1) = P0; 
H_all = zeros(ns,ns,N,maxIter); H_all(:,:,:,1) = H0; 
P_bar = zeros(ns,L,N,maxIter);
H_bar = zeros(ns,ns,N,maxIter);
P_g = zeros(ns,L,N,maxIter);
H_g = zeros(ns,ns,N,maxIter);
H_tilde = zeros(ns,ns,N,maxIter);
cost = zeros(1,maxIter);
for i = 1:maxIter
    if rem(i,50) == 0
        fprintf('Iter = %d\n',i)
        delta = delta*ratio;
    end
    
    % Smoothing
    if i == 1
        tt(i) = (1+sqrt(4*t0^2+1))/2;
        P_bar(:,:,:,i) = P_all(:,:,:,i);
        H_bar(:,:,:,i) = H_all(:,:,:,i);
    else
        tt(i) = (1+sqrt(4*tt(i-1)^2+1))/2;
        P_bar(:,:,:,i) = P_all(:,:,:,i) + (tt(i-1)-1)/tt(i)*(P_all(:,:,:,i)-P_all(:,:,:,i-1));
        H_bar(:,:,:,i) = H_all(:,:,:,i) + (tt(i-1)-1)/tt(i)*(H_all(:,:,:,i)-H_all(:,:,:,i-1));
    end
    
    % Gradient descent
    P_g(:,:,:,i) = P_bar(:,:,:,i) + delta*Lstar_X;
    H_g(:,:,:,i) = H_bar(:,:,:,i);
    
    % Proximal mapping
    for n = 1:N
        for j = 1:ns-1
            offdiagsum = sum(diag(H_g(:,:,n,i),j));
            H_tilde(:,:,n,i) = H_tilde(:,:,n,i) + diag(diag(H_g(:,:,n,i),j) - offdiagsum/(ns-j),j);
        end
        H_tilde(:,:,n,i) = H_tilde(:,:,n,i) + H_tilde(:,:,n,i)';         
    end
    diagsum = 0;
    for n = 1:N
        diagsum = diagsum + sum(diag(H_g(:,:,n,i)));        
    end
    for n = 1:N
        H_tilde(:,:,n,i) = H_tilde(:,:,n,i) + diag(diag(H_g(:,:,n,i))/diagsum);
    end
    for n = 1:N
        Zn = [H_tilde(:,:,n,i) -P_g(:,:,n,i); -P_g(:,:,n,i)' eye(L,L)];
        [Vn, Sn] = eig(Zn); 
        Zn_tilde = Vn*diag(max(diag(Sn),0))*Vn';
        H_all(:,:,n,i+1) = Zn_tilde(1:ns,1:ns);
        P_all(:,:,n,i+1) = -Zn_tilde(1:ns,ns+1:end);
    end
    % Compute cost function
    L_P = zeros(ns,N);
    for n = 1:N
        L_P(:,n) = (Snw.*P_all(:,:,n,i+1))*ones(L,1);
    end
    f1 = -real(trace(L_P'*X));
    f2 = 0; f3 = 0;
    for n = 1:N
        f2 = f2 + sum(diag(H_all(:,:,n,i+1)));
        for j = 1: ns-1
            f3 = f3 + abs(sum(diag(H_all(:,:,n,i+1),j)))^2;
        end
    end   
    f2 = abs(f2-1)^2;
    cost(i) = f1 + f2 + f3;
 
end


% Compute dual solution 
Y = -P_all(:,:,:,end)/2;
fgrid=0:0.01:1;              % frequency grid

% Compute dual polynomial
Phi_fY = zeros(L,N,length(fgrid));
norm_Phi_fY_fast = zeros(1,length(fgrid));
for ff = 1:length(fgrid)
    f = fgrid(ff);
    a=exp(1j*2*pi*t*f);       % atoms
    for n = 1:N
        Phi_fY(:,n,ff) = Y(:,:,n)'*a;        
    end
    norm_Phi_fY_fast(ff) = norm(squeeze(Phi_fY(:,:,ff)),'fro');
end

norm_Phi_fY_fast = norm_Phi_fY_fast/max(norm_Phi_fY_fast);

eps = 0.7*max(norm_Phi_fY_fast);
indx = find(norm_Phi_fY_fast >= eps);
f_pts = fgrid(indx);
[idx, c] = kmeans(f_pts(:), 3);
f_c = sort(c)'

Af = exp(1j*2*pi*t*f_c);
alpha1 = pinv(Af)*squeeze(Lstar_X(:,:,1));
alpha1 = alpha1(1,:)/norm(alpha1(1,:));
alpha1 = alpha1.';

s_new = Af(:,1);
w_fastANM_DPSS = pinv(X)*(s_new.*(Snw*alpha1));
w_fastANM_DPSS = w_fastANM_DPSS/(asv(as)*w_fastANM_DPSS);

runtime_fastANM = toc

%% %%%%%%% SMI %%%%%%%%%%

w_SMI = pinv(X)*s;
w_SMI = w_SMI./norm(w_SMI);

%% compute the radiation pattern

thList = linspace(-180,180,pp); % angle list
AF_SMI = 20*log10(abs(asv(thList.')*w_SMI));
AF_fastANM_DPSS_SMI = 20*log10(abs(asv(thList.')*w_fastANM_DPSS)); % array factor/radiation pattern

%% plot dual polynomial and the radiation pattern
figure(2)
clf;
lw=3;      % linewidth
Fs=25;     % FontSize
hold on;
plot(fgrid,norm_Phi_fY_fast,'-','Linewidth',lw-1,'color',[0.4 0 1]);
xlabel('Frequency','FontSize',Fs);
ylabel('Dual polynomial','FontSize',Fs)
h = legend('IVDST+DPSS+SMI');
set(h,'FontSize',Fs-8,'location','northeast');
set(gca,'FontSize',Fs);
grid on;
box on;
set(gcf, 'Color', 'w'); 
axis([0 1 0 1.05])



figure(3)
clf;
lw=3;      % linewidth
Fs=25;     % FontSize
h1 = plot(thList,AF_SMI,'k:','Linewidth',lw-1);
hold on;
h4 = plot(thList,AF_fastANM_DPSS_SMI,'-','Linewidth',lw-1,'color',[0.4 0 1]);
hold on;
plot([as as],[-60 15],'b','Linewidth',lw+1.5); % direction to intended transmitter
for thi = ai.' % directions to interferers
    plot([thi thi],[-60 15],'r','Linewidth',lw+1.5);
end
hold off
axis([-90 90 -40 10])
xlabel('Direction (deg)','FontSize',Fs)
ylabel('Radiation pattern (dB)','FontSize',Fs)
h = legend([h1 h4],'SMI','IVDST+DPSS+SMI');
set(h,'FontSize',Fs-8,'location','northeast','numcolumn',2);
set(gca,'FontSize',Fs);
xticks([-90,-60,-20,0,20,90])
yticks([-40,-30,-20,-10,0,10])
set(gcf, 'Color', 'w');    
grid on;


