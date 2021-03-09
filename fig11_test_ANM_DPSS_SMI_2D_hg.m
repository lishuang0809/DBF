% Copyright 2021. All Rights Reserved
% Code by Shuang Li
% For Paper, "Digital Beamforming Robust to Time-Varying Carrier Frequency Offset"
% arXiv:2103.04948, 2021.
% by Shuang Li, Payam Nayeri, and Michael B. Wakin
% This code is used to reproduce Figure 11.

%% Generate data matrix X
clc; clear; close all;

% Array setup
N = 4; % number of elements
d = 0.5; % element separation in wavelengths
x = -(N-1)*d/2 : d : (N-1)*d/2; % element positions
% Array steering vector (theta in degrees):
asv = @(theta) exp(1i*2*pi*sind(theta)*x);

% Signal setup
ns = 15; % number of samples % ns=30 for zigzag ns=15 for others
ni = 2; % number of interferers
as = -20; % angle to the desired signal (degrees. 0=broadside)
ai = [-60 ; 20];  % angles to interferers
ia = [1 1]; % interference amplitude (1 -> same as intended transmitter)
is = 1;     % signal amplitude
t = (0:ns-1).'; % time samples

trial = 20;
AF_SMIt = zeros(ni,trial);
AF_2Dt = zeros(ni,trial);
for tr = 1:trial
    fprintf('trial = %d\n',tr)
    % === linear freq drift ===
    
%     slp = 6;
%     slp1 = 0.001*randi(slp);  slp2 = -0.001*randi(slp); slp3 = -0.001*randi(slp);
%     fs = slp1*t; %0.0005*rand(ns,1);  %fs = 0.0001*t; % intended signal frequency offset
%     fi = t*[slp2 slp3]; %rand(ns,1)*[0.0005 0.001];%fi = t*[0.0001 0.0002];%[1/100 1/200]; % interferer frequency offset
%     
%     fs = fs - mean(fs);
%     fi = fi - mean(fi);
    
    % === random freq drift ===
%     fis0 = 0.03+0.001*randn(ni+1,1);
%     fs = zeros(ns,1);    fs(1) = fis0(1); %0.01;
%     fi = zeros(ns,ni);   fi(1,:) = fis0(2:end);    %[0.01 0.01];
%     sigma = 2e-3;
%     for i = 2:ns
%         fs(i) = fs(i-1) + sigma*randn(1,1);
%         fi(i,:) = fi(i-1,:) + sigma*randn(1,2);
%     end
    
    % === constant freq drift ===
    
    slp = 6;
    slp1 = 0.01*randi(slp);  slp2 = -0.01*randi(slp); slp3 = 0.01*randi(slp);
    fs = slp1*ones(ns,1);
    fi = ones(ns,1)*[slp2 slp3];
    
    
    % === zigzag freq drift ===
%     slp = 6;
%     slp1 = -0.001*randi(slp);  slp2 = -0.001*randi(slp); slp3 = 0.001*randi(slp);
%     fs = [slp1*(0:9)';slp1*9-slp1*(0:9)';slp1*(0:9)']; %0.0005*rand(ns,1);  %fs = 0.0001*t; % intended signal frequency offset
%     fi = [slp2*(0:9)' slp3*(0:9)';slp2*9-slp2*(0:9)' slp3*9-slp3*(0:9)';slp2*(0:9)' slp3*(0:9)']; %0.0005*rand(ns,1);  %fs = 0.0001*t; % intended signal frequency offset
%     
%     fs = fs - mean(fs);
%     fi = fi - mean(fi);
    
    
    maxdrift = max(max(abs([fs fi])))
    
    fs0 = 0.2;  %1/5;
    fi0 = [0.7 0.7];%[0.25 0.35];  %[1/10 1/20];
    
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
    
    % A more general way to generate X (use DPSS)
    W = 1.1*max(max(abs([fs fi]))); %0.005*max(max(abs([fs fi])));
    L = 4; %round(2*ns*W)-47
    Snw = dpss(ns,ns*W,L);
     
    Err_freqdrift = norm([pds pdi]-Snw*Snw'*[pds pdi],'fro')/norm([pds pdi],'fro')   % should be close to 0
    
    %% %%%%%%% ANM+DPSS+SMI 2D %%%%%%%%%%
    
    % Given X, estimate Xv(thus f^o) with ANM+DPSS when si is assumed to contain complex exponentials
    tic;
    tau = 1e-2*norm(X,'fro');
    cvx_begin   quiet 
    cvx_precision low
    
    variable T(2*ns-1,2*L-1,N) complex;
    variable ST(ns*L,ns*L,N) hermitian;
    variable Xv(ns,L,N) complex;
    variable tt(1,N)
    dual variable Y0{N}
    obj = 0;
    for n = 1:N
        obj = obj + trace(squeeze(ST(:,:,n)));
    end
    
    minimize(0.5*sum(tt)/N+0.5*obj/ns/N)
    subject to
    for n = 1:N
        for i = 1:2*ns-1
            for r = 1:ns
                for c = 1:ns
                    if c-r == ns-i
                        ST((r-1)*L+1:r*L,(c-1)*L+1:c*L,n) == ...
                            toeplitz(transpose(squeeze(T(i,L:end,n))),fliplr(squeeze(T(i,1:L,n))));
                    end
                end
            end
        end
        Y0{n} : [squeeze(ST(:,:,n)), vec(squeeze(Xv(:,:,n)).');
            vec(squeeze(Xv(:,:,n)).')',     tt(n)] == semidefinite(ns*L+1,ns*L+1);
        norm((Snw.*Xv(:,:,n))*ones(L,1)-X(:,n)) <= tau;
    end
    cvx_end
    
    
    %% Compute dual solution
    Y = zeros(ns,L,N);
    for n = 1:N
        Y(:,:,n) = -2*reshape(Y0{n}(1:end-1,end),L,ns).';
    end
    
    % Compute dual polynomial
    fgrid = 0:0.01:1;     % frequency grid
    tgrid = -90:1:90;     % angle grid
    lenf = length(fgrid);
    lent = length(tgrid);
    QQt0 = zeros(lenf,lent);
    for i = 1:lenf
        for j = 1:lent
            temp0 = zeros(L,1);
            for n = 1:N
                a=exp(1j*2*pi*t*fgrid(i));
                temp0 = temp0 + Y(:,:,n)'*a*exp(1j*2*pi*sind(tgrid(j))*x(n));
            end
            QQt0(i,j) = norm(temp0);
        end
    end
    
    
    %%
    eps = 0.7*max(max(QQt0));
    QQt0T = QQt0.';
    [pks,locs] = findpeaks(QQt0T(:), 'MinPeakHeight',eps, 'MinPeakProminence',0.5);
    [r,c] = ind2sub(size(QQt0T), locs);
    figure(4); clf;
    mesh(tgrid, fgrid, QQt0)
    hold on
    plot3(tgrid(r),fgrid(c), pks, 'r*')
    hold off
    grid on
    
    
    f_pts = fgrid(c);
    t_pts = tgrid(r);
    
    [t_idx, c_t] = kmeans(t_pts(:), 3);
    t_c = sort(c_t)'
    
    c_f = zeros(size(t_c));
    for i = 1:3
        idx_cf = find(t_idx==i);
        c_f(i) = mean(f_pts(idx_cf));
    end
    c_f = c_f';
    
    
    f_c = sort(c_f)'
    
    %%
    
    Af = exp(1j*2*pi*t*f_c);
    alpha1 = pinv(Af)*squeeze(Xv(:,:,1));
    alpha1 = alpha1(1,:)/norm(alpha1(1,:));
    alpha1 = alpha1.';
    
    s_new = Af(:,1);
    w_ANM_DPSS = pinv(X)*(s_new.*(Snw*alpha1));
    w_ANM_DPSS = w_ANM_DPSS/(asv(as)*w_ANM_DPSS);
    runtime_SDP_DPSS = toc
    
    %% %%%%%%% SMI %%%%%%%%%%
    
    w_SMI = pinv(X)*s;
    w_SMI = w_SMI/(asv(as)*w_SMI);
    
    %% compute the radiation pattern
    
    thList = linspace(-180,180,pp); % angle list
    AF_SMI = 20*log10(abs(asv(thList.')*w_SMI));
    AF_ANM_DPSS_SMI_2D = 20*log10(abs(asv(thList.')*w_ANM_DPSS)); % array factor/radiation pattern
    
    AF_SMIt(:,tr) = 20*log10(abs(asv(ai)*w_SMI));
    AF_2Dt(:,tr) = 20*log10(abs(asv(ai)*w_ANM_DPSS));
    %% plot dual polynomial and the radiation pattern
    figure(2)
    clf;
    lw=3;      % linewidth
    Fs=25;     % FontSize
    ms = 15;
    hold on;
    imagesc(tgrid,fgrid,QQt0);    
    hold on;
    plot([as ai'],[fs0 fi0], 'r*','markersize',ms,'linewidth',lw);
    ylabel('Frequency','FontSize',Fs);
    xlabel('Direction (deg)','FontSize',Fs)
    xticks([-90,-60,-20,0,20,90])
    set(gca,'FontSize',Fs);
    set(gcf, 'Color', 'w');
    axis tight;
    
    %%
    
    figure(3)
    clf;
    lw=3;      % linewidth
    Fs=25;     % FontSize
    h1 = plot(thList,AF_SMI,'k:','Linewidth',lw-1);
    hold on;
    h3 = plot(thList,AF_ANM_DPSS_SMI_2D,'-','Linewidth',lw-1,'color',[0.5 0.8 0.9]);
    hold on;
    plot([as as],[-60 15],'b','Linewidth',lw+1.5); % direction to intended transmitter
    for thi = ai.' % directions to interferers
        plot([thi thi],[-60 15],'r','Linewidth',lw+1.5);
    end
    hold off
    axis([-90 90 -40 10])
    xlabel('Direction (deg)','FontSize',Fs)
    ylabel('Radiation pattern (dB)','FontSize',Fs)
    h = legend([h1 h3],'SMI','2D-ANM+DPSS+SMI');
    set(h,'FontSize',Fs-8,'location','northeast','numcolumn',2);
    set(gca,'FontSize',Fs);
    xticks([-90,-60,-20,0,20,90])
    yticks([-40,-30,-20,-10,0,10])
    set(gcf, 'Color', 'w');
    grid on;
        
    
end

%%
figure(4)
clf;
lw=3;      % linewidth
Fs=25;     % FontSize
subplot(2,1,1)
histogram(AF_SMIt(1,:)); hold on;
histogram(AF_2Dt(1,:));
set(gca,'FontSize',Fs);
h = legend('SMI','2D-ANM+DPSS+SMI')
set(h,'FontSize',Fs-5);
ylim([0,trial+2])

title('Radiation pattern (dB) at -60','FontSize',Fs-5)
subplot(2,1,2)
histogram(AF_SMIt(2,:)); hold on;
histogram(AF_2Dt(2,:));
set(gca,'FontSize',Fs);
set(gcf, 'Color', 'w');
title('Radiation pattern (dB) at 20','FontSize',Fs-5)
ylim([0,trial])







