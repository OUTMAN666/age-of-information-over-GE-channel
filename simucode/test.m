%This file is created by Guan Qixing on 10/11/2021
%It simlate the AoI performance of LGFC policy with instant feedback and GE
%channel

clc;
clear;
close all;

p_vec=0.1:0.05:0.5; % the transition probability in GE channel model
eta=1-2*p_vec;
lambda=1/3;% need to be less than 0.5 to ensure stable queue
k = 3;
iter = 100;
% limitN = 100;
BFavgAoI_vec=zeros(iter,length(p_vec));
BLavgAoI_vec=zeros(iter,length(p_vec));
pFavgAoI_vec=zeros(iter,length(p_vec));
pLavgAoI_vec=zeros(iter,length(p_vec));
% pFanaAoI = zeros(1,length(p_vec));
for i=1:length(p_vec)
    p=p_vec(i)
    for j = 1 : iter
         BFavgAoI_vec(j,i)=getAoI_FCFS_GE(p,p,lambda);
         BLavgAoI_vec(j,i)=getAoI_LGFS_GE(p,p,lambda);
         pFavgAoI_vec(j,i)=regularArrival_FCFS(p,p,k);
         pLavgAoI_vec(j,i)=regularArrival_LGFS(p,p,k);
    end
    % pFanaAoI(i) = anaReguFCFS(p, p, limitN, k);
end
BFAoIexpectation = mean(BFavgAoI_vec,1);
BLAoIexpectation = mean(BLavgAoI_vec,1);
pFAoIexpectation = mean(pFavgAoI_vec,1);
pLAoIexpectation = mean(pLavgAoI_vec,1);
BFanaAoI = 1/lambda + 1 + 4 * lambda^2 / (1 - 2 * lambda) + ...
    (1 - 2 * lambda) * eta ./ (1- (1 - 2 * lambda) * eta) + ...
    4*lambda^2/(1-2*lambda)*(1./(1-eta)./(1-(1-2*lambda)*eta).^2-1);
BLanaAoI = 1/lambda + 1/2./p_vec;
pLanaAoI = (k + 1) / 2 + 1/2./p_vec;
figure;
plot(eta,BFanaAoI,'r-');
hold on;
plot(eta,BFAoIexpectation,'rs','MarkerFaceColor','r');
plot(eta,BLanaAoI,'b-');
plot(eta,BLAoIexpectation,'bo','MarkerFaceColor','b');
% plot(eta,pFanaAoI,'--m');
plot(eta,pFAoIexpectation,'mV','MarkerFaceColor','m');
plot(eta,pLanaAoI,'--k');
plot(eta,pLAoIexpectation,'k^','MarkerFaceColor','k');
hold off;
grid on;
xlabel('Channel Memory');
ylabel('Average AoI');
% title('\lambda=0.4 symmetric GE channel');
legend('Ber-FCFS-ana(\lambda=1/3)', 'Ber-FCFS-simu(\lambda=1/3)', ...
    'Ber-LGFS-ana(\lambda=1/3)', 'Ber-LGFS-simu(\lambda=1/3)', ...
    'per-FCFS-simu(k=3)', ...
    'per-pLGFS-ana(k=3)', 'per-pLGFS-simu(k=3)');
save memoryAoI.mat