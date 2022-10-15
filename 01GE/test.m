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
iter = 1000;
BFavgAoI_vec=zeros(iter,length(p_vec));
BLavgAoI_vec=zeros(iter,length(p_vec));
pFavgAoI_vec=zeros(iter,length(p_vec));
pLavgAoI_vec=zeros(iter,length(p_vec));
pFanaAoI = zeros(1,length(p_vec));
willavgAoI_vec = zeros(iter, length(p_vec));
willavgAoIpLGFS_vec = zeros(iter, length(p_vec));
for i=1:length(p_vec)
    p=p_vec(i)
    for j = 1 : iter
         BFavgAoI_vec(j,i)=getAoI_FCFS_GE(p,p,lambda);
         BLavgAoI_vec(j,i)=getAoI_LGFS_GE(p,p,lambda);
         pFavgAoI_vec(j,i)=regularArrival_FCFS(p,p,k);
         pLavgAoI_vec(j,i)=regularArrival_LGFS(p,p,k);
         willavgAoI_vec(j,i) = getAoI_will_GE(p,p);
         willavgAoIpLGFS_vec(j,i) = getAoI_will_GE_pLGFS(p,p);
    end
    pFanaAoI(i) = anaReguFCFS(p, p, 3);
end
BFAoIexpectation = mean(BFavgAoI_vec,1);
BLAoIexpectation = mean(BLavgAoI_vec,1);
pFAoIexpectation = mean(pFavgAoI_vec,1);
pLAoIexpectation = mean(pLavgAoI_vec,1);
willAoIexpectation = mean(willavgAoI_vec,1);
willpAoIexpectation = mean(willavgAoIpLGFS_vec,1);
BFanaAoI = 1/lambda + 1./(1-eta) + 4*lambda^2/(1-2*lambda)./(1-eta)./(1-eta+2*lambda*eta);
willanaAoI = 2 + 1/2./p_vec;
willpLGFSanaAoI = 1 + 1/2./p_vec;
BLanaAoI = 1/lambda + 1/2./p_vec;
pLanaAoI = (k + 1) / 2 + 1/2./p_vec;
figure;
plot(eta,BFanaAoI,'r-');
hold on;
plot(eta,pFanaAoI,'k-');
plot(eta,willanaAoI,'g-');

plot(eta,BLanaAoI,'--r');
plot(eta,pLanaAoI,'--k');
plot(eta,willpLGFSanaAoI,'--g');

plot(eta,BFAoIexpectation,'ro','MarkerFaceColor','r');
plot(eta,pFAoIexpectation,'ko','MarkerFaceColor','k');
plot(eta,willAoIexpectation,'go','MarkerFaceColor','g');

plot(eta,BLAoIexpectation,'r>');

plot(eta,pLAoIexpectation,'k>');

plot(eta,willpAoIexpectation,'g>');
hold off;
grid on;
xlabel('Channel Memory');
ylabel('Average AoI');
lgd = legend('Ber-FCFS-ana(\lambda=1/3)', 'per-FCFS-ana(K=3)', ...
    'will-FCFS-ana','Ber-pLGFS-ana(\lambda=1/3)', ...
     'per-pLGFS-ana(K=3)','will-pLGFS-ana', ...
    'Ber-FCFS-simu(\lambda=1/3)','per-FCFS-simu(K=3)', ...
    'will-FCFS-simu','Ber-pLGFS-simu(\lambda=1/3)', ...
     'per-pLGFS-simu(K=3)','will-pLGFS-simu');
lgd.NumColumns = 2;
save memoryAoI.mat
