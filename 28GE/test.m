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
% pFanaAoI = zeros(1,length(p_vec));
willavgAoI_vec = zeros(iter, length(p_vec));
willavgAoIp_vec = zeros(iter, length(p_vec));
for i=1:length(p_vec)
    p=p_vec(i)
    for j = 1 : iter
         BFavgAoI_vec(j,i)=getAoI_FCFS_GE(p,p,lambda);
         BLavgAoI_vec(j,i)=getAoI_LGFS_GE(p,p,lambda);
         pFavgAoI_vec(j,i)=regularArrival_FCFS(p,p,k);
         pLavgAoI_vec(j,i)=regularArrival_LGFS(p,p,k);
         willavgAoI_vec(j,i) = getAoI_will_GE(p,p);
         willavgAoIp_vec(j,i) = getAoI_will_GE_pLGFS(p,p);
    end
%    pFanaAoI(i) = anaReguFCFS(p, p, 3);
end
BFAoIexpectation = mean(BFavgAoI_vec,1);
BLAoIexpectation = mean(BLavgAoI_vec,1);
pFAoIexpectation = mean(pFavgAoI_vec,1);
pLAoIexpectation = mean(pLavgAoI_vec,1);
willAoIexpectation = mean(willavgAoI_vec,1);
willpAoIexpectation = mean(willavgAoIp_vec,1);
figure;
% plot(eta,BFanaAoI,'r-');
% plot(eta,willanaAoI,'--k');
plot(eta,BFAoIexpectation,'ro-','MarkerFaceColor','r');

hold on;
% plot(eta,BLanaAoI,'b-');
% plot(eta,pFanaAoI,'--m');
% plot(eta,pLanaAoI,'--k');
plot(eta,pFAoIexpectation,'ko-','MarkerFaceColor','k');
plot(eta,willAoIexpectation,'go-','MarkerFaceColor','g');
plot(eta,BLAoIexpectation,'--r>','LineWidth',1.5);
plot(eta,pLAoIexpectation,'--k>','LineWidth',1.5);
plot(eta,willpAoIexpectation,'--g>','LineWidth',1.5);
hold off;
grid on;
xlabel('Channel Memory');
ylabel('Average AoI');
lgd = legend('Ber-FCFS(\lambda=1/3)','per-FCFS(K=3)','will-FCFS', ...
    'Ber-pLGFS(\lambda=1/3)','per-pLGFS(K=3)','will-pLGFS');
lgd.NumColumns = 1;
save memoryAoI.mat