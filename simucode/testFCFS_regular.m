%This file is created by Guan Qixing on 04/04/2022
%It simlate the AoI performance of FCFC policy with GE channel

clc;
clear;
close all;

p_vec=0.1:0.05:0.5; % the transition probability in GE channel model
eta=1-2*p_vec;
k = 5;% need to be less than 0.5 to ensure stable queue

iter = 50;
avgAoI_vec=zeros(iter,length(p_vec));
for i=1:length(p_vec)
    p=p_vec(i)
    for j = 1 : iter
         avgAoI_vec(j,i)=regularArrival_FCFS(p,p,k);
%         [simu_delta(j,i),D_Ygeq1(j,i),Y(j,i),avgAoI_vec_rela(j,i)] = relationBetFCFS_LCFS(p,p,lambda);
    end
end
AoIexpectation = mean(avgAoI_vec,1);
% anaAoI = 1/lambda + 1 + 4 * lambda^2 / (1 - 2 * lambda) + ...
%     (1 - 2 * lambda) * eta ./ (1- (1 - 2 * lambda) * eta) + ...
%     4*lambda^2/(1-2*lambda)*(1./(1-eta)./(1-(1-2*lambda)*eta).^2-1);
figure;
% plot(eta,AoIexpectation,'b-');
plot(eta,AoIexpectation,'ro','MarkerFaceColor','r');
hold on;
% plot(eta,anaAoI,'b-');
hold off;
grid on;
xlabel('Channel Memory');
ylabel('AoI');
title('k=5 symmetric GE channel');
legend('PerArrival');
save memoryLatency_regular_FCFS.mat