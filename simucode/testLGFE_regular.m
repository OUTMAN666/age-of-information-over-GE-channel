%This file is created by Guan Qixing on 10/11/2021
%It simlate the AoI performance of LGFC policy with instant feedback and GE
%channel

clc;
clear;
close all;

p_vec=0.1:0.05:0.5; % the transition probability in GE channel model
eta=1-2*p_vec;
k = 10;% need to be less than 0.5 to ensure stable queue
iter = 50;
avgAoI_vec=zeros(iter,length(p_vec));
for i=1:length(p_vec)
    p=p_vec(i)
    for j = 1 : iter
         avgAoI_vec(j,i)=regularArrival_LGFS(p,p,k);
    end
end
AoIexpectation = mean(avgAoI_vec,1);
anaAoI = (k + 1) / 2 + 1/2./p_vec;
figure;
plot(eta,AoIexpectation,'b-');
%plot(eta,coefficient,'ro','MarkerFaceColor','r');
hold on;
plot(eta,anaAoI,'bo','MarkerFaceColor','b');
hold off;
grid on;
xlabel('Channel Memory');
ylabel('AoI');
title('k=10 symmetric GE channel');
legend('simu','ana');
save memoryLatency_LGFS.mat