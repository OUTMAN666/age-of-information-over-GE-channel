function [delta,ED_Ygeq1,EY,avgAoI] = relationBetFCFS_LCFS(p,r,k)
%This file is created by Xu Xiaoli on 23/09/2021
%It simlate the relationship of AoI performance between FCFS and LCFS

N =	10000; %the total number of time slot


%Simulate the packet arrival
PacketArrive = zeros(1,N);

PacketArrive(1 : k : N) = 1;

%Channel Parameters

% PG=0; %erasure probability when the channel is in Good State
% PB=1; %erasure probability when the channel is in Bad State

%Simulate the channel, Assume start from random state
GoodState=zeros(1,N); %Trace the channel state
Transition=(rand(1,N)<p); %whether a transition will occur at certain time
GoodState(1)=(rand<(r/(r+p))); %whether the initial state is G or B

NumTransitions=cumsum(Transition);%total number of Transitions 
flag=mod(NumTransitions,2);
GoodState(2:end)=mod(GoodState(1)+flag(1:end-1),2); 

%================Get the distribution of the inter-arrival time=============
PacketArriveTS=find(PacketArrive==1);%The time slot with new packet arrival
numPackets=length(PacketArriveTS);


%=============Get the distribution of the service time============

firstAttempt=zeros(1,numPackets); % to record the first attempt until successfull delivery
succTransmission=zeros(1,numPackets); %to record the time slot until it is successfully delivered
freshness=zeros(1,N); %to record the latest packet received by the receiver
%The first packet is transmitted immediately when it is generated, and it
%is delivered when the channel becomes good
firstAttempt(1)=PacketArriveTS(1);
succTransmission(1)=firstAttempt(1)+find(GoodState(firstAttempt(1):end)==1,1,'first');
%The remaining packets are served after it is arrived and the previous one is complete 
deltaTime = zeros(1,numPackets);
deltaTime(1) = find(PacketArrive(PacketArriveTS(1):succTransmission(1)-1)==1,1,'last')-1;
for i=2:numPackets
    %each packet is transmitted immediately 
    firstAttempt(i)=max(PacketArriveTS(i),succTransmission(i-1));
    totalAttempts=find(GoodState(firstAttempt(i):end)==1,1,'first');
    if ~isempty(totalAttempts)
        totalDelivered=i;
        succTransmission(i)=firstAttempt(i)+find(GoodState(firstAttempt(i):end)==1,1,'first');
        freshness(succTransmission(i-1):succTransmission(i)-1)=PacketArriveTS(i-1); %the latest packet between two successful transmission is the previous one
    else
        %cannot deliver until N time slots
        freshness(succTransmission(i-1):N)=PacketArriveTS(i-1); %the previous delivered information is the latest one
        totalDelivered=i-1;
        break;
    end
    % record the interval from the first arrival pkt to the last arrival pkt
    % among the  living time of this pkt
    deltaTime(i) = find(PacketArrive(PacketArriveTS(i):succTransmission(i)-1)==1,1,'last')-1;
end
successpkt = find(diff(succTransmission(1:end))>0,1,'last');
successinter = diff(succTransmission(1:successpkt));
successarrival = deltaTime(1:successpkt-1);
EY = mean(successarrival);
mulDY = successarrival.*successinter;
EDY = mean(mulDY);
ED_Ygeq1 = mean(successinter(successarrival>0));
EDmulEY = EY*ED_Ygeq1;
delta = EDY*lambda;
if totalDelivered==numPackets
    %all the packets are delivered
   freshness(succTransmission(i):N)=PacketArriveTS(i);
end
latency = succTransmission-PacketArriveTS;
T = latency(1:successpkt-1);
avgAoI=mean((1:N)-freshness)-delta;
anaTemp = ((2:20)*lambda-1+(1-lambda).^(2:20))/lambda;
simuTemp = zeros(1,19);
for i = 2:20
    simuTemp(i-1) = mean(successarrival(T==(i)));
end
% YT2 = successarrival(T==2)
% YdistT2 = hist(YT2,0:max(YT2))/length(YT2);
% YT3 = successarrival(T==3)
% YdistT3 = hist(YT3,0:max(YT3))/length(YT3);
anaidstYT10 = zeros(1,20);
anaidstYT10(1) = (1-lambda)^19;
anaidstYT10(2:20) = lambda*(1-lambda).^(19-(1:19));
YT10 = successarrival(T==20);
YdistT10 = hist(YT10,0:max(YT10))/length(YT10);

G0 = r/(p+r)*(1-p*lambda/r/(1-lambda));
D = p*lambda*G0/(p*lambda+(1-p-r)*lambda*(1-lambda));
B0 = p/(p+r)*(r*(1-lambda)-p*lambda)/(1-lambda)/(r+(1-p-r)*lambda);
beta = (p*lambda+(1-p-r)*lambda*(1-lambda))/(1-lambda)/(r+(1-p-r)*lambda);
alpha = beta/(lambda-lambda*beta+beta);
anaPt = zeros(1,max(T));
% anaPt(1) = G0-D+lambda/(lambda-lambda*beta+beta)*(B0+D);
anaPt(2:max(T)) = lambda/(lambda-lambda*beta+beta)*(B0+D)*alpha.^((2:max(T))-1);
anaPt(1) = 1-alpha+(G0-D)*alpha;
simuPt = hist(T,1:max(T))/length(T);


% figure;
% plot(1:max(T),anaPt,'bo-');
% hold on;
% plot(1:max(T),simuPt,'rs-');
% % plot(0:max(YT4),YdistT4,'k*');
% legend('Distribution of T ana','Distribution of T simu');



% interdistGivenArrival1=successinter(successarrival==1);
% interdist_conditioned_Arrival1=hist(interdistGivenArrival1,1:max(interdistGivenArrival1)) ...
%     /length(interdistGivenArrival1);
% 
% interdistGivenArrival2=successinter(successarrival==2);
% interdist_conditioned_Arrival2=hist(interdistGivenArrival2,1:max(interdistGivenArrival2)) ...
%     /length(interdistGivenArrival2); 
% 
% 
% interdistGivenArrival3=successinter(successarrival==3);
% interdist_conditioned_Arrival3=hist(interdistGivenArrival3,1:max(interdistGivenArrival3)) ...
%     /length(interdistGivenArrival3);
% 
% anadistDgivenY(1) = 1-p;
% anadistDgivenY(2:20) = p*r*(1-r).^((2:20)-2);
% 
% figure;
% plot(interdist_conditioned_Arrival1,'bo');
% hold on;
% plot(interdist_conditioned_Arrival2,'rs');
% plot(interdist_conditioned_Arrival2,'k*');
% plot(anadistDgivenY,'gd--');
% xlabel('D');
% ylabel('probability');
% legend('Distribution of interdeliver condition on interarrival=1', ...
%     'Distribution of interdeliver condition on interarrival=2', ...
%     'Distribution of interdeliver condition on interarrival=3', ...
%     'ana-dist-interdeliver');

end