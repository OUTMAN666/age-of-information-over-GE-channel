function avgAoI = regularArrival_LGFS(p,r,k)
%This file is created by Guan Qixing on 10/10/2021
%It simlate the AoI performance of LGFC policy with instantaneous feedback and GE
%channel

N=10000; %the total number of time slot


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
% the pkts that have been preempted can be seen as delivered at the same
% time when the last preemptive pkt was delivered
firstAttempt=zeros(1,numPackets); % to record the first attempt until successfull delivery
succTransmission=zeros(1,numPackets); %to record the time slot until it is successfully delivered
freshness=zeros(1,N); %to record the latest packet received by the receiver
%The first packet is transmitted immediately when it is generated, and it
%is delivered when the channel becomes good
firstAttempt(1)=PacketArriveTS(1);
succTransmission(1)=firstAttempt(1)+find(GoodState(firstAttempt(1):end)==1,1,'first');
%The remaining packets are served after it is arrived and the previous one is complete 
for i=2:numPackets
    %each packet is transmitted immediately 
    firstAttempt(i) = PacketArriveTS(i);
    totalAttempts=find(GoodState(firstAttempt(i):end)==1,1,'first');
    if ~isempty(totalAttempts)
        totalDelivered=i;
        %actually, there is only the last pkt that was been delivered
        %successfully
        succTransmission(i)=firstAttempt(i)+find(GoodState(firstAttempt(i):end)==1,1,'first');
        %the latest packet between two successful transmission is the previous one
        %fressness denote the time slot that the latest received pkt was generated
        freshness(succTransmission(i-1):succTransmission(i)-1)=PacketArriveTS(i-1); 
    else
        %cannot deliver until N time slots
        freshness(succTransmission(i-1):N)=PacketArriveTS(i-1); %the previous delivered information is the latest one
        totalDelivered=i-1;
        break;
    end
end

if totalDelivered==numPackets
    %all the packets are delivered
   freshness(succTransmission(i):N)=PacketArriveTS(i);
end
age = (1:N)-freshness;
avgAoI=mean(age);
% B = age(PacketArriveTS);
% EB = mean(B);
% serviceTime = succTransmission - firstAttempt;
% Sdist=hist(serviceTime,1:max(serviceTime));
% 
% X_vec=diff(PacketArriveTS);
% Xdist=hist(X_vec,1:max(X_vec)); %return the distribution of inter-arrival time 
% EX = mean(X_vec);
% EX2 = mean(X_vec.^2);
% 
% 
% 
% serviceTimedist_allPackets = Sdist/(numPackets);
% interTimedist_allPackets = Xdist/(numPackets-1);
% jointSX = interTimedist_allPackets'*serviceTimedist_allPackets;
% minSX = zeros(length(Xdist),length(Sdist));
% pr = 0;
% for i = 1:length(Xdist)
%     for j = 1:length(Sdist)
%         minSX(i,j) = min(i,j)*jointSX(i,j);
%         if j <= i
%             pr = pr+jointSX(i,j);
%         end
%     end
% end
% minSXexpectation = sum(sum(minSX));
% R = 1/2*(X_vec.^2-X_vec)+B(2:end).*min(X_vec,serviceTime(2:end));
% simuA = mean(R)/EX;
% A = 1/2*EX2/EX+EB/EX*minSXexpectation-1/2;
% % analysis results of serviceTime distribution
% anaServiceDist = zeros(1,max(serviceTime));
% anaServiceDist(1) = r/(p+r);
% 
% anaServiceDist(2:max(serviceTime)) = p/(p+r)*r*(1-r).^((2:max(serviceTime))-2);



% verify the independence of B_i with X_i & S_i
% Bdist = hist(B,1:max(B));
% Bdist_allPackets = Bdist/(numPackets);
% B_ip = B(2:end); 
% 
% BdistGivenX=B(X_vec==1);
% Bdist_conditioned_X=hist(BdistGivenX,1:max(BdistGivenX)); %return the distribution of inter-arrival time 
% Bdist_conditionedAllPkt_X=Bdist_conditioned_X/length(BdistGivenX);
% 
% 
% BdistGivenS=B_ip(serviceTime==1);
% Bdist_conditioned_S=hist(BdistGivenS,1:max(BdistGivenS)); %return the distribution of inter-arrival time 
% Bdist_conditionedAllPkt_S=Bdist_conditioned_S/length(BdistGivenS);

% figure;
% plot(Bdist_allPackets,'bo-');
% hold on;
% plot(Bdist_conditionedAllPkt_X,'rs-');
% plot(Bdist_conditionedAllPkt_S,'gs-');
% legend('Distribution of B','distribution of B condition on X=1','distribution of B condition on S=1');

end

