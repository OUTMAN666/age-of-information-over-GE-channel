function avgAoI = getAoI_will_GE_pLGFS(p,r)
N=10000; %the total number of time slot


%Simulate the packet arrival

PacketArrive = ones(1,N);

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
end

