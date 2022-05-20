function [AoI] = anaReguFCFS(p, r, N, k)

pG = r / (p + r);
pB = p / (p + r);

G_pre = pG;
B_pre = pB;

pa_b_G = zeros(N + 1, N + 1);
pa_b_B = zeros(N + 1, N + 1);

pa_b_G(1,1) = G_pre;
pa_b_B(1,1) = B_pre;

for i = 1 : N
    % G_cur = (1- p) * z * G_pre + r * B_pre;
    % B_cur = p * z * G_pre + (1 - r) * B_pre;
    pa_b_G(i + 1, 1) = r * pa_b_B(i, 1);
    pa_b_G(i + 1, 2 : i + 1) = (1 - p) * pa_b_G(i, 1 : i) + r * pa_b_B(i, 2 : i + 1);
    pa_b_B(i + 1, 1) = (1 - r) * pa_b_B(i, 1);
    pa_b_B(i + 1, 2 : i + 1) = p * pa_b_G(i, 1 : i) + (1 - r) * pa_b_B(i, 2 : i + 1);
end

pa_b = pa_b_G  + pa_b_B;


% �μӵ�����ֻ��p_0 - p_number ���� number + 1��
number = N;
% ��ʼ�����Ӿ��ȷֲ�
% steadyPΪ���ݰ����﷢�Ͷ�ʱ���г��ȵ���̬�ֲ�
steadyP = zeros(1, number + k);
steadyP(number + 1 : -1 : 1) = 0.5.^(1 : number + 1);

tempP = steadyP;




%==========================================
% Ӧ����ε����������г��ȵ���̬�ֲ���
%==========================================
while 1
    
    % ����p_1 : p_n
    
    for i = 1 : number
        temp = 0;
        for j = i - 1 : k + i - 1
            temp = temp + steadyP(j + 1) * pa_b(k + 1, j + 2 - i);
        end
        steadyP(i + 1) = temp;
    end
    
     % ����p_0����steadyP(1)
    temp = 0;
    
    for i = 0 : k - 1
        temp = temp + steadyP(i + 1) * sum(pa_b(k + 1, i + 2 : k + 1));
    end
    steadyP(1) = temp;
    
    % ������ֹ������
    if(norm(tempP - steadyP, 2) < 1.0e-4)
        break;
    end
    
    tempP = steadyP;
end

latency = 0;

for i = 0 : N
    temp = 0;
    for j = i : N
        temp = temp + (j + 1) * pa_b_G(j + 1, i + 1);
    end
    latency = latency + steadyP(i + 1) * temp;
end

AoI = latency + (k - 1) / 2;

end

