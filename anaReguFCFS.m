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


pa_b_GG = zeros(k + 1, k + 1);
% pa_b_GG(2) = p*(1-r)*r;
% pa_b_GG(3) = 2*p*r*(1-p);
% pa_b_GG(4) = (1-p)^3;
pa_b_GG(2, 2) = (1-p);


pa_b_GB = zeros(k + 1, k + 1);
% pa_b_GB(2) = p*(1-r)^2;
% pa_b_GB(3) = p*r*p+(1-p)*p*(1-r);
% pa_b_GB(4) = (1-p)^2*p;
pa_b_GB(2, 2) = p;



pa_b_BG = zeros(k + 1, k + 1);
% pa_b_BG(1) = r*(1-r)^2;
% pa_b_BG(2) = r*r*p+(1-p)*r*(1-r);
% pa_b_BG(3) = r*(1-p)^2;
pa_b_BG(2, 1) = r;



pa_b_BB = zeros(k + 1, k + 1);
% pa_b_BB(1) = (1-r)^3;
% pa_b_BB(2) = 2*p*r*(1-r);
% pa_b_BB(3) = r*(1-p)*p;
pa_b_BB(2, 1) = (1-r);

for i = 2 : k
    pa_b_BG(i + 1, 1) = r * pa_b_BB(i, 1);
    pa_b_BB(i + 1, 1) = (1-r) * pa_b_BB(i, 1);
    for j = 2 : i + 1
        pa_b_GG(i + 1, j) = (1-p) * pa_b_GG(i, j - 1) + r * pa_b_GB(i, j);
        pa_b_GB(i + 1, j) = p * pa_b_GG(i, j - 1) + (1-r) * pa_b_GB(i, j);
        pa_b_BG(i + 1, j) = (1-p) * pa_b_BG(i, j - 1) + r * pa_b_BB(i, j);
        pa_b_BB(i + 1, j) = p * pa_b_BG(i, j - 1) + (1-r) * pa_b_BB(i, j);
    end
end

% 参加迭代的只有p_0 - p_number 共计 number + 1个
m = N;
% steadyP为数据包到达发送端时队列长度的稳态分布

%======================================================
% 应该如何迭代的求解队列长度的稳态分布？解一个m+1阶的方阵
%======================================================
% while 1
%     
%      % 计算p_0，即steadyP(1)
%     temp = 0;
%     
%     for i = 0 : k - 1
%         temp = temp + steadyPG(i + 1) * sum(pa_b_GG(i + 2 : k + 1)) + ...
%             steadyPB(i + 1) * sum(pa_b_BG(i + 2 : k + 1));
%     end
%     steadyPG(1) = temp;
%     
%     temp = 0;
%     
%     for i = 0 : k - 1
%         temp = temp + steadyPG(i + 1) * sum(pa_b_GB(i + 2 : k + 1)) + ...
%             steadyPB(i + 1) * sum(pa_b_BB(i + 2 : k + 1));
%     end
%     steadyPB(1) = temp;
%     
%     % 计算p_1 : p_n
%     
%     for i = 1 : m
%         temp = 0;
%         for j = i - 1 : k + i - 1
%             temp = temp + steadyPG(j + 1) * pa_b_GG(j + 2 - i) + ...
%                 steadyPB(j + 1) * pa_b_BG(j + 2 - i);
%         end
%         steadyPG(i + 1) = temp;
%     end
%     
%     for i = 1 : m
%         temp = 0;
%         for j = i - 1 : k + i - 1
%             temp = temp + steadyPG(j + 1) * pa_b_GB(j + 2 - i) + ...
%                 steadyPB(j + 1) * pa_b_BB(j + 2 - i);
%         end
%         steadyPB(i + 1) = temp;
%     end
%     
%     steadyPG(m+2) = 1 - sum(steadyPG(1:m+1)) - sum(steadyPB(1:m+1));
%     
%     % 迭代终止的条件
%     if(norm(tempPG - steadyPG, 2) < 1.0e-5 && norm(tempPB - steadyPB, 2) < 1.0e-5)
%         break;
%     end
%     
%     tempPG = steadyPG;
%     tempPB = steadyPB;
% end


A = zeros(2 * m + 4, 2 * m + 4);
b = zeros(2 * m + 4, 1);
b(2 * m + 3, 1) = 1;
% p0
% p_0系数
A(1,1) = sum(pa_b_GG(k + 1, 2 : k + 1)) - 1;
A(1,2) = sum(pa_b_BG(k + 1, 2 : k + 1));
A(2,1) = sum(pa_b_GB(k + 1, 2 : k + 1));
A(2,2) = sum(pa_b_BB(k + 1, 2 : k + 1)) - 1;
% p_1:p_{k-1}系数
for j = 1 : k - 1
    A(1, 2 * j + 1) = sum(pa_b_GG(k + 1, j + 2 : k + 1));
    A(1, 2 * j + 2) = sum(pa_b_BG(k + 1, j + 2 : k + 1));
    A(2, 2 * j + 1) = sum(pa_b_GB(k + 1, j + 2 : k + 1));
    A(2, 2 * j + 2) = sum(pa_b_BB(k + 1, j + 2 : k + 1));
end
% p1~~p(m+2-k)
for i = 1 : m + 2 - k
    for j = i - 1 : i - 1 + k
        A(2 * i + 1, 2 * j + 1) = pa_b_GG(k + 1, j + 2 - i);
        A(2 * i + 1, 2 * j + 2) = pa_b_BG(k + 1, j + 2 - i);
        A(2 * i + 2, 2 * j + 1) = pa_b_GB(k + 1, j + 2 - i);
        A(2 * i + 2, 2 * j + 2) = pa_b_BB(k + 1, j + 2 - i);
    end
    A(2 * i + 1, 2 * i + 1) = A(2 * i + 1, 2 * i + 1) - 1;
    A(2 * i + 2, 2 * i + 2) = A(2 * i + 2, 2 * i + 2) - 1;
end
for i = m + 3 - k : m
    for j = i - 1 : m + 1
        A(2 * i + 1, 2 * j + 1) = pa_b_GG(k + 1, j + 2 - i);
        A(2 * i + 1, 2 * j + 2) = pa_b_BG(k + 1, j + 2 - i);
        A(2 * i + 2, 2 * j + 1) = pa_b_GB(k + 1, j + 2 - i);
        A(2 * i + 2, 2 * j + 2) = pa_b_BB(k + 1, j + 2 - i);
    end
    A(2 * i + 1, 2 * i + 1) = A(2 * i + 1, 2 * i + 1) - 1;
    A(2 * i + 2, 2 * i + 2) = A(2 * i + 2, 2 * i + 2) - 1;
end
A(2 * m + 3, :) = 1;
A(2 * m + 4, 2 * m + 4) = 1;
% A(2 * m + 4, 2 * m + 4) = 0;
steadyP_solve = pinv(A)*b;

steadyP = zeros(1, m + 2);
for i = 0 : m + 1
    steadyP(i + 1) = steadyP_solve(2 * i + 1) + steadyP_solve(2 * i + 2);
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

