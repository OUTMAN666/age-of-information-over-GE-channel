function [AoI] = getPn2(p, r, N, k)

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

% pa_b = pa_b_G  + pa_b_B;


% 初始化服从均匀分布
% steadyP为数据包到达发送端时队列长度的稳态分布
GN = zeros(1, N + 1);
BN = zeros(1, N + 1);

% 用随机到达情况的稳态分布替代，在相同到达速率的情况下
lambda = 1/(k);
G0 = 1/2*(1-lambda/(1-lambda));
GN(1) = G0;
B0 = 1/2*(r*(1-lambda)-p*lambda)/(1-lambda)/(r+(1-p-r)*lambda);
BN(1) = B0;
beta = (p*lambda+(1-p-r)*lambda*(1-lambda))/(1-lambda)/(r+(1-p-r)*lambda);
G1 = p*lambda*G0*beta/(p*lambda+(1-p-r)*lambda*(1-lambda));
GN(2) = G1;
BN(2) = B0*beta;

for i = 3 : N + 1
    GN(i) = beta * GN(i - 1);
    BN(i) = beta * BN(i - 1);
end

steadyP = GN + BN;

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

