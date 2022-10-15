function [final] = anaReguFCFS(p, r, k)
% for simplicity, we solve for K = 3, which could be extened for general K
% to solve beta in the paper easily
syms x



pa_b_GG = zeros(k + 1, k + 1);
pa_b_GG(2, 2) = (1-p);


pa_b_GB = zeros(k + 1, k + 1);
pa_b_GB(2, 2) = p;



pa_b_BG = zeros(k + 1, k + 1);
pa_b_BG(2, 1) = r;



pa_b_BB = zeros(k + 1, k + 1);
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


f=(x-pa_b_BB(4,1)-pa_b_BB(4,2)*x-pa_b_BB(4,3)*x^2-pa_b_BB(4,4)*x^3) ...
    *(x-pa_b_GG(4,1)-pa_b_GG(4,2)*x-pa_b_GG(4,3)*x^2-pa_b_GG(4,4)*x^3) ...
    -(pa_b_GB(4,1)+pa_b_GB(4,2)*x+pa_b_GB(4,3)*x^2+pa_b_GB(4,4)*x^3) ...
    *(pa_b_BG(4,1)+pa_b_BG(4,2)*x+pa_b_BG(4,3)*x^2+pa_b_BG(4,4)*x^3);

solution = solve(f==0,x);
beta0 = solution(solution>0);
beta = beta0(beta0<1);

M = zeros(3,3);


M(1,1) = 1-sum(pa_b_GG(4,2:4));
M(1,2) = -beta^0*sum(pa_b_BG(4,2:4))-beta^1*sum(pa_b_BG(4,3:4))-beta^2*sum(pa_b_BG(4,4:4));
M(1,3) = -beta^0*sum(pa_b_GG(4,3:4))-beta^1*sum(pa_b_GG(4,4:4));

M(2,1) = sum(pa_b_GB(4,2:4));
M(2,2) = beta^0*sum(pa_b_BB(4,2:4))+beta^1*sum(pa_b_BB(4,3:4))+beta^2*sum(pa_b_BB(4,4:4))-1;
M(2,3) = beta^0*sum(pa_b_GB(4,3:4))+beta^1*sum(pa_b_GB(4,4:4));

M(3,1) = 1;
M(3,2) = 1/(1-beta);
M(3,3) = 1/(1-beta);

right = zeros(3,1);
right(3,1) = 1;

res = pinv(M)*right;
final = res(1)+(1/(1-beta)+(p+r)/(r*(1-beta)^2))*res(3) ...
    +(1/(1-beta)+(r+beta*p)/(1-beta)^2)*res(2)/r + (k-1)/2;
end

