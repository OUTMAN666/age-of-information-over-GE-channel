clc;
clear;
k = 5;
lambda = 1/(k);
p=0.1;
G0 = 1/2*(1-lambda/(1-lambda));
B0 = 1/2*(p*(1-lambda)-p*lambda)/(1-lambda)/(p+(1-p-p)*lambda);

G0+B0