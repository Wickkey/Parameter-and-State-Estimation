%% README
% run the sections separately.
clc
clear
%% Q1
FUN = @(x,y) exp(-x./y).*exp(-y)./y;
Q = integral2(FUN,0,1,0.2,0.4);

fprintf("\n\nthe probability for the given range is %f \n",Q)
%% Q2
rng(1,'twister');

N1 = 1000;
mu = 1;
stddev = sqrt(2);

X = (stddev*randn(1,N1))+ mu;
Y = (3*X.^2 + 5*X)';
covtheoretical = [];
covestimate = [];

for i = 10:10:1000
    sigmamat = cov(X(1:i)',Y(1:i));
    covtheoretical = [covtheoretical sigmamat(1,2)];
    covestimate = [covestimate samplecov(i,X(1:i),Y(1:i))];
end

X_axis = [10:10:1000];

figure
plot(X_axis,covtheoretical, X_axis,covestimate)
legend('theoretical covariance','estimated covariance')
xlabel('No. of sample points');
ylabel('Covariance value')

fprintf("\n\nAt 1000 sample points,\nEstimated Covariance = %f",covestimate(100))
fprintf("\nTheoretical Covariance = %f \n",covtheoretical(100))
fprintf("Absolute difference = %f \n",abs(covtheoretical(100)-covestimate(100)))

%% Q4
%a)
rng(8,'twister');

N = 20000;
v = 10;
y = chi2rnd(v,1,N);
x_hat = [8:0.001:10];
L = length(x_hat);
Jvals = zeros(1,L);

for i = 1:L
    Jvals(i) = norm((y-x_hat(i)),1)/N;
end    

[M,index] = min(Jvals);
x_hat_opt = x_hat(index);

fprintf("\n\nOptimal MAE value is %f ",x_hat_opt);
fprintf("\nAbsolute average error is %f \n ",M);

%b)
probab1 = chi2cdf(1.1*x_hat_opt,v) - chi2cdf(0.9*x_hat_opt,v);
probab2 = chi2cdf(1.1*v,v) - chi2cdf(0.9*v,v);
fprintf("\nPr(0.9X*< X <1.1X*) = %f",probab1)
fprintf("\nPr(0.9mu< X <1.1mu) = %f \n",probab2)
%% functions

function sigmahat = samplecov(N,x,y)
    xmean = mean(x);
    ymean = mean(y);
    sigmahat = (x-xmean)*(y-ymean)/N;
end