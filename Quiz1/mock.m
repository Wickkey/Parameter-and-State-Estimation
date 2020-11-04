N = 1000;R=400;tol=1e-3;  %change n and r later
y = unifrnd(1,2,N*R,1);  %its a unif distribution, so works the same if we put in one vector.
% yk = mean(y');
% yk = yk';
med = median(y);
arr = abs(y-med);
sigma_estimator = median(arr);
sigma_true = var(y)^0.5;

if abs(sigma_true-sigma_estimator)<tol
    bflag = 0;
else
    bflag = 1;
end

%Q2
a = linspace(0,2,1000);
x = a*sigma_estimator;
[M,I] = min((x-sigma_true).^2);
alphaopt = a(I);