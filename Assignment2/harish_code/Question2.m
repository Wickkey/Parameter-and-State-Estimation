N = 500;  %Number of observations
R = 1000; %Realizations
sig_u = sqrt(2);
sig_e = sqrt(1.0666);
f10 = 0.5;
b20 = 2;

yk = zeros(N,R);
ek = zeros(N,R);
uk = zeros(N,R);

% Assuming y[0]=0, we generate data
ek(1,:) = sig_e*randn(1,R);
yk(1,:) = ek(1,:);
uk(1,:) = sig_u*randn(1,R);

for i=2:N
    uk(i,:) = sig_u*randn(1,R);
    ek(i,:) = sig_e*randn(1,R);
    if i>2
       yk(i,:) = -f10*yk(i-1,:) + b20*uk(i-2,:) + f10*ek(i-1,:) + ek(i,:);
    end
    if i<3
       yk(i,:) = -f10*yk(i-1,:) + f10*ek(i-1,:) + ek(i,:);
    end
end

autocov = zeros(R,2); % Helper matrix

for j=1:R
    [aut,lags] = xcov(yk(:,j),'Unbiased');
    id1 = find(lags==1);
    id0 = find(lags==0);
    autocov(j,1) = aut(id0);
    autocov(j,2) = aut(id1);
end

autocov = mean(autocov);

sigma_y2 = autocov(1);  % Variance of y
sigma_yy1 = autocov(2); % Autocovariance of y at lag 1

crosscov = zeros(R,2);

for i=1:R
        [cc, lags] = xcov(yk(:,i),uk(:,i),'Unbiased');
        id1 = find(lags==1);
        id2 = find(lags==2);
        crosscov(i,1) = cc(id1);
        crosscov(i,2) = cc(id2);
end

crosscov = mean(crosscov);

sigma_yu1 = crosscov(1);  % Cross covariance at lag 1
sigma_yu2 = crosscov(2);  % Cross covariance at lag 2

tsigma_y2 = sig_e^2 + (b20^2)*(sig_u^2)/(1 - f10^2); % Theoretical variance of y
tsigma_yy1 = -(f10)*(b20^2)*(sig_u^2)/(1-f10^2); % Theoretical autocovariance at lag 1

tsigma_yu1 = 0;  % Theoretical cross covariance at lag 1
tsigma_yu2 = b20*(sig_u^2);  % Theoretical cross covariance at lag 2

