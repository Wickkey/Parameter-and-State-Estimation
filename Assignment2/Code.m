%% README
% run sections separately.
clear
clc
%% Q1
clear
clc
K = 100; R = 100; %R = No_of_realization,k = given parameter in question.
mu = 0; std = 1;
e = mu+ std.*randn(K,R);
v = zeros(R,1);
v_i = 0;
for i = 1:R
    v(i) = v_i + sum(e(:,i));
end    

Var_v = var(v);
fprintf("As you can see, Var(v[k]) (%.2f) is proportional to K (%d) \n",Var_v,K)
%% Q2
clear
clc
rng('default');

sigma2u = 2;f = 0.5;b = 2;N = 502;

u = sqrt(sigma2u)*randn(1,N);
ystar = zeros(1, 500);

for i = 3:502
   ystar(i) = -f*ystar(i-1)+b*u(i-2);
end

ystar = ystar(1, 3:502); u = u(1, 3:502);

sigma2e = var(ystar)/10;
e = sqrt(sigma2e)*randn(1,500);

y = ystar + e;

[yy, yy_lags] = xcov(y, 'Unbiased');
[yu, yu_lags] = xcov(y, u, 'Unbiased');

zeropos = find(yy_lags == 0);
yy0_num = yy(zeropos);
onepos = zeropos + 1;
yy1_num = yy(onepos);
yu1_num = yu(onepos);
twopos = onepos + 1;
yu_2_num = yu(twopos);

yy0_th = sigma2e + (b*b*sigma2u)/(1-f*f);
yy1_th = -f*b*b*sigma2u/(1-f*f);
yu1_th = 0;
yu2_th = b*sigma2u;
fprintf("The values obtained via numerical method is fairly good on comparing with the theoretical estimates. \n");

%% Q3
clear
clc
a = load('a2_q3.mat');

N = size(a.vk,1);

%Plotting given data
figure
plot((1:N),a.vk)
xlabel('t'); ylabel('Val')
title('Plot of data points given')
box off;

%Plotting ACF
figure
subplot(211)
autocorr(a.vk,'NumLags',20);
ylabel('Sample ACF');
xlabel('')
box off

%Plotting PACF
subplot(212)
parcorr(a.vk,'NumLags',20);
ylabel('Sample PACF');
box off

% 3-b
vdk = diff(a.vk);

%Plotting ACF
figure
subplot(211)
autocorr(vdk,'NumLags',20);
ylabel('Sample ACF');
box off

%Plotting PACF
subplot(212)
parcorr(vdk,'NumLags',20);
ylabel('Sample PACF');
box off

% Estimating parameters

mod1 = arima(2,1,1);
mod1.Constant=0;
mod1est = estimate(mod1,a.vk);

res_mod1 = infer(mod1est,a.vk);
[ht1, pval1] = lbqtest(res_mod1)
r1 = summarize(mod1est)

mod2 = arima(3,1,0);
mod2.Constant=0;
mod2est = estimate(mod2, a.vk);

res_mod2 = infer(mod2est,a.vk);
[ht2, pval2] = lbqtest(res_mod2)
r2 = summarize(mod2est)

mod3 = arima(1,1,2);
mod3.Constant = 0;
mod3est=estimate(mod3,a.vk);

res_mod3 = infer(mod3est,a.vk);
[ht3, pval3] = lbqtest(res_mod3)
r3 = summarize(mod3est)

r1.AIC
r2.AIC
r3.AIC


%% Q4
clear

% Given data (True values)
N = 100;
a = 2;
b = 3;
sigma_e = 1;
mu_e = 0;

e = mu_e + sigma_e*randn(N,1);  % Gaussian data for epsilon
x = rand(N,1);   % Uniformly distributed data for x

% Calculating values of y[k] with given true values for parameters
y = a*x + b + e;

% Different values of a and b to get the most optimal estimate
a_vec = -5:0.1:10;
b_vec = -5:0.1:10;

Ll = zeros(size(b_vec,2), size(a_vec,2)); %Likelihood function

max = -inf;
for j=1:size(a_vec,2)
    for i=1:size(b_vec,2)
        Ll(i,j) = -(1/2)*log(2*pi*sigma_e^2) - (1/2)*((y-a_vec(j)*x-b_vec(i))'*(y-a_vec(j)*x-b_vec(i)))/sigma_e^2;
        
        if Ll(i,j)>max      % Conditional statement to update graphical maximum
            max = Ll(i,j);
            a_max = a_vec(j);
            b_max = b_vec(i);
        end
        
    end
end

figure
mesh(a_vec, b_vec, Ll);
hold on;
%patch([a_max+5; a_max+5; a_max-5; a_max-5], [b_max-5; b_max+5; b_max+5; b_max-5], [max; max; max; max], 'm');
plot3(a_max,b_max,max,'.r','markersize',20);
hold on;
plot3(a,b,max,'.b','markersize',20);
set(gca,'fontsize',10,'fontweight','bold');
xlabel('a','fontsize',12,'fontweight','bold');
ylabel('b','fontsize',12,'fontweight','bold');
zlabel('Log likelihood function','fontsize',12,'fontweight','bold');
title('3-D plot of Log-likelihood function','fontsize',14,'fontweight','bold');
hold off;
box off;
legend('Likelihood function','Graphical estimate','True value', 'Location' ,'NorthEast');

% Analytical Maximum Likelihood estimates

B = [x'*y; sum(y)];
A = [ x'*x sum(x); sum(x) N];
theta_hat = A^-1*B;

