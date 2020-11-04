%% Question 3-a
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

%% 3-b
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

%% Estimating parameters

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
