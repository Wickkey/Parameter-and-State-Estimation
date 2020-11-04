%% Q-1 (b)

N = 1000;   % Number of observations
R = 1000;  % Number of realizations
v = zeros(N,R);
v(1,:)=0;

for i=2:N
    e = randn(1,R);  %Standard Gaussian White Noise
    v(i,:) = v(i-1,:) + e;  % Random Signal
end

kvar = var(v,1,2);

plot((1:N),kvar,'linewidth',2);
ylabel('Variance','fontsize',13,'fontweight','bold');
xlabel('Time','fontsize',13,'fontweight','bold');
set(gca, 'fontsize', 12, 'fontweight', 'bold');
title('Variance vs Time','fontsize',14,'fontweight','bold');