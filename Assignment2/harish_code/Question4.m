%% Question 4 (b)

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

%% Analytical Maximum Likelihood estimates

B = [x'*y; sum(y)];
A = [ x'*x sum(x); sum(x) N];
theta_hat = A^-1*B;