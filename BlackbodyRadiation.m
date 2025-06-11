%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Blackbody Radiation

% Author:   ZhiXiang Wei
% Email:    3011860885@qq.com
% Version:  Matlab2022b
% Config:   AMD Ryzen7 6800H && RTX3060
% Usage:    1.no need to change anything, just run.
%           2.T: float, Temperature.
%           3.t_len: int, Length of runtime.
%           4.mu_min, mu_max, dmu: float, Frequency range.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;
T = 5000; % [K]
t_len = 5e4; % [s]
mu_min = 2e13; mu_max = 1e15; dmu = 1e3; % [Hz]
CompleteProcess(T, mu_min, mu_max, dmu, t_len);

function CompleteProcess(T, mu_min, mu_max, dmu, t_len)
arguments
T double = 5000;
mu_min double = 2e13;
mu_max double = 1e15;
dmu double = 1e3;
t_len int32 = 5e4;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
h = 6.62607015*1e-34; % [Js]
c = 299792458; % [m/s]
k = 1.3806505*1e-23; % [J/K]
% Parameters
dt = 1; % [s]
mu = linspace(mu_min,mu_max,dmu); % [Hz]
n = 1 ./ (exp((h*mu)/(k*T)) - 1); % Bose-Einstein Distribution
Np = poissrnd(n); % Binomial Distribution is approximately Poisson Distribution
Ep = zeros(size(mu));
gamma = 10^-(floor(max(log10(n)))+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bose-Einstein Distribution And Initialization of Energy Level Distribution
figure()
subplot(1,2,1)
plot(mu, n*gamma, "LineStyle","none", "Marker","o")
title("玻色爱因斯坦分布")
xlabel("频率")
ylabel("粒子数")
subplot(1,2,2)
bar(mu, Np)
title("初始化能级粒子分布")
xlabel("频率")
ylabel("粒子数")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Blackbody Radiation
figure()
B1 = bar(mu, h*mu.^3.*Ep ./ (mu(2)-mu(1)));
xlabel("频率")
ylabel('辐射强度')
ylim([0,6])
for t = 1:dt:t_len
    for i = 1:length(mu)
        n_i = Np(i);
        % Spontaneous transition probability
        P_emit = gamma*n(i);
        % Stimulated transition probability
        P_absorb = gamma*n(i);
        if rand() < P_emit && n_i > 0
            Np(i) = n_i - 1;
            Ep(i) = Ep(i) + 1;
        end
        if rand() < P_absorb
            Np(i) = n_i + 1;
        end
    end
    if rem(t, 500) < 1
    set(B1, "YData",h*mu.^3.*Ep ./ (mu(2)-mu(1)));
    drawnow limitrate
    end
end
drawnow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
hold on
len = length(Ep);
for i = 1:len
    if i < 200
        color = [5e-3*(201-i),0.1+2.5e-3*(201-i),0.8-2.5e-3*(201-i)];
    end
h1 = plot([mu(i),mu(i)],[0,Ep(i)], "Color",color);
end
h2 = plot(mu, n*double(t_len)*gamma,"red","LineWidth",1.5);
xlabel('频率');
ylabel('总发射光子数');
legend([h1,h2],["仿真值","理论值"]);
set(gca,"TickDir","out")
figure()
hold on
I = h* mu.^3.*Ep ./ (mu(2)-mu(1));
I_t = (2*h*mu.^3/c^2) ./ (exp( (h*mu)/(k*T) ) - 1)*2.25e8;
for i = 1:len
    color = [1e-3*(len-i),0.1+5e-4*(len-i),0.8-5e-4*(len-i)];
    h1 = plot([mu(i),mu(i)],[0,I(i)], "Color",color);
end
h2 = plot(mu, I_t,"red","LineWidth",2);
xlabel('频率');
ylabel('辐射强度');
legend([h1,h2],["仿真值","理论值"]);
set(gca,"TickDir","out")
% xlim([min(mu),max(mu)*1e-3])
end
