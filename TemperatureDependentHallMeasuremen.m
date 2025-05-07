%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temperature-Dependent Hall Measuremen

% Author:   ZhiXiang Wei
% Email:    3011860885@qq.com
% Version:  Matlab2022b
% Config:   AMD Ryzen7 6800H && RTX3060
% Usage:    1.no need to change anything, just run.
%           2.atom_size: int, Control the number of atom.
%           3.time: int, Control the time of optical pumping.
%           4.voltage_type: string, Control the shape of voltage, there are two
%           types: "square" and "triangle".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;
% 数据处理
UB = readmatrix("D:\100-Users\160-Desktop\变温霍尔效应\2025417\数据.xlsx", ...
    "Sheet","霍尔电阻电压磁场","Range",'B2:F1004');
UI = readmatrix("D:\100-Users\160-Desktop\变温霍尔效应\2025417\数据.xlsx", ...
    "Sheet","伏安特性","Range",'B2:E1071');
UT = readmatrix("D:\100-Users\160-Desktop\变温霍尔效应\2025417\数据.xlsx", ...
    "Sheet","霍尔电阻电压温度","Range",'B2:F1121');
UB_u = unique(UB(:,5));
UB_mean= zeros(length(UB_u),5);
for i = 1:length(UB_u)
   UB_id = UB(:,5) == UB_u(i);
   UB_mean(i,:) = sum(UB(UB_id,:),1) / sum(UB_id,1);
end

figure()
hold on
h1 = plot(UB_mean(:,5),UB_mean(:,2),"LineWidth",1.5,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(UB_mean(:,5),UB_mean(:,4),"LineWidth",1.5,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
title("霍尔/电阻电势-磁场曲线","FontSize",15,"FontWeight","bold")
xlabel("B(mT)","Interpreter","tex")
ylabel("E(V)","Interpreter","tex")
legend([h1,h2],["电阻电势","霍尔电势"],"Location","north")
ylim([-0.2,0.6])
grid on
xtickformat('%.1f')
ytickformat('%.1f')
figure()
hold on
h1 = plot(UI(:,3),UI(:,2),"LineWidth",1.5,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
title("伏安特性曲线","FontSize",15,"FontWeight","bold")
xlabel("I(mA)","Interpreter","tex")
ylabel("E(V)","Interpreter","tex")
legend(h1,"电阻电势","Location","north")
grid on
xtickformat('%.1f')
ytickformat('%.1f')
figure()
hold on
h1 = plot(UT(:,1),UT(:,2),"LineWidth",1.5,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(UT(1:2:end,1),UT(1:2:end,4),"LineWidth",1.5,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
plot(UT(2:2:end,1),UT(2:2:end,4),"LineWidth",1.5,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
title("霍尔/电阻电势-温度曲线","FontSize",15,"FontWeight","bold")
xlabel("T(K)","Interpreter","tex")
ylabel("E(V)","Interpreter","tex")
legend([h1,h2],["电阻电势","霍尔电势"],"Location","north")
ylim([-1,4])
grid on
xtickformat('%.1f')
ytickformat('%.1f')

T = UT(:,1);
V_Hp = UT(1:2:end,4);% V
V_Hn = UT(2:2:end,4);
V_sigma = UT(:,2);% V
I = UT(:,3)*1e-3;% A
B_p = UT(1:2:end,5)*1e-3;% T
B_n = UT(2:2:end,5)*1e-3;

% Constant
% l_NC = 0.127*1e-3; long = 1.0*1e-3; a = 1.0*1e-3; d = 0.1*1e-3;
l_NC = 6.0*1e-3; long = 6.0*1e-3; a = 4.0*1e-3; d = 0.2*1e-3;
% R
R_p = (V_Hp .* a) ./ (I(1:2:end) .* B_p);
R_n = (V_Hn .* a) ./ (I(2:2:end) .* B_n);
figure()
hold on
h1 = plot(T(1:2:end),abs(R_p),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(T(2:2:end),abs(R_n),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
title("霍尔系数-温度曲线","FontSize",15,"FontWeight","bold")
xlabel("T(K)","Interpreter","tex")
ylabel("kRH(m^{3}/c)","Interpreter","tex")
legend([h1,h2],["正向霍尔系数","反向霍尔系数"],"Location","north")
% ylim([-9,10])
% xlim([60,322])
grid on
xtickformat('%.1f')
ytickformat('%.1f')
% Rho
sigma = (I * l_NC) ./ (V_sigma * a * d);
rho = 1 ./ sigma;
figure()
h1 = plot(T,rho,"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
title("电阻率-温度曲线","FontSize",15,"FontWeight","bold")
xlabel("T(K)","Interpreter","tex")
ylabel("\rho(\Omega/m)","Interpreter","tex")
legend(h1,"电阻率","Location","north")
% ylim([-9,10])
% xlim([60,322])
grid on
xtickformat('%.1f')
ytickformat('%.1f')
% n
e = 1.602176634;% C
n = 1.93 ./ (R_n * e);
p = 1.93 ./ (R_p * e);
% n = 1 ./ (R_n * e);
% p = 1 ./ (R_p * e);
figure()
hold on
h1 = plot(T(1:2:end),abs(p),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(T(2:2:end),abs(n),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
title("载流子浓度-温度曲线","FontSize",15,"FontWeight","bold")
xlabel("T(K)","Interpreter","tex")
ylabel("n(10^{19}/m^{3})","Interpreter","tex")
legend([h1,h2],["正向载流子浓度","反向载流子浓度"],"Location","north")
% ylim([-9,10])
% xlim([60,322])
grid on
xtickformat('%.1f')
ytickformat('%.1f')
% Mu_H
mu_Hp = R_p .* sigma(1:2:end);
mu_Hn = R_n .* sigma(2:2:end);
figure()
hold on
h1 = plot(T(1:2:end),abs(mu_Hp),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(T(2:2:end),abs(mu_Hn),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
title("霍尔迁移率-温度曲线","FontSize",15,"FontWeight","bold")
xlabel("T(K)","Interpreter","tex")
ylabel("\mu_{H}(m^{2}/Vs)","Interpreter","tex")
legend([h1,h2],["正向霍尔迁移率","反向霍尔迁移率"],"Location","north")
% ylim([-9,10])
% xlim([60,322])
grid on
xtickformat('%.1f')
ytickformat('%.1f')

% ScatteringMechanism 
figure()
subplot(1,2,1)
hold on
h1 = plot(1./T(1:2:end),abs(R_p),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(1./T(2:2:end),abs(R_n),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
title("霍尔系数-温度倒数曲线","FontSize",15,"FontWeight","bold")
xlabel("1/T(1/K)","Interpreter","tex")
ylabel("kRH(m^{3}/c)","Interpreter","tex")
legend([h1,h2],["正向霍尔系数","反向霍尔系数"],"Location","north")
grid on
xtickformat('%.1f')
ytickformat('%.1f')
subplot(1,2,2)
hold on
h1 = plot(1./T(1:2:end),log(abs(R_p)),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(1./T(2:2:end),log(abs(R_n)),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
title("对数霍尔系数-温度倒数曲线","FontSize",15,"FontWeight","bold")
xlabel("1/T(1/K)","Interpreter","tex")
ylabel("ln(kRH)(m^{3}/c)","Interpreter","tex")
legend([h1,h2],["正向霍尔系数","反向霍尔系数"],"Location","north")
grid on
xtickformat('%.1f')
ytickformat('%.1f')
% mu_p/n
mu_p = mu_Hp ./ (e * p .* R_p);
mu_n = mu_Hn ./ (e * n .* R_n);
figure()
hold on
h1 = plot(T(1:2:end),mu_p,"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(T(2:2:end),mu_n,"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
title("载流子迁移率-温度曲线","FontSize",15,"FontWeight","bold")
xlabel("1/T(1/K)","Interpreter","tex")
ylabel("\mu_{n}(m^{2}/Vs)","Interpreter","tex")
legend([h1,h2],["正向载流子迁移率","反向载流子迁移率"],"Location","north")
grid on
xtickformat('%.1f')
ytickformat('%.1f')
% sigma_p/n
figure()
subplot(1,2,1)
hold on
h1 = plot(1./T(1:2:end),abs(sigma(1:2:end)),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(1./T(2:2:end),abs(sigma(2:2:end)),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
title("电导率-温度倒数曲线","FontSize",15,"FontWeight","bold")
xlabel("1/T(1/K)","Interpreter","tex")
ylabel("\sigma_{n}(S/m)","Interpreter","tex")
legend([h1,h2],["正向电导率","反向电导率"],"Location","north")
grid on
xtickformat('%.1f')
ytickformat('%.1f')
subplot(1,2,2)
hold on
h1 = plot(1./T(1:2:end),log(abs(sigma(1:2:end))),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(1./T(2:2:end),log(abs(sigma(2:2:end))),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
title("对数电导率-温度倒数曲线","FontSize",15,"FontWeight","bold")
xlabel("1/T(1/K)","Interpreter","tex")
ylabel("ln(\sigma_{n})(S/m)","Interpreter","tex")
legend([h1,h2],["正向电导率","反向电导率"],"Location","north")
grid on
xtickformat('%.1f')
ytickformat('%.1f')

% BandGap
B1 = 10; A1 = 100;
B2 = 10; A2 = 100;
k = 8.617333262145*1e-5;% eV/K
T_p = T(1:2:end);
T_n = T(2:2:end);
figure()
subplot(1,2,1)
logR_A1fit = fit(1./T_p(A1:end),log(abs(R_p(A1:end) .* T_p(A1:end).^(1.5))),"poly1");
logR_B1fit = fit(1./T_p(1:B1),log(abs(R_p(1:B1) .* T_p(1:B1).^(0.75))),"poly1");
hold on
h1 = plot(1./T_p(B1:A1),log(abs(R_p(B1:A1))),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(1./T_p(A1:end),log(abs(R_p(A1:end) .* T_p(A1:end).^(1.5))),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(2));
h3 = plot(1./T_p(1:B1),log(abs(R_p(1:B1) .* T_p(1:B1).^(0.75))),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
plot(1./T_p(A1:end),logR_A1fit.p1*1./T_p(A1:end)+logR_A1fit.p2,"LineWidth",1.2,"LineStyle","-","Marker","none","MarkerSize",6,"Color",color(2))
plot(1./T_p(1:B1),logR_B1fit.p1*1./T_p(1:B1)+logR_B1fit.p2,"LineWidth",1.2,"LineStyle","-","Marker","none","MarkerSize",6,"Color",color(3))
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
xlabel("1/T(1/K)","Interpreter","tex")
ylabel("ln(kRH\timesT^{3/k})(m^{3}/c)","Interpreter","tex")
title("正向霍尔系数","FontWeight","bold","FontSize",15)
legend([h1,h2,h3],["AB","A","B"],"Location","north")
grid on
xtickformat('%.1f')
ytickformat('%.1f')
subplot(1,2,2)
logRA2 = [1./T_n(A2:end),log(abs(R_n(A2:end) .* T_p(A2:end).^(1.5)))];
logRA2 = logRA2((logRA2(:,2)>-5*1e5 & logRA2(:,2)<5*1e5),:);
logR_A2fit = fit(logRA2(:,1),logRA2(:,2),"poly1");
logR_B2fit = fit(1./T_n(1:B2),log(abs(R_n(1:B2) .* T_p(1:B2).^(0.75))),"poly1");
hold on
h1 = plot(1./T_n(B2:A2),log(abs(R_n(B2:A2))),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(1./T_n(A2:end),log(abs(R_n(A2:end) .* T_n(A2:end).^(1.5))),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(2));
h3 = plot(1./T_n(1:B2),log(abs(R_n(1:B2) .* T_n(1:B2).^(0.75))),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
plot(1./T_n(A2:end),logR_A2fit.p1*1./T_n(A2:end)+logR_A2fit.p2,"LineWidth",1.2,"LineStyle","-","Marker","none","MarkerSize",6,"Color",color(2))
plot(1./T_n(1:B2),logR_B2fit.p1*1./T_n(1:B2)+logR_B2fit.p2,"LineWidth",1.2,"LineStyle","-","Marker","none","MarkerSize",6,"Color",color(3))
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
xlabel("1/T(1/K)","Interpreter","tex")
ylabel("ln(kRH\timesT^{3/k})(m^{3}/c)","Interpreter","tex")
title("反向霍尔系数","FontWeight","bold","FontSize",15)
legend([h1,h2,h3],["AB","A","B"],"Location","north")
grid on
xtickformat('%.1f')
ytickformat('%.1f')
sgtitle("对数霍尔系数-温度倒数曲线","FontSize",15,"FontWeight","bold")
% A
logR_p = log(abs(R_p(A1:end).*T_p(A1:end).^(1.5)));
logR_n = log(abs(R_n(A2:end).*T_n(A2:end).^(1.5)));
logR_p = [T_p(A1:end),logR_p];
logR_n = [T_n(A2:end),logR_n];
E_gp = 2*k * diff(logR_p(:,2)) ./ diff(1./logR_p(:,1));
E_gn = 2*k * diff(logR_n(:,2)) ./ diff(1./logR_n(:,1));
E_gp = [T_p(A1:end-1),E_gp];
E_gn = [T_n(A2:end-1),E_gn];
E_gp = E_gp(E_gp(:,2)<40 & E_gp(:,2)>-50,:);
E_gn = E_gn(E_gn(:,2)<40 & E_gn(:,2)>-50,:);
% B
logR_p = log(abs(R_p(1:B1).*T_p(1:B1).^(0.75)));
logR_n = log(abs(R_n(1:B2).*T_n(1:B2).^(0.75)));
logR_p = [T_p(1:B1),logR_p];
logR_n = [T_n(1:B2),logR_n];
E_ip = 2*k * diff(logR_p(:,2)) ./ diff(1./logR_p(:,1));
E_in = 2*k * diff(logR_n(:,2)) ./ diff(1./logR_n(:,1));
E_ip = [T_p(1:B1-1),E_ip];
E_in = [T_n(1:B2-1),E_in];
E_ip = E_ip(E_ip(:,2)<40 & E_ip(:,2)>-50,:);
E_in = E_in(E_in(:,2)<40 & E_in(:,2)>-50,:);
figure()
subplot(1,2,1)
hold on
h1 = plot(1./E_gp(:,1),E_gp(:,2),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(1./E_gn(:,1),E_gn(:,2),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
plot(1./E_gp(:,1),2*k * logR_A1fit.p1 * ones(length(E_gp(:,1))),"LineWidth",1.2,"LineStyle","-","Marker","none","MarkerSize",6,"Color",color(1));
plot(1./E_gn(:,1),2*k * logR_A2fit.p1 * ones(length(E_gn(:,1))),"LineWidth",1.2,"LineStyle","-","Marker","none","MarkerSize",6,"Color",color(3));
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
title("禁带宽度-温度曲线","FontSize",15,"FontWeight","bold")
xlabel("1/T(K)","Interpreter","tex")
ylabel("E_{g}(eV)","Interpreter","tex")
legend([h1,h2],["正向禁带宽度","反向禁带宽度"],"Location","north")
grid on
xtickformat('%.1f')
ytickformat('%.1f')
subplot(1,2,2)
hold on
h1 = plot(1./E_ip(:,1),E_ip(:,2),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(1./E_in(:,1),E_in(:,2),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
plot(1./E_ip(:,1),2*k * logR_B1fit.p1 * ones(length(E_ip(:,1))),"LineWidth",1.2,"LineStyle","-","Marker","none","MarkerSize",6,"Color",color(1));
plot(1./E_in(:,1),2*k * logR_B2fit.p1 * ones(length(E_in(:,1))),"LineWidth",1.2,"LineStyle","-","Marker","none","MarkerSize",6,"Color",color(3));
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
title("杂质电离能-温度曲线","FontSize",15,"FontWeight","bold")
xlabel("1/T(K)","Interpreter","tex")
ylabel("E_{i}(eV)","Interpreter","tex")
legend([h1,h2],["正向杂质电离能","反向杂质电离能"],"Location","north")
grid on
xtickformat('%.4f')
ytickformat('%.1f')
sgtitle("基于ln(R\timesT^{3/k})/(1/T)","FontSize",15,"FontWeight","bold")

% BandGap
B1 = 40; A1 = 65;
B2 = 40; A2 = 65;
k = 8.617333262145*1e-5;% eV/K
T_p = T(1:2:end);
T_n = T(2:2:end);
sigma_p = sigma(1:2:end);
sigma_n = sigma(2:2:end);

figure()
subplot(1,2,1)
logSigma_A1fit = fit(1./T_p(A1:end),log(sigma_p(A1:end)),"poly1");
logSigma_B1fit = fit(1./T_p(1:B1),log(sigma_p(1:B1)),"poly1");
hold on
h1 = plot(1./T_p,log(sigma_p),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(1./T_p(A1:end),log(sigma_p(A1:end)),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(2));
h3 = plot(1./T_p(1:B1),log(sigma_p(1:B1)),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
plot(1./T_p(A1:end),logSigma_A1fit.p1*1./T_p(A1:end)+logSigma_A1fit.p2,"LineWidth",1.2,"LineStyle","-","Marker","none","MarkerSize",6,"Color",color(2))
plot(1./T_p(1:B1),logSigma_B1fit.p1*1./T_p(1:B1)+logSigma_B1fit.p2,"LineWidth",1.2,"LineStyle","-","Marker","none","MarkerSize",6,"Color",color(3))
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
xlabel("1/T(1/K)","Interpreter","tex")
ylabel("ln\sigma_{n}(S/m)","Interpreter","tex")
title("正向电导率","FontWeight","bold","FontSize",15)
legend([h1,h2,h3],["AB","A","B"],"Location","north")
grid on
xtickformat('%.1f')
ytickformat('%.1f')
subplot(1,2,2)
logSigma_A2fit = fit(1./T_n(A2:end),log(sigma_n(A2:end)),"poly1");
logSigma_B2fit = fit(1./T_n(1:B2),log(sigma_n(1:B2)),"poly1");
hold on
h1 = plot(1./T_n,log(sigma_n),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(1./T_n(A2:end),log(sigma_n(A2:end)),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(2));
h3 = plot(1./T_n(1:B2),log(sigma_n(1:B2)),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
plot(1./T_n(A2:end),logSigma_A2fit.p1*1./T_n(A2:end)+logSigma_A2fit.p2,"LineWidth",1.2,"LineStyle","-","Marker","none","MarkerSize",6,"Color",color(2))
plot(1./T_n(1:B2),logSigma_B2fit.p1*1./T_n(1:B2)+logSigma_B2fit.p2,"LineWidth",1.2,"LineStyle","-","Marker","none","MarkerSize",6,"Color",color(3))
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
xlabel("1/T(1/K)","Interpreter","tex")
ylabel("ln\sigma_{n}(S/m)","Interpreter","tex")
title("反向电导率","FontWeight","bold","FontSize",15)
legend([h1,h2,h3],["AB","A","B"],"Location","north")
grid on
xtickformat('%.1f')
ytickformat('%.1f')
sgtitle("对数电导率-温度倒数曲线","FontSize",15,"FontWeight","bold")
% A
logSigma_p = log(sigma_p(A1:end).*T_p(A1:end).^(-6));
logSigma_n = log(sigma_n(A2:end).*T_n(A2:end).^(-6));
logSigma_p = [T_p(A1:end),logSigma_p];
logSigma_n = [T_n(A2:end),logSigma_n];
E_gp = -2*k * diff(logSigma_p(:,2)) ./ diff(1./logSigma_p(:,1));
E_gn = -2*k * diff(logSigma_n(:,2)) ./ diff(1./logSigma_n(:,1));
E_gp = [T_p(A1:end-1),E_gp];
E_gn = [T_n(A2:end-1),E_gn];
E_gp = E_gp(E_gp(:,2)<5e5 & E_gp(:,2)>-5e5,:);
E_gn = E_gn(E_gn(:,2)<5e5 & E_gn(:,2)>-5e5,:);
% B
logSigma_p = log(sigma_p(1:B1));
logSigma_n = log(sigma_n(1:B2));
logSigma_p = [T_p(1:B1),logSigma_p];
logSigma_n = [T_n(1:B2),logSigma_n];
E_ip = -2*k * diff(logSigma_p(:,2)) ./ diff(1./logSigma_p(:,1));
E_in = -2*k * diff(logSigma_n(:,2)) ./ diff(1./logSigma_n(:,1));
E_ip = [T_p(1:B1-1),E_ip];
E_in = [T_n(1:B2-1),E_in];
E_ip = E_ip(E_ip(:,2)<5e5 & E_ip(:,2)>-5e5,:);
E_in = E_in(E_in(:,2)<5e5 & E_in(:,2)>-5e5,:);
figure()
subplot(1,2,1)
hold on
h1 = plot(1./E_gp(:,1),E_gp(:,2),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(1./E_gn(:,1),E_gn(:,2),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
plot(1./E_gp(:,1),-2*k * logSigma_A1fit.p1 * ones(length(E_gp(:,1))),"LineWidth",1.2,"LineStyle","-","Marker","none","MarkerSize",6,"Color",color(1));
plot(1./E_gn(:,1),-2*k * logSigma_A2fit.p1 * ones(length(E_gn(:,1))),"LineWidth",1.2,"LineStyle","-","Marker","none","MarkerSize",6,"Color",color(3));
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
title("禁带宽度-温度曲线","FontSize",15,"FontWeight","bold")
xlabel("1/T(K)","Interpreter","tex")
ylabel("E_{g}(eV)","Interpreter","tex")
legend([h1,h2],["正向禁带宽度","反向禁带宽度"],"Location","north")
grid on
xtickformat('%.1f')
ytickformat('%.1f')
subplot(1,2,2)
hold on
h1 = plot(1./E_ip(:,1),E_ip(:,2),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(1));
h2 = plot(1./E_in(:,1),E_in(:,2),"LineWidth",1.2,"LineStyle","none","Marker",".","MarkerSize",6,"Color",color(3));
plot(1./E_ip(:,1),-2*k * logSigma_B1fit.p1 * ones(length(E_ip(:,1))),"LineWidth",1.2,"LineStyle","-","Marker","none","MarkerSize",6,"Color",color(1));
plot(1./E_in(:,1),-2*k * logSigma_B2fit.p1 * ones(length(E_in(:,1))),"LineWidth",1.2,"LineStyle","-","Marker","none","MarkerSize",6,"Color",color(3));
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5)
title("杂质电离能-温度曲线","FontSize",15,"FontWeight","bold")
xlabel("1/T(K)","Interpreter","tex")
ylabel("E_{i}(eV)","Interpreter","tex")
legend([h1,h2],["正向杂质电离能","反向杂质电离能"],"Location","north")
grid on
xtickformat('%.1f')
ytickformat('%.1f')
sgtitle("基于ln\sigma_{n}/(1/T)","FontSize",15,"FontWeight","bold")
