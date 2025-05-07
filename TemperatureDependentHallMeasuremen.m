%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temperature-Dependent Hall Measuremen

% Author:   ZhiXiang Wei
% Email:    3011860885@qq.com
% Version:  Matlab2022b
% Config:   AMD Ryzen7 6800H && RTX3060
% Usage:    1.no need to change anything, just run.
%           2.l_NC: float, poles distance.
%           3.l: float, length of Hall sample.
%           4.a: float, width of Hall sample.
%           5.d: float, thickness of Hall sample.
%           6.corr: float, Correction value based on Scattering Mechanism.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;
% 常数
l_NC = 0.127*1e-3;
long = 1.0*1e-3;
a = 1.0*1e-3;
d = 0.1*1e-3;
corr = 1;
% 数据导入
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

T = UT(:,1);
V_Hp = UT(1:2:end,4);% V
V_Hn = UT(2:2:end,4);
V_sigma = UT(:,2);% V
I = UT(:,3)*1e-3;% A
B_p = UT(1:2:end,5)*1e-3;% T
B_n = UT(2:2:end,5)*1e-3;

% R
R_p = (V_Hp .* a) ./ (I(1:2:end) .* B_p);
R_n = (V_Hn .* a) ./ (I(2:2:end) .* B_n);

% Rho
sigma = (I * l_NC) ./ (V_sigma * a * d);
rho = 1 ./ sigma;

% n (corr为修正因子)
e = 1.602176634;% C
n = corr ./ (R_n * e);
p = corr ./ (R_p * e);

% Mu_H
mu_Hp = R_p .* sigma(1:2:end);
mu_Hn = R_n .* sigma(2:2:end);

% mu_n
mu_p = mu_Hp ./ (e * p .* R_p);
mu_n = mu_Hn ./ (e * n .* R_n);

