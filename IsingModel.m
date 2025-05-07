%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ising Model

% Author:   ZhiXiang Wei
% Email:    3011860885@qq.com
% Version:  Matlab2022b
% Config:   AMD Ryzen7 6800H && RTX3060
% Usage:    1.no need to change anything, just run.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;
%磁场初始化
magnet_size = 30;
if magnet_size < 2
    error("求解尺度须为二维")
end
%自旋初始化
spin = 2*randi([0,1],magnet_size,magnet_size)-1;
%求解器初始化
Temp = 10; %K
time_length = 1; %s
%记录器初始化
record_magnet = zeros(time_length,1);
record_energy = zeros(time_length,1);
record_mean_mag = zeros(time_length,1);
record_mean_eng = zeros(time_length,1);
record_cv = zeros(time_length,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iger = 0;
for time = 1:time_length
    for Temp = [10:-0.05:0,0:0.05:10]
        iger = iger + 1;
        beta = 1/Temp;
        [spin,magnet,energy,deltaE] = MCstep(magnet_size,spin,beta);
        record_magnet(iger) = mean(magnet);
        record_energy(iger) = mean(energy);
        record_mean_mag(iger) = magnet(end); 
        record_mean_eng(iger) = energy(end);
        image(spin,"CDataMapping",'scaled');
        text(0,-1,['Temperature: ',num2str(Temp)],"FontSize",15,"FontWeight","bold")
        drawnow()
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
plot(record_mean_mag);

figure()
half = ceil(length(record_mean_eng)/2);
plot(half:-1:1,record_mean_eng(1:half),'b.');
hold on
plot(1:half+1,record_mean_eng(half:end),'r.');

figure
hold on
flag = 0;
temp = [10:-0.05:0,0:0.05:10];
for i = 1:length(record_mean_eng)
    if i < half
        plot(temp(i),record_mean_eng(i),'b.')
        xlabel('T')
    else
        plot(temp(i),record_mean_eng(i),'r.')
        xlabel('T')
    end
    drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%平均磁化率
Tc = 2.269;
Temp = linspace(0.001,5,500);
mean_chi = (1-1./sinh(2./Temp).^4).^(1/8);
mean_chi(Temp > Tc) = 0;
figure()
plot(Temp,mean_chi);
title('平均磁化率')
ylabel('\chi',Interpreter='tex'),xlabel('T',Interpreter='tex')

beta=0.1:0.01:1;
for i=1:1:size(beta,2)
lambda=2.*sinh(2.*beta(i))./cosh(2.*beta(i)).^2;
phi=linspace(0,pi/2,1000);
b=1./sqrt(1-lambda.^2.*sin(phi).^2);
B=trapz(phi,b);
A=2*tanh(2.*beta(i)).^2-1;
e_bar(i)=-coth(2.*beta(i)).*(1+2/pi.*A.*B);
end
beta1=beta(2:end)/2+beta(1:end-1)/2;
cv=-beta1.^2.*diff(e_bar)./diff(beta);

figure
plot(1./beta1,cv,'r','LineWidth',2);
title('比热容')
ylabel('C',Interpreter='tex'),xlabel('T',Interpreter='tex')


function [spin,record_magnet,record_energy,deltaE] = MCstep(magnet_size,spin,beta)
for step = 1:magnet_size*magnet_size
    flip_xy = unidrnd(magnet_size,1,2);

    center_x = rem(flip_xy(1),magnet_size)+1; center_y = rem(flip_xy(2),magnet_size)+1;
    center_left = rem(flip_xy(1)-1,magnet_size)+1; center_right = rem(flip_xy(1)+1,magnet_size)+1;
    center_up = rem(flip_xy(2)-1,magnet_size)+1; center_down = rem(flip_xy(2)+1,magnet_size)+1;

    state = spin(center_x,center_y);
    state1 = spin(center_x,center_up);
    state2 = spin(center_x,center_down);
    state3 = spin(center_left,center_y);
    state4 = spin(center_right,center_y);

    deltaE = 2*state*(state1+state2+state3+state4);
    if rand > exp(-beta*deltaE)
    else
        spin(center_x,center_y) = -state;
    end
    record_magnet(step) = measure_magnet(spin);
    record_energy(step) = measure_energy(spin);
end
end

function mean_magnet = measure_magnet(spin)
mean_magnet = sum(spin,"all");
end
function mean_energy = measure_energy(spin) 
State1=spin;State1(1,:)=[];State1=[State1; spin(1,:)];
State2=spin;State2(:,1)=[];State2=[State2, spin(:,1)];
mean_energy = -mean(spin.*State1+spin.*State2,"all");
end
