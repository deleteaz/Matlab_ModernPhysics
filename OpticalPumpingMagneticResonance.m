%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optical pumping magnetic resonance

% Author:   Deleteaz
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
atom_size = 100;
time = 30;
voltagetype = "square";

[x_SF1, x_SF2, y_SF1, y_SF2] = Initialize(atom_size);
AtomDraw(atom_size, x_SF1, x_SF2, y_SF1, y_SF2);
[x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2] = ...
    ExcitedTransition(atom_size, x_SF1, x_SF2, y_SF1, y_SF2);
AtomDraw(atom_size, x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2);
[x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2] = ...
    DeExcitedTransition(x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2);
AtomDraw(atom_size, x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2);

OpticalPumping(atom_size, time)

OpticalPumpMagneticResonance(atom_size, time, voltagetype)

function [x_SF1, x_SF2, y_SF1, y_SF2] = Initialize(atom_size) 
%原子初始化
arguments
    atom_size int32 = 100;
end
% J-I耦合
x_SF = [1:atom_size]'; 
y_SF = randi([1,2],atom_size,1);
con_SF1 = (y_SF>1); 
con_SF2 = (y_SF<2);
len_S1 = length(y_SF(con_SF1)); 
len_S2 = length(y_SF(con_SF2));
% 塞曼效应
x_SF1 = x_SF(con_SF1);  
x_SF2 = x_SF(con_SF2);
y_SF1 = randi([-1,1],len_S1,1);  
y_SF2 = randi([-2,2],len_S2,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2, light_I] = ...
    ExcitedTransition(atom_size ,x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2) 
%原子受激跃迁 dF = 0,+-1 || dM = +1
arguments
    atom_size int32 = 100
    x_SF1 (:,1) {mustBeNumeric} = []
    x_SF2 (:,1) {mustBeNumeric} = []
    y_SF1 (:,1) {mustBeNumeric} = []
    y_SF2 (:,1) {mustBeNumeric} = []
    x_PF1 (:,1) {mustBeNumeric} = []
    x_PF2 (:,1) {mustBeNumeric} = []
    y_PF1 (:,1) {mustBeNumeric} = []
    y_PF2 (:,1) {mustBeNumeric} = []
end
%原子跃迁条件
dF_S1_2_P1 = y_SF1 < 1; F_S1_2_P1 = randi([0,1],length(y_SF1),1);
dF_S1_2_P2 = y_SF1 < 2; F_S1_2_P2 = ~F_S1_2_P1;
dF_S2_2_P1 = y_SF2 < 1; F_S2_2_P1 = randi([0,1],length(y_SF2),1);
dF_S2_2_P2 = y_SF2 < 2; F_S2_2_P2 = ~F_S2_2_P1;
S1_2_P1 = F_S1_2_P1 & dF_S1_2_P1;
S1_2_P2 = F_S1_2_P2 & dF_S1_2_P2;
S2_2_P1 = F_S2_2_P1 & dF_S2_2_P1;
S2_2_P2 = F_S2_2_P2 & dF_S2_2_P2;
len_S1_2_P1 = sum(S1_2_P1); len_S1_2_P2 = sum(S1_2_P2); 
len_S2_2_P1 = sum(S2_2_P1); len_S2_2_P2 = sum(S2_2_P2); 
len_P1 = len_S1_2_P1 + len_S2_2_P1;
len_P2 = len_S1_2_P2 + len_S2_2_P2;
light_I = atom_size - (len_P1 + len_P2);

%跃迁原子
y_PF1 = [y_PF1; [y_SF1(S1_2_P1); y_SF2(S2_2_P1)] + 1]; 
y_PF2 = [y_PF2; [y_SF1(S1_2_P2); y_SF2(S2_2_P2)] + 1];
x_PF1 = [x_PF1; x_SF1(S1_2_P1); x_SF2(S2_2_P1)]; 
x_PF2 = [x_PF2; x_SF1(S1_2_P2); x_SF2(S2_2_P2)];
%未跃迁原子
y_SF1 = [y_SF1(~S1_2_P1 & ~S1_2_P2)]; 
y_SF2 = [y_SF2(~S2_2_P1 & ~S2_2_P2)];
x_SF1 = [x_SF1(~S1_2_P1 & ~S1_2_P2)]; 
x_SF2 = [x_SF2(~S2_2_P1 & ~S2_2_P2)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2, light_I] = ...
    DeExcitedTransition(x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2)
%原子退激跃迁 dF = 0,+-1 || dM = 0,+-1
arguments
    x_SF1 (:,1) {mustBeNumeric} = []
    x_SF2 (:,1) {mustBeNumeric} = []
    y_SF1 (:,1) {mustBeNumeric} = []
    y_SF2 (:,1) {mustBeNumeric} = []
    x_PF1 (:,1) {mustBeNumeric} = []
    x_PF2 (:,1) {mustBeNumeric} = []
    y_PF1 (:,1) {mustBeNumeric} = []
    y_PF2 (:,1) {mustBeNumeric} = []
end
%原子跃迁条件
dF_P1_2_S1 = y_PF1 < 1 & y_PF1 > -1; F_P1_2_S1 = randi([0,1],length(y_PF1),1);
dF_P1_2_S2 = y_PF1 < 2 & y_PF1 > -2; F_P1_2_S2 = ~F_P1_2_S1;
dF_P2_2_S1 = y_PF2 < 1 & y_PF2 > -1; F_P2_2_S1 = randi([0,1],length(y_PF2),1);
dF_P2_2_S2 = y_PF2 < 2 & y_PF2 > -2; F_P2_2_S2 = ~F_P2_2_S1;
dFm_P1_2_S1 = y_PF1 == 1;
dFm_P1_2_S2 = y_PF1 == 2;
dFm_P2_2_S1 = y_PF2 == 1;
dFm_P2_2_S2 = y_PF2 == 2; 
P1_2_S1 = F_P1_2_S1 & dF_P1_2_S1; 
P1_2_S2 = F_P1_2_S2 & dF_P1_2_S2;
P2_2_S1 = F_P2_2_S1 & dF_P2_2_S1;
P2_2_S2 = F_P2_2_S2 & dF_P2_2_S2;
mP1_2_S1 = F_P1_2_S1 & dFm_P1_2_S1; 
mP1_2_S2 = F_P1_2_S2 & dFm_P1_2_S2;
mP2_2_S1 = F_P2_2_S1 & dFm_P2_2_S1;
mP2_2_S2 = F_P2_2_S2 & dFm_P2_2_S2;
len_P1_2_S1 = sum(P1_2_S1); len_P1_2_S2 = sum(P1_2_S2); 
len_P2_2_S1 = sum(P2_2_S1); len_P2_2_S2 = sum(P2_2_S2); 
mlen_P1_2_S1 = sum(mP1_2_S1); mlen_P1_2_S2 = sum(mP1_2_S2); 
mlen_P2_2_S1 = sum(mP2_2_S1); mlen_P2_2_S2 = sum(mP2_2_S2); 
len_S1 = len_P1_2_S1 + len_P2_2_S1; mlen_S1 = mlen_P1_2_S1 + mlen_P2_2_S1;
len_S2 = len_P1_2_S2 + len_P2_2_S2; mlen_S2 = mlen_P1_2_S2 + mlen_P2_2_S2;
light_I = (len_S1 + len_S2 + mlen_S1 + mlen_S2);

%跃迁原子
y_SF1 = [y_SF1; [y_PF1(P1_2_S1); y_PF2(P2_2_S1)] + randi([-1,1],len_S1,1); [y_PF1(mP1_2_S1); y_PF2(mP2_2_S1)] + randi([-1,0],mlen_S1,1)]; 
y_SF2 = [y_SF2; [y_PF1(P1_2_S2); y_PF2(P2_2_S2)] + randi([-1,1],len_S2,1); [y_PF1(mP1_2_S2); y_PF2(mP2_2_S2)] + randi([-1,0],mlen_S2,1)];
x_SF1 = [x_SF1; x_PF1(P1_2_S1); x_PF2(P2_2_S1); x_PF1(mP1_2_S1); x_PF2(mP2_2_S1)]; 
x_SF2 = [x_SF2; x_PF1(P1_2_S2); x_PF2(P2_2_S2); x_PF1(mP1_2_S2); x_PF2(mP2_2_S2)];
%未跃迁原子
y_PF1 = [y_PF1(~P1_2_S1 & ~P1_2_S2 & ~mP1_2_S1 & ~mP1_2_S2)]; 
y_PF2 = [y_PF2(~P2_2_S1 & ~P2_2_S2 & ~mP2_2_S1 & ~mP2_2_S2)];
x_PF1 = [x_PF1(~P1_2_S1 & ~P1_2_S2 & ~mP1_2_S1 & ~mP1_2_S2)]; 
x_PF2 = [x_PF2(~P2_2_S1 & ~P2_2_S2 & ~mP2_2_S1 & ~mP2_2_S2)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hS1, hS2, hP1, hP2] = ...
    AtomDraw(atom_size, x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2) 
%原子跃迁图像
arguments
    atom_size int32 = 100
    x_SF1 (:,1) = []
    x_SF2 (:,1) = []
    y_SF1 (:,1) = []
    y_SF2 (:,1) = []
    x_PF1 (:,1) = []
    x_PF2 (:,1) = []
    y_PF1 (:,1) = []
    y_PF2 (:,1) = []
end
figure()
hold on
hS1 = plot(x_SF1,y_SF1-2,"LineStyle","none","Marker","o","Color",[0.3,0.5,1]);
hS2 = plot(x_SF2,y_SF2-7,"LineStyle","none","Marker","x","Color",[1,0.5,0.3]);
hP1 = plot(x_PF1,y_PF1+2,"LineStyle","none","Marker","o","Color",[0.3,0.5,1]);
hP2 = plot(x_PF2,y_PF2+7,"LineStyle","none","Marker","x","Color",[1,0.5,0.3]);

F = string(["-2","-1","0","+1","+2",nan,"-1","0","+1",nan,"-1","0","+1",nan,"-2","-1",0,"+1","+2"]);
for i = [-9,-8,-7,-6,-5,-3,-2,-1,1,2,3,5,6,7,8,9]
    text(atom_size,i,F(i+10))
end
text(-5,4,"P"); text(-3,-7,"F=2"); text(-3,-2,"F=1");
text(-5,-4,"S"); text(-3,7,"F=2"); text(-3,2,"F=1");
axis off
yline(0,"LineStyle","-","LineWidth",2)
yline([-9,-8,-7,-6,-5,-3,-2,-1,1,2,3,5,6,7,8,9],"LineStyle",":")
ylim([-10,10])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OpticalPumping(atom_size, time)
%光抽运
arguments
    atom_size int32 = 100;
    time int32 = 30;
end
%原子初始化
[x_SF1, x_SF2, y_SF1, y_SF2] = Initialize(atom_size);
[x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2, light_I(1,1)] = ...
    ExcitedTransition(atom_size, x_SF1, x_SF2, y_SF1, y_SF2);
[x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2, light_I(1,2)] = ...
    DeExcitedTransition(x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2);
[hS1, hS2, hP1, hP2] = AtomDraw(atom_size, x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2);
delete([hS1,hS2,hP1,hP2]);
%原子光抽运
for i = 1:time
    [x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2, light_I(i+1,1)] = ...
        ExcitedTransition(atom_size, x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2);
    [x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2, light_I(i+1,2)] = ...
        DeExcitedTransition(x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2);
    hS1 = plot(x_SF1,y_SF1-2,"LineStyle","none","Marker","o","Color",[0.3,0.5,1]);
    hS2 = plot(x_SF2,y_SF2-7,"LineStyle","none","Marker","x","Color",[1,0.5,0.3]);
    hP1 = plot(x_PF1,y_PF1+2,"LineStyle","none","Marker","o","Color",[0.3,0.5,1]);
    hP2 = plot(x_PF2,y_PF2+7,"LineStyle","none","Marker","x","Color",[1,0.5,0.3]);
    drawnow()
    pause(0.1)
    delete([hS1,hS2,hP1,hP2]);
end
hS1 = plot(x_SF1,y_SF1-2,"LineStyle","none","Marker","o","Color",[0.3,0.5,1]);
hS2 = plot(x_SF2,y_SF2-7,"LineStyle","none","Marker","x","Color",[1,0.5,0.3]);
hP1 = plot(x_PF1,y_PF1+2,"LineStyle","none","Marker","o","Color",[0.3,0.5,1]);
hP2 = plot(x_PF2,y_PF2+7,"LineStyle","none","Marker","x","Color",[1,0.5,0.3]);
figure()
hold on
plot(light_I,"LineWidth",2)
yline(0)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OpticalPumpMagneticResonance(atom_size, time, voltage_type)
%光泵磁共振
arguments
    atom_size int32 = 100;
    time double = 30;
    voltage_type = "square"
end
%原子初始化
[x_SF1, x_SF2, y_SF1, y_SF2] = Initialize(atom_size);
[x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2, light_I(1,1)] = ...
    ExcitedTransition(atom_size, x_SF1, x_SF2, y_SF1, y_SF2);
[x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2, light_I(1,2)] = ...
    DeExcitedTransition(x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2);
[hS1, hS2, hP1, hP2] = AtomDraw(atom_size, x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2);
delete([hS1,hS2,hP1,hP2]);
%原子光抽运
for i = 1:2*time
    [x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2, light_I(i+1,1)] = ...
        ExcitedTransition(atom_size, x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2);
    [x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2, light_I(i+1,2)] = ...
        DeExcitedTransition(x_SF1, x_SF2, y_SF1, y_SF2, x_PF1, x_PF2, y_PF1, y_PF2);
    %电压变换
    if (voltage_type == "square")
        %方波电压变换
        if i == time
            magnet = length([y_SF1;y_SF2]);
            y_SF = randi([1,2],magnet,1);
            x_SF = [1:atom_size]'; 
            con_SF1 = (y_SF>1);
            con_SF2 = (y_SF<2);
            len_S1 = length(y_SF(con_SF1));
            len_S2 = length(y_SF(con_SF2));
            y_SF1 = randi([-1,1],len_S1,1);
            y_SF2 = randi([-2,2],len_S2,1);
            x_SF1 = x_SF(con_SF1);
            x_SF2 = x_SF(con_SF2);
        end
    else
        %三角波电压变换
        k = 1/15;
        if i < time/2
            angle_shiftF1 = randperm(length(y_SF1), floor((1-k*i)*length(y_SF1)));
            angle_shiftF2 = randperm(length(y_SF2), floor((1-k*i)*length(y_SF2)));
        elseif i >= time/2 && i < time
            angle_shiftF1 = randperm(length(y_SF1), floor((-1+k*i)*length(y_SF1)));
            angle_shiftF2 = randperm(length(y_SF2), floor((-1+k*i)*length(y_SF2)));
        elseif i >= time && i < 3/2*time
            angle_shiftF1 = randperm(length(y_SF1), floor((3-k*i)*length(y_SF1)));
            angle_shiftF2 = randperm(length(y_SF2), floor((3-k*i)*length(y_SF2)));
        end
        if ~isempty(angle_shiftF1)
            y_SF1(angle_shiftF1) = randi([-1,1],length(angle_shiftF1),1);
        end
        if ~isempty(angle_shiftF2)
            y_SF2(angle_shiftF2) = randi([-2,2],length(angle_shiftF2),1);
        end
    end
    hS1 = plot(x_SF1,y_SF1-2,"LineStyle","none","Marker","o","Color",[0.3,0.5,1]);
    hS2 = plot(x_SF2,y_SF2-7,"LineStyle","none","Marker","x","Color",[1,0.5,0.3]);
    hP1 = plot(x_PF1,y_PF1+2,"LineStyle","none","Marker","o","Color",[0.3,0.5,1]);
    hP2 = plot(x_PF2,y_PF2+7,"LineStyle","none","Marker","x","Color",[1,0.5,0.3]);
    drawnow()
    pause(0.1)
    delete([hS1,hS2,hP1,hP2]);
end
hS1 = plot(x_SF1,y_SF1-2,"LineStyle","none","Marker","o","Color",[0.3,0.5,1]);
hS2 = plot(x_SF2,y_SF2-7,"LineStyle","none","Marker","x","Color",[1,0.5,0.3]);
hP1 = plot(x_PF1,y_PF1+2,"LineStyle","none","Marker","o","Color",[0.3,0.5,1]);
hP2 = plot(x_PF2,y_PF2+7,"LineStyle","none","Marker","x","Color",[1,0.5,0.3]);
figure()
hold on
x = [1:size(light_I,1)-1]';
plot(light_I(:,1),"LineWidth",2)
if (voltage_type == "square")
    plot(x, 5 * [ones(size(x,1)/2, 1); -ones(size(x,1)/2, 1)],"Color", [1,0.5,0.3], "LineWidth", 2)
else
    plot(x, 15 * [k*x(1:time/2); 2-k*x(time/2+1:time); 2-k*x(time+1:3/2*time); k*x(3/2*time+1:end)-4],"Color", [1,0.5,0.3], "LineWidth", 2)
end
end
