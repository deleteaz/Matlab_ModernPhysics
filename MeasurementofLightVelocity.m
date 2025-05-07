%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Measurement of Light Velocity 

% Author:   ZhiXiang Wei
% Email:    3011860885@qq.com
% Version:  Matlab2022b
% Config:   AMD Ryzen7 6800H && RTX3060
% Usage:    1.no need to change anything, just run.
%           2.n1: float, Refractive index of medium 1.
%           3.n2: float, Refractive index of medium 2.
%           4.degree: float, The Angle of the incident light with respect
%           to the X-axis.
%           5.time: int, The time the light travels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;
%折射率
n1 = 1;
n2 = 1.01;
degree = 0.1;
time = 15;

figure()
hold on
ylim([0,0.1])
%折射光
[x1, y1, x2, t2, y2, y3] = Initialize_i(degree, n1, n2);
for i = 1:time
[x1, t1, y1, x2, t2, y2, y3] = Loop_i(x2, t2, y2, n1, n2);
block = min(x2);
xline(block)
LightDraw_i(x1, t1, y1, x2, t2, time)

%二次反射光
[x1o, t1o, y1o, x2o, t2o, y2o, y3o] = Loop_o(x1, t1, y3, n1, n2);
LightDraw_o(x1o, t1o, y1o, x2o, t2o, time)
drawnow()
pause(0.1)
end

function [x1, y1, x2, t2, y2, y3] = Initialize_i(degree, n1, n2)
arguments
    degree double = 1
    n1 double = 1
    n2 double = 1.5
end
x1 = [0:1]';
x2 = max(x1) + [0:length(x1)-1]';
t1 = pi/180 * (degree);
t2 = asin(n1/n2 * sin(t1));
y1 = tan(t1) * x1;
y2 = tan(t2) * (x2 - x2(1)) + y1(end);
y3 = -tan(t1) * (x1 - x1(end)) + y1(end);
end

function [x1, y1, x2, t2, y2, y3] = Initialize_o(degree, n1, n2)
arguments
    degree double = 1
    n1 double = 1
    n2 double = 1.5
end
x1 = [0:1]';
x2 = min(x1) + [-length(x1)+1:0]';
t1 = pi/180 * (-degree);
t2 = asin(n1/n2 * sin(t1));
y1 = tan(t1) * x1;
y2 = tan(t2) * (x2 - x2(end)) + y1(1);
y3 = -tan(t1) * (x1 - x1(1)) + y1(1);
end

function [x1, t1, y1, x2, t2, y2, y3] = Loop_i(x2, t2, y2, n1, n2)
arguments
    x2 double
    t2 double
    y2 double
    n1 double = 1
    n2 double = 1.5
end
x1 = x2; x2 = max(x1) + [0:length(x1)-1]';
t1 = t2; t2 = asin(n1/n2 * sin(t1));
y1 = y2;
y2 = tan(t2) * (x2 - x2(1)) + y1(end);
y3 = -tan(t1) * (x1 - x1(end)) + y1(end);
end

function [x1, t1, y1, x2, t2, y2, y3] = Loop_o(x2, t2, y2, n1, n2)
arguments
    x2 double
    t2 double
    y2 double
    n1 double = 1
    n2 double = 1.5
end
x1 = x2; x2 = min(x1) + [-length(x1)+1:0]';
t1 = t2; t2 = asin(n2/n1 * sin(t1));
y1 = y2;
y2 = tan(t2) * (x2 - x2(end)) + y1(1);
y3 = -tan(t1) * (x1 - x1(1)) + y1(1);
end

function LightDraw_i(x1, t1, y1, x2, t2, time)
arguments
    x1 double
    t1 double
    y1 double
    x2 double
    t2 double
    time double = 15
end
x2long = [x2; [max(x2):time+3]'];
x3long = [1: max(x1)]';
plot(x1,y1,"LineWidth",2,"Color",[0.3,0.5,1])
plot(x2long, tan(t2) * (x2long - x2long(1)) + y1(end),"Color",[0.3,0.5,1])
plot(x3long, -tan(t1) * (x3long - x1(end)) + y1(end),"Color",[1,0.5,0.3])
end

function LightDraw_o(x1, t1, y1, x2, t2, time)
arguments
    x1 double
    t1 double
    y1 double
    x2 double
    t2 double
    time double = 15
end
x2long = [1:max(x2)]';
x3long = [min(x1): time+3]';
plot(x1,y1,"LineWidth",1,"Color",[1,0.7,0.5])
plot(x2long, -tan(t2) * (x2long - x2long(end)) + y1(1),"Color",[1,0.7,0.5])
plot(x3long, tan(t1) * (x3long - x1(1)) + y1(1),"Color",[0.5,0.7,1])
end
