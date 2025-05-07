%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface Plasmon Resonance

% Author:   ZhiXiang Wei
% Email:    3011860885@qq.com
% Version:  Matlab2022b
% Config:   AMD Ryzen7 6800H && RTX3060
% Usage:    1.no need to change anything, just run.
%           2.t_step: float, the step of time.
%           3.t_end: int, Control the time of single electron.
%           4.display: int, which figure displayed.
%           5.omega: int, the angular velocity of wave.
%           6.fai: int, the phase position of wave.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;
t_step = 0.5;
t_end = 25;
display = 0;
electron_single(t_step, t_end, display);

t_step = 0.05;
omega = 5;
fai = 4;
display = 1;
electron_group(t_step, omega, fai, display);

display = 1;
CompleteProcess(display);

function electron_single(t_step, t_end, display)
arguments
    t_step double = 0.5;
    t_end int32 = 25;
    display = 1;
end
% Drude模型(电子在电磁场与库仑力的耦合作用下进行振荡行为)
% dp/dt = -qE- qvxB - p/tao;
% m(y'') = -qE - q(y')xB - m(y')/tao
% y = [y y'] => [y' y''] = [y(2); -qE/m-q[y(2)xB]-y(2)/tao]
figure()
hold on
ylim([-0.7,0])
tspan = [0, t_step];
y0 = [0, 0.01];
E = tspan(end);
% tao = 10;
ball = @(t, y) [y(2); -sin(2*E)-y(2)];
oset = odeset("MaxStep",10);
[t_ball, y_ball] = ode45(ball, tspan, y0, oset);
for i = 1:t_end
    tspan = t_ball(end) + [0; t_step];
    y0 = [y_ball(end,1), y_ball(end,2)];
    E = tspan(end);
    ball = @(t, y) [y(2); -sin(2*E)-y(2)];
    [t_ball, y_ball] = ode45(ball, tspan, y0, oset);
    x_ball = zeros(length(t_ball),1);
    for j = 1:length(t_ball)
        if display == 1
            x_draw = x_ball(j);
            y_draw = y_ball(j);
            h1 = plot(x_draw,y_draw,"LineStyle","none","Marker","o","Color",[1,0.5,0.3]);
            drawnow()
            pause(0.1)
            delete(h1)
        else
            xlim([0,t_end*t_step])
            x_draw = t_ball(j);
            y_draw = y_ball(j);
            h1 = plot(x_draw,y_draw,"LineStyle","none","Marker","o","Color",[1,0.5,0.3]);
            drawnow()
            pause(0.1)
        end
    end
end
end

function electron_group(t_step, omega, fai, display)
arguments
    t_step double = 0.05;
    omega double = 5;
    fai double = 4;
    display = 1;
end
% 电子集体振荡行为与光的耦合
x0 = [1:t_step:pi]'; x1 = [pi+t_step:t_step:2*pi]'; x2 = [2*pi+t_step:t_step:3*pi]';
x = [x0; x1; x2];
y_elect = [0*x0; sin(omega*x1); 0*x2];
y_light = sin(omega*x+fai);

figure()
hold on
if display == 1
    ylim([-2.1,2.1])
    h1 = plot(x,y_elect,"LineWidth",1,"Color",[0.6,0.8,1]);
else
    xlim([0,51])
    ylim([0,81])
end
light_I = sum(abs((y_elect(length(x0) + (1:length(x2))) - y_light(length(x0) + (1:length(x2))))));
for i = 1:50
    y_light = sin(omega*x+fai-0.1*i);
    light_I(i,1) = sum(abs((y_elect(length(x0) + (1:length(x2))) - ...
        y_light(length(x0) + (1:length(x2))))));
    if display == 1

        h2 = plot(x,y_light,"LineWidth",1,"Color",[1,0.5,0.3]);
        h3 = plot(x,y_light+y_elect,"LineWidth",2,"Color",[0.3,0.7,0.3]);
        drawnow()
        pause(0.1)
        delete(h2)
        delete(h3)
    else
        h4 = plot(light_I,"LineWidth",2,"Color",[0.6,0.8,1]);
        drawnow()
        pause(0.1)
    end
end
end

function CompleteProcess(display)
arguments
    display = 1;
end
% Citation：
% 表面等离子共振效应中传统近似理论与薄膜光学理论
% doi:１０．３７８８／ｇｚｘｂ２０１０３９０７．１２１６
% Warning: I have changed the parameters.
if display == 1
    lambda = 632.8 * 1e-9; % Wavelength
    epsilon0 = 1.516^2; % Prism
    epsilon1 = -8.5 + 0.7i; % Metal
    epsilon2 = 1.33^2; % Medium
    d = 50 * 1e-9; % Distance from metal to medium

    theta_test = [50:0.005:95]';
    R_theta = zeros(size(theta_test));

    for i = 1:length(theta_test)
        theta = pi/180 * theta_test(i); %incident angle（°）

        k0x = (2*pi/lambda)*sqrt(epsilon0)*sin(theta);
        k0z = sqrt((2*pi/lambda)^2*epsilon0 - k0x^2);
        k1z = sqrt((2*pi/lambda)^2*epsilon1 - k0x^2);
        k2z = sqrt((2*pi/lambda)^2*epsilon2 - k0x^2);

        r01 = (epsilon1*k0z - epsilon0*k1z) / (epsilon1*k0z + epsilon0*k1z);
        r12 = (epsilon2*k1z - epsilon1*k2z) / (epsilon2*k1z + epsilon1*k2z);

        r012 = (r01 + r12*exp(2i*k1z*d)) / (1 + r01*r12*exp(2i*k1z*d));
        R_theta(i,1) = abs(r012);
    end
    figure()
    hold on
    plot(theta_test,R_theta,"LineWidth",2,"Color",[0.3,0.5,1]);
    grid on
    xlabel("\theta(°)")
    ylabel("Reflectivity")
    ylim([0.1,1])
else
    % 未能复现
    % Warning: it may has some bugs.
    lambda = 632.8 * 1e-9;
    epsilon1 = -18 + 0.7i;
    d = 50 * 1e-9;
    N0 = 3.24;
    syms n1 k1
    eqn1 = real(epsilon1) == n1^2 - k1^2;
    eqn2 = imag(epsilon1) == 2*n1*k1;
    [n1, k1] = solve([eqn1,eqn2]);
    n1 = double(n1(1));
    k1 = double(k1(1));
    N1 = n1 - k1*1i;
    N2 = 1;

    theta_test = [34:0.01:35.4]';
    R_thetav = zeros(size(theta_test));

    for i = 1:length(theta_test)
        theta0 = pi/180 * theta_test(i);
        theta1 = asin(N0/N1*sin(theta0));
        theta2 = asin(N1/N2*sin(theta1));

        itar0 = N0/cos(theta0);
        itar1 = N1/cos(theta1);
        itar2 = N2/cos(theta2);
        sigma1 = 2*pi*N1*d*cos(theta1)/lambda;
        A = [[cos(sigma1), 1i/itar1*sin(sigma1)]; [1i*itar1*sin(sigma1), cos(sigma1)]] * [1; itar2];
        B = A(1,1); C = A(2,1);

        r = (itar0 - C/B) / (itar0 + C/B);
        R_thetav(i,1) = abs(r)^2;
    end
    figure()
    plot(theta_test,R_thetav,"LineWidth",2,"Color",[0.3,0.5,1]);
end
end
