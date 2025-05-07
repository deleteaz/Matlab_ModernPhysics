%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ellipsometry

% Author:   ZhiXiang Wei
% Email:    3011860885@qq.com
% Version:  Matlab2022b
% Config:   AMD Ryzen7 6800H && RTX3060
% Usage:    1.no need to change anything, just run.
%           2.color: str, color of light.
%           3.t_start: int, when ligth in.
%           4.t_end: int, when light out.
%           5.dt: float, step of time.
%           6.Projective_display: int, the figure displayed.
%           7.omega_x, omega_y, omega_z: flaot, angular velocity of polarized light.
%           8.ds: float, thickness of Polarizer and Waveplate.
%           9.n1, n2: complex, refractive index of medium.
%           10.theta_min, theta_max: float, the incident angle of light from min to max.
%           11.i_max: int, runtime of light.
%           12.n_max: float, Maximum refractive index of medium.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;
color = ["red","green","blue","cyan","yellow"];
t_start = 1;
t_end = 20;
dt = 0.25;
Projective_display = 2;
lightPlot1(color,t_start,t_end,dt,Projective_display);

omega_x = 1; omega_y = 1; omega_z = 1;
ds = 0.25;
lightPlot2(color,t_start,t_end,dt,omega_x,omega_y,omega_z,ds) ;

n1 = 1;
n2 = 3.882 - 1i*0.019;
theta_min = -65;
theta_max = -65;
i_max = 25;
n_max = 1;
CompleteProcess(color,t_start,t_end,dt,omega_x,omega_y,omega_z,n1,n2,ds,theta_min,theta_max,i_max,n_max);

function lightPlot1(color,t_start,t_end,dt,Projective_display)
figure()
subplot(2, 6, [1:4,7:10])
view(3)
axis equal
t = (t_start:dt:t_end)';
x0 = zeros(length(t),1);
y0 = zeros(length(t),1);
z0 = zeros(length(t),1);
x = 2*sin(t);
y = 2*sin(1.5*t);
z = 2*sin(t);
hold on
% 3D-Axes
quiver3(t(1), 0, 0, (max(t)-min(t)), 0, 0, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3], "MaxHeadSize", 0.1)
quiver3(t(1), min(y), 0, 0, (max(y)-min(y)), 0, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3])
quiver3(t(1), 0, min(z), 0, 0, (max(z)-min(z)), "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3])
% Light
plot3(t, y0, z, "Color", color(1), "LineWidth", 1, "LineStyle", "--")
quiver3(t, y0, z0, x0, y0, z, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(1), "MaxHeadSize", 0.025)
plot3(t, y, z0, "Color", color(2), "LineWidth", 1, "LineStyle", "--")
quiver3(t, y0, z0, x0, y, z0, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(2), "MaxHeadSize", 0.025)

if Projective_display > 0
    % Projective area
    rect_pos = floor(length(t)/4);
    rect_bottom = min(z); rect_top = max(z); rect_left = min(y); rect_right = max(y);
    patch([t(rect_pos), t(rect_pos), t(rect_pos), t(rect_pos)], [rect_left,rect_left,rect_right,rect_right], [rect_bottom,rect_top,rect_top,rect_bottom], [0.5,0.7,0.9], "FaceAlpha", 0.3)
    % Light projection
    quiver3(t(rect_pos), y(rect_pos), 0, 0, 0, z(rect_pos), "off", "ShowArrowHead", "on", "LineWidth", 2, "Color", color(3))
    quiver3(t(rect_pos), 0, z(rect_pos), 0, y(rect_pos), 0, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", color(3))
    quiver3(t(rect_pos), 0, 0, 0, y(rect_pos), z(rect_pos), "off", "ShowArrowHead", "on", "LineWidth", 1.5, "Color", color(3))

    if Projective_display > 1
        % Compound light
        plot3(t, y, z, "Color", color(3), "LineWidth", 1.5, "LineStyle", "-")
        quiver3(t, y0, z0, x0, y, z, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(3), "MaxHeadSize", 0.05)

        subplot(2, 6, [5,6])
        hold on
        % 2D-Axes
        quiver(rect_right, 0, (rect_left-rect_right), 0, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3])
        quiver(0, rect_bottom, 0, (rect_top-rect_bottom), "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3])
        % Light projection
        quiver(0, 0, -y(rect_pos), z(rect_pos), "off", "ShowArrowHead", "on", "LineWidth", 2, "Color", color(3))
        quiver(0, z(rect_pos), -y(rect_pos), 0, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", color(1))
        quiver(-y(rect_pos), 0, 0, z(rect_pos), "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", color(2))

        subplot(2, 6, [11,12])
        hold on
        view(90,0)
        % 3D-Axes
        quiver3(t(1), 0, 0, (max(t)-min(t)), 0, 0, "off", "ShowArrowHead", "on", "LineWidth", 2, "Color", [0.3,0.3,0.3], "MaxHeadSize", 0.1)
        quiver3(t(1), min(y), 0, 0, (max(y)-min(y)), 0, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3])
        quiver3(t(1), 0, min(z), 0, 0, (max(z)-min(z)), "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3])
        % Compound light
        plot3(t, y, z, "Color", color(3), "LineWidth", 1.5, "LineStyle", "-")
    end
end
end

function lightPlot2(color,t_start,t_end,dt,omega_x,omega_y,omega_z,ds) 
df = 1/(omega_y*2*pi);
dlambda = 1/df;
figure()
for i = 0:25
    subplot(2, 6, [1:4,7:10])
    view(3)
    axis equal
    hold on
    % Light_SP1
    t1 = (t_start:dt:floor(t_end/2))';
    x1 = zeros(length(t1),1);
    y1 = zeros(length(t1),1);
    z1 = zeros(length(t1),1);
    x_1 = 2*sin(omega_x*t1 - i*dt);
    y_1 = 2*sin(omega_y*t1 - i*dt);
    z_1 = 2*sin(omega_z*t1 - i*dt);
    hp1_1 = plot3(t1, y1, z_1, "Color", color(1), "LineWidth", 1, "LineStyle", "--");
    hp1_2 = quiver3(t1, y1, z1, x1, y1, z_1, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(1), "MaxHeadSize", 0.025);
    hs1_1 = plot3(t1, y_1, z1, "Color", color(2), "LineWidth", 1, "LineStyle", "--");
    hs1_2 = quiver3(t1, y1, z1, x1, y_1, z1, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(2), "MaxHeadSize", 0.025);
    % Compound light_SP1
    hsp1_1 = plot3(t1, y_1, z_1, "Color", color(3), "LineWidth", 1.5, "LineStyle", "-");
    hsp1_2 = quiver3(t1, y1, z1, x1, y_1, z_1, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(3), "MaxHeadSize", 0.05);
    % Waveplates
    % ds = 0.25;
    WP_pos = length(t1);
    WP_bottom = min(z_1); WP_top = max(z_1); WP_left = min(y_1); WP_right = max(y_1);
    hw_1 = patch([t1(WP_pos), t1(WP_pos), t1(WP_pos), t1(WP_pos)], [WP_left,WP_left,WP_right,WP_right], [WP_bottom,WP_top,WP_top,WP_bottom], [0.5,0.7,0.9], "FaceAlpha", 0.3);
    hw_2 = patch([t1(WP_pos)+ds, t1(WP_pos)+ds, t1(WP_pos)+ds, t1(WP_pos)+ds], [WP_left,WP_left,WP_right,WP_right], [WP_bottom,WP_top,WP_top,WP_bottom], [0.5,0.7,0.9], "FaceAlpha", 0.3);
    % Light_SP2
    t2 = (floor(t_end/2)+ds:dt:t_end)';
    x2 = zeros(length(t2),1);
    y2 = zeros(length(t2),1);
    z2 = zeros(length(t2),1);
    x_2 = 2*sin(omega_x*t2 - i*dt);
    y_2 = 2*sin(omega_y*t2 - i*dt + dlambda/4);
    z_2 = 2*sin(omega_z*t2 - i*dt);
    hp2_1 = plot3(t2, y2, z_2, "Color", color(1), "LineWidth", 1, "LineStyle", "--");
    hp2_2 = quiver3(t2, y2, z2, x2, y2, z_2, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(1), "MaxHeadSize", 0.025);
    hs2_1 = plot3(t2, y_2, z2, "Color", color(2), "LineWidth", 1, "LineStyle", "--");
    hs2_2 = quiver3(t2, y2, z2, x2, y_2, z2, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(2), "MaxHeadSize", 0.025);
    % Compound light_SP2
    hsp2_1 = plot3(t2, y_2, z_2, "Color", color(3), "LineWidth", 1.5, "LineStyle", "-");
    hsp2_2 = quiver3(t2, y2, z2, x2, y_2, z_2, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(3), "MaxHeadSize", 0.05);
    % Light_WP
    tWP = (floor(t_end/2):dt:floor(t_end/2)+ds)';
    xWP = zeros(length(tWP),1);
    yWP = zeros(length(tWP),1);
    zWP = zeros(length(tWP),1);
    y_WP = (y_1(end)-y_2(1)) ./ (tWP(1)-tWP(end)) * (tWP - tWP(1)) + y_1(end);
    z_WP = (z_1(end)-z_2(1)) ./ (tWP(1)-tWP(end)) * (tWP - tWP(1)) + z_1(end);
    hw_3 = plot3(tWP, yWP, z_WP, "Color", color(1), "LineWidth", 1, "LineStyle", "--");
    hw_4 = plot3(tWP, y_WP, zWP, "Color", color(2), "LineWidth", 1, "LineStyle", "--");
    hw_5 = plot3(tWP, y_WP, z_WP, "Color", color(3), "LineWidth", 1.5, "LineStyle", "-");
    % 3D-Axes
    h_1 = quiver3(t1(1), 0, 0, (max(t2)-min(t1)), 0, 0, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3], "MaxHeadSize", 0.1);
    h_2 = quiver3(t1(1), min(y_1), 0, 0, (max(y_1)-min(y_1)), 0, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3]);
    h_3 = quiver3(t1(1), 0, min(z_1), 0, 0, (max(z_1)-min(z_1)), "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3]);

    subplot(2, 6, [5,6])
    axis equal
    hold on
    % 2D-Axes
    h2d1_1 = quiver(WP_right, 0, (WP_left-WP_right), 0, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3]);
    h2d1_2 = quiver(0, WP_bottom, 0, (WP_top-WP_bottom), "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3]);
    % Light projection
    h2d1_3 = quiver(0, 0, -y_1(end), z_1(end), "off", "ShowArrowHead", "on", "LineWidth", 2, "Color", color(3));
    h2d1_4 = quiver(0, z_1(end), -y_1(end), 0, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", color(1));
    h2d1_5 = quiver(-y_1(end), 0, 0, z_1(end), "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", color(2));

    subplot(2, 6, [11,12])
    axis equal
    hold on
    % 2D-Axes
    h2d2_1 = quiver(WP_right, 0, (WP_left-WP_right), 0, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3]);
    h2d2_2 = quiver(0, WP_bottom, 0, (WP_top-WP_bottom), "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3]);
    % Light projection
    h2d2_3 = quiver(0, 0, -y_2(1), z_2(1), "off", "ShowArrowHead", "on", "LineWidth", 2, "Color", color(3));
    h2d2_4 = quiver(0, z_2(1), -y_2(1), 0, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", color(1));
    h2d2_5 = quiver(-y_2(1), 0, 0, z_2(1), "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", color(2));

    drawnow()
    pause(0.1)
    if i < 25
        for j = 1:3
            delete(eval(['h_',num2str(j)]));
        end
        for j = 1:5
            delete(eval(['hw_',num2str(j)]));
        end
        for j = 1:2
            delete(eval(['hs1_',num2str(j)]));
            delete(eval(['hp1_',num2str(j)]));
            delete(eval(['hsp1_',num2str(j)]));
            delete(eval(['hs2_',num2str(j)]));
            delete(eval(['hp2_',num2str(j)]));
            delete(eval(['hsp2_',num2str(j)]));
        end
        for j = 1:5
            delete(eval(['h2d1_',num2str(j)]));
            delete(eval(['h2d2_',num2str(j)]));
        end
    end
end
end

function CompleteProcess(color,t_start,t_end,dt,omega_x,omega_y,omega_z,n1,n2,ds,theta_min,theta_max,i_max,n_max) 
df = 1/(omega_y*2*pi);
dlambda = 1/df;

figure()
% theta_min = -65;
% theta_max = -65;
% i_max = 25;
% n_max = 1;
for n1 = 1:0.05:n_max
for theta = (theta_min:-1:theta_max)*pi/180
    for i = 0:i_max
        subplot(2, 6, [1:4,7:10])
        hold on
        view(45,20)
        axis equal
        % Light_SP1
        t1 = (t_start:dt:floor(t_end/2))';
        x1 = zeros(length(t1),1);
        y1 = zeros(length(t1),1);
        z1 = zeros(length(t1),1);
        x_1 = 1.5*sin(omega_x*t1 - i*dt);
        y_1 = 1.5*sin(omega_y*t1 - i*dt);
        z_1 = 1.5*sin(omega_z*t1 - i*dt);
        [Rt1, Ry1, Rz1] = Rot([t1,y1,z_1] ,theta);
        hp1_1 = plot3(Rt1, Ry1, Rz1, "Color", color(1), "LineWidth", 1, "LineStyle", "--");
        [Rt1, Ry1, Rz1, Rvx1, Rvy1, Rvz1] = Rotv([t1,y1,z1], [x1,y1,z_1] ,theta);
        hp1_2 = quiver3(Rt1, Ry1, Rz1, Rvx1, Rvy1, Rvz1, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(1), "MaxHeadSize", 0.025);
        [Rt1, Ry1, Rz1] = Rot([t1,y_1,z1] ,theta);
        hs1_1 = plot3(Rt1, Ry1, Rz1, "Color", color(2), "LineWidth", 1, "LineStyle", "--");
        [Rt1, Ry1, Rz1, Rvx1, Rvy1, Rvz1] = Rotv([t1,y1,z1], [x1,y_1,z1] ,theta);
        hs1_2 = quiver3(Rt1, Ry1, Rz1, Rvx1, Rvy1, Rvz1, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(2), "MaxHeadSize", 0.025);
        % Compound light_SP1
        [Rt1, Ry1, Rz1] = Rot([t1,y_1,z_1] ,theta);
        hsp1_1 = plot3(Rt1, Ry1, Rz1, "Color", color(3), "LineWidth", 1.5, "LineStyle", "-");
        [Rt1, Ry1, Rz1, Rvx1, Rvy1, Rvz1] = Rotv([t1,y1,z1], [x1,y_1,z_1] ,theta);
        hsp1_2 = quiver3(Rt1, Ry1, Rz1, Rvx1, Rvy1, Rvz1, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(3), "MaxHeadSize", 0.05);
        % Light_SP2
        t2 = (floor(t_end/2)+ds:dt:t_end)';
        x2 = zeros(length(t2),1);
        y2 = zeros(length(t2),1);
        z2 = zeros(length(t2),1);
        x_2 = 1.5*sin(omega_x*t2 - i*dt);
        y_2 = 1.5*sin(omega_y*t2 - i*dt + dlambda/4);
        z_2 = 1.5*sin(omega_z*t2 - i*dt);
        [Rt2, Ry2, Rz2] = Rot([t2, y2, z_2], theta);
        hp2_1 = plot3(Rt2, Ry2, Rz2, "Color", color(1), "LineWidth", 1, "LineStyle", "--");
        [Rt2, Ry2, Rz2, Rvx2, Rvy2, Rvz2] = Rotv([t2,y2,z2], [x2,y2,z_2] ,theta);
        hp2_2 = quiver3(Rt2, Ry2, Rz2, Rvx2, Rvy2, Rvz2, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(1), "MaxHeadSize", 0.025);
        [Rt2, Ry2, Rz2] = Rot([t2, y_2, z2], theta);
        hs2_1 = plot3(Rt2, Ry2, Rz2, "Color", color(2), "LineWidth", 1, "LineStyle", "--");
        [Rt2, Ry2, Rz2, Rvx2, Rvy2, Rvz2] = Rotv([t2,y2,z2], [x2,y_2,z2] ,theta);
        hs2_2 = quiver3(Rt2, Ry2, Rz2, Rvx2, Rvy2, Rvz2, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(2), "MaxHeadSize", 0.025);
        % Compound light_SP2
        [Rt2, Ry2, Rz2] = Rot([t2, y_2, z_2], theta);
        hsp2_1 = plot3(Rt2, Ry2, Rz2, "Color", color(3), "LineWidth", 1.5, "LineStyle", "-");
        [Rt2, Ry2, Rz2, Rvx2, Rvy2, Rvz2] = Rotv([t2,y2,z2], [x2,y_2,z_2] ,theta);
        hsp2_2 = quiver3(Rt2, Ry2, Rz2, Rvx2, Rvy2, Rvz2, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(3), "MaxHeadSize", 0.05);
        % Waveplates_1
        WP_pos = length(t1);
        WP_bottom = min(z_1); WP_top = max(z_1); WP_left = min(y_1); WP_right = max(y_1);
        [Rt1_1, WP_y1, WP_z1] = Rot([t1(WP_pos) , WP_left, WP_bottom], theta);
        [Rt1_2, WP_y2, WP_z2] = Rot([t1(WP_pos) , WP_left, WP_top], theta);
        [Rt1_3, WP_y3, WP_z3] = Rot([t1(WP_pos) , WP_right, WP_top], theta);
        [Rt1_4, WP_y4, WP_z4] = Rot([t1(WP_pos) , WP_right, WP_bottom], theta);
        hw1_1 = patch([Rt1_1, Rt1_2, Rt1_3, Rt1_4], [WP_y1,WP_y2,WP_y3,WP_y4], [WP_z1,WP_z2,WP_z3,WP_z4], [0.5,0.7,0.9], "FaceAlpha", 0.3);
        [Rt1_1, WP_y1, WP_z1] = Rot([t1(WP_pos)+ds , WP_left, WP_bottom], theta);
        [Rt1_2, WP_y2, WP_z2] = Rot([t1(WP_pos)+ds , WP_left, WP_top], theta);
        [Rt1_3, WP_y3, WP_z3] = Rot([t1(WP_pos)+ds , WP_right, WP_top], theta);
        [Rt1_4, WP_y4, WP_z4] = Rot([t1(WP_pos)+ds , WP_right, WP_bottom], theta);
        hw1_2 = patch([Rt1_1, Rt1_2, Rt1_3, Rt1_4], [WP_y1,WP_y2,WP_y3,WP_y4], [WP_z1,WP_z2,WP_z3,WP_z4], [0.5,0.7,0.9], "FaceAlpha", 0.3);
        % Light_WP1
        tWP1 = (floor(t_end/2):dt:floor(t_end/2)+ds)';
        yWP1 = zeros(length(tWP1),1);
        zWP1 = zeros(length(tWP1),1);
        y_WP1 = (y_1(end)-y_2(1)) ./ (tWP1(1)-tWP1(end)) * (tWP1 - tWP1(1)) + y_1(end);
        z_WP1 = (z_1(end)-z_2(1)) ./ (tWP1(1)-tWP1(end)) * (tWP1 - tWP1(1)) + z_1(end);
        [RtWP1, RyWP1, RzWP1] = Rot([tWP1,yWP1,z_WP1] ,theta);
        hw1_3 = plot3(RtWP1, RyWP1, RzWP1, "Color", color(1), "LineWidth", 1, "LineStyle", "--");
        [RtWP1, RyWP1, RzWP1] = Rot([tWP1,y_WP1,zWP1] ,theta);
        hw1_4 = plot3(RtWP1, RyWP1, RzWP1, "Color", color(2), "LineWidth", 1, "LineStyle", "--");
        [RtWP1, RyWP1, RzWP1] = Rot([tWP1,y_WP1,z_WP1] ,theta);
        hw1_5 = plot3(RtWP1, RyWP1, RzWP1, "Color", color(3), "LineWidth", 1.5, "LineStyle", "-");
        % 3D-Axes1
        [Axes_x, Axes_y, Axes_z, Axes_vx, Axes_vy, Axes_vz] = Rotv([t1(1), 0, 0], [t2(end)-t1(1), 0, 0],theta);
        h_1 = quiver3(Axes_x, Axes_y, Axes_z, Axes_vx, Axes_vy, Axes_vz, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3], "MaxHeadSize", 0.1);
        [Axes_x, Axes_y, Axes_z, Axes_vx, Axes_vy, Axes_vz] = Rotv([t1(1), min(y_1), 0], [0, (max(y_1)-min(y_1)), 0] ,theta);
        h_2 = quiver3(Axes_x, Axes_y, Axes_z, Axes_vx, Axes_vy, Axes_vz, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3], "MaxHeadSize", 0.1);
        [Axes_x, Axes_y, Axes_z, Axes_vx, Axes_vy, Axes_vz] = Rotv([t1(1), 0, min(z_1)], [0, 0, (max(z_1)-min(z_1))] ,theta);
        h_3 = quiver3(Axes_x, Axes_y, Axes_z, Axes_vx, Axes_vy, Axes_vz, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3], "MaxHeadSize", 0.1);
        % Refractivity
        theta_r = -theta;
        theta_t = asin(n1/n2*sin(theta));
        rrp = (n2*cos(theta)-n1*cos(theta_t)) / (n2*cos(theta)+n1*cos(theta_t));
        rrp_abs = abs(rrp); rrp_ang = angle(rrp);
        rrs = (n1*cos(theta)-n2*cos(theta_t)) / (n1*cos(theta)+n2*cos(theta_t));
        rrs_abs = abs(rrs); rrs_ang = angle(rrs);
        rtp = 2*n1*cos(theta) / (n2*cos(theta)+n1*cos(theta_t));
        rtp_abs = abs(rtp); rtp_ang = angle(rtp);
        rts = 2*n1*cos(theta) / (n1*cos(theta)+n2*cos(theta_t));
        rts_abs = abs(rts); rts_ang = angle(rts);
        theta_t = real(theta_t);
        % Light_SP3
        t3 = (t_end:dt:3*floor(t_end/2)-ds)';
        x3 = zeros(length(t3),1);
        y3 = zeros(length(t3),1);
        z3 = zeros(length(t3),1);
        x_3 = 1.5*sin(omega_x*t3 - i*dt);
        y_3 = rrs_abs*sin(omega_y*t3 - i*dt + rrs_ang);
        z_3 = rrp_abs*sin(omega_z*t3 - i*dt + rrp_ang);
        [Rt3, Ry3, Rz3] = Rot([t3, y3, z_3], theta_r);
        hp3_1 = plot3(Rt3, Ry3, Rz3+2*Rz2(end), "Color", color(1), "LineWidth", 1, "LineStyle", "--");
        [Rt3, Ry3, Rz3, Rvx3, Rvy3, Rvz3] = Rotv([t3,y3,z3], [x3,y3,z_3] ,theta_r);
        hp3_2 = quiver3(Rt3, Ry3, Rz3+2*Rz2(end), Rvx3, Rvy3, Rvz3, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(1), "MaxHeadSize", 0.025);
        [Rt3, Ry3, Rz3] = Rot([t3, y_3, z3], theta_r);
        hs3_1 = plot3(Rt3, Ry3, Rz3+2*Rz2(end), "Color", color(2), "LineWidth", 1, "LineStyle", "--");
        [Rt3, Ry3, Rz3, Rvx3, Rvy3, Rvz3] = Rotv([t3,y3,z3], [x3,y_3,z3] ,theta_r);
        hs3_2 = quiver3(Rt3, Ry3, Rz3+2*Rz2(end), Rvx3, Rvy3, Rvz3, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(2), "MaxHeadSize", 0.025);
        % Compound light_SP3
        [Rt3, Ry3, Rz3] = Rot([t3, y_3, z_3], theta_r);
        hsp3_1 = plot3(Rt3, Ry3, Rz3+2*Rz2(end), "Color", color(3), "LineWidth", 1.5, "LineStyle", "-");
        [Rt3, Ry3, Rz3, Rvx3, Rvy3, Rvz3] = Rotv([t3,y3,z3], [x3,y_3,z_3] ,theta_r);
        hsp3_2 = quiver3(Rt3, Ry3, Rz3+2*Rz2(end), Rvx3, Rvy3, Rvz3, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(3), "MaxHeadSize", 0.05);
        % Light_SP4
        t4 = (3*floor(t_end/2):dt:2*t_end)';
        x4 = zeros(length(t4),1);
        y4 = zeros(length(t4),1);
        z4 = zeros(length(t4),1);
        x_4 = 1.5*sin(omega_x*t4 - i*dt);
        y_4 = rrs_abs*sin(omega_y*t4 - i*dt - rrs_ang + dlambda/4);
        z_4 = rrp_abs*sin(omega_z*t4 - i*dt - rrp_ang);
        [Rt4, Ry4, Rz4] = Rot([t4, y4, z_4], theta_r);
        hp4_1 = plot3(Rt4, Ry4, Rz4+2*Rz2(end), "Color", color(1), "LineWidth", 1, "LineStyle", "--");
        [Rt4, Ry4, Rz4, Rvx4, Rvy4, Rvz4] = Rotv([t4,y4,z4], [x4,y4,z_4] ,theta_r);
        hp4_2 = quiver3(Rt4, Ry4, Rz4+2*Rz2(end), Rvx4, Rvy4, Rvz4, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(1), "MaxHeadSize", 0.025);
        [Rt4, Ry4, Rz4] = Rot([t4, y_4, z4], theta_r);
        hs4_1 = plot3(Rt4, Ry4, Rz4+2*Rz2(end), "Color", color(2), "LineWidth", 1, "LineStyle", "--");
        [Rt4, Ry4, Rz4, Rvx4, Rvy4, Rvz4] = Rotv([t4,y4,z4], [x4,y_4,z4] ,theta_r);
        hs4_2 = quiver3(Rt4, Ry4, Rz4+2*Rz2(end), Rvx4, Rvy4, Rvz4, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(2), "MaxHeadSize", 0.025);
        % Compound light_SP4
        [Rt4, Ry4, Rz4] = Rot([t4, y_4, z_4], theta_r);
        hsp4_1 = plot3(Rt4, Ry4, Rz4+2*Rz2(end), "Color", color(3), "LineWidth", 1.5, "LineStyle", "-");
        [Rt4, Ry4, Rz4, Rvx4, Rvy4, Rvz4] = Rotv([t4,y4,z4], [x4,y_4,z_4] ,theta_r);
        hsp4_2 = quiver3(Rt4, Ry4, Rz4+2*Rz2(end), Rvx4, Rvy4, Rvz4, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(3), "MaxHeadSize", 0.05);
        % % Waveplates_2
        WP_pos = length(t3);
        WP_bottom = min(z_1); WP_top = max(z_1); WP_left = min(y_1); WP_right = max(y_1);
        [Rt3_1, WP_y1, WP_z1] = Rot([t3(WP_pos) , WP_left, WP_bottom], theta_r);
        [Rt3_2, WP_y2, WP_z2] = Rot([t3(WP_pos) , WP_left, WP_top], theta_r);
        [Rt3_3, WP_y3, WP_z3] = Rot([t3(WP_pos) , WP_right, WP_top], theta_r);
        [Rt3_4, WP_y4, WP_z4] = Rot([t3(WP_pos) , WP_right, WP_bottom], theta_r);
        hw2_1 = patch([Rt3_1, Rt3_2, Rt3_3, Rt3_4], [WP_y1,WP_y2,WP_y3,WP_y4], [WP_z1,WP_z2,WP_z3,WP_z4]+2*Rz2(end), [0.5,0.7,0.9], "FaceAlpha", 0.3);
        [Rt3_1, WP_y1, WP_z1] = Rot([t3(WP_pos)+ds , WP_left, WP_bottom], theta_r);
        [Rt3_2, WP_y2, WP_z2] = Rot([t3(WP_pos)+ds , WP_left, WP_top], theta_r);
        [Rt3_3, WP_y3, WP_z3] = Rot([t3(WP_pos)+ds , WP_right, WP_top], theta_r);
        [Rt3_4, WP_y4, WP_z4] = Rot([t3(WP_pos)+ds , WP_right, WP_bottom], theta_r);
        hw2_2 = patch([Rt3_1, Rt3_2, Rt3_3, Rt3_4], [WP_y1,WP_y2,WP_y3,WP_y4], [WP_z1,WP_z2,WP_z3,WP_z4]+2*Rz2(end), [0.5,0.7,0.9], "FaceAlpha", 0.3);
        % Light_WP2
        tWP2 = (3*floor(t_end/2)-ds:dt:3*floor(t_end/2))';
        yWP2 = zeros(length(tWP2),1);
        zWP2 = zeros(length(tWP2),1);
        y_WP2 = (y_3(end)-y_4(1)) ./ (tWP2(1)-tWP2(end)) * (tWP2 - tWP2(1)) + y_3(end);
        z_WP2 = (z_3(end)-z_4(1)) ./ (tWP2(1)-tWP2(end)) * (tWP2 - tWP2(1)) + z_3(end);
        [RtWP2, RyWP2, RzWP2] = Rot([tWP2,yWP2,z_WP2] ,theta_r);
        hw2_3 = plot3(RtWP2, RyWP2, RzWP2+2*Rz2(end), "Color", color(1), "LineWidth", 1, "LineStyle", "--");
        [RtWP2, RyWP2, RzWP2] = Rot([tWP2,y_WP2,zWP2] ,theta_r);
        hw2_4 = plot3(RtWP2, RyWP2, RzWP2+2*Rz2(end), "Color", color(2), "LineWidth", 1, "LineStyle", "--");
        [RtWP2, RyWP2, RzWP2] = Rot([tWP2,y_WP2,z_WP2] ,theta_r);
        hw2_5 = plot3(RtWP2, RyWP2, RzWP2+2*Rz2(end), "Color", color(3), "LineWidth", 1.5, "LineStyle", "-");
        % 3D-Axes-r
        [Axes_x, Axes_y, Axes_z, Axes_vx, Axes_vy, Axes_vz] = Rotv([t3(1), 0, 0], [t4(end)-t3(1), 0, 0],theta_r);
        hr_1 = quiver3(Axes_x, Axes_y, Axes_z+2*Rz2(end), Axes_vx, Axes_vy, Axes_vz, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3], "MaxHeadSize", 0.1);
        [Axes_x, Axes_y, Axes_z, Axes_vx, Axes_vy, Axes_vz] = Rotv([t4(end), min(y_4), 0], [0, (max(y_4)-min(y_4)), 0] ,theta_r);
        hr_2 = quiver3(Axes_x, Axes_y, Axes_z+2*Rz2(end), Axes_vx, Axes_vy, Axes_vz, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3], "MaxHeadSize", 0.1);
        [Axes_x, Axes_y, Axes_z, Axes_vx, Axes_vy, Axes_vz] = Rotv([t4(end), 0, min(z_4)], [0, 0, (max(z_4)-min(z_4))] ,theta_r);
        hr_3 = quiver3(Axes_x, Axes_y, Axes_z+2*Rz2(end), Axes_vx, Axes_vy, Axes_vz, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3], "MaxHeadSize", 0.1);
        % Light_SP5
        t5 = (0:dt:floor(t_end/2))';
        x5 = zeros(length(t5),1);
        y5 = zeros(length(t5),1);
        z5 = zeros(length(t5),1);
        x_5 = 1.5*sin(omega_x*t5 - i*dt);
        y_5 = rts_abs*sin(omega_y*t5 - i*dt - rts_ang);
        z_5 = rtp_abs*sin(omega_z*t5 - i*dt - rtp_ang);
        [Rt5, Ry5, Rz5] = Rot([t5, y5, z_5], theta_t);
        hp5_1 = plot3(Rt5+Rt2(end), Ry5, Rz5+Rz2(end), "Color", color(1), "LineWidth", 1, "LineStyle", "--");
        [Rt5, Ry5, Rz5, Rvx5, Rvy5, Rvz5] = Rotv([t5,y5,z5], [x5,y5,z_5] ,theta_t);
        hp5_2 = quiver3(Rt5+Rt2(end), Ry5, Rz5+Rz2(end), Rvx5, Rvy5, Rvz5, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(1), "MaxHeadSize", 0.025);
        [Rt5, Ry5, Rz5] = Rot([t5, y_5, z5], theta_t);
        hs5_1 = plot3(Rt5+Rt2(end), Ry5, Rz5+Rz2(end), "Color", color(2), "LineWidth", 1, "LineStyle", "--");
        [Rt5, Ry5, Rz5, Rvx5, Rvy5, Rvz5] = Rotv([t5,y5,z5], [x5,y_5,z5] ,theta_t);
        hs5_2 = quiver3(Rt5+Rt2(end), Ry5, Rz5+Rz2(end), Rvx5, Rvy5, Rvz5, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(2), "MaxHeadSize", 0.025);
        % Compound light_SP53
        [Rt5, Ry5, Rz5] = Rot([t5, y_5, z_5], theta_t);
        hsp5_1 = plot3(Rt5+Rt2(end), Ry5, Rz5+Rz2(end), "Color", color(3), "LineWidth", 1.5, "LineStyle", "-");
        [Rt5, Ry5, Rz5, Rvx5, Rvy5, Rvz5] = Rotv([t5,y5,z5], [x5,y_5,z_5] ,theta_t);
        hsp5_2 = quiver3(Rt5+Rt2(end), Ry5, Rz5+Rz2(end), Rvx5, Rvy5, Rvz5, "off", "ShowArrowHead", "on", "LineWidth", 0.5, "Color", color(3), "MaxHeadSize", 0.05);
        % 3D-Axes-t
        [Axes_x, Axes_y, Axes_z, Axes_vx, Axes_vy, Axes_vz] = Rotv([t5(1), 0, 0], [t5(end)-t5(1), 0, 0],theta_t);
        ht_1 = quiver3(Axes_x+Rt2(end), Axes_y, Axes_z+Rz2(end), Axes_vx, Axes_vy, Axes_vz, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3], "MaxHeadSize", 0.1);
        [Axes_x, Axes_y, Axes_z, Axes_vx, Axes_vy, Axes_vz] = Rotv([t5(end), min(y_5), 0], [0, (max(y_5)-min(y_5)), 0] ,theta_t);
        ht_2 = quiver3(Axes_x+Rt2(end), Axes_y, Axes_z+Rz2(end), Axes_vx, Axes_vy, Axes_vz, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3], "MaxHeadSize", 0.1);
        [Axes_x, Axes_y, Axes_z, Axes_vx, Axes_vy, Axes_vz] = Rotv([t5(end), 0, min(z_5)], [0, 0, (max(z_5)-min(z_5))] ,theta_t);
        ht_3 = quiver3(Axes_x+Rt2(end), Axes_y, Axes_z+Rz2(end), Axes_vx, Axes_vy, Axes_vz, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3], "MaxHeadSize", 0.1);
        % Platform
        hP_1 = patch([Rt1(1), Rt4(end), Rt4(end), Rt1(1)], [min(y_1), min(y_1), max(y_1), max(y_1)], Rz2(end)*ones(1,4), [0.9,0.9,0.1], "FaceAlpha", 0.1);
        hP_2 = patch([Rt2(end), Rt2(end), Rt2(end), Rt2(end)], [min(y_1), max(y_1), max(y_1), min(y_1)], [Rz1(1), Rz1(1), Rz5(end)+Rz2(end), Rz5(end)+Rz2(end)], [0.3,0.9,0.3], "FaceAlpha", 0.1);

        subplot(2, 6, [5,6])
        axis equal
        hold on
        % 2D-Axes
        h2d1_1 = quiver(WP_right, 0, (WP_left-WP_right), 0, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3]);
        h2d1_2 = quiver(0, WP_bottom, 0, (WP_top-WP_bottom), "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3]);
        % Light projection
        h2d1_3 = quiver(0, 0, -y_1(1), z_1(1), "off", "ShowArrowHead", "on", "LineWidth", 2, "Color", color(3));
        h2d1_4 = quiver(0, z_1(1), -y_1(1), 0, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", color(1));
        h2d1_5 = quiver(-y_1(1), 0, 0, z_1(1), "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", color(2));

        subplot(2, 6, [11,12])
        axis equal
        hold on
        % 2D-Axes
        h2d2_1 = quiver(WP_right, 0, (WP_left-WP_right), 0, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3]);
        h2d2_2 = quiver(0, WP_bottom, 0, (WP_top-WP_bottom), "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", [0.3,0.3,0.3]);
        % Light projection
        h2d2_3 = quiver(0, 0, -y_4(end), z_4(end), "off", "ShowArrowHead", "on", "LineWidth", 2, "Color", color(3));
        h2d2_4 = quiver(0, z_4(end), -y_4(end), 0, "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", color(1));
        h2d2_5 = quiver(-y_4(end), 0, 0, z_4(end), "off", "ShowArrowHead", "on", "LineWidth", 1, "Color", color(2));

        drawnow()
        pause(0.1)

        for j = 1:3
            delete(eval(['h_',num2str(j)]));
            delete(eval(['hr_',num2str(j)]));
            delete(eval(['ht_',num2str(j)]));
        end
        for j = 1:5
            delete(eval(['hw1_',num2str(j)]));
            delete(eval(['hw2_',num2str(j)]));
            delete(eval(['h2d1_',num2str(j)]));
            delete(eval(['h2d2_',num2str(j)]));
        end
        for j = 1:2
            delete(eval(['hs1_',num2str(j)]));
            delete(eval(['hp1_',num2str(j)]));
            delete(eval(['hsp1_',num2str(j)]));
            delete(eval(['hs2_',num2str(j)]));
            delete(eval(['hp2_',num2str(j)]));
            delete(eval(['hsp2_',num2str(j)]));
            delete(eval(['hs3_',num2str(j)]));
            delete(eval(['hp3_',num2str(j)]));
            delete(eval(['hsp3_',num2str(j)]));
            delete(eval(['hs4_',num2str(j)]));
            delete(eval(['hp4_',num2str(j)]));
            delete(eval(['hsp4_',num2str(j)]));
            delete(eval(['hs5_',num2str(j)]));
            delete(eval(['hp5_',num2str(j)]));
            delete(eval(['hsp5_',num2str(j)]));
            delete(eval(['hP_',num2str(j)]));
        end

    end
end
end
end

function [R_x, R_y, R_z] = Rot(M, theta)
% M = [row,row,row]
R = [[cos(theta) 0 sin(theta)]; [0 1 0]; [-sin(theta) 0 cos(theta)]];
trans = M*R;
R_x = trans(:,1);
R_y = trans(:,2);
R_z = trans(:,3);
end

function [R_x, R_y, R_z, R_vx, R_vy, R_vz] = Rotv(Mx, Mv, theta)
% Mx = Mv = [row,row,row]
R = [[cos(theta) 0 sin(theta)]; [0 1 0]; [-sin(theta) 0 cos(theta)]];
trans = Mx*R;
R_x = trans(:,1);
R_y = trans(:,2);
R_z = trans(:,3);
trans = Mv*R;
R_vx = trans(:,1);
R_vy = trans(:,2);
R_vz = trans(:,3);
end
