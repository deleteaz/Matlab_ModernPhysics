%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electron Diffraction
%
% Author:   ZhiXiang Wei
% Email:    3011860885@qq.com
% Version:  Matlab2022b
% Config:   AMD Ryzen7 6800H && RTX3060
% Usage:    1.no need to change anything, just run.
%           2.crystal_size: int, Size of crystal.
%           3.lattice_dx, lattice_dy: float, Step of the distance between lattices.
%           4.theta_max, dtheta: float, The angle and the precision for researching.
%           5.omega: float, Angular frequency.
%           6.A: float, Amplitude.
%           7.type: int, 2D(2) and 1D(1).
%           8.midpoint: int, on(1) and off(0) the midpoint.
%           9.N: int, Division accuracy.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;
crystal_size = 3;
lattice_dx = 1;
lattice_dy = 1;
theta_max = 75;
dtheta = 0.5;
omega = 19.05;
A = 0.25;
type = 2;
midpoint = 1;
N = 1500;
[theta_range, total, total_m] = CompleteProcess(crystal_size, lattice_dx, lattice_dy, theta_max, dtheta, omega, A, type, midpoint, N);
XRDplot(theta_range, theta_max, total, total_m, midpoint);

function [theta_range, total, total_m] = CompleteProcess(crystal_size, lattice_dx, lattice_dy, theta_max, dtheta, omega, A, type, midpoint, N)
arguments
crystal_size int32 = 3;
lattice_dx double = 1;
lattice_dy double = 1;
theta_max double = 75;
dtheta double = 0.5;
omega double = 19.05;
A double = 0.25;
type int32 = 2;
midpoint int32 = 1;
N int32 = 1500;
end
figure()
hold on
color = ["red","green","blue","yellow","magenta"];
color_atom = [0.5,0.5,0.5]; linestyle_atom = "none"; linewidth_atom = 2; Marker_atom = "o"; MarekerSize_atom = 20; MarkerFaceColor_atom = color;
axis equal
crystal_size = double(crystal_size);
xlim([0,2*crystal_size])
ylim([0,2*crystal_size])
x1 = (1:lattice_dx:crystal_size)';
y1 = (1:lattice_dy:crystal_size)';
switch type
    case 2
        x1 = (1:lattice_dx:crystal_size)';
        y1 = (1:lattice_dy:crystal_size)';
        len_x1 = length(x1); len_y1 = length(y1);
        [X1, Y1] = meshgrid(x1, y1);
        R_x1 = zeros(100,len_x1*len_y1); R_x2 = zeros(N,len_x1*len_y1);
        R_y1 = zeros(100,len_x1*len_y1); R_y2 = zeros(N,len_x1*len_y1);
        aR_x2 = zeros(N, len_x1*len_y1);
        aR_y2 = zeros(N, len_x1*len_y1);
        endpointx = [0;0];
        endpointy = [0;0];
        total = zeros(1, len_x1*len_y1); total_m = 0;

        plot(X1, Y1, "LineStyle",linestyle_atom,"Marker",Marker_atom,"MarkerSize",MarekerSize_atom,"MarkerFaceColor",MarkerFaceColor_atom(3));
        h1 = plot(R_x1,R_y1,"Color",color(1));
        h2 = plot(R_x2,R_y2,"Color",color(2));
        % h3 = plot(aR_x2,aR_y2,"Color",color(5));
        he = plot(endpointx, endpointy, "Color",color(5));
        if midpoint > 0
            x1_m = (1.5:lattice_dx:crystal_size-0.5)';
            y1_m = (1.5:lattice_dy:crystal_size-0.5)';
            len_x1_m = length(x1_m); len_y1_m = length(y1_m);
            [X1_m, Y1_m] = meshgrid(x1_m, y1_m);
            R_x1_m = zeros(100,len_x1_m*len_y1_m); R_x2_m = zeros(N,len_x1_m*len_y1_m);
            R_y1_m = zeros(100,len_x1_m*len_y1_m); R_y2_m = zeros(N,len_x1_m*len_y1_m);
            aR_x2_m = zeros(N, len_x1_m*len_y1_m);
            aR_y2_m = zeros(N, len_x1_m*len_y1_m);
            total_m = zeros(1, len_x1_m*len_y1_m);
            plot(X1_m, Y1_m, "LineStyle",linestyle_atom,"Marker",Marker_atom,"MarkerSize",MarekerSize_atom,"MarkerFaceColor",MarkerFaceColor_atom(1));
            h1_m = plot(R_x1_m,R_y1_m,"Color",color(3));
            h2_m = plot(R_x2_m,R_y2_m,"Color",color(4));
            % h3_m = plot(aR_x2_m,aR_y2_m,"Color",color(5));
        end

        theta_range = (0:-dtheta:-theta_max) * pi/180;
        save = 0;
        for theta = theta_range
            save = save+1;
            id = 0;
            for i = 1:len_x1
                for j = 1:len_y1
                    id = id + 1;
                    t1 = linspace(0, x1(i), 100)';
                    yc1 = A*sin(omega*t1) + y1(j);
                    [R_x1(:,id), R_y1(:,id)] = Rot([t1, yc1, ones(size(t1))], theta, [x1(i),y1(j)]);
                    t2 = linspace(x1(i), (max(x1) + crystal_size - (j-1)*lattice_dy*sin(-theta)) - (crystal_size-i)*lattice_dx*2e-4*(rad2deg(-theta))^1.905, N)';
                    yc2 = A*sin(omega*t2) + y1(j);
                    [R_x2(:,id), R_y2(:,id)] = Rot([t2, yc2, ones(size(t2))], -theta, [x1(i),y1(j)]);
                end
            end
            id = 0;
            for i = 1:len_x1
                for j = 1:len_y1
                    id = id + 1;
                    [aR_x2(:,id), aR_y2(:,id)] = Rot([R_x2(:,id),R_y2(:,id),ones(size(R_x2(:,id)))], theta, [x1(i),y1(j)]);
                    aR_y2(:,id) = aR_y2(:,id) - y1(j);
                end
            end
            total(:,save) = sum(aR_y2(end,:));
            [endpointx,endpointy] = Rot([(2*crystal_size) *ones(2,1), [20;-20], ones(2,1)], -theta, [x1(end),y1(1)]);
            for i = 1:len_x1*len_y1
                set(h1(i), "XData",R_x1(:,i), "YData",R_y1(:,i));
                set(h2(i), "XData",R_x2(:,i), "YData",R_y2(:,i));
                % set(h3(i), "XData",aR_x2(:,i), "YData",aR_y2(:,i));
            end
            set(he, "XData",endpointx, "YData",endpointy);

            if midpoint > 0
                id = 0;
                for i = 1:len_x1_m
                    for j = 1:len_y1_m
                        id = id + 1;
                        t1 = linspace(0, x1_m(i), 100)';
                        yc1 = A*sin(omega*t1) + y1_m(j);
                        [R_x1_m(:,id), R_y1_m(:,id)] = Rot([t1, yc1, ones(size(t1))], theta, [x1_m(i),y1_m(j)]);
                        t2 = linspace(x1_m(i), (max(x1_m) + crystal_size + 0.5 - (j-0.5)*lattice_dy*sin(-theta)) - (crystal_size - 0.5 -i)*lattice_dx*2e-4*(rad2deg(-theta))^1.905, N)';
                        yc2 = A*sin(omega*t2) + y1_m(j);
                        [R_x2_m(:,id), R_y2_m(:,id)] = Rot([t2, yc2, ones(size(t2))], -theta, [x1_m(i),y1_m(j)]);
                    end
                end
                id = 0;
                for i = 1:len_x1_m
                    for j = 1:len_y1_m
                        id = id + 1;
                        [aR_x2_m(:,id), aR_y2_m(:,id)] = Rot([R_x2_m(:,id),R_y2_m(:,id),ones(size(R_x2_m(:,id)))], theta, [x1_m(i),y1_m(j)]);
                        aR_y2_m(:,id) = aR_y2_m(:,id) - y1_m(j);
                    end
                end
                total_m(:,save) = sum(aR_y2_m(end,:));
                for i = 1:len_x1_m*len_y1_m
                    set(h1_m(i), "XData",R_x1_m(:,i), "YData",R_y1_m(:,i));
                    set(h2_m(i), "XData",R_x2_m(:,i), "YData",R_y2_m(:,i));
                    % set(h3_m(i), "XData",aR_x2_m(:,i), "YData",aR_y2_m(:,i));
                end
            end

            drawnow
        end
        drawnow

    case 1
        x1 = (crystal_size:lattice_dx:crystal_size)';
        y1 = (1:lattice_dy:crystal_size)';
        len_x1 = length(x1); len_y1 = length(y1);
        [X1, Y1] = meshgrid(x1, y1);
        R_x1 = zeros(100,len_x1*len_y1); R_x2 = zeros(N,len_x1*len_y1);
        R_y1 = zeros(100,len_x1*len_y1); R_y2 = zeros(N,len_x1*len_y1);
        aR_x2 = zeros(N, len_x1*len_y1);
        aR_y2 = zeros(N, len_x1*len_y1);
        endpointx = [0;0];
        endpointy = [0;0];
        total = zeros(1, len_x1*len_y1); total_m = 0;
        
        plot(X1, Y1, "LineStyle",linestyle_atom,"Marker",Marker_atom,"MarkerSize",MarekerSize_atom,"MarkerFaceColor",MarkerFaceColor_atom(3));
        h1 = plot(R_x1,R_y1,"Color",color(1));
        h2 = plot(R_x2,R_y2,"Color",color(2));
        % h3 = plot(aR_x2,aR_y2,"Color",color(5));
        he = plot(endpointx, endpointy, "Color",color(5));
        if midpoint > 0
            x1_m = (crystal_size:lattice_dx:crystal_size)';
            y1_m = (1.5:lattice_dy:crystal_size-0.5)';
            len_x1_m = length(x1_m); len_y1_m = length(y1_m);
            [X1_m, Y1_m] = meshgrid(x1_m, y1_m);
            R_x1_m = zeros(100,len_x1_m*len_y1_m); R_x2_m = zeros(N,len_x1_m*len_y1_m);
            R_y1_m = zeros(100,len_x1_m*len_y1_m); R_y2_m = zeros(N,len_x1_m*len_y1_m);
            aR_x2_m = zeros(N, len_x1_m*len_y1_m);
            aR_y2_m = zeros(N, len_x1_m*len_y1_m);
            total_m = zeros(1, len_x1_m*len_y1_m);
            plot(X1_m, Y1_m, "LineStyle",linestyle_atom,"Marker",Marker_atom,"MarkerSize",MarekerSize_atom,"MarkerFaceColor",MarkerFaceColor_atom(2));
            h1_m = plot(R_x1_m,R_y1_m,"Color",color(3));
            h2_m = plot(R_x2_m,R_y2_m,"Color",color(4));
            % h3_m = plot(aR_x2_m,aR_y2_m,"Color",color(5));
        end

        theta_range = (0:-dtheta:-theta_max) * pi/180;
        save = 0;
        for theta = theta_range
            save = save+1;
            id = 0;
            for i = 1:len_x1
                for j = 1:len_y1
                    id = id + 1;
                    t1 = linspace(0, x1(i), 100)';
                    yc1 = A*sin(omega*t1) + y1(j);
                    [R_x1(:,id), R_y1(:,id)] = Rot([t1, yc1, ones(size(t1))], theta, [x1(i),y1(j)]);
                    t2 = linspace(x1(i), (max(x1) + crystal_size - (j-1)*lattice_dy*sin(-theta)), N)';
                    yc2 = A*sin(omega*t2) + y1(j);
                    [R_x2(:,id), R_y2(:,id)] = Rot([t2, yc2, ones(size(t2))], -theta, [x1(i),y1(j)]);
                end
            end
            id = 0;
            for i = 1:len_x1
                for j = 1:len_y1
                    id = id + 1;
                    [aR_x2(:,id), aR_y2(:,id)] = Rot([R_x2(:,id),R_y2(:,id),ones(size(R_x2(:,id)))], theta, [x1(i),y1(j)]);
                    aR_y2(:,id) = aR_y2(:,id) - y1(j);
                end
            end
            total(:,save) = sum(aR_y2(end,:));
            [endpointx,endpointy] = Rot([(2*crystal_size) *ones(2,1), [20;-20], ones(2,1)], -theta, [x1(end),y1(1)]);
            for i = 1:len_x1*len_y1
                set(h1(i), "XData",R_x1(:,i), "YData",R_y1(:,i));
                set(h2(i), "XData",R_x2(:,i), "YData",R_y2(:,i));
                % set(h3(i), "XData",aR_x2(:,i), "YData",aR_y2(:,i));
            end
            set(he, "XData",endpointx, "YData",endpointy);

            if midpoint > 0
                id = 0;
                for i = 1:len_x1_m
                    for j = 1:len_y1_m
                        id = id + 1;
                        t1 = linspace(0, x1_m(i), 100)';
                        yc1 = A*sin(omega*t1) + y1_m(j);
                        [R_x1_m(:,id), R_y1_m(:,id)] = Rot([t1, yc1, ones(size(t1))], theta, [x1_m(i),y1_m(j)]);
                        t2 = linspace(x1_m(i), (max(x1_m) + crystal_size - (j-0.5)*lattice_dy*sin(-theta)), N)';
                        yc2 = A*sin(omega*t2) + y1_m(j);
                        [R_x2_m(:,id), R_y2_m(:,id)] = Rot([t2, yc2, ones(size(t2))], -theta, [x1_m(i),y1_m(j)]);
                    end
                end
                id = 0;
                for i = 1:len_x1_m
                    for j = 1:len_y1_m
                        id = id + 1;
                        [aR_x2_m(:,id), aR_y2_m(:,id)] = Rot([R_x2_m(:,id),R_y2_m(:,id),ones(size(R_x2_m(:,id)))], theta, [x1_m(i),y1_m(j)]);
                        aR_y2_m(:,id) = aR_y2_m(:,id) - y1_m(j);
                    end
                end
                total_m(:,save) = sum(aR_y2_m(end,:));
                for i = 1:len_x1_m*len_y1_m
                    set(h1_m(i), "XData",R_x1_m(:,i), "YData",R_y1_m(:,i));
                    set(h2_m(i), "XData",R_x2_m(:,i), "YData",R_y2_m(:,i));
                    % set(h3_m(i), "XData",aR_x2_m(:,i), "YData",aR_y2_m(:,i));
                end
            end

            drawnow
        end
        drawnow
end
end

function XRDplot(theta_range, theta_max, total, total_m, midpoint)
arguments
    theta_range (1,:)
    theta_max double
    total (1,:)
    total_m (1,:)
    midpoint int32 = 0;
end
% figure()
% hold on
% maxIndices = islocalmax(total.^2,"MinProminence",0.01);
% for i = 1:5:(length(total)-1)
%     plot(2*abs(rad2deg(theta_range(i:i+5))), (total(i:i+5)).^2,"LineWidth",2,"Color",color(1))
%     plot(2*abs(rad2deg(theta_range(i:i+5))).*maxIndices(i:i+5), total(i:i+5).^2.*maxIndices(i:i+5),"Color",[1,0.3,0.1],"Marker","|","Markersize",12,"LineStyle","none","LineWidth",1.5)
%     xticks(0:20:2*theta_max)
%     xlim([0,2*theta_max])
%     xlabel("角度2\theta/(°)")
%     ylabel("相对强度I/(f)")
%     set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
%     set(gca,"GridLineStyle",":","GridLineWidth",1.5,"TickDir","out")
%     drawnow()
%     pause(0.1)
% end
figure()
hold on
plot(2*abs(rad2deg(theta_range)), total.^2,"LineWidth",2)
maxIndices = islocalmax(total.^2,"MinProminence",0.01);
plot(2*abs(rad2deg(theta_range)).*maxIndices, total.^2.*maxIndices,"Color",[1,0.3,0.1],"Marker","|","Markersize",12,"LineStyle","none","LineWidth",1.5)
xticks(0:20:2*theta_max)
xlim([0,2*theta_max])
xlabel("角度2\theta/(°)")
ylabel("相对强度I/(f)")
title("基本晶格框架衍射谱")
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5,"TickDir","out")
if midpoint > 0
    figure()
    hold on
    plot(2*abs(rad2deg(theta_range)), (total_m).^2, "LineWidth",2)
    maxIndices = islocalmax(total_m.^2,"MinProminence",0.01);
    plot(2*abs(rad2deg(theta_range)).*maxIndices, total_m.^2.*maxIndices,"Color",[1,0.3,0.1],"Marker","|","Markersize",12,"LineStyle","none","LineWidth",1.5)
    xticks(0:20:2*theta_max)
    xlim([30,2*theta_max])
    xlabel("角度2\theta/(°)")
    ylabel("相对强度I/(f)")
    title("点心晶格框架衍射谱")
    set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
    set(gca,"GridLineStyle",":","GridLineWidth",1.5,"TickDir","out")
    figure()
    hold on
    plot(2*abs(rad2deg(theta_range)), (total+total_m).^2, "LineWidth",2)
    maxIndices = islocalmax((total+total_m).^2,"MinProminence",0.01);
    plot(2*abs(rad2deg(theta_range)).*maxIndices, (total+total_m).^2.*maxIndices,"Color",[1,0.3,0.1],"Marker","|","Markersize",12,"LineStyle","none","LineWidth",1.5)
    xticks(0:20:2*theta_max)
    xlim([30,2*theta_max])
    xlabel("角度2\theta/(°)")
    ylabel("相对强度I/(f)")
    title("总衍射谱")
    set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
    set(gca,"GridLineStyle",":","GridLineWidth",1.5,"TickDir","out")
end
end
function [R_x, R_y] = Rot(M, theta, pos)
% M = [row,row,row]
x = pos(1); y = pos(2);
R = [[cos(theta) -sin(theta) x*(1-cos(theta))+y*sin(theta)]; [sin(theta) cos(theta) y*(1-cos(theta))-x*sin(theta)]; [0 0 1]];
trans = M*R';
R_x = trans(:,1);
R_y = trans(:,2);
end
