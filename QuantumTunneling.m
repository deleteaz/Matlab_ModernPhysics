%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quantum Tunneling

% Author:   ZhiXiang Wei
% Email:    3011860885@qq.com
% Version:  Matlab2022b
% Config:   AMD Ryzen7 6800H && RTX3060
% Usage:    1.no need to change anything, just run.
%           2.d_max: float, Maximum of Barrier width.
%           3.U_max: float, Maximum of Barrier Height.
%           4.E_max: float, Maximum of Energy you want to research.
%           5.delta: float, Step of variable[d,U,E].
%           6.tol: float, Tolerance of error.
%           7.Base on you researching problem, selecting the type of Range.
%           8.y0: float, Bias of plotting.
%           9.diplay: Plot(1) or not(0).
%           10.N: int, Number of barrier divisions.
%           11.type: int, Type of barrier.
%           Square-barrier(1)/moving-Square-barrier(2)/Trapezoidal-barrier(3)/Trigonometric-barrier(4)/Gauss-barrier(5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc;clear;close all;
% Basis
d_max = 4;
U_max = 4;
E_max = 5.5;
delta = 1e-2;
tol = 1e-9;

% Range
d_limit = d_max;
% d_limit = (0.1:d_max*delta:d_max)';
U_limit = U_max;
% U_limit = (0:U_max*delta:U_max)';
% E_limit = E_max;
E_limit = ([tol,E_max*delta:E_max*delta:E_max])';

% Plotting
y0 = 1;
display = 1;
N = 500;
type = 1;

[R_sim, T_theory] = CompleteProcess(d_limit, U_limit, E_limit, N, tol, d_max, U_max, E_max, type, y0, display);
Transmission(R_sim, T_theory, d_limit, U_limit, E_limit, type);


function [R_sim, T_theory] = CompleteProcess(d_limit, U_limit, E_limit, N, tol, d_max, U_max, E_max, type, y0, display)
arguments
d_limit (:,1) = 4;
U_limit (:,1) = 4;
E_limit (:,1) = 4;
N int32 = 500;
tol double = 1e-9;
d_max double = 4;
U_max double = 4;
E_max double = 5.5;
type int32 = 1;
y0 double = 1;
display int32 = 0;
end
% Constants
m = 1;
hbar = 1;
e = 1;
% Arguments
J = zeros(length(E_limit), length(U_limit), length(d_limit));
J_D = zeros(size(J));
J_R = zeros(size(J));
R_sim = zeros(size(J));
T_theory = zeros(size(J));

color = ["red","green","blue","yellow","cyan"];
figure()
hold on
h1 = plot(0, 0, "Color",color(1),"LineWidth",3);
h1_1 = plot(0, 0, "Color",color(1),"LineWidth",1.5,"LineStyle","--");
h1_2 = plot(0, 0, "Color",color(3),"LineWidth",1.5,"LineStyle","--");
h2 = plot(0, 0, "Color",color(4),"LineWidth",3);
h2_1 = plot(0, 0, "Color",color(4),"LineWidth",1.5,"LineStyle","--");
h2_2 = plot(0, 0, "Color",color(2),"LineWidth",1.5,"LineStyle","--");
h3 = plot(0, 0, "Color",color(5),"LineWidth",3);
hl1 = xline([0, 0], "LineWidth",2,"LineStyle","--","Color",[0.3,0.3,0.3]);
hl2 = yline(0, "LineWidth",2,"Color",[0.3,0.3,0.3]);
hl3 = plot(0, 0, "LineWidth",2,"Color",[0.3,0.3,0.3]);
txt1 = text(0, 0, sprintf("U = %.2f to %.2f \nE = %.2f \nU > E", 0, 0, 0), ...
    "FontSize",12,"FontWeight","bold","Interpreter","tex");
txt2 = text(0, 0, sprintf("d = %.2f",0), ...
    "FontSize",12,"FontWeight","bold","Interpreter","tex");
ylim([-0.04*U_max,(1+0.03)*U_max])
xlabel("x")
ylabel("U(x)")
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5,"TickDir","out")
for i1 = 1:length(d_limit)
    d0 = d_limit(length(d_limit)+1-i1);
    spacePoint = [-10, 0, d_max, 15];
    % Site
    x1 = linspace(spacePoint(1), spacePoint(2), 100)';
    x2 = linspace(spacePoint(2), spacePoint(3), N)';
    x3 = linspace(spacePoint(3), spacePoint(4), 100)';
    for i2 = 1:length(U_limit)
        U0 = U_limit(i2);
        switch type
            case 1
                U = U0*ones(N,1); % 方势垒
            case 2
                U = U0*zeros(N,1); U(270-2*i1:400-2*i1) = U0; % 移动方势垒
            case 3
                U = U0-x2/d_max*i1*0.04; % 梯形势垒
            case 4
                U = U0*sin(x2); % 三角函数势垒
            case 5
                U = 3*exp( -((x2-d_max/2)/0.5/d0).^2 ); % 高斯势垒
        end
        U0max = max(U); U0min = min(U);
        for i3 = 1:length(E_limit)
            E0 = E_limit(i3);
            U(abs(E0 - U) < tol) = U(abs(E0 - U) < tol) + tol;

            M_shift = zeros(2,2,N);
            k2 = zeros(N,1);
            % Transfer
            j = 1;
            k1 = sqrt( 2*m*(E0) )/hbar;
            k2(j) = sqrt( 2*m*(E0-U(j)) )/hbar;
            eta = k1/k2(j);
            M_shift(:,:,j) = 0.5 * ...
                [[(1+eta)*exp(1i*(k1-k2(j))*x2(j)), (1-eta)*exp(-1i*(k1+k2(j))*x2(j))]; ...
                [(1-eta)*exp(1i*(k1+k2(j))*x2(j)), (1+eta)*exp(-1i*(k1-k2(j))*x2(j))]];
            for j = 2:N-1
                k1 = sqrt( 2*m*(E0-U(j-1)) )/hbar;
                k2(j) = sqrt( 2*m*(E0-U(j)) )/hbar;
                eta = k1/k2(j);
                M_shift(:,:,j) = 0.5 * ...
                    [[(1+eta)*exp(1i*(k1-k2(j))*x2(j)), (1-eta)*exp(-1i*(k1+k2(j))*x2(j))]; ...
                    [(1-eta)*exp(1i*(k1+k2(j))*x2(j)), (1+eta)*exp(-1i*(k1-k2(j))*x2(j))]] * ...
                    M_shift(:,:,j-1);
            end
            j = N;
            k1 = sqrt( 2*m*(E0-U(end)) )/hbar;
            k2(j) = sqrt( 2*m*(E0) )/hbar;
            eta = k1/k2(j);
            M_shift(:,:,j) = 0.5 * ...
                [[(1+eta)*exp(1i*(k1-k2(j))*x2(j)), (1-eta)*exp(-1i*(k1+k2(j))*x2(j))]; ...
                [(1-eta)*exp(1i*(k1+k2(j))*x2(j)), (1+eta)*exp(-1i*(k1-k2(j))*x2(j))]] * ...
                M_shift(:,:,j-1);

            % Simulated reflectivity
            R_sim(i3,i2,i1) = abs(M_shift(2,1,end)/M_shift(2,2,end)).^2;
            % Theoretical reflectivity（仅方势垒！）
            if type < 2
                T_theory(i3,i2,i1) = (4*k2(end)^2*k2(end-1)^2) / ((k2(end)^2-k2(end-1)^2)^2*sin(k2(end-1)*d0)^2 + 4*k2(end)^2*k2(end-1)^2);
            end

            if display < 1
            else
                % x1:Incident wave and Reflected wave
                A1 = 0.2;
                A2 = -M_shift(2,1,end)/M_shift(2,2,end) * A1; A2 = squeeze(A2);
                psi1_1 = A1 .* exp(1i*k2(end).*x1); % 入射波 (Right)
                psi1_2 = A2 .* exp(-1i*k2(end).*x1); % 反射波 (Left)
                % x2:Attenuation wave and Growth wave
                B1 = M_shift(1,1,:)*A1 + M_shift(1,2,:)*A2; B1 = squeeze(B1);
                B2 = M_shift(2,1,:)*A1 + M_shift(2,2,:)*A2; B2 = squeeze(B2);
                psi2_1 = B1 .* exp(1i*k2.*x2);  % 衰减波（Right）
                psi2_2 = B2 .* exp(-1i*k2.*x2); % 增长波（Left）
                % x3:Transmission wave
                C1 = M_shift(1,1,end)*A1 + M_shift(1,2,end)*A2; C1 = squeeze(C1);
                psi3_1 = C1 .* exp(1i*k2(end).*x3); % 透射波（Right）

                set(h1, "XData", x1, "YData",(real(psi1_1 + psi1_2) + y0).^2);
                set(h1_1, "XData", x1, "YData",(real(psi1_1) + y0).^2);
                set(h1_2, "XData", x1, "YData",(real(psi1_2) + y0).^2);
                set(h2, "XData", x2, "YData",(real(psi2_1 + psi2_2) + y0).^2);
                set(h2_1, "XData", x2, "YData",(real(psi2_1) + y0).^2);
                set(h2_2, "XData", x2, "YData",(real(psi2_2) + y0).^2);
                set(h3, "XData", x3, "YData",(real(psi3_1) + y0).^2);
                set(hl1(1), "Value", x2(1)); set(hl1(2), "Value", x2(end));
                set(hl2, "Value", 0);
                set(hl3, "XData", [x2(1); x2; x2(end)], "YData", [0; U; 0]);
                if E0 < U_max
                    set(txt1, "Position", [x1(1)+1, 1.25*(real(A1+A2)+y0).^2], "String", sprintf("U = %.2f to %.2f \nE = %.2f \nU > E", U0min, U0max, E0));
                else
                    set(txt1, "Position", [x1(1)+1, 1.25*(real(A1+A2)+y0).^2], "String", sprintf("U = %.2f to %.2f \nE = %.2f \nU < E", U0min, U0max, E0));
                end
                set(txt2, "Position", [x2(1), 0.9*gca().YLim(2)], "String", sprintf("d = %.2f", d0));
                drawnow
            end
        end
    end
end
drawnow
end

function Transmission(R_sim, T_theory, d_limit, U_limit, E_limit, type)
arguments
    R_sim (:,:,:)
    T_theory (:,:,:)
    d_limit (:,:,:)
    U_limit (:,:,:)
    E_limit (:,:,:)
    type int32 = 2;
end
R_sim = squeeze(R_sim);
if length(E_limit) > 1
    figure()
    color = ["red","green","blue","yellow","cyan"];
    hold on
    plot(E_limit./max(U_limit), 1-R_sim, "Color",color(1),"LineWidth",3)
    if type < 2
        plot(E_limit./max(U_limit), T_theory, "Color",color(3),"LineWidth",2,"LineStyle","--")
    end
    % If you want to use the following code, Please add your own code after [A1,A2,B1,B2,C1], which is easy.
    % plot(E_limit./max(U_limit), -J_R./J, "Color",color(5),"LineWidth",3)
    % plot(E_limit./max(U_limit), J, "Color",color(4),"LineWidth",1.5)
    % plot(E_limit./max(U_limit), J_D, "Color",color(2),"LineWidth",1.5)
    % plot(E_limit./max(U_limit), -J_R, "Color",color(3),"LineWidth",1.5)
    xlabel("能量比E/U")
    ylabel("透射系数T")
    legend(["仿真值","理论值"],"Location","northwest")
    set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
    set(gca,"GridLineStyle",":","GridLineWidth",1.5,"TickDir","out")
end
if length(d_limit) > 1
    figure()
    hold on
    plot(d_limit, 1-R_sim, "Color",color(1),"LineWidth",3)
    % plot(d, -J_R./J, "Color",color(1),"LineWidth",3)
    xlabel("高斯势垒方差衰减d")
    ylabel("透射系数T")
    set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on","Box","on")
    set(gca,"GridLineStyle",":","GridLineWidth",1.5,"TickDir","out")
end
end

