%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Laser Mode

% Author:   ZhiXiang Wei
% Email:    3011860885@qq.com
% Version:  Matlab2022b
% Config:   AMD Ryzen7 6800H && RTX3060
% Usage:    1.no need to change anything, just run.
%           2.Len_min, Len_max, dLen: float, Range of the researching Cavity length.
%           3.Lam: float, Wavelength of laser.
%           4.r: float, Reflectivity of Lens.
%           5.N: int, The number of laser reflections.
%           6.n: int, Constant.
%           7.q: int, Number of mode.
%           8.Apt: float, r of lens.
%           9.Nx, Ny: int, Division precision.
%           10.display: string, which for plotting.
%           11.m, n: int, Transverse mode order.
%           12.type: int, Square-cavity-mirror(1) or Circular-cavity-mirror(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;
Len_min = 6.000e-2; Len_max = 6.000e-2; dLen = 1e-6;
Lam = 632.8e-9; % [m]

r = 0.9; % [%]
N = 100;
n = 1;
q = 10000;
mode_1(Len_min, Len_max, dLen, r, N, n, q);

Apt = 5e-3; % [m]
Nx = 250;
Nz = 250;
display = "phase"; % "intensity"与"phase"
mode_2(Lam, Len_min, Len_max, dLen, r, N, Apt, Nx, Nz, display);

m = 1;
n = 1;
type = 1;
mode_3(Lam ,Len_min ,Len_max, dLen, m, n, type);

function mode_1(Len_min, Len_max, dLen, r, N, n, q)
arguments
    Len_min double = 6e-3;
    Len_max double = 6e-3;
    dLen double = 1e-6;
    r double = 0.9;
    N int32 = 100;
    n int32 = 1;
    q int32 = 10000;
end
figure()
hold on
% Constant
c = 3e8; % [m]
N = double(N); n = double(n); q = double(q);
Len_range = Len_min:dLen:Len_max; % [m]
Len_l = length(Len_range);

h1 = plot(0, 0, "LineWidth", 2, "Color", [1*1/Len_l,0.2*1/Len_l,1-1*1/Len_l]);
txt1 = text(0.05, 0.9, sprintf("腔长 L = %.2f m", 1), "FontSize", 12, "FontWeight","bold", "Units", "normalized");
txt2 = text(0.05, 0.8, sprintf("FSR = %.3f GHz", 1), "FontSize", 12, "FontWeight","bold", "Units", "normalized");
txt3 = text(0.05, 0.7, sprintf("反射率 r = %.2f", r), "FontSize", 12, "FontWeight","bold", "Units", "normalized");
txt4 = text(0.05, 0.6, sprintf("往返次数 = %d", 0), "FontSize", 12, "FontWeight","bold", "Units", "normalized");
title("激光谐振腔频率响应");
xlabel("频率 (Hz)");
ylabel("相对光强");
set(gca,"FontSize",15,"Box","on","LineWidth",2);
for iLen = 1:Len_l
    A_total = 0;
    Len = Len_range(iLen);
    FSR = c/(2*n*Len);
    center_freq = q * FSR;
    mu = center_freq + (-2.5*FSR:0.005*FSR:2.5*FSR)';
    phi = 4*pi*n*Len*mu/c;
    for k = 1:N
        % A_add = 0.05 + 0.05*exp(-1i*(phi-5.9*k));
        A_add = 0.1;
        A_total = r * A_total .* exp(-1i*phi);
        A_total = A_total + A_add;
        I = abs(A_total).^2;
        set(h1, "XData",mu, "YData",I, "Color",[1*iLen/Len_l,0.2*iLen/Len_l,1-1*iLen/Len_l]);
        set(txt1, "String", sprintf("腔长 L = %.2f m", Len))
        set(txt2, "String", sprintf("FSR = %.3f GHz", FSR/1e9))
        set(txt3, "String", sprintf("反射率 r = %.2f", r))
        set(txt4, "String", sprintf("往返次数 = %d", k))
        pause(0.05)
        drawnow
    end
    % h1 = plot(mu, I, "LineWidth", 2, "Color", [1*i/Len_l,0.2*i/Len_l,1-1*i/Len_l]);
    % drawnow
end
drawnow

figure("Position", [0, 0, 700, 400]);
mu = [center_freq, center_freq + 0.5*FSR];
A_total = zeros(1,length(mu));
I = zeros(N,length(mu));

sg1 = sgtitle(sprintf("往返次数: %d", 0),"FontSize",15);
subplot(1, 2, 1);
hold on
h1 = plot(real(A_total(1)), imag(A_total(1)), "Color",[0.3,0.5,1],"Marker","o","LineWidth",2,"MarkerSize",15);
q1 = quiver(0, 0, real(A_total(1)), imag(A_total(1)), "Color",[0.3,0.5,1],"ShowArrowHead","off","LineWidth",2);
tit1 = title(sprintf("谐振频率: %.3f THz", mu(1)/1e12));
txt1 = text(0.05, 0.85, sprintf("总振幅: %.2f", abs(A_total(1))), "FontSize", 12, "Units", "normalized");
xlim([-1,1]); ylim([-1,1])
grid on
set(gca,"FontSize",13,"Box","on","LineWidth",2,"Xtick",-1:0.5:1);
subplot(1, 2, 2);
hold on
h2 = plot(real(A_total(2)), imag(A_total(2)), "Color",[1,0.5,0.3],"Marker","o","LineWidth",2,"MarkerSize",15);
q2 = quiver(0, 0, real(A_total(2)), imag(A_total(2)), "Color",[1,0.5,0.3],"ShowArrowHead","off","LineWidth",2);
tit2 = title(sprintf("失谐频率: %.3f THz", mu(2)/1e12));
txt2 = text(0.05, 0.85, sprintf("总振幅: %.2f", abs(A_total(2))), "FontSize", 12, "Units", "normalized");
xlim([-0.1,0.1]); ylim([-0.1,0.1])
grid on
set(gca,"FontSize",13,"Box","on","LineWidth",2,"Xtick",-1:0.5:1);
for k = 1:N
    phi = 4*pi*n*Len*mu/c;
    A_new = 0.1;
    A_total = r * A_total .* exp(-1i*phi) + A_new;
    I(k,:) = abs(A_total).^2;
    set(sg1, "String", sprintf("往返次数: %d", k))
    set(h1, "XData",real(A_total(1)), "YData",imag(A_total(1)));
    set(q1, "UData",real(A_total(1)), "VData",imag(A_total(1)));
    set(txt1, "String", sprintf("总振幅: %.2f", abs(A_total(1))));
    set(h2, "XData",real(A_total(2)), "YData",imag(A_total(2)));
    set(q2, "UData",real(A_total(2)), "VData",imag(A_total(2)));
    set(txt2, "String", sprintf("总振幅: %.2f", abs(A_total(2))));
    drawnow;
end
drawnow

figure("Position", [0, 0, 700, 400]);
subplot(1,2,1)
plot(1:N, I(:,1), "Color", [0.3,0.5,1], "LineWidth", 2);
title("谐振频率");
xlabel("往返次数");
ylabel("相对光强");
grid on
set(gca,"FontSize",15,"Box","on","LineWidth",2);
subplot(1,2,2)
plot(1:N, I(:,2), "Color", [1,0.5,0.3], "LineWidth", 2);
title("失谐频率");
xlabel("往返次数");
ylabel("相对光强");
grid on
set(gca,"FontSize",15,"Box","on","LineWidth",2);
end

function mode_2(Lam, Len_min, Len_max, dLen, r, N, Apt, Nx, Nz, display)
arguments
    Lam double = 632.8e-9;
    Len_min double = 6e-3;
    Len_max double = 6e-3;
    dLen double = 1e-6;
    r double = 0.9;
    N int32 = 100;
    Apt double = 5e-3;
    Nx int32 = 250;
    Nz int32 = 250;
    display string = "phase";
end
Len_range = Len_min:dLen:Len_max; % [m]
N = double(N); Nx = double(Nx); Nz = double(Nz);
for iLen = 1:length(Len_range)
    Len = Len_range(iLen);
    x = linspace(-Apt*3, Apt*3, Nx);
    z = linspace(0, Len, Nz);
    [X, Z] = meshgrid(x, z);

    Apt_width = Apt/2;
    k = 2*pi/Lam;
    zR = pi*Apt_width^2/Lam;

    A = exp(-X.^2/Apt_width^2);
    R2 = X.^2;
    R2z = R2 + zR^2;
    phi = k*Z - atan(Z/zR) + k*R2./(2*R2z);

    if display == "intensity"
        field = A .* exp(1i*phi);
    else
        field = phi;
    end

    figure();
    hold on;
    h_img = imagesc(x*1e3, z, abs(field).^2);
    plot(x*1e3, 0, 'k-', 'LineWidth', 2);
    plot(x*1e3, 0 + Len, 'k-', 'LineWidth',2);
    rectangle('Position', [-Apt*1e3, -0.05*Len, 2*Apt*1e3, 0.1*Len],...
        'FaceColor', [0.7,0.7,0.7], 'EdgeColor', 'none');
    rectangle('Position', [-Apt*1e3, Len-0.05*Len, 2*Apt*1e3, 0.1*Len],...
        'FaceColor', [0.7,0.7,0.7], 'EdgeColor', 'none');
    xlabel('横向位置(mm)');
    ylabel('纵向位置(m)');
    title('谐振腔内光场传播');
    set(gca, "XTick",[],"YTick",[])
    daspect([500,1,1])

    dx = x(2)-x(1);
    dz = z(2)-z(1);
    fx = (-Nx/2:Nx/2-1)/(Nx*dx);
    H = exp(-1i*pi*Lam*fx.^2*dz);
    H = fftshift(H);
    mask_R = double(abs(x) < Apt);

    for t = 1:N
        for i = 2:Nz
            ifield = field(i-1, :);
            F = fft(ifield);
            F = F .* H;
            field(i, :) = ifft(F);
            field(i, :) = field(i, :) .* exp(1i*k*dz);
        end
        field(end, :) = r * field(end, :) .* mask_R;
        for i = (Nz-1):-1:1
            ifield = field(i+1, :);
            F = fft(ifield);
            F = F .* H;
            field(i, :) = ifft(F);
            field(i, :) = field(i, :) .* exp(1i*k*dz);
        end

        field(1, :) = r * field(1, :) .* mask_R;
        % field = field / max(abs(field(:)));
        if display == "intensity"
            colormap("colorcube");
            intensity = abs(field).^2;
            set(h_img, 'CData', intensity);
        else
            colormap(jet(256))
            phase_data = angle(field);
            set(h_img, 'CData', phase_data);
        end
        drawnow
    end
    drawnow
end
end

function mode_3(Lam ,Len_min ,Len_max, dLen, m, n, type)
arguments
    Lam double = 632.8e-9;
    Len_min double = 6e-3;
    Len_max double = 6e-3;
    dLen double = 1e-6;
    m int32 = 1;
    n int32 = 1;
    type int32 = 1;
end
m = double(m); n = double(n);
Len_range = Len_min:dLen:Len_max; % [m]
for iLen = 1:length(Len_range)
    Len = Len_range(iLen); % [m]
    if type < 2
        alpha = sqrt((2*pi)/(Len*Lam));
        omega = Len*Lam/pi;

        x = linspace(-3*sqrt(omega),3*sqrt(omega),250);
        y = x;

        Hx = HEM(alpha*x,m);
        Hy = HEM(alpha*y,n);
        Fx = exp(-x.^2/omega);
        Fy = exp(-y.^2/omega);
        if (n-m) < 1
            Con = 0.6^(m+n);
        else
            Con = (0.6-0.2*(m+n))^(m+n);
        end
        fig1 = figure();
        sgtitle(sprintf("方形腔镜TEM_{%d%d}模光强分布", m,n))
        subplot(2,2,1)
        title("H_{m}(x)")
        yyaxis("left")
        plot(x,Hx, "LineWidth",3,"LineStyle","-")
        set(gca,"FontSize",10,"Box","on","LineWidth",2)
        yyaxis("right")
        plot(y,Hy, "LineWidth",2,"LineStyle","--")
        set(gca,"FontSize",10,"Box","on","LineWidth",2);
        subplot(2,2,2)
        title("F(x)")
        yyaxis("left")
        plot(x,Fx, "LineWidth",3,"LineStyle","-")
        set(gca,"FontSize",10,"Box","on","LineWidth",2);
        yyaxis("right")
        plot(y,Fy, "LineWidth",2,"LineStyle","--")
        set(gca,"FontSize",10,"Box","on","LineWidth",2);
        subplot(2,2,3)
        title("F(x)H_{m}(x)")
        yyaxis("left")
        plot(x,Hx.*Fx, "LineWidth",3,"LineStyle","-")
        set(gca,"FontSize",10,"Box","on","LineWidth",2);
        yyaxis("right")
        plot(y,Hy.*Fy, "LineWidth",2,"LineStyle","--")
        set(gca,"FontSize",10,"Box","on","LineWidth",2);
        subplot(2,2,4)
        title("F^{2}(x)H_{m}^{2}(x)")
        yyaxis("left")
        plot(x,(Hx.*Fx).^2, "LineWidth",3,"LineStyle","-")
        set(gca,"FontSize",10,"Box","on","LineWidth",2);
        yyaxis("right")
        plot(y,(Hy.*Fy).^2, "LineWidth",2,"LineStyle","--")
        set(gca,"FontSize",10,"Box","on","LineWidth",2);

        [X,Y] = meshgrid(x,y);

        fig2 = figure();
        sgtitle(sprintf("对称共焦腔方形腔镜TEM_{%d%d}模光强分布",m,n))
        sp1 = subplot(2,2,1);
        HX = HEM(alpha*X,m);
        HY = HEM(alpha*Y,n);
        FX= exp(-X.^2/omega);
        FY = exp(-Y.^2/omega);
        U = Con*HX.*FX.*HY.*FY;
        imagesc((U).^2)
        title("横模光场分布", "FontSize",15,"FontWeight","bold")
        set(gca,"FontSize",12,"FontWeight","bold","Box","on","LineWidth",2)
        set(gca,"XTick",[],"YTick",[])
        axis image
        sp2 = subplot(2,2,3);
        mesh(U)
        axis off
        CM1 = Colormap_3(1,12);
        set(sp1, "Colormap",CM1)
        set(sp2, "Colormap",CM1)

        sp3 = subplot(2,2,2);
        U_total = zeros(size(X));
        FX = exp(-X.^2/omega);
        FY = exp(-Y.^2/omega);
        for i = 0:m
            for j = 0:n
                HX = HEM(alpha*X,i);
                HY = HEM(alpha*Y,j);
                if abs(i-j) < 1
                    Con = 0.6^(i+j);
                else
                    Con = (0.6-0.2*(i+j))^(i+j);
                end
                U_total = U_total + (Con*HX.*FX.*HY.*FY).^2;
            end
        end
        imagesc(U_total)
        title("总光场分布", "FontSize",15,"FontWeight","bold")
        set(gca,"FontSize",12,"FontWeight","bold","Box","on","LineWidth",2)
        set(gca,"XTick",[],"YTick",[])
        axis image
        sp4 = subplot(2,2,4);
        mesh(U_total)
        axis off
        CM2 = Colormap_1(4,12);
        set(sp3, "Colormap",CM2)
        set(sp4, "Colormap",CM2)
    else
        omega = sqrt(Len*Lam/pi);

        x = linspace(-3*omega, 3*omega, 250);
        y = x;

        phi = linspace(0,2*pi);
        zeta = x/omega;
        agr = (sqrt(2)*zeta).^m .* cos(m*phi(1));
        agp = (sqrt(2)*zeta(1)).^m .* cos(m*phi);
        Fr = exp(-zeta.^2);
        Fp = exp(-zeta(1).^2);
        ALr = ALG(2*zeta.^2,m,n);
        ALp = ALG(2*zeta(1).^2,m,n);

        fig1 = figure();
        sgtitle(sprintf("圆形腔镜TEM_{%d%d}模光强分布", m,n))
        subplot(2,4,1)
        plot(x,Fr, "LineWidth",3,"LineStyle","-","Color",[0.3,0.5,1])
        title("F(r) 与 AL(r)")
        set(gca,"FontSize",8,"Box","on","LineWidth",1.5);
        subplot(2,4,2)
        plot(x,agr, "LineWidth",3,"LineStyle","-","Color",[0.3,0.5,1])
        title("ag_{m}(r,phi)")
        set(gca,"FontSize",8,"Box","on","LineWidth",1.5)
        subplot(2,4,3)
        plot(x,agr.*ALr.*Fr, "LineWidth",3,"LineStyle","-","Color",[0.3,0.5,1])
        title("F(r)ag_{m}(r,phi)AL_{m}(r)")
        set(gca,"FontSize",8,"Box","on","LineWidth",1.5);
        subplot(2,4,4)
        plot(x,(agr.*ALr.*Fr).^2, "LineWidth",3,"LineStyle","-","Color",[0.3,0.5,1])
        title("F^{2}(r)ag_{m}^{2}(r,phi)AL_{m}^{2}(r)")
        set(gca,"FontSize",8,"Box","on","LineWidth",1.5);
        subplot(2,4,5)
        plot(x,ALr,"LineWidth",3,"LineStyle","-","Color",[1,0.5,0.3])
        set(gca,"FontSize",8,"Box","on","LineWidth",1.5);
        subplot(2,4,6)
        plot(phi,agp, "LineWidth",3,"LineStyle","-","Color",[1,0.5,0.3])
        set(gca,"FontSize",8,"Box","on","LineWidth",1.5);
        subplot(2,4,7)
        plot(phi,agp.*ALp.*Fp, "LineWidth",3,"LineStyle","-","Color",[1,0.5,0.3])
        set(gca,"FontSize",8,"Box","on","LineWidth",1.5);
        subplot(2,4,8)
        plot(phi,(agp.*ALp.*Fp).^2, "LineWidth",3,"LineStyle","-","Color",[1,0.5,0.3])
        set(gca,"FontSize",8,"Box","on","LineWidth",1.5);

        [X, Y] = meshgrid(x, y);
        R = sqrt(X.^2 + Y.^2); R(R < 1e-16) = 1e-16;
        phi = atan2(Y, X);
        zeta = R/omega;

        fig2 = figure();
        sgtitle(sprintf("对称共焦腔圆形腔镜TEM_{%d%d}模光强分布",m,n))
        sp1 = subplot(2,2,1);
        ag = (sqrt(2)*zeta).^m .* cos(m*phi);
        F = exp(-zeta.^2);
        AL = ALG(2*zeta.^2,m,n);
        if (n-m) < 1
            Con = 0.7^(m+n);
        else
            Con = (0.7-0.2*(m+n))^(m+n);
        end
        U = Con*ag.*F.*AL;
        imagesc((U).^2)
        title("横模光场分布", "FontSize",15,"FontWeight","bold")
        set(gca,"FontSize",12,"FontWeight","bold","Box","on","LineWidth",2)
        set(gca,"XTick",[],"YTick",[])
        axis image
        sp2 = subplot(2,2,3);
        mesh(U)
        axis off
        CM1 = Colormap_3(1,12);
        set(sp1, "Colormap",CM1)
        set(sp2, "Colormap",CM1)

        sp3 = subplot(2,2,2);
        U_total = zeros(size(R));
        F = exp(-zeta.^2);
        for i = 0:m
            for j = 0:n
                ag = (sqrt(2)*zeta).^i .* cos(i*phi);
                AL = ALG(2*zeta.^2,i,j);
                if abs(i-j) < 1
                    Con = 0.7^(i+j);
                else
                    Con = (0.7-0.2*(i+j))^(i+j);
                end
                U_total = U_total + (Con*ag.*F.*AL).^2;
            end
        end
        imagesc(U_total)
        title("总光场分布", "FontSize",15,"FontWeight","bold")
        set(gca,"FontSize",12,"FontWeight","bold","Box","on","LineWidth",2)
        set(gca,"XTick",[],"YTick",[])
        axis image
        sp4 = subplot(2,2,4);
        mesh(U_total)
        axis off
        CM2 = Colormap_1(4,12);
        set(sp3, "Colormap",CM2)
        set(sp4, "Colormap",CM2)
    end
end
end

function H = HEM(z,n)
% Hermite
H = 0;
for k = 0:floor(n*0.5)
    H = H + ((-1).^k.*factorial(n).*(2.*z).^(n-2*k)) ./ (factorial(k).*factorial(n-2*k));
end
end

function AL = ALG(r,m,n)
% Associated Laguerre L_n^m
AL = 0;
for k = 0:n
    AL = AL + (factorial(n+m)*(-r).^k) ./ (factorial(m+k)*factorial(k)*factorial(n-k));
end
end
