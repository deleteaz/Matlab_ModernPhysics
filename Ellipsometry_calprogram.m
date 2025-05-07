%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ellipsometry

% Author:   ZhiXiang Wei
% Email:    3011860885@qq.com
% Version:  Matlab2022b
% Config:   AMD Ryzen7 6800H && RTX3060
% Usage:    1.Please change P0_1,P0_2,A0_1,A0_2 before run.
%           2.P0_1, P0_2, A0_1, A0_2: float, extinction angle / polarization angle.
%           3.phi1: float, incident angle.
%           4.n1, n3: complex, refractive index of medium.
%           5.d, n2: the final result.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;
% Measurement Value
% Extinction angle / polarization angle
P0_1 = 177.33; A0_1 = 50.43;
P0_2 = 86.83; A0_2 = 141.73;
phi1 = 69.87;
n1 = 1;
n3 = 3.882 - 1i*0.019;

[psi_measure, Delta_measure] = Epcal(A0_1,P0_1,A0_2,P0_2);

AreaPlot(n1,n3,phi1);

result = solvefit(psi_measure,Delta_measure,phi1);
d = result(1)
n2 = result(2)

function [pop_gbest]= solvefit(psi_measure,Delta_measure,phi1)
GER = 200;
IGER = 1;
POP_SIZE = 100;
DIM = 2;
%x = [d,  n2]
LB = [100,1];
UB = [120,2];
muCR = 0.5;
muSF = 0.5;
pop = LB+(UB-LB).*rand(POP_SIZE,DIM);
pop_gbest = zeros(1,DIM);
fitness = inf*ones(POP_SIZE,1);
fitness_gbest = inf;
fitness_param_gbest = ones(1,DIM);
record = zeros(GER-1,1);
while GER > IGER
    r1 = randi([1,POP_SIZE],POP_SIZE,1);
    r2 = randi([1,POP_SIZE],POP_SIZE,1);
    r3 = randi([1,POP_SIZE],POP_SIZE,1);
    for i = 1:POP_SIZE
        while r1(i) == i
            r1(i) = randi([1,POP_SIZE],1);
        end
        while any([r2(i) == i,r2(i) == r1(i)])
            r2(i) = randi([1,POP_SIZE],1);
        end
        while any([r3(i) == i,r3(i) == r1(i),r3(i) == r2(i)])
            r3(i) = randi([1,POP_SIZE],1);
        end
    end
    pop_mutat = pop(r1,:) + muSF.*(pop(r2,:)-pop(r3,:));

    j = [1:DIM].*ones(POP_SIZE,DIM);
    jrand = randi([0,DIM],POP_SIZE,DIM);
    rc = rand(POP_SIZE,DIM);
    cross_id = (j == jrand) + (rc < muCR);
    pop_cross = cross_id.*pop_mutat + ~cross_id.*pop;

    LB_id = pop_cross < LB;
    UB_id = pop_cross > UB;
    LB_reset = (LB+pop)/2;
    UB_reset = (UB+pop)/2;
    pop_cross = LB_id.*LB_reset + ~LB_id.*pop_cross;
    pop_cross = UB_id.*UB_reset + ~UB_id.*pop_cross;

    for i = 1:POP_SIZE
        [fitness_param, fitness_cross] = fobj(pop_cross(i,:),[psi_measure,Delta_measure],phi1);
        if  fitness_cross < fitness(i)
            pop(i,:) = pop_cross(i,:);
            fitness(i) = fitness_cross;
            if  fitness_cross < fitness_gbest
                fitness_gbest = fitness_cross;
                fitness_param_gbest = fitness_param;
                pop_gbest = pop_cross(i,:);
            end
        end
    end
    [fitness,best_id] = sort(fitness);
    pop = pop(best_id,:);

    record(IGER) = fitness_gbest;
    IGER = IGER+1;
end

fprintf(['Measurement:\n psi: %.4f \tDelta: %.4f\n\n'...
    'Calculation(%.2fdeg):\n psi: %.4f \tDelta: %.4f\n\n' ...
    'd(膜厚): %.4fnm \tn2(薄膜折射率): %.4f'], ...
    psi_measure, Delta_measure, phi1, fitness_param_gbest(1), fitness_param_gbest(2), ...
    pop_gbest(1), pop_gbest(2));
figure()
plot(record)
xlabel("Iterations")
ylabel("RMSE")
yline(record(GER-1),"Color",[0.5,0.5,0.5],"LineWidth",1,"LineStyle","--")
text(GER-GER*0.1,record(GER-1),num2str(record(GER-1)))
end

function [y_param, y] = fobj(param,yTrue,phi1)
% Calculate Value
lambda = 632.8;
d = param(1);
n1 = 1;
n2 = param(2);
n3 = 3.882 - 1i*0.019;

phi2 = asind(n1/n2*sind(phi1));
phi3 = asind(n1/n3*sind(phi1));
r1p = (n2*cosd(phi1)-n1*cosd(phi2)) / (n2*cosd(phi1)+n1*cosd(phi2));
r1s = (n1*cosd(phi1)-n2*cosd(phi2)) / (n1*cosd(phi1)+n2*cosd(phi2));

r2p = (n3*cosd(phi2)-n2*cosd(phi3)) / (n3*cosd(phi2)+n2*cosd(phi3));
r2s = (n2*cosd(phi2)-n3*cosd(phi3)) / (n2*cosd(phi2)+n3*cosd(phi3));

delta = 2*pi*d*n2*cosd(phi2)/lambda;

G = (r1p + r2p*exp(-1i*2*delta)) / (1 + r1p*r2p*exp(-1i*2*delta)) * (1 + r1s*r2s*exp(-1i*2*delta)) / (r1s + r2s*exp(-1i*2*delta));
x = real(G); y = imag(G);
psi = atand(sqrt(x^2 + y^2));
if x > 0
    Delta = atand(y/x);
elseif x == 0 && y > 0
    Delta = 90;
elseif x == 0 && y < 0
    Delta = -90;
elseif x < 0 && y >= 0
    Delta = atand(y/x) + 180;
elseif x < 0 && y < 0
    Delta = atand(y/x) - 180;
end
y = sqrt((yTrue(1)-psi)^2 + (yTrue(2)-Delta)^2);
y_param = [psi, Delta];
end

function AreaPlot(n1,n3,phi1)
% Why Can't I use Newton's method
step = 200;
d_list = linspace(100,200,step);
n2_list = linspace(1,2,step);
psi_list = zeros(step,step);
Delta_list = zeros(step,step);
G_list = zeros(step,step);
for i = 1:step
    for j = 1:step
        d = d_list(i);
        lambda = 632.8;
        n2 = n2_list(j);

        phi2 = asind(n1/n2*sind(phi1));
        phi3 = asind(n1/n3*sind(phi1));
        r1p = (n2*cosd(phi1)-n1*cosd(phi2)) / (n2*cosd(phi1)+n1*cosd(phi2));
        r1s = (n1*cosd(phi1)-n2*cosd(phi2)) / (n1*cosd(phi1)+n2*cosd(phi2));

        r2p = (n3*cosd(phi2)-n2*cosd(phi3)) / (n3*cosd(phi2)+n2*cosd(phi3));
        r2s = (n2*cosd(phi2)-n3*cosd(phi3)) / (n2*cosd(phi2)+n3*cosd(phi3));

        delta = 2*pi*d*n2*cosd(phi2)/lambda;

        G = (r1p + r2p*exp(-1i*2*delta)) / (1 + r1p*r2p*exp(-1i*2*delta)) ...
            * (1 + r1s*r2s*exp(-1i*2*delta)) / (r1s + r2s*exp(-1i*2*delta));
        G_list(i,j) = G;
        x = real(G); y = imag(G);
        psi_list(i,j) = atand(sqrt(x^2 + y^2));
        if x > 0
            Delta_list(i,j) = atand(y/x);
        elseif x == 0 && y > 0
            Delta_list(i,j) = 90;
        elseif x == 0 && y < 0
            Delta_list(i,j) = -90;
        elseif x < 0 && y >= 0
            Delta_list(i,j) = atand(y/x) + 180;
        elseif x < 0 && y < 0
            Delta_list(i,j) = atand(y/x) - 180;
        end
    end
end

[x_plot, y_plot] = meshgrid(d_list, n2_list);
figure()
subplot(2,2,1)
h = mesh(x_plot, y_plot, psi_list);
xlabel("d（膜厚/nm）")
ylabel("n2（薄膜折射率）")
zlabel("\psi（椭偏角1）")
set(gca,"FontSize",15,"FontWeight","bold")
subplot(2,2,2)
h = mesh(x_plot, y_plot, Delta_list);
xlabel("d（膜厚/nm）")
ylabel("n2（薄膜折射率）")
zlabel("\Delta（椭偏角2）")
set(gca,"FontSize",15,"FontWeight","bold")
subplot(2,2,[3,4])
h = mesh(x_plot, y_plot, abs(G_list));
xlabel("d（膜厚/nm）")
ylabel("n2（薄膜折射率）")
zlabel("G（反射系数比）")
set(gca,"FontSize",15,"FontWeight","bold")
end

function [psi_measure, Delta_measure] = Epcal(A0_1,P0_1,A0_2,P0_2)
% Ellipsoidal parameter
if A0_1 > 90
    psi_measure_1 = 180 - A0_1; Delta_measure_1 = 2*P0_1 - 90;
elseif A0_1 <= 90
    psi_measure_1 = A0_1;
    if 2*P0_1 > 90
        Delta_measure_1 = 2*P0_1 - 270;
    elseif 2*P0_1 <= 90
        Delta_measure_1 = 2*P0_1 + 90;
    end
end
if A0_2 > 90
    psi_measure_2 = 180 - A0_2; Delta_measure_2 = 2*P0_2 - 90;
elseif A0_2 <= 90
    psi_measure_2 = A0_2;
    if 2*P0_2 > 90
        Delta_measure_2 = 2*P0_2 - 270;
    elseif 2*P0_2 <= 90
        Delta_measure_2 = 2*P0_2 + 90;
    end
end
% Eliminate eccentricity error
psi_measure = (psi_measure_1 + psi_measure_2) / 2;
Delta_measure = (Delta_measure_1 + Delta_measure_2) / 2;
fprintf('Measurement(deg):\n psi:\t%.4f \n Delta:\t%.4f',psi_measure,Delta_measure)
end
