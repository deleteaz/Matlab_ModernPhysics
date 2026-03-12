%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lattice Simulator
%
% Author:   ZhiXiang Wei
% Email:    3011860885@qq.com
% Version:  Matlab2022b
% Config:   AMD Ryzen7 6800H && RTX3060
% Usage:    1.no need to change anything, just run.
%           2.material: the material using in database.
%           3.Conventional_size: atomic size.
%           4.Conventional_angle: atomic angle.
%           5.atom_number: atomic number.
%           6.P, I, C, F: atomic type.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;
%% Paramters
% 晶格结构参数（a、b、c || α、β、γ）

% 参考数据库参数
% material = "Au";
% material_info = Materials_Database(material);
% Conventional_size = material_info{1};
% Conventional_angle = material_info{2};
% atom_number = material_info{3};
% material_type = material_info{4};
% P = material_type(1);
% I = material_type(2); 
% C = material_type(3); 
% F = material_type(4);

% 可调参数
Conventional_size = 3.5524 * ones(3,1);
Conventional_angle = 60 * ones(3,1);
atom_number = 28 * ones(4,1);
[P, I, C, F] = deal(1, 0, 0, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Conventional_angle = deg2rad(Conventional_angle);
Crystal = Lattice_Calculation(Conventional_size,Conventional_angle,[P, I, C, F]);
atom_position = Crystal{3};
Primitive_size = Crystal{4};
Primitive_angle = Crystal{5};
atom_cell = Crystal{6};

[a_Con, b_Con, c_Con] = deal(Conventional_size(1), Conventional_size(2), Conventional_size(3));
[alpha_Con, beta_Con, gamma_Con] = deal(Conventional_angle(1), Conventional_angle(2), Conventional_angle(3));
[a_Pri, b_Pri, c_Pri] = deal(Primitive_size(1), Primitive_size(2), Primitive_size(3));
[alpha_Pri, beta_Pri, gamma_Pri] = deal(Primitive_angle(1), Primitive_angle(2), Primitive_angle(3));
supercell_size = [2; 2; 2];

% 实空间
% 基矢量
a1 = a_Con*[1; 0; 0];
a2 = b_Con*[cos(alpha_Con); sin(beta_Con); 0];
a3 = c_Con*[cos(beta_Con);...
    (cos(alpha_Con)-cos(beta_Con)*cos(gamma_Con))/sin(gamma_Con);...
    sqrt(1 - cos(beta_Con)^2 - ((cos(alpha_Con)-cos(beta_Con)*cos(gamma_Con))/sin(gamma_Con))^2)];
[a1, a2, a3] = deal(round(a1,9), round(a2,9), round(a3,9));
Omega_R = dot(a1, cross(a2,a3));
% a1_Pri = a_Pri*[1; 0; 0];
% a2_Pri = b_Pri*[cos(alpha_Pri); sin(beta_Pri); 0];
% a3_Pri = c_Pri*[cos(beta_Pri);...
%     (cos(alpha_Pri)-cos(beta_Pri)*cos(gamma_Pri))/sin(gamma_Pri);...
%     sqrt(1 - cos(beta_Pri)^2 - ((cos(alpha_Pri)-cos(beta_Pri)*cos(gamma_Pri))/sin(gamma_Pri))^2)];
% Omega_cell = dot(a1_Pri, cross(a2_Pri, a3_Pri))
% 倒空间
b1 = 2*pi * cross(a2,a3) / Omega_R;
b2 = 2*pi * cross(a3,a1) / Omega_R;
b3 = 2*pi * cross(a1,a2) / Omega_R;
Omega_I = dot(b1, cross(b2,b3));
% Omega_I = (2*pi)^3 / Omega_R
% dot([a1,a2,a3],[b1,b2,b3])

% primitiveCell
figure()
hold on
view(3)
P_iarray = [0 0 0 0 1 1 1 1];
P_jarray = [0 0 1 1 1 1 0 0];
P_karray = [0 1 0 1 0 1 0 1];

if I > 0
    I_iarray = 0.5;
    I_jarray = 0.5;
    I_karray = 0.5;
else
    I_iarray = 0;
    I_jarray = 0;
    I_karray = 0;
end

if C > 0
    FC_iarray = [0.5 0 0 0.5 0 0];
    FC_jarray = [0.5 0 0 0.5 0 0];
    FC_karray = [1   0 0 0   0 0];

else
    if F > 0
        FC_iarray = [0.5 0   0.5 0.5 1   0.5];
        FC_jarray = [0.5 0.5 1   0.5 0.5 0  ];
        FC_karray = [1   0.5 0.5 0   0.5 0.5];
    else
        FC_iarray = [0 0 0 0 0 0];
        FC_jarray = [0 0 0 0 0 0];
        FC_karray = [0 0 0 0 0 0];
    end
end
iarray = [P_iarray, I_iarray, FC_iarray];
jarray = [P_jarray, I_jarray, FC_jarray];
karray = [P_karray, I_karray, FC_karray];

atom = zeros(size(Conventional_size,1),length(iarray));

for t = 1:8
    i = iarray(t);
    j = jarray(t);
    k = karray(t);
    if mod(t,2) > 0
        if j < 1
            atom(:,t) = [i*a_Con; j*b_Con; k*c_Con];
        else
            atom(:,t) = [i*a1(1) + j*a2(1); j*a2(2); k*c_Con];
        end
    else
        if j < 1
            atom(:,t) = [i*a1(1) + k*a3(1); j*b_Con + k*a3(2); k*a3(3)];
        else
            atom(:,t) = [i*a1(1) + j*a2(1) + k*a3(1); j*a2(2) + k*a3(2); k*a3(3)];
        end
    end
end
for t = 9:size(atom,2)
    i = iarray(t); j = jarray(t); k = karray(t);
    atom(:,t) = [
        i*a1(1) + j*a2(1) + k*a3(1); j*a2(2) + k*a3(2); k*a3(3)];
end

color = ["red","green","blue","cyan","yellow"];
color_line = [0.5,0.5,0.5]; linestyle_line = "-"; linewidth_line = 2;
color_atom = [0.5,0.5,0.5]; linestyle_atom = "none"; linewidth_atom = 2; Marker_atom = "o"; MarekerSize_atom = 25; MarkerFaceColor_atom = color;

line1 = 1:2:8;
line2 = 2:2:8;
for t = 1:length(line1)
    plot3([atom(1,line1(t)),atom(1,line1(t)+1)], ...
        [atom(2,line1(t)),atom(2,line1(t)+1)], ...
        [atom(3,line1(t)),atom(3,line1(t)+1)], ...
        "LineStyle",linestyle_line,"LineWidth",linewidth_line,"Color",color_line);
    if mod(t,4) > 0
        plot3([atom(1,line1(t)),atom(1,line1(t)+2)], ...
            [atom(2,line1(t)),atom(2,line1(t)+2)], ...
            [atom(3,line1(t)),atom(3,line1(t)+2)], ...
            "LineStyle",linestyle_line,"LineWidth",linewidth_line,"Color",color_line);
        plot3([atom(1,line2(t)),atom(1,line2(t)+2)], ...
            [atom(2,line2(t)),atom(2,line2(t)+2)], ...
            [atom(3,line2(t)),atom(3,line2(t)+2)], ...
            "LineStyle",linestyle_line,"LineWidth",linewidth_line,"Color",color_line);
    else
        plot3([atom(1,line1(t)),atom(1,line1(t-3))], ...
            [atom(2,line1(t)),atom(2,line1(t-3))], ...
            [atom(3,line1(t)),atom(3,line1(t-3))], ...
            "LineStyle",linestyle_line,"LineWidth",linewidth_line,"Color",color_line);
        plot3([atom(1,line2(t)),atom(1,line2(t-3))], ...
            [atom(2,line2(t)),atom(2,line2(t-3))], ...
            [atom(3,line2(t)),atom(3,line2(t-3))], ...
            "LineStyle",linestyle_line,"LineWidth",linewidth_line,"Color",color_line);
    end
end
if I > 0
    h2 = plot3(atom(1,9), ...
        atom(2,9), ...
        atom(3,9), ...
        "LineStyle",linestyle_atom,"Marker",Marker_atom,"MarkerSize",MarekerSize_atom,"MarkerFaceColor",MarkerFaceColor_atom(1));
end
if F > 0 || C > 0
    h3 = plot3(atom(1, 10:15), ...
        atom(2, 10:15), ...
        atom(3, 10:15), ...
        "LineStyle",linestyle_atom,"Marker",Marker_atom,"MarkerSize",MarekerSize_atom,"MarkerFaceColor",MarkerFaceColor_atom(5));
end
h1 = plot3(atom(1,1:8), ...
    atom(2,1:8), ...
    atom(3,1:8), ...
    "LineStyle",linestyle_atom,"Marker",Marker_atom,"MarkerSize",MarekerSize_atom,"MarkerFaceColor",MarkerFaceColor_atom(3));
% quiver3([0,0,0] ,[0,0,0] ,[0,0,0] ,atom_cell(1,:) ,atom_cell(2,:) ,atom_cell(3,:) ,"Color",[0,0,0],"LineWidth",2);
axis equal
box on

% superCell
figure()
hold on
view(3)

cell_atom_size = size(atom,2);
cell_atom_num = prod(supercell_size);
supercell = zeros(3,cell_atom_size*cell_atom_num);

flag = 0;
for i = 0:supercell_size(1) - 1
    for j = 0:supercell_size(2) - 1
        for k = 0:supercell_size(3) - 1
            shift = i*a1 + j*a2 + k*a3;
            supercell(:,flag+1:flag+cell_atom_size) = atom + shift;
            flag = flag + cell_atom_size;
        end
    end
end

color_line = [0.5,0.5,0.5]; linestyle_line = "-"; linewidth_line = 2;
color_atom = [0.5,0.5,0.5]; linestyle_atom = "none"; linewidth_atom = 2; Marker_atom = "o"; MarekerSize_atom = 25; MarkerFaceColor_atom = color;

point = 1:size(supercell,2);
lattice_point = point;
line1 = [];
line2 = [];
for i = 1:cell_atom_size:size(supercell,2)
    line1 = [line1,i:2:i+7];
    line2 = [line2,i+1:2:i+7];
end
for t = 1:length(line1)
    plot3([supercell(1,line1(t)),supercell(1,line1(t)+1)], ...
        [supercell(2,line1(t)),supercell(2,line1(t)+1)], ...
        [supercell(3,line1(t)),supercell(3,line1(t)+1)], ...
        "LineStyle",linestyle_line,"LineWidth",linewidth_line,"Color",color_line);
    if mod(t,4) > 0
        plot3([supercell(1,line1(t)),supercell(1,line1(t)+2)], ...
            [supercell(2,line1(t)),supercell(2,line1(t)+2)], ...
            [supercell(3,line1(t)),supercell(3,line1(t)+2)], ...
            "LineStyle",linestyle_line,"LineWidth",linewidth_line,"Color",color_line);
        plot3([supercell(1,line2(t)),supercell(1,line2(t)+2)], ...
            [supercell(2,line2(t)),supercell(2,line2(t)+2)], ...
            [supercell(3,line2(t)),supercell(3,line2(t)+2)], ...
            "LineStyle",linestyle_line,"LineWidth",linewidth_line,"Color",color_line);
    else
        plot3([supercell(1,line1(t)),supercell(1,line1(t-3))], ...
            [supercell(2,line1(t)),supercell(2,line1(t-3))], ...
            [supercell(3,line1(t)),supercell(3,line1(t-3))], ...
            "LineStyle",linestyle_line,"LineWidth",linewidth_line,"Color",color_line);
        plot3([supercell(1,line2(t)),supercell(1,line2(t-3))], ...
            [supercell(2,line2(t)),supercell(2,line2(t-3))], ...
            [supercell(3,line2(t)),supercell(3,line2(t-3))], ...
            "LineStyle",linestyle_line,"LineWidth",linewidth_line,"Color",color_line);
    end
end

P_point = []
I_point = [];
FC_point = [];

if I > 0
    for i = 1:cell_atom_size:size(supercell,2)
        I_point = [I_point, i+8];
    end
    h2 = plot3(supercell(1, I_point), ...
        supercell(2, I_point), ...
        supercell(3, I_point), ...
        "LineStyle",linestyle_atom,"Marker",Marker_atom,"MarkerSize",MarekerSize_atom,"MarkerFaceColor",MarkerFaceColor_atom(1));
end
if F > 0 || C > 0
    for i = 1:cell_atom_size:size(supercell,2)
        FC_point = [FC_point, i+9:i+14];
    end
    h3 = plot3(supercell(1, FC_point), ...
        supercell(2, FC_point), ...
        supercell(3, FC_point), ...
        "LineStyle",linestyle_atom,"Marker",Marker_atom,"MarkerSize",MarekerSize_atom,"MarkerFaceColor",MarkerFaceColor_atom(5));
end
for i = 1:cell_atom_size:size(supercell,2)
    P_point = [P_point, i:i+7];
end
h1 = plot3(supercell(1,P_point), ...
    supercell(2,P_point), ...
    supercell(3,P_point), ...
    "LineStyle",linestyle_atom,"Marker",Marker_atom,"MarkerSize",MarekerSize_atom,"MarkerFaceColor",MarkerFaceColor_atom(3));

supercell = round(supercell,12);
supercell = unique(supercell',"rows","stable")';
axis equal
box on

% XRD
radiation_source = "Cu";
switch radiation_source %Å
    case "Cu"
        lambda = 1.54056;
    case "Cr"
        lambda = 2.28962;
    case "Mo"
        lambda = 0.70926;
    case "Fe"
        lambda = 1.93597;
    case "Co"
        lambda = 1.78896;
    case "Ag"
        lambda = 0.5594;
    otherwise
        error("数据库没有相关材料数据")
end
hkl_max = 4;
hkl_list = zeros((hkl_max+1)^2*hkl_max,3);
flag = 0;
for h = 1:hkl_max
    for k = 0:hkl_max
        for l = 0:hkl_max
            flag = flag + 1;
            hkl_list(flag,:) = [h,k,l];
        end
    end
end
d = zeros(size(hkl_list,1),1);
theta_bragg = zeros(size(hkl_list,1),1);
intensity = zeros(size(hkl_list,1),1);
for t = 1:size(hkl_list,1)
    [h, k, l] = deal(hkl_list(t,1), hkl_list(t,2), hkl_list(t,3));
    G = h*b1 + k*b2 + l*b3;
    % d(t) = a*b*c * ...
    %     sqrt(1 - cos(alpha)^2 - cos(beta)^2 - cos(gamma)^2 + 2*cos(alpha)*cos(beta)*cos(gamma)) /...
    %     sqrt((b*c*h*sin(alpha))^2 + ...
    %          (a*c*k*sin(beta))^2 + ...
    %          (a*b*l*sin(gamma))^2 +...
    %          2*a^2*b*c*k*l*(cos(beta)*cos(gamma)-cos(alpha)) + ...
    %          2*a*b^2*c*h*l*(cos(alpha)*cos(gamma)-cos(beta)) + ...
    %          2*a*b*c^2*h*k*(cos(alpha)*cos(beta)-cos(gamma)));
    d(t) = 2*pi/norm(G);
    theta_bragg(t) = asind(lambda/(2*d(t)));
    F = 0;
    s = cos(2*theta_bragg(t));
    factor = 41.78217 * s^2 * (2.388*exp(-s^2*42.866) + 4.226*exp(-s^2*9.743) + 2.689*exp(-s^2*2.264));
    % factor = s^2*exp(-0.001*s^2*theta_bragg(t));
    % factor = 1;
    for i = 1:size(atom_position,2)
        [x, y, z] = deal(atom_position(1,i), atom_position(2,i), atom_position(3,i));
        F = F + (atom_number(i) - factor) * exp(1i * 2*pi * (h*x + k*y + l*z));
    end
    intensity(t) = abs(F)^2;
    theta_bragg(t) = round(2*real(theta_bragg(t)),9);
end
figure()
hold on
plot([0,180],[0,0],"LineWidth",4,"LineStyle","-","Color",color(1))
stem(theta_bragg,intensity,"LineWidth",2,"LineStyle","-","Marker","none","Color",color(1))
theta_bragg_unique = unique(theta_bragg);
intensity_max = max(intensity(theta_bragg<180));
for t = 1:length(theta_bragg_unique)
    idx = find(theta_bragg == theta_bragg_unique(t));
    if intensity(idx(1)) > 10
    text(theta_bragg_unique(t)-3, intensity(idx(1))+0.03*intensity_max*size(string(num2str(hkl_list(idx,:))),1), ...
        string(num2str(hkl_list(idx,:)) ));
    end
end
xticks(0:20:180)
xlim([-1,181])
ylim([0,1.2*intensity_max])
xlabel("角度2\theta/(°)")
ylabel("相对强度I/(f)")
title("XRD")
set(gca,"FontWeight","bold","FontSize",12.5,"LineWidth",1.5,"XMinorTick","on","YMinorTick","on")
set(gca,"GridLineStyle",":","GridLineWidth",1.5,"TickDir","out")

function y = guass(mu, sigma, x)
y = 1/(sqrt(2*pi)*sigma) * exp(-(x-mu).^2/(2*sigma^2));
end

%%参考文献
% http://staff.ustc.edu.cn/~weng/Teaching/SolidStatePhysics/SSP_RealCrystal.pdf
% 二种求算晶体平面间距d（hkl）的数学方法：doi:10.3866/pku.DXHX20150486
