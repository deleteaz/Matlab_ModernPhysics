%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lorentz Bimodal Fitting Algorithm

% Author:   ZhiXiang Wei
% Email:    3011860885@qq.com
% Version:  Matlab2022b
% Config:   AMD Ryzen7 6800H && RTX3060
% Usage:    1.Please change dk,I before run.
%           2.dk: array, Raman Shift.
%           3.I: array, light intensity.
%           4.pos: float, median Raman Shift of two peaks.
%           5.dk1, dk2: the final result(Raman Shift).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
dk = T(c1:c2,1);
I = T(c1:c2,2);
pos = 556
x = dk; yTrue = I;
result = solvefit(x,yTrue,pos);
dk1 = result(2)
dk2 = result(5)

function pop_gbest = solvefit(x,yTrue,pos)
GER = 500;
IGER = 1;
POP_SIZE = 200;
%x = [pa1,pb1,pc1,pa2,pb2,pc2,pd1]
LB = [0,pos-2,0,  0,pos-2,0, -5];
UB = [5,pos+2,0.8,5,pos+2,0.8,5];
DIM = length(LB);
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
        [fitness_param, fitness_cross] = fobj(pop_cross(i,:),yTrue,x);
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
figure()
plot(record,"Color",[0.3,0.5,1],"LineWidth",2,"LineStyle","-")
xlabel("Iterations")
ylabel("RMSE")
yline(record(GER-1),"Color",[0.5,0.5,0.5],"LineWidth",1,"LineStyle","--")
text(GER-GER*0.1,record(GER-1),num2str(record(GER-1)),"FontSize",10,"FontWeight","bold")
set(gca,"FontSize",15,"FontWeight","bold")
figure()
hold on
xfit = linspace(min(x),max(x),1000);
y1 = pop_gbest(1)./((xfit-pop_gbest(2)).^2+pop_gbest(3));
y2 = pop_gbest(4)./((xfit-pop_gbest(5)).^2+pop_gbest(6));
yT = yTrue*(ymax - ymin) + ymin;
y0Pred = (pop_gbest(7) + y1 + y2)*(ymax - ymin) + ymin;
y1 = y1*(ymax - ymin) + ymin;
y2 = y2*(ymax - ymin) + ymin;
plot(x,yT,"LineWidth",1.2,"LineStyle","-","Marker",".")
plot(xfit,y0Pred,"LineWidth",1.8,"LineStyle","-")
plot(xfit,y1,"LineWidth",1.5,"LineStyle","-")
plot(xfit,y2,"LineWidth",1.5,"LineStyle","-")
[y1_max,y1_id] = max(y1);
[y2_max,y2_id] = max(y2);
xline([xfit(y1_id),xfit(y2_id)],"Color",[0.5,0.5,0.5],"LineWidth",1,"LineStyle","--");
text(xfit(y1_id),y1_max,num2str(xfit(y1_id)),"FontSize",12,"FontWeight","bold");
text(xfit(y2_id),y2_max,num2str(xfit(y2_id)),"FontSize",12,"FontWeight","bold");
xlim([547,560])
set(gca,"FontSize",15,"FontWeight","bold")
end

function [y_param, y] = fobj(param,yTrue,x)
pa1 = param(1);
pb1 = param(2);
pc1 = param(3);
pa2 = param(4);
pb2 = param(5);
pc2 = param(6);
pd1 = param(7);
yPred = pd1 + pa1./((x-pb1).^2+pc1) + pa2./((x-pb2).^2+pc2);
y = sqrt(sum((yTrue-yPred).^2)/length(yTrue));
y_param = [pa1,pb1,pc1,pa2,pb2,pc2,pd1];
end
