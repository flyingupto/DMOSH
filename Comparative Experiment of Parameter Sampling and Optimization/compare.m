clc
clear
clf
close all
addpath('../src/');
load('initialn100.mat');
y=y_noisy;
[n, p] = size(A);
format shortG
e_golden=0.0001;
mu=0.95;%外迭代中 近似参数mu每次的下降比例
nu=1;%\nu的初始值
nu_end=0.01;
t=0.5;%分位数值

regular=["DP","MCP","SCAD","DPopt"];
%regular=["DP","DPopt"];
switch n
    case 80
        switch t
            case 0.25
                gamma=[0,0.345,0.33];
            case 0.5
                gamma=[0,0.35,0.3475];
            case 0.75
                gamma=[0,0.3425,0.345];
        end
    case 100
        switch t
            case 0.25
                gamma=[0,0.325,0.3375,0.305];
            case 0.5
                gamma=[0,0.375,0.3725,0.3];
            case 0.75
                gamma=[0,0.37,0.335, 0.30125];
        end
    case 150
        switch t
            case 0.25
                gamma=[0,0.3175,0.3325,0.2025];
            case 0.5
                gamma=[0,0.29,0.305,0.20125];
            case 0.75
                gamma=[0,0.34,0.3225, 0.20375];
        end

    case 220
        switch t
            case 0.25
                gamma=[0,0.343751,0.195,0.1125];
            case 0.5
                gamma=[0,0.345,0.17, 0.105];
            case 0.75
                gamma=[0,0.336,0.2025, 0.11];
        end
end


repeat=200;
x_compare=cell(1, length(regular));
C_compare=cell(1, length(regular));
L_compare=cell(1, length(regular));
F_compare=cell(1, length(regular));

for i=1:length(regular)
    for j=1:repeat
        [x,C,L,F] = TMFVB(A, y,e_golden,xtrue,mu,nu,nu_end,t,gamma(i),regular(i));
        x_compare{i}=[x_compare{i},x(:,end)];
        C_compare{i}=[C_compare{i},C];
        L_compare{i}=[L_compare{i},L(:,end)];
        F_compare{i}=[F_compare{i};F];
    end
end

x_mean=zeros(p, length(regular));
C_mean=zeros(8, length(regular));
L_mean=zeros(p, length(regular));
x_std=zeros(p, length(regular));
C_std=zeros(8, length(regular));
L_std=zeros(p, length(regular));
for i=1:length(regular)
    x_mean(:,i)=mean(x_compare{i},2);
    C_mean(:,i)=mean(C_compare{i},2);
    L_mean(:,i)=mean(L_compare{i},2);
end
for i=1:length(regular)
    x_std(:,i)=std(x_compare{i},0,2);
    C_std(:,i)=std(C_compare{i},0,2);
    L_std(:,i)=std(L_compare{i},0,2);
end
filename=['test',num2str(n),'_',num2str(t*100)];
save(filename);
last_iter=100;
figure
for i=1:length(regular)
    semilogy(L_mean(1:last_iter,i),'LineWidth',1.2)
    %plot(L_mean(1:last_iter,i));
    hold on;
    %ylim([0,111]);
end
legend(regular,'location','BEST')
xlabel('i','Fontsize',13);
ylabel('mean regularization parameter','Fontsize',13);
print([filename,'L.eps'],'-depsc','-r1200');

figure;
plot(xtrue,xtrue,'--','LineWidth',1.2,'HandleVisibility','off');
hold on;
for i=1:length(regular)
    plot(x_mean(:,i),xtrue,'o','LineWidth',1.2);
    hold on;
end
legend(regular,'Fontsize', 11,'location','BEST');
xlabel('$\widehat{x}$','interpreter','latex','Fontsize', 20);
ylabel('$x^{*}$','interpreter','latex','Fontsize', 20);
print([filename,'x.eps'],'-depsc','-r1200');

% for i=1:length(regular)
% figure
% M=C_compare{i}./max(C_compare{i},2);
% boxplot(M')
% end
