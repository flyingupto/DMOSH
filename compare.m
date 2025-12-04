clc
clear
clf
close all
addpath('../src/');
load('initialn200.mat');
y=y_noisy;
[n, p] = size(A);

e_golden=0.0001;
mu=0.95;%外迭代中 近似参数mu每次的下降比例
nu=1;%\nu的初始值
nu_end=0.01;
t=0.75;%分位数值

regular=["DP","LASSO"];
switch n
    case 200
        switch t
            case 0.25
                gamma=0.375;
            case 0.5
                gamma=0.43;
            case 0.75
                gamma=0.41;
        end
    case 500
        switch t
            case 0.25
                gamma=0.35;
            case 0.5
                gamma=0.39;
            case 0.75
                gamma=0.36;
        end
    case 1000
        switch t
            case 0.25
                gamma=0.375;
            case 0.5
                gamma=0.4;
            case 0.75
                gamma=0.37;
        end
end
%e0=0.05;%内迭代终止误差

repeat=50;
x_compare=cell(1, length(regular));
C_compare=cell(1, length(regular));
L_compare=cell(1, length(regular));
F_compare=cell(1, length(regular));

for i=1:length(regular)
    for j=1:repeat
        [x,C,L,F] = TMFVB(A, y,e_golden,xtrue,mu,nu,nu_end,t,gamma,regular(i));
        x_compare{i}=[x_compare{i},x(:,end)];
        C_compare{i}=[C_compare{i},C];
        L_compare{i}=[L_compare{i},L(:,end)];
        F_compare{i}=[F_compare{i};F];
    end
end
% for i=1:length(regular)
%     j=1;
%     while j<=repeat
%         %for j=1:repeat
%         [x,C,L,F] = TMFVB(A, y,e_golden,xtrue,mu,nu,nu_end,t,gamma,regular(i));
%         if C(1)==0
%             x_compare{i}=[x_compare{i},x(:,end)];
%             C_compare{i}=[C_compare{i},C];
%             L_compare{i}=[L_compare{i},L(:,end)];
%             F_compare{i}=[F_compare{i};F];
%             j=j+1;
%         else
%         end
%     end
% end

x_mean=zeros(p, length(regular));
C_mean=zeros(8, length(regular));
L_mean=zeros(p, length(regular));
x_std=zeros(p, length(regular));
C_std=zeros(8, length(regular));
L_std=zeros(p, length(regular));
regular=["DP","LASSO"];
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
xlabel('iterative number');
ylabel('$\Lambda$','interpreter','latex','Fontsize',18);
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