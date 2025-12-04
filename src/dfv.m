function y=dfv(x,tau,nu,n)%t分位数的光滑近似函数的导数 u是光滑近似参数 t是分位数 输出梯度向量
%nu=1;
y=zeros(n,1);
for i=1:n
    y(i)=(x(i)<(tau-1)*nu)*(tau-1)+x(i).*(x(i)>=(tau-1)*nu).*(x(i)<=tau*nu)/nu+tau*(tau*nu<x(i));
end
