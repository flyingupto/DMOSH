%t分位数的光滑近似函数 u是光滑近似参数 t是分位数
function y=fv(x,t,nu,n)
%nu=1;
y=zeros(n,1);
for i=1:n
    y(i)=(x(i)<(t-1)*nu).*((t-1)*x(i)-(t-1)^2*nu/2)+(x(i)>t*nu).*(t*x(i)-t^2*nu/2)+(x(i)>=(t-1)*nu & x(i)<=t*nu).*(x(i).*x(i)*0.5/nu);
end
