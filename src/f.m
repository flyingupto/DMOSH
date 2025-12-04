%t分位数函数
function y=f(x,t,n)
y=zeros(n,1);
for i=1:n
    if x(i)<0
        y(i)=x(i)*(t-1);
    else
        y(i)=x(i)*t;
    end
end
