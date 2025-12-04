function y=fd_doublepareto(x,c,a,p)%MCPµÄµ¼Êý
t=abs(x);
y=zeros(p,1);
for i=1:p
    y(i)=c/(a+abs(x(i)));
end