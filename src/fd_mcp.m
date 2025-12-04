function y=fd_mcp(x,lambda,gamma,p)%MCPµÄµ¼Êý
t=abs(x)./gamma;
y=zeros(p,1);
for i=1:p
    y(i)=(t(i)<=lambda).*(lambda-sign(x(i)).*t(i));
end