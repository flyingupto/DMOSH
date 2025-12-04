function y=fd_scad(x,lambda,gamma,p)%SCADµÄµ¼Êý
t=abs(x);
y=zeros(p,1);
for i=1:p
    y(i)=(t(i)<=lambda)*(lambda)+(t(i)>lambda & t(i)<gamma*lambda).*(gamma*lambda-t(i))/(gamma-1);
end