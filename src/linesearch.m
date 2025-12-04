%进退法确定线搜索区间 其中目标函数为近似的check函数
function [a,b]=linesearch(alpha0,h0,x,y,A,d,t,nu,n)%B是当前点 即[B_in,B] alpha0是当前步长 h0是步骤搜索长度 d是方向
%step1
f0=sum(fv(y-A*(x+alpha0*d),t,nu,n));
k=0;
%nu=1;
while(1)
    %step2
    alpha1=alpha0+h0;
    f1=sum(fv(y-A*(x+alpha1*d),t,nu,n));
    
    if f1<f0
        %step3
        h1=2*h0;
        Alpha=alpha0;
        alpha0=alpha1;
        f0=f1;
        k=k+1;
        %k+1后变量的更新
        alpha0=alpha1;
        h0=h1;
    else
        %step4
        if k==0
            h1=-1*h0;
            Alpha=alpha0;
            alpha0=alpha1;
            f1=f0;
            k=1;
            %k+1后变量的更新
            alpha0=alpha1;
            h0=h1;
        else
            a=min(Alpha,alpha1);
            b=max(Alpha,alpha1);
            break;
        end
    end
end