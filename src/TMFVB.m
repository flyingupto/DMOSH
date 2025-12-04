function [x_OP,C,L,F] = TMFVB(A, y,e_golden,xtrue,mu,nu,nu_end,t,gamma,regular)

[n, p] = size(A);

x=ones(p,1);
outiter=0;%外迭代次数
iter=0;%总迭代次数

zeta2= gamrnd(1, abs(x)+1);
tr=trace(A*A')/n;%计算利普希茨常数用
x_OP=[];
K=0;
tic;
%Q=eye(p);
%e0=0.05;
while nu>nu_end
    Q=eye(p);
    nu=mu*nu;
    k=0;
    e=10;
    maxk=8;
    while k<maxk %修改之前不乘mu
        switch regular
            case 'DP'%DMOSH
                %invlambda2=sampleInverseGaussian(x.^2, zeta2*sigma2,p);
                ee=1e-4;
                x(x<ee)=ee;
                invlambda2=iG(sqrt(zeta2./(x.^2)), zeta2,p);
                zeta2=exprnd(invlambda2./2);
                Lambda=sqrt(zeta2);
            case 'ADALASSO'
                %invlambda2=sampleInverseGaussian(x.^2, zeta2*sigma2,p);
                ee=1e-4;
                x(x<ee)=ee;
                Lambda=gamma./x;
            case 'LASSO'
                Lambda=gamma.*ones(p,1);
            case 'MCP'
                Lambda=gamma.*fd_mcp(x,0.25,4,p);
            case 'SCAD'
                Lambda=gamma.*fd_scad(x,0.25,4,p);
            case 'DPopt'
                switch n
                    case 220
                        r=2;
                        delta=5;
                    otherwise
                        r=2;
                        delta=9;
                end
                Lambda=gamma.*fd_doublepareto(x,r,delta,p);
            otherwise
                Lambda=gamma;
        end
        gk=-A'*dfv(y-A*x,t,nu,n)/n;
        xk=x;
        d=-Q*gk;

        %%%%%%%%%%进退法选择步长区间%%%%%%%%%%%%
        [S_a,S_b]=linesearch(0,nu/tr,x,y,A,d,t,nu,n);

        %%%%%%%%%黄金分割法确定步长%%%%%%%%%%%%%
        alpha=goldensection(S_a,S_b,e_golden,x,y,A,d,t,nu,n);

        z=xk+alpha*d;%参数更新forward
        x=sign(z).*max(abs(z)-alpha*Lambda,0);%参数更新backward
        gk1=-A'*dfv(y-A*x,t,nu,n)/n;
        xk1=x;
        m=gk1-gk; %由于g1和g2差别不大 所以迭代到后期 m会趋于0 导致QNaN
        if m==0
            cc=gamma
            break;
        end
        %
        %         if k>100
        %             c=1
        %             break;
        %         end
        s=xk1-xk;
        Qm=Q*m;
        Q=Q+s*s'/(s'*m)-Qm*Qm'/(m'*Qm);
        k=k+1;
        iter=iter+1;
        L(:,iter)=Lambda;
        %e=s'*s;
    end
    if m==0
        break;
    end
    outiter=outiter+1;
    K=K+k;
    FF=sum(f(y-A*x,t,n));
    F(outiter)=FF/n+sum(Lambda.*abs(x));%非近似的目标函数

    %%%%%%%%%%%%%%%%%%%%%
    x_OP(:,outiter)=x;
end
t2=toc;

%Sparse=sum(abs(sign(x)))/p;%%%%<1>稀疏度
FP=0;%%%%<2>False Positive假阳性
FN=0;%%%%<3>False Negative假阴性
for q=1:p
    if x(q)==0
        if xtrue(q)~=0
            FN=FN+1;
        end
    else
        if xtrue(q)==0
            FP=FP+1;
        end
    end
end
FN=FN/p;
FP=FP/p;
C=[];

C(1,1)=FN;
C(2,1)=sum(f(y-A*x,t,n))/n;%MQE
C(3,1)=(y-A*x)'*(y-A*x)/n;%MSE
C(4,1)=sum(abs(x(1:20)-xtrue(1:20)))/sum(abs(xtrue(1:20)));%RENC
C(5,1)=sum(abs(x-xtrue))/sum(abs(xtrue));%REC
C(6,1)=t2;%Time
C(7,1)=FP;
C(8,1)=K/outiter;%K


toc
end
