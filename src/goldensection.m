%%%%%黄金分割法确定步长%%%%%%
function l=goldensection(S_a,S_b,e_golden,x,y,A,d,t,nu,n)

while(S_b-S_a>e_golden)
    b=S_b-0.618*(S_b-S_a);
    a=S_a+0.618*(S_b-S_a);
    fb=sum(fv(y-A*(x+b*d),t,nu,n));
    fa=sum(fv(y-A*(x+a*d),t,nu,n));
    if fb<=fa
        S_b=a;
    else
        S_a=b;
    end
end
l=(S_a+S_b)/2;
    

