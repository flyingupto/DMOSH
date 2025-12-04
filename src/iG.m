function Tau2 = iG(mu,lambda,p)
% 初始化输出向量

Tau2 = zeros(p, 1);

% 对每个 mu 和 lambda 生成样本
for i = 1:p
    % 生成标准正态分布的随机数
    v = randn;
    y = v^2;
    % 计算中间值
    x = mu(i) + (mu(i)^2 * y) / (2 * lambda(i)) - ...
        (mu(i) / (2 * lambda(i))) * sqrt(4 * mu(i) * lambda(i) * y + mu(i)^2 * y^2);
    
    % 生成均匀分布的随机数
    u = rand;
    
    % 根据条件选择返回值
    if u <= mu(i) / (mu(i) + x)
        Tau2(i) = x;
    else
        Tau2(i) = mu(i)^2 / x;
    end
end
e=1e-4;
Tau2(Tau2<e)=e;
end