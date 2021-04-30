x = randn(100,1);
noise = randn(size(x))*.3;
 
p = [randi(5) rand*5];
y = x.^p(1) + p(2) + noise;
 
 
objFunction = @integerfun;
 
params0 = [2 2];
P  = fmincon(@(params) objFunction(params, y, x), params0);  
P2 = lsqnonlin(@(params) objFunction(params, y, x), params0);
 
 
lb = [0 0];
ub = [100 100];
%opts = optimoptions('surrogateopt','PlotFcn',"surrogateoptplot","ConstraintTolerance",1e-6);
opts = optimoptions('surrogateopt','PlotFcn',[], "ConstraintTolerance",1e-6);
 
intcon =   1;
P3 = surrogateopt(@(params) objFunction(params, y, x),lb,ub,intcon, opts);
 
disp([p; P; P2; P3])
figure(1), bar([p; P; P2; P3]');

function [err, pred] = integerfun(params, y, x)
 
n = round(params(1));
k = params(2);
pred = x.^n + k;
err = sum((y-pred).^2);
 
end