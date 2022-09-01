%% clear environment variables
clear; close all; clc
%% Example Data Process
mu = [1.175,1.19,0.396,1.12,0.346,0.679,0.089,0.73,0.481,1.08]';
Sigma = [
0.40755159 0.03175842 0.05183923 0.05663904 0.0330226 0.00827775 0.02165938 0.01332419 0.0343476 0.02249903
0.03175842 0.9063047 0.03136385 0.02687256 0.01917172 0.00934384 0.02495043 0.00761036 0.02874874 0.01336866
0.05183923 0.03136385 0.19490901 0.04408485 0.03006772 0.01322738 0.03525971 0.0115493 0.0427563 0.02057303
0.05663904 0.02687256 0.04408485 0.19528471 0.02777345 0.00526665 0.01375808 0.00780878 0.02914176 0.01640377
0.0330226 0.01917172 0.03006772 0.02777345 0.34059105 0.00777055 0.02067844 0.00736409 0.02542657 0.01284075
0.00827775 0.00934384 0.01322738 0.00526665 0.00777055 0.15983874 0.02105575 0.00518686 0.01723737 0.00723779
0.02165938 0.02495043 0.03525971 0.01375808 0.02067844 0.02105575 0.68056711 0.01377882 0.04627027 0.01926088
0.01332419 0.00761036 0.0115493 0.00780878 0.00736409 0.00518686 0.01377882 0.95526918 0.0106553 0.00760955
0.0343476 0.02874874 0.0427563 0.02914176 0.02542657 0.01723737 0.04627027 0.0106553 0.31681584 0.01854318
0.02249903 0.01336866 0.02057303 0.01640377 0.01284075 0.00723779 0.01926088 0.00760955 0.01854318 0.11079287];
n = 10;
% create correlation matrix
P = zeros(n,n);
for i=1:n
    for j=1:i
        P(i,j) = Sigma(i,j)/sqrt(Sigma(i,i)*Sigma(j,j)) ;
    end 
end
P = P + triu(P',1); % heatmap(P);

% create new correlation matrix
mu_change = (mu(3)+mu(4))/2;
mu(3) = mu_change - 0.001;
mu(4) = mu_change + 0.001;
mu;
Sigma_change = sqrt(sqrt(Sigma(3,3)* Sigma(4,4)));
Sigma(3,3) = (Sigma_change - 0.001)^2;
Sigma(4,4) = (Sigma_change + 0.001)^2;

P(3,1:2) = P(4,1:2);
P(3,5:10) = P(4,5:10);
P(1:2,3) = P(1:2,4);
P(5:10,3) = P(5:10,4);
P(3,4) = 0.9999;
P(4,3) = 0.9999; 
[Z,L] = eig(P); % heatmap(P);
% create new Sigma
Sigma_new = zeros(n,n);
for i=1:n
    for j=1:i
        Sigma_new(i,j) = P(i,j) * sqrt(Sigma(i,i)*Sigma(j,j)) ;
    end 
end
Sigma_new= Sigma_new + triu(Sigma_new',1);
Sigma = Sigma_new;

lambda_set = [(0:0.1:1) (1:50)];
%% Parameter Determination
% uncertainty matrix for quadratic: Omega
Omega = zeros(n,n);
Omega(logical(eye(n)))= diag(Sigma);
sqrt_Omega = sqrt(Omega);
% uncertainty level for quadratic: k
SR = mu./sqrt(diag(Sigma));
k = mean(SR)/2; % k =  0.688359105674851
% generate quadratic uncertainty set
N = 10000;
rng(0);
d_largest = (-1 + 2 * rand(n, N));    
size_set = [];
scale = 0.2574;
d = scale * d_largest; 
d_quad_set = [];
true_mu_set = [];
for i=1:N 
    if d(:,i)' * inv(Omega) * d(:,i) <= k*k
        d_quad_set = [d_quad_set d(:,i)];
    end
end
size_quad = size(d_quad_set,2) % 1000
true_mu_quad_set = mu + d_quad_set;


% uncertainty level for box: xi
sd = sqrt(diag(Sigma));
xi = sqrt( k*k / ((ones(n,1) ./sd)' * (ones(n,1) ./sd)) ); % 0.112111180274257
% generate box uncertainty set
scale = xi;
d = scale * d_largest; 
d_box_set = [];
true_mu_box_set = [];
for i=1:size_quad
    if abs(d(:,i)) <= xi
        d_box_set = [d_box_set d(:,i)];
    end
end
size_box = size(d_box_set,2); % 1000
true_mu_box_set = mu + d_box_set;
%% MVO: Markowitz Model
MVO.W = [];
MVO.ExpRet = [];
MVO.stdev = [];
for lambda = lambda_set 
    cvx_begin quiet
      variable w(n);  
      maximize(mu'*w - lambda/2*w'*Sigma*w);
      subject to
        ones(1,n)*w == 1;
        w >= 0;
    cvx_end
    MVO.W = [MVO.W w];
    MVO.ExpRet = [MVO.ExpRet mu'*w];
    MVO.stdev = [MVO.stdev sqrt(w' * Sigma * w)];
end

% calculate lambda of CML and set it as Lambda
[sharperatio, num] = max(MVO.ExpRet ./ MVO.stdev); % get optimal SR ==> Lambda
Lambda = lambda_set(num); % 3.799291773095309 18

%% ROquad1: RO with Quadratic Uncertainty Set (original)
ROquad1.W = [];
ROquad1.stdev = [];
ROquad1.ExpRet = [];
ROquad1.ExpRet_original = [];
for lambda = lambda_set 
    cvx_begin quiet
      variable w(n);  
      maximize(mu' * w - k * norm(sqrt_Omega * w) - lambda/2 * w' * Sigma * w);
      subject to
        ones(1,n)*w == 1;
        w >= 0;
    cvx_end
    ROquad1.W = [ROquad1.W w];
    ROquad1.ExpRet = [ROquad1.ExpRet mu' * w - k * norm(sqrt_Omega * w)];
    ROquad1.ExpRet_original = [ROquad1.ExpRet_original mu'*w];
    ROquad1.stdev = [ROquad1.stdev sqrt(w' * Sigma * w)];
end

%%  RObox1: RO with Box Uncertainty Set (original)
RObox1.W = [];
RObox1.stdev = [];
RObox1.ExpRet = [];
RObox1.ExpRet_original = [];
for lambda = lambda_set 
    cvx_begin quiet
      variable w(n);  
      maximize ((mu' * w) - xi * sum(abs(w)) - lambda/2 * w' * Sigma * w);
      subject to
        ones(1,n)*w == 1;
        w >= 0;
    cvx_end
    RObox1.W = [RObox1.W w];
    RObox1.ExpRet = [RObox1.ExpRet (mu' * w) - xi * sum(abs(w))];
    RObox1.ExpRet_original = [RObox1.ExpRet_original mu' * w];
    RObox1.stdev = [RObox1.stdev sqrt(w' * Sigma * w)];
end

%% ROquad3: RO with Quadratic Uncertainty Set (benchmark)
benchmark = ones(1,n)' / n;
% MVO with benchmark 
MVO3.W = [];
MVO3.stdev = [];
MVO3.ExpRet = [];
for lambda = lambda_set 
    cvx_begin quiet
      variable w(n);  
      maximize(mu' * (w-benchmark) - lambda/2 * (w-benchmark)' * Sigma * (w-benchmark));
      subject to
        ones(1,n)*w == 1;
        w >= 0;
    cvx_end
    MVO3.W = [MVO3.W w];
    MVO3.ExpRet = [MVO3.ExpRet mu' * (w-benchmark)];
    MVO3.stdev = [MVO3.stdev sqrt((w-benchmark)' * Sigma * (w-benchmark))];
end
% calculate lambda of CML and set it as Lambda
[sharperatio3, num] = max(MVO3.ExpRet ./ MVO3.stdev); % get optimal SR ==> Lambda
Lambda3 = lambda_set(num); % 2.1765


% ROquad3
ROquad3.W = [];
ROquad3.stdev = [];
ROquad3.ExpRet = [];
ROquad3.ExpRet_original = [];
for lambda = lambda_set 
    cvx_begin quiet
      variable w(n);  
      maximize(mu' * (w-benchmark) - k * norm(sqrt_Omega * (w-benchmark)) - lambda/2 * (w-benchmark)' * Sigma * (w-benchmark));
      subject to
        ones(1,n)*w == 1;
        w >= 0;
    cvx_end
    ROquad3.W = [ROquad3.W w];
    ROquad3.ExpRet = [ROquad3.ExpRet mu' * (w-benchmark) - k * norm(sqrt_Omega * (w-benchmark))];
    ROquad3.ExpRet_original = [ROquad3.ExpRet_original mu' * (w-benchmark)];
    ROquad3.stdev = [ROquad3.stdev sqrt((w-benchmark)' * Sigma * (w-benchmark))];
end

%%  RObox3: RO with Box Uncertainty Set (benchmark)
RObox3.W = [];
RObox3.stdev = [];
RObox3.ExpRet = [];
RObox3.ExpRet_orignal = [];
for lambda = lambda_set 
    cvx_begin quiet
      variable w(n);  
      maximize ((mu' * (w-benchmark)) - xi * sum(abs((w-benchmark))) - lambda/2 * (w-benchmark)' * Sigma * (w-benchmark));
      subject to
        ones(1,n)*w == 1;
        w >= 0;
    cvx_end
    RObox3.W = [RObox3.W w];
    RObox3.ExpRet = [RObox3.ExpRet (mu' * (w-benchmark)) - xi * sum(abs((w-benchmark)))];
    RObox3.ExpRet_orignal = [RObox3.ExpRet mu' * (w-benchmark)];
    RObox3.stdev = [RObox3.stdev sqrt((w-benchmark)' * Sigma * (w-benchmark))];
end

%% ROquad2: RO with Quadratic Uncertainty Set (zero net adjustment)
ROquad2.W = [];
ROquad2.stdev = [];
ROquad2.ExpRet = [];
ROquad2.ExpRet_original = [];
I = eye(n);
D = I;
e = ones(n,1);
for lambda = lambda_set 
    cvx_begin quiet
      variable w(n);  
      sqrt_Omega_modified = Omega - 1/(e'*D*Omega*D'*e)*Omega*D'*e*e'*D*Sigma;
      maximize(mu' * w - k * norm(sqrt_Omega_modified * w) - lambda/2 * w' * Sigma * w);
      subject to
        ones(1,n)*w == 1;
        w >= 0;
    cvx_end
    ROquad2.W = [ROquad2.W w];
    ROquad2.ExpRet = [ROquad2.ExpRet mu' * w - k * norm(sqrt_Omega_modified * w)];
    ROquad2.ExpRet_orignal = [ROquad2.ExpRet_original mu'* w];
    ROquad2.stdev = [ROquad2.stdev sqrt(w' * Sigma * w)];
end

%%  RObox2: RO with Box Uncertainty Set (zero net adjustment)
% find w0 in MVO or ROquad2
lambda = Lambda;
cvx_begin quiet
  variable w(n);  
  maximize(mu' * w - lambda/2 * w' * Sigma * w);
  subject to
    ones(1,n)*w == 1;
    w >= 0;
cvx_end
w0 = w;


RObox2.W = [];
RObox2.stdev = [];
RObox2.ExpRet = [];
RObox2.ExpRet_original = [];
RObox2.exitflag = [];
for lambda = lambda_set 
    w = zeros(n,1);
    cvx_begin quiet
      variable w(n);  
      maximize(mu' * w - k * norm(sqrt_Omega * w) - lambda/2 * w' * Sigma * w);
      subject to
        ones(1,n)*w == 1;
        w >= 0;
    cvx_end
    w0 = w;

    w = zeros(n,1);
    fun = @(w)objectivefcn1(w, n, mu, xi, D, lambda, Sigma);
        lb = zeros(n,1);
        Aeq = e';
        beq = 1;
        A = [];
        b = [];
        ub = [];
        nonlcon = []; 
        options = optimoptions('fmincon','Algorithm','sqp');
    [w,fval,exitflag,output] = fmincon(fun,w0,A,b,Aeq,beq,lb,ub, nonlcon,options);

    if exitflag == 0 || exitflag == -1 || exitflag == -2
        RObox2.W = [RObox2.W w];
        RObox2.ExpRet = [RObox2.ExpRet   NaN];
        RObox2.ExpRet_original = [RObox2.ExpRet_original  NaN];
        RObox2.stdev = [RObox2.stdev NaN];
        RObox2.exitflag = [RObox2.exitflag exitflag];
    else
        RObox2.W = [RObox2.W w];
        RObox2.ExpRet = [RObox2.ExpRet   -fval + lambda/2 * w' * Sigma * w];
        RObox2.ExpRet_original = [RObox2.ExpRet_original mu' * w];
        RObox2.stdev = [RObox2.stdev sqrt(w' * Sigma * w)];
        RObox2.exitflag = [RObox2.exitflag exitflag];
    end
end

% RObox2
w = zeros(n,1);
cvx_begin quiet
  variable w(n);  
  maximize(mu' * w - k * norm(sqrt_Omega * w) - Lambda/2 * w' * Sigma * w);
  subject to
    ones(1,n)*w == 1;
    w >= 0;
cvx_end
w0 = w;

w = zeros(n,1);
fun = @(w)objectivefcn1(w, n, mu, xi, D, lambda, Sigma);
    lb = zeros(n,1);
    Aeq = e';
    beq = 1;
    A = [];
    b = [];
    ub = [];
    nonlcon = []; 
    options = optimoptions('fmincon','Algorithm','sqp');
[w,fval,exitflag,output] = fmincon(fun,w0,A,b,Aeq,beq,lb,ub, nonlcon,options);

RObox2.w = w;

RObox2.upper = max(true_mu_box_set' * RObox2.w) 
RObox2.lower = min(true_mu_box_set' * RObox2.w) 
RObox2.interval = max(true_mu_box_set' * RObox2.w) - min(true_mu_box_set' * RObox2.w) 

%% calculate the range of expected utility
cvx_begin
variable w(n) ; 
maximize(mu' * w - Lambda/2*w'*Sigma*w);
  subject to
    ones(1,n)*w == 1;
    w >= 0;
cvx_end
MVO.w = w;
MVO.exp = mu' * w;

MVO.upperquad = max(true_mu_quad_set' * MVO.w) 
MVO.lowerquad = min(true_mu_quad_set' * MVO.w) 
MVO.intervalquad = max(true_mu_quad_set' * MVO.w) - min(true_mu_quad_set' * MVO.w) 

MVO.upperquadutility = max(true_mu_quad_set' * MVO.w - Lambda/2*MVO.w'*Sigma*MVO.w) 
MVO.lowerquadutility = min(true_mu_quad_set' * MVO.w - Lambda/2*MVO.w'*Sigma*MVO.w)
MVO.intervalquadutility = MVO.upperquadutility - MVO.lowerquadutility

MVO.upperbox = max(true_mu_box_set' * MVO.w) 
MVO.lowerbox = min(true_mu_box_set' * MVO.w) 
MVO.intervalbox = max(true_mu_box_set' * MVO.w) - min(true_mu_box_set' * MVO.w) 

MVO.upperboxutility = max(true_mu_box_set' * MVO.w - Lambda/2*MVO.w'*Sigma*MVO.w) 
MVO.lowerboxutility = min(true_mu_box_set' * MVO.w - Lambda/2*MVO.w'*Sigma*MVO.w)  
MVO.intervalboxutility = MVO.upperbox - MVO.lowerbox
% ROquad1
cvx_begin quiet
  variable w(n);  
  maximize(mu' * w - k * norm(sqrt_Omega * w) - Lambda/2 * w' * Sigma * w);
  subject to
    ones(1,n)*w == 1;
    w >= 0;
cvx_end
ROquad1.w = w;

ROquad1.upper = max(true_mu_quad_set' * ROquad1.w) 
ROquad1.lower = min(true_mu_quad_set' * ROquad1.w) 
ROquad1.interval = max(true_mu_quad_set' * ROquad1.w) - min(true_mu_quad_set' * ROquad1.w) 

ROquad1.upperutility = max(true_mu_quad_set' * ROquad1.w - Lambda/2 * ROquad1.w' * Sigma * ROquad1.w) 
ROquad1.lowerutility = min(true_mu_quad_set' * ROquad1.w - Lambda/2 * ROquad1.w' * Sigma * ROquad1.w) 
ROquad1.intervalutility = ROquad1.upperutility - ROquad1.lowerutility;

% ROquad2
cvx_begin quiet
  variable w(n);  
  sqrt_Omega_modified = Omega - 1/(e'*D*Omega*D'*e)*Omega*D'*e*e'*D*Sigma;
  maximize(mu' * w - k * norm(sqrt_Omega_modified * w) - Lambda/2 * w' * Sigma * w);
  subject to
    ones(1,n)*w == 1;
    w >= 0;
cvx_end
ROquad2.w = w;

ROquad2.upper = max(true_mu_quad_set' * ROquad2.w) 
ROquad2.lower = min(true_mu_quad_set' * ROquad2.w) 
ROquad2.interval = max(true_mu_quad_set' * ROquad2.w) - min(true_mu_quad_set' * ROquad2.w) 

ROquad2.upperutility = max(true_mu_quad_set' * ROquad2.w - Lambda/2 * ROquad2.w' * Sigma * ROquad2.w) 
ROquad2.lowerutility = min(true_mu_quad_set' * ROquad2.w - Lambda/2 * ROquad2.w' * Sigma * ROquad2.w) 
ROquad2.intervalutility = ROquad2.upperutility - ROquad2.lowerutility


% MVO3
cvx_begin
variable w(n) ; 
maximize(mu' * (w-benchmark) - Lambda3/2*(w-benchmark)'*Sigma*(w-benchmark));
  subject to
    ones(1,n)*w == 1;
    w >= 0;
cvx_end
MVO3.w = w;

MVO3.upperquad = max(true_mu_quad_set' * MVO3.w) 
MVO3.lowerquad = min(true_mu_quad_set' * MVO3.w) 
MVO3.intervalquad = max(true_mu_quad_set' * MVO3.w) - min(true_mu_quad_set' * MVO3.w) 

MVO3.upperquadutility = max(true_mu_quad_set' * (MVO3.w-benchmark) - Lambda3/2*(MVO3.w-benchmark)'*Sigma*(MVO3.w-benchmark))
MVO3.lowerquadutility = min(true_mu_quad_set' * (MVO3.w-benchmark) - Lambda3/2*(MVO3.w-benchmark)'*Sigma*(MVO3.w-benchmark))
MVO3.intervalquadutility = MVO3.upperquadutility - MVO3.lowerquadutility


MVO3.upperbox = max(true_mu_box_set' * MVO3.w) 
MVO3.lowerbox = min(true_mu_box_set' * MVO3.w) 
MVO3.intervalbox = max(true_mu_box_set' * MVO3.w) - min(true_mu_box_set' * MVO3.w) 

MVO3.upperboxutility = max(true_mu_box_set' * (MVO3.w-benchmark) - Lambda3/2*(MVO3.w-benchmark)'*Sigma*(MVO3.w-benchmark))
MVO3.lowerboxutility = min(true_mu_box_set' * (MVO3.w-benchmark) - Lambda3/2*(MVO3.w-benchmark)'*Sigma*(MVO3.w-benchmark))
MVO3.intervalboxutility = MVO3.upperboxutility - MVO3.lowerboxutility


% ROquad3
cvx_begin quiet
  variable w(n);  
  maximize(mu' * (w-benchmark) - k * norm(sqrt_Omega * (w-benchmark)) - Lambda3/2 * (w-benchmark)' * Sigma * (w-benchmark));
  subject to
    ones(1,n)*w == 1;
    w >= 0;
cvx_end
ROquad3.w = w;

ROquad3.upper = max(true_mu_quad_set' * ROquad3.w)
ROquad3.lower = min(true_mu_quad_set' * ROquad3.w) 
ROquad3.interval = max(true_mu_quad_set' * ROquad3.w) - min(true_mu_quad_set' * ROquad3.w) 

ROquad3.upperutility = max(true_mu_quad_set' * (ROquad3.w-benchmark)  - Lambda3/2 * (ROquad3.w-benchmark)' * Sigma * (ROquad3.w-benchmark))
ROquad3.lowerutility = min(true_mu_quad_set' * (ROquad3.w-benchmark)  - Lambda3/2 * (ROquad3.w-benchmark)' * Sigma * (ROquad3.w-benchmark))
ROquad3.intervalutility = ROquad3.upperutility - ROquad3.lowerutility


% RObox1
cvx_begin
  variable w(n);  
  maximize((mu' * w) - xi * sum(abs(w)) - Lambda/2 * w' * Sigma * w);
  subject to
    ones(1,n)*w == 1;
    w >= 0;
cvx_end
RObox1.w = w;

RObox1.upper = max(true_mu_box_set' * RObox1.w) 
RObox1.lower = min(true_mu_box_set' * RObox1.w) 
RObox1.interval = max(true_mu_box_set' * RObox1.w) - min(true_mu_box_set' * RObox1.w) 

RObox1.upperutility = max(true_mu_box_set' * RObox1.w - Lambda/2 * RObox1.w' * Sigma * RObox1.w) 
RObox1.lowerutility = min(true_mu_box_set' * RObox1.w - Lambda/2 * RObox1.w' * Sigma * RObox1.w) 
RObox1.intervalutility = RObox1.upperutility - RObox1.lowerutility


% RObox3
cvx_begin quiet
  variable w(n);  
  maximize ((mu' * (w-benchmark)) - xi * sum(abs((w-benchmark))) - Lambda3/2 * (w-benchmark)' * Sigma * (w-benchmark));
  subject to
    ones(1,n)*w == 1;
    w >= 0;
cvx_end
RObox3.w = w;

RObox3.upper = max(true_mu_box_set' * RObox3.w) 
RObox3.lower = min(true_mu_box_set' * RObox3.w) 
RObox3.interval = max(true_mu_box_set' * RObox3.w) - min(true_mu_box_set' * RObox3.w) 

RObox3.upperutility = max(true_mu_box_set' * (RObox3.w-benchmark) - Lambda3/2 * (RObox3.w-benchmark)' * Sigma * (RObox3.w-benchmark))
RObox3.lowerutility = min(true_mu_box_set' * (RObox3.w-benchmark) - Lambda3/2 * (RObox3.w-benchmark)' * Sigma * (RObox3.w-benchmark)) 
RObox3.intervalutility = RObox3.upperutility - RObox3.lowerutility

% RObox2
w = zeros(n,1);
cvx_begin quiet
  variable w(n);  
  maximize(mu' * w - k * norm(sqrt_Omega * w) - Lambda/2 * w' * Sigma * w);
  subject to
    ones(1,n)*w == 1;
    w >= 0;
cvx_end
w0 = w;

w = zeros(n,1);
fun = @(w)objectivefcn1(w, n, mu, xi, D, lambda, Sigma);
    lb = zeros(n,1);
    Aeq = e';
    beq = 1;
    A = [];
    b = [];
    ub = [];
    nonlcon = []; 
    options = optimoptions('fmincon','Algorithm','sqp');
[w,fval,exitflag,output] = fmincon(fun,w0,A,b,Aeq,beq,lb,ub, nonlcon,options);

RObox2.w = w;

RObox2.upper = max(true_mu_box_set' * RObox2.w) 
RObox2.lower = min(true_mu_box_set' * RObox2.w) 
RObox2.interval = max(true_mu_box_set' * RObox2.w) - min(true_mu_box_set' * RObox2.w) 

RObox2.upperutility = max(true_mu_box_set' * RObox2.w - Lambda/2 * w' * Sigma * w)
RObox2.lowertility = min(true_mu_box_set' * RObox2.w - Lambda/2 * w' * Sigma * w) 
RObox2.intervaltility = RObox2.upperutility - RObox2.lowertility 

%% plots
% MVO ROquad1 ROquad2
figure(1);
subplot(1,3,1);
f = area(lambda_set(1:61),MVO.W(:,1:61)');
title('MVO Composition');
xlabel('Risk Aversion');
ylabel('Optimal Weight of Each Asset');
legend('1','2','3','4','5','6','7','8','9','10');
xlim([0,50]);
ylim([0 1]);
alpha(0.48);
for i = 1:10
    f(i).EdgeColor = 'none';
end
subplot(1,3,2);
f = area(lambda_set(1:61),ROquad1.W(:,1:61)');
title('ROquad1 Composition');
xlabel('Risk Aversion');
ylabel('Optimal Weight of Each Asset');
legend('1','2','3','4','5','6','7','8','9','10');
xlim([0,50]);
ylim([0 1]);
alpha(0.48);
for i = 1:10
    f(i).EdgeColor = 'none';
end
subplot(1,3,3);
f = area(lambda_set(1:61),ROquad2.W(:,1:61)');
title('ROquad2 Composition');
xlabel('Risk Aversion');
ylabel('Optimal Weight of Each Asset');
legend('1','2','3','4','5','6','7','8','9','10');
xlim([0,50]);
ylim([0 1]);
alpha(0.48);
for i = 1:10
    f(i).EdgeColor = 'none';
end

figure(2);
plot(MVO.stdev,MVO.ExpRet,'-', 'LineWidth',2, ...
        'Color',[0.439215686274510   0.662745098039216   0.858823529411765]);hold on;
plot(ROquad1.stdev,ROquad1.ExpRet,'-', 'LineWidth',2, ...
        'Color',[0.890196078431372   0.380392156862745   0.19607843137254]);hold on;
plot(ROquad2.stdev,ROquad2.ExpRet,'-','LineWidth',2, ...
        'Color',[0.929411764705882   0.643137254901961   0.392156862745098]);
title('Efficient Frontier');
%ylim([0.7,1.2])
xlabel('Standard Deviation of Portfolio');
ylabel('Expected Return of Portfolio');
legend("MVO","ROquad1","ROquad2");

% MVO3 ROquad3
figure(3);
subplot(1,2,1);
f = area(lambda_set(1:61),MVO3.W(:,1:61)');
title('MVO3 Composition');
xlabel('Risk Aversion');
ylabel('Optimal Weight of Each Asset');
legend('1','2','3','4','5','6','7','8','9','10');
xlim([0,50]);
ylim([0 1]);
alpha(0.48);
for i = 1:10
    f(i).EdgeColor = 'none';
end
subplot(1,2,2);
f = area(lambda_set(1:61),ROquad3.W(:,1:61)');
title('ROquad3 Composition');
xlabel('Risk Aversion');
ylabel('Optimal Weight of Each Asset');
legend('1','2','3','4','5','6','7','8','9','10');
xlim([0,50]);
ylim([0 1]);
alpha(0.48);
for i = 1:10
    f(i).EdgeColor = 'none';
end

figure(4);
plot(MVO3.stdev,MVO3.ExpRet,'-', 'LineWidth',2, ...
        'Color',[0.439215686274510   0.662745098039216   0.858823529411765]);hold on;
plot(ROquad3.stdev,ROquad3.ExpRet,'-', 'LineWidth',2, ...
        'Color',[0.929411764705882   0.643137254901961   0.392156862745098]);
title('Efficient Frontier');
xlabel('Standard Deviation of Portfolio');
ylabel('Expected Return of Portfolio');
legend("MVO3","ROquad3");

% MVO RObox1 RObox2
figure(5);
subplot(1,3,1);
f = area(lambda_set(1:61),MVO.W(:,1:61)');
title('MVO Composition');
xlabel('Risk Aversion');
ylabel('Optimal Weight of Each Asset');
legend('1','2','3','4','5','6','7','8','9','10');
xlim([0,50]);
ylim([0 1]);
alpha(0.48);
for i = 1:10
    f(i).EdgeColor = 'none';
end
subplot(1,3,2);
f = area(lambda_set(1:61),RObox1.W(:,1:61)');
title('RObox1 Composition');
xlabel('Risk Aversion');
ylabel('Optimal Weight of Each Asset');
legend('1','2','3','4','5','6','7','8','9','10');
xlim([0,50]);
ylim([0 1]);
alpha(0.48);
for i = 1:10
    f(i).EdgeColor = 'none';
end
subplot(1,3,3);
f = area(lambda_set(1:61),RObox2.W(:,1:61)');
title('RObox2 Composition');
xlabel('Risk Aversion');
ylabel('Optimal Weight of Each Asset');
legend('1','2','3','4','5','6','7','8','9','10');
xlim([0,50]);
ylim([0 1]);
alpha(0.48);
for i = 1:10
    f(i).EdgeColor = 'none';
end

figure(6);
plot(MVO.stdev,MVO.ExpRet,'-', 'LineWidth',2, ...
        'Color',[0.439215686274510   0.662745098039216   0.858823529411765]);hold on;
plot(RObox1.stdev,RObox1.ExpRet,'-', 'LineWidth',2, ...
        'Color',[0.890196078431372   0.380392156862745   0.19607843137254]);hold on;
plot(RObox2.stdev,RObox2.ExpRet,'-','LineWidth',2, ...
        'Color',[0.929411764705882   0.643137254901961   0.392156862745098]);
title('Efficient Frontier');
xlabel('Standard Deviation of Portfolio');
ylabel('Expected Return of Portfolio');
legend("MVO","RObox1","RObox2");

% MVO3 RObox3
figure(7);
subplot(1,2,1);
f = area(lambda_set(1:61),MVO3.W(:,1:61)');
title('MVO3 Composition');
xlabel('Risk Aversion');
ylabel('Optimal Weight of Each Asset');
legend('1','2','3','4','5','6','7','8','9','10');
xlim([0,50]);
ylim([0 1]);
alpha(0.48);
for i = 1:10
    f(i).EdgeColor = 'none';
end
subplot(1,2,2);
f = area(lambda_set(1:61),RObox3.W(:,1:61)');
title('RObox3 Composition');
xlabel('Risk Aversion');
ylabel('Optimal Weight of Each Asset');
legend('1','2','3','4','5','6','7','8','9','10');
xlim([0,50]);
ylim([0 1]);
alpha(0.48);
for i = 1:10
    f(i).EdgeColor = 'none';
end


figure(8);
plot(MVO3.stdev,MVO3.ExpRet,'-', 'LineWidth',2, ...
        'Color',[0.439215686274510   0.662745098039216   0.858823529411765]);hold on;
plot(RObox3.stdev,RObox3.ExpRet,'-', 'LineWidth',2, ...
        'Color',[0.929411764705882   0.643137254901961   0.392156862745098]);
title('Efficient Frontier');
xlabel('Standard Deviation of Portfolio');
ylabel('Expected Return of Portfolio');
legend("MVO3","RObox3");


% comparision MVO ROquad1 RObox1
figure(9);
plot(MVO.stdev,MVO.ExpRet,'-', 'LineWidth',2, ...
        'Color',[0.439215686274510   0.662745098039216   0.858823529411765]);hold on;
plot(ROquad1.stdev,ROquad1.ExpRet,'-', 'LineWidth',2, ...
        'Color',[0.890196078431372   0.380392156862745   0.19607843137254]);hold on;
plot(RObox1.stdev,RObox1.ExpRet,'-','LineWidth',2, ...
        'Color',[0.929411764705882   0.643137254901961   0.392156862745098]);
title('Efficient Frontier');
xlabel('Standard deviation of portfolio');
ylabel('Expected return of portfolio');


% comparision MVO ROquad2 RObox2
figure(10);
plot(MVO.stdev,MVO.ExpRet,'-', 'LineWidth',2, ...
        'Color',[0.439215686274510   0.662745098039216   0.858823529411765]);hold on;
plot(ROquad2.stdev,ROquad2.ExpRet,'-', 'LineWidth',2, ...
        'Color',[0.890196078431372   0.380392156862745   0.19607843137254]);hold on;
plot(RObox2.stdev,RObox2.ExpRet,'-','LineWidth',2, ...
        'Color',[0.929411764705882   0.643137254901961   0.392156862745098]);
title('Efficient Frontier');
xlabel('Standard deviation of portfolio');
ylabel('Expected return of portfolio');

% comparision MVO3 ROquad3 RObox3
figure(11);
plot(MVO3.stdev,MVO3.ExpRet,'-', 'LineWidth',2, ...
        'Color',[0.439215686274510   0.662745098039216   0.858823529411765]);hold on;
plot(ROquad3.stdev,ROquad3.ExpRet,'-', 'LineWidth',2, ...
        'Color',[0.890196078431372   0.380392156862745   0.19607843137254]);hold on;
plot(RObox3.stdev,RObox3.ExpRet,'-','LineWidth',2, ...
        'Color',[0.929411764705882   0.643137254901961   0.392156862745098]);
title('Efficient Frontier');
xlabel('Standard deviation of portfolio');
ylabel('Expected return of portfolio');
