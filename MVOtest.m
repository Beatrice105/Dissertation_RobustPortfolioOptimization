%% Data Simulation
rng(30);
marketData = 100 - 4  + 8 * normrnd(0,1,12,4); % R = normrnd(MU,SIGMA,m); [6,14]

% transfer marketData to ret
m = size(marketData,1); % 12
n = size(marketData,2); % 3
returns = zeros(m,n);
for j=1:n
    for i=1:(m-1)
        returns(m-i,j) = marketData(i,j)/marketData(i+1,j);
    end
end
returns;
returns(m,:)=[];
avg = mean(returns);
mu = geo_mean(returns)';
mu

avg = mean(returns);
mu = geo_mean(returns)';
% calculate covariance matrix (Sigma)
Sigma = zeros(n,n);
T = m - 1;
for i=1:n
    for j=1:i
        Sigma(i,j) = ((returns(:,i) - avg(i))' * (returns(:,j) - avg(j)))/T ;
    end 
end
% Now flip the matrix and add
Sigma = Sigma + triu(Sigma',1);
Sigma
% corr = corrcoef(returns);

%% Markwitz Model
mu = geo_mean(returns)';
W = [];
ExpRet = [];
stdev = [];
lambda_set = [(0:0.5:20)];
for lambda = lambda_set 
    cvx_begin quiet
      variable w(n);  
      maximize(mu'*w - lambda/2*w'*Sigma*w);
      subject to
        ones(1,n)*w == 1;
        w >= 0;
    cvx_end
    W = [W w];
    ExpRet = [ExpRet mu'*w]; % true
    stdev = [stdev sqrt(w' * Sigma * w)];
end
W_mvo0 = W;
ExpRet_mvo0 = ExpRet;
stdev_mvo0 = stdev;
% calculate lambda of CML
% get optimal SR ==> Lambda
[Lambda, num] = max(ExpRet_mvo0 ./ stdev_mvo0);

%% sensitivity
mu_1 = [mu(1:2)*(1.01) ; mu(3:4)*(0.99)];
mu_2 = [mu(1:2)*(0.99) ; mu(3:4)*(1.01)];

mu = mu_1;
W = [];
ExpRet = [];
stdev = [];
ExpRet_actual = [];
for lambda = lambda_set 
    cvx_begin quiet
      variable w(n);  
      maximize(mu'*w - lambda/2*w'*Sigma*w);
      subject to
        ones(1,n)*w == 1;
        w >= 0;
    cvx_end
    W = [W w];
    ExpRet = [ExpRet mu'*w]; % estimated
    stdev = [stdev sqrt(w' * Sigma * w)];

    ExpRet_actual = [ExpRet_actual geo_mean(returns)*w]; % actual
end
W_mvo1 = W;
ExpRet_mvo1 = ExpRet;
stdev_mvo1 = stdev;
ExpRet_mvo1_actual = ExpRet_actual;

mu = mu_2;
W = [];
ExpRet = [];
stdev = [];
ExpRet_actual = [];
for lambda = lambda_set 
    cvx_begin quiet
      variable w(n);  
      maximize(mu'*w - lambda/2*w'*Sigma*w);
      subject to
        ones(1,n)*w == 1;
        w >= 0;
    cvx_end
    W = [W w];
    ExpRet = [ExpRet mu'*w]; % estimated
    stdev = [stdev sqrt(w' * Sigma * w)];

    ExpRet_actual = [ExpRet_actual geo_mean(returns)*w]; % actual
end
W_mvo2 = W;
ExpRet_mvo2 = ExpRet;
stdev_mvo2 = stdev;
ExpRet_mvo2_actual = ExpRet_actual;

figure(1);
t = tiledlayout(1,3);

nexttile
f2 = area(lambda_set,W_mvo0');
title('Portfolio O Composition');
xlabel('Risk Aversion');
ylabel('Optimal Weight of Each Asset');
legend('Asset 1', 'Asset 2', 'Asset 3','Asset 4');
alpha(0.5);
f2(1).FaceColor = [0.901960784313726   0.509803921568627   0.172549019607843];
f2(2).FaceColor = [0.929411764705882   0.643137254901961   0.392156862745098];
f2(3).FaceColor = [0.439215686274510   0.662745098039216   0.858823529411765];
f2(4).FaceColor = [0.129411764705882   0.494117647058824   0.811764705882353];
f2(1).EdgeColor = [0.901960784313726   0.509803921568627   0.172549019607843];
f2(2).EdgeColor = [0.929411764705882   0.643137254901961   0.392156862745098];
f2(3).EdgeColor = [0.439215686274510   0.662745098039216   0.858823529411765];
f2(4).EdgeColor = [0.129411764705882   0.494117647058824   0.811764705882353];

nexttile
f4 = area(lambda_set,W_mvo1');
title('Portfolio A Composition');
xlabel('Risk Aversion');
ylabel('Optimal Weight of Each Asset');
legend('Asset 1', 'Asset 2', 'Asset 3','Asset 4');
alpha(0.5);
f4(1).FaceColor = [0.901960784313726   0.509803921568627   0.172549019607843];
f4(2).FaceColor = [0.929411764705882   0.643137254901961   0.392156862745098];
f4(3).FaceColor = [0.439215686274510   0.662745098039216   0.858823529411765];
f4(4).FaceColor = [0.129411764705882   0.494117647058824   0.811764705882353];
f4(1).EdgeColor = [0.901960784313726   0.509803921568627   0.172549019607843];
f4(2).EdgeColor = [0.929411764705882   0.643137254901961   0.392156862745098];
f4(3).EdgeColor = [0.439215686274510   0.662745098039216   0.858823529411765];
f4(4).EdgeColor = [0.129411764705882   0.494117647058824   0.811764705882353];

nexttile
f6 = area(lambda_set,W_mvo2');
title('Portfolio B Composition');
xlabel('Risk Aversion');
ylabel('Optimal Weight of Each Asset');
legend('Asset 1', 'Asset 2', 'Asset 3','Asset 4');
alpha(0.5);
f6(1).FaceColor = [0.901960784313726   0.509803921568627   0.172549019607843];
f6(2).FaceColor = [0.929411764705882   0.643137254901961   0.392156862745098];
f6(3).FaceColor = [0.439215686274510   0.662745098039216   0.858823529411765];
f6(4).FaceColor = [0.129411764705882   0.494117647058824   0.811764705882353];
f6(1).EdgeColor = [0.901960784313726   0.509803921568627   0.172549019607843];
f6(2).EdgeColor = [0.929411764705882   0.643137254901961   0.392156862745098];
f6(3).EdgeColor = [0.439215686274510   0.662745098039216   0.858823529411765];
f6(4).EdgeColor = [0.129411764705882   0.494117647058824   0.811764705882353];

% efficient frontier
figure(2);
t = tiledlayout(1,2);

nexttile
plot(stdev_mvo0,ExpRet_mvo0,'--', 'LineWidth',2, ...
        'Color',[0.439215686274510   0.662745098039216   0.858823529411765]);hold on;
plot(stdev_mvo1,ExpRet_mvo1,'-', 'LineWidth',2, ...
        'Color',[0.890196078431372   0.380392156862745   0.19607843137254]);hold on;
plot(stdev_mvo1,ExpRet_mvo1_actual,'-','LineWidth',2, ...
        'Color',[0.929411764705882   0.643137254901961   0.392156862745098]);
title('Efficient Frontier');
xlabel('Standard Deviation of Portfolio');
ylabel('Expected Return of Portfolio');
legend("True Frontier","Estimated MVO Frontier","Actual MVO Frontier");

nexttile
plot(stdev_mvo0,ExpRet_mvo0,'--', 'LineWidth',2, ...
        'Color',[0.439215686274510   0.662745098039216   0.858823529411765]);hold on;
plot(stdev_mvo2,ExpRet_mvo2,'-', 'LineWidth',2, ...
        'Color',[0.890196078431372   0.380392156862745   0.19607843137254]);hold on;
plot(stdev_mvo2,ExpRet_mvo2_actual,'-','LineWidth',2, ...
        'Color',[0.929411764705882   0.643137254901961   0.392156862745098]);
title('Efficient Frontier');
xlabel('Standard Deviation of Portfolio');
ylabel('Expected Return of Portfolio');
legend("True Frontier","Estimated MVO Frontier","Actual MVO Frontier");

[Lambda, num] = max(ExpRet_mvo0 ./ stdev_mvo0);
Lambda

Lambda = 5

mu = geo_mean(returns)';
cvx_begin quiet
  variable w(n);  
  maximize(mu'*w - Lambda/2*w'*Sigma*w);
  subject to
    ones(1,n)*w == 1;
    w >= 0;
cvx_end
w0 = w;


mu = mu_1;
cvx_begin quiet
  variable w(n);  
  maximize(mu'*w - Lambda/2*w'*Sigma*w);
  subject to
    ones(1,n)*w == 1;
    w >= 0;
cvx_end
w1 = w;

mu = mu_2;
cvx_begin quiet
  variable w(n);  
  maximize(mu'*w - Lambda/2*w'*Sigma*w);
  subject to
    ones(1,n)*w == 1;
    w >= 0;
cvx_end
w2 = w;

[w0 w1 w2];
x=round([w0 w1 w2],4) % 有区别 说明Lambda不同

[geo_mean(returns)*w0, mu_1' * w1, mu_2' * w2]


% 法四 **********************************
mu = geo_mean(returns)';
mu_1 = [mu(1:2)*(1.01) ; mu(3:4)*(0.99)];
mu_2 = [mu(1:2)*(0.99) ; mu(3:4)*(1.01)];
mu_1
mu_best = mu_1;
mu_worst = [mu_1(1:2)-0.05 ; mu_1(3:4)+0.05];


Lambda = 5
mu = mu_best;
cvx_begin quiet
  variable w(n);  
  maximize(mu'*w - Lambda/2*w'*Sigma*w);
  subject to
    ones(1,n)*w == 1;
    w >= 0;
cvx_end
w1 = w;

% best
cvx_begin
  variable mu(4); 
  variable a;
  variable b;
  maximize(mu'*w1);
  subject to
  for i = 1:4
    % convex combinition
    mu == a * mu_best + b * mu_worst;
    a + b == 1;
    a>=0 ;
    b>=0 ;
  end
cvx_end
round(mu,6) == round(mu_best,6)

cvx_begin
  variable mu(4); 
  variable a;
  variable b;
  minimize(mu'*w1);
  subject to
  for i = 1:4
    % convex combinition
    mu == a * mu_best + b * mu_worst;
    a + b == 1;
    a>=0 ;
    b>=0 ;
  end
cvx_end
round(mu,6) == round(mu_worst,6)
% maximize --- mu_best
% minimize --- mu_worst

mu = mu_best;
cvx_begin quiet
  variable w(n);  
  maximize(mu'*w - Lambda/2*w'*Sigma*w);
  subject to
    ones(1,n)*w == 1;
    w >= 0;
cvx_end
w_best = w;


mu = mu_worst;
cvx_begin quiet
  variable w(n);  
  maximize(mu'*w - Lambda/2*w'*Sigma*w);
  subject to
    ones(1,n)*w == 1;
    w >= 0;
cvx_end
w_worst = w;

w_best
w_worst
% best case
mu_best' * w_best
mu_best' * w_worst

% worse case
mu_worst' * w_best
mu_worst' * w_worst

mu_best' * w_best-mu_best' * w_worst
mu_worst' * w_best-mu_worst' * w_worst



%% MVO with other constraints
w(3)+w(4) <= 0.6;
w(2)+w(3) >= 0.5;
w(3)+w(4) <= 0.6;
w(1) >= 0.2;
w(4) >= 0.2;

%% Markwitz Model
mu = geo_mean(returns)';
W = [];
ExpRet = [];
stdev = [];
lambda_set = [(0:0.5:20)];
for lambda = lambda_set 
    cvx_begin quiet
      variable w(n);  
      maximize(mu'*w - lambda/2*w'*Sigma*w);
      subject to
        ones(1,n)*w == 1;
        w >= 0;
        w(3)+w(4) <= 0.6;
        w(2)+w(3) >= 0.5;
        w(3)+w(4) <= 0.6;
        w(1) >= 0.2;
        w(4) >= 0.2;
    cvx_end
    W = [W w];
    ExpRet = [ExpRet mu'*w];
    stdev = [stdev sqrt(w' * Sigma * w)];
end
W_mvo0 = W;
ExpRet_mvo0 = ExpRet;
stdev_mvo0 = stdev;
% calculate lambda of CML
% get optimal SR ==> Lambda
[Lambda, num] = max(ExpRet_mvo0 ./ stdev_mvo0);
% 161.7181297577673

%% sensitivity
mu_1 = [mu(1:2)*(1.01) ; mu(3:4)*(0.99)]
mu_2 = [mu(1:2)*(0.99) ; mu(3:4)*(1.01)]

mu = mu_1
W = [];
ExpRet = [];
stdev = [];
lambda_set = [(0:0.5:20)];
for lambda = lambda_set 
    cvx_begin quiet
      variable w(n);  
      maximize(mu'*w - lambda/2*w'*Sigma*w);
      subject to
        ones(1,n)*w == 1;
        w >= 0;
        w(3)+w(4) <= 0.6;
        w(2)+w(3) >= 0.5;
        w(3)+w(4) <= 0.6;
        w(1) >= 0.2;
        w(4) >= 0.2;
    cvx_end
    W = [W w];
    ExpRet = [ExpRet mu'*w];
    stdev = [stdev sqrt(w' * Sigma * w)];
end
W_mvo1 = W;
ExpRet_mvo1 = ExpRet;
stdev_mvo1 = stdev;

mu = mu_2
W = [];
ExpRet = [];
stdev = [];
lambda_set = [(0:0.5:20)];
for lambda = lambda_set 
    cvx_begin quiet
      variable w(n);  
      maximize(mu'*w - lambda/2*w'*Sigma*w);
      subject to
        ones(1,n)*w == 1;
        w >= 0;
        w(3)+w(4) <= 0.6;
        w(2)+w(3) >= 0.5;
        w(3)+w(4) <= 0.6;
        w(1) >= 0.2;
        w(4) >= 0.2;
    cvx_end
    W = [W w];
    ExpRet = [ExpRet mu'*w];
    stdev = [stdev sqrt(w' * Sigma * w)];
end
W_mvo2 = W;
ExpRet_mvo2 = ExpRet;
stdev_mvo2 = stdev;

figure(7);
t = tiledlayout(1,3);

nexttile
f2 = area(lambda_set,W_mvo0');
title('Portfolio O Composition');
xlabel('Risk Aversion');
ylabel('Optimal Weight of Each Asset');
legend('Asset 1', 'Asset 2', 'Asset 3','Asset 4');
alpha(0.5);
f2(1).FaceColor = [0.901960784313726   0.509803921568627   0.172549019607843];
f2(2).FaceColor = [0.929411764705882   0.643137254901961   0.392156862745098];
f2(3).FaceColor = [0.439215686274510   0.662745098039216   0.858823529411765];
f2(4).FaceColor = [0.129411764705882   0.494117647058824   0.811764705882353];
f2(1).EdgeColor = [0.901960784313726   0.509803921568627   0.172549019607843];
f2(2).EdgeColor = [0.929411764705882   0.643137254901961   0.392156862745098];
f2(3).EdgeColor = [0.439215686274510   0.662745098039216   0.858823529411765];
f2(4).EdgeColor = [0.129411764705882   0.494117647058824   0.811764705882353];

nexttile
f4 = area(lambda_set,W_mvo1');
title('Portfolio A Composition');
xlabel('Risk Aversion');
ylabel('Optimal Weight of Each Asset');
legend('Asset 1', 'Asset 2', 'Asset 3','Asset 4');
alpha(0.5);
f4(1).FaceColor = [0.901960784313726   0.509803921568627   0.172549019607843];
f4(2).FaceColor = [0.929411764705882   0.643137254901961   0.392156862745098];
f4(3).FaceColor = [0.439215686274510   0.662745098039216   0.858823529411765];
f4(4).FaceColor = [0.129411764705882   0.494117647058824   0.811764705882353];
f4(1).EdgeColor = [0.901960784313726   0.509803921568627   0.172549019607843];
f4(2).EdgeColor = [0.929411764705882   0.643137254901961   0.392156862745098];
f4(3).EdgeColor = [0.439215686274510   0.662745098039216   0.858823529411765];
f4(4).EdgeColor = [0.129411764705882   0.494117647058824   0.811764705882353];

nexttile
f6 = area(lambda_set,W_mvo2');
title('Portfolio B Composition');
xlabel('Risk Aversion');
ylabel('Optimal Weight of Each Asset');
legend('Asset 1', 'Asset 2', 'Asset 3','Asset 4');
alpha(0.5);
f6(1).FaceColor = [0.901960784313726   0.509803921568627   0.172549019607843];
f6(2).FaceColor = [0.929411764705882   0.643137254901961   0.392156862745098];
f6(3).FaceColor = [0.439215686274510   0.662745098039216   0.858823529411765];
f6(4).FaceColor = [0.129411764705882   0.494117647058824   0.811764705882353];
f6(1).EdgeColor = [0.901960784313726   0.509803921568627   0.172549019607843];
f6(2).EdgeColor = [0.929411764705882   0.643137254901961   0.392156862745098];
f6(3).EdgeColor = [0.439215686274510   0.662745098039216   0.858823529411765];
f6(4).EdgeColor = [0.129411764705882   0.494117647058824   0.811764705882353];




mu = geo_mean(returns)';
cvx_begin quiet
  variable w(n);  
  maximize(mu'*w - Lambda/2*w'*Sigma*w);
  subject to
    ones(1,n)*w == 1;
    w >= 0;
        w(3)+w(4) <= 0.6;
        w(2)+w(3) >= 0.5;
        w(3)+w(4) <= 0.6;
        w(1) >= 0.2;
        w(4) >= 0.2;
cvx_end
w0 = w;


mu_1 = [mu(1:2)*(1.01) ; mu(3:4)*(0.99)]
mu_2 = [mu(1:2)*(0.99) ; mu(3:4)*(1.01)]

mu = mu_1;
cvx_begin quiet
  variable w(n);  
  maximize(mu'*w - Lambda/2*w'*Sigma*w);
  subject to
    ones(1,n)*w == 1;
    w >= 0;
        w(3)+w(4) <= 0.6;
        w(2)+w(3) >= 0.5;
        w(3)+w(4) <= 0.6;
        w(1) >= 0.2;
        w(4) >= 0.2;
cvx_end
w1 = w;

mu = mu_2;
cvx_begin quiet
  variable w(n);  
  maximize(mu'*w - Lambda/2*w'*Sigma*w);
  subject to
    ones(1,n)*w == 1;
    w >= 0;
        w(3)+w(4) <= 0.6;
        w(2)+w(3) >= 0.5;
        w(3)+w(4) <= 0.6;
        w(1) >= 0.2;
        w(4) >= 0.2;
cvx_end
w2 = w;

[w0 w1 w2]
x=round([w0 w1 w2],4)



mu_worst = min(returns);
mu_worst * w1 
mu_worst * w2



mu_best = max(returns);
mu_best * w1 
mu_best * w2

(mu_worst * w1 < mu_worst * w2) == (mu_best * w1 < mu_best * w2)

