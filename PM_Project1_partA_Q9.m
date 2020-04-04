%*************************************************************************
%*************************************************************************

% 9) portfolio of 10 industries: Food, Smoke, Clths, Steel, Oil, Agric, 
%    Drugs, Txtls, Autos, Chips.

%*************************************************************************
%*************************************************************************


clear;clc;

%**********************************************

% 1) "mean-variance locus" graph (without the risk-free asset)

data= xlsread('10_Industry_Portfolios(Q9).xlsm');
dateline = data(:,1);
dataset=data(1:end,2:11)./100;
names ={'Agric';'Food';'Smoke';'Clths';'Drugs';'Txtls';'Steel';'Autos';'Oil';'Chips'};

%[T,M]=size(data1);
%T_begin=T-12*5;
%T_end=T;
%dataset=data1(T_begin:T_end,:);
R=0.02/12;
sigma=cov(dataset);
isigma=inv(sigma); % invert sigma(variance-covariance matrice)
z_bar=(mean(dataset))';
n=length(z_bar);

Vone=ones(length(z_bar),1); %vector of ones
A=Vone'*isigma*Vone;
B=Vone'*isigma*z_bar;
C=z_bar'*isigma*z_bar;
delta=A*C-B^2;

mu=NaN(101,1);
lamda=NaN(101,1);
gamma=NaN(101,1);
sigma_sqr=NaN(101,1);
sharpe=NaN(101,1);

k=101;
ul=0.09;
ll=-0.01;

for i=1:k;
   mu(i,1)=ll+(ul-ll)/(k-1)*(i-1);
   lamda(i,1)=(C-mu(i,1)*B)/delta;
   gamma(i,1)=(mu(i,1)*A-B)/delta;
   weight=(lamda(i,1)*isigma*Vone)+(gamma(i,1)*isigma*z_bar);
   sigma2(i,1)=((A*(mu(i,1)^2))-2*B*mu(i,1)+C)/delta;
   sigma_sqr(i,1)=sqrt(sigma2(i,1));
   sharpe(i,1)=(mu(i,1)-R)/sqrt(sigma2(i,1));
end;

mumin=B/A;
sigmamin=sqrt(((A*(mumin^2))-2*B*mumin+C)/delta);
ind=(mu>mumin); %Indicates efficient horizon
ind2=(mu<mumin); %Indicates locus below efficient horizon

% prepare data to draw the asymptotic line
sigma_sqr2=[sigma_sqr;0]; 
mu1=B/A+sigma_sqr2*sqrt(delta/A);
mu2=B/A-sigma_sqr2*sqrt(delta/A);

% Creat plot

figure(1)
f2=plot(sigma_sqr2,mu1,':',sigma_sqr2,mu2,':');
set(f2,'color','black');
hold on;

for i=1:n
    hold on
    plot(sqrt(sigma(i,i)),z_bar(i),'k^')
        text(sqrt(sigma(i,i)),z_bar(i),['  ' names{i}])
end

f1=plot(sigma_sqr(ind), mu(ind),'-',sigma_sqr(ind2),mu(ind2),'--',...
    sigmamin,mumin,'.');
%
set(f1(1:2),'linewidth',2);
set(f1(1:2),'color','blue');
set(f1(3),'markersize',40);
set(f1(3),'color','red');
%
xlabel('Standard deviation of return','fontsize',14);
ylabel('Expected return','fontsize',14);
set(gca,'xlim',[0.00,0.2]);
set(gca,'ylim',[-0.01,0.03]);
grid;
title({'Mean-variance locus ( without the risk-free asset )',...
    '10 industries choosed: Food, Smoke, Clths, Steel, Oil, Agric, Drugs, Txtls, Autos, Chips'},...
    'fontsize',14);
hold off;


%Optimization 
%fun = (1/2)*w'*sigma*w;
%[w,sigma]= fmincon(@fun,var,A,b,Aeq,beq,lb,ub);

%**********************************************

% 2) "mean-variance locus" graph (with the risk-free asset)

mu_r=NaN(101,1);
gamma_r=NaN(101,1);
sigma_sqr_r=NaN(101,1);
sharpe_r=NaN(101,1);

for i=1:k;
   mu_r(i,1)=ll+(ul-ll)/(k-1)*(i-1);
   gamma_r(i,1)=(mu_r(i,1)-R)*inv(C-2*R*B+R^2*A);
   weight_r=(gamma_r(i,1)*isigma*(z_bar-R*Vone));
   sigmar2_r(i,1)=((mu_r(i,1)-R)^2)*inv(C-2*R*B+R^2*A);
   sigma_sqr_r(i,1)=sqrt(sigmar2_r(i,1));
   sharpe_r(i,1)=(mu_r(i,1)-R)/sigma_sqr_r(i,1);
end;

mu_ra=R+sigma_sqr_r*sqrt(C-2*R*B+R^2*A);
mumin_r=R;
ind_r=(mu_r>=mumin_r); %Indicates efficient horizon
ind2_r=(mu_r<mumin_r); %Indicates locus below efficient horizon

% Creat plot
figure(2)

for i=1:n
    hold on
    plot(sqrt(sigma(i,i)),z_bar(i),'k^')
        text(sqrt(sigma(i,i)),z_bar(i),['  ' names{i}])
end

f3=plot(sigma_sqr_r(ind_r), mu_r(ind_r),'-',sigma_sqr_r(ind2_r),mu_r(ind2_r),':');
set(f3,'linewidth',2);
set(f3,'color','blue');
xlabel('Standard deviation of return','fontsize',14);
ylabel('Expected return','fontsize',14);
set(gca,'xlim',[0.00,0.2]);
set(gca,'ylim',[-0.01,0.03]);
grid;
title({'Mean-variance locus ( with the risk-free asset )',...
    '10 industries choosed: Food, Smoke, Clths, Steel, Oil, Agric, Drugs, Txtls, Autos, Chips'},...
    'fontsize',14);

%**********************************************

% 3) tangency portfolio and sharpe ratio

weight_0t=0;
weight_t=(isigma*(z_bar-R*Vone))/(B-A*R);  %weight of tangency portfolio
mu_t=(C-B*R)/(B-A*R);                      %mean of tangency portfolio
sigma_t2=(C-2*R*B+R^2*A)/((B-A*R)^2);      %variance of tangency portfolio
sharpe_t=(mu_t-R)/sqrt(sigma_t2);          %sharpe ratio of tangency portfolio

%sharpesort=NaN(101,1);
%[sharpesort,index]=sort(sharpe,'descend');

figure(3)
plot(sigma_sqr,sharpe,'--');
hold on;
scatter(sqrt(sigma_t2),sharpe_t,'r');
xlabel('Standard deviation of return','fontsize',14);
ylabel('Sharpe ratio','fontsize',14);
title('Tengency portfolio,sharpe ratio','fontsize',14);
hold off;

%**********************************************

% 4) "mean-variance locus" graph (without the risk-free asset) with the
% short sale constraints

% Minimum variance portfolio

f = zeros(n, 1);       % There is no constant
A_4 = -1*eye(n,n);     % Besides the returns we forbid short selling
b_4 = zeros(n,1);     % Required return and weights greater/eqauls 0
Aeq_4m = ones(1,n);
beq_4m = 1;
w_4m(:,1) = quadprog(sigma,zeros(n,1),A_4,b_4,Aeq_4m,beq_4m);
mu_4m(1) = z_bar'*w_4m(:,1);
sigma_4m(1) = sqrt(w_4m(:,1)'*sigma*w_4m(:,1));

% Mean-variance portfolio locus
ngrid=1000;
mu_4 = linspace(min(z_bar),max(z_bar),ngrid);
Aeq_4 =[ones(1,n);z_bar'];   % All weights should sum up
sigma_4=mu_4;
w_4=zeros(n,ngrid);

% solve the optimization
for i = 1:ngrid;
   w_4(:,i)= quadprog(sigma,f,A_4,b_4,Aeq_4,[1;mu_4(i)],zeros(n,1));
   sigma_4(i) = sqrt(w_4(:,i)'*sigma*w_4(:,i));
   sharpe_4(i)=(mu_4(i)-R)/sigma_4(i);
end ;

%scatter(sigma_4,mu_4,'.','b')

ind_4=(mu_4>mu_4m); %Indicates efficient horizon
ind2_4=(mu_4<mu_4m);

figure(4)

for i=1:n
    hold on
    plot(sqrt(sigma(i,i)),z_bar(i),'k^')
        text(sqrt(sigma(i,i)),z_bar(i),['  ' names{i}])
end

f4=plot(sigma_4(ind_4), mu_4(ind_4),':',sigma_4(ind2_4),mu_4(ind2_4),'-',...
    sigma_4m,mu_4m,'.');
%
set(f4(1:2),'linewidth',2);
set(f4(1:2),'color','blue');
set(f4(3),'markersize',40);
set(f4(3),'color','red');
%
xlabel('Standard deviation of return','fontsize',14);
ylabel('Expected return','fontsize',14);
grid;
title({'Mean-variance locus ( without short selling )',...
    '5 industries choosed: Food, Smoke, Clths, Steel, Oil'},...
    'fontsize',14);
hold off;

%**********************************************

% 5) "mean-variance locus" graph (with risk-free asset) with the
% short sale constraints

k=500;
A_5 = -1*eye(n,n);
b_5 = zeros(n,1);

mu_5=linspace(min(z_bar),max(z_bar),k);
for i = 1:k;
    Aeq_5 = z_bar'-R*Vone';
    beq_5 = mu_5(i)-R;
    [w_5(:,i)] = quadprog(sigma,zeros(1,n),A_5,b_5,Aeq_5,beq_5,[],[],[]);
    sigma_5(i) =sqrt(w_5(:,i)'*sigma*w_5(:,i));
    w0_5(:,i) = 1-Vone'*w_5(:,i);
    w0_5a(:,i) = abs(w0_5(:,i));
end

ind = find(sigma_5==min(sigma_5));
sigma_5m = sigma_5(1,ind);
mu_5m = mu_5(1,ind);
w_5m = w_5(n,ind);
ind_5=(mu_5>mu_5m); %Indicates efficient horizon
ind2_5=(mu_5<mu_5m);

figure(5)
f5=plot(sigma_5(ind_5), mu_5(ind_5),':',sigma_5(ind2_5),mu_5(ind2_5),'-',...
    sigma_5m(1),mu_5m(1),'.');
%
set(f5(1:2),'linewidth',2);
set(f5(1:2),'color','blue');
set(f5(3),'markersize',40);
set(f5(3),'color','red');
%
xlabel('Standard deviation of return','fontsize',14);
ylabel('Expected return','fontsize',14);
grid;
title({'Mean-variance locus ( with the risk-free asset ) with short-sale constriants',...
    '5 industries choosed: Food, Smoke, Clths, Steel, Oil'},...
    'fontsize',14);
hold on

for i=1:n
    hold on
    plot(sqrt(sigma(i,i)),z_bar(i),'k^')
        text(sqrt(sigma(i,i)),z_bar(i),['  ' names{i}])
end
hold off;


figure(6)
scatter(sigma_4,mu_4,'.','b')
hold on
scatter(sigma_5,mu_5,'.','r')
xlabel('Standard deviation of return','fontsize',14);
ylabel('Return','fontsize',14);
title({'Tangency portfolio'},'fontsize',14)
hold on

for i=1:n
    hold on
    plot(sqrt(sigma(i,i)),z_bar(i),'k^')
        text(sqrt(sigma(i,i)),z_bar(i),['  ' names{i}])
end
hold off;




%**********************************************

% 6) Tangent portfolio and sharpe ratio

mu_5_2=linspace(mu_5m,max(z_bar),k);
for i = 2:k;
    Aeq_5 = z_bar'-R*Vone';
    beq_5_2 = mu_5_2(i)-R;
    [w_5_2(:,i)] = quadprog(sigma,zeros(1,n),A_5,b_5,Aeq_5,beq_5_2,[],[],[]);
    sigma_5_2(i) =sqrt(w_5_2(:,i)'*sigma*w_5_2(:,i));
    w0_5_2(:,i) = 1-Vone'*w_5_2(:,i);
    w0_5a_2(:,i) = abs(w0_5_2(:,i));
end
mu_5_2(1)=mu_5m;
sigma_5_2(1)=sigma_5m;
w_5_2(:,1)=w_5m;
w0_5_2(:,1) = 1-Vone'*w_5_2(:,1);
w0_5a_2(:,1) = abs(w0_5_2(:,1));
min(w0_5a_2)
ind=find(w0_5a_2==min(w0_5a_2));                 
weight6_0t = w0_5a_2(1,ind);                      % w0 of tangency portfolio
weight6_t = w_5_2(:,ind);                         % weight of tangency portfolio
mu6_t = (z_bar'-R*Vone')*weight6_t+R;             % mean of tangency portfolio
sigma6_t = sqrt(weight6_t'*sigma*weight6_t);      % variance of tangency portfolio
sharpe6_t=(mu6_t-R)/sigma6_t;                     % sharpe ratio of tangency portfolio

max(sharpe_4)
sharpe6_t

figure(7)
plot(sigma_4,sharpe_4,'--');
hold on;
scatter(sigma6_t,sharpe6_t,'r');
xlabel('Standard deviation of return','fontsize',14);
ylabel('Sharpe ratio','fontsize',14);
title('Tengency portfolio,sharpe ratio','fontsize',14);
hold off;




%**********************************************

% 7) the uncertainty around the mean-variance efficient locus

% bootstrap for every asset
rng default  % For reproducibility
[bootstat_1,bootsam_1] = bootstrp(1000,@(x)[mean(x) std(x)],dataset(:,1)');
[bootstat_2,bootsam_2] = bootstrp(1000,@(x)[mean(x) std(x)],dataset(:,2)');
[bootstat_3,bootsam_3] = bootstrp(1000,@(x)[mean(x) std(x)],dataset(:,3)');
[bootstat_4,bootsam_4] = bootstrp(1000,@(x)[mean(x) std(x)],dataset(:,4)');
[bootstat_5,bootsam_5] = bootstrp(1000,@(x)[mean(x) std(x)],dataset(:,5)');
[bootstat_6,bootsam_6] = bootstrp(1000,@(x)[mean(x) std(x)],dataset(:,6)');
[bootstat_7,bootsam_7] = bootstrp(1000,@(x)[mean(x) std(x)],dataset(:,7)');
[bootstat_8,bootsam_8] = bootstrp(1000,@(x)[mean(x) std(x)],dataset(:,8)');
[bootstat_9,bootsam_9] = bootstrp(1000,@(x)[mean(x) std(x)],dataset(:,9)');
[bootstat_10,bootsam_10] = bootstrp(1000,@(x)[mean(x) std(x)],dataset(:,10)');



mean_std_1 = sort(bootstat_1);
mean_std_2 = sort(bootstat_2);
mean_std_3 = sort(bootstat_3);
mean_std_4 = sort(bootstat_4);
mean_std_5 = sort(bootstat_5);
mean_std_6 = sort(bootstat_6);
mean_std_7 = sort(bootstat_7);
mean_std_8 = sort(bootstat_8);
mean_std_9 = sort(bootstat_9);
mean_std_10 = sort(bootstat_10);

zbar_l(1,1) = mean_std_1(200,1);
zbar_l(2,1) = mean_std_2(200,1);
zbar_l(3,1) = mean_std_3(200,1);
zbar_l(4,1) = mean_std_4(200,1);
zbar_l(5,1) = mean_std_5(200,1);
zbar_l(6,1) = mean_std_6(200,1);
zbar_l(7,1) = mean_std_7(200,1);
zbar_l(8,1) = mean_std_8(200,1);
zbar_l(9,1) = mean_std_9(200,1);
zbar_l(10,1) = mean_std_10(200,1);
zbar_u(1,1) = mean_std_1(800,1);
zbar_u(2,1) = mean_std_2(800,1);
zbar_u(3,1) = mean_std_3(800,1);
zbar_u(4,1) = mean_std_4(800,1);
zbar_u(5,1) = mean_std_5(800,1);
zbar_u(6,1) = mean_std_6(800,1);
zbar_u(7,1) = mean_std_7(800,1);
zbar_u(8,1) = mean_std_8(800,1);
zbar_u(9,1) = mean_std_9(800,1);
zbar_u(10,1) = mean_std_10(800,1);

sigma_ll(1,1) = mean_std_1(400,2);
sigma_ll(1,2) = mean_std_2(400,2);
sigma_ll(1,3) = mean_std_3(400,2);
sigma_ll(1,4) = mean_std_4(400,2);
sigma_ll(1,5) = mean_std_5(400,2);
sigma_ll(1,6) = mean_std_6(400,2);
sigma_ll(1,7) = mean_std_7(400,2);
sigma_ll(1,8) = mean_std_8(400,2);
sigma_ll(1,9) = mean_std_9(400,2);
sigma_ll(1,10) = mean_std_10(400,2);
sigma_l = sigma_ll'*sigma_ll;
sigma_uu(1,1) = mean_std_1(600,2);
sigma_uu(1,2) = mean_std_2(600,2);
sigma_uu(1,3) = mean_std_3(600,2);
sigma_uu(1,4) = mean_std_4(600,2);
sigma_uu(1,5) = mean_std_5(600,2);
sigma_uu(1,6) = mean_std_6(600,2);
sigma_uu(1,7) = mean_std_7(600,2);
sigma_uu(1,8) = mean_std_8(600,2);
sigma_uu(1,9) = mean_std_9(600,2);
sigma_uu(1,10) = mean_std_10(600,2);
sigma_u = sigma_uu'*sigma_uu;



%*************************************

% 7___1) "mean-variance locus" graph (without the risk-free asset)

%efficient frontier for 'zbar, sigma'
% Minimum variance portfolio
zbar = (mean(dataset))';
sigma = cov(dataset);
R = 0.02/12;
n = length(zbar);
Vone=ones(n,1);        %vector of ones
f = zeros(n, 1);       % There is no constant
Aeq_1m = ones(1,n);
beq_1m = 1;
w_1m(:,1) = quadprog(sigma,f,[],[],Aeq_1m,beq_1m);
mu_1m(1) = zbar'*w_1m(:,1);
sigma_1m(1) = sqrt(w_1m(:,1)'*sigma*w_1m(:,1));

% Mean-variance portfolio locus
ngrid=500;
mu_1 = linspace(mu_1m(1), max(zbar), ngrid);
Aeq_1 =[ones(1,10);zbar'];   % All weights should sum up
sigma_1=mu_1;
w_1=zeros(10,ngrid);


% solve the optimization
w_1 = NaN(10,500);
sigma_1 = NaN(1,500);
sharpe_1 = NaN(1,500);
for i = 1:ngrid;
   w_1(:,i)= quadprog(sigma,f,[],[],Aeq_1,[1;mu_1(i)]);
   sigma_1(i) = sqrt(w_1(:,i)'*sigma*w_1(:,i));
   sharpe_1(i)=(mu_1(i)-R)/sigma_1(i);
end ;


%efficient frontier for 'zbar_l, sigma'
% Minimum variance portfolio
w_1mzl(:,1) = quadprog(sigma,zeros(1,n),[],[],Aeq_1m,beq_1m);
mu_1mzl(1) = zbar_l'*w_1mzl(:,1);
sigma_1mzl(1) = sqrt(w_1mzl(:,1)'*sigma*w_1mzl(:,1));

% Mean-variance portfolio locus
mu_1zl = linspace(mu_1mzl(1), max(zbar_l), ngrid);
Aeq_1zl =[ones(1,10);zbar_l'];   % All weights should sum up
sigma_1zl=mu_1zl;
w_1zl=zeros(10,ngrid);


% solve the optimization
w_1zl = NaN(10,500);
sigma_1zl = NaN(1,500);
sharpe_1zl = NaN(1,500);
for i = 1:ngrid;
   w_1zl(:,i)= quadprog(sigma,f,[],[],Aeq_1zl,[1;mu_1zl(i)]);
   sigma_1zl(i) = sqrt(w_1zl(:,i)'*sigma*w_1zl(:,i));
   sharpe_1zl(i)=(mu_1zl(i)-R)/sigma_1zl(i);
end ;


%efficient frontier for 'zbar_u, sigma'
% Minimum variance portfolio
w_1mzu(:,1) = quadprog(sigma,zeros(1,n),[],[],Aeq_1m,beq_1m);
mu_1mzu(1) = zbar_u'*w_1mzu(:,1);
sigma_1mzu(1) = sqrt(w_1mzu(:,1)'*sigma*w_1mzu(:,1));

% Mean-variance portfolio locus
mu_1zu = linspace(mu_1mzu(1), max(zbar_u), ngrid);
Aeq_1zu =[ones(1,10);zbar_u'];   % All weights should sum up
sigma_1zu=mu_1zu;
w_1zu=zeros(10,ngrid);


% solve the optimization
w_1zu = NaN(10,500);
sigma_1zu = NaN(1,500);
sharpe_1zu = NaN(1,500);
for i = 1:ngrid;
   w_1zu(:,i)= quadprog(sigma,f,[],[],Aeq_1zu,[1;mu_1zu(i)]);
   sigma_1zu(i) = sqrt(w_1zu(:,i)'*sigma*w_1zu(:,i));
   sharpe_1zu(i)=(mu_1zu(i)-R)/sigma_1zu(i);
end ;


%create plot

figure(8)

f1 = plot(sigma_1, mu_1,'r.-');
hold on;
f2 = plot(sigma_1zl, mu_1zl,'b.-');
hold on;
f3 = plot(sigma_1zu, mu_1zu,'k.-');
hold on;
xlabel('Standard deviation of return','fontsize',20);
ylabel('Expected return','fontsize',20);
legend('mu','mu_l','mu_u');
title('sensitivity of efficient frontier to zbar');
grid;


%efficient frontier for 'zbar, sigma_l'
% Minimum variance portfolio
w_1msl(:,1) = quadprog(sigma_l,zeros(1,n),[],[],Aeq_1m,beq_1m);
mu_1msl(1) = zbar'*w_1msl(:,1);
sigma_1msl(1) = sqrt(w_1msl(:,1)'*sigma_l*w_1msl(:,1));

% Mean-variance portfolio locus
mu_1sl = linspace(mu_1msl(1), max(zbar), ngrid);
Aeq_1sl =[ones(1,10);zbar'];   % All weights should sum up
sigma_1sl=mu_1sl;
w_1sl=zeros(10,ngrid);


% solve the optimization
w_1sl = NaN(10,500);
sigma_1sl = NaN(1,500);
sharpe_1sl = NaN(1,500);
for i = 1:ngrid;
   w_1sl(:,i)= quadprog(sigma_l,f,[],[],Aeq_1sl,[1;mu_1sl(i)],zeros(10,1));
   sigma_1sl(i) = sqrt(w_1sl(:,i)'*sigma_l*w_1sl(:,i));
   sharpe_1sl(i)=(mu_1sl(i)-R)/sigma_1sl(i);
end ;


%efficient frontier for 'zbar, sigma_u'
% Minimum variance portfolio
w_1msu(:,1) = quadprog(sigma_u,f,[],[],Aeq_1m,beq_1m);
mu_1msu(1) = zbar'*w_1msu(:,1);
sigma_1msu(1) = sqrt(w_1msu(:,1)'*sigma_u*w_1msu(:,1));

% Mean-variance portfolio locus
mu_1su = linspace(mu_1msu(1), max(zbar), ngrid);
Aeq_1su =[ones(1,10);zbar'];   % All weights should sum up
sigma_1su=mu_1su;
w_1su=zeros(10,ngrid);

% solve the optimization
w_1su = NaN(10,500);
sigma_1su = NaN(1,500);
sharpe_1su = NaN(1,500);
for i = 1:ngrid;
   w_1su(:,i)= quadprog(sigma_u,f,[],[],Aeq_1su,[1;mu_1su(i)]);
   sigma_1su(i) = sqrt(w_1su(:,i)'*sigma_u*w_1su(:,i));
   sharpe_1su(i)=(mu_1su(i)-R)/sigma_1su(i);
end ;


%create plot
figure(9)
f1 = plot(sigma_1, mu_1,'r.-');
hold on;
f2 = plot(sigma_1sl, mu_1sl,'b.-');
hold on;
xlabel('Standard deviation of return','fontsize',20);
ylabel('Expected return','fontsize',20);
set(gca,'xlim',[0.03,0.04]);
legend('var-cov','var-cov_l');
title('sensitivity of efficient frontier to var-cov matrix');
grid;
hold off;

figure(10)
f3 = plot(sigma_1su, mu_1su,'k.-');
xlabel('Standard deviation of return','fontsize',20);
ylabel('Expected return','fontsize',20);
legend('var-cov_u');
title('sensitivity of efficient frontier to var-cov matrix');
grid;



%*************************************

% 7___2) "mean-variance locus" graph (with the risk-free asset)


% efficient frontier for 'zbar, sigma'
w_2 = NaN(10,500);
sigma_2 = NaN(1,500);
sharpe_2 = NaN(1,500);
w0_2 = NaN(1,500);
w0_2a = NaN(1,500);
mu_2 = linspace(R,max(zbar),ngrid);
for i = 1:ngrid;
    Aeq_2 = zbar'-R*Vone';
    beq_2 = mu_2(i)-R;
    [w_2(:,i)] = quadprog(sigma,f,[],[],Aeq_2,beq_2,[],[],[]);
    sigma_2(i) =sqrt(w_2(:,i)'*sigma*w_2(:,i));
    sharpe_2(i)=(mu_2(i)-R)/sigma_2(i);
    w0_2(:,i) = 1-Vone'*w_2(:,i);
    w0_2a(:,i) = abs(w0_2(:,i));
end


% efficient frontier for 'zbar_l, sigma'
w_2zl = NaN(10,500);
sigma_2zl = NaN(1,500);
sharpe_2zl = NaN(1,500);
w0_2zl = NaN(1,500);
w0_2azl = NaN(1,500);
mu_2zl = linspace(R,max(zbar_l),ngrid);
for i = 1:ngrid;
    Aeq_2zl = zbar_l'-R*Vone';
    beq_2zl = mu_2zl(i)-R;
    [w_2zl(:,i)] = quadprog(sigma,f,[],[],Aeq_2zl,beq_2zl,[],[],[]);
    sigma_2zl(i) =sqrt(w_2zl(:,i)'*sigma*w_2zl(:,i));
    sharpe_2zl(i)=(mu_2zl(i)-R)/sigma_2zl(i);
    w0_2zl(:,i) = 1-Vone'*w_2zl(:,i);
    w0_2azl(:,i) = abs(w0_2zl(:,i));
end



% efficient frontier for 'zbar_u, sigma'
w_2zu = NaN(10,500);
sigma_2zu = NaN(1,500);
sharpe_2zu = NaN(1,500);
w0_2zu = NaN(1,500);
w0_2azu = NaN(1,500);
mu_2zu = linspace(R,max(zbar_u),ngrid);
for i = 1:ngrid;
    Aeq_2zu = zbar_u'-R*Vone';
    beq_2zu = mu_2zu(i)-R;
    [w_2zu(:,i)] = quadprog(sigma,f,[],[],Aeq_2zu,beq_2zu,[],[],[]);
    sigma_2zu(i) =sqrt(w_2zu(:,i)'*sigma*w_2zu(:,i));
    sharpe_2zu(i)=(mu_2zu(i)-R)/sigma_2zu(i);
    w0_2zu(:,i) = 1-Vone'*w_2zu(:,i);
    w0_2azu(:,i) = abs(w0_2zu(:,i));
end

%create plot

figure(11)

f1 = plot(sigma_2, mu_2,'r.-');
hold on;
f2 = plot(sigma_2zl, mu_2zl,'b.-');
hold on;
f3 = plot(sigma_2zu, mu_2zu,'k.-');
hold on;
xlabel('Standard deviation of return','fontsize',20);
ylabel('Expected return','fontsize',20);
legend('mu','mu_l','mu_u');
title('sensitivity of efficient frontier to zbar');
grid;


% efficient frontier for 'zbar, sigma_l'
w_2sl = NaN(10,500);
sigma_2sl = NaN(1,500);
sharpe_2sl = NaN(1,500);
w0_2sl = NaN(1,500);
w0_2asl = NaN(1,500);
mu_2sl = linspace(R,max(zbar),ngrid);
for i = 1:ngrid;
    Aeq_2sl = zbar'-R*Vone';
    beq_2sl = mu_2sl(i)-R;
    [w_2sl(:,i)] = quadprog(sigma_l,f,[],[],Aeq_2sl,beq_2sl,[],[],[]);
    sigma_2sl(i) =sqrt(w_2sl(:,i)'*sigma_l*w_2sl(:,i));
    sharpe_2sl(i)=(mu_2sl(i)-R)/sigma_2sl(i);
    w0_2sl(:,i) = 1-Vone'*w_2sl(:,i);
    w0_2asl(:,i) = abs(w0_2sl(:,i));
end


% efficient frontier for 'zbar, sigma_u'
w_2su = NaN(10,500);
sigma_2su = NaN(1,500);
sharpe_2su = NaN(1,500);
w0_2su = NaN(1,500);
w0_2asu = NaN(1,500);
mu_2su = linspace(R,max(zbar),ngrid);
for i = 1:ngrid;
    Aeq_2su = zbar'-R*Vone';
    beq_2su = mu_2su(i)-R;
    [w_2su(:,i)] = quadprog(sigma_u,f,[],[],Aeq_2su,beq_2su,[],[],[]);
    sigma_2su(i) =sqrt(w_2su(:,i)'*sigma_u*w_2su(:,i));
    sharpe_2su(i)=(mu_2su(i)-R)/sigma_2su(i);
    w0_2su(:,i) = 1-Vone'*w_2su(:,i);
    w0_2asu(:,i) = abs(w0_2su(:,i));
end


% Creat plot
figure(12)
f1 = plot(sigma_2, mu_2,'r.-');
xlabel('Standard deviation of return','fontsize',20);
ylabel('Expected return','fontsize',20);
legend('var-cov');
title('sensitivity of efficient frontier to var-cov matrix');
grid;

figure(13)
f2 = plot(sigma_2sl, mu_2sl,'b.-');
xlabel('Standard deviation of return','fontsize',20);
ylabel('Expected return','fontsize',20);
legend('var-cov_l');
title('sensitivity of efficient frontier to var-cov matrix');
grid;

figure(14)
f3 = plot(sigma_2su, mu_2su,'k.-');
xlabel('Standard deviation of return','fontsize',20);
ylabel('Expected return','fontsize',20);
legend('var-cov_u');
title('sensitivity of efficient frontier to var-cov matrix');
grid;


%*****************************************

% 7___3) tangency portfolio and sharpe ratio


%tangency portfolio for 'zbar, sigma'
ind=find(w0_2a==min(w0_2a));                 
weight2_0t = w0_2a(1,ind);                      % w0 of tangency portfolio
weight2_t = w_2(:,ind);                         % weight of tangency portfolio
mu2_t = (zbar'-R*Vone')*weight2_t+R;             % mean of tangency portfolio
sigma2_t = sqrt(weight2_t'*sigma*weight2_t);      % variance of tangency portfolio
sharpe2_t=(mu2_t-R)/sigma2_t;                     % sharpe ratio of tangency portfolio


%tangency portfolio for 'zbar_l, sigma'
ind=find(w0_2azl==min(w0_2azl));                 
weight2zl_0t = w0_2azl(1,ind);                      % w0 of tangency portfolio
weight2zl_t = w_2zl(:,ind);                         % weight of tangency portfolio
mu2zl_t = (zbar_l'-R*Vone')*weight2zl_t+R;             % mean of tangency portfolio
sigma2zl_t = sqrt(weight2zl_t'*sigma*weight2zl_t);      % variance of tangency portfolio
sharpe2zl_t=(mu2zl_t-R)/sigma2zl_t;                     % sharpe ratio of tangency portfolio


%tangency portfolio for 'zbar_u, sigma'
ind=find(w0_2azu==min(w0_2azu));                 
weight2zu_0t = w0_2azu(1,ind);                      % w0 of tangency portfolio
weight2zu_t = w_2zu(:,ind);                         % weight of tangency portfolio
mu2zu_t = (zbar_u'-R*Vone')*weight2zu_t+R;             % mean of tangency portfolio
sigma2zu_t = sqrt(weight2zu_t'*sigma*weight2zu_t);      % variance of tangency portfolio
sharpe2zu_t=(mu2zu_t-R)/sigma2zu_t;                     % sharpe ratio of tangency portfolio


% create plot
figure(15)
plot(sigma_1zl,sharpe_1zl,'.');
hold on;
scatter(sigma2zl_t,sharpe2zl_t,'r');
xlabel('Standard deviation of return','fontsize',20);
ylabel('Sharpe ratio','fontsize',20);
hold off;

figure(16)
plot(sigma_1zu,sharpe_1zu,'--');
hold on;
scatter(sigma2zu_t,sharpe2zu_t,'r');
xlabel('Standard deviation of return','fontsize',20);
ylabel('Sharpe ratio','fontsize',20);
legend('mu_u');
title('tangency portfolio of 80% zbar');
hold off;


%tangency portfolio for 'zbar, sigma_l'

ind=find(w0_2asl==min(w0_2asl));                 
weight2sl_0t = w0_2asl(1,ind);                      % w0 of tangency portfolio
weight2sl_t = w_2sl(:,ind);                         % weight of tangency portfolio
mu2sl_t = (zbar'-R*Vone')*weight2sl_t+R;             % mean of tangency portfolio
sigma2sl_t = sqrt(weight2sl_t'*sigma_l*weight2sl_t);      % variance of tangency portfolio
sharpe2sl_t = (mu2sl_t-R)/sigma2sl_t;                     % sharpe ratio of tangency portfolio


%tangency portfolio for 'zbar, sigma_u'
ind=find(w0_2asu==min(w0_2asu));                 
weight2su_0t = w0_2asu(1,ind);                      % w0 of tangency portfolio
weight2su_t = w_2su(:,ind);                         % weight of tangency portfolio
mu2su_t = (zbar'-R*Vone')*weight2su_t+R;             % mean of tangency portfolio
sigma2su_t = sqrt(weight2su_t'*sigma_u*weight2su_t);      % variance of tangency portfolio
sharpe2su_t=(mu2su_t-R)/sigma2su_t;                     % sharpe ratio of tangency portfolio

%create plot
figure(17)
plot(sigma_1sl,sharpe_1sl,'.');
hold on;
scatter(sigma2sl_t,sharpe2sl_t,'r');
xlabel('Standard deviation of return','fontsize',20);
ylabel('Sharpe ratio','fontsize',20);
hold off;
figure(18)
plot(sigma_1su,sharpe_1su,'--');
hold on;
scatter(sigma2su_t,sharpe2su_t,'r');
xlabel('Standard deviation of return','fontsize',20);
ylabel('Sharpe ratio','fontsize',20);
hold off;


%*****************************************

% 7___4) "mean-variance locus" graph (without the risk-free asset) with the
% short sale constraints


%efficient frontier for 'zbar, sigma'

% Minimum variance portfolio
A_4 = -1*eye(n,n);     % Besides the returns we forbid short selling
b_4 = zeros(10, 1);     % Required return and weights greater/eqauls 0
Aeq_4m = ones(1,n);
beq_4m = 1;
w_4m(:,1) = quadprog(sigma,zeros(1,n),A_4,b_4,Aeq_4m,beq_4m);
mu_4m(1) = zbar'*w_4m(:,1);
sigma_4m(1) = sqrt(w_4m(:,1)'*sigma*w_4m(:,1));

% Mean-variance portfolio locus
mu_4 = linspace(mu_4m(1), max(zbar), ngrid);
Aeq_4 =[ones(1,10);zbar'];   % All weights should sum up
sigma_4=mu_4;
w_4=zeros(10,ngrid);

% solve the optimization
w_4 = NaN(10,500);
sigma_4 = NaN(1,500);
sharpe_4 = NaN(1,500);
for i = 1:ngrid;
   w_4(:,i)= quadprog(sigma,f,A_4,b_4,Aeq_4,[1;mu_4(i)],zeros(10,1));
   sigma_4(i) = sqrt(w_4(:,i)'*sigma*w_4(:,i));
   sharpe_4(i)=(mu_4(i)-R)/sigma_4(i);
end ;


%efficient frontier for 'zbar_l, sigma'

% Minimum variance portfolio
w_4mzl(:,1) = quadprog(sigma,zeros(1,n),A_4,b_4,Aeq_4m,beq_4m);
mu_4mzl(1) = zbar_l'*w_4mzl(:,1);
sigma_4mzl(1) = sqrt(w_4mzl(:,1)'*sigma*w_4mzl(:,1));

% Mean-variance portfolio locus
mu_4zl = linspace(mu_4mzl(1), max(zbar_l), ngrid);
Aeq_4zl =[ones(1,10);zbar_l'];   % All weights should sum up
sigma_4zl=mu_4zl;
w_4zl=zeros(10,ngrid);

% solve the optimization
w_4zl = NaN(10,500);
sigma_4zl = NaN(1,500);
sharpe_4zl = NaN(1,500);
for i = 1:ngrid;
   w_4zl(:,i)= quadprog(sigma,f,A_4,b_4,Aeq_4zl,[1;mu_4zl(i)],zeros(10,1));
   sigma_4zl(i) = sqrt(w_4zl(:,i)'*sigma*w_4zl(:,i));
   sharpe_4zl(i)=(mu_4zl(i)-R)/sigma_4zl(i);
end ;


%efficient frontier for 'zbar_u, sigma'

% Minimum variance portfolio
w_4mzu(:,1) = quadprog(sigma,zeros(1,n),A_4,b_4,Aeq_4m,beq_4m);
mu_4mzu(1) = zbar_u'*w_4mzu(:,1);
sigma_4mzu(1) = sqrt(w_4mzu(:,1)'*sigma*w_4mzu(:,1));

% Mean-variance portfolio locus
mu_4zu = linspace(mu_4mzu(1), max(zbar_u), ngrid);
Aeq_4zu =[ones(1,10);zbar_u'];   % All weights should sum up
sigma_4zu=mu_4zu;
w_4zu=zeros(10,ngrid);

% solve the optimization
w_4zu = NaN(10,500);
sigma_4zu = NaN(1,500);
sharpe_4zu = NaN(1,500);
for i = 1:ngrid;
   w_4zu(:,i)= quadprog(sigma,f,A_4,b_4,Aeq_4zu,[1;mu_4zu(i)],zeros(10,1));
   sigma_4zu(i) = sqrt(w_4zu(:,i)'*sigma*w_4zu(:,i));
   sharpe_4zu(i)=(mu_4zu(i)-R)/sigma_4zu(i);
end ;


%create plot

figure(19)
f1 = plot(sigma_4, mu_4,'r.-');
hold on;
f2 = plot(sigma_4zl, mu_4zl,'b.-');
hold on;
f3 = plot(sigma_4zu, mu_4zu,'k.-');
hold on;
xlabel('Standard deviation of return','fontsize',20);
ylabel('Expected return','fontsize',20);
set(gca,'xlim',[0.03,0.05]);
set(gca,'ylim',[-2*10^(-3),15*10^(-3)]);
legend('mu','mu_l','mu_u');
title('sensitivity of efficient frontier(with short-sale constraints) to zbar');
grid;
hold off;


%efficient frontier for 'zbar, sigma_l'

% Minimum variance portfolio
w_4msl(:,1) = quadprog(sigma_l,zeros(1,n),A_4,b_4,Aeq_4m,beq_4m);
mu_4msl(1) = zbar'*w_4msl(:,1);
sigma_4msl(1) = sqrt(w_4msl(:,1)'*sigma_l*w_4msl(:,1));

% Mean-variance portfolio locus
mu_4sl = linspace(mu_4msl(1), max(zbar), ngrid);
Aeq_4sl =[ones(1,10);zbar'];   % All weights should sum up
sigma_4sl=mu_4sl;
w_4sl=zeros(10,ngrid);

% solve the optimization
w_4sl = NaN(10,500);
sigma_4sl = NaN(1,500);
sharpe_4sl = NaN(1,500);
for i = 1:ngrid;
   w_4sl(:,i)= quadprog(sigma_l,f,A_4,b_4,Aeq_4sl,[1;mu_4sl(i)],zeros(10,1));
   sigma_4sl(i) = sqrt(w_4sl(:,i)'*sigma_l*w_4sl(:,i));
   sharpe_4sl(i)=(mu_4sl(i)-R)/sigma_4sl(i);
end ;



%efficient frontier for 'zbar, sigma_u'

% Minimum variance portfolio
w_4msu(:,1) = quadprog(sigma_u,zeros(1,n),A_4,b_4,Aeq_4m,beq_4m);
mu_4msu(1) = zbar'*w_4msu(:,1);
sigma_4msu(1) = sqrt(w_4msu(:,1)'*sigma_u*w_4msu(:,1));

% Mean-variance portfolio locus
mu_4su = linspace(mu_4msu(1), max(zbar), ngrid);
Aeq_4su =[ones(1,10);zbar'];   % All weights should sum up
sigma_4su=mu_4su;
w_4su=zeros(10,ngrid);

% solve the optimization
w_4su = NaN(10,500);
sigma_4su = NaN(1,500);
sharpe_4su = NaN(1,500);
for i = 1:ngrid;
   w_4su(:,i)= quadprog(sigma_u,f,A_4,b_4,Aeq_4su,[1;mu_4su(i)],zeros(10,1));
   sigma_4su(i) = sqrt(w_4su(:,i)'*sigma_u*w_4su(:,i));
   sharpe_4su(i)=(mu_4su(i)-R)/sigma_4su(i);
end ;


%create plot

figure(20)
f1 = plot(sigma_4, mu_4,'r.-');
hold on;
f2 = plot(sigma_4sl, mu_4sl,'b.-');
hold on;
f3 = plot(sigma_4su, mu_4su,'k.-');
hold on;
xlabel('Standard deviation of return','fontsize',20);
ylabel('Expected return','fontsize',20);
legend('var-cov','var-cov_l','var-cov_u');
title('sensitivity of efficient frontier(with short-sale constraints) to variance-covariance matrix');
grid;
hold off;


%*****************************************

% 7___5) "mean-variance locus" graph (with risk-free asset) with the
% short sale constraints

%efficient frontier for 'zbar, sigma'
A_5 = -1*eye(n,n);
b_5 = zeros(n,1);
w_5 = NaN(10,500);
sigma_5 = NaN(1,500);
sharpe_5 = NaN(1,500);
w0_5 = NaN(1,500);
w0_5a = NaN(1,500);
mu_5 = linspace(min(zbar),max(zbar),ngrid);
for i = 1:ngrid;
    Aeq_5 = zbar'-R*Vone';
    beq_5 = mu_5(i)-R;
    [w_5(:,i)] = quadprog(sigma,zeros(1,n),A_5,b_5,Aeq_5,beq_5,[],[],[]);
    sigma_5(i) =sqrt(w_5(:,i)'*sigma*w_5(:,i));
    sharpe_5(i)=(mu_5(i)-R)/sigma_5(i);
    w0_5(:,i) = 1-Vone'*w_5(:,i);
    w0_5a(:,i) = abs(w0_5(:,i));
end


%efficient frontier for 'zbar_l, sigma'
w_5zl = NaN(10,500);
sigma_5zl = NaN(1,500);
sharpe_5zl = NaN(1,500);
w0_5zl = NaN(1,500);
w0_5azl = NaN(1,500);
mu_5zl = linspace(min(zbar_l),max(zbar_l),ngrid);
for i = 1:ngrid;
    Aeq_5zl = zbar_l'-R*Vone';
    beq_5zl = mu_5zl(i)-R;
    [w_5zl(:,i)] = quadprog(sigma,f,A_5,b_5,Aeq_5zl,beq_5zl,[],[],[]);
    sigma_5zl(i) =sqrt(w_5zl(:,i)'*sigma*w_5zl(:,i));
    sharpe_5zl(i)=(mu_5zl(i)-R)/sigma_5zl(i);
    w0_5zl(:,i) = 1-Vone'*w_5zl(:,i);
    w0_5azl(:,i) = abs(w0_5zl(:,i));
end



%efficient frontier for 'zbar_u, sigma'
w_5zu = NaN(10,500);
sigma_5zu = NaN(1,500);
sharpe_5zu = NaN(1,500);
w0_5zu = NaN(1,500);
w0_5azu = NaN(1,500);
mu_5zu = linspace(min(zbar_u),max(zbar_u),ngrid);
for i = 1:ngrid;
    Aeq_5zu = zbar_u'-R*Vone';
    beq_5zu = mu_5zu(i)-R;
    [w_5zu(:,i)] = quadprog(sigma,f,A_5,b_5,Aeq_5zu,beq_5zu,[],[],[]);
    sigma_5zu(i) =sqrt(w_5zu(:,i)'*sigma*w_5zu(:,i));
    sharpe_5zu(i)=(mu_5zu(i)-R)/sigma_5zu(i);
    w0_5zu(:,i) = 1-Vone'*w_5zu(:,i);
    w0_5azu(:,i) = abs(w0_5zu(:,i));
end


%create plot
figure(21)

f1 = plot(sigma_5, mu_5,'r.-');
hold on;
f2 = plot(sigma_5zl, mu_5zl,'b.-');
hold on;
f3 = plot(sigma_5zu, mu_5zu,'k.-');
hold on;
xlabel('Standard deviation of return','fontsize',20);
ylabel('Expected return','fontsize',20);
legend('mu','mu_l','mu_u');
title('sensitivity of efficient frontier(with short-sale constraints) to zbar');
grid;
hold off;


%efficient frontier for 'zbar, sigma_l'
w_5sl = NaN(10,500);
sigma_5sl = NaN(1,500);
sharpe_5sl = NaN(1,500);
w0_5sl = NaN(1,500);
w0_5asl = NaN(1,500);
mu_5sl = linspace(min(zbar),max(zbar),ngrid);
for i = 1:ngrid;
    Aeq_5sl = zbar'-R*Vone';
    beq_5sl = mu_5sl(i)-R;
    [w_5sl(:,i)] = quadprog(sigma_l,f,A_5,b_5,Aeq_5sl,beq_5sl,[],[],[]);
    sigma_5sl(i) =sqrt(w_5sl(:,i)'*sigma_l*w_5sl(:,i));
    sharpe_5sl(i)=(mu_5sl(i)-R)/sigma_5sl(i);
    w0_5sl(:,i) = 1-Vone'*w_5sl(:,i);
    w0_5asl(:,i) = abs(w0_5sl(:,i));
end



%efficient frontier for 'zbar, sigma_u'
w_5su = NaN(10,500);
sigma_5su = NaN(1,500);
sharpe_5su = NaN(1,500);
w0_5su = NaN(1,500);
w0_5asu = NaN(1,500);
mu_5su = linspace(min(zbar),max(zbar),ngrid);
for i = 1:ngrid;
    Aeq_5su = zbar'-R*Vone';
    beq_5su = mu_5su(i)-R;
    [w_5su(:,i)] = quadprog(sigma_u,f,A_5,b_5,Aeq_5su,beq_5su,[],[],[]);
    sigma_5su(i) =sqrt(w_5su(:,i)'*sigma_u*w_5su(:,i));
    sharpe_5su(i)=(mu_5su(i)-R)/sigma_5su(i);
    w0_5su(:,i) = 1-Vone'*w_5su(:,i);
    w0_5asu(:,i) = abs(w0_5su(:,i));
end


%create plot

figure(22)

f1 = plot(sigma_5, mu_5,'r.-');
hold on;
f2 = plot(sigma_5sl, mu_5sl,'b.-');
hold on;
f3 = plot(sigma_5su, mu_5su,'k.-');
hold on;
xlabel('Standard deviation of return','fontsize',20);
ylabel('Expected return','fontsize',20);
legend('var-cov','var-cov_l','var-cov_u');
title('sensitivity of efficient frontier(with short-sale constraints) to variance-covariance matrix');
grid;
hold off;



%*****************************************

% 7___6) Tangent portfolio and sharpe ratio


%tangency portfolio for 'zbar, sigma'
ind=find(w0_5a==min(w0_5a));                 
weight5_0t = w0_5a(1,ind);                      % w0 of tangency portfolio
weight5_t = w_5(:,ind);                         % weight of tangency portfolio
mu5_t = (zbar'-R*Vone')*weight5_t+R;             % mean of tangency portfolio
sigma5_t = sqrt(weight5_t'*sigma*weight5_t);      % variance of tangency portfolio
sharpe5_t=(mu5_t-R)/sigma5_t;                     % sharpe ratio of tangency portfolio


%tangency portfolio for 'zbar_l, sigma'
ind=find(w0_5azl==min(w0_5azl));                 
weight5zl_0t = w0_5azl(1,ind);                      % w0 of tangency portfolio
weight5zl_t = w_5zl(:,ind);                         % weight of tangency portfolio
mu5zl_t = (zbar_l'-R*Vone')*weight5zl_t+R;             % mean of tangency portfolio
sigma5zl_t = sqrt(weight5zl_t'*sigma*weight5zl_t);      % variance of tangency portfolio
sharpe5zl_t=(mu5zl_t-R)/sigma5zl_t;                     % sharpe ratio of tangency portfolio


%tangency portfolio for 'zbar_u, sigma'
ind=find(w0_5azu==min(w0_5azu));                 
weight5zu_0t = w0_5azu(1,ind);                      % w0 of tangency portfolio
weight5zu_t = w_5zu(:,ind);                         % weight of tangency portfolio
mu5zu_t = (zbar_u'-R*Vone')*weight5zu_t+R;             % mean of tangency portfolio
sigma5zu_t = sqrt(weight5zu_t'*sigma*weight5zu_t);      % variance of tangency portfolio
sharpe5zu_t=(mu5zu_t-R)/sigma5zu_t;                     % sharpe ratio of tangency portfolio


% create plot
figure(23)
plot(sigma_4zl,sharpe_4zl,'b.-');
hold on;
scatter(sigma5zl_t,sharpe5zl_t,'r');
xlabel('Standard deviation of return','fontsize',20);
ylabel('Sharpe ratio','fontsize',20);
legend('mu_l','highest sharperatio_l');
title('tangency portfolio of 20% zbar');
hold off;
figure(24)
plot(sigma_4zu,sharpe_4zu,'k.-');
hold on;
scatter(sigma5zu_t,sharpe5zu_t,'r');
xlabel('Standard deviation of return','fontsize',20);
ylabel('Sharpe ratio','fontsize',20);
legend('mu_u','highest sharperatio_u');
title('tangency portfolio of 80% zbar');
hold off;



%tangency portfolio for 'zbar, sigma_l'
ind=find(w0_5asl==min(w0_5asl));                 
weight5sl_0t = w0_5asl(1,ind);                      % w0 of tangency portfolio
weight5sl_t = w_5sl(:,ind);                         % weight of tangency portfolio
mu5sl_t = (zbar'-R*Vone')*weight5sl_t+R;             % mean of tangency portfolio
sigma5sl_t = sqrt(weight5sl_t'*sigma_l*weight5sl_t);      % variance of tangency portfolio
sharpe5sl_t=(mu5sl_t-R)/sigma5sl_t;                     % sharpe ratio of tangency portfolio


%tangency portfolio for 'zbar, sigma_u'
ind=find(w0_5asu==min(w0_5asu));                 
weight5su_0t = w0_5asu(1,ind);                      % w0 of tangency portfolio
weight5su_t = w_5su(:,ind);                         % weight of tangency portfolio
mu5su_t = (zbar'-R*Vone')*weight5su_t+R;             % mean of tangency portfolio
sigma5su_t = sqrt(weight5su_t'*sigma_u*weight5su_t);      % variance of tangency portfolio
sharpe5su_t=(mu5su_t-R)/sigma5su_t;                     % sharpe ratio of tangency portfolio

%create plot
figure(25)
plot(sigma_4sl,sharpe_4sl,'b.-');
hold on;
scatter(sigma5sl_t,sharpe5sl_t,'r');
hold on;
plot(sigma_4su,sharpe_4su,'k.-');
hold on;
scatter(sigma5su_t,sharpe5su_t,'g');
xlabel('Standard deviation of return','fontsize',20);
ylabel('Sharpe ratio','fontsize',20);
legend('var-cov_l','highest sharperatio_l','var-cov_u','highest sharperatio_u');
title('tangency portfolio of 40% and 60% variance-covariance matrix');
hold off;



%**********************************************

% 8) Invest in up to 10 industries (maximum) in the universe of 48 
%    industries which maximum sharpe ratio.

data8= xlsread('Q8_48_Industry_Portfolios.xlsm');
dateline8 = data8(:,1);
dataset8=data8(1:end,2:end)./100;

%% Benchemark with 48 industrie portfoilo

% 48 industrie portfolio without short-sell constraint

R=0.02/12;
sigma8=cov(dataset8);
z_bar8=(mean(dataset8))';
n=length(z_bar8);

f = zeros(n, 1);       % There is no constant
Aeq_8m = ones(1,n);
beq_8m = 1;
w_8m(:,1) = quadprog(sigma8,zeros(1,n),[],[],Aeq_8m,beq_8m);
mu_8m(1) = z_bar8'*w_8m(:,1);
sigma_8m(1) = sqrt(w_8m(:,1)'*sigma8*w_8m(:,1));

% Mean-variance portfolio locus
ngrid=500;
mu_8 = linspace(mu_8m,max(z_bar8),ngrid);
Aeq_8 =[ones(1,48);z_bar8'];   % All weights should sum up
w_8=zeros(48,ngrid);
% solve the optimization
for i = 1:ngrid;
   w_8(:,i)= quadprog(sigma8,f,[],[],Aeq_8,[1;mu_8(i)]);
   sigma_82(i) = sqrt(w_8(:,i)'*sigma8*w_8(:,i));
   sharpe_8(i)=(mu_8(i)-R)/sigma_82(i);
end ;
%scatter(sigma_82,mu_8,'.','b')
sharp8_48=max(sharpe_8);
sharp8_48  %1.2532


% Bendchemark: 48 industrie portfolio with short-sell constraint

A_8 = -1*eye(n,n);     % Besides the returns we forbid short selling
b_8 = zeros(48, 1);     % Required return and weights greater/eqauls 0

w_8ms(:,1) = quadprog(sigma8,f,A_8,b_8,Aeq_8m,beq_8m);
mu_8ms(1) = z_bar8'*w_8ms(:,1);
sigma_8ms(1) = sqrt(w_8ms(:,1)'*sigma8*w_8ms(:,1));
mu_8s = linspace(mu_8ms,max(z_bar8),ngrid);

% solve the optimization
for i = 1:ngrid;
   w_8s(:,i)= quadprog(sigma8,f,A_8,b_8,Aeq_8,[1;mu_8(i)]);
   sigma_8s(i) = sqrt(w_8s(:,i)'*sigma8*w_8s(:,i));
   sharpe_8s(i)=(mu_8s(i)-R)/sigma_8s(i);
end ;
%scatter(sigma_8s,mu_8s,'.','b')
sharp8_48s=max(sharpe_8s);
sharp8_48s  %0.2910

% Bendchemark: 48 industrie portfolio equal weight

w8e=(1/n).*(ones(n,1));
sharpe_8e=((z_bar8'*w8e)-R)/sqrt((w8e'*sigma8*w8e));
sharpe_8e   %0.1038


%% Try 1: 10 industries with largest sharpe ratio from 48

% 10 strong without short-sell constraint

sharpe_8_10=NaN(1,48);
stddataset8=std(dataset8);
for i=1:48;
sharpe_8_10f(i)=(abs(z_bar8(i))-R)/(stddataset8(i));
end
[Bsharpe_8_10f,I] = sort(sharpe_8_10f,'descend'); 
I
strong5=NaN(60,5);
strong5(:,1)=dataset8(:,12); % MedEq 12
strong5(:,2)=dataset8(:,43); % Meals 43
strong5(:,3)=dataset8(:,45); % Insur 45
strong5(:,4)=dataset8(:,26); % Guns  26
strong5(:,5)=dataset8(:,34); % BusSv 34


portret(:,1)=dataset8(:,12);

data47=dataset8;
data47(:,12)=[];
p1=Portfolio('NumAssets',2);
p1.RiskFreeRate=R;
p1=setBudget(p1,1,1);
p1=setBounds(p1,-10,10);
[m,nn]=size(data47);
sr47=NaN(nn,1);
for i=1:nn
    ret=[portret,data47(:,i)];
    p1=p1.estimateAssetMoments(ret);
    wgt=p1.estimateMaxSharpeRatio;
    [prsk,pret]=p1.estimatePortMoments(wgt);
    sr47(i,1)=(pret-R)/prsk;
end
maxsr47=max(sr47);
ind1=find(sr47==maxsr47);  %% var21 of data47 is selected 
portret(:,2)=data47(:,ind1);
data46=data47;
data46(:,ind1)=[];


p2=Portfolio('NumAssets',3);
p2.RiskFreeRate=R;
p2=setBudget(p2,1,1);
p2=setBounds(p2,-10,10);
[m,n1]=size(data46);
sr46=NaN(n1,1);
for i=1:n1
    ret=[portret,data46(:,i)];
    p2=p2.estimateAssetMoments(ret);
    wgt=p2.estimateMaxSharpeRatio;
    [prsk,pret]=p2.estimatePortMoments(wgt);
    sr46(i,1)=(pret-R)/prsk;
end
maxsr46=max(sr46);
ind2=find(sr46==maxsr46); %% var20 of data46 is selected
portret(:,3)=data46(:,ind2);
data45=data46;
data45(:,ind2)=[];

p3=Portfolio('NumAssets',4);
p3.RiskFreeRate=R;
p3=setBudget(p3,1,1);
p3=setBounds(p3,-10,10);
[m,n2]=size(data45);
sr45=NaN(n2,1);
for i=1:n2
    ret=[portret,data45(:,i)];
    p3=p3.estimateAssetMoments(ret);
    wgt=p3.estimateMaxSharpeRatio;
    [prsk,pret]=p3.estimatePortMoments(wgt);
    sr45(i,1)=(pret-R)/prsk;
end
maxsr45=max(sr45);
ind3=find(sr45==maxsr45);
portret(:,4)=data45(:,ind3);
data44=data45;
data44(:,ind3)=[];

p4=Portfolio('NumAssets',5);
p4.RiskFreeRate=R;
p4=setBudget(p4,1,1);
p4=setBounds(p4,-10,10);
[m,n3]=size(data44);
sr44=NaN(n3,1);
for i=1:n3
    ret=[portret,data44(:,i)];
    p4=p4.estimateAssetMoments(ret);
    wgt=p4.estimateMaxSharpeRatio;
    [prsk,pret]=p4.estimatePortMoments(wgt);
    sr44(i,1)=(pret-R)/prsk;
end
maxsr44=max(sr44);
ind4=find(sr44==maxsr44);
portret(:,5)=data44(:,ind4);
data43=data44;
data43(:,ind4)=[];

p5=Portfolio('NumAssets',6);
p5.RiskFreeRate=R;
p5=setBudget(p5,1,1);
p5=setBounds(p5,-10,10);
[m,n4]=size(data43);
sr43=NaN(n4,1);
for i=1:n4
    ret=[portret,data43(:,i)];
    p5=p5.estimateAssetMoments(ret);
    wgt=p5.estimateMaxSharpeRatio;
    [prsk,pret]=p5.estimatePortMoments(wgt);
    sr43(i,1)=(pret-R)/prsk;
end
maxsr43=max(sr43);
ind5=find(sr43==maxsr43);
portret(:,6)=data43(:,ind5);
data42=data43;
data42(:,ind5)=[];

p6=Portfolio('NumAssets',7);
p6.RiskFreeRate=R;
p6=setBudget(p6,1,1);
p6=setBounds(p6,-10,10);
[m,n5]=size(data42);
sr42=NaN(n5,1);
for i=1:n5
    ret=[portret,data42(:,i)];
    p6=p6.estimateAssetMoments(ret);
    wgt=p6.estimateMaxSharpeRatio;
    [prsk,pret]=p6.estimatePortMoments(wgt);
    sr42(i,1)=(pret-R)/prsk;
end
maxsr42=max(sr42);
ind6=find(sr42==maxsr42);
portret(:,7)=data42(:,ind6);
data41=data42;
data41(:,ind6)=[];

p7=Portfolio('NumAssets',8);
p7.RiskFreeRate=R;
p7=setBudget(p7,1,1);
p7=setBounds(p7,-10,10);
[m,n6]=size(data41);
sr41=NaN(n6,1);
for i=1:n6
    ret=[portret,data41(:,i)];
    p7=p7.estimateAssetMoments(ret);
    wgt=p7.estimateMaxSharpeRatio;
    [prsk,pret]=p7.estimatePortMoments(wgt);
    sr41(i,1)=(pret-R)/prsk;
end
maxsr41=max(sr41);
ind7=find(sr41==maxsr41);
portret(:,8)=data41(:,ind7);
data40=data41;
data40(:,ind7)=[];

p8=Portfolio('NumAssets',9);
p8.RiskFreeRate=R;
p8=setBudget(p8,1,1);
p8=setBounds(p8,-10,10);
[m,n7]=size(data40);
sr40=NaN(n7,1);
for i=1:n7
    ret=[portret,data40(:,i)];
    p8=p8.estimateAssetMoments(ret);
    wgt=p8.estimateMaxSharpeRatio;
    [prsk,pret]=p8.estimatePortMoments(wgt);
    sr40(i,1)=(pret-R)/prsk;
end
maxsr40=max(sr40);
ind8=find(sr40==maxsr40);
portret(:,9)=data40(:,ind8);
data39=data40;
data39(:,ind8)=[];

p9=Portfolio('NumAssets',10);
p9.RiskFreeRate=R;
p9=setBudget(p9,1,1);
p9=setBounds(p9,-10,10);
[m,n8]=size(data39);
sr39=NaN(n8,1);
for i=1:n8
    ret=[portret,data39(:,i)];
    p9=p9.estimateAssetMoments(ret);
    wgt=p9.estimateMaxSharpeRatio;
    [prsk,pret]=p9.estimatePortMoments(wgt);
    sr39(i,1)=(pret-R)/prsk;
end
maxsr39=max(sr39);
ind9=find(sr39==maxsr39);
portret(:,10)=data39(:,ind9);


% The maximum sharpe ratio is 0.9046



% 10 strong with short-sell constraint

data47_ssc=dataset8;
data47_ssc(:,12)=[];
portret_ssc(:,1)=dataset8(:,12);
p1_ssc=Portfolio('NumAssets',2);
p1_ssc=p1_ssc.setDefaultConstraints;
p1_ssc.RiskFreeRate=R;
[m,n_ssc]=size(data47_ssc);
sr47_ssc=NaN(n_ssc,1);

for i=1:n_ssc
    ret=[portret_ssc,data47_ssc(:,i)];
    p1_ssc=p1_ssc.estimateAssetMoments(ret);
    wgt=p1_ssc.estimateMaxSharpeRatio;
    [prsk,pret]=p1_ssc.estimatePortMoments(wgt);
    sr47_ssc(i,1)=(pret-R)/prsk;
end
maxsr47_ssc=max(sr47_ssc);
ind1_ssc=find(sr47_ssc==maxsr47_ssc);  
portret_ssc(:,2)=data47_ssc(:,ind1_ssc);
data46_ssc=data47_ssc;
data46_ssc(:,ind1_ssc)=[];


p2_ssc=Portfolio('NumAssets',3);
p2_ssc=p2_ssc.setDefaultConstraints;
p2_ssc.RiskFreeRate=R;
[m,n1_ssc]=size(data46_ssc);
sr46_ssc=NaN(n1_ssc,1);
for i=1:n1_ssc
    ret=[portret_ssc,data46_ssc(:,i)];
    p2_ssc=p2_ssc.estimateAssetMoments(ret);
    wgt=p2_ssc.estimateMaxSharpeRatio;
    [prsk,pret]=p2_ssc.estimatePortMoments(wgt);
    sr46_ssc(i,1)=(pret-R)/prsk;
end
maxsr46_ssc=max(sr46_ssc);
ind2_ssc=find(sr46_ssc==maxsr46_ssc); 
portret_ssc(:,3)=data46_ssc(:,ind2_ssc);
data45_ssc=data46_ssc;
data45_ssc(:,ind2_ssc)=[];

p3_ssc=Portfolio('NumAssets',4);
p3_ssc=p3_ssc.setDefaultConstraints;
p3_ssc.RiskFreeRate=R;
[m,n2_ssc]=size(data45_ssc);
sr45_ssc=NaN(n2_ssc,1);
for i=1:n2_ssc
    ret=[portret_ssc,data45_ssc(:,i)];
    p3_ssc=p3_ssc.estimateAssetMoments(ret);
    wgt=p3_ssc.estimateMaxSharpeRatio;
    [prsk,pret]=p3_ssc.estimatePortMoments(wgt);
    sr45_ssc(i,1)=(pret-R)/prsk;
end
maxsr45_ssc=max(sr45_ssc);
ind3_ssc=find(sr45_ssc==maxsr45_ssc);
portret_ssc(:,4)=data45_ssc(:,ind3_ssc);
data44_ssc=data45_ssc;
data44_ssc(:,ind3_ssc)=[];

p4_ssc=Portfolio('NumAssets',5);
p4_ssc=p4_ssc.setDefaultConstraints;
p4_ssc.RiskFreeRate=R;
[m,n3_ssc]=size(data44_ssc);
sr44_ssc=NaN(n3_ssc,1);
for i=1:n3_ssc
    ret=[portret_ssc,data44_ssc(:,i)];
    p4_ssc=p4_ssc.estimateAssetMoments(ret);
    wgt=p4_ssc.estimateMaxSharpeRatio;
    [prsk,pret]=p4_ssc.estimatePortMoments(wgt);
    sr44_ssc(i,1)=(pret-R)/prsk;
end
maxsr44_ssc=max(sr44_ssc);
ind4_ssc=find(sr44_ssc==maxsr44_ssc);
portret_ssc(:,5)=data44_ssc(:,ind4_ssc);
data43_ssc=data44_ssc;
data43_ssc(:,ind4_ssc)=[];

p5_ssc=Portfolio('NumAssets',6);
p5_ssc=p5_ssc.setDefaultConstraints;
p5_ssc.RiskFreeRate=R;
[m,n4_ssc]=size(data43_ssc);
sr43_ssc=NaN(n4_ssc,1);
for i=1:n4_ssc
    ret=[portret_ssc,data43_ssc(:,i)];
    p5_ssc=p5_ssc.estimateAssetMoments(ret);
    wgt=p5_ssc.estimateMaxSharpeRatio;
    [prsk,pret]=p5_ssc.estimatePortMoments(wgt);
    sr43_ssc(i,1)=(pret-R)/prsk;
end
maxsr43_ssc=max(sr43_ssc);
ind5_ssc=find(sr43_ssc==maxsr43_ssc);
portret_ssc(:,6)=data43_ssc(:,ind5_ssc);
data42_ssc=data43_ssc;
data42_ssc(:,ind5_ssc)=[];

p6_ssc=Portfolio('NumAssets',7);
p6_ssc=p6_ssc.setDefaultConstraints;
p6_ssc.RiskFreeRate=R;
[m,n5_ssc]=size(data42_ssc);
sr42_ssc=NaN(n5_ssc,1);
for i=1:n5_ssc
    ret=[portret_ssc,data42_ssc(:,i)];
    p6_ssc=p6_ssc.estimateAssetMoments(ret);
    wgt=p6_ssc.estimateMaxSharpeRatio;
    [prsk,pret]=p6_ssc.estimatePortMoments(wgt);
    sr42_ssc(i,1)=(pret-R)/prsk;
end
maxsr42_ssc=max(sr42_ssc);
ind6_ssc=find(sr42_ssc==maxsr42_ssc);
portret_ssc(:,7)=data42_ssc(:,ind6_ssc);
data41_ssc=data42_ssc;
data41_ssc(:,ind6_ssc)=[];

p7_ssc=Portfolio('NumAssets',8);
p7_ssc=p7_ssc.setDefaultConstraints;
p7_ssc.RiskFreeRate=R;
[m,n6_ssc]=size(data41_ssc);
sr41_ssc=NaN(n6_ssc,1);
for i=1:n6_ssc
    ret=[portret_ssc,data41_ssc(:,i)];
    p7_ssc=p7_ssc.estimateAssetMoments(ret);
    wgt=p7_ssc.estimateMaxSharpeRatio;
    [prsk,pret]=p7_ssc.estimatePortMoments(wgt);
    sr41_ssc(i,1)=(pret-R)/prsk;
end
maxsr41_ssc=max(sr41_ssc);
ind7_ssc=find(sr41_ssc==maxsr41_ssc);
portret_ssc(:,8)=data41_ssc(:,ind7_ssc);
data40_ssc=data41_ssc;
data40_ssc(:,ind7_ssc)=[];

p8_ssc=Portfolio('NumAssets',9);
p8_ssc=p8_ssc.setDefaultConstraints;
p8_ssc.RiskFreeRate=R;
[m,n7_ssc]=size(data40_ssc);
sr40_ssc=NaN(n7_ssc,1);
for i=1:n7_ssc
    ret=[portret_ssc,data40_ssc(:,i)];
    p8_ssc=p8_ssc.estimateAssetMoments(ret);
    wgt=p8_ssc.estimateMaxSharpeRatio;
    [prsk,pret]=p8_ssc.estimatePortMoments(wgt);
    sr40_ssc(i,1)=(pret-R)/prsk;
end
maxsr40_ssc=max(sr40_ssc);
ind8_ssc=find(sr40_ssc==maxsr40_ssc);
portret_ssc(:,9) = data40_ssc (:,ind8_ssc);
data39_ssc=data40_ssc;
data39_ssc(:,ind8_ssc)=[];

p9_ssc=Portfolio('NumAssets',10);
p9_ssc=p9_ssc.setDefaultConstraints;
p9_ssc.RiskFreeRate=R;
[m,n8_ssc]=size(data39_ssc);
sr39_ssc=NaN(n8_ssc,1);
for i=1:n8_ssc
    ret=[portret_ssc,data39_ssc(:,i)];
    p9_ssc=p9_ssc.estimateAssetMoments(ret);
    wgt=p9_ssc.estimateMaxSharpeRatio;
    [prsk,pret]=p9_ssc.estimatePortMoments(wgt);
    sr39_ssc(i,1)=(pret-R)/prsk;
end
maxsr39_ssc=max(sr39_ssc);
ind9_ssc=find(sr39_ssc==maxsr39_ssc);
portret_ssc(:,10)=data39_ssc(:,ind9_ssc);

% the maximum sharpe ratio is 0.3389




%% Try 2: 10 industries with highest weight in 48 industries portfolios

% Without short-sell constraint (use absolue weight)
id = find(sharpe_8==sharp8_48);
w_48=w_8(:,id);

[Bw_48,k] = sort(abs(w_48),'descend'); 
k
wh5=NaN(60,5);
wh5(:,1)=dataset8(:,4); % Beer 4
wh5(:,2)=dataset8(:,3); % soda 3
wh5(:,3)=dataset8(:,47); % Fin  47
wh5(:,4)=dataset8(:,15); % Rubbr 15
wh5(:,5)=dataset8(:,21); % Mach  21


portretw(:,1)=dataset8(:,4);
data47w=dataset8;
data47w(:,4)=[];
p1w=Portfolio('NumAssets',2);
p1w.RiskFreeRate=R;
p1w=setBudget(p1w,1,1);
p1w=setBounds(p1w,-10,10);
[m,nw]=size(data47w);
sr47w=NaN(nw,1);
for i=1:nw
    ret=[portretw,data47w(:,i)];
    p1w=p1w.estimateAssetMoments(ret);
    wgt=p1w.estimateMaxSharpeRatio;
    [prsk,pret]=p1w.estimatePortMoments(wgt);
    sr47w(i,1)=(pret-R)/prsk;
end
maxsr47w=max(sr47w);
ind1w=find(sr47w==maxsr47w);  %% Medeq12
portretw(:,2)=data47w(:,ind1w);
data46w=data47w;
data46w(:,ind1w)=[];


p2w=Portfolio('NumAssets',3);
p2w.RiskFreeRate=R;
p2w=setBudget(p2w,1,1);
p2w=setBounds(p2w,-10,10);
[m,n1w]=size(data46w);
sr46w=NaN(n1w,1);
for i=1:n1w
    ret=[portretw,data46w(:,i)];
    p2w=p2w.estimateAssetMoments(ret);
    wgt=p2w.estimateMaxSharpeRatio;
    [prsk,pret]=p2w.estimatePortMoments(wgt);
    sr46w(i,1)=(pret-R)/prsk;
end
maxsr46w=max(sr46w);
ind2w=find(sr46w==maxsr46w); %% RLEst46
portretw(:,3)=data46w(:,ind2w);
data45w=data46w;
data45w(:,ind2w)=[];

p3w=Portfolio('NumAssets',4);
p3w.RiskFreeRate=R;
p3w=setBudget(p3w,1,1);
p3w=setBounds(p3w,-10,10);
[m,n2w]=size(data45w);
sr45w=NaN(n2w,1);
for i=1:n2w
    ret=[portretw,data45w(:,i)];
    p3w=p3w.estimateAssetMoments(ret);
    wgt=p3w.estimateMaxSharpeRatio;
    [prsk,pret]=p3w.estimatePortMoments(wgt);
    sr45w(i,1)=(pret-R)/prsk;
end
maxsr45w=max(sr45w);
ind3w=find(sr45w==maxsr45w); % Insur45
portretw(:,4)=data45w(:,ind3w);
data44w=data45w;
data44w(:,ind3w)=[];

p4w=Portfolio('NumAssets',5);
p4w.RiskFreeRate=R;
p4w=setBudget(p4w,1,1);
p4w=setBounds(p4w,-10,10);
[m,n3w]=size(data44w);
sr44w=NaN(n3w,1);
for i=1:n3w
    ret=[portretw,data44w(:,i)];
    p4w=p4w.estimateAssetMoments(ret);
    wgt=p4w.estimateMaxSharpeRatio;
    [prsk,pret]=p4w.estimatePortMoments(wgt);
    sr44w(i,1)=(pret-R)/prsk;
end
maxsr44w=max(sr44w);
ind4w=find(sr44w==maxsr44w); % Toys6
portretw(:,5)=data45w(:,ind4w);
data43w=data44w;
data43w(:,ind4w)=[];

p5w=Portfolio('NumAssets',6);
p5w.RiskFreeRate=R;
p5w=setBudget(p5w,1,1);
p5w=setBounds(p5w,-10,10);
[m,n4w]=size(data43w);
sr43w=NaN(n4w,1);
for i=1:n4w
    ret=[portretw,data43w(:,i)];
    p5w=p5w.estimateAssetMoments(ret);
    wgt=p5w.estimateMaxSharpeRatio;
    [prsk,pret]=p5w.estimatePortMoments(wgt);
    sr43w(i,1)=(pret-R)/prsk;
end
maxsr43w=max(sr43w);
ind5w=find(sr43w==maxsr43w);
portretw(:,6)=data43w(:,ind5w);
data42w=data43w;
data42w(:,ind5w)=[];

p6w=Portfolio('NumAssets',7);
p6w.RiskFreeRate=R;
p6w=setBudget(p6w,1,1);
p6w=setBounds(p6w,-10,10);
[m,n5w]=size(data42w);
sr42w=NaN(n5w,1);
for i=1:n5w
    ret=[portretw,data42w(:,i)];
    p6w=p6w.estimateAssetMoments(ret);
    wgt=p6w.estimateMaxSharpeRatio;
    [prsk,pret]=p6w.estimatePortMoments(wgt);
    sr42w(i,1)=(pret-R)/prsk;
end
maxsr42w=max(sr42w);
ind6w=find(sr42w==maxsr42w);
portretw(:,7)=data42w(:,ind6w);
data41w=data42w;
data41w(:,ind6w)=[];

p7w=Portfolio('NumAssets',8);
p7w.RiskFreeRate=R;
p7w=setBudget(p7w,1,1);
p7w=setBounds(p7w,-10,10);
[m,n6w]=size(data41w);
sr41w=NaN(n6w,1);
for i=1:n6w
    ret=[portretw,data41w(:,i)];
    p7w=p7w.estimateAssetMoments(ret);
    wgt=p7w.estimateMaxSharpeRatio;
    [prsk,pret]=p7w.estimatePortMoments(wgt);
    sr41w(i,1)=(pret-R)/prsk;
end
maxsr41w=max(sr41w);
ind7w=find(sr41w==maxsr41w);
portretw(:,8)=data41w(:,ind7w);
data40w=data41w;
data40w(:,ind7w)=[];

p8w=Portfolio('NumAssets',9);
p8w.RiskFreeRate=R;
p8w=setBudget(p8w,1,1);
p8w=setBounds(p8w,-10,10);
[m,n7w]=size(data40w);
sr40w=NaN(n7w,1);
for i=1:n7w
    ret=[portretw,data40w(:,i)];
    p8w=p8w.estimateAssetMoments(ret);
    wgt=p8w.estimateMaxSharpeRatio;
    [prsk,pret]=p8w.estimatePortMoments(wgt);
    sr40w(i,1)=(pret-R)/prsk;
end
maxsr40w=max(sr40w);
ind8w=find(sr40w==maxsr40w);
portretw(:,9)=data40w(:,ind8w);
data39w=data40w;
data39w(:,ind8w)=[];

p9w=Portfolio('NumAssets',10);
p9w.RiskFreeRate=R;
p9w=setBudget(p9w,1,1);
p9w=setBounds(p9w,-10,10);
[m,n8w]=size(data39w);
sr39w=NaN(n8w,1);
for i=1:n8w
    ret=[portret,data39w(:,i)];
    p9w=p9w.estimateAssetMoments(ret);
    wgt=p9w.estimateMaxSharpeRatio;
    [prsk,pret]=p9w.estimatePortMoments(wgt);
    sr39w(i,1)=(pret-R)/prsk;
end
maxsr39w=max(sr39w);
ind9w=find(sr39w==maxsr39w);
portretw(:,10)=data39w(:,ind9w);
                
% sharpe ration is 0.5241


% with short-sell constraint

ids = find(sharpe_8s==sharp8_48s);
w_48s=w_8s(:,ids);
[Bw_48s,ks] = sort(w_48s,'descend'); 
ks

wh5s=NaN(60,5);
wh5s(:,1)=dataset8(:,12); % MedEq 12
wh5s(:,2)=dataset8(:,36); % Chips 36
wh5s(:,3)=dataset8(:,27); % Gold  27
wh5s(:,4)=dataset8(:,26); % Guns  26
wh5s(:,5)=dataset8(:,45); % Insurance  45


data47_sscw=dataset8;
data47_sscw(:,4)=[];
portret_sscw(:,1)=dataset8(:,4);
p1_sscw=Portfolio('NumAssets',2);
p1_sscw=p1_sscw.setDefaultConstraints;
p1_sscw.RiskFreeRate=R;
[m,n_sscw]=size(data47_sscw);
sr47_sscw=NaN(n_sscw,1);

for i=1:n_sscw
    ret=[portret_sscw,data47_sscw(:,i)];
    p1_sscw=p1_sscw.estimateAssetMoments(ret);
    wgt=p1_sscw.estimateMaxSharpeRatio;
    [prsk,pret]=p1_sscw.estimatePortMoments(wgt);
    sr47_sscw(i,1)=(pret-R)/prsk;
end
maxsr47_sscw=max(sr47_sscw);
ind1_sscw=find(sr47_sscw==maxsr47_sscw);  % MedEq12
portret_sscw(:,2)=data47_sscw(:,ind1_sscw);
data46_sscw=data47_sscw;
data46_sscw(:,ind1_sscw)=[];


p2_sscw=Portfolio('NumAssets',3);
p2_sscw=p2_sscw.setDefaultConstraints;
p2_sscw.RiskFreeRate=R;
[m,n1_sscw]=size(data46_sscw);
sr46_sscw=NaN(n1_sscw,1);
for i=1:n1_sscw
    ret=[portret_sscw,data46_sscw(:,i)];
    p2_sscw=p2_sscw.estimateAssetMoments(ret);
    wgt=p2_sscw.estimateMaxSharpeRatio;
    [prsk,pret]=p2_sscw.estimatePortMoments(wgt);
    sr46_sscw(i,1)=(pret-R)/prsk;
end
maxsr46_sscw=max(sr46_sscw);
ind2_sscw=find(sr46_sscw==maxsr46_sscw); %Meal43
portret_sscw(:,3)=data46_sscw(:,ind2_sscw);
data45_sscw=data46_sscw;
data45_sscw(:,ind2_sscw)=[];

p3_sscw=Portfolio('NumAssets',4);
p3_sscw=p3_sscw.setDefaultConstraints;
p3_sscw.RiskFreeRate=R;
[m,n2_sscw]=size(data45_sscw);
sr45_sscw=NaN(n2_sscw,1);
for i=1:n2_sscw
    ret=[portret_sscw,data45_sscw(:,i)];
    p3_sscw=p3_sscw.estimateAssetMoments(ret);
    wgt=p3_sscw.estimateMaxSharpeRatio;
    [prsk,pret]=p3_sscw.estimatePortMoments(wgt);
    sr45_sscw(i,1)=(pret-R)/prsk;
end
maxsr45_sscw=max(sr45_sscw);
ind3_sscw=find(sr45_sscw==maxsr45_sscw); %Gold27
portret_sscw(:,4)=data45_sscw(:,ind3_sscw);
data44_sscw=data45_sscw;
data44_sscw(:,ind3_sscw)=[];

p4_sscw=Portfolio('NumAssets',5);
p4_sscw=p4_sscw.setDefaultConstraints;
p4_sscw.RiskFreeRate=R;
[m,n3_sscw]=size(data44_sscw);
sr44w_sscw=NaN(n3_sscw,1);
for i=1:n3_sscw
    ret=[portret_sscw,data44_sscw(:,i)];
    p4_sscw=p4_sscw.estimateAssetMoments(ret);
    wgt=p4_sscw.estimateMaxSharpeRatio;
    [prsk,pret]=p4_sscw.estimatePortMoments(wgt);
    sr44_sscw(i,1)=(pret-R)/prsk;
end
maxsr44_sscw=max(sr44_sscw);
ind4_sscw=find(sr44_sscw==maxsr44_sscw); %guns27
portret_sscw(:,5)=data44_sscw(:,ind4_sscw);
data43_sscw=data44_sscw;
data43_sscw(:,ind4_sscw)=[];

p5_sscw=Portfolio('NumAssets',6);
p5_sscw=p5_sscw.setDefaultConstraints;
p5_sscw.RiskFreeRate=R;
[m,n4_sscw]=size(data43_sscw);
sr43_sscw=NaN(n4_sscw,1);
for i=1:n4_sscw
    ret=[portret_sscw,data43_sscw(:,i)];
    p5_sscw=p5_sscw.estimateAssetMoments(ret);
    wgt=p5_sscw.estimateMaxSharpeRatio;
    [prsk,pret]=p5_sscw.estimatePortMoments(wgt);
    sr43_sscw(i,1)=(pret-R)/prsk;
end
maxsr43_sscw=max(sr43_sscw);
ind5_sscw=find(sr43_sscw==maxsr43_sscw);
portret_sscw(:,6)=data43_sscw(:,ind5_sscw);
data42_sscw=data43_sscw;
data42_sscw(:,ind5_sscw)=[];

p6_sscw=Portfolio('NumAssets',7);
p6_sscw=p6_sscw.setDefaultConstraints;
p6_sscw.RiskFreeRate=R;
[m,n5_sscw]=size(data42_sscw);
sr42_sscw=NaN(n5_sscw,1);
for i=1:n5_sscw
    ret=[portret_sscw,data42_sscw(:,i)];
    p6_sscw=p6_sscw.estimateAssetMoments(ret);
    wgt=p6_sscw.estimateMaxSharpeRatio;
    [prsk,pret]=p6_sscw.estimatePortMoments(wgt);
    sr42_sscw(i,1)=(pret-R)/prsk;
end
maxsr42_sscw=max(sr42_sscw);
ind6_sscw=find(sr42_sscw==maxsr42_sscw);
portret_sscw(:,7)=data42_sscw(:,ind6_sscw);
data41_sscw=data42_sscw;
data41_sscw(:,ind6_sscw)=[];

p7_sscw=Portfolio('NumAssets',8);
p7_sscw=p7_sscw.setDefaultConstraints;
p7_sscw.RiskFreeRate=R;
[m,n6_sscw]=size(data41_sscw);
sr41_sscw=NaN(n6_sscw,1);
for i=1:n6_sscw
    ret=[portret_sscw,data41_sscw(:,i)];
    p7_sscw=p7_sscw.estimateAssetMoments(ret);
    wgt=p7_sscw.estimateMaxSharpeRatio;
    [prsk,pret]=p7_sscw.estimatePortMoments(wgt);
    sr41_sscw(i,1)=(pret-R)/prsk;
end
maxsr41_sscw=max(sr41_sscw);
ind7_sscw=find(sr41_sscw==maxsr41_sscw);
portret_sscw(:,8)=data41_sscw(:,ind7_sscw);
data40_sscw=data41_sscw;
data40_sscw(:,ind7_sscw)=[];

p8_sscw=Portfolio('NumAssets',9);
p8_sscw=p8_sscw.setDefaultConstraints;
p8_sscw.RiskFreeRate=R;
[m,n7_sscw]=size(data40_sscw);
sr40_sscw=NaN(n7_sscw,1);
for i=1:n7_sscw
    ret=[portret_sscw,data40_sscw(:,i)];
    p8_sscw=p8_sscw.estimateAssetMoments(ret);
    wgt=p8_sscw.estimateMaxSharpeRatio;
    [prsk,pret]=p8_sscw.estimatePortMoments(wgt);
    sr40_sscw(i,1)=(pret-R)/prsk;
end
maxsr40_sscw=max(sr40_sscw);
ind8_sscw=find(sr40_sscw==maxsr40_sscw);
portret_sscw(:,9)=data40_sscw(:,ind8_sscw);
data39_sscw=data40_sscw;
data39_sscw(:,ind8_sscw)=[];

p9_sscw=Portfolio('NumAssets',10);
p9_sscw=p9_sscw.setDefaultConstraints;
p9_sscw.RiskFreeRate=R;
[m,n8_sscw]=size(data39_sscw);
sr39_sscw=NaN(n8_sscw,1);
for i=1:n8_sscw
    ret=[portret_sscw,data39_sscw(:,i)];
    p9_sscw=p9_sscw.estimateAssetMoments(ret);
    wgt=p9_sscw.estimateMaxSharpeRatio;
    [prsk,pret]=p9_sscw.estimatePortMoments(wgt);
    sr39_sscw(i,1)=(pret-R)/prsk;
end
maxsr39_sscw=max(sr39_sscw);
ind9_sscw=find(sr39_sscw==maxsr39_sscw);
portret_sscw(:,10)=data39_sscw(:,ind9_sscw);

% the maximum sharpe ratio is 0.3341




%% Try 3: Linear regression

% Without short-sell constraint
yl=(dataset8*w_48-R/48)/sigma_82(id); % sharpe ratio from the tangency porfolio of 
% all 48 industries.

modl_r=LinearModel.fit([dataset8(:,1),dataset8(:,2),dataset8(:,3),dataset8(:,4),...
    dataset8(:,5),dataset8(:,6),dataset8(:,7),dataset8(:,8),dataset8(:,9),...
    dataset8(:,10),dataset8(:,11),dataset8(:,12),dataset8(:,13),dataset8(:,14),...
    dataset8(:,15),dataset8(:,16),dataset8(:,17),dataset8(:,18),dataset8(:,19),...
    dataset8(:,20),dataset8(:,21),dataset8(:,22),dataset8(:,23),dataset8(:,24),...
    dataset8(:,25),dataset8(:,26),dataset8(:,27),dataset8(:,28),dataset8(:,29),...
    dataset8(:,30),dataset8(:,31),dataset8(:,32),dataset8(:,33),dataset8(:,34),...
    dataset8(:,35),dataset8(:,36),dataset8(:,37),dataset8(:,38),dataset8(:,39),...
    dataset8(:,40),dataset8(:,41),dataset8(:,42),dataset8(:,43),dataset8(:,44),...
    dataset8(:,45),dataset8(:,46),dataset8(:,47),dataset8(:,48)],yl,'RobustOpts','on');
modl_rc=table2array(modl_r.Coefficients);
coefficient=modl_rc(2:49,1);
[coefficient_2,I1] = sort(coefficient,'descend'); 
I1

% all these industries are significatif at 5%, 10 industries with the hightest
% coefficient beta are : Beer 4, Fin 47, Rubbr 15, Mach 21, Autos 23, Util  31
% Clths 10, Other  48, Mines  28, Guns  26.

r10=NaN(60,10);
r10(:,1)=dataset8(:,4); % Beer 4
r10(:,2)=dataset8(:,47); % Fin 47
r10(:,3)=dataset8(:,15); % Rubbr  15
r10(:,4)=dataset8(:,21); % Mach 21
r10(:,5)=dataset8(:,23); % Autos  23
r10(:,6)=dataset8(:,31); % Util  31
r10(:,7)=dataset8(:,10); % Clths 10
r10(:,8)=dataset8(:,48); % Other  48
r10(:,9)=dataset8(:,28); % Mines  28
r10(:,10)=dataset8(:,26); % Guns  26

sigma8_10r=cov(r10);
z_bar8_10r=(mean(r10))';
Aeq_8_10rm = ones(1,10);
beq_8_10rm = 1;
% Minimum variance portfolio
w_8_10rm(:,1) = quadprog(sigma8_10r,zeros(1,10),[],[],Aeq_8_10rm,beq_8_10rm);
mu_8_10rm = z_bar8_10r'*w_8_10rm(:,1);
sigma_8_rm = sqrt(w_8_10rm(:,1)'*sigma8_10r*w_8_10rm(:,1));

% Mean-variance portfolio locus
mu_8_10r = linspace(mu_8_10rm,max(z_bar8_10r),ngrid);
Aeq_8_10r=[ones(1,10);z_bar8_10r'];   % All weights should sum up
w_8_10r=zeros(10,ngrid);
% solve the optimization
for i = 1:ngrid
   w_8_10r(:,i)= quadprog(sigma8_10r,zeros(10,1),[],[],Aeq_8_10r,[1;mu_8_10r(i)]);
   sigma_8_10r2(i) = sqrt(w_8_10r(:,i)'*sigma8_10r*w_8_10r(:,i));
   sharpe_8_10r(i)=(mu_8_10r(i)-R)/sigma_8_10r2(i);
end 
%scatter(sigma_8_10r2,mu_8_10r,'.','b')
sharp8_10r=max(sharpe_8_10r);
sharp8_10r %0.3554


% With short-sell constraint
yls=(dataset8*w_48s-R/48)/sigma_8s(ids); % sharpe ratio from the tangency porfolio of 
% all48 industries.

modls_r=LinearModel.fit([dataset8(:,1),dataset8(:,2),dataset8(:,3),dataset8(:,4),...
    dataset8(:,5),dataset8(:,6),dataset8(:,7),dataset8(:,8),dataset8(:,9),...
    dataset8(:,10),dataset8(:,11),dataset8(:,12),dataset8(:,13),dataset8(:,14),...
    dataset8(:,15),dataset8(:,16),dataset8(:,17),dataset8(:,18),dataset8(:,19),...
    dataset8(:,20),dataset8(:,21),dataset8(:,22),dataset8(:,23),dataset8(:,24),...
    dataset8(:,25),dataset8(:,26),dataset8(:,27),dataset8(:,28),dataset8(:,29),...
    dataset8(:,30),dataset8(:,31),dataset8(:,32),dataset8(:,33),dataset8(:,34),...
    dataset8(:,35),dataset8(:,36),dataset8(:,37),dataset8(:,38),dataset8(:,39),...
    dataset8(:,40),dataset8(:,41),dataset8(:,42),dataset8(:,43),dataset8(:,44),...
    dataset8(:,45),dataset8(:,46),dataset8(:,47),dataset8(:,48)],yls,'RobustOpts','on')

modls_rc=table2array(modls_r.Coefficients);
coefficient_s=modls_rc(2:49,1);
[coefficient_s2,I2] = sort(coefficient_s,'descend'); 
I2
% all these industries are significatif at 5%, 10 industries with the hightest
% coefficient beta are the same with the 10 industries with highest weight in 48
% industries portfolios. But just the 3 first's beta coefficient are quantitatively
% considerable. MedEq 12, Chips 36, Gold 27, Guns 26, Insurance 45, Meals 43,
% BusSv 34, Aero  24, % Beer  4, % LabEq  37.


%% try 4 cluster tree

% Builds a multilevel hierarchy of clusters by creating a cluster tree

% Compute 48 industries data using Ward linkage 
Z = linkage(dataset8','average','chebychev');
T = cluster(Z,'maxclust',10);
% T gives us index, each inustry grouped in 1~10

% Find industries grouped in each group. and take the mean of each row of
% all the composant in each goupe.
[G1ci]=find(T==1); %19
G1c=zeros(60,1);
for i=1:length(G1ci)
    G1c=cat(2,G1c,dataset8(:,G1ci(i)));
    G1c=mean(G1c(:,2:end),2);
end
[G2ci]=find(T==2); %25
 G2c=zeros(60,1);
for i=1:length(G2ci)
    G2c=cat(2,G2c,dataset8(:,G2ci(i)));
    G2c=mean(G2c(:,2:end),2);
end
[G3ci]=find(T==3); %1，7
G3c=zeros(60,1);
for i=1:length(G3ci)
    G3c=cat(2,G3c,dataset8(:,G3ci(i)));
    G3c=mean(G3c(:,2:end),2);
end
[G4ci]=find(T==4); 
%2;3;4;6;8;9;10;11;12;13;14;15;17;18;21;22;23;24;26;30;
%31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;
G4c=zeros(60,1);
for i=1:length(G4ci)
    G4c=cat(2,G4c,dataset8(:,G4ci(i)));
    G4c=mean(G4c(:,2:end),2);
end
[G5ci]=find(T==5); %28
G5c=zeros(60,1);
for i=1:length(G5ci)
    G5c=cat(2,G5c,dataset8(:,G5ci(i)));
    G5c=mean(G5c(:,2:end),2);
end
[G6ci]=find(T==6); %5
G6c=zeros(60,1);
for i=1:length(G6ci)
    G6c=cat(2,G6c,dataset8(:,G6ci(i)));
    G6c=mean(G6c(:,2:end),2);
end
[G7ci]=find(T==7); %16
G7c=zeros(60,1);
for i=1:length(G7ci)
    G7c=cat(2,G7c,dataset8(:,G7ci(i)));
    G7c=mean(G7c(:,2:end),2);
end
[G8ci]=find(T==8); %20
G8c=zeros(60,1);
for i=1:length(G8ci)
    G8c=cat(2,G8c,dataset8(:,G8ci(i)));
    G8c=mean(G8c(:,2:end),2);
end
[G9ci]=find(T==9); %27
G9c=zeros(60,1);
for i=1:length(G9ci)
    G9c=cat(2,G9c,dataset8(:,G9ci(i)));
    G9c=mean(G9c(:,2:end),2);
end
[G10ci]=find(T==10); %29
G10c=zeros(60,1);
for i=1:length(G10ci)
    G10c=cat(2,G10c,dataset8(:,G10ci(i)));
    G10c=mean(G10c(:,2:end),2);
end

G3c_cat=NaN(60,length(G3ci)+1);
G3c_cat(:,1)=G3c;
for i=1:length(G3ci)
    G3c_cat(:,i+1)=dataset8(:,G3ci(i));
end
G4c_cat=NaN(60,length(G4ci)+1);
G4c_cat(:,1)=G4c;
for i=1:length(G4ci)
    G4c_cat(:,i+1)=dataset8(:,G4ci(i));
end

% Choose the represent of G3,G4 by correlation of the mean of G2 and other
% composant.  As groupe 1,2,5,6,7,8,9,10 only have one industry composant, just
% choose themselves.

G3c_crr=corr(G3c_cat);
[corrG3c_2,ic3]=sort(G3c_crr(:,1),'descend');
ic3  %industrie 45

G4c_crr=corr(G4c_cat);
[corrG4c_2,ic4]=sort(G4c_crr(:,1),'descend');
ic4

% Prepare data to optimization
c10=NaN(60,10);
c10(:,1)=dataset8(:,19); %Steel 19
c10(:,2)=dataset8(:,25); %Ships 25
c10(:,3)=dataset8(:,7);  %Fun 7
c10(:,4)=dataset8(:,48); %Other 48
c10(:,5)=dataset8(:,28); %Mines 28
c10(:,6)=dataset8(:,5);  %Smoke 5
c10(:,7)=dataset8(:,16); %Txtls 16
c10(:,8)=dataset8(:,20); %FabPr 20
c10(:,9)=dataset8(:,27); %Gold 27
c10(:,10)=dataset8(:,29);%Coal 29

sigma8_10c=cov(c10);
z_bar8_10c=(mean(c10))';
Aeq_8_10cm = ones(1,10);
beq_8_10cm = 1;
% Minimum variance portfolio
w_8_10cm(:,1) = quadprog(sigma8_10c,[],[],[],Aeq_8_10cm,beq_8_10cm);
mu_8_10cm = z_bar8_10c'*w_8_10cm(:,1);
sigma_8_cm = sqrt(w_8_10cm(:,1)'*sigma8_10c*w_8_10cm(:,1));

% Mean-variance portfolio locus
mu_8_10c = linspace(mu_8_10cm,max(z_bar8_10c),ngrid);
Aeq_8_10c=[ones(1,10);z_bar8_10c'];   % All weights should sum up
w_8_10c=zeros(10,ngrid);
% solve the optimization
for i = 1:ngrid
   w_8_10c(:,i)= quadprog(sigma8_10c,[],[],[],Aeq_8_10c,[1;mu_8_10c(i)]);
   sigma_8_10c2(i) = sqrt(w_8_10c(:,i)'*sigma8_10c*w_8_10c(:,i));
   sharpe_8_10c(i)=(mu_8_10c(i)-R)/sigma_8_10c2(i);
end 
%scatter(sigma_8_10c2,mu_8_10c,'.','b')
sharp8_10c=max(sharpe_8_10c);
sharp8_10c  %0.3245


%% try 5 k-means
% Partitions data into k distinct clusters based on distance to the 
% centroid of a cluster

% NOTE: As in k-means method, each industry data is a spot in a
% multidimensional space. And method choose each time randomly a
% start point and calculate the distance with other point, and give 
% groups. As the start point will differ each time, groups change each time
% Code below just a one time k-means. I created another document named
% kmean. At the second part of the document, we use the 'kmeanf' function
% to do 10 to 50 times k-means method, and give us the largest Sharpe ratio.

K8=kmeans(dataset8',10);
% K gives us index, each inustries grouped in 1~5

% Find the industry grouped in each group. and take the mean of each row of
% all the composant in each goup.
[G1ki]=find(K8==1); 
G1k=zeros(60,1);
    for i=1:length(G1ki)
        G1k=cat(2,G1k,dataset8(:,G1ki(i)));
        G1k=mean(G1k(:,2:end),2);
    end
    [G2ki]=find(K8==2);
    G2k=zeros(60,1);
    for i=1:length(G2ki)
        G2k=cat(2,G2k,dataset8(:,G2ki(i)));
        G2k=mean(G2k(:,2:end),2);
    end
    [G3ki]=find(K8==3);
    G3k=zeros(60,1);
    for i=1:length(G3ki)
        G3k=cat(2,G3k,dataset8(:,G3ki(i)));
        G3k=mean(G3k(:,2:end),2);
    end
    [G4ki]=find(K8==4);
    G4k=zeros(60,1);
    for i=1:length(G4ki)
        G4k=cat(2,G4k,dataset8(:,G4ki(i)));
        G4k=mean(G4k(:,2:end),2);
    end
    [G5ki]=find(K8==5);
    G5k=zeros(60,1);
    for i=1:length(G5ki)
        G5k=cat(2,G5k,dataset8(:,G5ki(i)));
        G5k=mean(G5k(:,2:end),2);
    end
    [G6ki]=find(K8==6);
    G6k=zeros(60,1);
    for i=1:length(G6ki)
        G6k=cat(2,G6k,dataset8(:,G6ki(i)));
        G6k=mean(G6k(:,2:end),2);
    end
    [G7ki]=find(K8==7);
    G7k=zeros(60,1);
    for i=1:length(G7ki)
        G7k=cat(2,G7k,dataset8(:,G7ki(i)));
        G7k=mean(G7k(:,2:end),2);
    end
    [G8ki]=find(K8==8);
    G8k=zeros(60,1);
    for i=1:length(G8ki)
        G8k=cat(2,G8k,dataset8(:,G8ki(i)));
        G8k=mean(G8k(:,2:end),2);
    end
    [G9ki]=find(K8==9);
    G9k=zeros(60,1);
    for i=1:length(G9ki)
        G9k=cat(2,G9k,dataset8(:,G9ki(i)));
        G9k=mean(G9k(:,2:end),2);
    end
    [G10ki]=find(K8==10);
    G10k=zeros(60,1);
    for i=1:length(G10ki)
        G10k=cat(2,G10k,dataset8(:,G10ki(i)));
        G10k=mean(G10k(:,2:end),2);
    end
    
    G1k_cat=NaN(60,length(G1ki)+1);
    G1k_cat(:,1)=G1k;
    for i=1:length(G1ki)
        G1k_cat(:,i+1)=dataset8(:,G1ki(i));
    end
    G2k_cat=NaN(60,length(G2ki)+1);
    G2k_cat(:,1)=G2k;
    for i=1:length(G2ki)
        G2k_cat(:,i+1)=dataset8(:,G2ki(i));
    end
    
G3k_cat=NaN(60,length(G3ki)+1);
G3k_cat(:,1)=G3k;
for i=1:length(G3ki)
    G3k_cat(:,i+1)=dataset8(:,G3ki(i));
end
G4k_cat=NaN(60,length(G4ki)+1);
G4k_cat(:,1)=G4k;
for i=1:length(G4ki)
    G4k_cat(:,i+1)=dataset8(:,G4ki(i));
end
G5k_cat=NaN(60,length(G5ki)+1);
G5k_cat(:,1)=G5k;
for i=1:length(G5ki)
    G5k_cat(:,i+1)=dataset8(:,G5ki(i));
end
G6k_cat=NaN(60,length(G6ki)+1);
G6k_cat(:,1)=G6k;
for i=1:length(G6ki)
    G6k_cat(:,i+1)=dataset8(:,G6ki(i));
end
G7k_cat=NaN(60,length(G7ki)+1);
G7k_cat(:,1)=G7k;
for i=1:length(G7ki)
    G7k_cat(:,i+1)=dataset8(:,G7ki(i));
end
G8k_cat=NaN(60,length(G8ki)+1);
G8k_cat(:,1)=G8k;
for i=1:length(G8ki)
    G8k_cat(:,i+1)=dataset8(:,G8ki(i));
end
G9k_cat=NaN(60,length(G9ki)+1);
G9k_cat(:,1)=G9k;
for i=1:length(G9ki)
    G9k_cat(:,i+1)=dataset8(:,G9ki(i));
end
G10k_cat=NaN(60,length(G10ki)+1);
G10k_cat(:,1)=G10k;
for i=1:length(G10ki)
    G10k_cat(:,i+1)=dataset8(:,G10ki(i));
end

% choose the represent of G1~G10 by correlation of the mean of group and other
% composant.
G1k_crr=corr(G1k_cat);
[corrG1k_2,ik1]=sort(G1k_crr(:,1),'descend');
rep1=dataset8(:,G1ki(ik1(2)-1));

G2k_crr=corr(G2k_cat);
[corrG2k_2,ik2]=sort(G2k_crr(:,1),'descend');
rep2=dataset8(:,G2ki(ik2(2)-1));

G3k_crr=corr(G3k_cat);
[corrG3k_2,ik3]=sort(G3k_crr(:,1),'descend');
rep3=dataset8(:,G3ki(ik3(2)-1));

G4k_crr=corr(G4k_cat);
[corrG4k_2,ik4]=sort(G4k_crr(:,1),'descend');
rep4=dataset8(:,G4ki(ik4(2)-1));

G5k_crr=corr(G5k_cat);
[corrG5k_2,ik5]=sort(G5k_crr(:,1),'descend');
rep5=dataset8(:,G5ki(ik5(2)-1));

G6k_crr=corr(G6k_cat);
[corrG6k_2,ik6]=sort(G6k_crr(:,1),'descend');
rep6=dataset8(:,G6ki(ik6(2)-1));

G7k_crr=corr(G7k_cat);
[corrG7k_2,ik7]=sort(G7k_crr(:,1),'descend');
rep7=dataset8(:,G7ki(ik7(2)-1));

G8k_crr=corr(G8k_cat);
[corrG8k_2,ik8]=sort(G8k_crr(:,1),'descend');
rep8=dataset8(:,G8ki(ik8(2)-1));

G9k_crr=corr(G9k_cat);
[corrG9k_2,ik9]=sort(G9k_crr(:,1),'descend');
rep9=dataset8(:,G9ki(ik9(2)-1));

G10k_crr=corr(G10k_cat);
[corrG10k_2,ik10]=sort(G10k_crr(:,1),'descend');
rep10=dataset8(:,G10ki(ik10(2)-1));

% prepare data to optimization
k10=NaN(60,10);
k10(:,1)=rep1; 
k10(:,2)=rep2;
k10(:,3)=rep3; 
k10(:,4)=rep4; 
k10(:,5)=rep5;
k10(:,6)=rep6; 
k10(:,7)=rep7;
k10(:,8)=rep8; 
k10(:,9)=rep9; 
k10(:,10)=rep10;


sigma8_10k=cov(k10);
z_bar8_10k=(mean(k10))';
Aeq_8_10km = ones(1,10);
beq_8_10km = 1;
% Minimum variance portfolio
w_8_10km(:,1) = quadprog(sigma8_10k,zeros(1,10),[],[],Aeq_8_10km,beq_8_10km);
mu_8_10km = z_bar8_10k'*w_8_10km(:,1);
sigma_8_km = sqrt(w_8_10km(:,1)'*sigma8_10k*w_8_10km(:,1));

% Mean-variance portfolio locus
mu_8_10k = linspace(mu_8_10km,max(z_bar8_10k),ngrid);
Aeq_8_10k=[ones(1,10);z_bar8_10k'];   % All weights should sum up
w_8_10k=zeros(10,ngrid);
% solve the optimization
for i = 1:ngrid
   w_8_10k(:,i)= quadprog(sigma8_10k,zeros(10,1),[],[],Aeq_8_10k,[1;mu_8_10k(i)]);
   sigma_8_10k2(i) = sqrt(w_8_10k(:,i)'*sigma8_10k*w_8_10k(:,i));
   sharpe_8_10k(i)=(mu_8_10k(i)-R)/sigma_8_10k2(i);
end 
scatter(sigma_8_10k2,mu_8_10k,'.','b')
sharp8_10k=max(sharpe_8_10k);
sharp8_10k   %0.3144




%% try 6 

clear;clc;
%*****************************************

% 8) 
% As we only want 5 indusduries from 48 indusries, so we foresse 43
% loops(48~5),each loop we will exclude 1 industry.In the first loop
% we exclude 1 industry by ordre from 1 to 48, and we caculate 48 times
% the left 47 industries portfolio,and we will choose the one who's left 
% giving the max sharp ratio to exclude,our 48 industries become 47, then 
% we do second loop, same way as before, the one who's left giving the max 
% sharp ratio exculded...After 43 loop, we have 5 industries left.

%% Without short-sell constraint

data8= xlsread('Q8_48_Industry_Portfolios.xlsm');
dateline8 = data8(:,1);
dataset8=data8(1:end,2:end)./100;

dataset81=dataset8;
idset=NaN(1,43);

for j=1:44
   for i=1:49-j
       dataset8r=dataset81;
       dataset8r(:,i)=[];
       dataset8_1=dataset8r; 
       sp(i)=Q8f(dataset8_1);
   end 
id=find(sp==max(sp)); 
dataset81(:,id)=[]; 
idset(j)=id; %by idset we know which industry is out each time
sp=[];
end 

% 1Fun,2ships,3Telcom,4Banks,5smoke,6BusSv,7Coal,8Comps,9Books,10Guns,11FabPr,
% 12Agric,13BldMt,14Hlth,15Aero,16FIn,17Trans,18Gold,19Whlsl,20Util,21PerSv,
% 22Rtail,23LabEQ,24Other,25Paper,26Food,27Autos,28Meals,29MedEq,30Chems,
% 31Oil,32Cnstr?33Mines,34steel,35Boxes,36Hshld,37Rubbr,38Soda,39Beer,
% 40Txtls,41Clths,43Mach in turn out.

% 10 strong: 4Beer,6Toys,10Clths,13Drugs,16Txtls,21Mach,22ElcEq,36Chips,
% 45Insur,46RlEst
% 5 strong:6Toys,22ElcEq,36Chips,45Insur,46RlEst

D8_5=cat(2,dataset8(:,6),dataset8(:,22),dataset8(:,36),...
    dataset8(:,45),dataset8(:,46));

D8_5_sharp=Q8f(D8_5);
D8_5_sharp  %0.6374

D8_10=cat(2,dataset8(:,4),dataset8(:,6),dataset8(:,10),dataset8(:,13),...
    dataset8(:,16),dataset8(:,21),dataset8(:,22),dataset8(:,36),...
    dataset8(:,45),dataset8(:,46));

D8_10_sharp=Q8f(D8_10);
D8_10_sharp  %0.7921


%% Without short-sell constraint
clear;clc;
data8= xlsread('Q8_48_Industry_Portfolios.xlsm');
dateline8 = data8(:,1);
dataset8=data8(1:end,2:end)./100;

dataset82=dataset8;
idset=NaN(1,43);

for j=1:44
   for i=1:49-j
       dataset8r=dataset82;
       dataset8r(:,i)=[];
       dataset8_2=dataset8r; 
       sp(i)=Q8fsquad2(dataset8_2);
   end 
id=find(sp==max(sp)); 
dataset82(:,id)=[]; 
idset(j)=id;
sp=[];
end 

% 10 strong with short sell, 2Food, 12MedEq, 17BldMt,26Guns,27Gold,31Util,36Chips,40Trans,
% 43 Meals, 45Insur  sharp ratio %0.3389
% 5 strong 12MedEq,17BldMt,26Guns,31Util,40Trans,45Insur

D8_5s=cat(2,dataset8(:,12),dataset8(:,17),dataset8(:,26),dataset8(:,31),...
   dataset8(:,40),dataset8(:,45));

D8_5_sharps=Q8fsquad2(D8_5s);
D8_5_sharps  %0.3175


D8_10s=cat(2,dataset8(:,2),dataset8(:,12),dataset8(:,17),dataset8(:,26),...
    dataset8(:,27),dataset8(:,31),dataset8(:,36),dataset8(:,40),...
    dataset8(:,43),dataset8(:,45));

D8_10_sharps=Q8fsquad2(D8_10s);
D8_10_sharps  %0.3389
