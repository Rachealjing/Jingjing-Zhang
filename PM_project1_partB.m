%*************************************************************************

% PART B

clear; clc;

% DATA:
% Monthly data on 10 Industry Portfolios and risk-free data from Fama-French
% 3 factorsfrom Ken French's website at:
% http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html.
% Data from 1926.07 to 2018.12.


%*************************************************************************

% 1) the portfolio that maximizes the Sharpe ratio without short-sale
% constraints

data_B = xlsread('10_Industry_Portfolios(part b).xlsm')./100;
data_riskfree = xlsread('risk-free asset.xlsm')./100;
n = length(data_B);
p_ret = NaN(n-60,7);

A = ones(1,10);
b = 1;
for i = 1:(n-60);
    ret(i,:) = data_B(i+60,2:11);
    p_asset = data_B(i:i+59,2:11);
    RiskFreeRate = data_riskfree(i+59,2);
    p_1 = Portfolio('NumAssets',10,'UpperBound', 1);
    p_1 = p_1.estimateAssetMoments(p_asset);
    p_1 = p_1.setEquality(A, b);
    [p_1_wgt(:,i)] = p_1.estimateMaxSharpeRatio;
    p_ret(i,1) = ret(i,:)*p_1_wgt(:,i);
end


%*************************************

% 2) the portfolio that maximizes the Sharpe ratio with short-sale
% constraints

for i = 1:(n-60);
    ret(i,:) = data_B(i+60,2:11);
    p_asset = data_B(i:i+59,2:11);
    RiskFreeRate = data_riskfree(i+59,2);
    p_2 = Portfolio('NumAssets',10,'LowerBound', 0,'UpperBound', 1);
    p_2 = p_2.estimateAssetMoments(p_asset);
    p_2 = p_2.setDefaultConstraints;
    [p_2_wgt(:,i)] = p_2.estimateMaxSharpeRatio;
    p_ret(i,2) = ret(i,:)*p_2_wgt(:,i);
end


%*************************************

% 3) the portfolio where the weight of each asset is inversely related to
% its variance

for i = 1:(n-60);
    ret(i,:) = data_B(i+60,2:11);
    p_asset = data_B(i:i+59,2:11);
    var_3(i,:) = 1./(var(p_asset));
    p_3_wgt(:,i) = (var_3(i,:)./sum(var_3(i,:)))';
    p_ret(i,3) = ret(i,:)*p_3_wgt(:,i);
end


%*************************************

% 4) the portfolio where the weight of each asset is inversely related to
% its volatility

for i = 1:(n-60);
    ret(i,:) = data_B(i+60,2:11);
    p_asset = data_B(i:i+59,2:11);
    std_4(i,:) = 1./(std(p_asset));
    p_4_wgt(:,i) = (std_4(i,:)./sum(std_4(i,:)))';
    p_ret(i,4) = ret(i,:)*p_4_wgt(:,i);
end


%*************************************

% 5) the portfolio where assets have the same weight

for i = 1:(n-60);
    ret(i,:) = data_B(i+60,2:11);
    p_asset = data_B(i:i+59,2:11);
    p_5_wgt(:,i) = (1/10)*ones(10,1);
    p_ret(i,5) = ret(i,:)*p_5_wgt(:,i);
end


%*************************************

% 6) the portfolio where the weight of each is linearly related to
% its market capitalization

data_firmsize = xlsread('PartB_Firm_Size.xlsx');
data_nbfirm = xlsread('PartB_Nb_Firms.xlsm');

for i = 1:(n-60);
    ret(i,:) = data_B(i+60,2:11);
    p_asset = data_B(i:i+59,2:11);
    Market_capital(i,:) = data_firmsize(i+59,2:11).*data_nbfirm(i+59,2:11);
    p_6_wgt(:,i) = (Market_capital(i,:)./sum( Market_capital(i,:)))';
    p_ret(i,6) = ret(i,:)*p_6_wgt(:,i);
end
    

%*************************************

% 7) the portfolio with the minimum variance

for i = 1:(n-60);
    ret(i,:) = data_B(i+60,2:11);
    p_asset = data_B(i:i+59,2:11);
    RiskFreeRate = data_riskfree(i+59,2);
    p_7 = Portfolio('NumAssets',10,'UpperBound', 1);
    p_7 = p_7.estimateAssetMoments(p_asset);
    p_7 = p_7.setEquality(A, b);
    p_7_wgt(:,i)= estimateFrontierLimits(p_7,'Min');
    p_ret(i,7) = ret(i,:)*p_7_wgt(:,i);
end

save('p1_7_ret','p_ret');

%*************************************

%% Compare and explication

p_ret2=importdata('p1_7_ret.mat');
time=data_riskfree(:,1);
names_P ={'P1snss','P2sss','P3v','P4vol','P5eqw','P6mc','P7mv'};

% Period July 1931 - October 2017
ind_upp1=find(data_B(:,1)==193107);  
ind_low=find(data_B(:,1)==201710);  
p_ret_31_17=p_ret2((ind_upp1-60):(ind_low-60),:);
rf_31_17=mean(data_riskfree(ind_upp1:ind_low,2));
p_ret_31_17_average=mean(p_ret_31_17);
p_ret_31_17_total=sum(p_ret_31_17);
p_ret_31_17_sharpe=(p_ret_31_17_average-rf_31_17)./std(p_ret_31_17);

figure(9)
subplot(2,2,1);
plot(time(ind_upp1:ind_low),p_ret_31_17);
xlabel('Portfolio 1-7','fontsize',14);
ylabel('Return','fontsize',14);
legend(names_P,'fontsize', 4,'Location','SouthEast');
title('Return for 1931/07-2017/10','fontsize',14);
subplot(2,2,2);
bar(p_ret_31_17_average);
%plot(names_P,p_ret_31_17_average);
set(gca,'XTickLabel',names_P)
xlabel('Portfolio 1-7','fontsize',14);
ylabel('Average return','fontsize',14);
title('Average return for 1931/07-2017/10','fontsize',14);
subplot(2,2,3);
bar(p_ret_31_17_total);
%plot(names_P,p_ret_31_17_total);
set(gca,'XTickLabel',names_P)
xlabel('Portfolio 1-7','fontsize',14);
ylabel('Total return','fontsize',14);
title('Total return for 1931/07-2017/10','fontsize',14);
subplot(2,2,4);
bar(p_ret_31_17_sharpe);
%plot(names_P,p_ret_31_17_sharpe);
set(gca,'XTickLabel',names_P)
xlabel('Portfolio 1-7','fontsize',14);
ylabel('Sharpe ratio','fontsize',14);
title('Sharpe ratio for 1931/07-2017/10','fontsize',14);

% Period January 1990 - October 2017
ind_upp2=find(data_B(:,1)==199001);   
p_ret_90_17=p_ret2((ind_upp2-60):(ind_low-60),:);
rf_90_17=mean(data_riskfree(ind_upp2:ind_low,2));
p_ret_90_17_average=mean(p_ret_90_17);
p_ret_90_17_total=sum(p_ret_90_17);
p_ret_90_17_sharpe=(p_ret_90_17_average-rf_90_17)./std(p_ret_90_17);

figure(10)
subplot(2,2,1);
plot(time(ind_upp2:ind_low),p_ret_90_17);
xlabel('Portfolio 1-7','fontsize',14);
ylabel('Return','fontsize',14);
legend(names_P,'fontsize',4,'Location','SouthEast');
title('Return for 1990/07-2017/10','fontsize',14);
subplot(2,2,2);
bar(p_ret_90_17_average);
%plot(names_P,p_ret_90_17_average);
set(gca,'XTickLabel',names_P)
xlabel('Portfolio 1-7','fontsize',14);
ylabel('Average return','fontsize',14);
title('Average return for 1990/01-2017/10','fontsize',14);
subplot(2,2,3);
bar(p_ret_90_17_total);
%plot(names_P,p_ret_90_17_total);
set(gca,'XTickLabel',names_P)
xlabel('Portfolio 1-7','fontsize',14);
ylabel('Total return','fontsize',14);
title('Total return for 1990/01-2017/10','fontsize',14);
subplot(2,2,4);
bar(p_ret_90_17_sharpe);
%plot(names_P,p_ret_90_17_sharpe);
set(gca,'XTickLabel',names_P)
xlabel('Portfolio 1-7','fontsize',14);
ylabel('Sharpe ratio','fontsize',14);
title('Sharpe ratio for 1990/01-2017/10','fontsize',14);

% Period January 2000 - October 2017
ind_upp3=find(data_B(:,1)==200001);   
p_ret_00_17=p_ret2((ind_upp3-60):(ind_low-60),:);
rf_00_17=mean(data_riskfree(ind_upp3:ind_low,2));
p_ret_00_17_average=mean(p_ret_00_17);
p_ret_00_17_total=sum(p_ret_00_17);
p_ret_00_17_sharpe=(p_ret_00_17_average-rf_00_17)./std(p_ret_00_17);

figure(11)
subplot(2,2,1);
plot(time(ind_upp3:ind_low),p_ret_00_17);
xlabel('Portfolio 1-7','fontsize',14);
ylabel('Return','fontsize',14);
legend(names_P,'fontsize',4,'Location','SouthEast');
title('Return for 2000/01-2017/10','fontsize',14);
subplot(2,2,2);
bar(p_ret_00_17_average);
%plot(names_P,p_ret_00_17_average);
set(gca,'XTickLabel',names_P)
xlabel('Portfolio 1-7','fontsize',14);
ylabel('Average return','fontsize',14);
title('Average return for 2000/01-2017/10','fontsize',14);
subplot(2,2,3);
bar(p_ret_00_17_total);
%plot(names_P,p_ret_00_17_total);
set(gca,'XTickLabel',names_P)
xlabel('Portfolio 1-7','fontsize',14);
ylabel('Total return','fontsize',14);
title('Total return for 2000/01-2017/10','fontsize',14);
subplot(2,2,4);
bar(p_ret_00_17_sharpe);
%plot(names_P,p_ret_00_17_sharpe);
set(gca,'XTickLabel',names_P)
xlabel('Portfolio 1-7','fontsize',14);
ylabel('Sharpe ratio','fontsize',14);
title('Sharpe ratio for 2000/01-2017/10','fontsize',14);