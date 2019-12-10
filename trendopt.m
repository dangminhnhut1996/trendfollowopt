%% 1. The Trendfollowing (Optimize with hypothesis test) 

% Note 1: The section with number denote important section.

% Note 2: section 1 and 2 do not have semi colon at the end of the
% function as intended to let user see which lower bound(map to which 
% upper eaiser

% Note 3: Uncomment code in below section to run the test

% Note 4: If uncomment the code,
% Please only use "Run Section" to check for each section as it will
% take a long time to run all these code. Feel free to use "Run" if time is
% not a matter.


clc; clear;
Data=readtable('NASDAQ1985-2019.xlsx');
  
[mtb,mkb]=learnhp(Data,4000,length(Data.Date),500)

%% 2. The Trendfollowing (Optimize without hypothesis test) 

[mtbw,mkbw]=learn(Data,4000,length(Data.Date),500)
%% Test Hypothesis test 

%[~,~,a,b,~]=trendfollow(Data,40,60);
%[~,~,a2,b2,~]=trendfollow(Data,39,64);
%[rtest,k,k2,h1,h2]=hptest(a,a2,b,b2);
%% 3. Test trend following
Data=readtable('NASDAQ1985-2019.xlsx');
[totalprofit,totalpreturn,a,b,u1]=trendfollow(Data,40,60);
[totalprofit2,totalpreturn2,a2,b2,u2]=trendfollow(Data,39,64);
%avprofit=mean(a);
%avpretrun=mean(b)*100;
%avprofit2=mean(a2);
%avpretrun2=mean(b2)*100;

%% Test Optimum with hypothesis test

%[tb,kb]=opthp(Data);
%% Test Optimum without Hypothesis test

%[tbw,kbw]=opt(Data);

%% Graph for RSI and EMA(RSI)

%{
Data=readtable('NASDAQ1985-2019.xlsx');
index = rsindex(Data); %RSI calculation
window_size = 14;
inddata=index(15:end,1);
exp = movavg(inddata,'exponential',window_size); %exponential moving average
date=Data.Date(15:end);
subplot(2,1,1);
plot(Data.Date,index.RelativeStrengthIndex)
title('Relative Strength Index for Nasdaq')
subplot(2,1,2);
plot(date,inddata.RelativeStrengthIndex,date,exp.RelativeStrengthIndex)
title('Relative Strength Index for Nasdaq exponential')
%}
%% Hypothesis Test

function [rtest,k,k2,h1,h2]=hptest(a,ar,b,br)
rtest=0;  
k=vartest2(a,ar,'Alpha',0.1); %F test 
k2=vartest2(b,br,'Alpha',0.1); %F test
if k==1
    [h1,~] = ttest2(a,ar,'tail', 'right','Alpha',0.1,'Vartype','unequal'); %T test with unequal variance
 
else
    [h1,~] = ttest2(a,ar,'tail', 'right','Alpha',0.1);%T test with equal variance
 
end
if k2==1
    [h2,~] = ttest2(b,br,'tail', 'right','Alpha',0.1,'Vartype','unequal');%T test with unequal variance
else
    [h2,~] = ttest2(b,br,'tail', 'right','Alpha',0.1);%T test with equal variance
end
    if h1==1 && h2==1
    rtest=1;
    end
end



%% optimum with hypothesis test

function [tb,kb]=opthp(data) 

tb=40;
kb=60;
[prob,preb,a,b,~]=trendfollow(data,40,60);
for t=30:40 
    for k=60:70
        [pro,pre,a2,b2,~]=trendfollow(data,t,k);
        if pro>prob % is totalprofit better
            if pre>preb % is total%return better
            [hyptest,~,~,~,~]=hptest(a,a2,b,b2); %hypothesis test( above section)
            if hyptest==1
            a=a2;
            b=b2;
            prob=pro;
            preb=pre;
            tb=t;
            kb=k;
            end
        end
    end
end
end
end
%% learn with hypothesis test
function [mtb,mkb]=learnhp(Data,starts,ends,step) 
g=starts;
ntb=1;
nkb=1;
te=0;


if mod((ends-starts), step)==0
    te=1;
end

while g~=ends || g==ends && te==1  % taking into account the case ((ends-starts)/step) mod step =0
    if g<ends
        x=Data(1:g,:);
        [tb,kb]=opthp(x);
        mtb(ntb)=tb;
        mkb(nkb)=kb;
        ntb=nkb+1;
        nkb=nkb+1;
        g=g+step; 
    else    
    x=Data(1:ends,:);
    [tb,kb]=opthp(x);
    mtb(ntb)=tb;
    mkb(nkb)=kb;
    g=ends;
    te=0;
    end
end

 
end   
%% Optimum without hypothesis test

function [tb,kb]=opt(data) 
prob=0;
preb=0;
tb=0;
kb=0;
for t=30:40 
    for k=60:70
        [pro,pre]=trendfollow(data,t,k);
        if pro>prob
            if pre>preb
            prob=pro;
            preb=pre;
            tb=t;
            kb=k;
        end
    end
end
end
end
%% learn without hypothesis test
function [mtb,mkb]=learn(Data,starts,ends,step) 
g=starts; 
ntb=1; % initiate index number in lower bound array
nkb=1; % initiate index number in upper bound array
te=0; 


if mod((ends-starts), step)==0
    te=1;
end

while g~=ends || g==ends && te==1 % taking into account the case ((ends-starts)/step) mod step =0
    if g<ends
        x=Data(1:g,:);
        [tb,kb]=opt(x);
        mtb(ntb)=tb;
        mkb(nkb)=kb;
        ntb=nkb+1; 
        nkb=nkb+1; 
        g=g+step; 
    else    
    x=Data(1:ends,:);
    [tb,kb]=opt(x);
    mtb(ntb)=tb;
    mkb(nkb)=kb;
    g=ends;
    te=0;
    end
end

 
end   

%% Trendfollowing algorithm
function [totalprofit,totalpreturn,profit,preturn,bs]=trendfollow(Data,lowerbound,upperbound)
index = rsindex(Data); %RSI calculation
window_size = 14;
inddata=index(15:end,1);
exp = movavg(inddata,'exponential',window_size); %exponential moving average
date=Data.Date(15:end);
bs=zeros(length(date),2);
os=lowerbound;
ob=upperbound;
p=0;
q=0;
b=1; % initiate index number in buy array
s=1; % initiate index number in sell array
ab=0; %already buy
p1=1; %start long positon
p2=1; %close long position
q1=1; %start short positon
q2=1; %close short position
P=zeros(length(date),2);
Q=zeros(length(date),2);
for i=1:length(date)

if sign(Data.Day_DayBefore(i))== 1
    
    if inddata.RelativeStrengthIndex(i)> exp.RelativeStrengthIndex(i) && (exp.RelativeStrengthIndex(i)<os||exp.RelativeStrengthIndex(i)>ob)
     
     if p==0 && q==0 && ab==0   %    If no position opened

         p=1;                   %       open  a long position, P'
         P(p1,1)=i;
         p1=p1+1;
         
     else if q==1 && ab==1     %    else if short position opened
             q=0;               %       Close out short position, Q'
             Q(q2,1)=i;
             q2=q2+1;
             bs(s,2)=i;         %keep track sell
             s=s+1;             
             ab=0;              % reset already buy
         end
     end
    end


    

else if sign(Data.Day_DayBefore(i))== -1
        if inddata.RelativeStrengthIndex(i)< exp.RelativeStrengthIndex(i) && (exp.RelativeStrengthIndex(i)<os||exp.RelativeStrengthIndex(i)>ob)
            if p==0 && q==0 && ab==1    %    If no position opened
                q=1;                    %    open  a long position, Q'
                Q(q1,1)=i;
                q1=q1+1;
         
            else if p==1                %    else if long position opened
             p=0;                       %       Close out short position, P'
             P(p2,1)=i;
             p2=p2+1;
             bs(b,1)=i;                 % keep track buying
             b=b+1;                     
             ab=1;                      % turn already buy on
                end
            end
        end
    end


end
    if mod(i,390)==0             % If end of market (390 days)
    p=0;   q=0;                 % Close all opened position 
    end
end

j=1;
while bs(j,1)~=0 && bs(j,2)~=0
    profit(j)=Data.Close(bs(j,2)+14)-Data.Close(bs(j,1)+14); %sell price-buy price
    preturn(j)=(Data.Close(bs(j,2)+14)-Data.Close(bs(j,1)+14))/Data.Close(bs(j,1)+14);%(sell price-buy price)/buy price
    j=j+1;
    
end
totalprofit=sum(profit);
totalpreturn=sum(preturn)*100;
end

