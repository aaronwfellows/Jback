function [m,b,Sy,Sm,Sb,r2,YPRED,CONFINT,SE_slope,SE_int,tcrit,CI_2,Probab ] = myreg(X,Y)
% function [m,b,Sy,Sm,Sb,r2,YPRED,CONFINT,SE_slope,SE_int,tcrit ,CI_2 ,Probab]=MYREG(X,Y);

%had to change the name of this function since the former name of this
%(reg.m) is already used by matlab v6 and v7

% New Output: CI_2
% 
% This is a nx2 matrix with the standard error for the predicted Yi for a given value of Xi
% Formula from Box 14.2 p 472 Sokal and Rolf




%A little routine to calculate regression statistics on 2 arrays %
%Here are the values we calculate
%
%	1. Slope, m
%	2. intercept, b
%	3. std deviation of the residuals, Sy
%	4. std deviation of the slope, Sm
%	5. std deviation of the intercept, Sb
%	6. std deviation of analytical error, Sc (not included here)
%	7. r2,
%


%   These statistics come from "Principles of Instrumental Analysis", 4th Ed
%   by Skoog and Leary (1992). Saunders College Publishing PP A17 - A18.


% Standard Errors for regression coefficients are based on Moore and
% McCabe (1998) "Introduction to the Practice of Statistics" 3rd Ed. WH
% Freeman & Co, NY.Chapter 10



%make sure dimensions agree
% try
    if lnn(Y)~=0 || lnn(X)~=0,
    if size(Y,2)>size(Y,1),Y=Y';end
    if size(X,2)>size(X,1),Y=Y';end
    end
% catch
%     disp('L45');
% end

if lnn(X)==0,
    disp('The X series is totally full of NaNs, cannot calculate regression stats');
    %disp('Programed paused');
    %pause
    m=NaN;b=NaN;Sy=NaN;Sm=NaN;Sb=NaN;r2=NaN;CONFINT=NaN;YPRED=NaN*ones(1,length(X));
    SE_slope=NaN;SE_int=NaN;tcrit=NaN;CI_2=NaN;Probab=NaN;
    return
end

if lnn(Y)==0,
    disp('The Y series is totally full of NaNs, cannot calculate regression stats');
    %disp('Program paused');
    %pause
    m=NaN;b=NaN;Sy=NaN;Sm=NaN;Sb=NaN;r2=NaN;CONFINT=NaN;YPRED=NaN*ones(1,length(X));
    SE_slope=NaN;SE_int=NaN;tcrit=NaN;CI_2=NaN;Probab=NaN;
    return

end


Xold=X;
Yold=Y;
%get rid of NaNs
[X,Y]=eatNaN(X,Y);

%calculate preliminary statistics

meanx=mean(X);
meany=mean(Y);
N=length(X);
%NY=length(Y);%for diagnostics only
%sum of squares
devsX=meanx-X;
devsY=meany-Y;
if size(devsY)~=size(devsX),
    devsY=devsY';
end


Sxx=sum((devsX).^2);
Syy=sum((devsY).^2);
Sxy=sum((devsX).*(devsY));




% 1. Slope

m=Sxy/Sxx;

%	2. intercept, b

b=meany-m*meanx;

%	3. std deviation of the residuals, Sy

Sy=sqrt((Syy-m^2*Sxx)/(N-2));

%	4. std deviation of the slope, Sm

Sm=Sy/sqrt(Sxx);

%	5. std deviation of the intercept, Sb

sigma_x2=sum(X.^2);

sigmax_2=(sum(X))^2;

Sb=Sy.*sqrt(1/(N - sigmax_2/sigma_x2));


%6. correlation coeffiecient

r2=rsq(Y,X);

%   more features added on 12/5/04. Confidence Intervals on a regression
%   line.
%
%   From Moore and McCabe (1998) Introduction to the Practice of
%   Statistics. 3rd Ed. Freeman and Co. New York. pp 671-691

%The formula for a 95 percent CI for a mean response to a dependent
%variable is
%
%  95 CI   =   tcrit  * SE
%
%    where SE = s *  sqrt [ 1 + 1/n + (X - XBAR )^2 / SS(xi - xbar ) ]

% and tcrit is the value for the t(n-2) density curve C between -tcrit and
% tcrit
%
% and  s = sqrt ( SS(ydevs) / (n-2) )
%
%
%
%
%




Ypred = Xold * m +b;
Yres = Yold - Ypred;
SS_Yres = nansum ( Yres.^2 );

var_Yres = SS_Yres /(N-2);
sd_Yres = sqrt(var_Yres);

CI_level = 95; %express as percernt

alp = 0.5 * ( 1 - CI_level/100 );
df = N - 2;

tcrit = abs ( tinv ( alp , df ) );

SE_y = sd_Yres * sqrt ( 1 + 1/N + (Xold - meanx ).^2/Sxx); %THIS IS FOR A
% PREDICTING A FUTURE VALUE 
    SE_slope = sd_Yres/sqrt(Sxx);
SE_int = sd_Yres*sqrt( (1/N) + meanx^2/sqrt(Sxx));

ConfInt = tcrit * SE_y;%This is the CI for the mean response

%CI_upper = Ypred + ConfInt;
%CI_lower = Ypred - ConfInt;

YPRED = Ypred;
CONFINT = ConfInt;

CI_given_upper = Ypred + SE_y;
CI_given_lower = Ypred - SE_y;

CI_2 = [CI_given_upper CI_given_lower];

%workout P
%[bXX,bintXX,rXX,rintXX,statsXX] = my_regress(Y,X);
[~,~,~,~,statsXX] = my_regress(Y,X);
Probab = statsXX(3);


return


function [ACLEAN, BCLEAN]=eatNaN(A,B)

%make sure arrays are horizontal
if size(A,1)>size(A,2),A=A';end
if size(B,1)>size(B,2),B=B';end

%if arrays are different sizes trim one to the other
if length(A)>length(B),
    A=A(1:length(B));
elseif length(B)>length(A),
    B=B(1:length(A));
end






if lnn(A)==0,
    ACLEAN(1:length(A))=NaN;
    BCLEAN(1:length(A))=NaN;
    return
elseif lnn(B)==0,
    ACLEAN(1:length(A))=NaN;
    BCLEAN(1:length(A))=NaN;
    return
else
    % disp(['~isnan A is ' num2str(length(~isnan(A))) ' long']);
    % disp([' isnan B is ' num2str(length(~isnan(B))) ' long']);


    A1=A(~isnan(A));
    B1=B(~isnan(A));

    A2=A1(~isnan(B1));
    B2=B1(~isnan(B1));

end

ACLEAN=A2;
BCLEAN=B2;

return

function [len_isnotanan]=lnn(X)

len_isnotanan=length(find(~isnan(X)));

return