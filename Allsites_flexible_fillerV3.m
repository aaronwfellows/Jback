function Allsites_flexible_fillerV3(iSite)
%iSite=10;
%clear -iSite;
% load('C:\Analysis\California\Tower data\Sage_v3_2.mat');
% Startyear = 2006;
% Energyclose = 0.90;
% Good = DATA(:,1) >= 200;
% DATA = DATA(Good,:);
% 

% load('C:\Analysis\California\Tower data\Grass_v3_2.mat');
% Startyear = 2006;
% Energyclose = 0.90;
% 
% Good = DATA(:,1) >= 200;
% DATA = DATA(Good,:);

% 
% load('C:\Analysis\California\Tower data\JamesRes_v3_2.mat');
% Startyear = 2006;
% Energyclose = 0.73;
% 
% load('C:\Analysis\California\Tower data\Pinyon_v3_2.mat');
% Startyear = 2006;
% Energyclose = 0.93;
% 
% load('C:\Analysis\California\Tower data\PinyonBurn_v3_2.mat');
% Startyear = 2006;
% Energyclose = 0.75;
% 
% load('C:\Analysis\California\Tower data\LowDesert_v3_2.mat');
% Startyear = 2006;
% Energyclose = 0.67; 
% 
% load('C:\towerData\WebsiteData\SJER_v3_2.mat');
% Startyear = 2010;
% Energyclose = 0.80;
% 
if iSite == 10
    load('C:\towerData\WebsiteData\Soaproot_v3_2.mat'); %#ok<*MCABF>
    Energyclose = 0.73;
elseif iSite == 7
    load('C:\towerData\WebsiteData\P301_v3_2.mat');
    Energyclose = 0.8585;
elseif iSite == 8
    load('C:\towerData\WebsiteData\SJER_v3_2.mat');
    Energyclose = 0.7763;
elseif iSite == 5
    load('C:\towerData\WebsiteData\Sage_v3_2.mat');
    Energyclose = 0.8516;
elseif iSite == 4
    load('C:\towerData\WebsiteData\Grass_v3_2.mat');
    Energyclose = 0.8940;
else
    disp('No energy budget closure term for this site yet');
    return
end
% 
% load('C:\towerData\WebsiteData\Shorthair_v3_2.mat');
% Startyear = 2010;
% Energyclose = 0.78;  %Estimated



Ustarfilter = 0.2;
dimthresh = 200;
brightthresh = 400;



%Creat continous meteorology fields first
%use two filling approaches
%for short missing periods (<=4 missing intervals) just linear extrapolate
%for longer missing intervals fill in the mean of 25 day blocks

%Version 2 - In general, a simple (conservative) approach)
% R calculated by extrapolating u*>0.3 periods to darkness using nonlinear
% fit.  Includes cleaning to remove outliers and uses regression info from
% previous period if results are unreasonable (really small or large R)
% 
% GEE calcuted from same nonlinear fit, based on sobtracting R from NEE or 
% calculating the GEE fit by dropping outthe first (R) term
%
% NEE calculated similar to GEE
%
% Et claculated from linear regression between Et and K
% Et set to zero at night - this was needed since regrsssion yielded some
% high night Ets in summer





%create filled T first
%start by calculating length of missing period
Tconsensus = DATA(:,5);
bad = isnan(Tconsensus);
Tconsensus(bad) = DATA(bad,6);
%bad = isnan(Tconsensus);
Tunfill = Tconsensus;
Tfill = Tconsensus;
missing = isnan(Tunfill);
counter = 0; 

cummissing = NaN * missing;
for i = 1:length(missing)
    if missing(i) == 0;
       counter = 0 ;
    else
       counter = counter + 1;
    end
    cummissing(i) = counter;
end

%now create liner extrapolate
Tinterp = interp1(DATA(~missing,1),Tunfill(~missing),DATA(:,1),'linear');

%now creat series of mean diel curves for 25-day chnunks
binwid = 25;
start = min(DATA(:,1));
binnum = floor((max(DATA(:,1)- start)/binwid));
bincenter = (0:(binnum-1))* binwid + start;
Toutbin = zeros(binnum,48);
Hour = 24 * (DATA(:,1) - floor(DATA(:,1)));
good = isfinite(Tunfill);

for i = 1:binnum;
    focusday = (DATA(:,1) > (bincenter(i) - binwid/2)) & (DATA(:,1) <= (bincenter(i) + binwid/2));
    for j = 1:48;
        focushour = Hour > j/2 - 0.75 & Hour < j/2 - 0.25;
        Toutbin(i,j) = mean(Tunfill(good&focusday&focushour));
    end
end
Toutbin = [bincenter' Toutbin];

%Now do the filling
for i = 1:length(missing);
    if missing(i);
        if cummissing(i) < 5;
            Tfill(i) = Tinterp(i);
        else
            rownumb = floor((max(DATA(i,1)- start+binwid/2)/binwid)) + 1;
            colnumb = round(Hour(i)*2) + 2;
            if rownumb > binnum;
                rownumb = binnum;
            end
            Tfill(i) = Toutbin(rownumb,colnumb);
        end
    end
end


%Now create filled light
%start by calculating length of missing period
Kunfill = DATA(:,13);
missing = isnan(Kunfill);
Kunfill(missing) = DATA(missing,11).*0.52647;
Kfill = Kunfill;
missing = isnan(Kunfill);
counter = 0; 

for i = 1:length(missing);
    if missing(i) == 0;
       counter = 0 ;
    else
       counter = counter + 1;
    end
    cummissing(i) = counter;
end

%now create liner extrapolate
Kinterp = interp1(DATA(~missing,1),Kunfill(~missing),DATA(:,1),'linear');

%now creat series of mean diel curves for 25-day chnunks
binwid = 25;
start = min(DATA(:,1));
binnum = floor((max(DATA(:,1)- start)/binwid));
bincenter = (0:(binnum-1))* binwid + start;
Koutbin = zeros(binnum,48);
Hour = 24 * (DATA(:,1) - floor(DATA(:,1)));
good = isfinite(Kunfill);

for i = 1:binnum;
    focusday = (DATA(:,1) > (bincenter(i) - binwid/2)) & (DATA(:,1) <= (bincenter(i) + binwid/2));
    for j = 1:48;
        focushour = Hour > j/2 - 0.75 & Hour < j/2 - 0.25;
        Koutbin(i,j) = mean(Kunfill(good&focusday&focushour));
    end
end
Koutbin = [bincenter' Koutbin];

%Now do the filling
for i = 1:length(missing);
    if missing(i);
        if cummissing(i) < 5;
            Kfill(i) = Kinterp(i);
        else
            rownumb = floor((max(DATA(i,1)- start+binwid/2)/binwid)) + 1;
            colnumb = round(Hour(i)*2) + 2;
            if rownumb > binnum;
                rownumb = binnum;
            end
            Kfill(i) = Koutbin(rownumb,colnumb);
        end
    end
end



%Now create filled Rnet
%start by calculating length of missing period
Rnetunfill = DATA(:,15);
missing = isnan(Rnetunfill);
Rnetunfill(missing) = Kfill(missing)*0.70582-43.776;
Rnetfill = Rnetunfill;


% %Now do CO2 flux filling 
% %%% Now calculate Resp using light curves
NEEunfill = DATA(:,19);
calm = DATA(:,2) < Ustarfilter & isfinite(DATA(:,2));
NEEunfill(calm) = NaN;
bad = DATA(:,19) > 10 | DATA(:,19) < -35;
NEEunfill(bad) = NaN;

NEEfit = DATA(:,19) * NaN;
Respfit = DATA(:,19) * NaN;
%GEEfit = DATA(:,19) * NaN;

%now creat series of light curves for 10-day chnunks and calculate the flux
%at each irradaince

binwid = 10;
start = binwid/2+min(DATA(:,1));
binnum = floor((max(DATA(:,1)- start)/binwid));
bincenter = (0:(binnum-1))* binwid + start + binwid/2;
% Toutbin = zeros(binnum,48);
% Hour = 24 * (DATA(:,1) - floor(DATA(:,1)));
% good = isfinite(Kunfill);
GEEfit = NEEunfill * NaN;

parfitdimold = [-0.001  3];
parfitmediumold = [-0.001  3];
parfitbrightold = [-0.001  3];

for i = 1:binnum;
    focusday = (DATA(:,1) > (bincenter(i) - 5)) & (DATA(:,1) <= (bincenter(i) + 5));
%     par_day = Kfill(focusday & isfinite(NEEunfill));
%     nee_day = NEEunfill(focusday & isfinite(NEEunfill));
  
    dim = Kfill < dimthresh;
    medium = Kfill >= dimthresh & Kfill < brightthresh;
    bright = Kfill >= brightthresh;

    par_day = Kfill(focusday & isfinite(NEEunfill) & dim);
    nee_day = NEEunfill(focusday & isfinite(NEEunfill) & dim);

    if length(nee_day) > 5
        parfitdim = polyfit(par_day, nee_day, 1);
        parfitdimold = parfitdim;
    else 
        parfitdim = parfitdimold;   
    end
        
    par_day = Kfill(focusday & isfinite(NEEunfill) & medium);
    nee_day = NEEunfill(focusday & isfinite(NEEunfill) & medium);

    if length(nee_day) > 5
        parfitmedium = polyfit(par_day, nee_day, 1);
        parfitmediumold = parfitmedium;
    else 
        parfitmedium = parfitmediumold;   
    end
    
    par_day = Kfill(focusday & isfinite(NEEunfill) & bright);
    nee_day = NEEunfill(focusday & isfinite(NEEunfill) & bright);

    if length(nee_day) > 5
        parfitbright = polyfit(par_day, nee_day, 1);
        parfitbrightold = parfitbright;
    else 
        parfitbright = parfitmediumold;   
    end
    
    Respfit(focusday) = parfitdim(2);        
    NEEfit(focusday & dim) = parfitdim(2) + (parfitdim(1)*Kfill(focusday & dim));
    NEEfit(focusday & medium) = parfitmedium(2) + (parfitmedium(1)*Kfill(focusday & medium));
    NEEfit(focusday & bright) = parfitbright(2) + (parfitbright(1)*Kfill(focusday & bright));
    NEEfit(NEEfit < -40) = NaN;
    GEEfit(focusday) =  NEEfit(focusday) - Respfit(focusday);
    par_day = Kfill(focusday);

    %     plot(par_day,nee_day,'.')
    %     i
    %     b
    %     pause
    
end

% 
% 
% %Now calculate annual respiration
missing = isnan(Respfit);
Respfit(missing) = 0;
Respfill = Respfit;

% %%% Now calculate GEE
GEEunfill = DATA(:,19) - Respfill;
calm = DATA(:,2) < Ustarfilter & isfinite(DATA(:,2));
GEEunfill(calm) = NaN;
light = Kfill > 2;
GEEunfill(~light) = 0;
missing = isnan(GEEunfill);
GEEfill = GEEunfill;
GEEfill(missing) = GEEfit(missing);
missing = isnan(GEEfill);
GEEfill(missing) = 0;

% %Finally calculate filled NEE and NEP
NEEunfill = DATA(:,19);
calm = DATA(:,2) < Ustarfilter & isfinite(DATA(:,2));
NEEunfill(calm) = NaN;
NEEfill = NEEunfill;
missing = isnan(NEEfill);
NEEfill(missing) = NEEfit(missing);
missing = isnan(NEEfill);
NEEfill(missing) = 0;


%%%%%%
%%% Now do water and energy
%%%%%%
%% Now calculate E
Eunfill = DATA(:,20);
calm = DATA(:,2) < Ustarfilter & isfinite(DATA(:,2));
Eunfill(calm) = NaN;
%Efit = Eunfill * NaN;

%now creat series of light curves for 10-day chnunks and calculate the flux
%at each irradaince

binwid = 10;
start = binwid/2+min(DATA(:,1));
binnum = floor((max(DATA(:,1)- start)/binwid));
bincenter = (0:(binnum-1))* binwid + start + binwid/2;
% Toutbin = zeros(binnum,48);
% Hour = 24 * (DATA(:,1) - floor(DATA(:,1)));
% good = isfinite(Kunfill);
Efit = Eunfill * NaN;
parfitold = [0.3 0.006];

countthis=0;

for i = 1:binnum;
    focusday = (DATA(:,1) > (bincenter(i) - 5)) & (DATA(:,1) <= (bincenter(i) + 5));
    par_day = Kfill(focusday & isfinite(Eunfill));
    et_day = Eunfill(focusday & isfinite(Eunfill));
    
    X = [ones(size(par_day)) par_day];
    y = et_day;
    
    if length(et_day) > 20
        b = regress(y,X);
    else
        b = parfitold;
        countthis = countthis+1;
    end
    parfitold = b;
    
    Efit(focusday) = b(1) + b(2)*Kfill(focusday);
    
    %     plot(par_day,nee_day,'.')
    %     i
    %     b
    %     pause
    
end

%Now do the E filling
missing = isnan(Eunfill);
Efill = Eunfill;
Efill(missing) = Efit(missing);
missing = isnan(Efill);
Efill(missing) = 0;


%%%%%%
%%%%%%
%% Now calculate H
Hunfill = DATA(:,17);

calm = DATA(:,2) < Ustarfilter & isfinite(DATA(:,2));
Hunfill(calm) = NaN;

%Hfit = Hunfill * NaN;

%now creat series of light curves for 10-day chnunks and calculate the flux
%at each irradaince

binwid = 30;
start = binwid/2+min(DATA(:,1));
binnum = floor((max(DATA(:,1)- start)/binwid));
bincenter = (0:(binnum-1))* binwid + start + binwid/2;
Toutbin = zeros(binnum,48);
Hour = 24 * (DATA(:,1) - floor(DATA(:,1)));
good = isfinite(Kunfill);
Hfit = Hunfill * NaN;
parfitold = [-20 0.6];

for i = 1:binnum;
    focusday = (DATA(:,1) > (bincenter(i) - binwid/2)) & (DATA(:,1) <= (bincenter(i) + binwid/2));
    par_day = Kfill(focusday & isfinite(Hunfill));
    nee_day = Hunfill(focusday & isfinite(Hunfill));
    
    X = [ones(size(par_day)) par_day];
    y = nee_day;
    
    if length(nee_day) > 20
        b = regress(y,X);
    else
        b = parfitold;
    end
    parfitold = b;
    
    Hfit(focusday) = b(1) + b(2)*Kfill(focusday);
    
%         plot(par_day,nee_day,'.')
%         i
%         b
%         pause
    
end

%Now do the H filling
missing = isnan(Hunfill);
Hfill = Hunfill;
Hfill(missing) = Hfit(missing);
missing = isnan(Hfill);
Hfill(missing) = 0;

%%%%%%%%
%%%%%%%%
%%%%%%%%
%Now calculate calculate all of the sums and annual fluxes
% cumH = cumsum(Hfill)*(1800)/(1000*1000)./Energyclose;%0.90 accoutns for energy budget closure
% cumE = cumsum(Efill)*(18.02*1800)/(1000*1000)./Energyclose;%*0.90 accoutns for energy budget closure
% cumNEE = cumsum(NEEfill)/(2.31*2)./Energyclose;%*0.90 accoutns for energy budget closure
% cumGEE = cumsum(GEEfill)/(2.31*2)./Energyclose;%*0.90 accoutns for energy budget closure
% cumResp = cumsum(Respfill)/(2.31*2)./Energyclose;%*0.90 accoutns for energy budget closure

%Now calculate annual sums
Resp = [];
GEE  = [];
NEE  = [];
E    = [];
H    = [];
DATA(:,1) = DATA(:,1)+datenum(2006,1,0);
for yr = 2007:2013
    yii = DATA(:,1) > datenum(yr-1,10,1)-0.01 & DATA(:,1) < datenum(yr,10,1)+0.01;
    if sum(yii) > 0
        Resp = [Resp sum(Respfill(yii))]; %#ok<*AGROW>
        GEE  = [GEE  sum(GEEfill(yii))]; 
        NEE  = [NEE  sum(NEEfill(yii))];
        E    = [E    sum(Efill(yii))];
        H    = [H    sum(Hfill(yii))];
    else
        Resp = [Resp 0];
        GEE  = [GEE  0];
        NEE  = [NEE  0];
        E    = [E 0];
        H    = [H 0];        
    end
end

H    = H*1800*1e-6./Energyclose;
E    = E*1800*18.02*1e-6./Energyclose;
NEE  = NEE/(2.31*2)./Energyclose;
GEE  = GEE/(2.31*2)./Energyclose;
Resp = Resp/(2.31*2)./Energyclose;

% cumH = cumsum(Hfill)*(1800)/(1000*1000)./Energyclose;%0.90 accoutns for energy budget closure
% cumE = cumsum(Efill)*(18.02*1800)/(1000*1000)./Energyclose;%*0.90 accoutns for energy budget closure
% cumNEE = cumsum(NEEfill)/(2.31*2)./Energyclose;%*0.90 accoutns for energy budget closure
% cumGEE = cumsum(GEEfill)/(2.31*2)./Energyclose;%*0.90 accoutns for energy budget closure
% cumResp = cumsum(Respfill)/(2.31*2)./Energyclose;%*0.90 accoutns for energy budget closure

    
AnnC = [Resp' GEE' NEE'];
AnnFluxes = [Resp' GEE' NEE' E' H']


Filledup = [DATA(:,1) Kfill Rnetfill Tfill NEEfill GEEfill Respfill Efill Hfill];
% save('C:\Analysis\California\Tower data\Sage_Filled.mat','Filledup');
% save('C:\Analysis\California\Tower data\Grass_Filled.mat','Filledup');
% save('C:\Analysis\California\Tower data\Pinyon_Filled.mat','Filledup');
% save('C:\Analysis\California\Tower data\PFburn_Filled.mat','Filledup');
% save('C:\Analysis\California\Tower data\Desert_Filled.mat','Filledup');
% save('C:\Analysis\California\Tower data\James_Filled.mat','Filledup');
%save('C:\towerData\filled\SJER_Filled.mat','Filledup');
if iSite == 10
    save('C:\towerData\filled\Soaproot_Filled.mat','Filledup');
elseif iSite == 7
    save('C:\towerData\filled\P301_Filled.mat','Filledup');
elseif iSite == 8
    save('C:\towerData\filled\SJER_Filled.mat','Filledup');
elseif iSite == 5
    save('C:\towerData\filled\Sage_Filled.mat','Filledup');
elseif iSite == 4
    save('C:\towerData\filled\Grass_Filled.mat','Filledup');
else
    disp('did not save');
end
%save('C:\towerData\filled\Shorthair_Filled.mat','Filledup');




% 
% %Summing
% binwid = 10;
% start = binwid/2+min(DATA(:,1));
% binnum = floor((max(DATA(:,1)- start)/binwid));
% bincenter = (0:(binnum-1))* binwid + start + binwid/2;
% Sageoutbin = zeros(binnum,6);
% Sageoutbin(:,1) = bincenter;
% for i = 1:binnum;
%     focusday = (DATA(:,1) > (bincenter(i) - 0.5)) & (DATA(:,1) <= (bincenter(i) + 0.5));
%     Sageoutbin(i,2) = sum(FilledSage(focusday,5))./(2.31*2)./0.90;%*0.90 accoutns for energy budget closure
%     Sageoutbin(i,3) = sum(FilledSage(focusday,6))./(2.31*2)./0.90;%*0.90 accoutns for energy budget closure
%     Sageoutbin(i,4) = sum(FilledSage(focusday,7))./(2.31*2)./0.90;%*0.90 accoutns for energy budget closure
%     Sageoutbin(i,5) = sum(FilledSage(focusday,8))*(18.02*1800)/(1000*1000)./0.90;%*0.90 accoutns for energy budget closure
%     Sageoutbin(i,6) = sum(FilledSage(focusday,9))*(1800)/(1000*1000)./0.90./0.90;%*0.90 accoutns for energy budget closure
% end
% 
% 

return;
