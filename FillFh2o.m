function [DATA, Efill_regress, Efill_robust, Kfill] = FillFh2o(DATA, nfill)
%Fill the EC data

%% Now create filled light
% start by calculating length of missing period
Kunfill = DATA(:,13);
missing = isnan(Kunfill);
Kfill   = Kunfill;

counter = 0;

for i = 1:length(missing);
    if missing(i) == 0
       counter = 0 ;
    else
       counter = counter + 1;
    end
    cummissing(i) = counter;
end

%% now create liner extrapolate
Kinterp = interp1(DATA(~missing,1),Kunfill(~missing),DATA(:,1),'linear');

%% now create series of mean diel curves for 25-day chnunks
binwid = nfill;
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

%% Now do the filling
for i = 1:length(missing);
    if missing(i);
        if cummissing(i) < 5;
            Kfill(i) = Kinterp(i);
        else
            rownumb = floor((max(DATA(i,1)- start+12.5)/25)) + 1;
            colnumb = round(Hour(i)*2) + 2;
            if rownumb > binnum;
                rownumb = binnum;
            end
            Kfill(i) = Koutbin(rownumb,colnumb);
        end
    end
end

%%

%%% Now do water and energy
%%%%%%
%%% Now calculate E
Eunfill = DATA(:,20);
Efit_regress = Eunfill * NaN;
Efit_robust = Eunfill * NaN;

%now creat series of light curves and calculate the flux
%at each irradaince

binwid = nfill;
start = binwid/2+min(DATA(:,1));
binnum = floor((max(DATA(:,1)- start)/binwid));
bincenter = (0:(binnum-1))* binwid + start + binwid/2;
Toutbin = zeros(binnum,48);
Hour = 24 * (DATA(:,1) - floor(DATA(:,1)));
good = isfinite(Kunfill); 
parfitold = [0.8 0.006];

for i = 1:binnum;
    focusday = DATA(:,1) > (bincenter(i) - (binwid/2)) & DATA(:,1) <= (bincenter(i) + (binwid/2));
    par_day = Kfill(focusday & isfinite(Eunfill));
    et_day = Eunfill(focusday & isfinite(Eunfill));
    
    X = [ones(size(par_day)) par_day];
    y = et_day;
    
    if length(et_day) > 10
        b1 = regress(y,X);
        b2 = robustfit(X(:,2), y);        
    else 
        b1 = parfitold;
        b2 = parfitold;
    end
  
    %change parfitold to most recent period
    parfitold = b1;

    
    Efit_regress(focusday) = b1(1) + b1(2)*Kfill(focusday);
    Efit_robust(focusday) = b2(1) + b2(2)*Kfill(focusday);
         %plot(par_day,et_day,'b.') 
         %hold on
         %plot(par_day,b1(1) + b1(2)*par_day,'r.') 
         %hold on
         %plot(par_day,b2(1) + b2(2)*par_day,'k.') 
         %i
         %b1
         %b2
    
end

%Now do the E filling
missing = isnan(Eunfill);
Efill_regress = Eunfill;
Efill_regress(missing) = Efit_regress(missing);
missing = isnan(Efill_regress);
Efill_regress(missing) = 0;


missing = isnan(Eunfill);
Efill_robust = Eunfill;
Efill_robust(missing) = Efit_robust(missing);
missing = isnan(Efill_robust);
Efill_robust(missing) = 0;

