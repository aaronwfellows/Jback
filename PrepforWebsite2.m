function PrepforWebsite2(site)
%Formats cleaned data for Website

%This program 
%removes the leaf wetness and fuel moisture sensors 
%averages up the soil temperature sensors 
%renames the header
clear -site


out_root='C:\towerData\WebsiteData\';
%---------------------------------------------- -------------------------%
%grab the site data that interests you
%-----------------------------------------------------------------------%
if site == 1
    load('C:\towerData\combined\DC_Burn_cleaned.mat'); %This is the best data
    site_name = 'PinyonBurn'
    AddK =1;
    site_type = 1
    %get the soil heat storage
    load('C:\towerData\EnergyBudgetAnalysis\DC_Burn\SoilHeatStorage_D_BurnSoilHeat_2.mat')
elseif site == 2
    load('C:\towerData\combined\DC_LowDes_cleaned.mat')
    site_name = 'LowDesert'
    AddK =0;
    site_type = 1
    %get the soil heat storage
    load('C:\towerData\EnergyBudgetAnalysis\DC_LowDes\SoilHeatStorage_D_LowDesSoilHeat_2.mat')
elseif site == 3    
    load('C:\towerData\combined\DC_Pinyon_cleaned.mat')
    site_name = 'Pinyon'
    AddK =1;
    site_type = 1
    %get the soil heat storage
    load('C:\towerData\EnergyBudgetAnalysis\DC_Pinyon\SoilHeatStorage_D_PinyonSoilHeat_2.mat')
elseif site == 4     
    load('C:\towerData\combined\LR_Grass_cleaned.mat')
    site_name = 'Grass'
    AddK =1;
    site_type = 1
    %get the soil heat storage
    load('C:\towerData\EnergyBudgetAnalysis\LR_Grass\SoilHeatStorage_D_grassSoilHeat_2.mat')
elseif site == 5 
    load('C:\towerData\combined\LR_Sage_cleaned.mat')
    site_name = 'Sage'
    AddK =1;
    site_type = 1
    %get the soil heat storage
    load('C:\towerData\EnergyBudgetAnalysis\LR_Sage\SoilHeatStorage_D_sageSoilHeat_2.mat')
elseif site == 6 
    load('C:\towerData\combined\JamesRes_cleaned.mat')
    site_name = 'JamesRes'
    AddK =0;
    site_type = 1
    %get the soil heat storage
    load('C:\towerData\EnergyBudgetAnalysis\JamesRes\SoilHeatStorage_D_JamesResSoilHeat_2.mat')
elseif site == 7
    load('C:\towerData\combined\P301_cleaned.mat')
    site_name = 'P301'
    AddK =0;
    site_type = 3
    %get the soil heat storage
    load('C:\towerData\EnergyBudgetAnalysis\P301\SoilHeatStorage_D_P301SoilHeat_2.mat')
elseif site == 8
    load('C:\towerData\combined\SJER_cleaned.mat')
    site_name = 'SJER'
    AddK =0;
    site_type = 3
    %get the soil heat storage
    load('C:\towerData\EnergyBudgetAnalysis\SJER\SoilHeatStorage_D_SJERSoilHeat_2.mat')
elseif site == 9
    load('C:\towerData\combined\Shorthair_cleaned.mat')
    site_name = 'Shorthair'
    AddK =0;
    site_type = 3
    %get the soil heat storage
    load('C:\towerData\EnergyBudgetAnalysis\Shorthair\SoilHeatStorage_D_ShorthairSoilHeat_2.mat')
elseif site == 10
    load('C:\towerData\combined\Soaproot_cleaned.mat')
    site_name = 'Soaproot'
    AddK =0;
    site_type = 3
    %get the soil heat storage
    load('C:\towerData\EnergyBudgetAnalysis\Soaproot\SoilHeatStorage_D_SoaprootSoilHeat_2.mat')
else
    disp('you did not laod a site')
end

%-----------------------------------------------------------------------%
%Across site fills and calculation for final output
%-----------------------------------------------------------------------%
%soil surface temp - average the surface temp sensors
T_soil_surface = nanmean(Cleaned_D(:, 41:44),2); 

%Fill temp
%consider this if T107 does not agree with the hmp
isn=isnan(D(:,20)) ==1 | isnan(D(:,22)) == 1;
P=polyfit(D(~isn,20),D(~isn,22),1);
Tadj=(D(:,22)-P(1,2))./P(1,1); %consider this

%Fill Data
clear TempFill 
TempFill=Cleaned_D(:,20);
bad=isnan(TempFill(:,1));
TempFill(bad,1)=Cleaned_D(bad,22);%could use Tadj - check this first


missingT=isnan(Cleaned_D(:,20))==1 & isnan(TempFill(:,1))==0;

figure(1)
plot(Cleaned_D(:,1), Cleaned_D(:,20), 'k.')
hold on
plot(Cleaned_D(:,1), Cleaned_D(:,22), 'b.')
hold on
plot(Cleaned_D(missingT,1), TempFill(missingT,1), 'rp')
xlabel('time')
ylabel('temp C')
legend('hmp', 'T107', 'no hmp - yes T107')

pause
close(figure(1))

%do the filling
Cleaned_D(:,20) = TempFill(:,1);

%-----------------------------------------------------------------------%
%fill K radiation across datasets
%-----------------------------------------------------------------------%
if AddK ==1;
    
    %backup the orginal data --> Cleaned_D_back = site Cleaned_D
    Cleaned_D_back = Cleaned_D;
    Cleaned_Header_back = Cleaned_Header;
    
    %remove Cleaned_D to make way for an additional site data
    clear Cleaned_D Cleaned_Header
    
    %this should load a nearby site Cleaned_D
    if site == 1 %pinyon burn
        %fill missing solar in with Pinyon 
        load('C:\towerData\combined\DC_Pinyon_cleaned.mat')
        site_name_fillgaps = 'Pinyon';
    elseif site == 3 %pinyon 
        %this should load in DC_Pinyon data and call it Cleaned_D
        load('C:\towerData\combined\DC_Burn_cleaned.mat'); %This is the best data
        site_name_fillgaps = 'PinyonBurn';
    elseif site == 4 %grass
        %this should load in Sage data and call it Cleaned_D
        load('C:\towerData\combined\LR_Sage_cleaned.mat')
        site_name_fillgaps = 'Sage';
    elseif site == 5 %sage
        %this should load in Grass data and call it Cleaned_D
        load('C:\towerData\combined\LR_Grass_cleaned.mat')
        site_name_fillgaps = 'Grass';       
    end
 
%-----------------------------------------------------------------------%
    %Get common time between the 2 data sets
%-----------------------------------------------------------------------%
    % Where do the time stamps on the 2 datasets overlap?
    Timeback = round(Cleaned_D_back(:,1)*48); % convert to integers and round
    Time = round(Cleaned_D(:,1)*48);
    [commonValues,iTback,iT] = intersect(Timeback, Time);
    
    %Make a K row that looks like Cleaned_D_backup
    clear KMerge
    KMerge = Cleaned_D_back(:,1) *NaN;
    KMerge(iTback,1) = Cleaned_D(iT,18);
    
    %check - These should have strong agreement
    figure(1)
    plot(KMerge(:,1), Cleaned_D_back(:,18), 'k.')
    xlabel('K merge site')
    ylabel('K orginal site')
    
    figure(2)
    plot(Cleaned_D_back(:,1), Cleaned_D_back(:,18), 'ro')
    hold on
    plot(Cleaned_D_back(:,1), KMerge(:,1), 'k.')
    xlabel('time')
    ylabel('K')
    legend('K orginal site', 'K merge site')
    pause
    close(figure(1), figure(2))

%-----------------------------------------------------------------------%
    %do the filling
%-----------------------------------------------------------------------%
    clear Kfill
    KFill=Cleaned_D_back(:,18); %PFburn is the first choice
    missing=isnan(KFill(:,1));
    KFill(missing,1)=KMerge(missing,1);

    disp('Merged nearby site K to this site:')
    site_name
    
    Cleaned_D_back(:,18) = KFill(:,1);
    
%-----------------------------------------------------------------------%    
    %switch name back
%-----------------------------------------------------------------------%
    clear Cleaned_D Cleaned_Header
    
    Cleaned_D = Cleaned_D_back;
    Cleaned_Header = Cleaned_Header_back;
    
else
    
    disp('Added no K to site:')
    site_name
    
end
    
%--------------------------------------------------------------%
%check that the ground heat is lined up
disp('if these are not the exact same number')
size(G(:,1))
size(Cleaned_D(:,1))
soil_line_should_be_zero = nansum(G(:,1)-Cleaned_D(:,1))
doubleCheck_soil_line_should_be_zero = nansum(G(:,2)-Cleaned_D(:,15))

figure(1)
plot(Cleaned_D(:,1), G(:,1), 'k.')
pause
close(figure(1))

figure(2)
plot(Cleaned_D(:,1), G(:,1)-Cleaned_D(:,1), 'k.')
pause
close(figure(2))

figure(3)
plot(Cleaned_D(:,15), G(:,2), 'k.')
pause
close(figure(3))

figure(4)
plot(Cleaned_D(:,1), Cleaned_D(:,15)-G(:,2), 'k.')
pause
close(figure(4))

%--------------------------------------------------------------%

%-----------------------------------------------------------------------%
%Set-up the dataset for output
%-----------------------------------------------------------------------%
len=size(Cleaned_D,1);
wid=43;

DATA = ones(len, wid) * NaN;
HEADER =ones(wid, 1) * NaN;
HEADER = {HEADER};

%fill in the data - remove CO2 concentration, h2o concentration, leaf wetness, 
%fuel moisture, and simplify Temp

%Time
datai=1;
cleanedi=1;

%change time to start at 1/1/2006
if site == 7
    %P301 currently it starts at datenum(2008,1,0)
    %we need to add 2007,and2006
    time = Cleaned_D(:,cleanedi)+365+365;
    DATA(:,datai) = time;
elseif site == 8 %SJER
    %SJER currently it starts at datenum(2009,1,0)
    %we need to add 2008,2007,and2006
    time = Cleaned_D(:,cleanedi)+366+365+365;
    DATA(:,datai) = time;
elseif site == 9
    %Shorthair currently it starts at datenum(2009,1,0)
    %we need to add 2008,2007,and2006
    time = Cleaned_D(:,cleanedi)+ 366+365+365;
    DATA(:,datai) = time;
elseif site == 10 %Soaproot
    %Soaproot currently it starts at datenum(2010,1,0)
    %we need to add 2009,2008,2007,and2006
    time = Cleaned_D(:,cleanedi)+365+366+365+365;
    DATA(:,datai) = time;
else
    %all southern CA sites start at 1/1/2006 ==> no adjustment
    DATA(:,datai) = Cleaned_D(:,cleanedi);
end
    
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 
    
%Ust
datai=2;
cleanedi=4;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%wind direction
datai=3;
cleanedi=5;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%wind speed
datai=4;
cleanedi=6;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%Tair
datai=5;
cleanedi=20;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%sonic temp raw
datai=6;
cleanedi=2;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%sonic temp dried
datai=7;
cleanedi=3;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%rain
datai=8;
cleanedi=23;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%RH
datai=9;
cleanedi=21;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%H2o HMP
datai=10;
cleanedi=24;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%PAR in
datai=11;
cleanedi=16;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%PAR out
datai=12;
cleanedi=17;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%Solar In
datai=13;
cleanedi=18;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%Solar Out
datai=14;
cleanedi=19;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%R net
datai=15;
cleanedi=15;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%Ground heat
datai=16;
str = 'G';
gi=3;%soil heat storage
DATA(:,datai) = G(:,gi); % we will fill this later
HEADER(datai,1)={str}; 

%Sensible Heat
datai=17;
cleanedi=9;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%Latent Heat
datai=18;
cleanedi=10;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%Fco2
datai=19;
cleanedi=13;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%Fh2o
datai=20;
cleanedi=14;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%NDVI
datai=21;
cleanedi=30;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%SoilT
str=['T' '_' 'Soil' '_' 'Surface'];
datai=22;
DATA(:,datai) = T_soil_surface(:,1);
HEADER(datai,1)={str}; 

%Soil Moisture
str=['Soil' '_' 'Moisture' '_' '1'];
datai=23;
cleanedi=35;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)={str}; 

%Soil Moisture
str=['Soil' '_' 'Moisture' '_' '2'];
datai=24;
cleanedi=36;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)={str}; 

%Soil Moisture
str=['Soil' '_' 'Moisture' '_' '3'];
datai=25;
cleanedi=37;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)={str}; 

%Soil Moisture
str=['Soil' '_' 'Moisture' '_' '4'];
datai=26;
cleanedi=38;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)={str}; 

%Matric Potential Sensors have two set-ups - SoCal and Sierra
if site_type == 1; %2 has no sensors
    str=['StartT' '_' '5cm'];
    datai=27;
    cleanedi=48;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str}; 

    str=['StartT' '_' '10cm'];
    datai=28;
    cleanedi=49;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str};  

    str=['StartT' '_' '25cm'];
    datai=29;
    cleanedi=50;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str}; 

    str=['StartT' '_' '50cm'];
    datai=30;
    cleanedi=51;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str}; 

    str=['StartT' '_' '100cm'];
    datai=31;
    cleanedi=52;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str}; 

    str=['StartT' '_' '200cm'];
    datai=32;
    cleanedi=53;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str}; 

    str = ['DelT' '_' '5cm'];
    datai=33;
    cleanedi=54;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str}; 

    str = ['DelT' '_' '10cm'];
    datai=34;
    cleanedi=55;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str}; 

    str = ['DelT' '_' '25cm'];
    datai=35;
    cleanedi=56;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str};  

    str = ['DelT' '_' '50cm'];
    datai=36;
    cleanedi=57;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str};  

    str = ['DelT' '_' '100cm'];
    datai=37;
    cleanedi=58;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str};  

    str = ['DelT' '_' '200cm'];
    datai=38;
    cleanedi=59;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str};  

elseif site_type == 3; 
% these are orginized in the opposite direction --> big number is shallow
%there are 4 sensors - the 5 adn 25 cm depths are missing
    
    str=['StartT' '_' '5cm'];
    datai=27;
    cleanedi=53;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str}; 

    str=['StartT' '_' '10cm'];
    datai=28;
    cleanedi=51;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str};  

    str=['StartT' '_' '25cm'];
    datai=29;
    cleanedi=52;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str}; 

    str=['StartT' '_' '50cm'];
    datai=30;
    cleanedi=50;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str}; 

    str=['StartT' '_' '100cm'];
    datai=31;
    cleanedi=49;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str}; 

    str=['StartT' '_' '200cm'];
    datai=32;
    cleanedi=48;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str}; 

    str = ['DelT' '_' '5cm'];
    datai=33;
    cleanedi=59;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str}; 

    str = ['DelT' '_' '10cm'];
    datai=34;
    cleanedi=57;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str}; 

    str = ['DelT' '_' '25cm'];
    datai=35;
    cleanedi=58;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str};  

    str = ['DelT' '_' '50cm'];
    datai=36;
    cleanedi=56;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str};  

    str = ['DelT' '_' '100cm'];
    datai=37;
    cleanedi=55;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str};  

    str = ['DelT' '_' '200cm'];
    datai=38;
    cleanedi=54;
    DATA(:,datai) = Cleaned_D(:,cleanedi);
    HEADER(datai,1)={str};  
    
end

%source flags
%Ust source flags
datai=39;
cleanedi=25;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%H source flags
datai=40;
cleanedi=26;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%LE source flags
datai=41;
cleanedi=27;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%Fco2 source flags
datai=42;
cleanedi=28;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%Fh2o source flags
datai=43;
cleanedi=29;
DATA(:,datai) = Cleaned_D(:,cleanedi);
HEADER(datai,1)=Cleaned_Header(cleanedi,1); 

%-----------------------------------------------------------------------%
%save and output the data
%-----------------------------------------------------------------------%   
disp('saving outputs for web:')

filename = [out_root site_name '_' 'v3' '_' '2']


%save as a mat file
save(filename, 'DATA', 'HEADER');

%save as a csv file
%filename = [filename '.csv'];
%dlmwrite(filename,DATA, 'precision', '%.6f')

%-----------------------------------------------------------------------%
%Need to make a new Rn at the Sage and Pinyon burn sites
%-----------------------------------------------------------------------%
%this will write over the PinyonBurn file with a fixed up Rn sensor constructed
%from the grass Rn sensor
if site == 1
    Rn_patch_PFBurn()
end

%this will write over the Sage file with a fixed up Rn sensor constructed
%from the grass Rn sensor
if site == 5
    Rn_patch_Sage()
end

