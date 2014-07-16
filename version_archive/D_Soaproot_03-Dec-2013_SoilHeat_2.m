function SoilHeat_2(site)
%G for grass and sage
Date_data_compiled_on = date;

%What site?
%site = 6;

%check the soil output
soilcheck =1;

if site == 1 
    load('C:\towerData\combined\DC_Burn_cleaned.mat')
    site_name = 'D_Burn'
    out_root='C:\towerData\EnergyBudgetAnalysis\DC_Burn\';
    BulkDensity=1.3; %[Mg/m^3]
    site_setup=1; % 6 matric potential sensors, 5, 10, 25, 50, 100, 200 cm depths
elseif site == 2     
    load('C:\towerData\combined\DC_LowDes_cleaned.mat')
    site_name = 'D_LowDes'
    out_root='C:\towerData\EnergyBudgetAnalysis\DC_LowDes\';
    BulkDensity=1.3; %[Mg/m^3]
    site_setup=1; % 6 matric potential sensors, 5, 10, 25, 50, 100, 200 cm depths
elseif site == 3     
    load('C:\towerData\combined\DC_Pinyon_cleaned.mat')
    site_name = 'D_Pinyon'
    out_root='C:\towerData\EnergyBudgetAnalysis\DC_Pinyon\';
    BulkDensity=1.3; %[Mg/m^3]
    site_setup=1; % 6 matric potential sensors, 5, 10, 25, 50, 100, 200 cm depths
elseif site == 4     
    load('C:\towerData\combined\LR_Grass_cleaned.mat')
    site_name = 'D_grass'
    out_root='C:\towerData\EnergyBudgetAnalysis\LR_Grass\';
    BulkDensity=1.3; %[Mg/m^3]
    site_setup=1; % 6 matric potential sensors, 5, 10, 25, 50, 100, 200 cm depths
elseif site == 5 
    load('C:\towerData\combined\LR_Sage_cleaned.mat')
    site_name = 'D_sage'
    out_root='C:\towerData\EnergyBudgetAnalysis\LR_Sage\';
    BulkDensity=1.3; %[Mg/m^3]
    site_setup=1; % 6 matric potential sensors, 5, 10, 25, 50, 100, 200 cm depths
elseif site == 6 
    load('C:\towerData\combined\JamesRes_cleaned.mat')
    site_name = 'D_JamesRes'
    out_root='C:\towerData\EnergyBudgetAnalysis\JamesRes\';
    BulkDensity=1.3; %[Mg/m^3]
    site_setup=1; % 6 matric potential sensors, 5, 10, 25, 50, 100, 200 cm depths
elseif site == 7
    load('C:\towerData\combined\P301_cleaned.mat')
    site_name = 'D_P301';
    out_root='C:\towerData\EnergyBudgetAnalysis\P301\';
    BulkDensity=1.3; %[Mg/m^3]
    site_setup=3; % 4 matric potential sensors, 10, 50, 100, 200 cm depths
elseif site == 8
    load('C:\towerData\combined\SJER_cleaned.mat')
    site_name = 'D_SJER';
    out_root='C:\towerData\EnergyBudgetAnalysis\SJER\';
    BulkDensity=1.3; %[Mg/m^3]
    site_setup=3; % 4 matric potential sensors, 10, 50, 100, 200 cm depths
elseif site == 9
    load('C:\towerData\combined\Shorthair_cleaned.mat')
    site_name = 'D_Shorthair';
    out_root='C:\towerData\EnergyBudgetAnalysis\Shorthair\';
    BulkDensity=1.3; %[Mg/m^3]
    site_setup=3; % 4 matric potential sensors, 10, 50, 100, 200 cm depths
elseif site == 10
    load('C:\towerData\combined\Soaproot_cleaned.mat')
    site_name = 'D_Soaproot';
    out_root='C:\towerData\EnergyBudgetAnalysis\Soaproot\';
    BulkDensity=1.3; %[Mg/m^3]
    site_setup=3; % 4 matric potential sensors, 10, 50, 100, 200 cm depths
end  
%-----------------------------------------------------------------------%
%calulate the heat storage
%delQ=Cv*depth*(T(i)-T(i-1))/(time(i)-time(i-1)) [W/m^2]
%----%a closer look at units:
%----%Cv - [J/C/m^3 of soil space], depth - [m of soil space],...
%----%(T(i)-T(i-1)) - [C] ==> [J/C/m^3] * [m] * [C] = [J/m^2]
%determine change in Q over the 30 min = delQ(i)
%----%delQ(i) = (Q(i) - Q(i-1))/(time(i)-time(i-1))
%----%==> J/m^2/30 min ==> 30min/(30*60s) ==> [J/m^2/s] ==> [W/m^2]
%-----------------------------------------------------------------------%
%how many rows in Cleaned_D?
len = size(Cleaned_D,1);
o=ones(len,1);
%-----------------------------------------------------------------------%
%The problem: The diel variations in soil temperature are exaggerated at a couple periods during
%the day. ==> We need to remove the exaggerated temperatures to determine soil heat storage.

%The matric potential sensors use a thermocouple to measure soil temperature. 
%The thermocouples need a reference temperature at the data logger to determine the 
%soil temperature. The exaggerated soil temperature periods correspond to large changes 
%in the thermocouple's reference temperature.  
%The idea: the soil data logger warms up in an uneven way.  
%The temperature at the thermocouple data logger junction is not the same 
%temperature as the reference temperature.  
%==> The soil temperature is incorrect.

%Assumptions:
%The exaggeration in temperature is the same for all thermocouples
%The deepest thermocouple temperature does not change.
%Differencing each thermocouple with the deepest thermocouple will give the
%correct variation in soil temperature (not the correct absolute
%temperature).

%-----------------------------------------------------------------------%
%Fix up the temperature so it is usable (see above) [degree C]
if site_setup==1; 
    T5=Cleaned_D(:,48)-Cleaned_D(:,53);
    T10=Cleaned_D(:,49)-Cleaned_D(:,53);
    T25=Cleaned_D(:,50)-Cleaned_D(:,53);
    T50=Cleaned_D(:,51)-Cleaned_D(:,53);
    T100=Cleaned_D(:,52)-Cleaned_D(:,53);
    T200=Cleaned_D(:,53)-Cleaned_D(:,53);
elseif site_setup ==3;
    T10=Cleaned_D(:,51)-Cleaned_D(:,48);
    T50=Cleaned_D(:,50)-Cleaned_D(:,48);
    T100=Cleaned_D(:,49)-Cleaned_D(:,48);
    T200=Cleaned_D(:,48)-Cleaned_D(:,48);
end

%-----------------------------------------------------------------------%
%for the So Cal sites - 6 soil matric potential sensors [degree C]
%---the sensors are at 5, 10, 25, 50, 100, and 200 cm depths

%However, we to get a weighted T for the top 30cm that corresponds to the V soil water
%measurements
if site_setup==1; 
    
    %-----------------------------------------------------------------------%
    %T in the top 30cm of soil
    T30 = (T5 *(7.5/30)) + (T10 *(10/30)) + (T25 *(12.5/30)); %weight the T over the top 30cm of soil
    T34 = T25; %remaining depth of soil closest to the 25cm probe

    Ts= [T30, T34, T50, T100, T200];      

    %-----------------------------------------------------------------------%
    %delT nans [degree C]
    %(T(i)-T(i-1))
    delT30=ones(len,1)*NaN;
    delT34=ones(len,1)*NaN;
    delT50=ones(len,1)*NaN;
    delT100=ones(len,1)*NaN;
    delT200=ones(len,1)*NaN;

    %Units [C]
    for i = 2:len
        delT30(i)=T30(i,1)-T30((i-1),1);
        delT34(i)=T34(i,1)-T34((i-1),1);
        delT50(i)=T50(i,1)-T50((i-1),1);
        delT100(i)=T100(i,1)-T100((i-1),1);
        delT200(i)=T200(i,1)-T200((i-1),1);
    end

    %store these 
    delTs=[delT30, delT34, delT50, delT100, delT200];
    %-----------------------------------------------------------------------%
elseif site_setup == 3

    %-----------------------------------------------------------------------%
    T30 = T10;
    T34 = T10 .*NaN;
    
    Ts= [T30, T34, T50, T100, T200]; 
    %-----------------------------------------------------------------------%
    %Units [C]
    for i = 2:len
        delT30(i)=T30(i,1)-T30((i-1),1);
        delT50(i)=T50(i,1)-T50((i-1),1);
        delT100(i)=T100(i,1)-T100((i-1),1);
        delT200(i)=T200(i,1)-T200((i-1),1);
    end

    delT30=delT30';
    delT50=delT50';
    delT100=delT100';
    delT200=delT200';
        
    delT34=delT30.*NaN;
    
    %store these 
    delTs=[delT30, delT34, delT50, delT100, delT200];
    
end

%-----------------------------------------------------------------------%
%depth [m]
if site_setup==1; 
    d30=0.30*o; %[m of soil]
    d34=0.075*o; %[m of soil]
    d50=0.375*o; %[m of soil]
    d100=0.75*o; %[m of soil]
    d200=0.5*o; %[m of soil]
elseif site_setup==3;
    d30=0.30*o; %[m of soil]
    d34=0*o; %[m of soil]
    d50=0.45*o; %[m of soil]
    d100=0.75*o; %[m of soil]
    d200=0.5*o; %[m of soil]
end
%-----------------------------------------------------------------------%
%How much water is in the top 30 cm
Vwater=nanmean(Cleaned_D(:,35:38),2);  %fraction or [m^3 water/m^3 soil]

%store these 
Vs=Cleaned_D(:,35:38);

%We only have matric potential below 30cm.
%There are too many assumptions to get water volume below 30cm &
%temperature does not change much below 30cm ==> ignore this
%-----------------------------------------------------------------------%
%-----------------------------------------------------------------------%
%volumetric heat capacity, Cv 
%volumetric heat capacity of a soil is the sum of the heat capacities of the
%soil components (Campell and Norman, p. 117)
%---- assume water and soil only - ignore air contribution (should be small)
%-----------------------------------------------------------------------%

%Density of soil minerals (Campell and Norman, p. 118)
Cs_minerals = 870*o; %[J/kg/C]
BulkDensity=(BulkDensity*1000)*o; %[kg/m^3]

%Water info
dens_water = (1000)*o; % kg of water/m^3 (Campell and Norman, p. 118)
Cs_water = (4180)*o; % J/kg of water/K (Campell and Norman, p. 118)

%Units [J/kg mineral/C] * [kg mineral /m^3 soil] ==> J/m^3 of soil/C
%      [J/kg of water/C] * [kg of water/m^3 of water] * [m^3 of water/m^3 of soil] ==> J/m^3 of soil/C
Cv30= (Cs_minerals .* BulkDensity) + (Vwater .* dens_water .* Cs_water); %units J/C/m^3 of soil space

%We have no direct estimate of water storage in these lower layers
Cv34= (Cs_minerals .* BulkDensity); %units J/C/m^3 of soil space
Cv50= (Cs_minerals .* BulkDensity); %units J/C/m^3 of soil space
Cv100= (Cs_minerals .* BulkDensity); %units J/C/m^3 of soil space
Cv200= (Cs_minerals .* BulkDensity); %units J/C/m^3 of soil space

%store these 
Cvs=[Cv30, Cv34, Cv50, Cv100, Cv200];
%-----------------------------------------------------------------------%
%change in time 
%Units (30min * (60 s/1 min))
del_time = (30 * 60) * o;%[s]
%-----------------------------------------------------------------------%
%delQ=Cv*depth*(T(i)-T(i-1))/(time(i)-time(i-1)) [W/m^2]

%Units:
%[J/C/m^3 of soil space] * [m of soil] *[C] / [s] ==>
%J/m^2/s ==> W/m^2
%-----------------------------------------------------------------------%

%heat storage in each soil layer [W/m^2]
delQ30 =Cv30 .* d30 .*delT30 ./ del_time;
delQ34 =Cv34 .* d34 .*delT34 ./ del_time;
delQ50 =Cv50 .* d50 .* delT50 ./ del_time;
delQ100 =Cv100 .*d100 .* delT100 ./ del_time;
delQ200 =Cv200 .*d200 .* delT200 ./ del_time;

%store these
Qs=[delQ30, delQ34, delQ50, delQ100, delQ200];
%-----------------------------------------------------------------------%
%heat stored over the top 2m of soil
if site_setup ==1
    delQ = delQ30 + delQ34 + delQ50 + delQ100 + delQ200; % [W/m^2]
elseif site_setup ==3
    delQ = delQ30 + delQ50 + delQ100 + delQ200; % [W/m^2]
end
%-----------------------------------------------------------------------%
G=[Cleaned_D(:,1), Cleaned_D(:,15), delQ(:,1)]; %time; Rn; G
%save the G
save([out_root 'SoilHeatStorage_' site_name 'SoilHeat_2' '.mat'], 'Date_data_compiled_on', 'G');
%-----------------------------------------------------------------------%
%Check the work and eenrgy fluxes
if soilcheck ==1
    soilheatcheck_2(Cleaned_D, Ts, delTs, Vwater, Vs, Cvs, Qs, delQ, site_name, out_root);
else
    disp('you did not check the soil heat flux outputs')
end
%-----------------------------------------------------------------------%
%site energy budget closure
%-----------------------------------------------------------------------%
%energy budget closure - we want turbulent periods
turb=Cleaned_D(:,4) > 0.3 & isfinite(Cleaned_D(:,9)) & isfinite(Cleaned_D(:,10)) & isfinite(Cleaned_D(:,15)) & isfinite(delQ(:,1));
%Availible energy
ae = Cleaned_D(turb,15)-delQ(turb,1);
%Net Radiation
Rn = Cleaned_D(turb,15);
%sensible + latent heat
HLE=Cleaned_D(turb,9)+Cleaned_D(turb,10);

disp('%-----------------------------------------------------------------------%')
disp('Are we keeping a sufficient number of observations?')
disp('%-----------------------------------------------------------------------%')
%n_obs
n_obs = sum(turb) %the number of observations we are using, obs removed resulting from missing sensors 

possible_obs=Cleaned_D(:,4) > 0.3; %There are this many turbulent observations
len_poss = sum(possible_obs); %number of turbulent observations

len = size(turb,1) %length of entire matrix


disp('turbulent and finite energy fluxes /(obs with Cleaned_D(:,4) > 0.3)')
f_of_possible_obs_Ustfiltered = n_obs/len_poss

disp('turbulent and finite energy fluxes /(all obs)')
f_of_possible_obs_total = n_obs/len

disp('Cleaned_D(:,4) > 0.3 obs /(all obs)')
f_Ustfiltered_obs_total = len_poss/len
disp('%-----------------------------------------------------------------------%')
%-----------------------------------------------------------------------%
site_name

disp('%-----------------------------------------------------------------------%')
disp('Rn-soil heat storage(AE) vs H + LE')
disp('%-----------------------------------------------------------------------%')
%with zero intercept
energy_closure_zero = ae \ HLE

%using polyfit
energy_closure_polyfit = polyfit(ae, HLE, 1)

%robustfit
energy_closure_robustfit = robustfit(ae,HLE)
disp('%-----------------------------------------------------------------------%')
%-----------------------------------------------------------------------%
disp('%-----------------------------------------------------------------------%')
disp('Rn vs H + LE ---- No soil heat storage')
disp('%-----------------------------------------------------------------------%')
energy_closure_zero_no_delQ = Rn \ HLE

%using polyfit
energy_closure_polyfit_no_delQ = polyfit(Rn, HLE, 1)

%robustfit
energy_closure_robustfit_no_delQ = robustfit(Rn,HLE)
disp('%-----------------------------------------------------------------------%')
%-----------------------------------------------------------------------%
%save the energy budget information
EnergyBudgetInfo = [n_obs; len_poss; len;...
    f_of_possible_obs_Ustfiltered; f_of_possible_obs_total; f_Ustfiltered_obs_total;...
    energy_closure_zero;...
    energy_closure_polyfit(1,1); energy_closure_polyfit(1,2);...
    energy_closure_robustfit(2,1); energy_closure_robustfit(1,1);...
    energy_closure_zero_no_delQ;...
    energy_closure_polyfit_no_delQ(1,1); energy_closure_polyfit_no_delQ(1,2);...
    energy_closure_robustfit_no_delQ(2,1); energy_closure_robustfit_no_delQ(1,1)];

Header_EnergyBudgetInfo = {'n_obs_with_AE_Ust_filtered'; 'n_obs_Ust_filtered'; 'n_obs_all';...
    'n_obs_with_AE_Ust_filtered/n_obs_Ust_filtered'; 'n_obs_with_AE_Ust_filtered/n_obs_all'; 'n_obs_Ust_filtered/n_obs_all';...
    'energy_closure_zero_AE';...
    'energy_closure_polyfit_AE_slope'; 'energy_closure_polyfit_AE_intercept';...
    'energy_closure_robustfit_AE_slope'; 'energy_closure_robustfit_AE_intercept';...
    'energy_closure_zero_no_delQ';...
    'energy_closure_polyfit_no_delQ_slope'; 'energy_closure_polyfit_no_delQ_intercept';...
    'energy_closure_robustfit_no_delQ_slope'; 'energy_closure_robustfit_no_delQ_intercept'};

save([out_root 'EnergyBudgetInfo_' site_name '.mat'], 'Date_data_compiled_on', 'EnergyBudgetInfo', 'Header_EnergyBudgetInfo');
%-----------------------------------------------------------------------%
%check out the scatter plot
%Availible energy
figure(1)
plot(ae, HLE, 'k.')
hold on
plot(ae, energy_closure_zero*ae, 'r-')
hold on
plot(ae, energy_closure_polyfit(1,1) *ae + energy_closure_polyfit(1,2), 'g-')
hold on
plot(ae, energy_closure_robustfit(2,1)*ae + energy_closure_robustfit(1,1), 'b-')
title(site_name)
legend('ae v HLE', 'energy_closure_zero', 'energy_closure_polyfit', 'energy_closure_robustfit');
xlabel('Rn - soil heat storage (W m^-^2)')
ylabel('H + LE (W m^-^2)')

f='EnergyBudgetClosure_AE';
saveas(figure(1), [out_root f]);

%Net Radiation
figure(2)
plot(Rn, HLE, 'k.')
hold on
plot(Rn, energy_closure_zero_no_delQ*Rn, 'r-')
hold on
plot(Rn, energy_closure_polyfit_no_delQ(1,1) *Rn + energy_closure_polyfit_no_delQ(1,2), 'g-')
hold on
plot(Rn, energy_closure_robustfit_no_delQ(2,1)*Rn + energy_closure_robustfit_no_delQ(1,1), 'b-')
title(site_name)
legend('Rn v HLE', 'energy_closure_zero_no_delQ', 'energy_closure_polyfit_no_delQ', 'energy_closure_robustfit_no_delQ');
xlabel('Rn (W m^-^2)')
ylabel('H + LE (W m^-^2)')

f='EnergyBudgetClosure_Rn';
saveas(figure(1), [out_root f]);

copyfile('SoilHeat_2.m',['./version_archive/' site_name '_' date '_SoilHeat_2.m'],'f');
