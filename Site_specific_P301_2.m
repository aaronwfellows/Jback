function [D, TEMP, H2O, rho_fix, v_wind] = Site_specific_P301_2(HEADER, D, Press, R_mol, Tc, a)


%The following is the data from the SiteDetails xls spreadsheet.  
%i added a few dates that corresponded with some jumps in the irga
%calibration -- This should help. - awf
IRGAchanged = [9,4,2008; 10,1,2010; 10,6,2011; 6,28,2012; 9,26,2013]; %;9,29,2012]; %This needs to be updated 5/5/2012 - awf
    
%if you want to fix up the incoming radiation 
if a==1
%Time stamp was shifted to California time from day 528 or 529 to 677

timeshifted = find(D(1,:) >= 528.6 & D(1,:) <= 676.74); %column numbers
gototime=timeshifted+round(48*(6/24)); %daylight savings - GMT is ahead of California
D(2:end, gototime) = D(2:end, timeshifted);

%NaN out the days around the shift - we don't really know what appened
b_e = min(gototime);
b_e = D(1, b_e);
b_l = max(gototime);
b_l = D(1, b_l);

buffer = D(1,:) > (b_e - 0.5) &  D(1,:) < (b_e + 0.5);
buffer = buffer | D(1,:) > (b_l - 0.5) &  D(1,:) < (b_l + 0.5);
D(2:end, buffer)=NaN;

%This file is a site specific clean-up for DMERGE data at the P301 in the Sierra Nevada Mountains
%awf 2/9/2013

%OUTPUTS:
%D is the filtered DMERGE data file
%TEMP is a cleaned and filled temperature in deg C 
%H2O is the water mixing ratio [mmol h2o/ mol of dry air]
%rho_fix is bad temp and water vapor measurements
%----> rho_fix is a handle that points to fast calculated
%----> fluxes that used bad dry density estimates 
%v_wind is cleaned up bad wind components


%----------------------------------------------------------------------
%fix up the goes wind problem in the u direction
%----------------------------------------------------------------------

%use ones through the program
len=length(D(1,:));
o=ones(1,len);

%use time to remove bad data
time = D(1,:);

adj = time < 1005.74;
D(207,adj)=D(207,adj)*100;

   
    
    correct_shading = 0; %1 if you want to fill in the shade
    
    if correct_shading ==1
    %---------------------------------------------------------------------
    %% PAR_IN, PYRR_IN
    %---------------------------------------------------------------------

    % this is pretty crude, but it seems to work better than any modeling
    % thing I tried
    uniqdays = unique(floor(D(1,:)));
    % incoming data columns
    incoming = [88 90 214 216];
    for jj=1:4
        column = incoming(jj);
        bad = D(column,:);
        for ii=nanmin(uniqdays):nanmax(uniqdays)
            dayii  = floor(D(1,:)) == ii & D(1,:)-floor(D(1,:)) > 0.8;
            daymax = nanmax(D(column,dayii));
            if ~isnan(daymax)
                maxii  = find(D(column,:) == daymax & floor(D(1,:))==ii);
                % convert to percentage of max
                %bad(maxii-4:maxii-1) = D(column,maxii-4:maxii-1)./daymax;
                %bad(maxii-4:maxii-1)
                %bad(maxii-4:maxii-1) = bad(maxii-4:maxii-1) + corr';
                %bad(maxii-4:maxii-1)
                %pause
                %bad(maxii-4:maxii-1) = bad(maxii-4:maxii-1).*daymax;
                bad(maxii-5:maxii-1) = cos(pi*[5:-1:1]/24 - pi/48)' * daymax;
                %orr(column,maxii-4:maxii-1) = D(column,maxii-4:maxii-1)./bad';
            end
        end
        D(column,:) = bad;
    end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% By measurement
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %----------------------------------------------------------------------
    %% HMP_Temp 
    %----------------------------------------------------------------------
    %data logger
    badHMP_dl = isnan(D(97,:)) | isempty(D(97,:)) | ~isfinite(D(97,:));
    badHMP_dl = badHMP_dl | D(97,:)<-20;
    badHMP_dl = badHMP_dl | time > 592.124 & time < 592.126;
    badHMP_dl = badHMP_dl | time > 1095.435 & time < 1095.44;
    badHMP_dl = badHMP_dl | time > 1226.7 & time < 1226.85;
    D(97,badHMP_dl) = NaN;
    
    %goes
    badHMP_goes = isnan(D(223,:)) | isempty(D(223,:)) | ~isfinite(D(223,:));
    badHMP_goes = badHMP_goes | D(223,:)<-20;
    badHMP_goes = badHMP_goes | time < 249.231;
    badHMP_goes = badHMP_goes | time > 592.124 & time < 592.126;
    badHMP_goes = badHMP_goes | time > 528.66 & time < 528.8;
    badHMP_goes = badHMP_goes | time > 624.26 & time < 624.28;
    badHMP_goes = badHMP_goes | time > 729.42 & time < 729.46;
    badHMP_goes = badHMP_goes | time > 801.435 & time < 801.445;
    badHMP_goes = badHMP_goes | time > 1095.435 & time < 1095.44;
    badHMP_goes = badHMP_goes | time > 1776.43 & time < 1776.44;
    
    focus = time > 271.5 & time < 274;
    badHMP_goes = badHMP_goes | focus & D(223,:) > 22;
    
    D(223,badHMP_goes) = NaN;

    %build a consensus HMP_Temp [C]
    thmp=D(97,:);%HMP Temp
    bad=isnan(thmp);
    thmp(bad)=D(223,bad); %GOES TEMP
    
    %----------------------------------------------------------------------
    %% Tsonic
    %----------------------------------------------------------------------
    %clean mat
    bad_Ts_mat=isnan(D(8,:)) | isempty(D(8,:)) | ~isfinite(D(8,:));
    bad_Ts_mat = bad_Ts_mat | time > 307.65 & time < 308.88;
    bad_Ts_mat = bad_Ts_mat | time > 488.6 & time < 489.26;
    bad_Ts_mat = bad_Ts_mat | time > 528.6249 & time < 528.63;
    bad_Ts_mat = bad_Ts_mat | time > 722.15 & time < 723.2;
    bad_Ts_mat = bad_Ts_mat | time > 749.26 & time < 749.75;
    bad_Ts_mat = bad_Ts_mat | time > 755.43 & time < 755.74;
    bad_Ts_mat = bad_Ts_mat | time > 767.6 & time < 770.6;
    bad_Ts_mat = bad_Ts_mat | time > 782.3 & time < 782.5;
    bad_Ts_mat = bad_Ts_mat | time > 786.76 & time < 787.385;
    bad_Ts_mat = bad_Ts_mat | time > 798.8115 & time < 798.13;
    bad_Ts_mat = bad_Ts_mat | time > 820.87 & time < 821.45;
    bad_Ts_mat = bad_Ts_mat | time > 841.955 & time < 842.12;
    bad_Ts_mat = bad_Ts_mat | time > 877.62 & time < 877.89;
    bad_Ts_mat = bad_Ts_mat | time > 878.54 & time < 878.62;
    bad_Ts_mat = bad_Ts_mat | time > 946.995 & time < 947.002;
    bad_Ts_mat = bad_Ts_mat | time > 1027.49 & time < 1029.53;
    bad_Ts_mat = bad_Ts_mat | time > 1054.93 & time < 1062.4;
    bad_Ts_mat = bad_Ts_mat | time > 1071.135 & time < 1071.61;
    bad_Ts_mat = bad_Ts_mat | time > 1080.6 & time < 1082.18;
    bad_Ts_mat = bad_Ts_mat | time > 1082.66 & time < 1090;
    bad_Ts_mat = bad_Ts_mat | time > 1199.72 & time < 1200.45;
    bad_Ts_mat = bad_Ts_mat | time > 1224.66 & time < 1225.3;
    bad_Ts_mat = bad_Ts_mat | time > 1231.52 & time < 1231.525;
    bad_Ts_mat = bad_Ts_mat | time > 1245.5 & time < 1246.51;
    bad_Ts_mat = bad_Ts_mat | time > 1253.6 & time < 1254;
    bad_Ts_mat = bad_Ts_mat | time > 1276.7 & time < 1277.055;
    bad_Ts_mat = bad_Ts_mat | time > 1412.2 & time < 1414.35;
    bad_Ts_mat = bad_Ts_mat | time > 1481.88 & time < 1482.4;
    bad_Ts_mat = bad_Ts_mat | time > 1499.94 & time < 1500.595;
    bad_Ts_mat = bad_Ts_mat | time > 1503.87 & time < 1503.88;
    bad_Ts_mat = bad_Ts_mat | time > 1538.24 & time < 1538.4;
    bad_Ts_mat = bad_Ts_mat | time > 1557.025 & time < 1557.68;
    bad_Ts_mat = bad_Ts_mat | time > 1578.3 & time < 1579.1;
    D(8,bad_Ts_mat)=NaN;

    %clean dl
    bad_Ts_dl = isnan(D(84,:)) | isempty(D(84,:)) | ~isfinite(D(84,:));
    bad_Ts_dl = bad_Ts_dl | time > 488 & time < 489.26;
    bad_Ts_dl = bad_Ts_dl | time > 652.2 & time < 653.68;
    bad_Ts_dl = bad_Ts_dl | time > 692.03 & time < 692.18;
    bad_Ts_dl = bad_Ts_dl | time > 711.39 & time < 715.82;
    bad_Ts_dl = bad_Ts_dl | time > 722.15 & time < 723.2;
    bad_Ts_dl = bad_Ts_dl | time > 749.26 & time < 749.75;
    bad_Ts_dl = bad_Ts_dl | time > 755.43 & time < 755.74;
    bad_Ts_dl = bad_Ts_dl | time > 767.6 & time < 770.6;
    bad_Ts_dl = bad_Ts_dl | time > 782.3 & time < 782.5;
    bad_Ts_dl = bad_Ts_dl | time > 786.76 & time < 787.385;
    bad_Ts_dl = bad_Ts_dl | time > 798.8115 & time < 798.13;
    bad_Ts_dl = bad_Ts_dl | time > 820.87 & time < 821.45;
    bad_Ts_dl = bad_Ts_dl | time > 841.955 & time < 842.12;
    bad_Ts_dl = bad_Ts_dl | time > 849.2 & time < 850;
    bad_Ts_dl = bad_Ts_dl | time > 877.62 & time < 877.89;
    bad_Ts_dl = bad_Ts_dl | time > 878.54 & time < 878.62;
    bad_Ts_dl = bad_Ts_dl | time > 1027.49 & time < 1029.53;
    bad_Ts_dl = bad_Ts_dl | time > 1054.93 & time < 1062.4;
    bad_Ts_dl = bad_Ts_dl | time > 1071.135 & time < 1071.61;
    bad_Ts_dl = bad_Ts_dl | time > 1080.6 & time < 1082.18;
    bad_Ts_dl = bad_Ts_dl | time > 1082.66 & time < 1090;
    bad_Ts_dl = bad_Ts_dl | time > 1199.72 & time < 1200.45;
    bad_Ts_dl = bad_Ts_dl | time > 1224.66 & time < 1225.3;
    bad_Ts_dl = bad_Ts_dl | time > 1226.7 & time < 1226.85;
    bad_Ts_dl = bad_Ts_dl | time > 1231.52 & time < 1231.525;
    bad_Ts_dl = bad_Ts_dl | time > 1245.5 & time < 1246.51;
    bad_Ts_dl = bad_Ts_dl | time > 1253.6 & time < 1254;
    bad_Ts_dl = bad_Ts_dl | time > 1276.7 & time < 1277.055;
    bad_Ts_dl = bad_Ts_dl | time > 1412.2 & time < 1414.35;
    bad_Ts_dl = bad_Ts_dl | time > 1481.88 & time < 1482.4;
    bad_Ts_dl = bad_Ts_dl | time > 1499.94 & time < 1500.595;
    bad_Ts_dl = bad_Ts_dl | time > 1503.87 & time < 1503.88;
    bad_Ts_dl = bad_Ts_dl | time > 1538.24 & time < 1538.4;
    bad_Ts_dl = bad_Ts_dl | time > 1557.025 & time < 1557.68;
    bad_Ts_dl = bad_Ts_dl | time > 1578.3 & time < 1579.1;
    D(84,bad_Ts_dl)=NaN;

    %clean GOES
    bad_Ts_goes = isnan(D(210,:)) | isempty(D(210,:)) | ~isfinite(D(210,:));
    bad_Ts_goes = bad_Ts_goes | D(206,:) ==0;
    bad_Ts_goes = bad_Ts_goes | time > 272 & time < 274;
    bad_Ts_goes = bad_Ts_goes | time > 488 & time < 489.26;
    bad_Ts_goes = bad_Ts_goes | time > 624.26 & time < 624.28;
    bad_Ts_goes = bad_Ts_goes | time > 652.2 & time < 653.68;
    bad_Ts_goes = bad_Ts_goes | time > 692.03 & time < 692.18;
    bad_Ts_goes = bad_Ts_goes | time > 711.39 & time < 715.82;
    bad_Ts_goes = bad_Ts_goes | time > 722.15 & time < 723.2;
    bad_Ts_goes = bad_Ts_goes | time > 749.26 & time < 749.75;
    bad_Ts_goes = bad_Ts_goes | time > 755.43 & time < 755.74;
    bad_Ts_goes = bad_Ts_goes | time > 767.6 & time < 770.6;
    bad_Ts_goes = bad_Ts_goes | time > 782.3 & time < 782.5;
    bad_Ts_goes = bad_Ts_goes | time > 786.76 & time < 787.385;
    bad_Ts_goes = bad_Ts_goes | time > 798.8115 & time < 798.13;
    bad_Ts_goes = bad_Ts_goes | time > 820.87 & time < 821.45;
    bad_Ts_goes = bad_Ts_goes | time > 841.955 & time < 842.12;
    bad_Ts_goes = bad_Ts_goes | time > 849.2 & time < 850;
    bad_Ts_goes = bad_Ts_goes | time > 877.62 & time < 877.89;
    bad_Ts_goes = bad_Ts_goes | time > 878.54 & time < 878.62;
    bad_Ts_goes = bad_Ts_goes | time > 1027.49 & time < 1029.53;
    bad_Ts_goes = bad_Ts_goes | time > 1054.93 & time < 1062.4;
    bad_Ts_goes = bad_Ts_goes | time > 1071.135 & time < 1071.61;
    bad_Ts_goes = bad_Ts_goes | time > 1080.6 & time < 1082.18;
    bad_Ts_goes = bad_Ts_goes | time > 1082.66 & time < 1090;
    bad_Ts_goes = bad_Ts_goes | time > 1199.72 & time < 1200.45;
    bad_Ts_goes = bad_Ts_goes | time > 1224.66 & time < 1225.3;
    bad_Ts_goes = bad_Ts_goes | time > 1231.52 & time < 1231.525;
    bad_Ts_goes = bad_Ts_goes | time > 1245.5 & time < 1246.51;
    bad_Ts_goes = bad_Ts_goes | time > 1253.6 & time < 1254;
    bad_Ts_goes = bad_Ts_goes | time > 1276.7 & time < 1277.055;
    bad_Ts_goes = bad_Ts_goes | time > 1412.2 & time < 1414.35;
    bad_Ts_goes = bad_Ts_goes | time > 1481.88 & time < 1482.4;
    bad_Ts_goes = bad_Ts_goes | time > 1499.94 & time < 1500.595;
    bad_Ts_goes = bad_Ts_goes | time > 1503.87 & time < 1503.88;
    bad_Ts_goes = bad_Ts_goes | time > 1538.24 & time < 1538.4;
    bad_Ts_goes = bad_Ts_goes | time > 1557.025 & time < 1557.68;
    bad_Ts_goes = bad_Ts_goes | time > 1578.3 & time < 1579.1;
    bad_Ts_goes = bad_Ts_goes | time > 1784.5 & time < 1785.36;
    D(210,bad_Ts_goes)=NaN;

    %----------------------------------------------------------------------
    %%   Water Vapor Pressure (HMP Probe) - convert RH to mixing ratio
    %----------------------------------------------------------------------
    %we cleaned up the hmp temp 
    %==> fix up h2o_hmp calculation that were done in the data logger

    %Lowe ~1977
    A_0=6.107799961;
    A_1=0.4436518521;
    A_2=0.01428945805;
    A_3=0.0002650648471;
    A_4=0.000003031240396;
    A_5=0.00000002034080948;
    A_6=0.00000000006136820929;

    %need to use the hmp temp because this is where the humidity is being
    %measured
    esat=0.1*(A_0+thmp.*(A_1+thmp.*(A_2+thmp.*(A_3+thmp.*(A_4+thmp.*(A_5+thmp.*A_6)))))); %[kPa]
 
    %RH - clean-up relative humidity
    %----------------------------------------------------------------------
    %% HMP_RH
    %----------------------------------------------------------------------
    
    %------------------------------------------------
    %There may be some bad HMP_RH - ex day 782 to 784
    %and GOES RH < day 475????
    %------------------------------------------------
    
    %data logger
    bad= D(98,:) < 0.01;
    bad= bad | time < 465.9;
    bad= bad | badHMP_dl; %In this case this works
    D(98,bad) = NaN;

    %This is not the best way to correct these numbers and will introduce a
    %small bias but I could not tell if there was a gain problem or an
    %offset problem 
    adjust = D(98,:) >= 1;
    D(98,adjust) = 1;
    
    %goes
    bad= D(224,:) < 0.01;
    bad= bad | time < 465.9;
    bad= bad | D(224,:) > 10;
    bad= bad | badHMP_goes;
    D(224,bad) = NaN;

    %problem
    adjust = D(224,:) >= 1;
    D(224,adjust) = 1;
    
    e=D(98,:).*esat;%[kPa]
    mol_frac_dl = e/Press; %mol h2o/mol wet air
    mixing_dl=mol_frac_dl./(1-mol_frac_dl); %mol h2o/mol dry air
    D(100,:)= mixing_dl*1000; %mmol h2o/mol dry air

    e=D(224,:).*esat;%[kPa]
    mol_frac_goes = e/Press;  
    mixing_goes=mol_frac_goes./(1-mol_frac_goes); %mol h2o/mol dry air
    D(226,:)=mixing_goes*1000; %mmol h2o/mol dry air

    %----------------------------------------------------------------------
    %% IRGA H2O concentrations - mmol h20/mol dry air
    %carry bad IRGA H2O concentrations  through the water fluxes
    %----------------------------------------------------------------------
    %clean fast
    bad_H2Oc_mat = isnan(D(42,:)) | isempty(D(42,:)) | ~isfinite(D(42,:));
    bad_H2Oc_mat = bad_H2Oc_mat | D(42, :) > 29.8; %h2o conc are saturated at 30 mmol/mol 
    bad_H2Oc_mat = bad_H2Oc_mat | D(42, :) < 0; %h2o conc are saturated at 30 mmol/mol 
    
    focus = time > 2080 & time < 2160;
    bad_H2Oc_mat = bad_H2Oc_mat | (focus & D(42, :) < 1); %h2o conc are saturated at 30 mmol/mol 
    D(42, bad_H2Oc_mat) = NaN;

  
    %clean data logger
    bad_H2Oc_dl = isnan(D(86,:)) | isempty(D(86,:)) | ~isfinite(D(86,:));
    bad_H2Oc_dl = bad_H2Oc_dl | D(86, :) > 29.8 | D(86,:) < 0; %h2o conc are saturated at 30 mmol/mol 
    bad_H2Oc_dl = bad_H2Oc_dl | D(78, :) ==0; %no variabiltiy 
    bad_H2Oc_dl = bad_H2Oc_dl | time > 2092.8 & time < 2096.745; 
    bad_H2Oc_dl = bad_H2Oc_dl | time > 2121.38 & time < 2125.4; 
    
    D(86, bad_H2Oc_dl) = NaN;

   
    %convert to mixing ratio
    mixing_dl=D(86, :)./(1-(D(86, :)/1000)); %mol h2o/mol dry air
    D(86, :)= mixing_dl; %mmol h2o/mol dry air

    %clean GOES
    bad_H2Oc_goes = isnan(D(212,:)) | isempty(D(212,:)) | ~isfinite(D(212,:));
    bad_H2Oc_goes = bad_H2Oc_goes | D(212, :) > 29.8; %h2o conc are saturated at 30 mmol/mol 
    bad_H2Oc_goes = bad_H2Oc_goes | D(204, :) == 0; %no variabilty 
    bad_H2Oc_goes = bad_H2Oc_goes | D(212, :) < 0.1; %no variabilty 
    bad_H2Oc_goes = bad_H2Oc_goes | time < 465;
    bad_H2Oc_goes = bad_H2Oc_goes | time > 1005.78 & time < 1005.82;
    bad_H2Oc_goes = bad_H2Oc_goes | time > 2092.8 & time < 2096.745; 
    bad_H2Oc_goes = bad_H2Oc_goes | time > 2121.38 & time < 2125.4; 
    
    D(212, bad_H2Oc_goes) = NaN;

    %convert to mixing ratio
    mixing_goes=D(212, :)./(1-(D(212, :)/1000)); %mol h2o/mol dry air
    D(212, :)= mixing_goes; %mmol h2o/mol dry air

    %----------------------------------------------------------------------
    %build a consensus H2O mixing ratio [mmol h20/mol of dry air]
    H2O=D(100,:);%HMP 
    bad=isnan(H2O(1,:));
    H2O(1,bad)=D(226,bad); %GOES HMP
    bad=isnan(H2O(1,:));
    H2O(1,bad)=D(42,bad);%fast
    bad=isnan(H2O(1,:));
    H2O(1,bad)=D(86,bad);%dl
    bad=isnan(H2O(1,:));
    H2O(1,bad)=D(212,bad);%goes

    %a fill H2O - h2o has a small effect on density ~ few percent - so we
    %should fill in h2o for density corrections 
    medianh2o=nanmedian(H2O(1,:));
    h2ofill=medianh2o.*o;

    bad=isnan(H2O(1,:));
    H2O(1,bad)=h2ofill(1,bad);

    %----------------------------------------------------------------------
    %convert sonic temp over to actual temperature and fill
    %----------------------------------------------------------------------
    %water mixing ratio to mol fraction [mmol h20/mol of wet air]
    h2o_mol_frac = H2O./(1+(H2O/1000));
    %----------------------------------------------------------------------
    %Kaimal and Gaynor - 1991
    Tkact_mat   = (D(8,:)+Tc)  ./ (1 + 0.00032 * h2o_mol_frac); %[K]
    Tkact_dl   = (D(84,:)+Tc)  ./ (1 + 0.00032 * h2o_mol_frac); %[K]
    Tkact_goes   = (D(210,:)+Tc)  ./ (1 + 0.00032 * h2o_mol_frac); %[K]
    
    Tkact_mat   = Tkact_mat-Tc; %[C]
    Tkact_dl   = Tkact_dl-Tc; %[C]
    Tkact_goes   = Tkact_goes-Tc; %[C]
    %----------------------------------------------------------------------
    %build a consensus Temp
    TEMP=D(97,:);%HMP Temp
    bad=isnan(TEMP(1,:));
    TEMP(1,bad)=D(223,bad); %GOES TEMP
    bad=isnan(TEMP(1,:));
    TEMP(1,bad)=Tkact_mat(1,bad);%fast
    bad=isnan(TEMP(1,:));
    TEMP(1,bad)=Tkact_dl(1,bad);%dl
    bad=isnan(TEMP(1,:));
    TEMP(1,bad)=Tkact_goes(1,bad);%goes
    %----------------------------------------------------------------------

    %fast calculated a bad dry density if the water vapor was wrong or temp was
    %wrong - The error was probably small but we can adjust this with the ratio
    %of best_dry_density/dry_density_that_was_used ==> Rho Calc in combined ./ D(48,:)

    %When rho_fix ==1, we should apply the density correction
    rho_fix = bad_Ts_mat | bad_H2Oc_mat;
    %~49% of the observations will be adjusted
    %----------------------------------------------------------------------
    %winds
    %----------------------------------------------------------------------
    %verticle windspeed (w) is incorrect for all time < 740 for all data sources

    %component wind problems - variance is 0 is impossible
    %u_wind
    u_mat_var=D(9,:);
    Umat_isnan=isnan(D(5,:)) | isempty(D(5,:)) | ~isfinite(D(5,:));
    ubad_wind_mat=u_mat_var(1,:)==0;
    ubad_wind_mat = ubad_wind_mat | Umat_isnan(1,:);
    D(5,ubad_wind_mat)=NaN;
    
    %v_wind
    v_mat_var=D(10,:);
    Vmat_isnan=isnan(D(6,:)) | isempty(D(6,:)) | ~isfinite(D(6,:));
    vbad_wind_mat = v_mat_var(1,:)==0;
    vbad_wind_mat = vbad_wind_mat | Vmat_isnan(1,:)==1;
    D(6,vbad_wind_mat)=NaN;


    w_mat_var=D(11,:);
    Wmat_isnan=isnan(D(7,:)) | isempty(D(7,:)) | ~isfinite(D(7,:));
    wbad_wind_mat = w_mat_var(1,:)==0;
    wbad_wind_mat = wbad_wind_mat | Wmat_isnan(1,:)==1;
    D(7,wbad_wind_mat)=NaN;
    
    %clean dl
    %DATA LOGGER
    u_dl_var=D(66,:);%'Ux variance
    Udl_isnan=isnan(D(81,:)) | isempty(D(81,:)) | ~isfinite(D(81,:));
    ubad_wind_dl=u_dl_var(1,:)==0;
    ubad_wind_dl = ubad_wind_dl | Udl_isnan(1,:)==1;
    D(81,ubad_wind_dl)=NaN;

    v_dl_var=D(71,:);%'Uy variance
    Vdl_isnan=isnan(D(82,:))==1 | isempty(D(82,:)) ==1 | ~isfinite(D(82,:));
    vbad_wind_dl = v_dl_var(1,:)==0;
    vbad_wind_dl = vbad_wind_dl | Vdl_isnan(1,:)==1;
    D(82,vbad_wind_dl)=NaN;


    w_dl_var=D(60,:);% Uz variance
    Wdl_isnan=isnan(D(83,:)) | isempty(D(83,:)) | ~isfinite(D(83,:));
    wbad_wind_dl = w_dl_var(1,:)==0;
    wbad_wind_dl = wbad_wind_dl | Wdl_isnan(1,:)==1;
    D(83,wbad_wind_dl)=NaN;

    %clean goes
    u_goes_var = D(192,:);% Ux variance
    Ugoes_isnan=isnan(D(207,:)) | isempty(D(207,:)) | ~isfinite(D(207,:));
    ubad_wind_goes = u_goes_var(1,:) == 0;
    ubad_wind_goes = ubad_wind_goes | Ugoes_isnan ==1;
    D(207,ubad_wind_goes) = NaN;% Ux
    
    v_goes_var = D(197,:);% Uy variance 
    Vgoes_isnan=isnan(D(208,:)) | isempty(D(208,:)) | ~isfinite(D(208,:));
    vbad_wind_goes = v_goes_var==0;
    vbad_wind_goes = vbad_wind_goes | Vgoes_isnan ==1;
    D(208,vbad_wind_goes) = NaN;% Uy
    
    w_goes_var = D(186,:);% Uz variance
    Wgoes_isnan=isnan(D(209,:)) | isempty(D(209,:)) | ~isfinite(D(209,:));
    wbad_wind_goes = w_goes_var==0;
    wbad_wind_goes = wbad_wind_goes | Wgoes_isnan ==1;
    D(209,wbad_wind_goes) = NaN;% Uy

    %fill winds
    %u wind
    u=D(5,:);%u wind
    bad=isnan(u(1,:));
    u(1,bad)=D(81,bad); %u wind
    bad=isnan(u(1,:));
    u(1,bad)=D(207,bad);%u wind

    %v wind
    v=D(6,:);%v wind
    bad=isnan(v(1,:));
    v(1,bad)=D(82,bad); %v wind
    bad=isnan(v(1,:));
    v(1,bad)=D(208,bad);%v wind

    %w wind
    w=D(7,:);%w wind
    bad=isnan(w(1,:));
    w(1,bad)=D(83,bad); %w wind
    bad=isnan(w(1,:));
    w(1,bad)=D(209,bad);%w wind

    v_wind=[u;v;w];

    %----------------------------------------------------------------------
    %% Wind Direction
    %----------------------------------------------------------------------
    %mat wind dir is wrong when component winds are wrong
    bad_wind_mat = ubad_wind_mat | vbad_wind_mat | wbad_wind_mat;
    D(22,bad_wind_mat)=NaN;

    %----------------------------------------------------------------------
    %% Ustar
    %----------------------------------------------------------------------
    %clean FAST
    bad = D(9,:) == 0; % u variance == 0 
    bad = bad | isnan(D(9,:)) == 1; %missing u
    bad = bad | D(10,:) == 0; % u variance == 0 
    bad = bad | isnan(D(10,:)) == 1; %missing u
    bad = bad | D(11,:) == 0; %uz variability = 0
    bad = bad | isnan(D(11,:)) == 1; %missing uz
    D(21,bad)=NaN;

    %clean data logger
    %bad = bad_wind_dl;
    bad = D(60,:) == 0;
    bad = bad | isnan(D(60,:)); % is there uz uz data
    bad = bad | D(66,:) == 0; % sensor not working
    bad = bad | isnan(D(66,:)); % is there ux ux data
    bad = bad | time > 1783.915 & time < 1783.92;
    D(61,bad)=NaN;

    bad = D(60,:) == 0;
    bad = bad | isnan(D(60,:)); % is there uz uz data
    bad = bad | D(71,:) == 0; % sensor not working
    bad = bad | isnan(D(71,:)); % is there uy uy data
    bad = bad | time > 1783.915 & time < 1783.92;
    D(62,bad)=NaN;

    %clean GOES
    bad = D(186,:) == 0; % uz uz
    bad = bad | isnan(D(186,:)); % is there uz uz data
    bad = bad | D(192,:) == 0; % sensor not working
    bad = bad | isnan(D(192,:)); % is there ux ux data
    bad = bad | time > 1783.915 & time < 1783.92;
    D(187,bad)=NaN;

    bad = D(186,:) == 0; % uz uz
    bad = bad | isnan(D(186,:)); % is there uz uz data
    bad = bad | D(197,:) == 0; % sensor not working
    bad = bad | isnan(D(197,:)); % is there uy uy data
    bad = bad | time > 1783.915 & time < 1783.92;
    D(188,bad)=NaN;

    %----------------------------------------------------------------------
    %% Sensible heat
    %----------------------------------------------------------------------
    %clean fast - this is despiked in fastflux
    bad = D(12,:) == 0; % Ts var == 0 
    bad = bad | isnan(D(12,:))==1; %missing Ts
    bad = bad | D(11,:) == 0; %w variability = 0
    bad = bad | isnan(D(11,:))==1; %missing w
    D(29,bad)=NaN;

    bad = D(29,:) < -300;
    bad = bad | D(29,:) > 1000;
    D(29,bad)=NaN;

    %clean data logger
    %need a bad covariance term
    bad = D(60,:) == 0;
    bad = bad | isnan(D(60,:))==1; % is there uz uz data
    bad = bad | isnan(D(80,:))==1; % is there uy uy data
    bad = bad | D(80,:) == 0; % sensor not working
    bad = bad | D(65,:) < -15; % sensor not working
    bad = bad | D(65,:) > 5; % sensor not working
    D(65,bad)=NaN;

    %clean GOES 
    bad = D(186,:) == 0; % uz uz
    bad = bad | isnan(D(186,:))==1; % is there uz uz data
    bad = bad | isnan(D(206,:))==1; % is there ux ux data
    bad = bad | D(206,:) == 0; % sensor not working
    bad = bad | D(191,:) > 4; % sensor not working
    D(191, bad) = NaN;

    %----------------------------------------------------------------------
    %% CO2 concentration
    %----------------------------------------------------------------------
    %clean fast
    %D(35,bad_fastCO2)=NaN;

    %clean data logger
    bad_dlCO2 = D(85,:) < 250;
    D(85, bad_dlCO2) = NaN;

    D(85,:)=D(85,:)./(1-(h2o_mol_frac(1,:)/1000)); %co2/mol dry air
    
    %clean GOES
    bad_goesCO2 = D(211,:) < 250;
    D(211, bad_goesCO2) = NaN;
    
    D(211,:)=D(211,:)./(1-(h2o_mol_frac(1,:)/1000)); %co2/mol dry air
    %----------------------------------------------------------------------
    %% FCO2
    %----------------------------------------------------------------------
    %clean fast
    bad = D(11,:) == 0; %w variability = 0
    bad = bad | isnan(D(11,:))==1; %missing w
    bad = bad | D(36,:) == 0; %co2 variability = 0
    bad = bad | isnan(D(36,:))==1; %missing co2
    D(46, bad) = NaN;

    %out of reasonable range range
    bad = D(46,:) < -32;
    bad = bad | D(46,:) > 20;
    D(46, bad) = NaN;

    %clean data logger FCO2 by taking out uz_co2 covariance
    bad_fco2_dl = isnan(D(60,:)); % no uz uz
    bad_fco2_dl = bad_fco2_dl | D(60,:)==0;%no wind variability
    bad_fco2_dl = bad_fco2_dl | D(75,:)==0; % co2 variability
    bad_fco2_dl = bad_fco2_dl | D(75,:)> 3000; % co2 variability
    bad_fco2_dl = bad_fco2_dl | isnan(D(75,:)); % no data
    bad_fco2_dl = bad_fco2_dl | isnan(D(63,:)); % no data
    D(63, bad_fco2_dl) = NaN;
    
    %clean GOES
    bad_fco2_goes = isnan(D(186,:)); % sonic is not working no uz uz
    bad_fco2_goes = bad_fco2_goes | D(186,:)==0;%no wind variability
    bad_fco2_goes = bad_fco2_goes | D(201,:)<=0; % co2 variability
    bad_fco2_goes = bad_fco2_goes | D(75,:)> 3000; % fast co2 var catches GOES varibility
    bad_fco2_goes = bad_fco2_goes | isnan(D(201,:)); % sonic is not working no uz uz
    bad_fco2_goes = bad_fco2_goes | isnan(D(189,:)); % no data
    D(189, bad_fco2_goes) = NaN;

    %----------------------------------------------------------------------
    %% Adjust water fluxes using IRGA H2O variabilty compared with HMP variability
    %NOTE: IRGA H2O variabilty seems high compared with HMP variability.
    %Therefore, we over estimate FH2O.  This will reduce the water flux.

    %Make sure the water vapor concentrations are cleaned first.
    %----------------------------------------------------------------------
    %IRGA calibrations
    disp('need to include irga swap calibrations');
    Cnums=datenum(IRGAchanged(:,3), IRGAchanged(:,1), IRGAchanged(:,2));
    Cnumsince2009=Cnums-datenum(2008,1,0); % this is for the P301

    firstday = min(time);
    lastday = max(time);

    length_Cnumsince2009=length(Cnumsince2009(:,1));
    Cnumsince2009(1,1)=firstday-1;
    Cnumsince2009(length_Cnumsince2009+1,1)=lastday;
    
    %--------------------------%
    %[I added this - Aaron]
    %--------------------------%
    length_Cnumsince2009 = length_Cnumsince2009+1;
    
    % serviced figure
    f1=(Cnumsince2009./Cnumsince2009)*40;
    f2=(Cnumsince2009./Cnumsince2009)*-5;
    f3=(Cnumsince2009./Cnumsince2009)*40;
    F=[f1,f2,f3]';
    [m,n]=size(F);
    l=m*n;
    F=reshape(F,l,1);
    C09=[Cnumsince2009; Cnumsince2009; Cnumsince2009];
    C09=sort(C09);

    figure(2)
    plot(time, D(100,:), 'k.');
    hold on
    plot(time, D(42,:), 'b.');
    hold on
    plot(C09, F, 'r-');
    legend('black=hmp', 'blue=irga', 'irga calibration');

    Pmat=[];
    Pdl=[];
    Pgoes=[];


    for i=2:length_Cnumsince2009
        %proc data
        good_mat = time(1,:) < Cnumsince2009(i)-2 & time(1,:) > Cnumsince2009(i-1)+2 & isfinite(D(42,:)) & isfinite(D(100,:));
        irga=D(42,good_mat);
        hmp=D(100,good_mat);
        p = polyfit(hmp,irga,1) %p(1) is slope; p(2) is intercept
        c_tm = time(1,:) <= Cnumsince2009(i) & time(1,:) > Cnumsince2009(i-1); %adjustment period

        figure(1)
        plot(hmp,irga, 'k.')
        hold on
        plot(hmp, p(1,1)*hmp+p(1,2), 'r-')
        xlabel('hmp h2o')
        ylabel('dl irga h2o')
        pause 
        close(figure(1))
    
        %record so you can look later
        tt= median(time(1,good_mat)); 
        Pmat=[Pmat; tt, p];

        %The offset does not make sense here. Adjust water fluxes by dividing by the slope.
        outrageous_p = isnan(p)==1 | ~isfinite(p) | p(1) > 1.47 | p(1) < 0.8 | length(irga) < 30;
        
        if outrageous_p == 1
            disp('can not adjust h2o fluxes 1')
            tt
        else
            D(47,c_tm)=D(47,c_tm)./p(1); %Fh20
            D(30,c_tm)=D(30,c_tm)./p(1); %LH
        end
      
    
        %data logger
        good_dl = time(1,:) < Cnumsince2009(i)-2 & time(1,:) > Cnumsince2009(i-1)+2 & isfinite(D(86,:))==1 & isfinite(D(100,:))==1;
        irga=D(86,good_dl);
        hmp=D(100,good_dl);
        p = polyfit(hmp,irga,1) %p(1) is slope; p(2) is intercept
        c_tm = time(1,:) <= Cnumsince2009(i) & time(1,:) > Cnumsince2009(i-1); %adjustment period

        
        figure(1)
        plot(hmp,irga, 'k.')
        hold on
        plot(hmp, p(1,1)*hmp+p(1,2), 'r-')
        xlabel('hmp h2o')
        ylabel('dl irga h2o')
        pause 
        close(figure(1))
        
        
        %record so you can look later
        tt= median(time(1,good_dl)); 
        Pdl=[Pdl; tt, p];

        outrageous_p = isnan(p)==1 | ~isfinite(p) | p(1) > 1.47 | p(1) < 0.8 | length(irga) < 30;

        if outrageous_p == 1
            disp('can not adjust h2o fluxes 2')
            tt
        else
        %The offset does not make sense here. Adjust water fluxes by dividing by the slope.
            D(64,c_tm)=D(64,c_tm)./p(1); %cov h2o dl
        end
    
    
        %goes
        good_goes = time(1,:) < Cnumsince2009(i)-2 & time(1,:) > Cnumsince2009(i-1)+2 & isfinite(D(212,:))==1 & isfinite(D(226,:))==1;
        irga=D(212,good_goes);
        hmp=D(226,good_goes);
        p = polyfit(hmp,irga,1) %p(1) is slope; p(2) is intercept
        c_tm = time(1,:) <= Cnumsince2009(i) & time(1,:) > Cnumsince2009(i-1); %adjustment period

        %record so you can look later
        tt= median(time(1,good_goes)); 
        Pgoes=[Pgoes; tt, p];

        %The offset does not make sense here. Adjust water fluxes by dividing by the slope.
        outrageous_p = isnan(p) | ~isfinite(p) | p(1) > 1.47 | p(1) < 0.75 | length(irga) < 30;
        
        figure(1)
        plot(hmp,irga, 'k.')
        hold on
        plot(hmp, p(1,1)*hmp+p(1,2), 'r-')
        xlabel('hmp h2o')
        ylabel('dl irga h2o')
        pause 
        close(figure(1))
        
        if outrageous_p == 1
            disp('can not adjust h2o fluxes 3')
            tt
        else
            D(190,c_tm)=D(190,c_tm)./p(1); %cov h2o goes

        end 
    end
    %----------------------------------------------------------------------
    %% FH2O
    %----------------------------------------------------------------------
    %clean fast
    bad_mat_et = bad_H2Oc_mat; 
    bad_mat_et = bad_mat_et | D(43,:) == 0; %
    bad_mat_et = bad_mat_et | isnan(D(43,:)) == 1; %
    bad_mat_et = bad_mat_et | D(11,:) == 0; %uz variability = 0
    bad_mat_et = bad_mat_et | isnan(D(11,:)) == 1; %missing uz
    bad_mat_et = bad_mat_et | D(47,:) < -10;
    bad_mat_et = bad_mat_et | D(47,:) > 20;

    D(47, bad_mat_et) = NaN;
    
    %clean data logger
    bad_dl_et = bad_H2Oc_dl;
    bad_dl_et = bad_dl_et | isnan(D(60,:))==1; % sonic is not working no uz uz
    bad_dl_et = bad_dl_et | D(60,:)==0;%no wind variability
    bad_dl_et = bad_dl_et | D(78,:)==0; % h2o variability
    bad_dl_et = bad_dl_et | isnan(D(78,:))==1; % sonic is not working no uz uz
    bad_dl_et = bad_dl_et | isnan(D(64,:)) == 1; % no data

    D(64, bad_dl_et) = NaN;

    %clean GOES
    bad_goes_et = bad_H2Oc_goes;
    bad_goes_et = bad_goes_et | isnan(D(186,:))==1; % sonic is not working no uz uz
    bad_goes_et = bad_goes_et | D(186,:)==0;%no wind variability
    bad_goes_et = bad_goes_et | D(204,:)==0; % h2o variability
    bad_goes_et = bad_goes_et | isnan(D(204,:))==1; % 
    bad_goes_et = bad_goes_et | isnan(D(190,:)) == 1; % no data

    D(190, bad_goes_et) = NaN;

    %----------------------------------------------------------------------
    %% Latent heat
    %----------------------------------------------------------------------
    %fast
    %bad_le_mat = bad_mat_et | D(30,:) == 0;
    bad_le_mat = bad_mat_et;
    D(30,bad_le_mat)=NaN;

    %----------------------------------------------------------------------
    %% Rn
    %----------------------------------------------------------------------
    %clean data logger
    bad = time < 302.1;
    D(87,bad) = -1*D(87,bad);

    %clean goes
    %----------------------------------------------------------------------
    %% Incoming solar radiation - pyranometer = SOLAR_IN
    %----------------------------------------------------------------------
    %clean data logger

    
    %clean goes


    %----------------------------------------------------------------------
    %% Outgoing solar radiation - pyranometer = SOLAR_OUT
    %----------------------------------------------------------------------
    %clean data logger


    %clean goes
    bad = D(215,:) > 1400; 
    bad = bad | D(215,:) < -500; 
    D(215,bad) = NaN;

    %----------------------------------------------------------------------
    %% PAR_In
    %----------------------------------------------------------------------
    %clean data logger
    focus = time < 302.06;
    D(90,focus) = -1 * D(90,focus);
    %clean goes
    bad = D(216,:) > 3000;
    bad = bad | D(216,:) < -500;
    D(216,bad) = NaN;

    %----------------------------------------------------------------------
    %% PAR_Out
    %----------------------------------------------------------------------
    %clean data logger
    focus = time < 302.06;
    D(91,focus) = -1 * D(91,focus);
% 
    %clean goes
%     bad = zeros(size(D(90,:)));
%     D(217,bad) = NaN;

    %----------------------------------------------------------------------
    %% T_107
    %----------------------------------------------------------------------
    %dl
    bad = D(102,:) < -30;
    bad = bad | time > 1226.7 & time < 1226.84;
    D(102,bad) = NaN;

    %goes
    bad = D(228,:) < -30;
    bad = bad | time > 464.9 & time < 1005.74;
    
    focus = time > 240 & time < 280;
    focus = focus | time > 412 & time < 416;
    focus = focus | time > 1059 & time < 1060;
    focus = focus | time > 1095 & time < 1096;
    focus = focus | time > 1776 & time < 1777;
    bad = bad | focus & D(228,:) == 0;
    D(228,bad) = NaN;

    %----------------------------------------------------------------------
    %% Rain
    %----------------------------------------------------------------------
    %data logger
    bad=D(103,:)>5000;
    D(103,bad) = NaN;

    %goes
    bad=D(229,:)>200;
    D(229,bad) = NaN;
    %----------------------------------------------------------------------
    %% NDVI
    %----------------------------------------------------------------------
    bad = D(96,:) > 1;
    D(96,bad)=NaN;

    bad = D(222,:) > 1;
    D(222,bad)=NaN;

    %----------------------------------------------------------------------
    %Soil T
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    %Clean Soil Moisture
    %----------------------------------------------------------------------

    bad = D(282,:) > 0.6;
    D(282,bad) = NaN;
    
    %----------------------------------------------------------------------
    %Fuel_M
    %----------------------------------------------------------------------


    %----------------------------------------------------------------------
    %LWS
    %----------------------------------------------------------------------


    %----------------------------------------------------------------------
    %Clean Matric potential sensors 
    %----------------------------------------------------------------------


elseif a==2

    %Mop-up

    TEMP =[];
    H2O =[];
    rho_fix =[];
    v_wind =[];    

    %Sensible Heat
    %dl
    bad = D(339,:)<-300;
    bad = bad | D(339,:)> 1000;
    D(339,bad)=NaN;

    %goes
    bad = D(340,:)<-300;
    bad = bad | D(340,:)> 1000;
    D(340, bad) = NaN;

    %Fh2o
    %dl
    bad_dl_et = D(347,:) < -10;
    bad_dl_et = bad_dl_et | D(347,:) > 20;
    D(347, bad_dl_et) = NaN;

    %goes
    bad_goes_et = D(348,:) < -10;
    bad_goes_et = bad_goes_et | D(348,:) > 20;
    D(348, bad_goes_et) = NaN;

    %latent heat
    bad_le_dl = bad_dl_et;
    D(343,bad_le_dl)=NaN;

    %clean GOES 
    bad_le_goes = bad_goes_et;
    D(344,bad_le_goes)=NaN;

    %FCO2
    %dl
    bad = D(345,:) < -32;
    bad = bad | D(345,:) > 20;
    D(345, bad) = NaN;

    %goes
    bad = D(346,:) < -32;
    bad = bad | D(346,:) > 20;
    D(346, bad) = NaN;

    %Ustar
    %clean dl
    bad = D(322,:) > 3;
    D(322,bad)=NaN;

    %clean GOES
    bad = D(323,:) > 3;
    D(323,bad)=NaN;


end