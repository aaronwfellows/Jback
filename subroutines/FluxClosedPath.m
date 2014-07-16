function [CO2,H2O,RHOM,FCO2,FH2O,HLATENT,HSENSIBLE,GAINS, Ntemp, Nco2, Nh2o]=FluxClosedPath(w,Ts,SONDIAG,Xc,Xw,IRGADIAG,PFlux,P_estimate)

var_defs();
global Mw_da Tc R_mol CpH2Oc CpAirc
P = P_estimate;
qplot='n';
% May 11, 2012 - turned on maximizes covariance between w and T so relies
% solely on the shaw model to get time constant and delay (amsm)

% Oct 19, 2004 - no longer maximizes covariance between w and T so relies
% solely on the shaw model to get time constant and delay (amsm)

% Feb 25, 2003 - the number of samples to shift when looking for the
% maximum covariance modified to seconds instead of samples to acct
% for changes insampling frequency for this dataset. Also, changed
% the number of samples to shift for temperature (heat flux) to
% zero, since these should be lined up

% aug 15, 2002 - we don't filter the vertical velocity for making
% the covariance with degraded temperature and uncorrected co2 or
% h2o signal, because the effect should cancel when calculating the
% gain because the vertical velocity is both in the numerator and
% denominator (see Berger paper)

% nov 12, 2001 - calcuating air density here now.  The correct
% densisty to use for the flux is the dry air density at the
% sampling inlet, not the density of air at the irga.  Previous
% recent runs used the irga cell density and were giving fluxes
% 15 percent too low, which makes sense because the pressure in
% the cell is about 85 percent of atmospheric.  In calculating
% ambient air density use the nominal value for pressure 101.3
% kPa, since we do not have an extra pressure sensor, and since
% we need the closed path to function independent of the open path


%procedure

%1) shift arrays Co2, H2o
%2) filter temperature
%3) maximize covariance with w and T
%4) maximize covariance with w and Tf
%5) calculate the gain factor from unfiltered/filtered ratio, Gco2, Gh2o
%6) maximize covariance with w and co2, w and h2o
%7) adjust using the gain factors Gco2, Gh2o

%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Air Density
%%%%%%%%%%%%%%%%%%%%%%%%

% calculate  air molar density (moles moist air/m^3 moist air)
% rho = P * (1-Xw) / (R_u * Td )
% P (kpA)
% Xw = mole fraction h2o - convert with factor of 1000
% Td - K

% - Convert temperature if sent in Celcius
%try
if median(Ts)<100
    Ts=Ts+Tc;
end
%catch 
%    disp('line 61 Flux..');
%end

ixw = despike(Xw,8,-15,45,'Xw');
%ixc = despike(Xc,8,200,600,'Xc');
its = despike(Ts,8,223,353,'Ts');
%disp('despikeing in FCPC');
%pause
Td    = Ts./(ones(size(Ts))+0.32*Xw/1000);

rho_a =  1e3/R_mol * P * (1-mean(1e-3* Xw(find(its&ixw)))) / mean(Td(find(ixw&its))); %#ok<*FNDSB>
rho_w =  1e3/R_mol * P * (  mean(1e-3* Xw(find(its&ixw)))) / mean(Td(find(ixw&its)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIFIC HEAT CAPACITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dry air Cp - not a moist air Cp
Cpv = CpH2Oc;     %  J/Kg K

%Cp = CpAirc ;  %mean(CpAirc*(ones(size(Xw))-Xw) + Cpv.*Xw)
%Cp = mean(CpAirc*(ones(size(Xw))-Xw) + Cpv.*Xw);
Cp = mean(CpAirc*(ones(size(Xw))-(1e-3*Xw)) + Cpv.*(1e-3*Xw));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LATENT HEAT OF VAPORIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     Lv    - Latent heat of vaporization (J/Kg)
%                in J/kg from Stull p 641

% Feb 8 2001 - changed the index to iok, since the latent heat
% of vaporization is based on the dried temperature which requires
% the moisture

Lv = (2.501-0.00237*(mean(Td(find(its)))-Tc))*10^6;

%%%%%%%%%%%%%%%%%%%%
% Shift Arrays
%%%%%%%%%%%%%%%%%%%%

len   = length(Xc);
delay = PFlux.sf*floor((PFlux.delayco2+PFlux.delayh2o)/2);


% cut off end of  signals

Xc=Xc(1,delay+1:len);
Xw=Xw(1,delay+1:len);

IRGADIAG=IRGADIAG(:,delay+1:len);

% cut off beginning of undelayed sensors

w   = w(1:len-delay);
SONDIAG  = SONDIAG(:,1:len-delay);
Ts       = Ts(1:len-delay);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SENSIBLE HEAT FLUX (W/m^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% have to use the max covar because water vapor goes into drying
% Tsonic (think about this some).

%we want an estiamte of the delay from covarmax - awf
[wTd,Nwt] = covarmax( w(find(SONDIAG&IRGADIAG(1,:))) , Td(find(SONDIAG&IRGADIAG(1,:))),10*PFlux.sf );
Ntemp=Nwt;

%[wTd,Nwt] = covarmax( w(find(SONDIAG&IRGADIAG(1,:))) , Td(find(SONDIAG&IRGADIAG(1,:))),2*PFlux.sf );
HSENSIBLE =  Mw_da/1000*rho_a*Cp*wTd; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter temperature and LI7000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% account for 1 Hz sampling frequency, such that the filter
% normalized frequency has to be set to 1. This means we must be
% losing some iformation (think about it more)

wn=1/PFlux.tauco2/(PFlux.sf/2);
if wn>1; wn=0.99; end


[bco2,aco2]=butter(1,wn,'low');

wn=1/PFlux.tauh2o/(PFlux.sf/2);
if wn>1; wn=0.99; end
[bh2o,ah2o]=butter(1,wn,'low');

% first interpolate for bad points (added nov 13, 2001)

iok=find(SONDIAG);

Ts=interp1(iok,Ts(iok),1:length(Ts),'linear','extrap');
w=interp1(iok,w(iok),1:length(Ts),'linear','extrap');

Tsco2f = filter(bco2,aco2,Ts);
Tsh2of = filter(bh2o,ah2o,Ts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Truncate all time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Skip the first minute of the time series to allow
% for spinup of the filters - use h2o since it has a longer
% spinup

Ntrunc=floor(60*PFlux.sf) ;

len=length(Xc);

% cut off end of  signals

Xc=Xc(Ntrunc:len);
Xw=Xw(Ntrunc:len);

IRGADIAG=IRGADIAG(:,Ntrunc:len);

w   = w(:,Ntrunc:len);
SONDIAG  = SONDIAG(:,Ntrunc:len);
Ts       = Ts(Ntrunc:len);

Tsco2f   = Tsco2f(Ntrunc:len);
Tsh2of   = Tsh2of(Ntrunc:len);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET MAXIMUM COVARIANCES, UNFILTERED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%May 11, 2013 - modified to use covarmax and not cov - awf

% feb 25, 2003 - modify for no shift for temperuatre, and for just 2 second shift for
% co2 and water vapor

tmp = cov( w(find(SONDIAG)),Ts(find(SONDIAG))) ;
wTs = tmp(1,2);


%[wXc,NwXc] = covarmax(w(find(SONDIAG&IRGADIAG(1,:))),Xc(find(SONDIAG&IRGADIAG(1,:))), 2*PFlux.sf) ;
%[wXw,NwXw] = covarmax(w(find(SONDIAG&IRGADIAG(1,:))),Xw(find(SONDIAG&IRGADIAG(1,:))), 2*PFlux.sf) ;

%we want an estiamte of the delay from covarmax - awf
[wXc,NwXc] = covarmax(w(find(SONDIAG&IRGADIAG(1,:))),Xc(find(SONDIAG&IRGADIAG(1,:))), 10*PFlux.sf) ;
[wXw,NwXw] = covarmax(w(find(SONDIAG&IRGADIAG(1,:))),Xw(find(SONDIAG&IRGADIAG(1,:))), 10*PFlux.sf) ;
Nco2=NwXc;
Nh2o=NwXw;

%tmp=cov(w(find(SONDIAG&IRGADIAG(1,:))),Xc(find(SONDIAG&IRGADIAG(1,:)))) ;
%wXc=tmp(1,2);
%NwXc=0;

%tmp=cov(w(find(SONDIAG&IRGADIAG(1,:))),Xw(find(SONDIAG&IRGADIAG(1,:)))) ;
%wXw=tmp(1,2);
%NwXw=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET MAXIMUM COVARIANCES, FILTERED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%should the buttered Ts be coarmaxed??? - awf

% feb 25, 2002, removed maximizing covariance

[wTsco2f,NwTsco2]   = covarmax(w(find(SONDIAG)),Tsco2f(find(SONDIAG)), 2*PFlux.sf) ;
[wTsh2of,NwTsh2o]   = covarmax(w(find(SONDIAG)),Tsh2of(find(SONDIAG)), 2*PFlux.sf) ;

%tmp = cov(w(find(SONDIAG)),Tsco2f(find(SONDIAG))) ;

%wTsco2f = tmp(1,2);
%NwTsco2 = 0;

%tmp = cov(w(find(SONDIAG)),Tsh2of(find(SONDIAG))) ;

%wTsh2of = tmp(1,2);
%NwTsh2o = 0;


%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE GAIN FACTORS
%%%%%%%%%%%%%%%%%%%%%%%%

Gco2  = abs(wTs/wTsco2f);
Gh2o  = abs(wTs/wTsh2of);

%disp('pause');
%pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADJUST CLOSED PATH FLUXES BASED ON GAIN FACTORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wXc_adj =   wXc * Gco2;
wXw_adj =   wXw * Gh2o;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WATER VAPOR FLUX AND LATENT HEAT FLUX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  E = mean( rho_a ) * <w'q'>    [millimoles h2o/m2/s]
%
%  H_l = rho_a * Lv * cov(w,q)   [Watts/m^2]
%
%          H_l   - Latent heat flux (W/m^2)
%          rho_a - dry air density (moles dry air/m^3 moist air)

Ecorr = rho_a*wXw_adj;
HLcorr = 18/1000*Lv*Ecorr/1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CO2 FLUX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  F = mean( rho_a ) * <w'x'>    (umol/m^2/s)
%
%  rho_a - dry air density (moles dry air/m^3 moist air)
%      x - mole fraction of co2 (umol/mol)
%
%  wco2d_adj - covariance(w,co2) -->  co2 in (umol/mol dry air)
%                  including high frequency correction, Webb Corr

Fco2 = rho_a*wXc_adj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Format OUput
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=Xc(1,find(IRGADIAG(1,:)));
CO2 = [min(x);  max(x); median(x);  mean(x); std(x); skewness(x); kurtosis(x) ];

x=Xw(1,find(IRGADIAG(1,:)));
H2O = [min(x);  max(x); median(x);  mean(x); std(x); skewness(x); kurtosis(x) ];

FCO2 = Fco2;
FH2O = Ecorr;

% replaced third componetn with tower top density estimate
% 11/15/2001 for testing

RHOM      = [rho_a;rho_w];
HLATENT   = HLcorr;

GAINS = [Gco2; Gh2o];


%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING

if findstr(qplot,'y')

    iok=find(SONDIAG&IRGADIAG(1,:));

    size(iok)
    size(w)
    size(Xc)
    size(Xw)

    max(iok)

    %Pw = hspec_ss(w(iok),0,PFlux.sf);
    Pt = hspec_ss(Ts(iok),0,PFlux.sf);
    Ptfc = hspec_ss(Tsco2f(iok),0,PFlux.sf);
    Ptfh = hspec_ss(Tsh2of(iok),0,PFlux.sf);

    % co2

    PXc  = hspec_ss(Xc(1,iok),0,PFlux.sf);
    PXw  = hspec_ss(Xw(1,iok),0,PFlux.sf);

    % plot co2 power spectra with unfiltered and filtered temperature

    figure(7);
    clf
    plot([Ts;Tsco2f;Tsh2of]')

    figure(1);
    subplot(221)
    f=Pt(:,1);
    loglog(f,f.^(5/3).*Pt(:,2)/var(Ts),f,f.^(5/3).*Ptfc(:,2)/var(Tsco2f),f,f.^(5/3).*PXc(:,2)/var(Xc(iok)),[1e-4 10],[.01 .01],'k')
    set(gca,'xlim',[1/1800 1]);
    title('Temperature (meas, blue), (t with co2 filter, g), (co2 closed path, r)');
    ylabel('f^{(5/3)}*P','fontsize',12,'fontweight','bold');


    subplot(222)
    f=Pt(:,1);
    loglog(f,f.^(5/3).*Pt(:,2)/var(Ts),f,f.^(5/3).*Ptfh(:,2)/var(Tsh2of),f,f.^(5/3).*PXw(:,2)/var(Xw(iok)),[1e-4 10],[.01 .01],'k')
    set(gca,'xlim',[1/1800 1]);
    title('Temperature (meas, blue), (temp with h2o filter, g), (h2o closed path, r)');
    ylabel('f^{(5/3)}*P','fontsize',12,'fontweight','bold');


    %%%% CO2 Spectra and Cospectra
    %[FwXc,PwXc,sig2c] = cross_spectra(w(iok),Xc(1,iok),PFlux.sf,NwXc);
    %[FwXw,PwXw,sig2h] = cross_spectra(w(iok),Xw(1,iok),PFlux.sf,NwXw);
    [FwXc,PwXc,Wonky] = cross_spectra(w(iok),Xc(1,iok),PFlux.sf,NwXc);
    [FwXw,PwXw,Wonky] = cross_spectra(w(iok),Xw(1,iok),PFlux.sf,NwXw);

    OgiveXc=ogive(FwXc,PwXc(:,4));
    OgiveXw=ogive(FwXw,PwXw(:,4));

    % temperature

    %[Fwt,Pwt,sigt] = cross_spectra(w(iok),Ts(iok),PFlux.sf,Nwt);
    %[Fwtfco2,Pwtfco2,sigtfco2] = cross_spectra(w(iok),Tsco2f(iok),PFlux.sf,NwTsco2);
    %[Fwtfh2o,Pwtfh2o,sigtfh2o] = cross_spectra(w(iok),Tsh2of(iok),PFlux.sf,NwTsh2o);
    [Fwt,Pwt,Wonky] = cross_spectra(w(iok),Ts(iok),PFlux.sf,Nwt);
    [Fwtfco2,Pwtfco2,Wonky] = cross_spectra(w(iok),Tsco2f(iok),PFlux.sf,NwTsco2);
    [Fwtfh2o,Pwtfh2o,Wonky] = cross_spectra(w(iok),Tsh2of(iok),PFlux.sf,NwTsh2o);

    covwc=cov(w(iok),Xc(iok));covwc=covwc(1,2);
    covww=cov(w(iok),Xw(iok));covww=covww(1,2);
    covwt=cov(w(iok),Ts(iok));covwt=covwt(1,2);


    subplot(223);
    semilogx(Fwt,Fwt.*Pwt(:,4)/covwt,Fwtfco2,Fwtfco2.*Pwtfco2(:,4)/covwt,FwXc,FwXc.*PwXc(:,4)/covwc);
    ylabel('f*S(f)');
    set(gca,'xlim',[1/1800 1]);
    title('wt (blue),  wt (co2 filter, green),  wco2 (red)');

    subplot(224);
    semilogx(Fwt,Fwt.*Pwt(:,4)/covwt,Fwtfh2o,Fwtfh2o.*Pwtfh2o(:,4)/covwt,FwXw,FwXw.*PwXw(:,4)/covww);
    ylabel('f*S(f)');
    set(gca,'xlim',[1/1800 1]);
    title('wt (blue),  wt (h2o filter, green),  wh2o (red)');

    figure(3);

    subplot(5,4,9);
    hold on;
    loglog(FwXc,FwXc.*PwXc(:,2),'r--');
    ylabel('f*S(f)');
    set(gca,'xlim',[1/1800 1]);
    set(gca,'fontsize',8);
    title('P_{w} - closed calculation dashed','fontsize',8);
    hold off;

    subplot(5,4,10);
    hold on;
    loglog(Pt(:,1),Pt(:,1).*Pt(:,2),'r--');
    ylabel('f*S(f)');
    set(gca,'xlim',[1/1800 1]);
    title('P_{T} ','fontsize',8);
    hold off;
    set(gca,'fontsize',8);

    subplot(5,4,11);
    hold on;
    loglog(FwXc,FwXc.*PwXc(:,3),'r--');
    ylabel('f*S(f)');
    set(gca,'xlim',[1/1800 1]);
    title('P_{CO2}','fontsize',8);
    hold off;
    set(gca,'fontsize',8);

    subplot(5,4,12);
    hold on;
    loglog(FwXw,FwXw.*PwXw(:,3),'r--');
    ylabel('f*S(f)');
    set(gca,'xlim',[1/1800 1]);
    title('P_{H2O}','fontsize',8);
    hold off;
    set(gca,'fontsize',8);

    subplot(527);
    hold on;
    semilogx(FwXc,FwXc.*PwXc(:,4),'r--',[FwXc(1) FwXc(end)],[0 0],'k');
    ylabel('f*S(f)');
    set(gca,'xlim',[1/1800 1]);
    title('CO_{WCO2} - closed dashed','fontsize',8);
    hold off;
    set(gca,'fontsize',8);

    subplot(528);
    hold on
    semilogx(FwXc,OgiveXc,'r--');
    set(gca,'xlim',[1/1800 1]);
    ylabel('Ogive');
    title('OGIVE  WCO2','fontsize',8);
    hold off;
    set(gca,'fontsize',8);

    subplot(529);
    hold on;
    semilogx(FwXw,FwXw.*PwXw(:,4),'r--',[FwXw(1) FwXw(end)],[0 0],'k');
    set(gca,'xlim',[1/1800 1]);
    ylabel('f*S(f)');
    title('Co_{WH2O}','fontsize',8);
    hold off;
    set(gca,'fontsize',8);

    subplot(5,2,10);
    hold on;
    semilogx(FwXw,OgiveXw,'r--');
    set(gca,'xlim',[1/1800 1]);
    ylabel('Ogive');
    title('OGIVE H2O','fontsize',8);
    hold off;
    set(gca,'fontsize',8);

end
%pause
return




function [F,P,SIGOUT] = cross_spectra(a,b,sf,N)

if N>0

    a=a(1:end-N);size(a);
    b=b(N+1:end);size(b);

elseif N<0

    N=abs(N);

    a=a(N+1:end);size(a);
    b=b(1:end-N);size(b) ;

end

% normalize if necessary

if nargin > 5

    if FLAG==1

        P=hspec_ss(a,b,sf);

        covar=cov(a,b);

        P(:,2)=P(:,2)/std(a);
        P(:,3)=P(:,3)/std(b);
        P(:,4)=P(:,4)/covar(1,2);

        disp('Spectra/cospectrum normalized');

    end

else

    P=hspec_ss(a,b,sf);

end

F=P(:,1);
SIGOUT=[a;b];

return





% $$$     figure(4);
% $$$
% $$$     yunit=1/62;
% $$$     yheight=14*yunit;
% $$$
% $$$     xleft=.25;
% $$$     xlength=.5;
% $$$
% $$$     axes('position',[xleft yunit*40 xlength yheight]);
% $$$     set(gca,'nextplot','add')
% $$$
% $$$     loglog(FwXc,FwXc.*PwXc(:,3),'k--',[1/2.9 1/2.9],[1e-2 5],'k','linewidth',1.5);
% $$$     ylabel('fS(f)','fontsize',14);
% $$$     set(gca,'xlim',[1/1800 1],'ylim',[1e-2 5]);
% $$$     set(gca,'fontsize',12);
% $$$     set(gca,'xticklabel',[]);
% $$$     set(gca,'linewidth',2)
% $$$     text(1e-3,3,'a)','fontsize',14);
% $$$
% $$$     axes('position',[xleft yunit*25 xlength yheight]);
% $$$     set(gca,'nextplot','add')
% $$$
% $$$     semilogx(FwXc,FwXc.*PwXc(:,4),'k--',[FwXc(1) FwXc(end)],[0 0],'k',[1/2.9 1/2.9],[-.05 .2],'k','linewidth',1.5);
% $$$     ylabel('fCo_{wc}(f)','fontsize',14);
% $$$
% $$$
% $$$     set(gca,'xlim',[1/1800 1],'ylim',[-.05 .2]);
% $$$     set(gca,'fontsize',12);
% $$$     set(gca,'xticklabel',[]);
% $$$     set(gca,'linewidth',2)
% $$$     text(1e-3,.15,'b)','fontsize',14);
% $$$
% $$$     axes('position',[xleft yunit*10 xlength yheight]);
% $$$     set(gca,'nextplot','add')
% $$$
% $$$     semilogx(FwXc,OgiveXc,'k--',[1/2.9 1/2.9],[0 .4],'k','linewidth',1.5);
% $$$     set(gca,'xlim',[1/1800 1],'ylim',[0 .4]);
% $$$     ylabel('Ogive','fontsize',14);
% $$$     set(gca,'fontsize',12);
% $$$     set(gca,'linewidth',2)
% $$$     xlabel('Freq. (Hz)');
% $$$     text(1e-3,.35,'c)','fontsize',14);
