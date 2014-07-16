function [TAU,DELAY,RMIN] = taucalc(sig1,sig2,sf,P,method,PlotInfo)

var_defs();
% PlotInfo - Structure
% PlotInfo.figure  = fig #
% PlotInfo.sub1    = string to locate plot1
% PlotInfo.sub4    = string to locate plot2
% PlotInfo.sub3    = string to locate plot3
% PlotInfo.sub2    = string to locate plot4
% PlotInfo.title   = title
% PlotInfo.legend   = title

% INPUT parameters:

% P.tau0        =  2.0;
% P.taumin      =  1;
% P.taumax      =  5.0;
% P.delay0      = 10.0;
% P.delaymin    =  8.0;
% P.delaymax    = 15.0;
% P.FreqLoTau      =  1/300;
% P.FreqHiTau      =  .6;


plot_spec=0;

% 6/22/2001 - adding second estimate fo tau by degrading sig1 using
%             a first order filter, and comparing power spectra to
%             get a best fit

%  Method 1: Ad Hoc
%
%  A search for a first order filter time constant which is applied to
%  sig2 and gives the closest power spectrum as sig1. One concern here
%  is that the phase spectrum is not explicitly treated, only the
%  amplitude spectrum is.
%
%  There are several parameters need to allow enough flexibility for
%  the script to work% over a broader range of conditions.
%
%  Taus -   time constant parameters - for a first order filter, it
%           represents 1/freq where the signal has decayed to a certain
%           percentage of the input power - (.5 or something?)
%
%  taumin - the minimum realistic time constant.  For example, tau
%           should be greater than zero.
%  taumax - the maximum realistic time constant
%
%
%  Frequencies - the method for finding the 'best fit' time constant is
%              to calculate a normalized measure R of the error in the fit for
%              each value of tau, then to choose R tau that minimizes the
%              error.  The range of frequencies over which to calculate the
%              error in the fit is specified by
%
%         FreqLo - the lowest frequency for which the difference between
%                the measured spectrum for sig2 and the filtered spectrum
%                for sig1 should be computed.
%
%         FreqHi - highest frequency to compute the fit
%

% detrend

sig1 = detrend(sig1,1)/std(sig1);
sig2 = detrend(sig2,1)/std(sig2);

% set the length of sig1, sig2 to the nearest power of 2, or at least to an even number for now
disp('taucalc.m - consider setting signals going into hspec to an even multiple of 2 samples...');

if rem(length(sig1),2)==1
    sig1=sig1(1:length(sig1)-1);
    sig2=sig2(1:length(sig2)-1);
end

%sizesig1=size(sig1)
%sizesig2=size(sig2)

if strcmp(method,'adhoc')

    DELAY=[];

    taus = P.AdHoctaumin:(P.AdHoctaumax-P.AdHoctaumin)/20:P.AdHoctaumax;

    for itau = 1:length(taus)

        wn=1/taus(itau)/(sf/2);
        [b,a]=butter(1,wn,'low');

        % degrade temperature

        sig1f = filter(b,a,sig1);

        P2=hspec_ss(sig1f,sig2,sf);
        P1=hspec_ss(sig1,0,sf);
        f1=P2(:,1);

        ihi =  find(abs(f1-P.AdHocFreqHiTau)==min(abs(f1-P.AdHocFreqHiTau)));
        ilo =  find(abs(f1-P.AdHocFreqLoTau)==min(abs(f1-P.AdHocFreqLoTau)));

        % Get the fit
        R(itau)= sum( (log(P2(ilo:ihi,2))-log(P2(ilo:ihi,3))).^2);
        disp('itau tau R(itau)');
        disp(num2str([itau taus(itau) R(itau) ]));

        % ESTIMATE TIME DELAY FROM THE PHASE SPECTRUM
        %
        % The delay is approximately given by
        %
        % phi(f) = 2*pi*f*DT
        %
        %        1  d phi(f)
        % DT = ---- -------
        %      2*pi   df

        ihi =  find(abs(f1-P.AdHocFreqHiDelay)==min(abs(f1-P.AdHocFreqHiDelay)));
        ilo =  find(abs(f1-P.AdHocFreqLoDelay)==min(abs(f1-P.AdHocFreqLoDelay)));

        phi   = unwrap(atan2(P2(:,5),P2(:,4)));
        DELAY = (1/(2*pi)) * (phi(ihi)-phi(ilo) ) / (f1(ihi)-f1(ilo))  ;


        figure(1); clf
        subplot(211);
        loglog(f1,f1.*P2(:,2),f1,f1.*P2(:,3));
        title('sig1 blue, sig2 green');
        subplot(212);
        semilogx(f1,phi,[f1(ihi) f1(ilo)],[phi(ihi) phi(ilo)],'ro');
        title('phase angle');

        %disp('pause')
        %pause

    end


    TAU = taus(R==min(R));
    %%% TAU = taus(find(R==min(R)));
    RMIN = min(R);


    if nargin>5

        wn=1/TAU/(sf/2);
        [b,a]=butter(1,wn,'low');

        % degrade temperature
        sig1f = filter(b,a,sig1);
        P2=hspec_ss(sig1f,sig2,sf);

        %%FIGURE subplot
        eval([ 'figure(' int2str(PlotInfo.figure) ');']);
        eval([ 'subplot' PlotInfo.sub1 ';']);

        plot([sig2' sig1' sig1f'])
        set(gca,'xlim',[0 length(sig1)])
        set(gca,'fontsize',8);
        title( [ PlotInfo.title1 '  ' ]);

        %%FIGURE subplot
        eval([ 'subplot' PlotInfo.sub2 ';']);
        loglog(f1,f1.*P2(:,3),f1,f1.*P1(:,2),f1,f1.*P2(:,2),'linewidth',1.75);

        set(gca,'xlim',[10^(-3.5) 4]);
        set(gca,'ylim',[10^(-5)   1]);

        ylims=get(gca,'ylim');
        hold on;
        loglog([P.AdHocFreqLoTau P.AdHocFreqLoTau],ylims,'r',[P.AdHocFreqHiTau P.AdHocFreqHiTau],ylims,'r');
        hold off

        grid on

        title( [ PlotInfo.title2 '  '  num2str(TAU,4)  ]);

        if ~isempty(char(PlotInfo.legend))
            b=legend(char(PlotInfo.legend),3);
            set(b,'fontsize',10);
        end

        xlabel('frequency (Hz)')
        set(gca,'fontsize',10);
        hold off

        drawnow
        %pause
    end


elseif strcmp(method,'shaw')

    RMIN=[];

    disp('using Shaw model');

    %  Method 2:  Apply a model to the measured phase spectrum between the
    %         reference sig1 and degraded sig2 using a non-linear least squares
    %         regression.  The model is in Shaw et al and models the phase
    %         spectrum as the sum of a time shift and a smearing or first
    %         order filter.  A constant time delay causes a phase that is
    %         linear in frequency.  The first order filter has an arctangent
    %         of frequency.  Practically, for our long tubes, the time delay
    %         seems to dominate the phase spectrum so that we don't get a
    %         good estimate of the time constant.  Seems to work better for
    %         shorter tubing.
    %
    %  Parameters:

    % calculate from slope of phase angle plot (note use ORIGINAL timeseries, not shifted)

    % 2/4/99 - use matlab csd to calculate the phase angle, not hspec...

    if plot_spec==1
        figure(105);
        clf;
        subplot(211);
        plot([sig1;sig2]');
    end
    window   = hanning(256*4);
    noverlap = 0;

    % power spectra - only for comparison
    %%% [P12,f] = spectrum.psd(sig1,256*4,sf,window,noverlap,'none');
    %[P12,f] = psd(sig1,256*4,sf,window,noverlap,'none');  %#ok<FDEPR>
    [C12,f] = csd(sig1,sig2,256*4,sf,window,noverlap,'none'); 
    %%% [C12,f] = csd(sig1,sig2,256*4,sf,window,noverlap,'none');
    Co12    = real(C12);
    Q12     = imag(C12);
    phi     = unwrap(atan2(Q12,Co12));
    if plot_spec==1
        subplot(223);semilogx(f,f.*Co12,f,f.*Q12);
        title('co-spectrum (blue), quad-spectrum (green)');
        set(gca,'xlim',[.004 2]);
        subplot(224);semilogx(f,phi);
        title('phase-spectrum');
        set(gca,'xlim',[.004 2]);
    end
    phioff =  median(phi(4:8));

    if phioff > 3*pi/2
        phioff = 2*pi;
    elseif phioff > pi/2 && phioff <= 3*pi/2
        phioff = pi;
    elseif phioff < -pi/2 && phioff >= -3*pi/2
        phioff = -pi;
    elseif phioff < -3*pi/2
        phioff = -2*pi;
    else
        phioff =0;
    end

    % example function from Shaw paper (using our delay of 1.5 sec, and
    % shaw's acetone time constant of .1 s;
    % restrict the frequency range over which to fit the model between
    % 2 and ~80 seconds. (high frequency 'noise' causes significant
    % deviations of data from model when a fit is forced - this frequency
    % range is between index

    % find frequency limits closest

    iHighBeg =  find(abs(f-P.ShawFreqHiTau)==min(abs(f-P.ShawFreqHiTau)));
    iLow     =  find(abs(f-P.ShawFreqLoTau)==min(abs(f-P.ShawFreqLoTau)));

    options=optimset('display','off');

    RNbest    = [];
    taubest   = [];
    delaybest = [];
    ibest=[];

    icount=0;
    ihi=iHighBeg-1;

    tautmp=[];
    ftmp=[];
    deltmp=[];
    rtmp=[];

    % while icount<20

    itnum = 0;
    for ijk=iHighBeg:3:length(f)

        itnum=itnum+1;

        %ihi is the index of the highest freq in the range being tested
        ihi = ihi+1;

        freq_stt = f(iLow);
        freq_end = f(ihi);

        Initial_Guesses = [P.Shawtau0 P.Shawdelay0];
        Min_values = [P.Shawtaumin P.Shawdelaymin];
        Max_values = [P.Shawtaumax P.Shawdelaymax];
        
        %Initial_Guesses = [.3 2];
        %Min_values = [0 0];
        %Max_values = [10 10];
        
        %tries to fit model on increasing ranges of frequency
%         [x2,resnorm] =
%         lsqcurvefit('phasemodel',Initial_Guesses,f(iLow:ihi),phi(iLow:ihi)-phioff,[P.Shawtaumin P.Shawdelaymin],[P.Shawtaumax P.Shawdelaymax],options);
        [x2,resnorm] = lsqcurvefit('phasemodel',Initial_Guesses,f,phi-phioff,Min_values,Max_values,options);



        %RNtest is the diagnostic to see whether the fit to the model is
        %better in this freqency range
        RNtest = sqrt(resnorm)/(ihi-iHighBeg+1);


        if ihi == iHighBeg   % first time through (initial freq range
            RNbest = sqrt(resnorm) ;
            ibest=1;
            icount=0;
            taubest   = x2(1);
            delaybest = x2(2);
            fsttbest = f(ihi);
        elseif RNtest < RNbest
            RNbest = RNtest;
            ibest = ihi-iHighBeg+1;
            fsttbest = f(ihi);
            icount=0;    % stop counting if improvement in fit occurs
            taubest   = x2(1);
            delaybest = x2(2);
        else    % no improvement in fit
            icount=icount+1;
        end

%%%         if ibest==39
%%%             disp('p');
%%%         end
        
        
        %         ftmp is the starting frequency in the range being tested
        ftmp=[ftmp f(ihi)]; %#ok<*AGROW>
        tautmp=[tautmp x2(1)];
        deltmp=[deltmp x2(2)];
        rtmp=[rtmp sqrt(resnorm)/(ihi-iHighBeg+1)];

        %disp(num2str([icount ihi sqrt(resnorm)/(ihi-iHighBeg+1) ibest RNbest taubest delaybest ]));
%         disp(['it #' num2str(itnum) ...
%             ' # ' num2str(icount) ...
%             ', ihi= ' num2str(ihi) ...
%             ', freqs: ' num2str(freq_stt) '-' num2str(freq_end) ...
%             ', fit curr= ' num2str(RNtest) ...
%             ', fit best= ' num2str(RNbest) ...
%             ', ibest= ' num2str(ibest) ...
%             ', fsttbest= ' num2str(fsttbest) ...
%             ', taubest=  ' num2str(taubest) ...
%             ', delaybest= ' num2str(delaybest) ]);

        if ihi==length(f)
            icount=20;
        end

        %new phi fit array
         
        Phi_Fit = [itnum icount ihi freq_stt freq_end RNtest RNbest ibest fsttbest taubest delaybest];
        
        PHI_FIT(ijk,1:length(Phi_Fit)) =  Phi_Fit;
        
    end   % of trying to fit across increasingly higher frequency ranges


    % ok new change here:
    
    %     take the average of the parameters for the 3 best fits 
    % get the rows when an improvement in the fit occurred
    %[TMP_FITS,sortix] = sort(PHI_FIT(:,6));
    [Wonky,sortix] = sort(PHI_FIT(:,6));

    %get three top fits
    Three_Top_Fits_ix = sortix(end-2:end);

    PHI_FITS_TOP3 = PHI_FIT(Three_Top_Fits_ix,:);


    HIFREQ_BEST = median(PHI_FITS_TOP3(:,5));
    TAUBEST = median(PHI_FITS_TOP3(:,10));
    DELAYBEST = median(PHI_FITS_TOP3(:,11));
   
    %get stats on tau and delay
%     tautmp_avg = nanmean(tautmp);
%     tautmp_std = nanstd(tautmp);
%     tautmp_len = lnn(tautmp);
%     tautmp_med = nanmedian(tautmp);
%     tautmp_min = nanmin(tautmp);
%     tautmp_max = nanmax(tautmp);
% 
%     deltmp_avg = nanmean(deltmp);
%     deltmp_std = nanstd(deltmp);
%     deltmp_len = lnn(deltmp);
%     deltmp_med = nanmedian(deltmp);
%     deltmp_min = nanmin(deltmp);
%     deltmp_max = nanmax(deltmp);

    %get medians of good numbers
    %     iok=find(tautmp>.2 & tautmp<2 & deltmp>1.5 & deltmp<5);
    tautmp_iok_bix= tautmp>0 & tautmp<3 ;  %tautmp_iok_ix=find(tautmp_iok_bix);
    deltmp_iok_bix= deltmp>0.25 & deltmp<5; %deltmp_iok_ix=find(deltmp_iok_bix);

    iok = find(tautmp_iok_bix & deltmp_iok_bix);
    % iok = 1:length(tautmp);

    if isempty(iok)
        medtau = 0;
        meddelay = 0;
    else
        medtau = median(tautmp(iok));
        meddelay = median(deltmp(iok));
    end

    if plot_spec == 1
        figure(102);clf
        subplot(411);plot(ftmp);title('ftmp')
        subplot(423);plot(tautmp);
        subplot(425);plot(deltmp);

        subplot(424);plot(ftmp,tautmp,ftmp([1 end]),[ medtau   medtau  ],'r');title('tautmp')
        subplot(426);plot(ftmp,deltmp,ftmp([1 end]),[ meddelay meddelay],'r');title('deltmp')
        subplot(414);plot(rtmp);title('rtmp')
    end
    
    % commented out this original code since I want to see whether the tau
    % and delay with the best fit (smallest RN value) are better to use 
% % % %     TAU   = medtau;   % taubest;
% % % %     DELAY = meddelay; % delaybest;
% % % %     
    
    DELAY = DELAYBEST; 
    TAU = TAUBEST;
    
    
    %F     = f(ibest);
    %RN    = RNbest;

    if nargin > 5

        figure(106);
        clf

        Phasemod1 = phioff + atan(-2*pi*f*TAU)-2*pi*f*DELAY;

        semilogx(f(iLow:ibest),Phasemod1(iLow:ibest),'c','linewidth',6)
        hold on;
        semilogx(f,phi,'k.',f,Phasemod1,'r')
        set(gca,'xlim',[.0001 10]);
        set(gca,'ylim',[-60 20]);
        grid on

        if nargin>5
            title(['Shaw (blue): ',num2str(TAU,4)],'fontsize',14)
        else
            title(['Shaw (blue): ',num2str(TAU,4)],'fontsize',14)
        end

        legend('Modelled Values Best','Phi values from CoSpec','All modelled Values','location','South');

        xlabel('frequency (Hz)')
        ylabel('Phase Angle (rad)')
        hold off

        drawnow
        
        %pause
    end   % of whether to plot or not

end   % of method choice


return

