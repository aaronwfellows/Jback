function [UVWROT,SONDIAG,UVWTMEAN,THETA,UVWTVAR,COVUVWT,UVWTSKEW,UVWTKURT,USTAR,HBUOYANT,TRANSPORT]=fluxcsat3(uvwt,diagnostic,dataset,rotation,params)

% aug 13, 2002 - modify so that the despike routine executes
% after the  sonic diagnostic routine. before, if the sonic
% diagnostic was present, the despiking wasnt done, but this
% would overlook NaNs for example, and its probably better to do
% despiking anyway

% aug 5, 2002 - modify to accept the diagnostic data from the
% Canada data.  The program returns 4 diagnostic bits in
% locations 0,1,2,3, instead of what is output by the CSAT3,
% which has the diagnostic bits between 12, 13, 14, 15

if nargin < 4
    rotation='double';
elseif isempty(rotation)
    rotation='double';
    params=[];
end

if nargin < 3
    dataset = 'brazil';
end


%
% 10-15-2001 - modified to send SONDIAG to coordrot. That routine
% was set up for it already, but I wasn't using it. I'm not sure
% how that effected the earlier fluxes.
%
%     disp('9/23/2001: adding approximate buoyant heat flux to output');
%
% processes the measured sonic outputs from the campbell CSAT 3
%
% INPUTS:
%
% uvwt - NX4 array
%    ROW 1: measured sonic u component
%    ROW 2: measured sonic v componen
%    ROW 3: measured sonic w component
%    ROW 4: measured sonic t component
%
% OUTPUTS:
%
% UVWROT - NX3 array - rotated sonic components such that mean(v) = mean(w) = 0
%    ROW 1: sonic component rotated into the mean wind direction
%    ROW 2: sonic cross-wind component
%    ROW 3: sonic w component
%
% SONDIAG - NX1 -diagnostic vector for sonic for each sample contains a 1 if the
%    measurement was good or a 0 if there was a spike
%
% UVWTMEAN - 4X1 - mean values for (despiked) sonic measurements in measured (not rotated) coordinates
%    ROW 1: mean measured u component
%    ROW 2: mean measured v componentUVWTMEAN_
%    ROW 3: mean measured w component
%    ROW 4: mean measured sonic temperature
%
% THETA: - 1X1 - meteorological mean wind angle - it is the compass angle in degrees that
%        the wind is blowing FROM (0 = North, 90 = east, etc)
%
% UVWTVAR - 4X1 -  variances of ROTATED wind components and the sonic temperature
%    ROW 1: along-wind velocity variance
%    ROW 2: cross-wind velocity variance
%    ROW 3: vertical-wind velocity variance
%    ROW 4: sonic temperature variance
%
% COVUVWT - 6X1 - covariances of ROTATED wind components and the sonic temperature
%    ROW 1: uw co-variance
%    ROW 2: vw co-variance
%    ROW 3: uv co-variance
%    ROW 4: ut co-variance
%    ROW 5: vt co-variance
%    ROW 6: wt co-variance
%
% USTAR - NX1 friction velocity (m/s)

if ~isnan(median(uvwt(4,:))) && ~isempty(median(uvwt(4,:)))
    if median(uvwt(4,:)) > 100
        uvwt(4,:) = uvwt(4,:)-273.15;
    end
end

[iu] = despike(uvwt(1,:),8,-20,20,'U');  uvwt(1,find(~iu)) = NaN*ones(size(find(~iu))); %#ok<*FNDSB>
[iv] = despike(uvwt(2,:),8,-20,20,'V');  uvwt(2,find(~iv)) = NaN*ones(size(find(~iv)));
[iw] = despike(uvwt(3,:),8,-20,20,'W');  uvwt(3,find(~iw)) = NaN*ones(size(find(~iw)));
[it] = despike(uvwt(4,:),8,-50,50,'T');  uvwt(4,find(~it)) = NaN*ones(size(find(~it)));

SONDESPIKE = ones(1,length(iu));

SONDESPIKE( find( sum([iu;iv;iw;it]) < 4 ) ) = zeros(1, length(find( sum([iu;iv;iw;it]) < 4 )) );

if nargin > 1
    SONDIAG = SONDESPIKE & diagword(floor(abs(diagnostic)),dataset);
end

%  Process the data only if there are enough good records
%disp(['Good Points:  ' int2str(length(find(SONDIAG))) ',    Bad Points:  ' int2str(length(find(~SONDIAG)))]);

if length(find(SONDIAG))>300

    UVWTMEAN =  [ mean(uvwt(1,find(SONDIAG))); mean(uvwt(2,find(SONDIAG))); mean(uvwt(3,find(SONDIAG))); mean(uvwt(4,find(SONDIAG))) ];

    THETA = (-atan2(UVWTMEAN(2),UVWTMEAN(1)) + pi/2) *180/pi;
    if THETA < 0; THETA = THETA+360; end

    % calculate statistics that require all channels are despiked
    % rotate the measured wind vector.

    if findstr(rotation,'double')

        UVWROT  = coordrot(uvwt(1:3,:),0,SONDIAG);

        covs    = cov( [ UVWROT(:,find(SONDIAG) ); uvwt(4,find(SONDIAG)) ]' );
        UVWTROT = [ UVWROT(:,find(SONDIAG) ); uvwt(4,find(SONDIAG)) ] ;
        qsqr    =.5*( sum( UVWROT(:,find(SONDIAG) ).^2 ) );

        iok   = find(SONDIAG);
        [hs,Wonky] = covarmax( UVWROT(3,iok) , uvwt(4,iok) , 16 );

    elseif findstr(rotation,'planar')

        % Find Wind Sector: will there be a problem with absolute sonic
        % angle due to orientation?
        % fix problem with mean between 359 and 360

        THETA = (-atan2(UVWTMEAN(2),UVWTMEAN(1)) + pi/2) *180/pi;
        if THETA < 0; THETA = THETA+360; end


        % The row is the mean of the

        params(361,1:3)=params(1,1:3);
        params=interp1([0:360]',params,THETA); %#ok<NBRAK>

        %pause
        % Find correct row


        [UVWROT,gamma]  = coordrot(uvwt(1:3,:),2,SONDIAG,params);

        % rotate back by gamma about the z-axis.  This puts the
        % measurements roughly in the measurement frame, although
        % rotated into the planar fit coordinates.  For UVWROT, the mean
        % v-component is zero, but w-mean is not necessarily zero.  When
        % we rotate around gamma, v-mean will be the north to south wind
        % component - it should be similar to the measured v-bar...
        %
        % trans.m does the coordinate transformation, and the FLAG=1
        % means that the vector is rotated back to the measurement frame

        uvwrot = trans( mean(UVWROT(:,find(SONDIAG))')',[0 0 gamma]',1); %#ok<UDIM>

        UVWTMEAN =  [ UVWTMEAN [ uvwrot; UVWTMEAN(4) ]];

        %pause

        UVWTROT =   [ UVWROT(:,find(SONDIAG) ); uvwt(4,find(SONDIAG)) ]'  ;
        covs = cov( UVWTROT );

        qsqr    =  .5*sum( UVWROT(:,find(SONDIAG) ).^2 );

        iok   = find(SONDIAG);
        [hs,Wonky] = covarmax( UVWROT(3,iok), uvwt(4,iok), 16 );

    end

    % calculate stresses

    %disp('Sonic signals linear detrended before calculating fluxes');

    % add here to rotate horizontally if an angle is sent
    COVUVWT = [ covs(1,3); covs(2,3); covs(1,2); covs(1,4); covs(2,4); covs(3,4)];

    %            u'w'       v'.w'       u'.v'       u't'      v't'       w't'   
    
    UVWTVAR = diag(covs);

    % calculate skewness

    UVWTSKEW = skewness(UVWTROT')';
    UVWTKURT = kurtosis(UVWTROT')';

    USTAR = sqrt( sqrt(covs(1,3)^2 + covs(2,3)^2) );

    % calculate turbulent transport term

    TRANSPORT = mean(UVWROT(3,find(SONDIAG)).*qsqr);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BUOYANCY  FLUX , approximate (W/m^2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    HBUOYANT =  29/1000*38.6*1004*hs;

    % $$$ figure(3);clf
    % $$$
    % $$$ subplot(521)
    % $$$ plot(UVWROT(:,find(SONDIAG))');
    % $$$ set(gca,'xlim',[0 size(UVWROT,2)]);
    % $$$ title('U (bl), V (gr), W (r)')
    % $$$
    % $$$ subplot(522)
    % $$$ plot(uvwt(4,find(SONDIAG))');
    % $$$ set(gca,'xlim',[0 size(UVWROT,2)]);
    % $$$ title(['Ts']);

    %drawnow

else    % case of not enough good points

    UVWROT    = NaN*ones(3,size(uvwt,2));
    THETA     = NaN;
    UVWTVAR   = NaN*ones(4,1);
    COVUVWT   = NaN*ones(6,1);
    UVWTSKEW  = NaN*ones(4,1);
    UVWTKURT  = NaN*ones(4,1);
    USTAR     = NaN;
    HBUOYANT  = NaN;
    TRANSPORT = NaN;


    if findstr(rotation,'planar')
        UVWTMEAN  = NaN*ones(4,2);
    else
        UVWTMEAN  = NaN*ones(4,1);
    end


end

return


function SONDIAG = diagword(diagnostic,dataset)

% set to absoulte value of diagnostic, because there was a negative
% in one of the diagnostic words, which shouldnt be possible

if strcmp(dataset,'brazil')

    err_amp_low  = bitand(diagnostic,2^12)==2^12;
    err_amp_high = bitand(diagnostic,2^13)==2^13;
    err_no_lock  = bitand(diagnostic,2^14)==2^14;
    err_sos      = bitand(diagnostic,2^15)==2^15;

elseif strcmp(dataset,'canada')

    err_amp_low  = bitand(diagnostic,2^0)==2^0;
    err_amp_high = bitand(diagnostic,2^1)==2^1;
    err_no_lock  = bitand(diagnostic,2^2)==2^2;
    err_sos      = bitand(diagnostic,2^3)==2^3;

end

%figure(2);
%subplot(411);
%plot([err_amp_low'])
%title('Amp low');
%subplot(412);
%plot([err_amp_high'])
%title('Amp high');
%subplot(413);
%plot([err_no_lock'])
%title('No lock');
%subplot(414);
%plot([err_sos'])
%title('speed of sound');


SONDIAG = ~(err_amp_low | err_amp_high | err_no_lock | err_sos);


return




