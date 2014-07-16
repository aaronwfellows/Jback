function [N_obs, headerN_obs, D_obs_out] = combinedinfo(siteName, D, D_obs_out, N_obs, a)
%checks the merge and outputs info on what you removed
figout ='..\combined\CombinedInfo\';
figout = [figout siteName '\'];

if a ==1
    disp('calculating how many obs you have')
    n_row = size(D,1);
    n_col = size(D,2);
    N_obs =ones(n_row,2)*NaN;
    headerN_obs = {'N_obs'};
    
    for row = 1:n_row
        good = isfinite(D(row,:)) == 1 & isnan(D(row,:)) == 0;
        D_obs_out(row,:)=good(1,:); 
        n_measure = sum(good);
        N_obs(row,1) = n_col;
        N_obs(row,2) = n_measure;
    end
    
elseif a ==2
    
    %D grew because we added new rows - calculations
    addNaN=ones(34,2)*NaN;
    N_obs=[N_obs; addNaN];
    
    disp('calculating how many obs you removed')
    n_row = size(D,1);
    N_obs2 =ones(n_row,1)*NaN;
    
    for row = 1:n_row
        good = isfinite(D(row,:)) == 1 & isnan(D(row,:)) == 0;
        prefilter_good=D_obs_out(row,:);
        D_obs_out(row,:) = prefilter_good(1,:) - good(1,:);
        n_measure = sum(good);
        N_obs2(row,1) = n_measure;
    end
       
    removed_obs = N_obs(:,2) - N_obs2(:,1);
    f_obs = removed_obs./N_obs(:,2);
    
    N_obs=[N_obs(:,1), N_obs(:,2), N_obs2, removed_obs, f_obs];
    headerN_obs = {'Length_D', 'N_obs', 'N_obs after filter', 'n removed obs', 'fraction obs removed'};
    
    %--------------------------------------------------------------%
    %combined ust
    Ust(1,:)=D(21,:);
    bad=isnan(Ust(1,:));
    Ust(1,bad)=D(322,bad);
    bad=isnan(Ust(1,:));
    Ust(1,bad)=D(323,bad);
    %we will use just ust filtered observations 
    turb = Ust(1,:) >0.3;
    
    %Gain and delay settings
    co2_delay_s = nanmean(D(55,:))
    h2o_delay_s = nanmean(D(57,:))
    co2_gain = nanmean(D(58,turb))
    h2o_gain = nanmean(D(59,turb))
    
    Gain_Delay = [co2_delay_s; h2o_delay_s; co2_gain; h2o_gain];
    headerGain_Delay ={'co2_delay_s'; 'h2o_delay_s'; 'co2_gain (ust>0.3)'; 'h2o_gain (ust>0.3)'};
    %date;
    Date_data_compiled_on = date;
    %--------------------------------------------------------------%
    %save this information
    save([figout 'CombinedInfo_' siteName '.mat'], 'Date_data_compiled_on', 'headerGain_Delay', 'Gain_Delay', 'N_obs', 'headerN_obs', 'D_obs_out');
    %--------------------------------------------------------------%
    %--------------------------------------------------------------%
    disp('you entered the combinedinfo subroutine')
    disp('- This will check what you are going to combined and output info.')
    disp('There is a pause between plots.  Hit enter to move to the next variable.')
    %--------------------------------------------------------------%
    %--------------------------------------------------------------%
    f='fraction obs removed';
    figure(1)
    bar(N_obs(:,5), 'k')
    xlabel('row pointing to DMERGE variable (see HMERGE)')
    ylabel(f)
    
    saveas(figure(1), [figout f]);
    
    pause
    close(figure(1))
    %--------------------------------------------------------------%
    %plot Tau and plt delay
    f='co2 delay (s)';
    figure(1)
    plot(D(1,turb) ,D(55,turb), 'k.')
    xlabel('tower time (day)')
    ylabel(f)
    
    saveas(figure(1), [figout f]);
    
    pause
    close(figure(1)) 
    %--------------------------------------------------------------%
    f='h2o delay (s)';
    figure(1)
    plot(D(1,turb) ,D(57,turb), 'k.')
    xlabel('tower time (day)')
    ylabel(f)
    
    saveas(figure(1), [figout f]);
    
    pause
    close(figure(1)) 
    %--------------------------------------------------------------%
    f='co2 gain';
    figure(1)
    plot(D(1,turb) ,D(58,turb), 'k.')
    xlabel('tower time (day)')
    ylabel(f)
    
    co2_gain = nanmean(D(58,turb))
    
    saveas(figure(1), [figout f]);
    
    pause
    close(figure(1)) 
    %--------------------------------------------------------------%
    f='h2o gain';
    figure(1)
    plot(D(1,:) ,D(59,:), 'k.')
    xlabel('tower time (day)')
    ylabel(f)
    
    h2o_gain = nanmean(D(59,turb))
    
    saveas(figure(1), [figout f]);
    
    pause
    close(figure(1)) 
    % ------------------------------------------------------------------%
    %sonic T
    f='sonic T(^oC)';
    figure(1)
    %fast
    plot(D(1,:), D(8,:), 'ko')
    hold on

    %dl T
    plot(D(1,:), D(84,:), 'bx')
    hold on

    %goes T
    plot(D(1,:), D(210,:), 'r.')
    legend('k = mat', 'b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel(f)

    f2='sonic T(^oC) xy';
    figure(2)
    plot(D(8,:), D(84,:), 'bo')
    hold on
    plot(D(8,:), D(210,:), 'r.')
    legend('b=mat v. dl', 'r=mat v. goes')
    xlabel(' T (^oC)')
    ylabel(' T (^oC)')

    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2)) 
    
    % ------------------------------------------------------------------%
    %T actual
    f= 'T actual (^oC)';
    figure(1)
    %fast
    plot(D(1,:), D(330,:), 'ko')
    hold on

    %dl T
    plot(D(1,:), D(331,:), 'bx')
    hold on

    %goes T
    plot(D(1,:), D(332,:), 'r.')
    legend('k = mat', 'b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel(f)

    f2 = 'T actual (^oC) xy';
    figure(2)
    plot(D(330,:), D(331,:), 'bo')
    hold on
    plot(D(330,:), D(332,:), 'r.')
    legend('b=mat v. dl', 'r=mat v. goes')
    xlabel(' T (^oK)')
    ylabel(' T (^oK)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2)) 
    
    % ------------------------------------------------------------------%
    %Ustar
    disp('P1 is the slope and intercept of mat and dl')
    disp('P2 is the slope and intercept of mat and goes')
    
    isn=isnan(D(21,:)) ==1 | isnan(D(322,:)) == 1;
    P1=polyfit(D(21,~isn),D(322,~isn),1)
    B1 = robustfit(D(21,~isn),D(322,~isn))
    
    isn=isnan(D(21,:)) ==1 | isnan(D(323,:)) == 1;
    P2=polyfit(D(21,~isn),D(323,~isn),1)
    B2 = robustfit(D(21,~isn),D(323,~isn))
    
    f='Friction velocity (m_s-1)';
    figure(1)
    %fast
    plot(D(1,:), D(21,:), 'ko')
    hold on

    %dl T
    plot(D(1,:), D(322,:), 'bx')
    hold on

    %goes T
    plot(D(1,:), D(323,:), 'r.')
    legend('k = mat', 'b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel(f)

    f2='Friction velocity (m_s-1) mat_dl';
    figure(2)
    plot(D(21,:), D(322,:), 'b.')
    hold on
    plot(D(21,:), P1(1,1) * D(21,:) + P1(1,2), 'g-')
    hold on
    plot(D(21,:), B1(2,1) * D(21,:) + B1(1,1), 'k-')
    legend('b=mat v. dl', 'polyfit', 'robustfit')
    xlabel('mat Friction velocity (m/s)')
    ylabel('dl Friction velocity (m/s)')
    
    f3='Friction velocity (m_s-1) mat_goes';
    figure(3)
    plot(D(21,:), D(323,:), 'r.')
    hold on
    plot(D(21,:), P2(1,1) * D(21,:) + P2(1,2), 'g-')
    hold on
    plot(D(21,:), B2(2,1) * D(21,:) + B2(1,1), 'k-')
    legend('r=mat v. goes', 'polyfit', 'robustfit')
    xlabel('mat Friction velocity (m/s)')
    ylabel('goes Friction velocity (m/s)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    saveas(figure(3), [figout f3]);
    close(figure(1),figure(2),figure(3)) 
    % ------------------------------------------------------------------%
    %wind dir
    f='winddir';
    figure(1)
    %fast
    plot(D(1,:), D(324,:), 'ko')
    hold on

    %dl T
    plot(D(1,:), D(325,:), 'bx')
    hold on

    %goes T
    plot(D(1,:), D(326,:), 'r.')
    legend('k = mat', 'b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel(f)

    f2='winddir xy';
    figure(2)
    plot(D(324,:), D(325,:), 'b.')
    hold on
    plot(D(324,:), D(326,:), 'r.')
    legend('b=mat v. dl', 'r=mat v. goes')
    xlabel('wind dir (deg)')
    ylabel('wind dir (deg)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2)) 
    
    % ------------------------------------------------------------------%
    % wind speed
    disp('P1 is the slope and intercept of mat and dl')
    disp('P2 is the slope and intercept of mat and goes')
    
    isn=isnan(D(327,:)) ==1 | isnan(D(328,:)) == 1;
    P1=polyfit(D(327,~isn),D(328,~isn),1)
    
    isn=isnan(D(327,:)) ==1 | isnan(D(329,:)) == 1;
    P2=polyfit(D(327,~isn),D(329,~isn),1)

    f='wind sp(m_s-1)';
    figure(1)
    %fast
    plot(D(1,:), D(327,:), 'ko')
    hold on

    %dl T
    plot(D(1,:), D(328,:), 'bx')
    hold on

    %goes T
    plot(D(1,:), D(329,:), 'r.')
    legend('k = mat', 'b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel(f)

    f2='wind sp(m_s-1) mat_dl';
    figure(2)
    plot(D(327,:), D(328,:), 'b.')
    hold on
    plot(D(327,:), P1(1,1) * D(327,:) + P1(1,2), 'g-')
    legend('b=mat v. dl')
    xlabel('mat wind sp(m/s)')
    ylabel('dl wind sp(m/s)')
    
    f3='wind sp(m_s-1) mat_goes';
    figure(3)
    plot(D(327,:), D(329,:), 'r.')
    hold on
    plot(D(327,:), P2(1,1) * D(327,:) + P2(1,2), 'g-')
    legend('r=mat v. goes')
    xlabel('mat wind sp (m/s)')
    ylabel('goes wind sp (m/s)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    saveas(figure(3), [figout f3]);
    close(figure(1),figure(2),figure(3)) 
    % ------------------------------------------------------------------%
    % dry density
    disp('P1 is the slope and intercept of mat and dl')
    disp('P2 is the slope and intercept of mat and goes')
    
    isn=isnan(D(336,:)) ==1 | isnan(D(337,:)) == 1;
    P1=polyfit(D(336,~isn),D(337,~isn),1)
    
    isn=isnan(D(336,:)) ==1 | isnan(D(338,:)) == 1;
    P2=polyfit(D(336,~isn),D(338,~isn),1)

    f='dry density (mol_m^-3)';
    figure(1)
    %fast
    plot(D(1,:), D(336,:), 'ko')
    hold on

    %dl T
    plot(D(1,:), D(337,:), 'bx')
    hold on

    %goes T
    plot(D(1,:), D(338,:), 'r.')
    legend('k = mat', 'b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel(f)

    f2='dry density (mol_m^-3) mat_dl';
    figure(2)
    plot(D(336,:), D(337,:), 'b.')
    hold on
    plot(D(336,:), P1(1,1) * D(336,:) + P1(1,2), 'g-')
    legend('b=mat v. dl')
    xlabel('mat dry density (mol/m^3)')
    ylabel('dl dry density (mol/m^3)')
    
    f3='dry density (mol_m^-3) mat_goes';
    figure(3)
    plot(D(336,:), D(338,:), 'r.')
    hold on
    plot(D(336,:), P2(1,1) * D(336,:) + P2(1,2), 'g-')
    legend('r=mat v. goes')
    xlabel('mat dry density (mol/m^3)')
    ylabel('goes dry density (mol/m^3)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    saveas(figure(3), [figout f3]);
    close(figure(1),figure(2),figure(3)) 
    % ------------------------------------------------------------------%
    % moist density
    disp('P1 is the slope and intercept of mat and dl')
    disp('P2 is the slope and intercept of mat and goes')
    
    isn=isnan(D(333,:)) ==1 | isnan(D(334,:)) == 1;
    P1=polyfit(D(333,~isn),D(334,~isn),1)
    
    isn=isnan(D(333,:)) ==1 | isnan(D(338,:)) == 1;
    P2=polyfit(D(333,~isn),D(335,~isn),1)

    f='moist density (mol_m^-3)';
    figure(1)
    %fast
    plot(D(1,:), D(333,:), 'ko')
    hold on

    %dl T
    plot(D(1,:), D(334,:), 'bx')
    hold on

    %goes T
    plot(D(1,:), D(335,:), 'r.')
    legend('k = mat', 'b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel(f)

    f2='moist density (mol_m^-3) mat_dl';
    figure(2)
    plot(D(333,:), D(334,:), 'b.')
    hold on
    plot(D(333,:), P1(1,1) * D(333,:) + P1(1,2), 'g-')
    legend('b=mat v. dl')
    xlabel('mat moist density (mol / m^3)')
    ylabel('dl moist density (mol / m^3)')
    
    f3='moist density (mol_m^-3) mat_goes';
    figure(3)
    plot(D(333,:), D(335,:), 'r.')
    hold on
    plot(D(333,:), P2(1,1) * D(333,:) + P2(1,2), 'g-')
    legend('r=mat v. goes')
    xlabel('mat moist density (mol / m^3)')
    ylabel('goes moist density (mol / m^3)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    saveas(figure(3), [figout f3]);
    close(figure(1),figure(2),figure(3)) 
    % ------------------------------------------------------------------%
    %sensible heat
    disp('P1 is the slope and intercept of mat and dl')
    disp('P2 is the slope and intercept of mat and goes')
    
    isn=isnan(D(29,:)) ==1 | isnan(D(339,:)) == 1;
    P1=polyfit(D(29,~isn),D(339,~isn),1)
    
    isn=isnan(D(29,:)) ==1 | isnan(D(340,:)) == 1;
    P2=polyfit(D(29,~isn),D(340,~isn),1)

    f='sens heat (W_m^-2)';
    figure(1)
    %fast
    plot(D(1,:), D(29,:), 'ko')
    hold on

    %dl T
    plot(D(1,:), D(339,:), 'bx')
    hold on

    %goes T
    plot(D(1,:), D(340,:), 'r.')
    legend('k = mat', 'b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel(f)

    f2 = 'sens heat (W_m^-2) mat_dl';
    figure(2)
    plot(D(29,:), D(339,:), 'b.')
    hold on
    plot(D(29,:), P1(1,1) * D(29,:) + P1(1,2), 'g-')
    legend('b=mat v. dl')
    xlabel('mat sens heat (W/m^2)')
    ylabel('dl sens heat (W/m^2)')
    
    f3 = 'sens heat (W_m^-2) mat_goes';
    figure(3)
    plot(D(29,:), D(340,:), 'r.')
    hold on
    plot(D(29,:), P2(1,1) * D(29,:) + P2(1,2), 'g-')
    legend('r=mat v. goes')
    xlabel('mat sens heat (W/m^2)')
    ylabel('goes sens heat (W/m^2)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    saveas(figure(3), [figout f3]);
    close(figure(1),figure(2),figure(3))    
    % ------------------------------------------------------------------%
    %latent heat
    disp('P1 is the slope and intercept of mat and dl')
    disp('P2 is the slope and intercept of mat and goes')
    
    isn=isnan(D(30,:)) ==1 | isnan(D(343,:)) == 1;
    P1=polyfit(D(30,~isn),D(343,~isn),1)
    B1 = robustfit(D(30,~isn),D(344,~isn))
    
    isn=isnan(D(30,:)) ==1 | isnan(D(344,:)) == 1;
    P2=polyfit(D(30,~isn),D(344,~isn),1)
    B2 = robustfit(D(30,~isn),D(344,~isn))

    f= 'LE (W_m^-2)';
    figure(1)
    %fast
    plot(D(1,:), D(30,:), 'ko')
    hold on

    %dl T
    plot(D(1,:), D(343,:), 'bx')
    hold on

    %goes T
    plot(D(1,:), D(344,:), 'r.')
    legend('k = mat', 'b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('LE (W/m^2)')

    f2= 'LE (W_m^-2) mat_dl';
    figure(2)
    plot(D(30,:), D(343,:), 'b.')
    hold on
    plot(D(30,:), P1(1,1) * D(30,:) + P1(1,2), 'g-')
    hold on
    plot(D(30,:), B1(2,1) * D(30,:) + B1(1,1), 'k-')
    legend('b=mat v. dl', 'polyfit', 'robustfit')
    xlabel('mat LE (W/m^2)')
    ylabel('dl LE (W/m^2)')
    
    f3= 'LE (W_m^-2) mat_goes';
    figure(3)
    plot(D(30,:), D(344,:), 'r.')
    hold on
    plot(D(30,:), P2(1,1) * D(30,:) + P2(1,2), 'g-')
    hold on
    plot(D(30,:), B2(2,1) * D(30,:) + B2(1,1), 'k-')
    legend('r=mat v. goes', 'polyfit', 'robustfit')
    xlabel('mat LE (W/m^2)')
    ylabel('goes LE (W/m^2)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    saveas(figure(3), [figout f3]);
    close(figure(1),figure(2),figure(3))
    % ------------------------------------------------------------------%
    %co2 conc
    disp('P1 is the slope and intercept of mat and dl')
    disp('P2 is the slope and intercept of mat and goes')
    
    isn=isnan(D(35,:)) ==1 | isnan(D(85,:)) == 1;
    P1=polyfit(D(35,~isn),D(85,~isn),1)
    
    isn=isnan(D(35,:)) ==1 | isnan(D(211,:)) == 1;
    P2=polyfit(D(35,~isn),D(211,~isn),1)

    f= 'co2 mixing ratio';
    figure(1)
    %fast
    plot(D(1,:), D(35,:), 'ko')
    hold on

    %dl T
    plot(D(1,:), D(85,:), 'bx')
    hold on

    %goes T
    plot(D(1,:), D(211,:), 'r.')
    legend('k = mat', 'b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('co2 mixing ratio (micromol/mol dry air')

    f2= 'co2 mixing ratio mat_dl';
    figure(2)
    plot(D(35,:), D(85,:), 'b.')
    hold on
    plot(D(35,:), P1(1,1) * D(35,:) + P1(1,2), 'g-')
    legend('b=mat v. dl')
    xlabel('mat co2 mixing ratio')
    ylabel('dl co2 mixing ratio')
    
    f3= 'co2 mixing ratio mat_goes';
    figure(3)
    plot(D(35,:), D(211,:), 'r.')
    hold on
    plot(D(35,:), P2(1,1) * D(35,:) + P2(1,2), 'g-')
    legend('r=mat v. goes')
    xlabel('mat co2 mixing ratio')
    ylabel('goes co2 mixing ratio')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    saveas(figure(3), [figout f3]);
    close(figure(1),figure(2),figure(3))
    % ------------------------------------------------------------------%
    %h2o conc
    disp('P1 is the slope and intercept of mat and dl')
    disp('P2 is the slope and intercept of mat and goes')
    
    isn=isnan(D(42,:)) ==1 | isnan(D(86,:)) == 1;
    P1=polyfit(D(42,~isn),D(86,~isn),1)
    
    isn=isnan(D(42,:)) ==1 | isnan(D(212,:)) == 1;
    P2=polyfit(D(42,~isn),D(212,~isn),1)

    f='h2o conc';
    figure(1)
    %fast
    plot(D(1,:), D(42,:), 'ko')
    hold on

    %dl T
    plot(D(1,:), D(86,:), 'bx')
    hold on

    %goes T
    plot(D(1,:), D(212,:), 'r.')
    legend('k = mat', 'b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('h2o conc (mmol/mol dry air)')

    f2='h2o conc mat_dl';
    figure(2)
    plot(D(42,:), D(86,:), 'b.')
    hold on
    plot(D(42,:), P1(1,1) * D(42,:) + P1(1,2), 'g-')
    legend('b=mat v. dl')
    xlabel('mat h2o conc')
    ylabel('dl h2o conc')
    
    f3='h2o conc mat_goes';
    figure(3)
    plot(D(42,:), D(212,:), 'r.')
    hold on
    plot(D(42,:), P2(1,1) * D(42,:) + P2(1,2), 'g-')
    legend('r=mat v. goes')
    xlabel('h2o conc (W/m^2)')
    ylabel('h2o conc (W/m^2)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    saveas(figure(3), [figout f3]);
    close(figure(1),figure(2),figure(3))
    % ------------------------------------------------------------------%
    %Fco2
    
    %dl
    isn=isnan(D(46,:)) ==1 | isnan(D(345,:)) == 1;
    reg = D(345,:) <= 0 & D(46,:) <= 0 & turb ==1 & ~isn;
    A=D(46,reg)';
    B=D(345,reg)';
    P11=A\B % matlab gives you the least squares solution to the A P = B
    neg1 = D(46,:) <= 0;
    B11 = robustfit(D(46,reg),D(345,reg))
    
    reg = D(345,:) > 0 & D(46,:) > 0 & turb ==1 & ~isn;
    A=D(46,reg)';
    B=D(345,reg)';
    P12=A\B % matlab gives you the least squares solution to the A P = B
    pos1 = D(46,:) > 0;
    B12 = robustfit(D(46,reg),D(345,reg))

    %goes
    isn=isnan(D(46,:)) ==1 | isnan(D(346,:)) == 1;
    reg = D(346,:) <= 0 & D(46,:) <= 0 & turb ==1  & ~isn;
    A=D(46,reg)';
    B=D(346,reg)';
    P21=A\B % matlab gives you the least squares solution to the A P = B
    neg2 = D(46,:) <= 0;
    B21 = robustfit(D(46,reg),D(346,reg))
    
    reg = D(346,:) > 0 & D(46,:) > 0 & turb ==1  & ~isn;
    A=D(46,reg)';
    B=D(346,reg)';
    P22=A\B % matlab gives you the least squares solution to the A P = B
    pos2 = D(46,:) > 0;
    B22 = robustfit(D(46,reg),D(346,reg))
    
    f='Fco2 (micromol_m-2_s-1)';
    figure(1)
    %fast
    plot(D(1,:), D(46,:), 'ko')
    hold on

    %dl T
    plot(D(1,:), D(345,:), 'bx')
    hold on

    %goes T
    plot(D(1,:), D(346,:), 'r.')
    legend('k = mat', 'b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('Fco2 (micromol/m2/s)')

    f2='Fco2 (micromol_m-2_s-1) mat_dl';
    figure(2)
    plot(D(46,:), D(345,:), 'b.')
    hold on
    plot(D(46,neg1), P11(1,1) * D(46,neg1), 'g-')
    hold on
    plot(D(46,neg1), B11(2,1) * D(46,neg1), 'k-')
    hold on
    plot(D(46,pos1), P12(1,1) * D(46,pos1), 'g-')
    hold on
    plot(D(46,pos1), B12(2,1) * D(46,pos1), 'k-')
    legend('b=mat v. dl', 'least square forced through 0', 'robustfit slope only')
    xlabel('mat Fco2 (micromol/m2/s)')
    ylabel('dl Fco2 (micromol/m2/s)')
    
    f3='Fco2 (micromol_m-2_s-1) mat_goes';
    figure(3)
    plot(D(46,:), D(346,:), 'b.')
    hold on
    plot(D(46,neg2), P21(1,1) * D(46,neg2), 'g-')
    hold on
    plot(D(46,neg2), B21(2,1) * D(46,neg2), 'k-')
    hold on
    plot(D(46,pos2), P22(1,1) * D(46,pos2), 'g-')
    hold on
    plot(D(46,pos2), B22(2,1) * D(46,pos2), 'k-')
    legend('b=mat v. dl', 'least square forced through 0', 'robustfit slope only')
    xlabel('mat Fco2 (micromol/m2/s)')
    ylabel('dl Fco2 (micromol/m2/s)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    saveas(figure(3), [figout f3]);
    close(figure(1),figure(2),figure(3))   
    % ------------------------------------------------------------------%
    %Fh2o   
    disp('P1 is the slope and intercept of mat and dl')
    disp('P2 is the slope and intercept of mat and goes')
    
    isn=isnan(D(47,:)) ==1 | isnan(D(347,:)) == 1;
    P1=polyfit(D(47,~isn),D(347,~isn),1)
    B1 = robustfit(D(47,~isn),D(347,~isn))
    
    isn=isnan(D(47,:)) ==1 | isnan(D(348,:)) == 1;
    P2=polyfit(D(47,~isn),D(348,~isn),1)
    B2 = robustfit(D(47,~isn),D(348,~isn))
    
    %dl
    f='Fh2o (millimol h2o_m-2_s-1)';
    figure(1)
    %fast
    plot(D(1,:), D(47,:), 'ko')
    hold on

    %dl T
    plot(D(1,:), D(347,:), 'bx')
    hold on

    %goes T
    plot(D(1,:), D(348,:), 'r.')
    legend('k = mat', 'b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('Fh2o (millimol h2o/m2/s)')

    f2='Fh2o (millimol h2o_m-2_s-1) mat_dl';
    figure(2)
    plot(D(47,:), D(347,:), 'b.')
    hold on
    plot(D(47,:), P1(1,1) * D(47,:) + P1(1,2), 'g-')
    hold on
    plot(D(47,:), B1(2,1) * D(47,:) + B1(1,1), 'k-')
    legend('b=mat v. dl', 'polyfit', 'robustfit')
    xlabel('mat Fh2o (millimol h2o/m2/s)')
    ylabel('dl Fh2o (millimol h2o/m2/s)')
    
    f3='Fh2o (millimol h2o_m-2_s-1) mat_goes';
    figure(3)
    plot(D(47,:), D(348,:), 'b.')
    hold on
    plot(D(47,:), P2(1,1) * D(47,:) + P2(1,2), 'g-')
    hold on
    plot(D(47,:), B2(2,1) * D(47,:) + B2(1,1), 'k-')
    legend('b=mat v. dl', 'polyfit', 'robustfit')
    xlabel('mat Fh2o (millimol h2o/m2/s)')
    ylabel('dl Fh2o (millimol h2o/m2/s)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    saveas(figure(3), [figout f3]);
    close(figure(1),figure(2),figure(3)) 
    % ------------------------------------------------------------------%
    %Rn  
    %dl
    figure(1)
    plot(D(1,:), D(87,:), 'bx')
    hold on

    %goes T
    f='Rn (W_m^-2)';
    plot(D(1,:), D(213,:), 'r.')
    legend('b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('Rn (W/m^2)')

    f2='Rn (W_m^-2) dl_goes';
    figure(2)
    plot(D(87,:), D(213,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl Rn (W/m^2)')
    ylabel('goes Rn (W/m^2)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2))  
    
    % ------------------------------------------------------------------%
    %PAR_In  
    %dl
    f='PAR_In (micromol photons_m-2_s-1)';
    figure(1)
    plot(D(1,:), D(90,:), 'bx')
    hold on

    %goes T
    plot(D(1,:), D(216,:), 'r.')
    legend('b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('PAR_In (micromol photons/m2/s)')

    f2='PAR_In (micromol photons_m-2_s-1) dl_goes';
    figure(2)
    plot(D(90,:), D(216,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl PAR_In (micromol photons/m2/s)')
    ylabel('goes PAR_In (micromol photons/m2/s)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2))  
    
    % ------------------------------------------------------------------%
    %PAR_Out  
    %dl
    f = 'PAR_Out (micromol photons_m-2_s-1)';
    figure(1)
    plot(D(1,:), D(91,:), 'bx')
    hold on

    %goes T
    plot(D(1,:), D(217,:), 'r.')
    legend('b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('PAR_Out (micromol photons/m2/s)')

    f2 = 'PAR_Out (micromol photons_m-2_s-1) dl_goes';
    figure(2)
    plot(D(91,:), D(217,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl PAR_Out (micromol photons/m2/s)')
    ylabel('goes PAR_Out (micromol photons/m2/s)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2)) 
    
    % ------------------------------------------------------------------%
    %SOLAR_IN  
    %dl
    f = 'K (W_m-2)';
    figure(1)
    plot(D(1,:), D(88,:), 'bx')
    hold on

    %goes 
    plot(D(1,:), D(214,:), 'r.')
    legend('b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('K (W/m2)')

    f2 = 'K (W_m-2) dl_goes';
    figure(2)
    plot(D(88,:), D(214,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl K (W/m2)')
    ylabel('goes K (W/m2)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2)) 
    
    % ------------------------------------------------------------------%
    %SOLAR_OUT 
    %dl
    f= 'K out (W_m-2)';
    figure(1)
    plot(D(1,:), D(89,:), 'bx')
    hold on

    %goes 
    plot(D(1,:), D(215,:), 'r.')
    legend('b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('K out (W/m2)')

    f2='K out (W_m-2) dl_goes';
    figure(2)
    plot(D(89,:), D(215,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl K out (W/m2)')
    ylabel('goes K out (W/m2)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2)) 
    
    % ------------------------------------------------------------------%
    %compare radiation
    disp('compare different radiation sensors - did not save these plots')
    
    
    figure(1)
    plot(D(88,:), D(87,:), 'k.') 
    xlabel('K in')
    ylabel('Rn')
    
    pause
    f='Kin v Rn dl';
    saveas(figure(1), [figout f]);
    close(figure(1))
    
    figure(1)
    plot(D(88,:), D(213,:), 'k.')
    xlabel('K in')
    ylabel('Rn')

    pause
    f='Kin v Rn goes';
    saveas(figure(1), [figout f]);
    close(figure(1))
    
    figure(1)
    plot(D(88,:), D(214,:), 'k.') 
    xlabel('K in')
    ylabel('goes pyrr in')

    pause
    f='Kin v goes pyrr in';
    saveas(figure(1), [figout f]);
    close(figure(1))
    
    figure(1)
    plot(D(88,:), D(89,:), 'k.') 
    xlabel('K in')
    ylabel('K out')

    pause
    f='Kin v K out dl';
    saveas(figure(1), [figout f]);
    close(figure(1))
    
    
    figure(1)
    plot(D(88,:), D(215,:), 'k.') 
    xlabel('K in')
    ylabel('K out')

    pause
    f='Kin v K out goes';
    saveas(figure(1), [figout f]);
    close(figure(1))
    
    figure(1)
    plot(D(88,:), D(90,:), 'k.') 
    xlabel('K in')
    ylabel('PAR in')

    pause
    f='Kin v PAR in dl';
    saveas(figure(1), [figout f]);
    close(figure(1))
    
    figure(1)
    plot(D(88,:), D(216,:), 'k.') 
    xlabel('K in')
    ylabel('PAR in')

    pause
    f='Kin v PAR in goes';
    saveas(figure(1), [figout f]);
    close(figure(1))
    
    figure(1)
    plot(D(88,:), D(91,:), 'k.') 
    xlabel('K in')
    ylabel('PAR out')

    pause
    f='Kin v PAR out dl';
    saveas(figure(1), [figout f]);
    close(figure(1))
    
    figure(1)
    plot(D(88,:), D(217,:), 'k.') 
    xlabel('K in')
    ylabel('PAR out')

    pause
    f='Kin v PAR out goes';
    saveas(figure(1), [figout f]);
    close(figure(1))
    % ------------------------------------------------------------------%
    %HMP_Temp 
    %dl
    f='HMP_Temp (W_m-2)';
    figure(1)
    plot(D(1,:), D(97,:), 'bx')
    hold on

    %goes 
    plot(D(1,:), D(223,:), 'r.')
    legend('b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('HMP_Temp (W/m2)')

    
    f2='HMP_Temp (W_m-2) dl_goes';
    figure(2)
    plot(D(97,:), D(223,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl HMP_Temp (C)')
    ylabel('goes HMP_Temp (C)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2))
    
    % ------------------------------------------------------------------%
    %HMP_RH
    %dl
    f= 'HMP_RH';
    figure(1)
    plot(D(1,:), D(98,:), 'bx')
    hold on

    %goes 
    plot(D(1,:), D(224,:), 'r.')
    legend('b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('HMP_RH')

    f2= 'HMP_RH dl_goes';
    figure(2)
    plot(D(98,:), D(224,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl HMP_RH')
    ylabel('goes HMP_RH')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2))
    
    % ------------------------------------------------------------------%
    %T_107
    %dl
    f='T_107';
    figure(1)
    plot(D(1,:), D(102,:), 'bx')
    hold on

    %goes 
    plot(D(1,:), D(228,:), 'r.')
    legend('b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('T_107')

    f2='T_107 dl_goes';
    figure(2)
    plot(D(102,:), D(228,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl T_107')
    ylabel('goes T_107')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2))
    
    %------------------------------------------------------------------%
    %rain
    %dl
    f1='rain (mm)';
    figure(1)
    plot(D(1,:), D(103,:), 'bx')
    hold on

    %goes 
    plot(D(1,:), D(229,:), 'r.')
    legend('b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('rain (mm)')

    f2='rain (mm) dl_goes';
    figure(2)
    plot(D(103,:), D(229,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl rain (mm)')
    ylabel('goes rain (mm)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2))
    
    % ------------------------------------------------------------------%
    %h2o conc
    %dl
    f='hmp h2o conc (pp thousand)';
    figure(1)
    plot(D(1,:), D(100,:), 'bx')
    hold on

    %goes 
    plot(D(1,:), D(226,:), 'r.')
    legend('b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('hmp h2o conc (pp thousand)')

    f2='hmp h2o conc (pp thousand) dl_goes';
    figure(2)
    plot(D(100,:), D(226,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl hmp h2o conc (pp thousand)')
    ylabel('goes hmp h2o conc (pp thousand)')
    
    pause
    
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2))
    
    % ------------------------------------------------------------------% 
    %cross check
    figure(1)
    plot(D(100,:), D(42,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl hmp h2o conc (pp thousand)')
    ylabel('goes irga h2o conc (pp thousand)')
   
    pause 
    f='irga h2o conc (pp thousand) hmp_mat';
    saveas(figure(1), [figout f]);
    close(figure(1))
    
    figure(1)
    plot(D(100,:), D(86,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl hmp h2o conc (pp thousand)')
    ylabel('goes irga h2o conc (pp thousand)')
   
    pause
    f='irga h2o conc (pp thousand) hmp_dl';
    saveas(figure(1), [figout f]);
    close(figure(1))
    
    figure(1)
    plot(D(100,:), D(212,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl hmp h2o conc (pp thousand)')
    ylabel('goes irga h2o conc (pp thousand)')
   
    pause
    f='irga h2o conc (pp thousand) hmp_goes';
    saveas(figure(1), [figout f]);
    close(figure(1))
    % ------------------------------------------------------------------% 
    %NDVI
    %dl
    figure(1)
    plot(D(1,:), D(96,:), 'bx')
    hold on

    %goes 
    plot(D(1,:), D(222,:), 'r.')
    legend('b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('NDVI')

    figure(2)
    plot(D(96,:), D(222,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl NDVI')
    ylabel('goes NDVI')
    
    pause
    
    f='NDVI';
    f2='NDVI dl_goes';
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2))
    
    % ------------------------------------------------------------------%
    
    %RED IN
    %dl
    figure(1)
    plot(D(1,:), D(92,:), 'bx')
    hold on

    %goes 
    plot(D(1,:), D(218,:), 'r.')
    legend('b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('RED IN')

    figure(2)
    plot(D(92,:), D(218,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl RED IN')
    ylabel('goes RED IN')
    
    pause
    
    f='RED IN';
    f2='RED IN dl_goes';
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2))
    
    % ------------------------------------------------------------------%
    
    %RED OUT
    %dl
    figure(1)
    plot(D(1,:), D(94,:), 'bx')
    hold on

    %goes 
    plot(D(1,:), D(220,:), 'r.')
    legend('b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('RED IN')

    figure(2)
    plot(D(94,:), D(220,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl RED OUT')
    ylabel('goes RED OUT')
    
    pause
    
    f='RED OUT';
    f2='RED OUT dl_goes';
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2))
    
    % ------------------------------------------------------------------%
    
    %NIR In
    %dl
    figure(1)
    plot(D(1,:), D(93,:), 'bx')
    hold on

    %goes 
    plot(D(1,:), D(219,:), 'r.')
    legend('b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('NIR In')

    figure(2)
    plot(D(93,:), D(219,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl NIR In')
    ylabel('goes NIR In')
    
    pause
    
    f='NIR In';
    f2='NIR In dl_goes';
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2))
    
    % ------------------------------------------------------------------%
    
    %NIR Out
    %dl
    figure(1)
    plot(D(1,:), D(95,:), 'bx')
    hold on

    %goes 
    plot(D(1,:), D(221,:), 'r.')
    legend('b=dl', 'r = goes')
    xlabel('tower time (day)')
    ylabel('NIR Out')

    figure(2)
    plot(D(95,:), D(221,:), 'b.')
    legend('b=dl v. goes')
    xlabel('dl NIR Out')
    ylabel('goes NIR Out')
    
    pause
    
    f='NIR Out';
    f2='NIR Out dl_goes';
    saveas(figure(1), [figout f]);
    saveas(figure(2), [figout f2]);
    close(figure(1),figure(2))
      
    
    %------------------------------------------------------------------%
    %soil moisture    
    %dl
    figure(1)
    plot(D(1,:), D(282,:), 'rx')
    hold on
    plot(D(1,:), D(283,:), 'bx')
    hold on
    plot(D(1,:), D(284,:), 'gx')
    hold on
    plot(D(1,:), D(285,:), 'kx')
    xlabel('tower time (day)')
    ylabel('fraction Soil Moisture')
    
    pause
    
    f='SoilMoisture';
    saveas(figure(1), [figout f]);
    close(figure(1))
    
    %------------------------------------------------------------------%
    %soil temp   
    %dl
    figure(1)
    plot(D(1,:), D(288,:), 'rx')
    hold on
    plot(D(1,:), D(289,:), 'bx')
    hold on
    plot(D(1,:), D(290,:), 'gx')
    hold on
    plot(D(1,:), D(291,:), 'kx')
    xlabel('tower time (day)')
    ylabel('Soil Temperature')
    
    pause
    
    f='Soil Temperature';
    saveas(figure(1), [figout f]);
    close(figure(1))
    
    %------------------------------------------------------------------%
    %LWS 
    %dl
    figure(1)
    plot(D(1,:), D(293,:), 'rx')
    hold on
    plot(D(1,:), D(294,:), 'bx')
    hold on
    plot(D(1,:), D(295,:), 'gx')
    xlabel('tower time (day)')
    ylabel('LWS')
    
    pause
    
    f='LWS';
    saveas(figure(1), [figout f]);
    close(figure(1))
    
    %------------------------------------------------------------------%
    %Matric Potential 
    %dl
    figure(1)
    plot(D(1,:), D(296,:), 'rx')
    hold on
    plot(D(1,:), D(297,:), 'bx')
    hold on
    plot(D(1,:), D(298,:), 'gx')
    hold on
    plot(D(1,:), D(299,:), 'mx')
    hold on
    plot(D(1,:), D(300,:), 'yp')
    hold on
    plot(D(1,:), D(301,:), 'kx')
    legend('depth1', 'depth2', 'depth3', 'depth4', 'depth5', 'depth6')
    xlabel('tower time (day)')
    ylabel('Matric Start T')
    
    pause
    
    f='Matric Start T';
    saveas(figure(1), [figout f]);
    close(figure(1))
    
    %------------------------------------------------------------------%
    %DelT
    %dl
    figure(1)
    plot(D(1,:), D(314,:), 'rx')
    hold on
    plot(D(1,:), D(315,:), 'bx')
    hold on
    plot(D(1,:), D(316,:), 'gx')
    hold on
    plot(D(1,:), D(317,:), 'mx')
    hold on
    plot(D(1,:), D(318,:), 'yp')
    hold on
    plot(D(1,:), D(319,:), 'kx')
    legend('depth1', 'depth2', 'depth3', 'depth4', 'depth5', 'depth6')
    xlabel('tower time (day)')
    ylabel('DelT')
    
    pause
    
    f='DelT';
    saveas(figure(1), [figout f]);
    close(figure(1))
    %------------------------------------------------------------------%
    
else
    disp('broken combinedinfo')
    N_obs =[];
    headerN_obs =[];   
end