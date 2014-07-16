function [D, HMERGE, Cleaned_D, DataSource, Cleaned_Header] = shorthairCombined(D, HMERGE, Cleaned_D, DataSource, Cleaned_Header)


    %focus periods - 
    focus1 = D(1,:) < 735.4376; 
    focus2 = D(1,:) > 735.4376; 
    
    %----------------------------------------------------------------------
    %1.Time
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(1,:)= D(1,:);
    
    %Fil Header
    Cleaned_Header(1,1)={'TIME'};
    
    %----------------------------------------------------------------------
    %2.Tsonic
    %----------------------------------------------------------------------
    
    %Fill Data
    Cleaned_D(2,:)=D(8,:);%fast
    bad=isnan(Cleaned_D(2,:));
    Cleaned_D(2,bad)=D(84,bad);%dl
    bad=isnan(Cleaned_D(2,:));
    Cleaned_D(2,bad)=D(210,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(8,:));
    good_dl=(isnan(D(8,:)) & ~isnan(D(84,:)));
    good_goes=(isnan(D(8,:)) & isnan(D(84,:)) & ~isnan(D(210,:)));
    DataSource(2,good_mat)=1;
    DataSource(2,good_dl)=2;
    DataSource(2,good_goes)=3;
    
    %Fill Header
    Cleaned_Header(2,1)={'TSONIC'};
    
    %----------------------------------------------------------------------
    %3.Tactual
    %----------------------------------------------------------------------
    %****************************************
    %Fill Data
    Cleaned_D(3,:)=D(330,:);%fast
    bad=isnan(Cleaned_D(3,:));
    Cleaned_D(3,bad)=D(331,bad);%dl
    bad=isnan(Cleaned_D(3,:));
    Cleaned_D(3,bad)=D(332,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(330,:));
    good_dl=(isnan(D(330,:)) & ~isnan(D(331,:)));
    good_goes=(isnan(D(330,:)) & isnan(D(331,:)) & ~isnan(D(332,:)));
    DataSource(3,good_mat)=1;
    DataSource(3,good_dl)=2;
    DataSource(3,good_goes)=3;
    
    
    %Fill Header
    Cleaned_Header(3,1)={'TACTUAL'};
    %****************************************
    %----------------------------------------------------------------------
    %4.Ustar
    %----------------------------------------------------------------------
    %homogenize
    disp('homogenize Ustar with robustfit - shorthair using 2 periods')
    
    %----------------------------------------------------------------------
    %dl
    %----------------------------------------------------------------------
    %period 1
    isn=isnan(D(21,:)) ==1 | isnan(D(322,:)) == 1 | focus2(1,:) ==1; % we do not want period 2
    %P=polyfit(D(21,~isn),D(322,~isn),1);
    %Dadj=(D(322,:)-P(1,2))./P(1,1);
    B = robustfit(D(21,~isn),D(322,~isn)); % B1(2,1) is the slope and B1(1,1) is the intercept
    Dadj=(D(322,:)-B(1,1))./B(2,1);
    D(322,focus1)=Dadj(1, focus1);
    %----------------------------------------------------------------------
    %period 2
    isn=isnan(D(21,:)) ==1 | isnan(D(322,:)) == 1 | focus1(1,:) ==1; % we do not want period 2
    %P=polyfit(D(21,~isn),D(322,~isn),1);
    %Dadj=(D(322,:)-P(1,2))./P(1,1);
    B = robustfit(D(21,~isn),D(322,~isn)); % B1(2,1) is the slope and B1(1,1) is the intercept
    Dadj=(D(322,:)-B(1,1))./B(2,1);
    D(322,focus2)=Dadj(1, focus2);
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    %goes
    %----------------------------------------------------------------------
    %period 1
    isn=isnan(D(21,:)) ==1 | isnan(D(323,:)) == 1 | focus2(1,:) ==1;
    %P=polyfit(D(21,~isn),D(323,~isn),1);
    %Dadj=(D(323,:)-P(1,2))./P(1,1);
    B = robustfit(D(21,~isn),D(323,~isn)); % B1(2,1) is the slope and B1(1,1) is the intercept
    Dadj=(D(323,:)-B(1,1))./B(2,1);
    D(323,focus1)=Dadj(1,focus1);
    %----------------------------------------------------------------------
    %period 2
    isn=isnan(D(21,:)) ==1 | isnan(D(323,:)) == 1 | focus1(1,:) ==1;
    %P=polyfit(D(21,~isn),D(323,~isn),1);
    %Dadj=(D(323,:)-P(1,2))./P(1,1);
    B = robustfit(D(21,~isn),D(323,~isn)); % B1(2,1) is the slope and B1(1,1) is the intercept
    Dadj=(D(323,:)-B(1,1))./B(2,1);
    D(323,focus2)=Dadj(1,focus2);
    %----------------------------------------------------------------------
    
    %****************************************
    %Fill Data
    Cleaned_D(4,:)=D(21,:);
    bad=isnan(Cleaned_D(4,:));
    Cleaned_D(4,bad)=D(322,bad);
    bad=isnan(Cleaned_D(4,:));
    Cleaned_D(4,bad)=D(323,bad);
    
    %Where is the data from?
    good_mat=~isnan(D(21,:));
    good_dl=(isnan(D(21,:)) & ~isnan(D(322,:)));
    good_goes=(isnan(D(21,:)) & isnan(D(322,:)) & ~isnan(D(323,:)));
    DataSource(4,good_mat)=1;
    DataSource(4,good_dl)=2;
    DataSource(4,good_goes)=3;
    
    %put USTAR source into Cleaned_D
    Cleaned_D(25,:)=DataSource(4,:);
    bad=Cleaned_D(25,:)==0;
    Cleaned_D(25,bad)=NaN;
  
    %Fill Header
    Cleaned_Header(4,1)={'USTAR'};
    Cleaned_Header(25,1)={'USTAR_SOURCE_FLAG'};
    %****************************************
    %----------------------------------------------------------------------
    %5.Wind Direction
    %----------------------------------------------------------------------
    %****************************************
    %Fill Data
    Cleaned_D(5,:)=D(324,:);%fast
    bad=isnan(Cleaned_D(5,:));
    Cleaned_D(5,bad)=D(325,bad);%dl
    bad=isnan(Cleaned_D(5,:));
    Cleaned_D(5,bad)=D(326,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(324,:));
    good_dl=(isnan(D(324,:)) & ~isnan(D(325,:)));
    good_goes=(isnan(D(324,:)) & isnan(D(325,:)) & ~isnan(D(326,:)));
    DataSource(5,good_mat)=1;
    DataSource(5,good_dl)=2;
    DataSource(5,good_goes)=3;
    
    %Fill Header
    Cleaned_Header(5,1)={'WDIR'};
    %****************************************
    %----------------------------------------------------------------------
    %6.Wind Speed
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    %homogenize
    isn=isnan(D(327,:)) ==1 | isnan(D(328,:)) == 1 | focus2(1,:) ==1;
    P=polyfit(D(327,~isn),D(328,~isn),1);
    Dadj=(D(328,:)-P(1,2))./P(1,1);
    D(328,focus1)=Dadj(1,focus1);
    %----------------------------------------------------------------------
    isn=isnan(D(327,:)) ==1 | isnan(D(328,:)) == 1 | focus1(1,:) ==1;
    P=polyfit(D(327,~isn),D(328,~isn),1);
    Dadj=(D(328,:)-P(1,2))./P(1,1);
    D(328,focus2)=Dadj(1,focus2);
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    isn=isnan(D(327,:)) ==1 | isnan(D(329,:)) == 1 | focus2(1,:) ==1;
    P=polyfit(D(327,~isn),D(329,~isn),1);
    Dadj=(D(329,:)-P(1,2))./P(1,1);
    D(329,focus1)=Dadj(1,focus1);
    %----------------------------------------------------------------------
    isn=isnan(D(327,:)) ==1 | isnan(D(329,:)) == 1 | focus1(1,:) ==1;
    P=polyfit(D(327,~isn),D(329,~isn),1);
    Dadj=(D(329,:)-P(1,2))./P(1,1);
    D(329,focus2)=Dadj(1,focus2);
    %----------------------------------------------------------------------
    %****************************************
    %Fill Data
    Cleaned_D(6,:)=D(327,:);%fast
    bad=isnan(Cleaned_D(6,:));
    Cleaned_D(6,bad)=D(328,bad);%dl
    bad=isnan(Cleaned_D(6,:));
    Cleaned_D(6,bad)=D(329,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(327,:));
    good_dl=(isnan(D(327,:)) & ~isnan(D(328,:)));
    good_goes=(isnan(D(327,:)) & isnan(D(328,:)) & ~isnan(D(329,:)));
    DataSource(6,good_mat)=1;
    DataSource(6,good_dl)=2;
    DataSource(6,good_goes)=3;
    
    %Fill Header
    Cleaned_Header(6,1)={'WSPD'};
    %****************************************
    %----------------------------------------------------------------------
    %7. Dry Density
    %----------------------------------------------------------------------
    %homogenize
    isn=isnan(D(336,:)) ==1 | isnan(D(337,:)) == 1;
    P=polyfit(D(336,~isn),D(337,~isn),1);
    Dadj=(D(337,:)-P(1,2))./P(1,1);
    D(337,:)=Dadj;
    
    isn=isnan(D(336,:)) ==1 | isnan(D(338,:)) == 1;
    P=polyfit(D(336,~isn),D(338,~isn),1);
    Dadj=(D(338,:)-P(1,2))./P(1,1);
    D(338,:)=Dadj;
    %****************************************
    %Fill Data
    Cleaned_D(7,:)=D(336,:);%fast
    bad=isnan(Cleaned_D(7,:));
    Cleaned_D(7,bad)=D(337,bad);%dl
    bad=isnan(Cleaned_D(7,:));
    Cleaned_D(7,bad)=D(338,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(336,:));
    good_dl=(isnan(D(336,:)) & ~isnan(D(337,:)));
    good_goes=(isnan(D(336,:)) & isnan(D(337,:)) & ~isnan(D(338,:)));
    DataSource(7,good_mat)=1;
    DataSource(7,good_dl)=2;
    DataSource(7,good_goes)=3;
    
    %Fill Header
    Cleaned_Header(7,1)={'DENSITY_DRY'};
    %****************************************
    %----------------------------------------------------------------------
    %8. Moist Density
    %----------------------------------------------------------------------
    %homogenize
    isn=isnan(D(333,:)) ==1 | isnan(D(334,:)) == 1;
    P=polyfit(D(333,~isn),D(334,~isn),1);
    Dadj=(D(334,:)-P(1,2))./P(1,1);
    D(334,:)=Dadj;
    
    isn=isnan(D(333,:)) ==1 | isnan(D(335,:)) == 1;
    P=polyfit(D(333,~isn),D(335,~isn),1);
    Dadj=(D(335,:)-P(1,2))./P(1,1);
    D(335,:)=Dadj;
    %****************************************
    %Fill Data
    Cleaned_D(8,:)=D(333,:);%fast
    bad=isnan(Cleaned_D(8,:));
    Cleaned_D(8,bad)=D(334,bad);%dl
    bad=isnan(Cleaned_D(8,:));
    Cleaned_D(8,bad)=D(335,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(333,:));
    good_dl=(isnan(D(333,:)) & ~isnan(D(334,:)));
    good_goes=(isnan(D(333,:)) & isnan(D(334,:)) & ~isnan(D(335,:)));
    DataSource(8,good_mat)=1;
    DataSource(8,good_dl)=2;
    DataSource(8,good_goes)=3;
    
    %Fill Header
    Cleaned_Header(8,1)={'DENSITY_MOIST'};
    %----------------------------------------------------------------------
    %9. Sensible heat
    %----------------------------------------------------------------------
    %dl
    %----------------------------------------------------------------------
    %homogenize
    isn=isnan(D(29,:)) ==1 | isnan(D(339,:)) == 1 | focus2(1,:) ==1;
    P=polyfit(D(29,~isn),D(339,~isn),1);
    Dadj=(D(339,:)-P(1,2))./P(1,1);
    D(339,focus1)=Dadj(1,focus1);
    
    isn=isnan(D(29,:)) ==1 | isnan(D(339,:)) == 1 | focus1(1,:) ==1;
    P=polyfit(D(29,~isn),D(339,~isn),1);
    Dadj=(D(339,:)-P(1,2))./P(1,1);
    D(339,focus2)=Dadj(1, focus2);
    
    %----------------------------------------------------------------------
    isn=isnan(D(29,:)) ==1 | isnan(D(340,:)) == 1 | focus2(1,:) ==1;
    P=polyfit(D(29,~isn),D(340,~isn),1);
    Dadj=(D(340,:)-P(1,2))./P(1,1);
    D(340,focus1)=Dadj(1, focus1);
    
    isn=isnan(D(29,:)) ==1 | isnan(D(340,:)) == 1 | focus1(1,:) ==1;
    P=polyfit(D(29,~isn),D(340,~isn),1);
    Dadj=(D(340,:)-P(1,2))./P(1,1);
    D(340,focus2)=Dadj(1, focus2);
    %----------------------------------------------------------------------
    %****************************************
    
    %Fill Data
    Cleaned_D(9,:)=D(29,:);%fast
    bad=isnan(Cleaned_D(9,:));
    Cleaned_D(9,bad)=D(339,bad);%dl
    bad=isnan(Cleaned_D(9,:));
    Cleaned_D(9,bad)=D(340,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(29,:));
    good_dl=(isnan(D(29,:)) & ~isnan(D(339,:)));
    good_goes=(isnan(D(29,:)) & isnan(D(339,:)) & ~isnan(D(340,:)));
    DataSource(9,good_mat)=1;
    DataSource(9,good_dl)=2;
    DataSource(9,good_goes)=3;
    
    %put source into Cleaned_D
    Cleaned_D(26,:)=DataSource(9,:);
    bad=Cleaned_D(26,:)==0;
    Cleaned_D(26,bad)=NaN;
    
    %Fill Header
    Cleaned_Header(9,1)={'SENSIBLE_HEAT'};
    Cleaned_Header(26,1)={'SENSIBLE_HEAT_SOURCE_FLAG'};
   
    %----------------------------------------------------------------------
    %11. CO2 concentration
    %----------------------------------------------------------------------
    %homogenize
    isn=isnan(D(35,:)) ==1 | isnan(D(85,:)) == 1;
    P=polyfit(D(35,~isn),D(85,~isn),1);
    Dadj=(D(85,:)-P(1,2))./P(1,1);
    D(85,:)=Dadj;
    
    isn=isnan(D(35,:)) ==1 | isnan(D(211,:)) == 1;
    P=polyfit(D(35,~isn),D(211,~isn),1);
    Dadj=(D(211,:)-P(1,2))./P(1,1);
    D(211,:)=Dadj;
    %****************************************
    
    %Fill Data
    Cleaned_D(11,:)=D(35,:);%fast
    bad=isnan(Cleaned_D(11,:));
    Cleaned_D(11,bad)=D(85,bad);%dl
    bad=isnan(Cleaned_D(11,:));
    Cleaned_D(11,bad)=D(211,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(35,:));
    good_dl=(isnan(D(35,:)) & ~isnan(D(85,:)));
    good_goes=(isnan(D(35,:)) & isnan(D(85,:)) & ~isnan(D(211,:)));
    DataSource(11,good_mat)=1;
    DataSource(11,good_dl)=2;
    DataSource(11,good_goes)=3;
    
    %Fill Header
    Cleaned_Header(11,1)={'CO2_CONC'};
    %----------------------------------------------------------------------
    %12. H2O concentration - mixing ratio
    %----------------------------------------------------------------------
    %homogenize
    isn=isnan(D(42,:)) ==1 | isnan(D(86,:)) == 1;
    P=polyfit(D(42,~isn),D(86,~isn),1);
    Dadj=(D(86,:)-P(1,2))./P(1,1);
    D(86,:)=Dadj;
    
    isn=isnan(D(42,:)) ==1 | isnan(D(212,:)) == 1;
    P=polyfit(D(42,~isn),D(212,~isn),1);
    Dadj=(D(212,:)-P(1,2))./P(1,1);
    D(212,:)=Dadj;
    %****************************************
    
    %Fill Data
    Cleaned_D(12,:)=D(42,:);%fast
    bad=isnan(Cleaned_D(12,:));
    Cleaned_D(12,bad)=D(86,bad);%dl
    bad=isnan(Cleaned_D(12,:));
    Cleaned_D(12,bad)=D(212,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(42,:));
    good_dl=(isnan(D(42,:)) & ~isnan(D(86,:)));
    good_goes=(isnan(D(42,:)) & isnan(D(86,:)) & ~isnan(D(212,:)));
    DataSource(12,good_mat)=1;
    DataSource(12,good_dl)=2;
    DataSource(12,good_goes)=3;
    
    %Fill Header
    Cleaned_Header(12,1)={'H2O_CONC'};
    %----------------------------------------------------------------------
    %13. FCO2
    %----------------------------------------------------------------------
    %split into neg fluxes and pos fluxes before adjusting the gain
      
    %we will use just ust filtered observations 
    turb = Cleaned_D(4,:)>0.3;
    
    %homogenize
    Dadj(1,:)=D(345,:).*NaN;
    isn=isnan(D(46,:)) ==1 | isnan(D(345,:)) == 1;
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %dl
    %----------------------------------------------------------------------
    %period 1 & neg Fco2
    reg = D(345,:) <= 0 & D(46,:) <= 0 & turb ==1 & ~isn & focus1(1,:) ==1;
    %P=polyfit(D(46,reg),D(345,reg),1);
    %least squares regression with intercept through zero
    A=D(46,reg)';
    B=D(345,reg)';
    P=A\B; % matlab gives you the least squares solution to the A P = B
    neg = D(345,:) <= 0 & focus1(1,:) ==1;
    Dadj(1,neg)=D(345,neg)/P(1,1);
    
    %period 1 & pos Fco2
    reg = D(345,:) > 0 & D(46,:) > 0 & turb ==1 & ~isn & focus1(1,:) ==1;
    %P=polyfit(D(46,reg),D(345,reg),1);
    A=D(46,reg)';
    B=D(345,reg)';
    P=A\B; % matlab gives you the least squares solution to the A P = B
    pos = D(345,:) > 0 & focus1(1,:) ==1;
    Dadj(1,pos)=D(345,pos)./P(1,1);

    %----------------------------------------------------------------------
    %period 2 & neg Fco2
    reg = D(345,:) <= 0 & D(46,:) <= 0 & turb ==1 & ~isn & focus2(1,:) ==1;
    %P=polyfit(D(46,reg),D(345,reg),1);
    %least squares regression with intercept through zero
    A=D(46,reg)';
    B=D(345,reg)';
    P=A\B; % matlab gives you the least squares solution to the A P = B
    neg = D(345,:) <= 0 & focus2(1,:) ==1;
    Dadj(1,neg)=D(345,neg)/P(1,1);
    
    %period 2 & pos Fco2
    reg = D(345,:) > 0 & D(46,:) > 0 & turb ==1 & ~isn & focus2(1,:) ==1;
    %P=polyfit(D(46,reg),D(345,reg),1);
    A=D(46,reg)';
    B=D(345,reg)';
    P=A\B; % matlab gives you the least squares solution to the A P = B
    pos = D(345,:) > 0 & focus2(1,:) ==1;
    Dadj(1,pos)=D(345,pos)./P(1,1);
    
    D(345,:)=Dadj;
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %goes
    %----------------------------------------------------------------------
    Dadj(1,:)=D(346,:).*NaN;
    isn=isnan(D(46,:)) ==1 | isnan(D(346,:)) == 1;
    
    
    reg = D(346,:) <= 0 & D(46,:) <= 0 & turb ==1  & ~isn & focus1(1,:) ==1;
    %P=polyfit(D(46,reg),D(346,reg),1);
    A=D(46,reg)';
    B=D(346,reg)';
    P=A\B; % matlab gives you the least squares solution to the A P = B
    neg = D(346,:) <= 0 & focus1(1,:) ==1;
    Dadj(1,neg)=D(346,neg)/P(1,1);
    
    reg = D(346,:) > 0 & D(46,:) > 0 & turb ==1  & ~isn & focus1(1,:) ==1;
    %P=polyfit(D(46,reg),D(346,reg),1);
    A=D(46,reg)';
    B=D(346,reg)';
    P=A\B; % matlab gives you the least squares solution to the A P = B
    pos = D(346,:) > 0 & focus1(1,:) ==1;
    Dadj(1,pos)=D(346,pos)/P(1,1);

    %----------------------------------------------------------------------
    reg = D(346,:) <= 0 & D(46,:) <= 0 & turb ==1  & ~isn & focus2(1,:) ==1;
    %P=polyfit(D(46,reg),D(346,reg),1);
    A=D(46,reg)';
    B=D(346,reg)';
    P=A\B; % matlab gives you the least squares solution to the A P = B
    neg = D(346,:) <= 0 & focus2(1,:) ==1;
    Dadj(1,neg)=D(346,neg)/P(1,1);
    
    reg = D(346,:) > 0 & D(46,:) > 0 & turb ==1  & ~isn & focus2(1,:) ==1;
    %P=polyfit(D(46,reg),D(346,reg),1);
    A=D(46,reg)';
    B=D(346,reg)';
    P=A\B; % matlab gives you the least squares solution to the A P = B
    pos = D(346,:) > 0 & focus2(1,:) ==1;
    Dadj(1,pos)=D(346,pos)/P(1,1);
    D(346,:)=Dadj;
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    %****************************************
    
    %Fill Data
    Cleaned_D(13,:)=D(46,:);%fast
    bad=isnan(Cleaned_D(13,:));
    Cleaned_D(13,bad)=D(345,bad);%dl
    bad=isnan(Cleaned_D(13,:));
    Cleaned_D(13,bad)=D(346,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(46,:));
    good_dl=(isnan(D(46,:)) & ~isnan(D(345,:)));
    good_goes=(isnan(D(46,:)) & isnan(D(345,:)) & ~isnan(D(346,:)));
    DataSource(13,good_mat)=1;
    DataSource(13,good_dl)=2;
    DataSource(13,good_goes)=3;
    
    %put source into Cleaned_D
    Cleaned_D(28,:)=DataSource(13,:);
    bad=Cleaned_D(28,:)==0;
    Cleaned_D(28,bad)=NaN;
    
    %Fill Header
    Cleaned_Header(13,1)={'CO2_FLUX'};
    Cleaned_Header(28,1)={'CO2_FLUX_SOURCE_FLAG'};
    %----------------------------------------------------------------------
    %14. FH2O
    %----------------------------------------------------------------------
    %we will use just ust filtered observations 
    turb = Cleaned_D(4,:)>0.3;
    
    %homogenize
    %----------------------------------------------------------------------
    isn=isnan(D(47,:)) ==1 | isnan(D(347,:)) == 1;
    reg = isn == 0 & turb == 1 & focus1(1,:) ==1;
    %P=polyfit(D(47,~isn),D(347,~isn),1);
    %Dadj=(D(347,:)-P(1,2))./P(1,1);
    disp('using robustfit to homogenize Fh2o')
    B = robustfit(D(47,reg),D(347,reg)); % B1(2,1) is the slope and B1(1,1) is the intercept
    Dadj=(D(347,:)-B(1,1))./B(2,1);
    D(347,focus1)=Dadj(1,focus1);
    
    isn=isnan(D(47,:)) ==1 | isnan(D(347,:)) == 1;
    reg = isn == 0 & turb == 1 & focus2(1,:) ==1; 
    %P=polyfit(D(47,~isn),D(347,~isn),1);
    %Dadj=(D(347,:)-P(1,2))./P(1,1);
    disp('using robustfit to homogenize Fh2o')
    B = robustfit(D(47,reg),D(347,reg)); % B1(2,1) is the slope and B1(1,1) is the intercept
    Dadj=(D(347,:)-B(1,1))./B(2,1);
    D(347,focus2)=Dadj(1,focus2);
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    isn=isnan(D(47,:)) ==1 | isnan(D(348,:)) == 1;
    reg = isn == 0 & turb == 1 & focus1(1,:) ==1;
    %P=polyfit(D(47,~isn),D(348,~isn),1);
    %Dadj=(D(348,:)-P(1,2))./P(1,1);
    disp('using robustfit to homogenize Fh2o')
    B = robustfit(D(47,reg),D(348,reg)); % B1(2,1) is the slope and B1(1,1) is the intercept
    Dadj=(D(348,:)-B(1,1))./B(2,1);
    D(348,focus1)=Dadj(1, focus1);
    
    isn=isnan(D(47,:)) ==1 | isnan(D(348,:)) == 1;
    reg = isn == 0 & turb == 1 & focus2(1,:) ==1;
    %P=polyfit(D(47,~isn),D(348,~isn),1);
    %Dadj=(D(348,:)-P(1,2))./P(1,1);
    disp('using robustfit to homogenize Fh2o')
    B = robustfit(D(47,reg),D(348,reg)); % B1(2,1) is the slope and B1(1,1) is the intercept
    Dadj=(D(348,:)-B(1,1))./B(2,1);
    D(348,focus2)=Dadj(1, focus2);
    %----------------------------------------------------------------------
    %****************************************
    
    %Fill Data
    Cleaned_D(14,:)=D(47,:);%fast
    bad=isnan(Cleaned_D(14,:));
    Cleaned_D(14,bad)=D(347,bad);%dl
    bad=isnan(Cleaned_D(14,:));
    Cleaned_D(14,bad)=D(348,bad);%goes
    
    %Where is the data from?
    good_mat=~isnan(D(47,:));
    good_dl=(isnan(D(47,:)) & ~isnan(D(347,:)));
    good_goes=(isnan(D(47,:)) & isnan(D(347,:)) & ~isnan(D(348,:)));
    DataSource(14,good_mat)=1;
    DataSource(14,good_dl)=2;
    DataSource(14,good_goes)=3;
    
    %put source into Cleaned_D
    Cleaned_D(29,:)=DataSource(14,:);
    bad=Cleaned_D(29,:)==0;
    Cleaned_D(29,bad)=NaN;
    
    %Fill Header
    Cleaned_Header(14,1)={'H2O_FLUX'};
    Cleaned_Header(29,1)={'H2O_FLUX_SOURCE_FLAG'};
    
    %----------------------------------------------------------------------
    %10. Latent heat
    %----------------------------------------------------------------------
    %convert Fh2o to LE
    %calculate Latent heat of vaporization using the sonic temperature
    T = D(349,:);
    Lv   = (2.501-0.00237 .* T(1,:))*10^6; %units: J/kg
    
    %calculate latent heat using combined Fh2o
    %units: [g h2o/mol h2o * J/kg h2o * mmol h2o/m^2/s * (1/1000) * (1/1000)] = [W/m2]
    Cleaned_D(10,:) = 1e-6 * 18 .* Lv .* Cleaned_D(14,:); %[W/m2]
    
    %put source into Cleaned_D - the source is from Fh2o
    Cleaned_D(27,:)=Cleaned_D(29,:);
    
    %Fill Header
    Cleaned_Header(10,1)={'LATENT_HEAT'};
    Cleaned_Header(27,1)={'LATENT_HEAT_SOURCE_FLAG'};
    
    %----------------------------------------------------------------------
    %15. Rn
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(15,:)=D(87,:);%dl
    bad=isnan(Cleaned_D(15,:));
    Cleaned_D(15,bad)=D(213,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(87,:));
    good_goes=(isnan(D(87,:)) & ~isnan(D(213,:)));
    DataSource(15,good_dl)=2;
    DataSource(15,good_goes)=3;

    %Fill Header
    Cleaned_Header(15,1)={'RNET'};
    %----------------------------------------------------------------------
    %16. PAR_In
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(16,:)=D(90,:);%dl
    bad=isnan(Cleaned_D(16,:));
    Cleaned_D(16,bad)=D(216,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(90,:));
    good_goes=(isnan(D(90,:)) & ~isnan(D(216,:)));
    DataSource(16,good_dl)=2;
    DataSource(16,good_goes)=3;

    %Fill Header
    Cleaned_Header(16,1)={'PAR_IN'};
    %----------------------------------------------------------------------
    %17. PAR_OUT
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(17,:)=D(91,:);%dl
    bad=isnan(Cleaned_D(17,:));
    Cleaned_D(17,bad)=D(217,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(91,:));
    good_goes=(isnan(D(91,:)) & ~isnan(D(217,:)));
    DataSource(17,good_dl)=2;
    DataSource(17,good_goes)=3;

    %Fill Header
    Cleaned_Header(17,1)={'PAR_OUT'};

    %----------------------------------------------------------------------
    %18. Incoming solar radiation - pyranometer = SOLAR_IN
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(18,:)=D(88,:);%dl
    bad=isnan(Cleaned_D(18,:));
    Cleaned_D(18,bad)=D(214,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(88,:));
    good_goes=(isnan(D(88,:)) & ~isnan(D(214,:)));
    DataSource(18,good_dl)=2;
    DataSource(18,good_goes)=3;

    %Fill Header
    Cleaned_Header(18,1)={'SOLAR_IN'};

    %----------------------------------------------------------------------
    %19. Outgoing solar radiation - pyranometer = SOLAR_OUT
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(19,:)=D(89,:);%dl
    bad=isnan(Cleaned_D(19,:));
    Cleaned_D(19,bad)=D(215,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(89,:));
    good_goes=(isnan(D(89,:)) & ~isnan(D(215,:)));
    DataSource(19,good_dl)=2;
    DataSource(19,good_goes)=3;

    %Fill Header
    Cleaned_Header(19,1)={'SOLAR_OUT'};
    
    %----------------------------------------------------------------------
    %20. HMP_Temp 
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(20,:)=D(97,:);%dl
    bad=isnan(Cleaned_D(20,:));
    Cleaned_D(20,bad)=D(223,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(97,:));
    good_goes=(isnan(D(97,:)) & ~isnan(D(223,:)));
    DataSource(20,good_dl)=2;
    DataSource(20,good_goes)=3;

    %Fill Header
    Cleaned_Header(20,1)={'T_HMP'}; 
    
    %----------------------------------------------------------------------
    %21. HMP_RH
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(21,:)=D(98,:);%dl
    bad=isnan(Cleaned_D(21,:));
    Cleaned_D(21,bad)=D(224,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(98,:));
    good_goes=(isnan(D(98,:)) & ~isnan(D(224,:)));
    DataSource(21,good_dl)=2;
    DataSource(21,good_goes)=3;

    %Fill Header
    Cleaned_Header(21,1)={'RH'}; 
        
    %----------------------------------------------------------------------
    %22. T_107
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(22,:)=D(102,:);%dl
    bad=isnan(Cleaned_D(22,:));
    Cleaned_D(22,bad)=D(228,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(102,:));
    good_goes=(isnan(D(102,:)) & ~isnan(D(228,:)));
    DataSource(22,good_dl)=2;
    DataSource(22,good_goes)=3;

    %Fill Header
    Cleaned_Header(22,1)={'T_107'}; 
    
    %----------------------------------------------------------------------
    %23. Rain
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(23,:)=D(103,:);%dl
    bad=isnan(Cleaned_D(23,:));
    Cleaned_D(23,bad)=D(229,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(103,:));
    good_goes=(isnan(D(103,:)) & ~isnan(D(229,:)));
    DataSource(23,good_dl)=2;
    DataSource(23,good_goes)=3;

    %Fill Header
    Cleaned_Header(23,1)={'RAIN'}; 
    
    %----------------------------------------------------------------------
    %24. Water vapor concentration
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(24,:)=D(100,:);%dl
    bad=isnan(Cleaned_D(24,:));
    Cleaned_D(24,bad)=D(226,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(100,:));
    good_goes=(isnan(D(100,:)) & ~isnan(D(226,:)));
    DataSource(24,good_dl)=2;
    DataSource(24,good_goes)=3;

    %Fill Header
    Cleaned_Header(24,1)={'H2O_HMP'}; 
    
    %----------------------------------------------------------------------
    %30. NDVI
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(30,:)=D(96,:);%dl
    bad=isnan(Cleaned_D(30,:));
    Cleaned_D(30,bad)=D(222,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(96,:));
    good_goes=(isnan(D(96,:)) & ~isnan(D(222,:)));
    DataSource(30,good_dl)=2;
    DataSource(30,good_goes)=3;

    %Fill Header
    Cleaned_Header(30,1)={'NDVI'}; 
    
    %----------------------------------------------------------------------
    %31. RED IN
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(31,:)=D(92,:);%dl
    bad=isnan(Cleaned_D(31,:));
    Cleaned_D(31,bad)=D(218,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(92,:));
    good_goes=(isnan(D(92,:)) & ~isnan(D(218,:)));
    DataSource(31,good_dl)=2;
    DataSource(31,good_goes)=3;

    %Fill Header
    Cleaned_Header(31,1)={'Red In'}; 
    
    %----------------------------------------------------------------------
    %32. RED OUT
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(32,:)=D(94,:);%dl
    bad=isnan(Cleaned_D(32,:));
    Cleaned_D(32,bad)=D(220,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(94,:));
    good_goes=(isnan(D(94,:)) & ~isnan(D(220,:)));
    DataSource(32,good_dl)=2;
    DataSource(32,good_goes)=3;

    %Fill Header
    Cleaned_Header(32,1)={'Red Out'}; 
    
    %----------------------------------------------------------------------
    %33. NIR In
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(33,:)=D(93,:);%dl
    bad=isnan(Cleaned_D(33,:));
    Cleaned_D(33,bad)=D(219,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(93,:));
    good_goes=(isnan(D(93,:)) & ~isnan(D(219,:)));
    DataSource(33,good_dl)=2;
    DataSource(33,good_goes)=3;

    %Fill Header
    Cleaned_Header(33,1)={'NIR In'}; 
    
    %----------------------------------------------------------------------
    %34. NIR Out
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(34,:)=D(95,:);%dl
    bad=isnan(Cleaned_D(34,:));
    Cleaned_D(34,bad)=D(221,bad);%goes
    
    %Where is the data from?
    good_dl=~isnan(D(95,:));
    good_goes=(isnan(D(95,:)) & ~isnan(D(221,:)));
    DataSource(34,good_dl)=2;
    DataSource(34,good_goes)=3;

    %Fill Header
    Cleaned_Header(34,1)={'NIR Out'}; 
    
    %----------------------------------------------------------------------
    %35. Soil Temperature and Moisture 
    %----------------------------------------------------------------------
    %Fill Data
    Cleaned_D(35:44,:)=D(282:291,:);%dl
    Cleaned_D(45:53,:)=D(293:301,:);%dl
    Cleaned_D(54:59,:)=D(314:319,:);%dl
    Cleaned_D(60:61,:)=D(320:321,:);%dl
    
    %Where is the data from?
    DataSource(35:61,:)=5;

    %Fill Header
    Cleaned_Header(35,1)={'Soil_M(1)'}; 
    Cleaned_Header(36,1)={'Soil_M(2)'}; 
    Cleaned_Header(37,1)={'Soil_M(3)'}; 
    Cleaned_Header(38,1)={'Soil_M(4)'}; 
    Cleaned_Header(39,1)={'Fuel_M(1)'}; 
    Cleaned_Header(40,1)={'Fuel_M(2)'}; 
    Cleaned_Header(41,1)={'Soil_T(1)'}; 
    Cleaned_Header(42,1)={'Soil_T(2)'}; 
    Cleaned_Header(43,1)={'Soil_T(3)'}; 
    Cleaned_Header(44,1)={'Soil_T(4)'}; 
    Cleaned_Header(45,1)={'LWS(1)'}; 
    Cleaned_Header(46,1)={'LWS(2)'}; 
    Cleaned_Header(47,1)={'LWS(3)'}; 
    Cleaned_Header(48,1)={'StartT_C(1)'}; 
    Cleaned_Header(49,1)={'StartT_C(2)'}; 
    Cleaned_Header(50,1)={'StartT_C(3)'}; 
    Cleaned_Header(51,1)={'StartT_C(4)'}; 
    Cleaned_Header(52,1)={'StartT_C(5)'}; 
    Cleaned_Header(53,1)={'StartT_C(6)'}; 
    Cleaned_Header(54,1)={'DelT_C(1)'}; 
    Cleaned_Header(55,1)={'DelT_C(2)'}; 
    Cleaned_Header(56,1)={'DelT_C(3)'}; 
    Cleaned_Header(57,1)={'DelT_C(4)'}; 
    Cleaned_Header(58,1)={'DelT_C(5)'}; 
    Cleaned_Header(59,1)={'DelT_C(6)'}; 
    Cleaned_Header(60,1)={'AirT'}; 
    Cleaned_Header(61,1)={'Snow Depth'}; 
    
    %Flip to maintain format
    Cleaned_D=Cleaned_D';
    DataSource=DataSource';