%function [Urot,alpha,beta]=coordrot(U,IFLAG,iok,params)
%
% 1/21/2001 - modified to not consider NaNs when calculating means
%
% Routine to rotate the 3 X N wind vector matrix U=U(u,v,w) components .
% The rotation can be either in a the horizontal (x,y) plane, or 3D
% in (x,y,z), depending upon IFLAG.  
%
% IFLAG = 1
% rotate only in the horizontal plane, such that the mean crosswind
% velocity is zero
%
% IFLAG = 0 OR NO IFLAG
% perform rotation of all axes such that both the mean crosswind and
% vertical wind are zero

function [Urot,alpha,beta]=coordrot(U,IFLAG,iok,params)

%  disp('modified coordrot 6/18/2001');

  if nargin==1
    IFLAG=0;
  end

if nargin==3 && ~isempty(find(~iok, 1))
 
  U(:,~iok) = NaN*ones(size(U,1),length(find(~iok)));

end

   iunan=find(isnan(U(1,:))); %#ok<NASGU>

if IFLAG == 0 || IFLAG == 1
   
   ubar=mean( U(1,  ~isnan( U(1,:) )  ));
   vbar=mean( U(2,  ~isnan( U(2,:) )  ));
   wbar=mean( U(3,  ~isnan( U(3,:) )  ));

   alpha=atan2(vbar,ubar);

   if IFLAG		 	% rotate 3D

     beta=0;

   else

      uhor=sqrt(ubar^2+vbar^2);
      beta=atan2(wbar,uhor);
 
   end;

   Urot(1,:)= U(1,:)*cos(alpha)*cos(beta) + U(2,:)*sin(alpha)*cos(beta) + U(3,:)*sin(beta);
   Urot(2,:)=-U(1,:)*sin(alpha)           + U(2,:)*cos(alpha)                        ;
   Urot(3,:)=-U(1,:)*cos(alpha)*sin(beta) - U(2,:)*sin(alpha)*sin(beta) + U(3,:)*cos(beta);


elseif IFLAG==2
  
  disp('Planar Fit');
 
  c0=params(1);
  alpha=params(2);
  beta=params(3);
  
  % set b0, b1, b2
    
%  c0 = params(1);
%  b1 = params(2);
%  b2 = params(3);
  
%  p31=-b1/sqrt(b1^2+b2^2+1);
%  p32=-b2/sqrt(b1^2+b2^2+1);
%  p33=1/sqrt(b1^2+b2^2+1);

%  sinb=-p32/sqrt(p32^2+p33^2);
%  cosb=p33/sqrt(p32^2+p33^2);
%  sina=p31;
%  cosa=sqrt(p32^2+p33^2);

  DMAT = [  cos(alpha) 0 sin(alpha) ;
            0    1   0  ;
           -sin(alpha) 0 cos(alpha)];

  CMAT = [  1    0     0 ;
            0  cos(beta) -sin(beta);
            0  sin(beta)  cos(beta)];

   P=DMAT'*CMAT';

   Urot = P*(U - [0 0 c0]'*ones(1,size(U,2)));
  
   vbar = mean(Urot(2,:));
   ubar = mean(Urot(1,:));
   
   gamma=atan2(vbar,ubar);
   
   M=[ cos(gamma) sin(gamma) 0;
      -sin(gamma) cos(gamma) 0;
      0              0       1];

   Urot = M*Urot;
   
   % return gamma as alpha, the second output, so i can rotate back
   % in the horizontal plane to save components...

   alpha=gamma;
   
end


return
