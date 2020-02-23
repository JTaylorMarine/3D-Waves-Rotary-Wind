                       %% 3-D UNIFIED MODEL
                       %% 3-D UNIFIED MODEL
clear all; 
                   %% INPUT PARAMETERS 
                   %  INPUT PARAMETERS
                   
                   % MODEL DOMAIN
Lx=1000000;                % X-Size of the model domain (m);
Ly=1000000;                % Y-Size of the model domain (m);   
H0=100;                    % Depth (m)
dx=20000;                  % Spatial x grid step (m);
dy=20000;                  % Spatial y grid step (m)

                    % TOPOGRAPHY 1: SEAMOUNT
Hg=80;                     % Height of the bank  (m)
% XWidth=1500000;            % X-Width of the bank (m)
% YWidth=1500000;            % Y-Width of the bank (m)
% Xcentre=0;                 % X-positioin of the bank centre
% Ycentre=1000000;           % Y-positioin of the bank centre                       
XWidth=100000;             % X-Width of the bank (m)
YWidth=100000;             % Y-Width of the bank (m)
Xcentre=800000;            % X-positioin of the bank centre
Ycentre=500000;            % Y-positioin of the bank centre                       
                    % TOPOGRAPHY 3: SLOPE-SHELF
LB1=500000;                % LEFT  X-coordinate of the shelf break (m)(Must be between 0 and Lx)
LB2=900000;                % RIGHT X-coordinate of the shelf break (m)(Must be between 0 and Lx)
Alfa=0.3;                    % Angle of bottom inclination (deg)
BotSmooth=5;               % Smoothing parameter of the bottom (1,2,3,...)
                   % INCOMING WAVE
AmpZ=0;                    % Amplitude of the incident wave (m)
T=2;                       % Wave Period (hours)
                   % INITIAL FREE SURFACE ELEVATION (TSUNAMI)
AmpQuake=0;                % Initial wave amplitude in the epicentre (m)
XQuake=500000;             % X position of the epicentre (m)
YQuake=500000;             % Y position of the epicentre (m)
WidthQuake=100000;         % Half-Width of the earthquake (m)
                   % INITIAL FREE SURFACE ELEVATION (SEICHE)
AmpSeiche=1;               % Initial wave amplitude in the epicentre (m)
                   % NUMERICAL INPUT PARAMETERS
Lat=0;                     % Latitude (deg)
N=1000;                    % Number of temporal steps
nc=1;                      % Parameterof nonlinearity: must be either 0 or 1
filt=0.2;                  % Parameter for filtering high frequency (from 0 to 1)
nn=20;                     % Parameter for the movie to show every nn-th frame            

                  %% CALCULATED MODEL PREPARATION 
x=[0:dx:Lx];                % X-grid
y=[0:dy:Ly];                % Y-grid
J=length(x);                % Number of the OX grid steps (length of the domain) 
I=length(y);                % Number of the OY grid steps (width of the domain)
Omega=0.0000729;            % Angular speed of the Earth rotaion (rad/sec)              
f=2*Omega*sind(Lat);        % Coriolis parameter
g=9.81;                     % Acceleration due to the gravity (m/sec*sec)
dt=dx/sqrt(2*g*H0)/2;       % temporal step

                   %%Wind 

ra=1.2;              %Density of air (kg/m^3)
rwater=1028;         %Desnity fo water (kg/m^3)
Ca=1.3e-3;           %Drag coefficient
Um=10;               %Easterly wind (m/sec)
Vm=0;                %Notherly wind (m/sec)

for j=1:J
    for i=1:I
        Uwind(j,i)=Um*cos(pi*(i-1)*dy/Ly);
        Vwind(j,i)=0;
        tUwind(j,i)=ra*Ca*Uwind(j,i)*abs(Uwind(j,i));
        tVwind(j,i)=ra*Ca*Vwind(j,i)*abs(Vwind(j,i));
    end
end
            %% PARAMETERS OF THE INCIDENT WAVE

sig=2*pi/T/3600;                % Wave frequency (1/sec)
k=sqrt((sig^2-f^2)/g/H0);       % Wave number (1/m)
Lambda=2*pi/k;                  % Wavelength (m)
AmpU=AmpZ*g*k*sig/(sig^2-f^2);  % Amplitude of U velocity (m/sec)
AmpV=AmpZ*g*k*f/(sig^2-f^2);    % Amplitude of U velocity (m/sec)

                       %% TOPOGRAPHY 1: SEAMOUNT
                      
for j=1:J 
      for i=1:I
          H(j,i)=H0;
      end
end
Jcentre=Xcentre/dx;
Jwidth=XWidth/dx;
Icentre=Ycentre/dy;
Iwidth=YWidth/dy;
Jbank1=Jcentre-Jwidth;
Jbank2=Jcentre+Jwidth;
Ibank1=Icentre-Iwidth;
Ibank2=Icentre+Iwidth;
if Ibank1<=0
    Ibank11=1;
else
    Ibank11=Ibank1;
end

if Ibank2>=I
    Ibank22=I;
else
    Ibank22=Ibank2;
end
%%%%%%%%%%%%%%%%5
 if Jbank1<=0
    Jbank11=1;
else
    Jbank11=Jbank1;
end

if Jbank2>=J
    Jbank22=J;
else
    Jbank22=Jbank2;
end   

for j=Jbank11:Jbank22
      for i=Ibank11:Ibank22
           H(j,i)=H0-Hg*cos(pi*(j-Jcentre)/(Jbank2-Jbank1))^2*...
                  cos(pi*(i-Icentre)/(Ibank2-Ibank1))^2;           % TOPO 1
%          H(j,i)=H0-Hg;                                           % TOPO 2
         
     end
end
%                        %% TOPOGRAPHY 3: SHELF-SLOPE
%                        % SHELF-SLOPE                             % TOPO 3
                                                                   
for j=1:J 
      for i=1:I
          H(j,i)=H0;
      end
end
 

JB1=LB1/dx+1;
JB2=LB2/dx+1;

for i=1:I 
    JB=JB1+floor((i-1)/(I-1)*(JB2-JB1));
      for j=JB:-1:1
          H(j,i)=H0-(JB-j)*dx*tand(Alfa);
          if H(j,i)<H0-Hg
              H(j,i)=H0-Hg;
          else
              H(j,i)=H(j,i);
          end
      end
end

           HH=smooth2a(H,BotSmooth,BotSmooth);
           H(:,:)=HH(:,:);         
            


               %% INITIAL CONDITIONS 1 (Tsunami)             
        z(1,1:J,1:I)=0;        % Surface displacement (m)
        u(1,1:J,1:I)=0;        % X velocity (m/sec)
        v(1,1:J,1:I)=0;        % Y velocity (m/sec)
        
      %      Initial conditions for the tsunami wave

if AmpQuake>0
       for i=1:I
       for j=1:J 
          Dist=sqrt(((j-1)*dx-XQuake)^2+((i-1)*dy-YQuake)^2);  
         if Dist/WidthQuake <= 1
             z(1,j,i)=AmpQuake*(cos(0.5*pi*Dist/WidthQuake))^2;
         end
       end
       end
else 
                    %      INITIAL CONDITION 3 (SEICHE)
       for i=1:I
       for j=1:J 
             z(1,j,i)=AmpSeiche*(j-1)/(J-1);
         end
       end
end
                          %%  MAIN BLOCK: FTCS scheme for n=1
    n=1; 
     for j=2:J-1 
      for i=2:I-1
        u(n+1,j,i)=u(n,j,i)+f*dt*v(n,j,i)-g*dt*(z(n,j+1,i)-z(n,j-1,i))/2/dx;
        v(n+1,j,i)=v(n,j,i)-f*dt*u(n,j,i)-g*dt*(z(n,j,i+1)-z(n,j,i-1))/2/dy;
        z(n+1,j,i)=z(n,j,i)-dt*((H(j+1,i)+z(n,j+1,i))*u(n,j+1,i)-...
                   (H(j-1,i)+z(n,j-1,i))*u(n,j-1,i))/2/dx-...
                    dt*((H(j+1,i)+z(n,j,i+1))*v(n,j,i+1)-...
                    (H(j,i-1)+z(n,j,i-1))*v(n,j,i-1))/2/dy;
      end   
    end
                           % BOUNDARY CONDITIONS  
                          
                           % Open boundary (North)

        u(n+1,J,:)=0;                         % Wall
        z(n+1,J,:)= z(n+1,J-1,:);             % Wall
        v(n+1,J,:)= v(n+1,J-1,:);             % Wall
        
        if AmpZ>0 
             u(n+1,J,:)=-AmpU*sin(sig*dt*n);  % Wave
             z(n+1,J,:)=AmpZ*sin(sig*dt*n);   % Wave
             v(n+1,J,:)=-AmpV*cos(sig*dt*n);  % Wave
        end
                 
         if  sig*dt*n>pi
              u(n+1,J,:)=0;                   % Wall
              z(n+1,J,:)= z(n+1,J-1,:);       % Wall
              v(n+1,J,:)= v(n+1,J-1,:);       % Wall   
        end
                           
                           % Left i=1 and right i=I boundaries    
        u(n+1,:,1)=u(n+1,:,2);
        z(n+1,:,1)=z(n+1,:,2);
        v(n+1,:,1)=0;
        u(n+1,:,I)=u(n+1,:,I-1);
        z(n+1,:,I)=z(n+1,:,I-1);
        v(n+1,:,I)=0;
    
                             % Bottom j=1 boundary
                                 %  Shoreline (south)
        u(n+1,1,:)=0;             % Wall
        z(n+1,1,:)=z(n+1,2,:);    % dz/dx=0 
        v(n+1,1,:)=v(n+1,2,:);    % dv/dx=0       

                        % Corners
         u(n+1,1,1)=0;                    % Left-bottom 
         u(n+1,1,I)=0;                    % Right-bottom 
         u(n+1,J,I)=(u(n+1,J-1,I)+u(n+1,J,I-1))/2; % Right-Top
         u(n+1,J,1)=(u(n+1,J-1,1)+u(n+1,J,2))/2;   % Left-top    
       
         v(n+1,1,1)=0;                    % Left-bottom 
         v(n+1,1,I)=0;                    % Right-bottom 
         v(n+1,J,I)=0;                    % Right-Top
         v(n+1,J,1)=0;                    % Left-top
        
         z(n+1,1,1)=(z(n+1,1,2)+z(n+1,2,1))/2;     % Left-bottom 
         z(n+1,1,I)=(z(n+1,1,I-1)+z(n+1,2,I))/2;   % Right-bottom
         z(n+1,J,I)=(z(n+1,J-1,I)+z(n+1,J,I-1))/2; % Right-Top
         z(n+1,J,1)=(z(n+1,J-1,1)+z(n+1,J,2))/2;   % Left-top    


                                %%  MAIN BLOCK: CTCS scheme for n>1
gdtdx=g*dt/dx;                                
gdtdy=g*dt/dy; 
fdt=2*f*dt;
ncdtdx=nc*dt/dx;
ncdtdy=nc*dt/dy;
dtdx=dt/dx;
dtdy=dt/dy;

  for n=2:N-1                  
           for j=2:J-1  
               for i=2:I-1
      u(n+1,j,i)=u(n-1,j,i)+fdt*v(n,j,i)-gdtdx*(z(n,j+1,i)-z(n,j-1,i))...    
     -ncdtdx*u(n,j,i)*(u(n,j+1,i)-u(n,j-1,i))-ncdtdy*v(n,j,i)*(u(n,j,i+1)-u(n,j,i-1))...
             +2*dt*tUwind(j,i)/rwater/(H(j,i)+nc*z(n,j,i));
      v(n+1,j,i)=v(n-1,j,i)-fdt*u(n,j,i)-gdtdy*(z(n,j,i+1)-z(n,j,i-1))...
     -ncdtdx*u(n,j,i)*(v(n,j+1,i)-v(n,j-1,i))-ncdtdy*v(n,j,i)*(v(n,j,i+1)-v(n,j,i-1))... 
         +2*dt*tVwind(j,i)/rwater/(H(j,i)+nc*z(n,j,i));        
      z(n+1,j,i)=z(n-1,j,i)-dtdx*((H(j+1,i)+nc*z(n,j+1,i))*u(n,j+1,i)-(H(j-1,i)+nc*z(n,j-1,i))*u(n,j-1,i))...
     -dtdy*((H(j,i+1)+nc*z(n,j,i+1))*v(n,j,i+1)-(H(j,i-1)+nc*z(n,j,i-1))*v(n,j,i-1));               
               end   
           end
                           % BOUNDARY CONDITIONS  
                          
                           % Open boundary (North)

        u(n+1,J,:)=0;                         % Wall
        z(n+1,J,:)= z(n+1,J-1,:);             % Wall
        v(n+1,J,:)= v(n+1,J-1,:);             % Wall
        
        if AmpZ>0 
             u(n+1,J,:)=-AmpU*sin(sig*dt*n);  % Wave
             z(n+1,J,:)=AmpZ*sin(sig*dt*n);   % Wave
             v(n+1,J,:)=-AmpV*cos(sig*dt*n);  % Wave
        end
                 
         if  sig*dt*n>pi
              u(n+1,J,:)=0;                   % Wall
              z(n+1,J,:)= z(n+1,J-1,:);       % Wall
              v(n+1,J,:)= v(n+1,J-1,:);       % Wall   
        end
                           
                           % Left i=1 and right i=I boundaries    
        u(n+1,:,1)=u(n+1,:,2);
        z(n+1,:,1)=z(n+1,:,2);
        v(n+1,:,1)=0;
        u(n+1,:,I)=u(n+1,:,I-1);
        z(n+1,:,I)=z(n+1,:,I-1);
        v(n+1,:,I)=0;
    
                             % Bottom j=1 boundary
                                 %  Shoreline (south)
        u(n+1,1,:)=0;             % Wall
        z(n+1,1,:)=z(n+1,2,:);    % dz/dx=0 
        v(n+1,1,:)=v(n+1,2,:);    % dv/dx=0       

                        % Corners
         u(n+1,1,1)=0;                    % Left-bottom 
         u(n+1,1,I)=0;                    % Right-bottom 
         u(n+1,J,I)=(u(n+1,J-1,I)+u(n+1,J,I-1))/2; % Right-Top
         u(n+1,J,1)=(u(n+1,J-1,1)+u(n+1,J,2))/2;   % Left-top    
       
         v(n+1,1,1)=0;                    % Left-bottom 
         v(n+1,1,I)=0;                    % Right-bottom 
         v(n+1,J,I)=0;                    % Right-Top
         v(n+1,J,1)=0;                    % Left-top
        
         z(n+1,1,1)=(z(n+1,1,2)+z(n+1,2,1))/2;     % Left-bottom 
         z(n+1,1,I)=(z(n+1,1,I-1)+z(n+1,2,I))/2;   % Right-bottom
         z(n+1,J,I)=(z(n+1,J-1,I)+z(n+1,J,I-1))/2; % Right-Top
         z(n+1,J,1)=(z(n+1,J-1,1)+z(n+1,J,2))/2;   % Left-top    
 
%          filt=0;
    for j=2:J-1 
      for i=2:I-1
        zn(j,i)=z(n+1,j,i)*(1-filt)+0.25*filt*(z(n+1,j-1,i)+z(n+1,j+1,i)+z(n+1,j,i-1)+z(n+1,j,i+1));
      end   
    end
    
   for j=2:J-1 
      for i=2:I-1
        z(n+1,j,i)=zn(j,i);
      end   
    end
         
end        % End of the temporal loop
      
                    %%  VISUALISATION - WAVE
if AmpQuake>AmpSeiche
   Comp=AmpQuake;
else 
    Comp=AmpSeiche;
end

Scale=0;
if AmpZ>Comp
    Scale=AmpZ;
    koef=2.5;
else
     Scale=Comp;
     koef=1.5;
end


k=H0/5/Scale;                    
[Y,X]=meshgrid(y/1000,x/1000);        
fig1=figure(1);
 set(gcf,'position',[100 70 1700 880])
%nn=10;
numframes=N;
AA=moviein(numframes,fig1);
set(fig1,'NextPlot','replacechildren')
for i=1:nn:numframes
T=dt*i/3600;  
zfin(:,:)=z(i,:,:);   

surf(X,Y,zfin,'linestyle','none','FaceColor','interp',	'EdgeColor','none','FaceLighting','phong');
caxis([-Scale Scale])
alpha(0.85)
hold on
surf(X,Y,-H/k,'linestyle','none')
box on
hold off
zlim([0*Scale Scale*koef]);

view(-50,50)
camlight left
% Create ylabel
zlabel('Displacement (m)');
xlabel('X-Distance (km)','rotation',25)
ylabel('Y-Distance (km)','rotation',-15)

AA(:,i)=getframe(fig1);
title(['Time:' num2str(T) 'h']);
pause(0.2)
end

    %% Visualization of Wind Driven Circulation
ufin(:,:)=u(N-1,:,:);
vfin(:,:)=v(N-1,:,:);

[Y X]=meshgrid(y/1000,x/1000);

nq=3;
figure;
quiver(X(2:nq:J,2:nq:I),Y(2:nq:J,2:nq:I),ufin(2:nq:J,2:nq:I),...
vfin(2:nq:J,2:nq:I),'AutoScaleFactor',2,...
'LineWidth',1.5,'MaxHeadSize',5.5);
xlim([0 Lx/1000]);
ylim([0 Ly/1000]);
%Create xLabel
xlabel({'Distance (km)',''});
%Create yLabel
ylabel({'Distance (km)',''});
