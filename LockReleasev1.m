clc;clear all; close all; format short

% pwd = '/home/valera/ParGCCOM/TestCases/ShearLayer/SL_128x6x128/';
% ProbName = 'shearlayer';
% ProbSize = '128x6x128';
% 
% GridFileName =  [pwd ProbName '.' ProbSize '.dat' ];  ;
% 
% ProbSizeFileName = [pwd ProbName '.' ProbSize 'probsize.ser.dat' ];
% 
% INIfile = 'gcom_ini_DC_128x128x6.nc';
% 
% [x,y,z,x_r,y_r,z_r,IMax,JMax,KMax]=ReadGridASCII(ProbSizeFileName,GridFileName);

%Nodes:
IMax = 101; JMax = 101;
%Spacing:
h = 1/IMax;

%Cells:
m = IMax-1; n = JMax-1;

xlength = linspace(-16,16,n+1);
ylength = linspace(-1,1,m+1);

[X,Y] = ndgrid(xlength,ylength);
%%

uic = zeros(IMax ,JMax);
vic = zeros(IMax ,JMax);
pic = zeros(IMax ,JMax);
Tic = zeros(IMax, JMax);
Sic = zeros(IMax, JMax);


p = zeros(IMax+1,JMax+1);
u = zeros(IMax+1,JMax+1);
v = zeros(IMax+1,JMax+1);
T = zeros(IMax+1,JMax+1);
S = zeros(IMax+1,JMax+1);

xlength = linspace(-16,16,n);
ylength = linspace(-1,1,m);

[Xc,Yc] = ndgrid(xlength,ylength);

xlength = linspace(-16,16,n+1);
ylength = linspace(-1,1,m);

[Xu,Yu] = ndgrid(xlength,ylength);

xlength = linspace(-16,16,n);
ylength = linspace(-1,1,m+1);

[Xv,Yv] = ndgrid(xlength,ylength);

xlength = linspace(-16,16,n+2);
ylength = linspace(-1,1,m+2);

[Xp,Yp] = ndgrid([1 1.5 : 1 : JMax-0.5 JMax], [1 1.5 : 1 : IMax-0.5 IMax]);

Xp = Xp - JMax/2.0;
XMax = max(max(Xp));
Xp = Xp/XMax*16.0; 

Yp = Yp - IMax/2.0;
YMax = max(max(Yp));
Yp = Yp/YMax;               

Dx = X(1,2)-X(1,1);
Dy = Y(2,1)-Y(1,1);
dw = 80.0;
dp = 0.05;
%%

TMAX = 16.1324636; TMIN = 10.0;

T0 = TMAX*100;

for j=1:JMax
    for i=1:IMax

        Tic(i,j) = TMIN + (TMAX-TMIN)*(0.5-0.5*erf(Y(i,j)/0.01));
        
    
        end
end

figure(1);
u_v = squeeze(Tic);
pcolor(u_v);
title('T')
colorbar;
shading interp;

Tic = Tic/T0;

T = Tic;

%%
rho_ref = 1027.d0;
S_ref   = 35.d0;
T_ref   = 10.d0;
drho_dS = 0.781d0;
drho_dT = -0.1708d0;
LStar = 0.05;
UStar = 1.0;
gforce = 9.80;

for j=1:JMax
    for i=1:IMax

        dens(i,j) = rho_ref + drho_dT*(T(i,j)*T0-T_ref);  %EOS  
        
    end
end


DMax = max(max(dens));
DMin = min(min(dens));

rho_P0 = (DMax+DMin)/2.0;

figure(6);
u_v = squeeze(dens);
pcolor(u_v);
title('rho')
colorbar;
shading interp;

%%
buoy = zeros(IMax,JMax-1);

buoy(1,:) = gforce*LStar*(dens(1,1:n)-rho_P0)/(rho_P0*(UStar^2));

for j=1:JMax-1
    for i=1:IMax-1

        buoy(i,j) = gforce*LStar*(0.5*(dens(i,j)+dens(i+1,j))-rho_P0)/(rho_P0*(UStar^2));

    end
end

buoy(end,:) = buoy(end-1,:);

figure(6);
u_v = squeeze(buoy);
pcolor(u_v);
title('(rho-rho_0)/rho_0')
colorbar;
shading interp;

bic = buoy;


%%
% 
% 
% for j=1:JMax
%     for i=1:IMax
% 
%         uic(i,j) = X(i,j)/16;
%         vic(i,j) = Y(i,j);
%         
%     
%         end
% end




%%  Creating operator
k = 2; % Mimetic order of accuracy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dy = Xc(2,1)-Xc(1,1);
Dx = Dy;
D = div2D(k,m,Dx,n,Dy); %div2DCurv(k,X,Y);
G =  grad2D(k,m,Dx,n,Dy); % grad2DCurv(k,X,Y);

%L = D*G;  %lap2D(k,IMax-1,Dx,JMax-1,Dy);
L = lap2D(k,IMax-1,Dx,JMax-1,Dy);
L = L + robinBC2D(k,m,Dx,n,Dy,0,1);

idx1 = numel(uic);
idy1 = numel(vic);
idx2 = idx1 + idy1;


dt = 0.1;
rho = rho_P0;
D = rho/dt*D;
G = -dt/rho*G;
%% Preparing t=0 arrays

uic = zeros(IMax ,JMax);
vic = zeros(IMax ,JMax);


iatu = (X(1:IMax-1,:)+X(2:IMax,:))*0.50;
jatu = (Y(1:IMax-1,:)+Y(2:IMax,:))*0.50;

iatv = (X(:,1:JMax-1)+X(:,2:JMax))*0.50;
jatv = (Y(:,1:JMax-1)+Y(:,2:JMax))*0.50;

Fuic = griddedInterpolant(X, Y, uic); 
Fvic = griddedInterpolant(X, Y, vic); 
Fpic = griddedInterpolant(X, Y, pic); 
%FTic = griddedInterpolant(X, Y, Tatc); 


uatu = Fuic(iatu,jatu);
vatu = Fvic(iatu,jatu);
patu = Fpic(iatu,jatu);

uatv = Fuic(iatv,jatv);
vatv = Fvic(iatv,jatv);
patv = Fpic(iatv,jatv);

u = Fuic(Xp,Yp);
v = Fvic(Xp,Yp);
p = Fpic(Xp,Yp);
%T = Fpic(Xp,Yp);

T = Tic;
buoy = bic;

% 
   figure(1);
    u_v = squeeze(uatu(:,:));
    subplot(2,2,1)
    pcolor(u_v);
    title('u')
    colorbar;
    shading interp;
    subplot(2,2,2)
    u_2 = squeeze(p(:,:));
    pcolor(u_2);
    title('p')
    colorbar;
    shading interp;
    subplot(2,2,3)
    u_3 = squeeze(vatv(:,:));
    pcolor(u_3);
    title('v')
    colorbar;
    shading interp;
    subplot(2,2,4)
    u_3 = squeeze(T(:,:));
    pcolor(u_3);
    title('T')
    colorbar;
    shading interp;


% Reshaping


Vv = reshape(vatv',[], 1);    %Do we need to transpose now with ndgrid?
Uv = reshape(uatv',[], 1);
Pv = reshape(patv',[], 1);

Vu = reshape(vatu',[], 1);
Uu = reshape(uatu',[], 1);
Pu = reshape(patu',[], 1);

Vp = reshape(v',[], 1);
Up = reshape(u',[], 1);
Pp = reshape(p',[], 1);

cs = 1447.0;


tf = 1.7;
nf = tf/dt ;
beta = 1/cs;
Re = 750;
IRe = 1.0/Re;
nu = 0.0014;
time_steps = linspace(0,tf,nf);

[uatu,vatv,p] = applyboundaries2D(uatu,vatv,p); 
%uatu(end,:) = 1.0;

Fuatu = griddedInterpolant(iatu, jatu, uatu); 
Fvatv = griddedInterpolant(iatv, jatv, vatv); 
u = Fuatu(Xp,Yp);
v = Fvatv(Xp,Yp);

Fu = griddedInterpolant(Xp,Yp, u);
Fv = griddedInterpolant(Xp,Yp, v);

uatv = Fu(iatv,jatv);
vatu = Fv(iatu,jatu);

    figure(2);
    u_v = squeeze(uatv(:,:));
    subplot(2,2,1)
    pcolor(u_v);
    title('u')
    colorbar;
    shading interp;
    subplot(2,2,2)
    u_2 = squeeze(p(:,:));
    pcolor(u_2);
    title('p')
    colorbar;
    shading interp;
    subplot(2,2,3)
    u_3 = squeeze(vatu(:,:));
    pcolor(u_3);
    title('v')
    colorbar;
    shading interp;
%     subplot(2,2,4)
%     pcolor(Xp,Yp,cav)
%     hold on
%     quiver(Xp,Yp,u,v)
%     hold off
%     xlim([0 1])
%     ylim([0 1])
%     colorbar;
%     shading interp;
%     title(['vorticity Dt = ' num2str((time+1)*dt)])

%Dy = Xp(1,3)-Xp(1,2);
%Dx = Dy;

%Dx = 1.0; 
%Dy = 1.0; 



%% Boussinesq Pred-Corr formulation:
for time = 1:size(time_steps,2)

    %Predict u*

       %RHS = zeros(IMax+1,JMax+1);
    %u_star = zeros(IMax+1,JMax+1); %u;
    %v_star = zeros(IMax+1,JMax+1); %v;
    
    u_star = uatu;
    v_star = vatv;
    
   
    for j = 2:JMax-1
        for i = 2:IMax-2
    
            d2u_dy2 = (uatu(i-1,j)-2*uatu(i,j)+uatu(i+1,j))/(Dy^2);
            d2u_dx2 = (uatu(i,j-1)-2*uatu(i,j)+uatu(i,j+1))/(Dx^2);
                        
            du_dx = (uatu(i,j+1)-uatu(i,j-1))/(2*Dx);
            du_dy = (uatu(i+1,j)-uatu(i-1,j))/(2*Dy);
            udu_dx  = uatu(i,j)*du_dx;
            vdu_dy  = vatu(i,j)*du_dy;  

                      
            u_star(i,j) = uatu(i,j) - dt*(udu_dx+vdu_dy - IRe*(d2u_dy2+d2u_dx2));
        end
    end
    
    %Predict v*
     
    for j = 2:JMax-2
        for i = 2:IMax-1
    
            d2v_dy2 = (vatv(i-1,j)-2*vatv(i,j)+vatv(i+1,j))/(Dy^2);
            d2v_dx2 = (vatv(i,j-1)-2*vatv(i,j)+vatv(i,j+1))/(Dx^2);
                        
            dv_dx = (vatv(i,j+1)-vatv(i,j-1))/(2*Dx);
            dv_dy = (vatv(i+1,j)-vatv(i-1,j))/(2*Dy);       
            udv_dx  = uatv(i,j)*dv_dx;
            vdv_dy  = vatv(i,j)*dv_dy;
                              
            v_star(i,j) = vatv(i,j) - dt*( udv_dx+vdv_dy - IRe*(d2v_dy2+d2v_dx2));% - buoy(i,j) );
           
        end
    end

    v_star = v_star + dt*buoy;
    
    [u_star,v_star,~] = applyboundaries2Dstar(u_star,v_star,p);

    figure(4);
    u_v = squeeze(u_star(:,:));
    subplot(2,2,1)
    pcolor(u_v);
    title('u*')
    colorbar;
    shading interp;
    subplot(2,2,2)
    u_2 = squeeze(v_star(:,:));
    pcolor(u_2);
    title('v*')
    colorbar;
    shading interp;
    subplot(2,2,3)
    u_3 = squeeze(p(:,:));
    pcolor(u_3);
    title('p')
    colorbar;
    shading interp;
    subplot(2,2,4)
    u_3 = squeeze(T(:,:));
    pcolor(u_3);
    title('T')
    colorbar;
    shading interp;
    
    Uu = reshape(u_star',[], 1);
    Vv = reshape(v_star',[], 1);
   
    %Get p (continuity)
    
    R = [Uu ; Vv];  
    Pp = L\(D*R);

    Gp = G*[ Pp ]; 
    dpdx = Gp(1:(IMax*(JMax-1)));
    dpdy = Gp((IMax*(JMax-1)+1):end);
    
    %Correct u,v:
    
    Uu = Uu - dpdx;
    Vv = Vv - dpdy;
    
    uatu = reshape(Uu,[IMax JMax-1])' ; 
    vatv = reshape(Vv,[IMax-1 JMax])' ; 
    
    [uatu,vatv,p] = applyboundaries2D(uatu,vatv,p);   
    
    Fuatu = griddedInterpolant(iatu, jatu, uatu); 
    Fvatv = griddedInterpolant(iatv, jatv, vatv); 
    u = Fuatu(Xp,Yp);
    v = Fvatv(Xp,Yp);

    p = reshape(Pp,[IMax+1 JMax+1])';
    
    Fu = griddedInterpolant(Xp,Yp, u);
    Fv = griddedInterpolant(Xp,Yp, v);

    uatv = Fu(iatv,jatv);
    vatu = Fv(iatu,jatu);    
     
    uatc = Fu(X,Y);
    vatc = Fv(X,Y); 
    
    [~,cav] = curl(Xp,Yp,u,v);
              
    %Interpolate and refresh arrays:
  
    Pv = reshape(patv',[], 1);
    Pu = reshape(patu',[], 1);
    
    Vv = reshape(vatv',[], 1);
    Uv = reshape(uatv',[], 1);
 
    Vu = reshape(vatu',[], 1);
    Uu = reshape(uatu',[], 1);    
    
    for j = 2:JMax-2
      for i = 2:IMax-2
    
            d2T_dy2 = (T(i-1,j)-2*T(i,j)+T(i+1,j))/(Dy^2);
            d2T_dx2 = (T(i,j-1)-2*T(i,j)+T(i,j+1))/(Dx^2);
                        
            dT_dx = (T(i,j+1)-T(i,j-1))/(2*Dx);
            dT_dy = (T(i+1,j)-T(i-1,j))/(2*Dy);       
            udT_dx  = uatc(i,j)*dT_dx;
            vdT_dy  = vatc(i,j)*dT_dy;
                    
           
            T(i,j) = T(i,j) - dt*(udT_dx+vdT_dy - IRe*(d2T_dy2+d2T_dx2));
           
        end
    end
    
%     T(1,:) = 0.0 ; %T(2,:);
%     T(end,:) = 0.0 ; % T(end-1,:);    
% 
%     
%     T(:,1) = 0.0 ; % T(:,2);
%     T(:,end) = 0.0 ; %T(:,end-1);  
%     
     
    dens = rho_ref + drho_dT*(T*T0-T_ref);  %EOS  
    
    buoy(1,:) = gforce*LStar*(dens(1,1:n)-rho_P0)/(rho_P0*(UStar^2));

    for j=1:JMax-1
        for i=1:IMax-1

            buoy(i,j) = gforce*LStar*(0.5*(dens(i,j)+dens(i+1,j))-rho_P0)/(rho_P0*(UStar^2));

        end
    end
    
    buoy(end,:) = buoy(end-1,:);
    
    figure(2);
    u_v = squeeze(u(:,:));
    subplot(3,2,1)
    pcolor(u_v);
    title('u')
    colorbar;
    shading interp;
    subplot(3,2,2)
    u_3 = squeeze(v(:,:));
    pcolor(u_3);
    title('v')
    colorbar;
    shading interp;
    subplot(3,2,3)
    u_2 = squeeze(p(:,:));
    pcolor(u_2);
    title('p')
    colorbar;
    shading interp;
    subplot(3,2,4)
    u_3 = squeeze(T(:,:));
    pcolor(u_3);
    title(['T ' num2str((time+1)*dt)])
    colorbar;
    shading interp;
    subplot(3,2,5)
    u_3 = squeeze(buoy(:,:));
    pcolor(u_3);
    title(['buoy'])
    colorbar;
    shading interp;
    
    
    
%     pcolor(Xp,Yp,cav)
%     hold on
%     quiver(Xp,Yp,u,v)
%     hold off
%     xlim([0 1])
%     ylim([0 1])
%     colorbar;
%     shading interp;
%     title(['vorticity Dt = ' num2str((time+1)*dt)])


end


