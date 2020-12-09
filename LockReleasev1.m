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

[X,Y] = meshgrid(xlength,ylength);
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

[Xc,Yc] = meshgrid(xlength,ylength);

xlength = linspace(-16,16,n+1);
ylength = linspace(-1,1,m);

[Xu,Yu] = meshgrid(xlength,ylength);

xlength = linspace(-16,16,n);
ylength = linspace(-1,1,m+1);

[Xv,Yv] = meshgrid(xlength,ylength);

xlength = linspace(-16,16,n+2);
ylength = linspace(-1,1,m+2);

[Xp,Yp] = meshgrid([1 1.5 : 1 : JMax-0.5 JMax], [1 1.5 : 1 : IMax-0.5 IMax]);

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

        Tic(i,j) = TMIN + (TMAX-TMIN)*(0.5-0.5*erf(X(i,j)/0.01));
        
    
        end
end

figure(1);
u_v = squeeze(Tic);
pcolor(u_v);
title('T')
colorbar;
shading interp;

Tatc = Tic/T0;

%%
rho_ref = 1027.d0;
S_ref   = 35.d0;
T_ref   = 10.d0;
drho_dS = 0.781d0;
drho_dT = -0.1708d0;
LStar = 0.05;
UStar = 1.0;
gforce = 9.80;

dens = rho_ref + drho_dT*(Tatc*T0-T_ref);  %EOS  

DMax = max(max(dens));
DMin = min(min(dens));

rho_P0 = (DMax+DMin)/2.0;

buoy = dens;

buoy(1,:) = gforce*LStar*(dens(1,:)-rho_P0)/(rho_P0*(UStar^2));

for j=2:JMax
    for i=1:IMax-1

        buoy(i,j) = gforce*LStar*(0.5*(dens(i,j)+dens(i+1,j))-rho_P0)/(rho_P0*(UStar^2));

    end
end



%%  Creating operator
k = 2; % Mimetic order of accuracy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dy = Xc(1,2)-Xc(1,1);
Dx = Dy;
D = div2D(k,m,Dx,n,Dy); %div2DCurv(k,X,Y);
G =  grad2D(k,m,Dx,n,Dy); % grad2DCurv(k,X,Y);

L = D*G;  %lap2D(k,IMax-1,Dx,JMax-1,Dy);
L = L + robinBC2D(k,m,Dx,n,Dy,1,0);

idx1 = numel(uic);
idy1 = numel(vic);
idx2 = idx1 + idy1;


dt = 0.1;
rho = 1;
D = rho/dt*D;
G = -dt/rho*G;
%% Preparing t=0 arrays

iatu = (X(1:IMax-1,:)+X(2:IMax,:))*0.50;
jatu = (Y(1:IMax-1,:)+Y(2:IMax,:))*0.50;

iatv = (X(:,1:JMax-1)+X(:,2:JMax))*0.50;
jatv = (Y(:,1:JMax-1)+Y(:,2:JMax))*0.50;
    
uatu = interp2(X, Y, uic, iatu,jatu); %u-space
uatv = interp2(X, Y, uic, iatv,jatv);

vatv = interp2(X, Y, vic, iatv,jatv);
patv = interp2(X, Y, pic, iatv,jatv);
vatu = interp2(X, Y, vic, iatu,jatu);
patu = interp2(X, Y, pic, iatu,jatu);

uatp = interp2(X, Y, uic, Xp, Yp);
vatp = interp2(X, Y, vic, Xp, Yp);
patp = interp2(X, Y, pic, Xp, Yp);

u = uatp;
v = vatp;
p = patp;

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
    subplot(2,2,4)
    pcolor(Xp,Yp,cav)
    hold on
    quiver(Xp,Yp,u,v)
    hold off
    xlim([0 1])
    ylim([0 1])
    colorbar;
    shading interp;
    title(['vorticity Dt = ' num2str((time+1)*dt)])
%%
% Reshaping


Vv = reshape(vatv',[], 1);
Uv = reshape(uatv',[], 1);
Pv = reshape(patv',[], 1);

Vu = reshape(vatu',[], 1);
Uu = reshape(uatu',[], 1);
Pu = reshape(patu',[], 1);

Vp = reshape(v',[], 1);
Up = reshape(u',[], 1);
Pp = reshape(p',[], 1);

cs = 1447.0;


tf = 0.1;
nf = tf/dt ;
beta = 1/cs;
Re = 750;
IRe = 1.0/Re;
nu = 0.0014;
time_steps = linspace(0,tf,nf);

[uatu,vatv,p] = applyboundaries2D(uatu,vatv,p); 
%uatu(end,:) = 1.0;
%%
u = interp2(iatu, jatu, uatu, Xp, Yp); %scaterring interpolant
v = interp2(iatv, jatv, vatv, Xp, Yp);

u(end,:) = u(end-1,:); %Needs one more
u(1,:)   = u(2,:);
v(:,end) = v(:,end-1);
v(:,1)   = v(:,2);

uatv = interp2(Xp,Yp, u, iatv, jatv);
vatu = interp2(Xp,Yp, v, iatu, jatu);

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
    subplot(2,2,4)
    pcolor(Xp,Yp,cav)
    hold on
    quiver(Xp,Yp,u,v)
    hold off
    xlim([0 1])
    ylim([0 1])
    colorbar;
    shading interp;
    title(['vorticity Dt = ' num2str((time+1)*dt)])

%Dy = Xp(1,3)-Xp(1,2);
%Dx = Dy;

%Dx = 1.0; 
%Dy = 1.0; 



%% Boussinesq Pred-Corr formulation:
%for time = 1:size(time_steps,2)

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

            
%             du_dx = (uatu(i,j+1)-uatu(i,j-1))/(2*Dx);
%             du_dy = (vatu(i+1,j)-vatu(i-1,j))/(2*Dy);
%             udu_dx  = uatu(i,j)*du_dx;
%             vdu_dy  = uatu(i,j)*du_dy;  
%             
            %udu_dx = kawamura2D(u,u,Dx,i,j,'x');
            %udv_dy = kawamura2D(u,v,Dy,i,j,'y');
            
            %RHS(i,j) = -udu_dx-vdu_dy + nu*(d2u_dy2+d2u_dx2);
            
            u_star(i,j) = uatu(i,j) - dt*(udu_dx+vdu_dy - IRe*(d2u_dy2+d2u_dx2));
        end
    end
    
    %u_star(:,end) = u_star(:,end-1);
    %u_star(end-1:end,end) = u_star(end-2:end-1,end-1);
    
    
    
    %RHS = D*[ Uu.*Uu - txx ; Uv.*Vv - txy_v ];
%     RHS = D*[ Uu - txx ; Uv - txy_v ];
%     RHS = reshape(RHS,[IMax+1 JMax+1])';
%     RHS = applyRHSboundaries2D(RHS,IMax,JMax);
%     u_star = SSPRK102D(u,RHS,dt);
    
    %Predict v*
     
    for j = 2:JMax-2
        for i = 2:IMax-1
    
            d2v_dy2 = (vatv(i-1,j)-2*vatv(i,j)+vatv(i+1,j))/(Dy^2);
            d2v_dx2 = (vatv(i,j-1)-2*vatv(i,j)+vatv(i,j+1))/(Dx^2);
                        
            dv_dx = (vatv(i,j+1)-vatv(i,j-1))/(2*Dx);
            dv_dy = (vatv(i+1,j)-vatv(i-1,j))/(2*Dy);       
            udv_dx  = uatv(i,j)*dv_dx;
            vdv_dy  = vatv(i,j)*dv_dy;
                    
%             dv_dx = (uatv(i,j+1)-uatv(i,j-1))/(2*Dx);
%             dv_dy = (vatv(i+1,j)-vatv(i-1,j))/(2*Dy);       
%             udv_dx  = vatv(i,j)*dv_dx;
%             vdv_dy  = vatv(i,j)*dv_dy;
                        
            %udv_dx = kawamura2D(u,v,Dx,i,j,'x');
            %vdv_dy = kawamura2D(v,v,Dy,i,j,'y');
           
            %RHS(i,j) = -udv_dx-vdv_dy + nu*(d2v_dy2+d2v_dx2);
           
            v_star(i,j) = vatv(i,j) - dt*(udv_dx+vdv_dy - IRe*(d2v_dy2+d2v_dx2)) -buoy(i,j);
           
        end
    end


%     RHS = D*[ Uu.*Vu - txy_u  ; Vv.*Vv - tyy ];
%     RHS = D*[ Vu - txy_u  ; Vv - tyy ];
%     RHS = reshape(RHS,[IMax+1 JMax+1])';
%     RHS = applyRHSboundaries2D(RHS,IMax,JMax);
%     v_star = SSPRK102D(v,RHS,dt);
    
    [u_star,v_star,~] = applyboundaries2Dstar(u_star,v_star,p);
    %u_star(end,:) = 1.0 ;

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
    
%%
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
    
    [uatu,vatv,~] = applyboundaries2D(uatu,vatv,p);
    %uatu(end,:) = 1.0;
    
    u = interp2(iatu, jatu, uatu, Xp, Yp); %scattered interpolant
    v = interp2(iatv, jatv, vatv, Xp, Yp);
    p = reshape(Pp,[IMax+1 JMax+1])';
        
    u(end,:) = u(end-1,:);
    u(1,:)   = u(2,:);
    v(:,end) = v(:,end-1);
    v(:,1)   = v(:,2);

    uatv = interp2(Xp,Yp, u, iatv, jatv);
    vatu = interp2(Xp,Yp, v, iatu, jatu);
     
     
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
    
            d2T_dy2 = (Tatc(i-1,j)-2*Tatc(i,j)+Tatc(i+1,j))/(Dy^2);
            d2T_dx2 = (Tatc(i,j-1)-2*Tatc(i,j)+Tatc(i,j+1))/(Dx^2);
                        
            dT_dx = (Tatc(i,j+1)-Tatc(i,j-1))/(2*Dx);
            dT_dy = (Tatc(i+1,j)-Tatc(i-1,j))/(2*Dy);       
            udT_dx  = uatv(i,j)*dT_dx;
            vdT_dy  = vatv(i,j)*dT_dy;
                    
           
            Tatc(i,j) = Tatc(i,j) - dt*(udT_dx+vdT_dy - IRe*(d2T_dy2+d2T_dx2));
           
        end
    end
    
    Tatc(end,:) = Tatc(end-1,:);
    Tatc(1,:)   = Tatc(2,:);
    Tatc(:,end) = Tatc(:,end-1);
    Tatc(:,1)   = Tatc(:,2);
    
    dens = rho_ref + drho_dT*(Tatc*T0-T_ref);  %EOS  
    
    buoy(1,:) = gforce*LStar*(dens(1,:)-rho_P0)/(rho_P0*(UStar^2));

    for j=2:JMax
        for i=1:IMax-1

        buoy(i,j) = gforce*LStar*(0.5*(dens(i,j)+dens(i+1,j))-rho_P0)/(rho_P0*(UStar^2));

        end
    end
    
    figure(2);
    u_v = squeeze(u(:,:));
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
    u_3 = squeeze(v(:,:));
    pcolor(u_3);
    title('v')
    colorbar;
    shading interp;
    subplot(2,2,4)
    pcolor(Xp,Yp,cav)
    hold on
    quiver(Xp,Yp,u,v)
    hold off
    xlim([0 1])
    ylim([0 1])
    colorbar;
    shading interp;
    title(['vorticity Dt = ' num2str((time+1)*dt)])


%end


