% The initial condition employed in the code Dyniso

v0 = 27.1893;% reference velocity

l = 55.0; % reference length scale
% l = 50.8;

a0 =  0.56102E+01;   a(1) = -0.11236E+01;   a(2) = -0.30961E+00;
a(3) =  0.33172E+00;   a(4) = -0.10959E+00;   a(5) = -0.22320E-01; 
a(6) =  0.66575E-02;

c1 = 2.*pi/l;
c2 = c1/(v0^2);

k = 1:0.1:1000;
arg = a0+a(1)*log(c1*k)+a(2)*log(c1*k).^2+a(3)*log(c1*k).^3+...
    a(4)*log(c1*k).^4+a(5)*log(c1*k).^5+a(6)*log(c1*k).^6;
Ek = c2*exp(arg);

figure;loglog(k,Ek)

%%
% The experiment results of Comte-Bellot & Corrsin (1971)
uref = 27.1893;
lref = 55.0/2./pi;

load('CBC_exp.mat');
k42 = k_42*lref;
E42 = E_42/(uref.^2*lref);
k98 = k_98*lref;
E98 = E_98/(uref.^2*lref);
k171 = k_171*lref;
E171 = E_171/(uref.^2*lref);

figure;loglog(k42,E42,'m+',k98,E98,'ro',k171,E171,'bs','LineWidth',2);hold on;
axis([1,400,10^-5,0.2])
grid on

% fig=figure;
% p = plot(yp/ystar,sqrt(Rxxmean)/ustar*ystar,...
%     yp/ystar,sqrt(Ryymean)/ustar*ystar,yp/ystar,sqrt(Rzzmean)/ustar*ystar);
% %p = plot(zp,up,zp,up);
% %p = plot(x,cf_fp,xl,cf_rib_1280);
% %p = plot(x,cf_case3,x,cf_case3_5p);
% p(1).LineWidth =2;
% p(1).Color = 'k';
% 
% p(2).Color = 'b';
% p(2).LineWidth = 2;
% p(2).LineStyle = '--';
% 
% p(3).Color = 'r';
% p(3).LineWidth = 2;
% p(3).LineStyle = '-.';
% %
% axis([0 180 0 0.4]);
%axis([0.4 1000 0 50]);
%axis([300 940 0.0008 0.0062]);
%axis([300 940 0.0008 0.0078]);
grid on
set(gca,'FontSize',20);
yticks([10^-5 10^-3 10^-1])
t=xlabel('$$\kappa$$','Interpreter','Latex','FontSize',24);
t=ylabel('$$E(\kappa)$$','Interpreter','Latex','FontSize',24);
%%
%
%%
ek42 = load('liutex_spectra.0');
ek98 = load('liutex_spectra.100');
ek171 = load('liutex_spectra.230');
loglog(ek42(:,1),ek42(:,2),'m-','LineWidth',2);
loglog(ek98(1:20,1),ek98(1:20,2),'r','LineWidth',2);
loglog(ek171(1:20,1),ek171(1:20,2),'b','LineWidth',2);
%%
print(gcf,'liutex_les','-dpng')

%%
ek42 = load('smagorinsky_spectra.0');
ek98 = load('smagorinsky_spectra.100');
ek171 = load('smagorinsky_spectra.230');
loglog(ek42(:,1),ek42(:,2),'m','LineWidth',2);
loglog(ek98(1:20,1),ek98(1:20,2),'r','LineWidth',2);
loglog(ek171(1:20,1),ek171(1:20,2),'b','LineWidth',2);
load CBC_512_98_filt2_spec.mat;
loglog(k_mag(1:30),e_k_mag(1:30),'color',[0.5 0.5 0.5],'LineWidth',1.5)
load CBC_512_171_filt2_spec.mat;
loglog(k_mag(1:30),e_k_mag(1:30),'color',[0.5 0.5 0.5],'LineWidth',1.5)

%%
print(gcf,'smagorinsky_les','-dpng')

%%
% comparison with the results by Rozema's DSM subgrid model
load('DSMresults_Rozema2015_fig1.mat')
sp98(:,1)=sp98(:,1)/55.88*55/2/pi;
sp171(:,1)=sp171(:,1)/55.88*55/2/pi;
sp98(:,2)=sp98(:,2)*55.88/55*2*pi;
sp171(:,2)=sp171(:,2)*55.88/55*2*pi;
%
%figure;
%loglog(k42,E42,'*',k98,E98,'o',k171,E171,'*');hold on;
loglog(sp98(:,1),sp98(:,2),'g--','LineWidth',1.5);
loglog(sp171(:,1),sp171(:,2),'g--','LineWidth',1.5);


