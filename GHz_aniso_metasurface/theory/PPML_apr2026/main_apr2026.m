clearvars;
close all; 
clc;

% LiNbO3 pol rotating metasurface
addpath('PPML_2d_basic_general_part20260401')

% gemoetric parameters
a = 800;       % periodicity
rad_hole = 200; % hole radius
dslab = 300;    % slab thickness

a1 = [a,0];  % square lattice
a2 = [0,a];

extra.plotflag = 0;  % set to 1 if you want to show the computational Fourier lattice
extra.plotfignum = 40;
extra.plotmn = 1;
trunc.scheme = 'parall';   
trunc.rad1 = 4;     % half number of computational spatial harmonics in first lattice direction
trunc.rad2 = 4;     % "    "    "    computational            "         second    "     "

N1 = 256; N2 = 256; % real space discretization grid (for FTs)
lattice  = create_lattice(a1,a2,N1,N2,trunc,extra); clearvars extra;
%%

tmpth = linspace(0,360,100);       % create a polygon that describes the circular hole
polygon(:,1) = rad_hole*cosd(tmpth) + a/2;
polygon(:,2) = rad_hole*sind(tmpth) + a/2;
clearvars tmp*

epsLN = diag([2.21^2, 2.14^2, 2.21^2]); %in X'Y'Z' ref. syst. (see notes)
psi = 30; % rotation in the wafer plane (see notes) 
R = [cosd(psi), sind(psi),0; -sind(psi), cosd(psi),0; 0,0,1];
for i = 1:3
for j = 1:3
tmp = 0;
    for ip = 1:3
    for jp = 1:3
    tmp = tmp + epsLN(ip,jp)*R(i,ip)*R(j,jp);
    end
    end
    epsLNrot(i,j) = tmp;
end
end

epsA.xx = epsLNrot(1,1); epsA.xy = epsLNrot(1,2);  % rotated LiNbO3 is material A 
epsA.yx = epsLNrot(2,1); epsA.yy = epsLNrot(2,2);
epsA.z  = epsLNrot(3,3);

epsB.xx = 1; epsB.xy = 0; % air is material B
epsB.yx = 0; epsB.yy = 1;
epsB.z = 1;

epsmask = create_epsmask_from_polygon(lattice,polygon,epsA,epsB);  % create the real space permittivity mask

plotpars.plotflag = 1;
plotpars.plotfignum = 50;
plotpars.N1t = 3; plotpars.N2t = 3;
plot_epsmask(lattice,epsmask,plotpars);           clearvars plotpars;

fteps = create_fteps(lattice,epsmask);     % create the fourier transforms as needed in RCWA               
%%

%              superstr.    |  patt. slab    |  subst.       
struc.d       = [0            dslab               0    ]; % 
struc.fts     = {             fteps                    }; % material A 
struc.sigma   = [           0                0         ];    

struc.epssup = 1;                                
struc.epssub = 1;

lambdav = [1200:2:2000]; % use nanometers because the pattern dimensions were given in nm
for ilambda = 1:length(lambdav)
    ilambda
lambda = lambdav(ilambda);
theta  = 1.01; % never set zero
phi = 90;

globpars.k0   = 2*pi/lambda;         % wavevector

inlayer = 'sup';
switch inlayer
    case 'sup'  
    globpars.kparx = globpars.k0 * sqrt(struc.epssup) * sind(theta)*cosd(phi);
    globpars.kpary = globpars.k0 * sqrt(struc.epssup) * sind(theta)*sind(phi);
    case 'sub'
    globpars.kparx = globpars.k0 * sqrt(struc.epssub) * sind(theta)*cosd(phi);
    globpars.kpary = globpars.k0 * sqrt(struc.epssub) * sind(theta)*sind(phi);
end

wm = solve_2d_basic_general(lattice,struc,globpars); % calculates the working matrices

extra.plotflag = 0;
difdata = DJM_2d_basic_general(wm,lattice,struc,globpars,extra); % calculates the diffraction Jones matrices

T{ilambda} = difdata.Jsupsub{1,1}; % the Jones matrix for transmission from superstrate to substrate
end
clearvars wm;

%%

pols = [1;0]; % s-pol. if phi = 90, this is x-pol. 
polp = [0;1]; % p-pol. if phi = 90, this is y-pol. 
for ilambda = 1:length(lambdav)
Txx(ilambda) = abs( (pols')*(T{ilambda}*pols) )^2;  % inietto x e  analizzo x 
Tyy(ilambda) = abs( (polp')*(T{ilambda}*polp) )^2;  % inietto y e  analizzo y  
Txy(ilambda) = abs( (polp')*(T{ilambda}*pols) )^2;  % inietto x e  analizzo y
Tyx(ilambda) = abs( (pols')*(T{ilambda}*polp) )^2;  % inietto y e  analizzo x
end
figure(1)
subplot(211)
plot(lambdav,Txx,'k'); hold on;
plot(lambdav,Tyy,'b'); hold on;
legend("Txx","Tyy")
subplot(212)
plot(lambdav,Txy,'c'); hold on;
plot(lambdav,Tyx,'m--'); hold on;
legend("Txy","Tyx")