function [output] = main_LNO_super(aa,rr,ds,sigma,flag,nsc,lambdav)

% Same as main_LNO_func, but the unit cell is a supercell of nsc x nsc holes (nsc = 10 if not given): every hole is
% shifted from its equilibrium position by a gaussian noise of standard deviation sigma (nm), the same in x and y.
% sigma = 0 gives the perfect lattice. lambdav = wavelengths to simulate in nm (1520:0.5:1570 if not given). Heavy: the number of harmonics grows as (2*4*nsc+1)^2 (nsc = 3: 625, nsc = 5: 1681, nsc = 10: 6561).

% LiNbO3 pol rotating metasurface
addpath('PPML_2d_basic_general_part20260401')

% gemoetric parameters
a = aa;       % periodicity
rad_hole = rr; % hole radius
dslab = ds;    % slab thickness - etched
dslabh = 300 - dslab;   % slab thickness - not etched
if nargin < 6, nsc = 10; end
L = nsc*a;    % period of the supercell

a1 = [L,0];  % square supercell
a2 = [0,L];

extra.plotflag = 0;  % set to 1 if you want to show the computational Fourier lattice
extra.plotfignum = 40;
extra.plotmn = 1;
trunc.scheme = 'parall';   
trunc.rad1 = 4*nsc;     % half number of computational spatial harmonics in first lattice direction
trunc.rad2 = 4*nsc;     % "    "    "    computational            "         second    "     "

N1 = 2^nextpow2(100*nsc); N2 = N1; % real space discretization grid (for FTs)
lattice  = create_lattice(a1,a2,N1,N2,trunc,extra); clearvars extra;
%%

rng(1);                                   % the same defects at every call, so that results at different a, r can be compared
[ix, iy] = meshgrid(0:nsc-1);
centri = [(ix(:)+0.5)*a, (iy(:)+0.5)*a] + sigma*randn(nsc^2, 2);   % hole centres: equilibrium position + gaussian shift (nm)


epsLN = diag([2.21^2, 2.14^2, 2.21^2]); %in X'Y'Z' ref. syst. (see notes)
psi = 0; % rotation in the wafer plane (see notes) 
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

epsmask = maschera(lattice,centri,rad_hole,epsA,epsB);  % create the real space permittivity mask

plotpars.plotflag = 0;
plotpars.plotfignum = 50;
plotpars.N1t = 3; plotpars.N2t = 3;
%plot_epsmask(lattice,epsmask,plotpars);           clearvars plotpars;

fteps = create_fteps(lattice,epsmask);     % create the fourier transforms as needed in RCWA

epsmaskA = maschera(lattice,zeros(0,2),rad_hole,epsA,epsA);  % homogeneous mask: material A everywhere
ftepsA   = create_fteps(lattice,epsmaskA);                           % FTs of the unpatterned slab
%%

%              superstr.    |  patt. slab   unpatt. slab |  subst.       
struc.d       = [0            dslab            dslabh          0    ]; % 
struc.fts     = {             fteps            ftepsA               }; % material A
struc.sigma   = [           0                0          0           ];    

struc.epssup = 1;                                
struc.epssub = 1;

if nargin < 7, lambdav = 1520:0.5:1570; end % use nanometers because the pattern dimensions were given in nm

for ilambda = 1:length(lambdav)
    ilambda
lambda = lambdav(ilambda);
theta  = 0.01; % never set zero
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
if flag.plotS
figure(1)
subplot(211)
plot(lambdav,Txx,'k'); hold on;
plot(lambdav,Tyy,'b'); hold on;
legend("Txx","Tyy")
subplot(212)
plot(lambdav,Txy,'c'); hold on;
plot(lambdav,Tyx,'m--'); hold on;
legend("Txy","Tyx")
end
datasave = [lambdav',Txx',Tyy',Txy',Tyx'];
output = datasave;
end

function epsmask = maschera(lattice,centri,r,epsA,epsB)   % permittivity on the grid: material B inside the circles of radius r around centri, material A outside
dentro = false(size(lattice.P,1),1);
for i = 1:size(centri,1)
    dentro = dentro | (lattice.P(:,1)-centri(i,1)).^2 + (lattice.P(:,2)-centri(i,2)).^2 < r^2;
end
M = false(lattice.N1,lattice.N2);
M(sub2ind([lattice.N1 lattice.N2],lattice.mnDir(:,1),lattice.mnDir(:,2))) = dentro;
for f = {'xx','xy','yx','yy','z'}
    epsmask.(f{1}) = epsA.(f{1})*~M + epsB.(f{1})*M;
end
end
