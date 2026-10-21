function [Txx, Tyy, Txy, Tyx] = simula_membrana(a, rad_hole, dslab, lambdav)
% Same simulation as main_apr2026.m (LiNbO3 slab with a square lattice of circular holes), as a function of the
% geometry and without plots. All lengths in nm.
%   a        periodicity of the square lattice
%   rad_hole hole radius
%   dslab    slab thickness
%   lambdav  wavelengths (vector)
% Returns the transmissions Txx, Tyy (co-pol) and Txy, Tyx (cross-pol), row vectors with one value per wavelength.

addpath(fullfile(fileparts(mfilename('fullpath')), 'PPML_2d_basic_general_part20260401'));

a1 = [a,0];  % square lattice
a2 = [0,a];

extra.plotflag = 0;
trunc.scheme = 'parall';
trunc.rad1 = 4;     % half number of computational spatial harmonics in first lattice direction
trunc.rad2 = 4;
N1 = 256; N2 = 256; % real space discretization grid (for FTs)
lattice = create_lattice(a1,a2,N1,N2,trunc,extra);

tmpth = linspace(0,360,100);       % polygon that describes the circular hole
polygon = [rad_hole*cosd(tmpth(:)) + a/2, rad_hole*sind(tmpth(:)) + a/2];

epsLN = diag([2.21^2, 2.14^2, 2.21^2]); % in X'Y'Z' ref. syst.
psi = 30;                               % rotation in the wafer plane
R = [cosd(psi), sind(psi),0; -sind(psi), cosd(psi),0; 0,0,1];
epsLNrot = R * epsLN * R.';

epsA.xx = epsLNrot(1,1); epsA.xy = epsLNrot(1,2);  % rotated LiNbO3 is material A
epsA.yx = epsLNrot(2,1); epsA.yy = epsLNrot(2,2);
epsA.z  = epsLNrot(3,3);

epsB.xx = 1; epsB.xy = 0; % air is material B
epsB.yx = 0; epsB.yy = 1;
epsB.z = 1;

epsmask = create_epsmask_from_polygon(lattice,polygon,epsA,epsB);
fteps = create_fteps(lattice,epsmask);

%              superstr.    |  patt. slab    |  subst.
struc.d       = [0            dslab               0    ];
struc.fts     = {             fteps                    };
struc.sigma   = [           0                0         ];
struc.epssup = 1;
struc.epssub = 1;

theta = 1.01;   % never set zero
phi = 90;
pols = [1;0];   % s-pol. if phi = 90, this is x-pol.
polp = [0;1];   % p-pol. if phi = 90, this is y-pol.

Txx = zeros(1,numel(lambdav)); Tyy = Txx; Txy = Txx; Tyx = Txx;
for ilambda = 1:numel(lambdav)
    globpars.k0 = 2*pi/lambdav(ilambda);
    globpars.kparx = globpars.k0 * sqrt(struc.epssup) * sind(theta)*cosd(phi);
    globpars.kpary = globpars.k0 * sqrt(struc.epssup) * sind(theta)*sind(phi);
    wm = solve_2d_basic_general(lattice,struc,globpars);
    extra.plotflag = 0;
    difdata = DJM_2d_basic_general(wm,lattice,struc,globpars,extra);
    T = difdata.Jsupsub{1,1};  % Jones matrix for transmission from superstrate to substrate
    Txx(ilambda) = abs((pols')*(T*pols))^2;  % inject x, analyse x
    Tyy(ilambda) = abs((polp')*(T*polp))^2;  % inject y, analyse y
    Txy(ilambda) = abs((polp')*(T*pols))^2;  % inject x, analyse y
    Tyx(ilambda) = abs((pols')*(T*polp))^2;  % inject y, analyse x
end
end
