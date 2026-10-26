function plot_epsmask(lattice,epsmask,plotpars)

% ========= INPUT
%
% lattice: stucture as created by create_lattice(...)
% epsmask: is a structure containing fields  epsmask.xx
%                                           epsmask.xy
%                                            etc.
%                        each of which is a [Nx x Ny complex matrix] 
%                        containing the values of permittivity 
%                        at the real space grid points 
%
% plotpars.plotfignum: figure number
% plotpars.N1t: controls the number of repetitions of unit cell along a1 axis for visualization
% plotpars.N2t: controls the number of repetitions of unit cell along a2 axis for visualization


fnv = {'xx','xy','yx','yy','z'}; % list of names from tensor components of permittivity
spp = [1,2,4,5,3]; % list of subplot positions

for iprog = 1:size(lattice.P,1)
    X( lattice.mnDir(iprog,1), lattice.mnDir(iprog,2) ) = lattice.P( iprog, 1 );
    Y( lattice.mnDir(iprog,1), lattice.mnDir(iprog,2) ) = lattice.P( iprog, 2 );
end

figure(plotpars.plotfignum)
for ifn = 1:5
subplot(2,3,spp(ifn))


mx = 0;Mx = 0;
my = 0;My = 0;
for i1 = -plotpars.N1t:plotpars.N1t
    for i2 = -plotpars.N2t:plotpars.N2t
        tmp = i1*lattice.a1 + i2*lattice.a2;
        Xplot = X + tmp(1);    
        Yplot = Y + tmp(2);
        h = pcolor(Xplot,Yplot,abs(epsmask.(fnv{ifn}))); hold on;
        h.EdgeColor = 'none';
        
        if (i1 == 0 && i2 == 0) || (i1 == 0 && i2 == 1) 
        plot([tmp(1),lattice.a1(1)+tmp(1)],[tmp(2),lattice.a1(2)+tmp(2)],'k-','linewidth',2)
        end
        if (i1 == 0 && i2 == 0) || (i1 == 1 && i2 == 0)
        plot([tmp(1),lattice.a2(1)+tmp(1)],[tmp(2),lattice.a2(2)+tmp(2)],'k-','linewidth',2)
        end

mx = min(mx,tmp(1));
Mx = max(Mx,tmp(1));
my = min(my,tmp(2));
My = max(My,tmp(2));
        
    end
end

axis([mx Mx my My]);
axis equal;
colorbar;
title(['eps ',fnv{ifn}])
end




