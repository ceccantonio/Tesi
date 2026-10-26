function epsmask = create_epsmask_from_polygon(lattice,polygon,epsA,epsB)

% lattice: stucture as created by create_lattice(...)
% polygon: matrix [x1,y1;
%                  x2,y2;
%				   .....;
%				   xN,yN]
%          of the vertices of the polygon that contains material B.
%          To define multiply connected polygons, separate them by NaNs:
%                  [x1_1 ,y1_1;
%                   x1_2 ,y1_2;
%                   .....    ;
%                   x1_N1,y1_N1;
%					NaN  ,NaN ;
%					x2_1 ,y2_1;
%                   x2_2 ,y2_2;
%                   .....    ;
%                   x2_N2,y2_N2]
%
% epsA.xx , epsA.xy, ... = value of permittivity components of material A [complex]
% ..idem with B


fnv = {'xx','xy','yx','yy','z'}; % list of names from tensor components of permittivity
    
In =  inpolygon(lattice.P(:,1),lattice.P(:,2),polygon(:,1),polygon(:,2)); % finding single indices of points within polygon

for ifn = 1:5
tmpepsmask = ones(lattice.N1,lattice.N2)*epsA.(fnv{ifn}) ;
for iprog = 1:length(In)
    if In(iprog)
    tmpepsmask( lattice.mnDir(iprog,1), lattice.mnDir(iprog,2) ) = epsB.(fnv{ifn});
    end
end
epsmask.(fnv{ifn}) = tmpepsmask;
end

clearvars tmp*

