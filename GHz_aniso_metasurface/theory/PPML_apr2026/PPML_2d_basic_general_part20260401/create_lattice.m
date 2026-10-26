function lattice = create_lattice(a1,a2,N1,N2,trunc,extra)

% Builds the quntities relative to direct and reciprocal lattices.
%
%    ====== INPUT 
%
% a1, a2 -> direct space unit cell vectors [each 1x2 real]
% N1, N2 -> direct space unit cell discretization [each integer, better if it is a power of 2 for subsequent FFT]
% trunc -> truncation scheme details.
%             * for circular trunc scheme, set trunc.scheme = 'circ' 
%                                          and trunc.rad    = [an integer] which is the truncation radius
%             * for rectangular truncation scheme, set trunc scheme = 'parall'
%                                                  and trunc.rad1 = [integer]; trunc.rad2 = [integer];
% extra -> extra.plotflag = [1 or 0] to plot or not to plot
%          extra.plotfignum = [integer], the figure number                                                  
%          extra.plotmn = [1 or 0] to plot the (m,n) indices for each reciprocal space basis element
%
%
%    ====== OUTPUT
%
% lattice.a1 = a1;  % direct lattice basis vector. [1x2 real]
% lattice.a2 = a2;
% lattice.g1 = g1;  % rec lattice basis vector. [1x2 real]
% lattice.g2 = g2;  
% lattice.mn = mn;  % list of m&n indices for all points in truncated rec lattice. [NPW x 2 matrix] 
% lattice.G = G;    % list of G   vectors for all points in truncated rec lattice. [NPW x 2 matrix] 
% lattice.N         % the number of harmonics (a.k.a. npw)
% lattice.P = P;    % list of direct space grid points.  [(N1xN2) x 2 matrix]
% lattice.mnDir = mnDir; % list of m&n indices for direct space grid points
% lattice.N1 = N1;
% lattice.N2 = N2;


%%%%%%%%%% reciprocal lattice primitive vectors
om = dot( [a1,0] , cross( [a2,0] , [0,0,1] ) );
tmp = 2*pi/om*cross( [a2,0] , [0,0,1] ); 
g1 = tmp(1:2);
tmp = 2*pi/om*cross( [0,0,1], [a1,0]  );
g2 = tmp(1:2);
clearvars tmp*;



%%% reciprocal lattice basis construction

switch trunc.scheme

case 'circ'
trunc_ord = trunc.rad;

trunc_rad = sqrt(norm(g1)*norm(g2))*(trunc_ord+0.001); % truncation radius for PW basis

halfnpw1t = ceil( 1.1*trunc_rad*norm(g2)/norm(cross( [g1,0], [g2,0] )) );
halfnpw2t = ceil( 1.1*trunc_rad*norm(g1)/norm(cross( [g1,0], [g2,0] )) ); % temporary. 

% create a parallelogramic grid based on N1t, N2t truncation
% and chooses only the vectors that lie within the trunc_rad
iprog = 1;
tmpG = [];
tmpmn = [];
for i1 = -halfnpw1t:halfnpw1t
    for i2 = -halfnpw2t:halfnpw2t
        tmp = i1*g1 + i2*g2;
        if norm(tmp) < trunc_rad
        tmpG(iprog,:) = tmp;
        tmpmn(iprog,:) = [i1,i2];
        iprog = iprog+1;
        end
    end
end
clearvars tmp i* halfnpw1t halfnpw2t;

% sort the vectors
[~,I] = sort(vecnorm(tmpG,2,2));
for iprog = 1:length(I)
G(iprog,:) = tmpG(I(iprog),:);
mn(iprog,:) = tmpmn(I(iprog),:);
end
clearvars tmp* I* i* 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'parall'

halfnpw1 = trunc.rad1;
halfnpw2 = trunc.rad2;

iprog = 1;
tmpmn = [];
for i1 = -halfnpw1:halfnpw1
    for i2 = -halfnpw2:halfnpw2
        tmp = i1*g1 + i2*g2;
        tmpG(iprog,:) = tmp;
        tmpmn(iprog,:) = [i1,i2];
        iprog = iprog+1;
    end
end
clearvars i*

% sort the vectors
[~,I] = sort(vecnorm(tmpG,2,2));
for iprog = 1:length(I)
G(iprog,:) = tmpG(I(iprog),:);
mn(iprog,:) = tmpmn(I(iprog),:);
end
clearvars tmp* I* i* 


end

%%%%%% end reciprocal lattice construction






%%%%% begin direct lattice grid construction
P = zeros(N1*N2,2);
mnDir = zeros(N1*N2,2);
iprog = 1;
for i1 = 1:N1
    for i2 = 1:N2
    P(iprog,:) = a1*(i1-1)/N1 + a2*(i2-1)/N2;
    mnDir(iprog,:) = [i1,i2];
    iprog = iprog+1;
    end
end
clearvars i*


%%%% (optionally) plotting the rec. space points
if extra.plotflag
figure(extra.plotfignum)
scatter(G(:,1),G(:,2),'.'); hold on; axis equal
plot([0, g1(1)],[0,g1(2)],'k-','linewidth',2)
plot([0, g2(1)],[0,g2(2)],'k-','linewidth',2)
switch trunc.scheme
case 'circ'
    plot(cosd(0:5:360)*trunc_rad,sind(0:5:360)*trunc_rad,'k-')
case 'parall'

end
if extra.plotmn 
for i = 1:size(G,1)
text(G(i,1)+g1(1)/10+g2(1)/10,G(i,2)+g1(2)/10+g2(2)/10,[num2str(mn(i,1)),' ',num2str(mn(i,2))]);
text(G(i,1)+g1(1)/10+g2(1)/10,G(i,2)+g1(2)/3+g2(2)/3,['I = ',num2str(i)]);
end
end
end



%%% returning the values to output structure

lattice.a1 = a1;  % direct lattice  basis vector. 1x2 vector
lattice.a2 = a2;
lattice.g1 = g1;  % rec lattice basis vector. 1x2 vector 
lattice.g2 = g2;  
lattice.mn = mn;  % list of m&n indices for all points in truncated rec lattice. Npw x 2 matrix
lattice.G = G;    % list of G   vectors for all points in truncated rec lattice. Npw x 2 matrix
lattice.N = size(G,1); % the number of harmonic (a.k.a. npw)
lattice.P = P;    % list of direct space grid points.  (N1xN2) x 2
lattice.mnDir = mnDir;
lattice.N1 = N1;
lattice.N2 = N2;