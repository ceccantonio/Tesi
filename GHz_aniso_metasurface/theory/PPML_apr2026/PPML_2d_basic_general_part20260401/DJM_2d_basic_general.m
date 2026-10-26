%%
function difdata = DJM_2d_basic_general(wm,lattice,struc,globpars,extra)

L = size(wm.q,3)-2;
N = lattice.N;

epssup = struc.epssup;
epssub = struc.epssub;

k0 = globpars.k0;
kparx = globpars.kparx;
kpary = globpars.kpary;

kx = lattice.G(:,1) + kparx;
ky = lattice.G(:,2) + kpary;

%%%%%%% determining the diff orders into superstrate & substrate
kparG = [kx,ky];


%%%%%%% determining the open diff orders into superstrate & substrate
iprog = 1;
IIsup = [];
for I = 1:N
    if ( norm(kparG(I,:)) < k0*sqrt(epssup) )
        IIsup(iprog) = I;
        iprog = iprog+1;
    end
end
difdata.IIsup = IIsup;

iprog = 1;
IIsub = [];
for I = 1:N
    if ( norm(kparG(I,:)) < k0*sqrt(epssub) )
        IIsub(iprog) = I;
        iprog = iprog+1;
    end
end
difdata.IIsub = IIsub;

Nsup = length(IIsup);
Nsub = length(IIsub);
difdata.Nsup = Nsup;
difdata.Nsub = Nsub;

%%%%%%%%%%% plotting diffraction pattern
if extra.plotflag
figure(extra.plotfignum)

subplot(121) % diff pattern in superstrate
plot([-1,1]*k0*sqrt(epssup), [0,0],'k'); hold on;
plot([0,0],[-1,1]*k0*sqrt(epssup),'k'); hold on;
scatter(kparG(IIsup,1),kparG(IIsup,2),'.'); hold on; axis equal
plot(cosd(0:5:360)*k0*sqrt(epssup),sind(0:5:360)*k0*sqrt(epssup),'k-')

textspac = min(k0*sqrt(epssup)/5,abs(lattice.g2(2))/5);
for i = 1:Nsup
I = IIsup(i);
text(kparG(I,1),kparG(I,2)+textspac,...
    ['[',num2str(lattice.mn(I,1)),' ',num2str(lattice.mn(I,2)),']']);
text(kparG(I,1),kparG(I,2)+textspac*2,...
    ['I = ',num2str(I)]);
text(kparG(I,1),kparG(I,2)+textspac*3,...
    ['Isup = ',num2str(i)]);
end

subplot(122) % diff pattern in substrate
plot([-1,1]*k0*sqrt(epssub), [0,0],'k'); hold on;
plot([0,0],[-1,1]*k0*sqrt(epssub),'k'); hold on;
scatter(kparG(IIsub,1),kparG(IIsub,2),'.'); hold on; axis equal
plot(cosd(0:5:360)*k0*sqrt(epssub),sind(0:5:360)*k0*sqrt(epssub),'k-')

textspac = min(k0*sqrt(epssub)/5,abs(lattice.g2(2))/5);
for i = 1:Nsub
I = IIsub(i);
text(kparG(I,1),kparG(I,2)+textspac,...
    ['[',num2str(lattice.mn(I,1)),' ',num2str(lattice.mn(I,2)),']']);
text(kparG(I,1),kparG(I,2)+textspac*2,...
    ['I = ',num2str(I)]);
text(kparG(I,1),kparG(I,2)+textspac*3,...
    ['Isub = ',num2str(i)]);
end
end


%%%%%%%%%%%%%%%%%%%%%% local polarization basis & matrices Uin, Uout
Uin  = zeros(2*(Nsup+Nsub));
Uout = zeros(2*(Nsup+Nsub));

%%%%% superstrate
for i = 1:Nsup
I = IIsup(i);

%%% cosine of the diffraction angle
difdata.abscosthsup(i) =  abs( dot([kparG(I,:), wm.q(I,I,1)] , [0,0,1]) / norm([kparG(I,:), wm.q(I,I,1)]) );
 
%%% input waves
shat = -cross([ kparG(I,:),0] , [0,0,1] );
shat = shat/norm(shat);

phat = cross( [kparG(I,:), wm.q(I,I,1)] ,shat);
phat = phat/norm(phat);

U(1:2,1) = shat(1:2); 
U(1:2,2) = phat(1:2); 
Uin((2*i-1):(2*i),(2*i-1):(2*i)) = U;
Uintilde((2*i-1):(2*i),(2*i-1):(2*i)) = U/( sqrt(sqrt(epssup))*sqrt( difdata.abscosthsup(i)) );

%%% output waves
shat = -cross([ kparG(I,:),0] , [0,0,1] );
shat = shat/norm(shat);

phat = cross( [kparG(I,:), -wm.q(I,I,1)] ,shat);
phat = phat/norm(phat);

U(1:2,1) = shat(1:2); 
U(1:2,2) = phat(1:2); 
Uout((2*i-1):(2*i),(2*i-1):(2*i)) = U;
Uouttilde((2*i-1):(2*i),(2*i-1):(2*i)) = U/( sqrt(sqrt(epssup))*sqrt( difdata.abscosthsup(i)) );

end 


%%%%% substrate
for i = 1:Nsub
I = IIsub(i);

%%% cosine of the diffraction angle
difdata.abscosthsub(i) = abs( dot([kparG(I,:), wm.q(I,I,end)] , [0,0,1]) / norm([kparG(I,:), wm.q(I,I,end)]) );

%%% input waves
shat = -cross([ kparG(I,:),0] , [0,0,1] );
shat = shat/norm(shat);

phat = cross( [kparG(I,:), -wm.q(I,I,end)] ,shat);
phat = phat/norm(phat);

U(1:2,1) = shat(1:2); 
U(1:2,2) = phat(1:2); 
Uin((2*i-1+2*Nsup):(2*i+2*Nsup), (2*i-1+2*Nsup):(2*i+2*Nsup)) = U;
Uintilde((2*i-1+2*Nsup):(2*i+2*Nsup), (2*i-1+2*Nsup):(2*i+2*Nsup)) = U/( sqrt(sqrt(epssub))*sqrt( difdata.abscosthsub(i) ) );

%%% output waves
shat = -cross([ kparG(I,:),0] , [0,0,1] );
shat = shat/norm(shat);

phat = cross( [kparG(I,:), wm.q(I,I,end)] ,shat);
phat = phat/norm(phat);

U(1:2,1) = shat(1:2); 
U(1:2,2) = phat(1:2); 
Uout((2*i-1+2*Nsup):(2*i+2*Nsup), (2*i-1+2*Nsup):(2*i+2*Nsup)) = U;
Uouttilde((2*i-1+2*Nsup):(2*i+2*Nsup), (2*i-1+2*Nsup):(2*i+2*Nsup)) = U/( sqrt(sqrt(epssub))*sqrt( difdata.abscosthsub(i) ) );

end 


% difdata.Uin  = Uin;
% difdata.Uout = Uout;


Win  = zeros(4*N);
Wout = zeros(4*N);
for I = 1:N
Win(2*(I-1)+1,   I) =  wm.A(N+I,  I,1) * exp(1i*wm.q(I,I,1)*struc.d(1));
Win(2*(I-1)+1, N+I) =  wm.A(N+I,N+I,1) * exp(1i*wm.q(I,I,1)*struc.d(1));
Win(2*(I-1)+2,   I) = -wm.A(  I,  I,1) * exp(1i*wm.q(I,I,1)*struc.d(1));
Win(2*(I-1)+2, N+I) = -wm.A(  I,N+I,1) * exp(1i*wm.q(I,I,1)*struc.d(1));

Win(2*N+2*(I-1)+1, 2*N+  I) = -wm.A(N+I,  I,L+2) * exp(1i*wm.q(I,I,L+2)*struc.d(L+2));
Win(2*N+2*(I-1)+1, 2*N+N+I) = -wm.A(N+I,N+I,L+2) * exp(1i*wm.q(I,I,L+2)*struc.d(L+2));
Win(2*N+2*(I-1)+2, 2*N+  I) =  wm.A(  I,  I,L+2) * exp(1i*wm.q(I,I,L+2)*struc.d(L+2));
Win(2*N+2*(I-1)+2, 2*N+N+I) =  wm.A(  I,N+I,L+2) * exp(1i*wm.q(I,I,L+2)*struc.d(L+2));

Wout(2*(I-1)+1, 2*N+  I) = -wm.A(N+I,  I,1) ;
Wout(2*(I-1)+1, 2*N+N+I) = -wm.A(N+I,N+I,1) ;
Wout(2*(I-1)+2, 2*N+  I) =  wm.A(  I,  I,1) ;
Wout(2*(I-1)+2, 2*N+N+I) =  wm.A(  I,N+I,1) ;

Wout(2*N+2*(I-1)+1,   I) =  wm.A(N+I,  I,L+2);
Wout(2*N+2*(I-1)+1, N+I) =  wm.A(N+I,N+I,L+2);
Wout(2*N+2*(I-1)+2,   I) = -wm.A(  I,  I,L+2);
Wout(2*N+2*(I-1)+2, N+I) = -wm.A(  I,N+I,L+2);
end

V = zeros(2*(Nsup+Nsub),4*N);
for Isup = 1:Nsup
V(2*(Isup-1)+1 , 2*(IIsup(Isup)-1)+1  ) = 1;
V(2*(Isup-1)+2 , 2*(IIsup(Isup)-1)+2  ) = 1;
end
for Isub = 1:Nsub
V(2*Nsup + 2*(Isub-1)+1 , 2*N + 2*(IIsub(Isub)-1)+1  ) = 1;
V(2*Nsup + 2*(Isub-1)+2 , 2*N + 2*(IIsub(Isub)-1)+2  ) = 1;
end

% difdata.Win = Win;
% difdata.Wout = Wout;
% difdata.V = V;

D  = Uout\( V*( (Wout*( ([wm.S1(:,:,L+2), wm.S2(:,:,L+2); wm.S3(:,:,1), wm.S4(:,:,1)]/Win )*(V') ) )*Uin ) );
Dt = Uouttilde\( V*( (Wout*( ([wm.S1(:,:,L+2), wm.S2(:,:,L+2); wm.S3(:,:,1), wm.S4(:,:,1)]/Win )*(V') ) )*Uintilde ) );
difdata.D = D;
difdata.Dt = Dt;

difdata.Jsupsup = {}; difdata.Jsupsub = {}; difdata.Jsubsup = {}; difdata.Jsubsub = {};
difdata.Jtsupsup = {}; difdata.Jtsupsub = {}; difdata.Jtsubsup = {}; difdata.Jtsubsub = {};
% sup -> sup Jones matrices
for Isup = 1:Nsup
for Ksup = 1:Nsup
    difdata.Jsupsup{Isup,Ksup} = D( 2*Ksup-1 : 2*Ksup , 2*Isup-1 : 2*Isup );
    difdata.Jtsupsup{Isup,Ksup} = Dt( 2*Ksup-1 : 2*Ksup , 2*Isup-1 : 2*Isup );
end
end

% sup -> sub Jones matrices
for Isup = 1:Nsup
for Ksub = 1:Nsub
    difdata.Jsupsub{Isup,Ksub} = D(2*Nsup+2*Ksub-1 : 2*Nsup+2*Ksub , 2*Isup-1 : 2*Isup );
    difdata.Jtsupsub{Isup,Ksub} = Dt(2*Nsup+2*Ksub-1 : 2*Nsup+2*Ksub , 2*Isup-1 : 2*Isup );
end
end

% sub -> sup Jones matrices
for Isub = 1:Nsub
for Ksup = 1:Nsup
    difdata.Jsubsup{Isub,Ksup} = D( 2*Ksup-1 : 2*Ksup , 2*Nsup+2*Isub-1 : 2*Nsup+2*Isub );
    difdata.Jtsubsup{Isub,Ksup} = Dt( 2*Ksup-1 : 2*Ksup , 2*Nsup+2*Isub-1 : 2*Nsup+2*Isub );
end
end

% sub -> sub Jones matrices
for Isub = 1:Nsub
for Ksub = 1:Nsub
    difdata.Jsubsub{Isub,Ksub} = D( 2*Nsup+2*Ksub-1 : 2*Nsup+2*Ksub , 2*Nsup+2*Isub-1 : 2*Nsup+2*Isub );
    difdata.Jtsubsub{Isub,Ksub} = Dt( 2*Nsup+2*Ksub-1 : 2*Nsup+2*Ksub , 2*Nsup+2*Isub-1 : 2*Nsup+2*Isub );
end
end