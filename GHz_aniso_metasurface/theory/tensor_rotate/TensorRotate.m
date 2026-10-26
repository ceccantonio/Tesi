clearvars; close all;

% rotation of a crystal in a fixed reference frame x-y-z
% p  and p4  are the tensors of the UNROTATED crystal
% pT and pT4 are the tensors of the ROTATED crystal.

p11 = 5.78; p12 = -1.01; p13 = -1.47; p33 = 5.02; p14 = -1.02; p44 = 17; % parametri di compliance del LiNbO3 ,
                                         %  imposti essere coincidenti con
                                         %  il sistema di rif. x-y-z.   

% loading the tensor in Voigt indexing
p = zeros(6);
p = [p11, p12, p13, p14, 0,0;
     p12, p11, p13, -p14, 0,0;
     p13, p13, p33, 0,0,0;
     p14, -p14, 0,p44, 0,0;
     0,0,0,0,p44,2*p14;
     0,0,0,0,2*p14,2*(p11-p12)];

% transforming the tensor into extended indexing
for i1 = 1:3
for i2 = 1:3
for i3 = 1:3
for i4 = 1:3
p4(i1,i2,i3,i4) = p(xyz2voigt(i1,i2),xyz2voigt(i3,i4));
end
end
end
end

% rotation in the plane defined by [100] and [001] 
theta = pi/4;
T = [cos(theta) , 0, sin(theta);
     0          , 1, 0         ;
     -sin(theta), 0, cos(theta)];


for i1 = 1:3
for i2 = 1:3
for i3 = 1:3
for i4 = 1:3
    tmp = 0;
    for ib1 = 1:3
    for ib2 = 1:3
    for ib3 = 1:3
    for ib4 = 1:3
    tmp = tmp + T(i1,ib1)*T(i2,ib2)*T(i3,ib3)*T(i4,ib4)*p4(ib1,ib2,ib3,ib4);
    end
    end
    end
    end
    pT4(i1,i2,i3,i4) = tmp;
end
end
end
end


for I1 = 1:6
for I2 = 1:6
    [i1,i2] = voigt2xyz(I1); 
    [i3,i4] = voigt2xyz(I2);
pT(I1,I2) = pT4(i1,i2,i3,i4);
end
end

% now pT and pT4 are the tensors of the rotated crystal 

