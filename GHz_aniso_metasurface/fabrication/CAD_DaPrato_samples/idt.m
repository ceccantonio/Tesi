function [y1,y2] = idt(pos,b,l1,l2,N,type,layer1,layer2)


switch type
    case 1
   c = -1;
    case 2
   c = 1;
end


% pos -> [x,y] 
%x          % vertical axis of IDT area
%y=y+c*7;   % top (bottom) of IDT area (depending on type = 1 or 2)
                  % + 7 um di tolleranza
                  
pos = pos + [0,c*7];
x = pos(1); y = pos(2);

% internal idt finger width
wi = b/5;
% external idt finger width
we = b/4;
eps = 0.5; % overlap


global gs

holec = cell(2*N,1);

k = 1;
for n = 1:N

holec{k} = [x-l2/2-eps,    y+c*((n-1)*b + (we-wi)/2);...
            x-l2/2-eps+l2, y+c*((n-1)*b + (we-wi)/2);...
            x-l2/2-eps+l2, y+c*((n-1)*b + (we-wi)/2)+c*wi;...
            x-l2/2-eps,    y+c*((n-1)*b + (we-wi)/2)+c*wi;...
            x-l2/2-eps,    y+c*((n-1)*b + (we-wi)/2)];
          
k = k+1;

holec{k} = [x-l2/2+eps,    y+c*((n-1)*b + b/2 + (we-wi)/2);...  
            x-l2/2+eps+l2, y+c*((n-1)*b + b/2 + (we-wi)/2);...  
            x-l2/2+eps+l2, y+c*((n-1)*b + b/2 + (we-wi)/2)+c*wi;...  
            x-l2/2+eps,    y+c*((n-1)*b + b/2 + (we-wi)/2)+c*wi;...  
            x-l2/2+eps,    y+c*((n-1)*b + b/2 + (we-wi)/2)]; 

k = k+1;
end
gs(end+1) = gds_element('boundary', 'xy',holec, 'layer',layer1);



holec = cell(2*N,1);

k = 1;
for n = 1:N

holec{k} = [x-l1/2,           y+c*(n-1)*b;...
            x-l1/2+(l1-l2)/2, y+c*(n-1)*b;...
            x-l1/2+(l1-l2)/2, y+c*(n-1)*b+c*we;...
            x-l1/2,           y+c*(n-1)*b+c*we;...
            x-l1/2,           y+c*(n-1)*b];
          
k = k+1;

holec{k} = [x+l2/2,           y+c*(n-1)*b + c*b/2;...
            x+l2/2+(l1-l2)/2, y+c*(n-1)*b + c*b/2;...
            x+l2/2+(l1-l2)/2, y+c*(n-1)*b + c*b/2+c*we;...
            x+l2/2,           y+c*(n-1)*b + c*b/2+c*we;...
            x+l2/2,           y+c*(n-1)*b + c*b/2]; 

k = k+1;
end
gs(end+1) = gds_element('boundary', 'xy',holec, 'layer',layer2);

y1 = y;
y2 = y1+c*N*b;
