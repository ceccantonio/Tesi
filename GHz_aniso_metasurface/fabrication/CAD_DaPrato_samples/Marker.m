function Marker(pos,type,layer)
global gs
% pos -> [x,y] position of marker center
% type -> 0 or 1 

% marker line width
w = 20;
% marker side length
l = 50;

switch type
    case 0
hole =  [0,0; l,0; l,w; w,w; w,l; 0,l; 0,0] + pos;
gs(end+1) = gds_element('boundary', 'xy',hole, 'layer',layer);
hole = -[0,0; l,0; l,w; w,w; w,l; 0,l; 0,0] + pos;
gs(end+1) = gds_element('boundary', 'xy',hole, 'layer',layer);
    case 1
hole =  [0,0; -l,0; -l,w; -w,w; -w,l; 0,l; 0,0] + pos;
gs(end+1) = gds_element('boundary', 'xy',hole, 'layer',layer);
hole = -[0,0; -l,0; -l,w; -w,w; -w,l; 0,l; 0,0] + pos;
gs(end+1) = gds_element('boundary', 'xy',hole, 'layer',layer);
end