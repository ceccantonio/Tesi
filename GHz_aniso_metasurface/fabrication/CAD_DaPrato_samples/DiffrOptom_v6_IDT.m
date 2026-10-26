clearvars; clc;

% Initialize Cell Parameters
% addpath(genpath('C:\Users\Simone\Documents\Z_Pisa2\Software\GDSII Library')); 
addpath(genpath('C:\Users\alesp\Documents\MATLAB\GDSII Library')); 

DiffrOptom_v6_PositionInit
load Positions


PhCsize = 100; %lato del PhC [um] 
a = 1.32;   % non toccare
f = 0.3;    % solo a scopo rappresentativo

Lambdav = linspace(2.5,2.85,12); % da dati su Berlino 1 

% create a structure to hold elements
global gs;
gs = gds_structure('DiffrOptom_v6_IDT');   

% chips layout
%for ichip = 1:size(poschip,1)
for ichip = 1:1
tmp1 = chipsize(ichip,1);
tmp2 = chipsize(ichip,2);
hole = poschip(ichip,:) + [0,0; 0,tmp2; tmp1,tmp2; tmp1,0; 0,0];
gs(end+1) = gds_element('boundary', 'xy',hole, 'layer',0);
end
clearvars tmp* i*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% chips
%for ichip = 1:size(poschip,1)
for ichip = 1:1
%%%% markers
Marker( poschip(ichip,:) + [1000,1000] ,0 , 25);
Marker( poschip(ichip,:) + [7000,1000] ,1 , 25);
Marker( poschip(ichip,:) + [1000,7000] ,1 , 25);
Marker( poschip(ichip,:) + [7000,7000] ,0 , 25);

postxt = poschip(ichip,:) + [1000,1200];
posstr = convertStringsToChars(num2str(ichip));
gs(end+1) = gdsii_ptext(posstr, postxt, 100, 25, 0);
postxt = poschip(ichip,:) + [7000,1200];
posstr = convertStringsToChars(num2str(ichip));
gs(end+1) = gdsii_ptext(posstr, postxt, 100, 25, 0);
postxt = poschip(ichip,:) + [1000,7200];
posstr = convertStringsToChars(num2str(ichip));
gs(end+1) = gdsii_ptext(posstr, postxt, 100, 25, 0);
postxt = poschip(ichip,:) + [7000,7200];
posstr = convertStringsToChars(num2str(ichip));
gs(end+1) = gdsii_ptext(posstr, postxt, 100, 25, 0);
 
    
for imem = 1:size(posmem,1)
%for imem = 1:1


%%%% bottom transducers
Lambda = Lambdav(imem);
Npairs = round(200/Lambda);
[~,~] = idt(posmem(imem,:) + poschip(ichip,:) + [0,-200], Lambda,210,160,Npairs,1,31,32);

%%%% top transducers
Lambda = Lambdav(imem+6);
Npairs = round(200/Lambda);
[~,~] = idt(posmem(imem,:) + poschip(ichip,:) + [0, 200], Lambda,210,160,Npairs,2,31,32);


%%%% PhCs
fx = f; 
fy = f;

PhC_rect( posmem(imem,:) + poschip(ichip,:),...
         a,fx,fy,PhCsize,20)


end

end 






%-------------------------------------------------------------------------%
% create a library to hold the structure
glib = gds_library('all', 'uunit',1e-6, 'dbunit',1e-9, gs);
% finally write the library to a file
write_gds_library(glib, '!DiffrOptom_v6_IDT_tmp.gds');

disp('End');