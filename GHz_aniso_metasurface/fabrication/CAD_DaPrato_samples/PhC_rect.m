function PhC_rect(pos,a,fx,fy,PhCsize,layer)

global gs

        nix = floor(PhCsize/a);
        niy = floor(PhCsize/a); 
        holec = cell(nix*niy,1);
        k = 1;
        for ix = 1:nix
        for iy = 1:niy
        k

        holec{k} = a * [0,0; fx,0; fx,fy; 0,fy; 0,0] + a*[ix-1,iy-1] ...
            - [1,1]*PhCsize/2  + pos;
        k = k+1;
        end
        end
        
    gs(end+1) = gds_element('boundary', 'xy',holec, 'layer',layer);