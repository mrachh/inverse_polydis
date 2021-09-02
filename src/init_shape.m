function [verts,xyin,angs,targs] = init_shape(nangs)

    verts = [14,14,22,21,27,25,27,23,22,18,20,17,14,11,8,10,6,5,1,3,1,7,6,13,13];
    verts = (verts-13.5)/6;
    verts2= [1,8,7,10,15,16,21,20,22,19,27,25,31,25,27,19,22,20,21,16,15,10,7,8,1];
    verts2= (verts2-15)/6;
    verts = [verts;verts2];
pi = atan(1)*4;    
angs  = (0:(nangs-1))/nangs*2*pi;
targs = 10*[cos(angs);sin(angs)];
angs  = angs';
xyin  = [0;0];
end

