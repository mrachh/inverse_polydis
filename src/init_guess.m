function [verts] = init_guess(nverts,xyin,r0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    pi = atan(1)*4;
    angs = (0:(nverts-1))/nverts*2*pi;
    xs = cos(angs)*r0+xyin(1);
    ys = sin(angs)*r0+xyin(2);
    verts = [xs;ys];
end

