addpath('../src');
verts = [0,1,1,0;0,0,1,1];
xyin = [0.3;0.4];
zk = 1.1;
angs = [0;pi/2];
targs = [3.1, 4.2, 4; -1.7,2.5,-3.6];
[u,chnkr,bd_sol,F,err1] = helm_dirichlet_solver(verts,zk,targs,angs,xyin);