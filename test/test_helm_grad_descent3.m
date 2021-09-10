if (1>0)
addpath('../src');
nangs = 40;
[verts,xyin,angs,targs] = init_shape(nangs);

nvs = 4;
r0  = 2;
[vs] = init_guess(nvs,xyin,r0);

nvs_out = nvs;
zk = 1.1;

chnkrs = [];
err_curr = 10000;
end

while (zk<10)
disp(zk)

tic
[u,chnkr,bd_sol,F,err1] = helm_dirichlet_solver(verts,zk,targs,angs,xyin);
toc
return;

xy_s = sum(vs')'/nvs;
[u_s,chnkr_s,bd_sol_s,F_s,err_s] = helm_dirichlet_solver(vs,zk,targs,angs,xy_s);
chnkrs = [chnkrs, chnkr_s];


err_curr = norm(u-u_s,'fro');
return

vs_in = vs;
[vs_out,e_new] = opt_sing_freq(vs_in,nvs,zk,u,targs,angs,...
         nangs);
%[vs_out,ier,e_new] = opt_sing_freq_min(vs_in,nvs,zk,u,targs,angs,...
%         nangs,rcut);
%if (ier == 0)
    vs = vs_out;
%    if (e_new < err_curr)
%        vs = vs_out;
        err_curr = e_new;
%    end    
%end     
%vs = vs_out;

iffixed = false;
while (iffixed == false)
[vs_out,nvs_out] = add_verts(vs,nvs,zk,u,targs,angs,nangs);

vs_in = vs_out;
nvs_in= nvs_out;
rcut  = 0;
[vs_out,ier,e_new] = opt_sing_freq_min(vs_in,nvs_in,zk,u,...
    targs,angs,nangs,rcut);
if (ier == 0)
    vs = vs_out;
    if (nvs_out == nvs)
        iffixed = true;
    end
    nvs = nvs_out;
else
    if (e_new < err_curr)
        vs = vs_out;
        err_curr = e_new;
        if (nvs_out == nvs)
            iffixed = true;
        end
        nvs = nvs_out;
    end    
end
end

zk = zk + 0.25;
end