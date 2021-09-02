function [dp] = helm_dprime(zk,chnkr,F,sigma)
    [~,~,u,v] = lege.exps(chnkr.k);
    sigma_use = reshape(sigma,[chnkr.k,chnkr.nch]);
    sigma_coefs = complex(zeros(size(sigma_use)));
    sigma_d_coefs = complex(zeros(size(sigma_use)));
    sigma_d = complex(zeros(size(sigma_use)));
    for i=1:chnkr.nch
        sigma_coefs(:,i) = u*sigma_use(:,i);
        sigma_d_coefs(1:(chnkr.k-1),i) = lege.derpol(sigma_coefs(:,i));
        sigma_d(:,i) = v*sigma_d_coefs(:,i);
        dsdt = sqrt(chnkr.d(1,:,i).^2 + chnkr.d(2,:,i).^2).*chnkr.h(i);
        dsdt = dsdt';
        sigma_d(:,i) = sigma_d(:,i)./dsdt;
    end
    rnorms = normals(chnkr);
    dens2 = squeeze(rnorms(1,:,:));
    dens3 = squeeze(rnorms(2,:,:));
    dens1 = reshape(sigma_d,[chnkr.k*chnkr.nch,1]);
    dens2 = reshape(dens2,[chnkr.k*chnkr.nch,1]);
    dens3 = reshape(dens3,[chnkr.k*chnkr.nch,1]);
    s1 = rskelf_mv(F,dens1);
    s1 = reshape(s1,[chnkr.k,chnkr.nch]);
    s2 = rskelf_mv(F,dens2.*sigma);
    s3 = rskelf_mv(F,dens3.*sigma);
    
    for i=1:chnkr.nch
        sigma_coefs(:,i) = u*s1(:,i);
        sigma_d_coefs(1:(chnkr.k-1),i) = lege.derpol(sigma_coefs(:,i));
        sigma_d(:,i) = v*sigma_d_coefs(:,i);
        dsdt = sqrt(chnkr.d(1,:,i).^2 + chnkr.d(2,:,i).^2).*chnkr.h(i);
        dsdt = dsdt';
        sigma_d(:,i) = sigma_d(:,i)./dsdt;
    end
    sigma_d = reshape(sigma_d,[chnkr.k*chnkr.nch,1]);
    dp = sigma_d + zk^2.*(s2.*dens2 + s3.*dens3);
end