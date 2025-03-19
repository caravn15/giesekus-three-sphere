function [Xt_s1,Xt_s2,Xt_s3,ft_s1,ft_s2,ft_s3,velX_s1,velX_s2,velX_s3] = ...
    SortData(Xt_all,ft_all,velX_all,ld_sims)

numt = 2501;

for Jj = ld_sims

    Xt   = Xt_all{Jj};
    ft   = ft_all{Jj};
    velX = velX_all{Jj};

    Xt_s1((Jj-1)*numt+2-Jj:Jj*numt-Jj+1,1) = squeeze(Xt(1,1,1:end-1));
    Xt_s2((Jj-1)*numt+2-Jj:Jj*numt-Jj+1,1) = squeeze(Xt(1,2,1:end-1));
    Xt_s3((Jj-1)*numt+2-Jj:Jj*numt-Jj+1,1) = squeeze(Xt(1,3,1:end-1));

    ft_s1((Jj-1)*numt+2-Jj:Jj*numt-Jj+1,1) = squeeze(ft(1,1,1:end));
    ft_s2((Jj-1)*numt+2-Jj:Jj*numt-Jj+1,1) = squeeze(ft(1,2,1:end));
    ft_s3((Jj-1)*numt+2-Jj:Jj*numt-Jj+1,1) = squeeze(ft(1,3,1:end));

    velX_s1((Jj-1)*numt+2-Jj:Jj*numt-Jj+1,1) = squeeze(velX(1,1,1:end));
    velX_s2((Jj-1)*numt+2-Jj:Jj*numt-Jj+1,1) = squeeze(velX(1,2,1:end));
    velX_s3((Jj-1)*numt+2-Jj:Jj*numt-Jj+1,1) = squeeze(velX(1,3,1:end));

end

end