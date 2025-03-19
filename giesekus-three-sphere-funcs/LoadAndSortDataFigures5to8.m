function LoadAndSortDataFigures5to8()

clear all; close all; clc;

%% setup

% sims to load
ld_sims = 1:5;

%% load data

for Ii=1:625

    clear Xt_all;

    for Jj = ld_sims

        load(['data/set_',num2str(Jj),'_sim_',num2str(Ii-1),'.mat'],'Xt',...
            'ft','velX');
        Xt_all{Jj}   = Xt;
        ft_all{Jj}   = ft;
        velX_all{Jj} = velX;

        % calculate Ubar(s)
        if Jj == 1
            X0 = Xt(1,2,1);
        end
        Ubar(Ii,Jj) = Xt(1,2,end-1) - X0;

    end

    [Xt_s1{Ii},Xt_s2{Ii},Xt_s3{Ii},ft_s1{Ii},ft_s2{Ii},ft_s3{Ii},...
        velX_s1{Ii},velX_s2{Ii},velX_s3{Ii}] = ...
        SortData(Xt_all,ft_all,velX_all,ld_sims);
    fprintf(['Loaded data for Ii = ',num2str(Ii-1),'\n'])
end

Xt_s1 = cell2mat(Xt_s1);
Xt_s2 = cell2mat(Xt_s2);
Xt_s3 = cell2mat(Xt_s3);

ft_s1 = cell2mat(ft_s1);
ft_s2 = cell2mat(ft_s2);
ft_s3 = cell2mat(ft_s3);

velX_s1 = cell2mat(velX_s1);
velX_s2 = cell2mat(velX_s2);
velX_s3 = cell2mat(velX_s3);

for Ii=1:625

    W = velX_s1(:,Ii)'*ft_s1(:,Ii) + velX_s2(:,Ii)'*ft_s2(:,Ii) + ...
        velX_s3(:,Ii)'*ft_s3(:,Ii);

    Wbar(Ii,1) = W/size(ft_s1,1);
    eff(Ii,1)  = Ubar(Ii,ld_sims(end))^2/Wbar(Ii,1);

end

%% save

save('data/All_Sims_Loaded_Data.mat');

end

