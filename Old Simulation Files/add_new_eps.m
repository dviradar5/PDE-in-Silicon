function [eps_new,Jnew,FCC] = add_new_eps(eps,osc,src_array,grid3d,obj_array,J)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
load('doping_parameters','SP','I','old_eps','FCC');

if sum(I(:))>0

    %% Calculated absorpt photon  

    alpha_z=sqrt((abs(old_eps)-real(old_eps))/2)*4*pi/(SP.lambda)*SP.dl(2)*SP.um;


    xy2=size(I);
    dif=round((size(eps{1})-size(I))/2);


    % special clean -(the unncessary indices)
    for jjj=1:3
        eps{jjj}(:,1)=eps{jjj}(:,2);
        ind0=find(mean(eps{jjj})>mean(eps{jjj}(:,1)));
        ind0=ind0(1)+1;
        eps{jjj}(:,ind0)=eps{jjj}(:,ind0+1); %this line
    end
    alpha_z=alpha_z(dif(1)+(1:xy2(1)),ind0-2+(1:xy2(2)));
    %figure; imagesc(alpha_z); title('alpha_z')

    Inphoton=I*SP.lambda*SP.E/SP.h_bar/SP.c0/(SP.dl(1)*SP.um*100)^3/SP.eV;


    nphotons=Inphoton.*(1-exp(-alpha_z));
    if sum(FCC(:))==0
        c_e=nphotons;
        c_h=nphotons;
    else
        c_e=nphotons+FCC;
        c_h=nphotons+FCC;
    end
    FCC=c_e;


    eps_new=eps;
    for jjj=1:3
        eps_new{jjj}(dif(1)+(1:xy2(1)),ind0-2+(1:xy2(2)))=SP.PDE_eps(c_e,c_h,eps{jjj}(dif(1)+(1:xy2(1)),ind0-2+(1:xy2(2))));
    end
else
    eps_new=eps;
    FCC=zeros(10);
end
if SP.gauss_src
%     [Jnew] = Gauss_Src(SP.lambda,src_array,grid3d,SP.um,SP.waist,SP.z0);
[Jnew] = Gauss_Src(SP.lambda,src_array,grid3d,SP.waist);%change 2025
else
    Jnew=J;
end

end

