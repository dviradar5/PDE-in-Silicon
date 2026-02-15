function [J1] = Gauss_Src(osc,src_array,grid3d,obj_array)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
ax=grid3d.bound;
xc=linspace(ax(1),ax(4),grid3d.N(1));
yc=linspace(ax(2),ax(5),grid3d.N(2));
zc=linspace(ax(3),ax(6),grid3d.N(3));
[X,Y,Z] = meshgrid(xc,zc,yc);
[~,r,z] = cart2pol(X,Y,Z);


%%%%%%%%%change 2025 
% disp(obj_array) % Check what obj_array contains
% if ~isstruct(obj_array)
%     error('obj_array is not a structure array.');
% end
% 
% if ~isfield(obj_array, 'material')
%     error('obj_array does not have a "material" field.');
% end
% 
% if ~isfield(obj_array(1,1).material, 'eps')
%     error('obj_array.material does not have an "eps" field.');
% end
%%%%%%%%%%%%%%%
% obj_array
% osc
n_index=sqrt(obj_array(1,1).material.eps(1,1));
lambda=osc.in_L0;        
k = 2*pi/lambda*n_index;       % wave number;

for a_axis=Axis.elems
    J1{a_axis}=zeros(grid3d.N);
    if src_array.polarization==a_axis
        w0=diff(src_array.shape.lprim{1});        
        zR = k*w0^2/2;  % Rayleigh distance;
        z0=src_array.l{2};
        z=z-z0;
        R=@(z) z+zR^2./z;
        w = w0*sqrt(1+(z/zR).^2);        
        ww=(z+r.^2./2./R(z));
        m_factor=z(1)+1;
        m_factor=[m_factor m_factor+lambda/n_index];
        G =@(z,r) 1./sqrt(1+z.^2/zR^2).*exp(-r.^2./w.^2).*exp(-0.5i*k*(z.*r.^2)./(z.^2+zR^2));
        %gg=(m_factor(1)<ww)-(m_factor(2)<ww);
        gg=1-(m_factor(2)<ww);
        ggg=abs(gg.*G(z,r));
        J1{a_axis}=squeeze(ggg);
%         figure;imagesc(squeeze(gg)); axis image
%         figure;imagesc(real(squeeze(G(z,r)))); axis image
    end
end
end

