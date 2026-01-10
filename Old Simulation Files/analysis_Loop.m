clc;clear;close all
exp_name='N1';
dir_exp=['results\',exp_name];
m=dir([dir_exp,'_LN*.mat']);
zz_loc=15;
Ir_loop=0;
I_loop=0;
marks_sets={'--','-','-.',':'};
figure(1); subplot(121); hold all; legend;
subplot(122); hold all; legend;
figure(2);
set_num=round(linspace(1,20,5));
count=0;
for ii=1:length(m)
    load([dir_exp,'_LN',num2str(ii),'.mat'])
    ex_sim=['Energy = ',num2str(SP.E*ii*1e9),'nJ'];
    xx=linspace(SP.sizes(1),SP.sizes(2),size(I,1));
    yy=linspace(SP.sizes(5),SP.sizes(6),size(I,2));
    [Z,X]=meshgrid(yy,xx);
    [~,ind]=min(abs(zz_loc-yy));
    Ir=create_radial_image(I(:,ind));
    figure; subplot(3,2,1);
    imagesc(xx,xx,Ir); axis image
    title(['Profile, z = ',num2str(zz_loc)])
    subplot(3,2,2);
    I_loop=I_loop+I;
    Ir_loop=Ir_loop+Ir;
    imagesc(xx,xx,Ir_loop); axis image
    title(['overall at z =',num2str(zz_loc)])
    subplot(3,2,3);
    %alpha_z=imag(SP.n_index)*4*pi/SP.lambda;
    alpha_z=0.0086*4*pi/SP.lambda;
    %imagesc(yy,xx,I.*exp(alpha_z*SP.um*Z)); title('I');
    imagesc(yy,xx,I); title('I');
    subplot(3,2,4);
    imagesc(yy,xx,I_loop); title('I');
    colorbar;
    axis image
    line(yy([ind, ind]),xx([1, end]),'Color','red','LineStyle','--')
    subplot(3,2,[5,6]);
    imagesc(yy,xx,FCC);title('FCC');
    colorbar;
    axis image
    sgtitle(ex_sim);
    figure(1);
    if any(ii==set_num)
        subplot(121);
        plot(xx,Ir(round(length(xx)/2),:),marks_sets{1+mod(ii,length(marks_sets))},'DisplayName',ex_sim);
        xlabel('x [\mum]'); ylabel('Intensity [a.u.]')
        subplot(122);
        plot(xx,Ir_loop(round(length(xx)/2),:),'DisplayName',ex_sim);
        xlabel('x [\mum]'); ylabel('Intensity [a.u.]')
        sgtitle(['z = ',num2str(zz_loc)])
    end
    figure(2);
    if any(ii==set_num)
        count=count+1;
        subplot(2,length(set_num),count)
        
        imagesc(xx,xx,Ir); axis image
        title(ex_sim,'Probe'); xlabel('x [\mum]');ylabel('y [\mum]');
        subplot(2,length(set_num),count+length(set_num))
        imagesc(xx,xx,Ir_loop); axis image
        title(ex_sim,'Pump'); xlabel('x [\mum]');ylabel('y [\mum]');
        sgtitle(['z = ',num2str(zz_loc),'\mum'])
    end

end
disp(sum(FCC(:)))
load([dir_exp,'_single_LN',num2str(2),'.mat'])
disp(sum(FCC(:)))