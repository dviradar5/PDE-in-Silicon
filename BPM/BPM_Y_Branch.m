% Beam Propagation Method 30 Marzo 2007
% Edgar Guevara Codina 
% Dispositivos Optoelectronicos 
% Propagacion de un pulso gaussiano por un acoplador Y
close all; clear all;       % cierra ventanas y borra variables
%----------- Declaracion de Variables ------------------------------- 
xa = -100e-6;               % Coordenada Inicial
xb = 100e-6;                % Coordenada Final
num_samples = 1000;         % Numero de muestras (par)
H = xb - xa;                % Tamaño del dominio espacial
del = H/num_samples;        % Espaciado de las muestras
x = linspace (xa, xb-del, num_samples);     % Dominio espacial
A = 1;                      % Amplitud del pulso
x0 = 0.1e-6;                % Desplazamiento
W0 = 0.8e-6;                % Radio de la cintura del pulso
nmax = 3.492;               % Indice de refraccion maximo 3.492 SiGe
nmin = 3.481;               % Indice de refraccion minimo 3.481 Si
lambda = 0.633e-6;          % Longitud de onda
k0=(2*pi/lambda);           % Numero de onda
deltaz = 1e-6;              % Intervalo de Propagacion
zmax = 900;                 % Distancia Maxima de Propagacion (um)
% ---------- Reticula para graficar pulso propagado -------------------
[xx,yy] = meshgrid ([xa:del:xb-del],[1:1:zmax]);
% ---------- Generacion del pulso gaussiano ---------------------------
modo = A*exp (-((x+x0)/W0).^2);  % Pulso Gaussiano
dftmodo = fix(fft(modo));   % Transformada Discreta de Fourier de dicho pulso
%figure(1);
% Grafica el pulso original
%subplot(3,1,1); plot(x,modo,'r.');grid on; xlim([xa xb]);
%title(sprintf('Pulso Gaussiano Original %de^{-(x+%g/%g)^2}',A,x0,W0));
% Grafica la magnitud del pulso original
%subplot(3,1,2); plot(abs(dftmodo),'g.');grid on; title('Magnitud');
% Grafica la fase del pulso original
%subplot(3,1,3); plot(angle(dftmodo),'b.');grid on; title('Fase');
% --------------- Obtencion del indice promedio -----------------------
nbar = (nmax + nmin)/2;
% ---------- Generacion del perfil de la guia de onda -----------------
zz = imread('ybranch.bmp','BMP');       %Cargar Imagen con el perfil
zz=+zz;                                 %Convertir a numeros reales
zz=zz.';                                %Transpuesta
%figure(2);
for j = 1:1:num_samples,                %Cambiar indices de refraccion
    for k = 1:1:zmax,
        if(zz(k,j)==1)
            n(k,j)=nmin;
        else
            n(k,j)=nmax;
        end
    end
end
%Variacion del indice para un DPS (nbar+1)
%for j = 600:1:900,
    %for k = 500:1:520,
        %if(n(j,k)==nmax)
            %n(j,k)=nbar+1;
            %end
            %end
            %end
%surf(xx,yy,n);             %Grafica del perfil de la guia de onda
%view (0,90);
%shading interp;
%xlim([xa xb]);
%plot(x,n,'b.'); title('Perfil de indice de refraccion'); xlim([xa xb]);
% ---------------- Calculo de la fase 1 normalizada -------------------
kx1 = linspace(0,num_samples/2 + 1,num_samples/2);
kx2 = linspace(-num_samples/2,-1,num_samples/2);
kx = (2*pi/H)*[kx1 kx2];
%figure(3);
%plot(kx,'r.'); title('Vector de onda Kx');
fase1 = exp((i*deltaz*kx.^2)./(nbar*k0 + sqrt(max(0,nbar^2*k0*2 - kx.^2))));
% --------------- Calculo de la atenuacion -----------------
for j = 1:1:num_samples,
    if (j >= 370) && (j <= 630)
        od(j) = 0;
    else
        od(j) = 7500;
    end
end
%figure(4);
%plot(x,od,'g.'); title('Atenuacion de la señal'); xlim([xa xb]);

% --------------- Ciclo principal del programa ------------------------
%figure(5);
%hold on;
for k = 1:1:zmax,
        %Calcular el desfase contribuido por el efecto de lente
        fase2 = exp(-(od + i*(n(k,:) - nbar)*k0)*deltaz);
        % Aplicamos la Transformada Inversa
        modo = ifft((fft(modo).*fase1)).*fase2;
        zz(k,:) = abs(modo);
 end
% Graficamos la magnitud de nuestros pulsos propagados
surf(xx,yy,zz);
colormap hsv;
shading interp;
grid on;
view (0,90);
xlim([xa/4 xb/4]);
title('Acoplador Y');
%figure(6);
%contour(xx,yy,zz,30);
%xlim([xa/4 xb/4]);
%title('Grafica de Contorno del acoplador Y');