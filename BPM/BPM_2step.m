% Beam Propagation Method 13 Marzo 2007
% Edgar Guevara Codina 
% Dispositivos Optoelectronicos 
% Propagacion de un pulso gaussiano por dos guias de onda con perfil rectangular
close all; clear all;       % cierra ventanas y borra variables
%----------- Declaracion de Variables ------------------------------- 
xa = -100e-6;               % Coordenada Inicial
xb = 100e-6;                % Coordenada Final
num_samples = 1000;         % Numero de muestras (par)
H = xb - xa;                % Tamaño del dominio espacial
del = H/num_samples;        % Espaciado de las muestras
x = linspace (xa, xb-del, num_samples);     % Dominio espacial
A = 1;                      % Amplitud del pulso
x0 = 7e-6;                  % Desplazamiento
W0 = 3e-6;                  % Radio de la cintura del pulso
nmax = 1.5;                 % Indice de refraccion maximo
nmin = 1.499;               % Indice de refraccion minimo
x1 = -12e-6;                % lado izquierdo del 1er escalon
x2 = -2e-6;                 % lado derecho del 1er escalon
x3 = 2e-6;                  % lado izquierdo del 1er escalon
x4 = 12e-6;                 % lado derecho del 1er escalon
%m1 = (nmax - nmin)/(-xmin);% pendiente positiva
%m2 = (nmin - nmax)/(xmax); % pendiente negativa
lambda = 1e-6;              % Longitud de onda
k0=(2*pi/lambda);           % Numero de onda
deltaz = 1e-6;              % Intervalo de Propagacion
zmax = 10000;               % Distancia Maxima de Propagacion (um)
% ---------- Reticula para graficar pulso propagado -------------------
[xx,yy] = meshgrid ([xa:del:xb-del],[0:10:zmax]);
% ---------- Generacion del pulso gaussiano ---------------------------
modo = A*exp (-((x+x0)/W0).^2);  % Pulso Gaussiano
dftmodo = fix(fft(modo));   % Transformada Discreta de Fourier de dicho pulso
figure(1);
% Grafica el pulso original
subplot(3,1,1); plot(x,modo,'r.');grid on; xlim([xa xb]);
title(sprintf('Pulso Gaussiano Original %de^{-(x+%g/%g)^2}',A,x0,W0));
% Grafica la magnitud del pulso original
subplot(3,1,2); plot(abs(dftmodo),'g.');grid on; title('Magnitud');
% Grafica la fase del pulso original
subplot(3,1,3); plot(angle(dftmodo),'b.');grid on; title('Fase');
% --------------- Obtencion del indice promedio -----------------------
nbar = (nmax + nmin)/2;
% ---------- Generacion del perfil de la guia de onda -----------------
for j = 1:1:num_samples
    if (x(j) >= x1) & (x(j) <= x2)
        n(j) = nmax;
    elseif (x(j) >= x3) & (x(j) <= x4)    
        n(j) = nmax;
    else
        n(j) = nmin;
    end
end
figure(2);
plot(x,n,'b.'); title('Perfil de indice de refraccion'); xlim([xa xb]);
% ---------------- Calculo de la fase 1 normalizada -------------------
kx1 = linspace(0,num_samples/2 + 1,num_samples/2);
kx2 = linspace(-num_samples/2,-1,num_samples/2);
kx = (2*pi/H)*[kx1 kx2];
figure(3);
plot(kx,'r.'); title('Vector de onda Kx');
fase1 = exp((i*deltaz*kx.^2)./(nbar*k0 + sqrt(max(0,nbar^2*k0*2 - kx.^2))));
% --------------- Calculo de la fase 2 con atenuacion -----------------
for j = 1:1:num_samples
    if (j >= 400) && (j <= 600)
        od(j) = 0;
    else
        od(j) = 3500;
    end
end
figure(4);
plot(x,od,'g.'); title('Atenuacion de la señal'); xlim([xa xb]);
fase2 = exp(-(od + i*(n - nbar)*k0)*deltaz);
% --------------- Ciclo principal del programa ------------------------
figure(5);
hold on;
for k = 0:1:zmax/10
    for j = 1:1:10
        % Aplicamos la Transformada Inversa
        modo = ifft((fft(modo).*fase1)).*fase2;
        zz(k+1,:) = abs(modo);
    end
end
% Graficamos la magnitud de nuestros pulsos propagados
surf(xx,yy,zz,'FaceColor','interp','edgecolor','none','FaceLighting','phong');
colormap jet;
camlight right;
grid on;
view (30,45);
xlim([xa xb]);
title('Magnitud del Pulso Gaussiano Propagado');
figure(6);
contour(xx,yy,zz,30);
title('Grafica de Contorno');