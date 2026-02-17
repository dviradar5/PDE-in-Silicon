% Beam Propagation Method 20 Febrero 2007
% Edgar Guevara Codina 
% Dispositivos Optoelectronicos 
% Genera un pulso gaussiano y lo propaga de 0 a 1000 um
% Grafica superficie
close;clear all;
% ---------- Generacion -----------------
x1 = -60;                           % Coordenada Inicial
x2 = 60;                            % Coordenada Final
num_samples = 240;                  % Numero de muestras (par)
H = x2-x1;                          % Tamaño del dominio espacial
del = H/num_samples;                        % Espaciado de las muestras
x = linspace (x1, x2-del, num_samples);     % Dominio espacial
A = 1;                              % Amplitud del pulso
W0 = 8;                             % Radio de la cintura del pulso
modo = A*exp (-(x/W0).^2);          % Pulso Gaussiano
dftmodo = fix(fft(modo));           % Transformada Discreta de Fourier de dicho pulso
figure(1);
subplot(3,1,1); plot(x,modo,'rx');grid on; title('Pulso Gaussiano Original Ae^{-(x/8)^2}');
% Grafica el pulso original
subplot(3,1,2); plot(abs(dftmodo),'gx');grid on; title('Magnitud');
% Grafica la magnitud del pulso original
subplot(3,1,3); plot(angle(dftmodo),'bx');grid on; title('Fase');
% Grafica la fase del pulso original
% ---------- Propagacion -----------------
s1 = 1:num_samples/2;               % Region 1 del desfase
s2 = num_samples/2+1:num_samples;   % Region 2 del desfase
s = [s1 s2];                        % Region total del desfase
n = 1;                              % Indice de refraccion del medio
lambda = 0.8;                       % Longitud de onda
%L=[0 500 1000];                    % Distancias de Propagacion
k0=n*(2*pi/(n*lambda));             % Numero de onda
z=linspace(0,0,num_samples);        % Desplazamiento en z
figure(2);
hold off;
% Generamos la reticula para graficar
[x1,y1] = meshgrid ([x1:del:x2-del],[0:10:990]); 
for L = 0:10:990
    desfase1 = n*k0*L*(1.*(1-0.5*(lambda^2*(s1-1).^2/H^2)));
    desfase2 = n*k0*L*(1.*(1-0.5*(lambda^2*(num_samples+1-s2).^2/H^2)));
    desfase =[desfase1 desfase2];
    nuevadft = abs(dftmodo).*exp(i*(angle(dftmodo)+(desfase)));
    % Aplicamos la Transformada Inversa
    nuevomodo=ifft(nuevadft);
    % Generamos la matriz para graficar
    z1(L/10+1,:) = abs(nuevomodo);
end
% Graficamos la magnitud de nuestros pulsos propagados
surf(x1,y1,z1);
shading interp;
colormap jet;
grid on;
title('Magnitud de los Pulsos Gaussianos Propagados');