archivo = load("guitarras.mat");    % Carga del archivo de sonido

fs = archivo.fs;    % Obtención de la frecuencia de muestreo del archivo
ts = archivo.ts;    % Obtención del delta de tiempo del archivo

% Creación del vector tiempo en base a ts y los datos del archivo
t = [0:ts:size(archivo.tom)/fs - ts];   
x = archivo.tom;    % Obtención de los datos de la muestra

% Reproducción del sonido de los datos
% sound(archivo.tom, fs)

% Gráfico del audio como referencia
figure
plot(t, x)
title('Referencia del audio')
xlabel('Tiempo (s)')

App = max(x) - min(x);  % Amplitud peak-to-peak

% Obtención de los limites para la sumatoria de la potencia en base al vector de tiempo
limiteInf = find(t==1.36);
limiteSup = find(t==1.37);

% Cálculo de la potencia
potencia = sum(abs(x(limiteInf:limiteSup)).^2) / (limiteSup - limiteInf)

% Creación del vector de tiempo y datos de la primera cuerda 
% (desde que comienza hasta 10ms después)
tCuerda = t(limiteInf:limiteSup);   
xCuerda = x(limiteInf:limiteSup);

% Gráfico de la primera cuerda (desde que comienza hasta 10ms después)
figure
plot(tCuerda, xCuerda)
title('Pulsación de la cuerda')
xlabel('Tiempo (s)')

% Cálculo de la transformada de Fourier. 
% (Código obtenido del ejemplo 4.0 entregado por el profesor en canvas)
N = 1024; % tamaño de la Fast Fourier Transform (FFT)
X = fft(xCuerda, N); % FFT es habitualmente compleja, X(1) contiene componente 
% continua. N es el tamaño de la secuencia. Usar potencia de 2 para que
% el algoritmo sea mas eficiente

df = fs/N; % resolución de frecuencia en Hz para el grafico
Index = 0:N-1; % para el grafico de la FFT
f = Index*df; % indice del eje x convertido a frecuencia

Y = fftshift(abs(X));
sf = ([0:N-1]-round((N-1)/2))*df;

% Gráfico en escala logarítmica de la transformada de Fourier con amplitudes en dB.
figure
semilogx(sf, 20*log10(Y));
title('Transformada X(f) usando fftshift')
ylabel('Amplitud (dB)')
xlabel('Hz')

% Ecuación de sintesis de la muestra completa en base a amplitudes de frecuencias 
% representativas (frecuencia fundamental: 43.1 Hz)
xSintesisdB = -15.4 - 15.9*cos(2*pi*43.1*t);
xSintesisdB = xSintesisdB - 17.1*cos(2*pi*86.1*t);
xSintesisdB = xSintesisdB - 19.2*cos(2*pi*130.2*t);
xSintesisdB = xSintesisdB - 23.1*cos(2*pi*172.3*t); 
xSintesisdB = xSintesisdB - 28*cos(2*pi*215.3*t); 
xSintesisdB = xSintesisdB - 28*cos(2*pi*861*t);
xSintesisdB = xSintesisdB - 28*cos(2*pi*904.4*t);

% Gráfico de la ecuación de sintesis como referencia
figure
plot(t, xSintesisdB)
title('Ecuación de síntesis')
xlabel('Tiempo (s)')

% Reproducción del sonido para comparar con la muestra
% sound(xSintesisdB, fs)
