% Limpiar el espacio de trabajo y la ventana de comandos
clearvars;
clear;
clc;
clf;

% Parámetros iniciales
nl = 1;            % Número de loops de cables
ds = 0.1;          % Diferencial de longitud
x = -5:ds:5;       % Vectores de -5 a 5 con un incremento de ds
y = -5:ds:5;       % Definen nuestro espacio de trabajo en 3D
z = x;
Lx = length(x);    % Tamaño de los vectores
Ly = length(y);
Lz = length(z);
rw = 0.2;          % Espesor del cable
I = 300;           % Corriente eléctrica
mo = 4*pi*1e-7;    % Permeabilidad del espacio libre [T*m/A]
km = mo * I / (4 * pi);  % Constante magnética [T]

N = 15;            % Número de puntos por loop
R = 1.5;           % Radio del cable
sz = 1;            % Loop de paso en la dirección del eje z
s = 1;             % Número del loop
dtheta = 2 * pi / N;  % Paso de ángulo diferencial [radianes]
dl = R * dtheta;   % Diferencial de longitud
ang = 0:dtheta:2*pi-dtheta;  % Vector de ángulos

% Cálculo de las coordenadas y componentes direccionales
for i = 1:nl
    Px(s:s+N-1) = R * cos(ang);  % Coordenadas X del loop circular
    Py(s:s+N-1) = R * sin(ang);  % Coordenadas Y del loop circular
    Pz(s:s+N-1) = -nl/2 * sz + (i-1) * sz;  % Coordenadas Z del loop circular
    dx(s:s+N-1) = -Py(s:s+N-1) * dtheta;  % Componentes X de los vectores direccionales
    dy(s:s+N-1) = Px(s:s+N-1) * dtheta;   % Componentes Y de los vectores direccionales
    s = s + N;  % Actualizar el índice para el siguiente loop
end

dz(1:N*nl) = 0;  % Componentes Z de los vectores direccionales en 0

% Visualización del cable de corrientes usando quiver3
quiver3(Px, Py, Pz, dx, dy, dz, 0.5, '-r', 'LineWidth', 2);
view(-34, 33);
figure(1);

% Inicialización de matrices para el campo magnético
dBx = zeros(Lx, Ly, Lz);  % Campo magnético en X
dBy = dBx;                % Campo magnético en Y
dBz = dBx;                % Campo magnético en Z

% Cálculo del campo magnético usando la ley de Biot-Savart
for I = 1:Lx
    for J = 1:Ly
        for K = 1:Lz
            for L = 1:nl*N
                rx = x(I) - Px(L);
                ry = y(J) - Py(L);
                rz = z(K) - Pz(L);
                r = sqrt(rx^2 + ry^2 + rz^2 + rw^2);
                r3 = r^3;

                dBx(I, J, K) = dBx(I, J, K) + km * dy(L) * rz / r3;
                dBy(I, J, K) = dBy(I, J, K) + km * dx(L) * rz / r3;
                dBz(I, J, K) = dBz(I, J, K) + km * (dx(L) * ry - dy(L) * rx) / r3;
            end
        end
    end
end

% Magnitud del campo magnético
Bmag = sqrt(dBx.^2 + dBy.^2 + dBz.^2);
centery = round(Ly / 2);
Bx_xz = squeeze(dBx(:, centery, :));
Bz_xz = squeeze(dBz(:, centery, :));
Bxz = squeeze(Bmag(:, centery, :));

% Visualización del campo magnético
figure(2);
hold on;
pcolor(x, z, (Bxz').^(1/3)); shading interp; colormap jet; colorbar;
h1 = streamslice(x, z, Bx_xz', Bz_xz', 3);
set(h1, 'Color', [0.8 1 0.9]);
xlabel('x');
ylabel('z');
title('Campo magnético de una corriente circular');
axis([-5 5, -5 5]);

% Parámetros del dipolo magnético
mag = 2000;  % Momento magnético [Am^2]
m = 0.004;   % Masa del imán [kg]
g = 9.81;    % Gravedad
w = m * -g;  % Peso del imán
zo = 5;      % Posición inicial del dipolo en z
dt = 0.05;   % Paso de tiempo

zm(1) = zo;      % Vector con las posiciones del imán
zmfree(1) = zo;  % Vector con posición del imán en caída libre
tt(1) = 0;       % Vector de tiempo
vz(1) = 0;       % Vector de la velocidad en z
vzfree(1) = 0;   % Vector de la velocidad en z en caída libre
cc = 1;          % Contador

% Animación de la caída del dipolo magnético
path = animatedline;
while zm(cc) > -5 % El bucle se ejecuta mientras la posición del dipolo sea mayor que -5

    addpoints(path, 0, zm(cc)); % Agrega puntos a la animación del camino del dipolo.
    drawnow; % Actualiza la figura para mostrar los cambios.
    head = scatter(0, zm(cc), 100, 'filled'); % Dibuja el dipolo en la posición actual en la figura.

    % Calcula la fuerza magnética sobre el dipolo en función de su posición actual zm((cc))
    Fm(cc) = (6 * mo * I * R^2 * mag * zm(cc)) / (4 * (zm(cc)^2 + R^2)^(5/2));
    F(cc) = Fm(cc) + w; % Calcula la fuerza total F sobre el dipolo sumando la fuerza magnética Fm y la fuerza gravitacional w.
    a = F(cc) / m; % Calcula la aceleración a del dipolo
    pause(0.01);
    zm(cc+1) = zm(cc) + vz(cc) * dt + 0.5 * a * dt^2; % Actualiza la posición del dipolo
    zmfree(cc+1) = zmfree(cc) + vzfree(cc) * dt - 0.5 * g * dt^2;
    vz(cc+1) = (zm(cc+1) - zm(cc)) / dt; % Calcula la velocidad del dipolo en la siguiente iteración.
    vzfree(cc+1) = (zmfree(cc+1) - zmfree(cc)) / dt;
    cc = cc + 1;
    delete(head); % Elimina el dipolo dibujado en la posición anterior para actualizar su posición en la animación.

end

% Gráficas de resultados
figure(4);
subplot(1, 2, 1);
hold on;
plot(zm(1:length(Fm)), 1000 * Fm, '-b', 'LineWidth', 2);
plot([0, 0], [-150, 150], '-.k', 'LineWidth', 2);
grid on;
xlabel('Posición z (m)');
ylabel('Fuerza magnética (mN)');
title('Fuerza magnética de un anillo de corriente sobre un imán que cae');
legend('Fuerza magnética en la dirección Z', 'Ubicación del loop actual', 'Location', 'southeast');

subplot(1, 2, 2);
hold on;
tt = 0:dt:(cc-1) * dt;
plot(tt, zm, '-r', 'LineWidth', 2);
plot(tt, zmfree, '--b', 'LineWidth', 2);
plot([0, 1.8], [0, 0], '-.k', 'LineWidth', 2);
grid on;
xlabel('Tiempo (s)');
ylabel('Posición z (m)');
title('Posición vs tiempo de un dipolo magnético cayendo a través de un anillo de corriente');
legend('Caída sobre un anillo de corriente', 'Caída libre (sin fuerza magnética)', 'Ubicación del loop actual', 'Location', 'southwest');
axis([0 1.8 -6 6]);

% Parámetros para el siguiente cálculo
clear;
clc;
clf;

mag = 500;        % Momento magnético
Rring = 0.5;      % Radio del anillo [m]
zo = 0.1;         % Posición inicial del imán
zring = 0;        % Posición del anillo en z
dt = 0.01;        % Paso de tiempo
t(1) = 0;         % Vector de tiempo
zm(1) = zo;       % Vector de posición en z
cc = 1;           % Contador
vz(1) = 0;        % Vector de velocidad
mo = 4*pi*1e-7;    % Permeabilidad del espacio libre [T*m/A]
figure(1);


while zm(cc) > 0.0162
    pause(0.001) % Pausa en cada iteración
    clf % Limpia las figuras del workspace
    
    % Calculamos el campo magnético y flujo magnético en la posición actual
    [x, y, phiB1, Bz] = B_due_M(zm(cc), mag, Rring);
    
    % Actualizamos la posición del imán en su caída libre
    zm(cc+1) = zm(cc) + vz(cc) * dt - 0.5 * 9.81 * dt^2;
    
    % Actualizamos la velocidad del imán en su caída libre
    vz(cc+1) = (zm(cc+1) - zm(cc)) / dt;
    
    % Calculamos el campo magnético y flujo magnético en la nueva posición
    [x, y, phiB2, Bz] = B_due_M(zm(cc+1), mag, Rring);
    
    % Calculamos la fuerza electromotriz
    fem(cc) = (phiB2 - phiB1) / dt;
    
    % Gráfica 1: Fuerza electromotriz vs Tiempo
    subplot(2, 2, 1)
    hold on
    grid on
    xlabel('time, s')
    ylabel('fem, mV')
    plot(t(1:cc), 100 * fem(1:cc), '-k', 'LineWidth', 1)
    plot(t(1:cc), 100 * fem(1:cc), '*r', 'LineWidth', 2)
    
    % Gráfica 2: Altura del imán vs Tiempo
    subplot(2, 2, 2)
    hold on
    axis([0 0.3 -10 10])
    grid on
    xlabel('time, s')
    ylabel('magnet height, cm')
    plot(t(1:cc), 100 * zm(1:cc), 'ob', 'LineWidth', 2)
    
    % Gráfica 3: Campo magnético en el espacio (color)
    subplot(2, 2, 3)
    hold on
    pcolor(x, y, zm(cc) / abs(zm(cc)) * abs(abs(0.005^2 * Bz)).^(1/3)) % Aqui se escala para graficar mejor
    shading interp
    colormap hot
    colorbar
    view(-45, -45)
    
    % Gráfica 4: Campo magnético en el plano del anillo (superficie)
    subplot(2, 2, 4)
    hold on
    mesh(x, y, zm(cc) / abs(zm(cc)) * abs(abs(10^2 * Bz)).^(1/3)) % Se escala para graficar mejor
    view(-30, -3)
    axis([-Rring Rring -Rring Rring -5 15])
    
    % Incrementamos el contador y el tiempo para la siguiente iteración
    cc = cc + 1;
    t(cc) = t(cc-1) + dt;
end

function [x, y, phiB, Bz] = B_due_M(z, mag, Rring)
    % Constante de permeabilidad del vacío
    mo = 4 * pi * 1e-7;
    
    % Incremento de longitud para calcular el espacio
    ds = 0.005;
    
    % Matriz de x que tiene valores del diámetro del anillo
    x = -Rring:ds:Rring;
    
    % Matriz de y que tiene valores del diámetro del anillo
    y = -Rring:ds:Rring;
    
    % Número de puntos en la matriz x
    Lx = length(x);
    
    % Número de puntos en la matriz y
    Ly = length(y);
    
    % Matriz vacía del mismo tamaño de nuestro espacio de trabajo para guardar
    % los valores del campo magnético
    Bz = zeros(Lx, Ly);
    
    % Valor inicial del flujo magnético
    phiB = 0;
    
    % Loop que recorre la matriz x
    for i = 1:Lx
        % Loop que recorre la matriz y
        for j = 1:Ly
            % Distancia del centro del anillo al punto en cuestión
            r = sqrt(x(i)^2 + y(j)^2);
            
            % Condicional que asegura que el punto que estamos calculando se
            % encuentra dentro del anillo
            if r < Rring
                % Cálculo del campo magnético en ese punto
                Bz(i, j) = mo / (4 * pi) * (3 * z * (mag * z) - mag * (x(i)^2 + y(j)^2 + z^2)) / ((x(i)^2 + y(j)^2 + z^2 + (ds / 10)^2)^(5/2));
                
                % Suma del valor del flujo magnético en ese punto al total
                phiB = phiB + ds^2 * Bz(i, j);
            end
        end
    end
end

