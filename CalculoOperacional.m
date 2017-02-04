%% Calculo operacional

% Trabajo realizado por los alumnos:
% Xiao Luo
% Jorge Calatayud Maeso
% Grupo: E4

%%  3.Circuito RC
clear

% Calculo mediante el comando dsolve:
sol1 = dsolve('1e3*1e-6*Dv + v = 10*heaviside(t)','v(0) = 0');

% Datos para el calculo mediante laplace:
syms t s Y;
Vg = 10*heaviside(t);
VG = laplace(Vg, t, s); 
v0 = 0;
Y1 = s*Y - v0;
R = 1e3;
C = 1e-6;

% Calculo mediante comandos laplace e ilaplace:
Sol2 = solve(R*C*Y1 + Y - VG, Y);
sol2 = ilaplace(Sol2, s, t);

% Representacion:
fprintf('3. Circuito RC\n')
fprintf('Solucion mediante dsolve: %s\n', char(sol1));
fprintf('Solucion mediante laplace e ilaplace: %s\n', char(sol2));

%%  4.Circuito RLC
clear

% Apartado a)

% Datos para laplace:
syms t s Y;
Vg = 10*heaviside(t);
VG = laplace(Vg, t, s);
v0 = 0;
dv0 = 0;
Y1 = s*Y - v0;
Y2 = s*Y1 - dv0;
R = 1e-1;
C = 1e-6;

% Caso 1:
L = 25e-8;

Sol1 = solve(L*C*Y2 + R*C*Y1 + Y - VG, Y);
sol1 = ilaplace(Sol1, s, t);
sol1 = vpa(sol1, 2);

% Caso 2:
L = 25e-10;

Sol2 = solve(L*C*Y2 + R*C*Y1 + Y - VG, Y);
sol2 = ilaplace(Sol2, s, t);
sol2 = vpa(sol2, 2);

% Caso 3:
L = 25e-11;

Sol3 = solve(L*C*Y2 + R*C*Y1 + Y - VG, Y);
sol3 = ilaplace(Sol3, s, t);
sol3 = vpa(sol3, 2);

% Representacion:
fprintf('\n4. Circuito LRC (Resueltos con laplace e ilaplace)\n');
fprintf('  a)\n');
fprintf('\tcaso 1:\n');
fprintf('\t\tSolucion; %s\n', char(sol1));
fprintf('\tcaso 2:\n');
fprintf('\t\tSolucion; %s\n', char(sol2));
fprintf('\tcaso 3:\n');
fprintf('\t\tSolucion; %s\n', char(sol3));

% Apartado b)
clear

% Datos para laplace:
syms t s Y;
Vg = 10*sin(2*pi*t*5e4);
VG = laplace(Vg, t, s);
v0 = 0;
dv0 = 0;
Y1 = s*Y - v0;
Y2 = s*Y1 - dv0;
R = 1e-1;
C = 1e-6;
L = 25e-11;

% Calculo:
Sol = solve(L*C*Y2 + R*C*Y1 + Y - VG, Y);
sol = ilaplace(Sol, s, t);
sol = vpa(sol,2);

% Representacion:
fprintf('  b)\n');
fprintf('\t\tSolucion; %s\n', char(sol));

%%  5.Otro circuito RLC
clear

% Datos:
syms t s Y1 Y2;
Vg = 10*heaviside(t);
VG = laplace(Vg, t, s);
R1 = 1e1;
R2 = 1e-1;
L = 1e-3;
C = 1e-3;
i10 = 10/(R1 + R2);
i20 = 10/(R1 + R2);
di10 = 10*R2/(L*(R1 + R2));
di20 = 10*R2/(L*(R1 + R2));
Y11 = s*Y1 - i10;
Y12 = s*Y11 - di10;
Y21 = s*Y2 - i20;
Y22 = s*Y21 - di20;

% Calculo:
[Sol1, Sol2] = solve(L*Y11 + R1*Y1 - L*Y21 - VG, L*C*Y12 + R2*C*Y21 + Y2 - L*C*Y22, Y1, Y2);
sol2 = ilaplace(Sol2, s, t);

v = int(sol2, 0, t)/C;
v = vpa(v,2);

% Representacion:
fprintf('\n5. Otro circuito RLC (Resuelto con laplace e ilaplace)\n');
fprintf('\tSolucion: v(t) = %s\n', char(v));

%%  6.Circuito RL
clear

% Datos:
syms t s Y;
Vg = 10*heaviside(t);
VG = laplace(Vg, t, s);
R = 1e2;
L = 1e-3;
i0 = 0;
Y1 = s*Y - i0;

% Apartado 6.1:
sol1 = dsolve('1e-5*DI + I = 1e-1*heaviside(t)','I(0) = 0');
sol1 = L*diff(sol1, 1);
sol1 = vpa(sol1,2);

% Apartado 6.2:
Sol2 = solve(L*Y1 + R*Y - VG, Y);
VL = L*s*Sol2;
vL = ilaplace(VL, s, t);

% Representacion:
fprintf('\n6. Circuito RL\n');
fprintf('Solucion mediante dsolve: %s\n', char(sol1));
fprintf('Solucion mediante laplace e ilaplace: %s\n', char(vL));

% Apartado 6.3:
% poniendo la ecuacion en su forma canonica, se ve que la variable tau
% tiene valor de L/R. La solucion homogenea tiene forma de A*exp(-t/tau).
% La solucion forzada sale de inmediato y es de valor V_G(t)/R. Teniendo en
% cuenta la condicion inicial de i(0) = 0, de la solucion general se saca
% la constante A, que es de valor -V_G(t)/R. Por tanto la solucion final
% es V_G(t)*(1-exp(-t/tau))/R. Sustituyendo los valores numericos en la
% ecuacion se ve que es exactamente lo mismo de lo que se obtiene en matlab


% Representacion grafica:
% Puesto que hacia aproximadamente 5*tau entra en regimen permanente, aqui
% se elige T = 10*tau para poder visualizarlo bien.
ezplot(vL, [0, 1e-4], 1);

%%  7.Circuito RLC2
clear

% Datos:
syms t s Y1 Y2;
R1 = 1e-1;
R2 = 1e3;
R3 = 1;
L = 25e-6;
C = 1e-6;
V = 10;
i10 = 0;
i20 = 0;
Y11 = s*Y1 - i10;
Y21 = s*Y2 - i20;

% Apartado 7.1:
Vg = V*heaviside(t);
VG = laplace(Vg, t, s);
[Sol1, Sol2] = solve(L*Y11 + (R1 + R2)*Y1 - R2*Y2 - VG, (R2 + R3)*C*Y21 + Y2 - R2*C*Y11, Y1, Y2);
sol2 = ilaplace(Sol2);
v1 = int(sol2, 0, t)/C;
v1 = vpa(v1, 2);
% El resultado de la integral es: exp(-22000*t)*(-1.08685e-6*sin(2e5*t)-9.88045e-6*cos(2e5*t))

% Apartado 7.2:
Vg = V*sin(2*pi*t*5e4)*heaviside(t);
VG = laplace(Vg, t, s);
[Sol1, Sol2] = solve(L*Y11 + (R1 + R2)*Y1 - R2*Y2 - VG, (R2 + R3)*C*Y21 + Y2 - R2*C*Y11, Y1, Y2);
sol2 = ilaplace(Sol2);
v2 = int(sol2, 0, t)/C;
v2 = vpa(v2, 2);

% Representacion:
fprintf('\n7. Circuito RLC2\n');
fprintf('  Apartado 7.1:\n');
fprintf('\tSolucion: v(t) = %s\n', char(v1));
fprintf('  Apartado 7.2:\n');
fprintf('\tSolucion: v(t) = %s\n', char(v2));

% Representacion grafica:
h1 = ezplot(char(v2), [0, (R2 + R3)*C*0.5]);
hold on;
h2 = ezplot(char(Vg), [0, (R2 + R3)*C*0.5]);
hold off;
set(h1, 'Color', 'red');
set(h2, 'Color', 'blue');
legend('Vc(t)','Vg(t)');


% Puesto que pasando la segunda ecuacion a su forma canonica se observa de
% inmediato que tau2 = (R2 + R3)*C = 1001*1e-6 ~= 1e-3 y a 5*tau2 se
% entra en regimen permanente. Por problema de densidad de la funcion en la
% grafica, aqui proponemos un T=0.5*tau2 para facilitar la observacion.

%%  8.Analisis del regimen permanente senoidal
clear

% Solucion del Apartado 7.2 del problema 7:
syms t s Y1 Y2;
R1 = 1e-1;
R2 = 1e3;
R3 = 1;
L = 25e-6;
C = 1e-6;
V = 10;
i10 = 0;
i20 = 0;
Y11 = s*Y1 - i10;
Y21 = s*Y2 - i20;
Vg = V*sin(2*pi*t*5e4)*heaviside(t);
VG = laplace(Vg, t, s);
[Sol1, Sol2] = solve(L*Y11 + (R1 + R2)*Y1 - R2*Y2 - VG, (R2 + R3)*C*Y21 + Y2 - R2*C*Y11, Y1, Y2);
sol2 = ilaplace(Sol2);
v2 = int(sol2, 0, t)/C;
v2 = vpa(v2, 2);

% Aqui se usa primero el codigo que proporciona el fichero word del
% enunciado

% Calculo de H(s) del circuito RCL 2
 
syms R1 R2 R3 L C s
syms I1 I2 VG VC
 
% Ecuaciones de las mallas
Eq1=(R1+R2+s*L)*I1-R2*I2-VG;
Eq2=(R3+R2+1/(s*C))*I2-R2*I1;
 
% Solucion de las ecuaciones
[I1, I2]=solve(Eq1,Eq2,I1,I2);
 
% Tensi車n de salida VC(s)
VC=(1/(s*C))*I2;
 
% Funci車n de transferencia H(s)
H=VC/VG;
[N, D]=numden(H);
D=collect(D,s);
H=N/D;
pretty(H)
 
% Sustituci車n de valores circuitales y para s=jw0
H=subs(H,{R1 R2 R3 L C s},{1e-1 1e3 1e0 25e-6 1e-6 1j*2*pi*5e4});

% Atenuaci車n para s=jw0
A=double(abs(H));
% Fase para s=jw0
Fi=double(angle(H));
disp(['|H(jw0)| = ' num2str(A) ' ;  <H(jw0) = ' num2str(Fi)])


% Comprobacion:
EcTeorica = 10*A*sin(2*pi*5e4*t + Fi);
EcTeorica = vpa(EcTeorica, 2);
fprintf('\n8. Analisis del regimen permanente senoidal\n');
fprintf('  Solucion calculada en el problema 7:\t%s\n', char(v2));
fprintf('  Solucion calculada en el problema 8:\t%s\n', char(EcTeorica));
fprintf('  Observando las graficas se ve que el Regimen Senoidal Permanente de las dos soluciones coinciden.\n')

g1 = ezplot(char(v2), [0, 1e-3*0.5]);
hold on;
g2 = ezplot(char(EcTeorica), [0, 1e-3*0.5]);
hold off;
set(g1, 'Color', 'red');
set(g2, 'Color', 'blue');
legend('P.7','P.8');
