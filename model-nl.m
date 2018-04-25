% Copyright (C) 2018 - Mykhayl Puzanov
%
% Date of creation: 24.04.2018
%
% The test of varivolt movements and control

T = 1; % sampling time

% X = [x1; x2] where x1 - is position, x2 - is velocity
% x1 has range 200 .. 1000mm 
% we cannot measure the position directly, but we have an indirect measurement by voltage
% U = kv * x1
% U has range [30V .. 95V]

% measurement 
Pmin = 200; % [mm] minimal position
Pmax = 1000; % [mm] maximal position
Umin = 30; % [V] minimal voltage
Umax = 95; % [V] maximal voltage

% a changing of position is doing by a motor wich can be represented as First-Order-Unit
% y(n+1) = a*y(n) + b*u(n+1)
% a = Tp / (Tp + T)
% b = Ku * T / (Tp + T)

nm = 1350; % [1/min] nominal speed

Tp = 2;

Qmax = 100;
Qlimp = Qmax;
Qlimn = -Qmax;

Ku = (nm) / (Qmax);

a2 = Tp / (Tp + T);
b2 = Ku * T / (Tp + T);

% a transformer translation factor depends from a lineer position of secondary coil
% the linear position {p} is an intergral from the motor speed {w}
% p(n+1) = p(n) + kw * w(n + 1)


Tmr = 120; % [s] a movement time from end to end

vx = (Pmax - Pmin) / Tmr; % nominal linear speed
kw = vx / nm; % speed translation factor
a1 = kw * T;
b1 = 0;

% A state-space matrix

A = [1, a1; 0, a2]
% For a control we use the "voltage" which is changing a motor speed
B = [0; b2]

% measurement matrix
kv = (Umax - Umin)/(Pmax - Pmin);
U0 = Umax - kv*Pmax;
% U = kv*P + U0
C = [kv, 0]

% feedworward
%D = [0; 0]
D = 0

% check a controlability
% n = 2, R = [B, A*B]
R = [B, A*B]
disp(rank(R))

% check observability
% n = 2 
O = [C; C*A]
disp(rank(O))

% opened-loop system
So = ss(A, B, C, D, T);

% closed-loop model
% complex poles
pr = -0.87; % is good
pi = 0.1;
p = [complex(-pr, +pi), complex(-pr, -pi)]
Ns = 14.898461 % scaling factor, acceptable only for these poles

K = place(A,B, p)

Ac = A - B*K

eig(A)
eig(Ac)

Bc = B * Ns

Sc = ss(Ac, Bc, C, D, T);

% Controller with an observer
% an observer character
op = [0.98; 0.67]
L = place(A', C', op)'

Nc = 6.9374759
Ao = A - L*C


% simulation

% simulate reference signals
t=0:T:T*360;
n = length(t);
ui = zeros(1, n); ui(1) = 1;
us = ones(1, n); us(1) = 0;
uf = abs(sin(0.02*t));


% assign signal
r = us * (Umax - Umin) * 0.5;

% initial state
x0 = [Pmin;0];

%opened-loop system
[y1, t1, x1] = lsim(So, r, t, x0);
[y2a, t2a, x2a] = lsim(Sc, r, t, x0);

% closed-loop system (ideal)
%x2 = ltitr(Ac, Bc, u, x0)
%y2 = C*x2 + D*u;
x2 = zeros(2, n); x2(:,1) = x0;
y2 = zeros(1, n);
u2 = r;
for i=1:n-1   
    ux = Ns * r(i) - K *x2(:, i); % error signal
    % include a limiter of manipulation signal
    if ux > Qlimp 
        ux = Qlimp; 
    end
    if ux < Qlimn
        ux = Qlimn; 
    end
    u2(i) = ux;
    % plant simulation
    y2(:, i) = C * x2(:, i);
    x2(:, i + 1) = A * x2(:, i) + B * u2(i);
end
y2(:, n) = C * x2(:, n);

% plotting

f1 = figure();
subplot(211)
plot(t, r, t2a, y2a)
subplot(212)
plot(t, r, t, y2)


