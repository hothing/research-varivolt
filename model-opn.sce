// Copyright (C) 2018 - Mykhayl Puzanov
//
// Date of creation: 24.04.2018
//
// The test of varivolt movements and control

T = 0.1 // sampling time

// X = [x1; x2] where x1 - is position, x2 - is velocity
// x1 has range 200 .. 1000mm 
// we cannot measure the position directly, but we have an indirect measurement by voltage
// U = kv * x1
// U has range [30V .. 95V]

// measurement 
Pmin = 0 // [mm] minimal position
Pmax = 1000 // [mm] maximal position
Umin = 30 // [V] minimal voltage
Umax = 95 // [V] maximal voltage

// a changing of position is doing by a motor wich can be represented as First-Order-Unit
// y(n+1) = a*y(n) + b*u(n+1)
// a = Tp / (Tp + T)
// b = Ku * T / (Tp + T)

nm = 1350 // [1/min] nominal speed
nms = nm / 60 // nominal speed [1/s]

Tp = 2
Qmax = 100
Qlimp = Qmax
Qlimn = -Qmax

Ku = (nms) / (Qmax)

a2 = Tp / (Tp + T)
b2 = Ku * T / (Tp + T)

// a transformer translation factor depends from a lineer position of secondary coil
// the linear position {p} is an intergral from the motor speed {w}
// p(n+1) = p(n) + kw * w(n + 1)


Tmr = 120 // [s] a movement time from end to end

vx = (Pmax - Pmin) / Tmr // nominal linear speed
kw = vx / nms // speed translation factor
a1 = kw * T
b1 = 0

// A state-space matrix

A = [1, a1; 0, a2]
// For a control we use the "voltage" which is changing a motor speed
B = [0; b2]

// measurement matrix
kv = (Umax - Umin)/(Pmax - Pmin)
U0 = Umax - kv*Pmax
// U = kv*P + U0
C = [kv, 0]

// feedworward
//D = [0; 0]
D = 0

// check a controlability
// n = 2, R = [B, A*B]
R = [B, A*B]
disp(rank(R))

// check observability
// n = 2 
O = [C; C*A]
disp(rank(O))

// opened-loop system
So = syslin(T, A, B, C, D)

// closed-loop model
// complex poles
pr = -0.87 // is good
pi = 0.1
p = [-pr + pi*%i, -pr - pi*%i]
Ns = 14.898461 // scaling factor, acceptable only for these poles

K = ppol(A,B, p)

Ac = A - B*K

spec(Ac)

Bc = B * Ns

Sc = syslin(T, Ac, Bc, C, D)

// Controller with an observer
// an observer character
op = [0.98; 0.67]
L = ppol(A', C', op)'

Nc = 6.9374759
Ao = A - L*C


// simulation

// simulate reference signals
Tmax = Tmr
t=0:T:Tmax
n = length(t)
ui = zeros(1, n); ui(1) = 1;
us = ones(1, n); us(1) = 0;
uf = abs(sin(0.02*t))


// assign signal
r = us * 100.0

// initial state
x0 = [Pmin;0]

//opened-loop system
x1 = ltitr(A, B, r, x0)
y1 = C*x1 + D*r;

// plotting

f1 = figure();
subplot(211)
plot([t', t'], [r', y1'])
subplot(212)
plot(t', x1')
