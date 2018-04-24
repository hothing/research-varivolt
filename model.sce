// The test of varivolt movements and control

T = 1 // sampling time

// X = [x1; x2] where x1 - is position, x2 - is velocity
// x1 has range 200 .. 1000mm 
// we cannot measure the position directly, but we have an indirect measurement by voltage
// U = kv * x1
// U has range [30V .. 95V]

// measurement 
Pmin = 200 // [mm] minimal position
Pmax = 1000 // [mm] maximal position
Umin = 30 // [V] minimal voltage
Umax = 95 // [V] maximal voltage

// a changing of position is doing by a motor wich can be represented as First-Order-Unit
// y(n+1) = a*y(n) + b*u(n+1)
// a = Tp / (Tp + T)
// b = Ku * T / (Tp + T)

Tp = 2
Ku = 1
a2 = Tp / (Tp + T)
b2 = Ku * T / (Tp + T)

// a transformer translation factor depends from a lineer position of secondary coil
// the linear position {p} is an intergral from the motor speed {w}
// p(n+1) = p(n) + kw * w(n + 1)

nm = 1350 // [1/min] nominal speed
Tmr = 120 // [s] a movement time from end to end

kw = ((Pmax - Pmin) / Tmr) / nm
a1 = kw * T
b1 = 0

// A state-space matrix

A = [1, a1; 0, a2]
// For a control we use the "voltage" which is changing a motor speed
B = [0; b2]

// measurement matrix
kv = (Umax - Umin)/(Pmax - Pmin)
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
//K = [22, 2]
p = [-3 , -1]

pr = -0.01 // -0.1 is good
pi = 0.0
p = [-pr + pi*%i, -pr - pi*%i]

K = ppol(A,B, p)

Ac = A - B*K
spec(Ac)
Ns = 7328.132  // scaling factor
Bc = B * Ns

Sc = syslin(T, Ac, Bc, C, D)

// an observer making
opr = 100 // 
opi = 0.0 //
op = [-opr + opi*%i, -opr - opi*%i]
L = ppol(A', C', op)'

At1 = A - B*K
At2 = B * K
At3 = zeros(A)
At4 = A - L*C

At = [At1, At2; At3, At4];

Bt = [Bc; zeros(B) ];

Ct = [C, zeros(C) ];

Sco = syslin(T, At, Bt, Ct, D)

// simulation

// use short impulse signal for 
ui = zeros(1, 200); ui(1) = 1;
us = ones(1, 200); ui(1) = 0;

// assign signal
u = us

// initial state
x0 = [0;0]

//opened-loop system
x1 = ltitr(A, B, u, x0)
y1 = C*x1 + D*u;

// closed-loop system (ideal))
x2 = ltitr(Ac, Bc, u, x0)
y2 = C*x2 + D*u;

// closed-loop system with the observer
xco0 = [0;0;0;0]
x3 = ltitr(At, Bt, u, xco0)
y3 = Ct*x3 + D*u;

subplot(311)
plot(y1)
subplot(312)
plot(y2)
subplot(313)
plot(y3)
