// The test of varivolt movements and control

T = 1 // sampling time

// X = [x1; x2] where x1 - is position, x2 - is velocity
// x1 has range 200 .. 1000mm 
// we cannot measure the position directly, but we have an indirect measurement by voltage
// U = kv * x1
// U has range [30V .. 95V]
A = [1, T; 0, 1]
// By fact, for control we can use only a speed setpoint, but!
// We cannot a motor can be represented with first-oder unit and it has a ramp function
// the simplest way to account this is to use an acceleration as control signal in our model
B = [T^2/2; T]
// B = [0; T]

// measurement 
x1min = 200 // [mm] minimal position
x1max = 1000 // [mm] maximal position
Umin = 30 // [V] minimal voltage
Umax = 95 // [V] maximal voltage

kv = (Umax - Umin)/(x1max - x1min)

//C = [kv, 0; 0, 0]
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
Ns = 12.1  // scaling factor
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
[y1, x1] = flts(u, So, x0)

// closed-loop system (ideal))
[y2, x2] = flts(u, Sc, x0); 

// closed-loop system with the observer
xco0 = [0;0;0;0]
[y3, x3] = flts(u, Sco, xco0); 

subplot(311)
plot(y1)
subplot(312)
plot(y2)
subplot(313)
plot(y3)
