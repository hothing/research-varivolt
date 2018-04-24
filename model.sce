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
B = [0; T]

// measurement 
x1min = 200 // [mm] minimal position
x1max = 1000 // [mm] maximal position
Umin = 30 // [V] minimal voltage
Umax = 95 // [V] maximal voltage

kv = (Umax - Umin)/(x1max - x1min)

C = [kv, 0; 0, 1]
// feedworward
D = [0; 0]

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

// simulation

// use short impulse signal for 
ui = zeros(1, 200); ui(1) = 1;
us = ones(1, 200); ui(1) = 0;

// assign signal
u = us

//opened-loop system
x1=dsimul(So,u)

// closed-loop system
x2=dsimul(Sc,u); 

subplot(121)
plot(x1(1,:))
subplot(122)
plot(x2(1,:))
