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

nm = 1350 // [1/min] nominal speed

Tp = 2
Ku = nm / 100
a2 = Tp / (Tp + T)
b2 = Ku * T / (Tp + T)

// a transformer translation factor depends from a lineer position of secondary coil
// the linear position {p} is an intergral from the motor speed {w}
// p(n+1) = p(n) + kw * w(n + 1)


Tmr = 120 // [s] a movement time from end to end

vx = (Pmax - Pmin) / Tmr // nominal linear speed
kw = vx / nm // speed translation factor
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

pr = -0.9 // -0.1 is good
pi = 0.0
p = [-pr + pi*%i, -pr - pi*%i]

K = ppol(A,B, p)

Ac = A - B*K
spec(Ac)
Ns = 5.5384616  // scaling factor
Bc = B * Ns

Sc = syslin(T, Ac, Bc, C, D)

// Controller with an observer
// an observer character
opr = -0.5 // 
opi = 0.0 //
op = [-opr + opi*%i, -opr - opi*%i]
L = ppol(A', C', op)'

Ao = A - L*C
Bo = B

// simulation

// use short impulse signal for 
ui = zeros(1, 200); ui(1) = 1;
us = ones(1, 200); us(1) = 0;

// assign signal
u = us * Umax

// initial state
x0 = [Pmin;0]

//opened-loop system
x1 = ltitr(A, B, u, x0)
y1 = C*x1 + D*u;

// closed-loop system (ideal)
x2 = ltitr(Ac, Bc, u, x0)
y2 = C*x2 + D*u;

// closed-loop system with the observer
Nco = 7.47

x3 = zeros(2, length(u)); x3(:,1) = x0
y3 = zeros(1,length(u))
x3hat = zeros(2, length(u)); //x3hat(:,1) = x0
y3hat = y3
m3 = u
n = length(u)
for i=1:n-1 do   
    m3(i) = Nco*u(i) - K *x3hat(:, i) // error signal
    // plant simulation
    y3(:, i) = C * x3(:, i);
    x3(:, i + 1) = A * x3(:, i) + B * m3(i);
    // controller simulation
    y3hat(:, i) = C * x3hat(:, i);
    ye = y3(:, i) - y3hat(:, i);
    x3hat(:, i + 1) = Ac * x3hat(:, i) + B * u(i) + L * ye;
end
y3(:, n) = C * x3(:, n);
y3hat(:, n) = C * x3hat(:, n);

m4 = - K *x3hat

f1 = figure();
subplot(221)
plot(y1)
subplot(222)
plot(y2)
subplot(223)
plot(y3)
subplot(224)
plot(y3hat)
f2 = figure();
subplot(221)
plot(x1')
subplot(222)
plot(x2')
subplot(223)
plot(x3')
subplot(224)
plot(x3hat')
