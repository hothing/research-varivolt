// Copyright (C) 2018 - Mykhayl Puzanov
//
// Date of creation: 24.04.2018
//
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
Qmin = 0
Qmax = 100
Ku = (nm - 0) / (Qmax -Qmin)
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
Ns = 14.898461 // acceptable only for these poles


p = [-pr + pi*%i, -pr - pi*%i]

K = ppol(A,B, p)

Ac = A - B*K

spec(Ac)
//Ns = 5.5384616  // scaling factor

Bc = B * Ns

Sc = syslin(T, Ac, Bc, C, D)

// simulation

// use short impulse signal for 
ui = zeros(1, 200); ui(1) = 1;
us = ones(1, 200); us(1) = 0;

// assign signal
r = us * Umax
n = length(r)

// initial state
x0 = [Pmin;0]

//opened-loop system
x1 = ltitr(A, B, r, x0)
y1 = C*x1 + D*r;

// closed-loop system (ideal)
//x2 = ltitr(Ac, Bc, u, x0)
//y2 = C*x2 + D*u;
x2 = zeros(2, n); x2(:,1) = x0
y2 = zeros(1, n)
u2 = r
for i=1:n-1 do   
    ux = Ns * r(i) - K *x2(:, i) // error signal
    // include a limiter of manipulation signal
    if ux > Qmax then ux = Qmax; end
    if ux < Qmin then ux = Qmin; end
    u2(i) = ux
    // plant simulation
    y2(:, i) = C * x2(:, i);
    x2(:, i + 1) = A * x2(:, i) + B * u2(i);
end
y2(:, n) = C * x2(:, n);

// closed-loop system with the observer (Kalman filter)

// Kalman filter weighting matrixes 
Usd = 16
Q = [0.1, 0; 0, 0.01]
S = Usd^2
Nf = 11.693973

t=0:T:T*3600
deff('u=ureference(t)','u=abs(sin(0.01*t))')

us = ureference();
r = us * 50
n = length(r)

x3 = zeros(2, n); x3(:,1) = x0
y3 = zeros(1, n)
x3hat = zeros(2, n); //x3hat(:,1) = x0
y3hat = y3
u3 = r
P = zeros(A)

ns3 = grand(1, n, "nor", 0, Usd)
y3n = y3
ye = y3

// prepare a Kalman filter gain matrix 
for j=1:600 do
    Pminus = A * P * A' + Q
    Kp = (C * Pminus * C' + S)^-1
    Kk = Pminus * C' * Kp // [n x 1]
    P = (eye(P) - Kk*C) * Pminus
end

for i=1:n-1 do   
    // plant simulation
    y3(:, i) = C * x3(:, i);    
    x3(:, i + 1) = A * x3(:, i) + B * u3(i);
    // simulate a noised measurement
    y3n(:,i) = y3(:,i) + ns3(:,i)
    
    // controller simulation
    y3hat(:, i) = C * x3hat(:, i);
    // Kalman filter
    x3hatminus = A * x3hat(:, i) + B * u3(:, i)
    //Pminus = A * P * A' + Q
    //Kp = (C * Pminus * C' + S)^-1
    //Kk = Pminus * C' * Kp // [n x 1]
    //P = (eye(P) - Kk*C) * Pminus
    // estimation error
    ey = y3n(:,i) - C * x3hatminus
    ye(:,i) = ey
    // X expectation
    x3hat(:, i + 1) = x3hatminus + Kk * ey

    // calculate manipulation 
    ux = Nf * r(i) - K *x3hat(:, i)
    //ux = r(i) // open the system
    // include a limiter of manipulation signal
    if ux > Qmax then ux = Qmax; end
    if ux < Qmin then ux = Qmin; end
    u3(i+1) = ux
end
y3(:, n) = C * x3(:, n);
y3hat(:, n) = C * x3hat(:, n);

xe3 = x3 - x3hat

// plotting

f1 = figure();
subplot(221)
plot(y1)
subplot(222)
plot(y2)
subplot(223)
plot(y3)
subplot(224)
plot(y3n)

f2 = figure();
subplot(221)
plot(x1')
subplot(222)
plot(x2')
subplot(223)
plot(x3')
subplot(224)
plot(x3hat')

f3 = figure();
subplot(121)
plot(u2)
subplot(122)
plot(u3)

f4 = figure();
plot(xe3')
