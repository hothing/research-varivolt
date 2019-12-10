// Copyright (C) 2018 - Mykhayl Puzanov
//
// Date of creation: 24.04.2018
//
// This is Scilab script.
// See https://www.scilab.org/download
//
// The script calculates the parameters for 
// the state-space controller and Kalman filter.
// The reason why the SS-controller with KF has been used is
// demand to very accurate voltage regulation
// without windings position encoders 
// and with very noised voltage measurement.

// measurement parameters
Umin = 37.6 // [V] minimal voltage
Umax = 107 // [V] maximal voltage

// windings drive parameters
Tp = 0.5 // [s] ]the drive lag time
nm = 1350 // [1/min] nominal motor speed
G = 0.523 // [mm/360 deg] translation factor
OutMax = 100 // [%] maximal manipulation signal
Xmin = 0 // [mm] minimal position
Xmax = 884 // [mm] maximal position

// Kalman filter parameters
Q = [0.5, 0.0; 0.0, 0.0] // Kalman filter weighting matrix 
Usd = 12 // standart deviation of the measured voltage signal
Fsd = 1.0 // noise magnitude factor to check sensivity of model

// SS-Controller parameters
// complex poles
pr = -0.95 // real part sets the reaction time and smoothness
// Note: PR must be in range (-1.0 .. 0.0)
// Small value means faster reaction
// Big value means slower and smoother reaction
// Practical range is (-0.98 .. -0.5)
pi = 0.0 // image part sets the oscillation compensation
// For this type of drive there is no oscillation part

// Options
show_trends = 9

// PLC cycle/sampling time
ST = 0.1 // [s] sampling time

////////////////////////////////////////////////////////////////////////////////
Qlimp = OutMax
Qlimn = -OutMax

// a changing of position is doing by a motor 
// which can be represented as First-Order-Unit
// y(n+1) = a*y(n) + b*u(n+1)
// a = Tp / (Tp + T)
// b = Ku * T / (Tp + T)

Ku = nm / (60 * Qlimp)
a2 = Tp / (Tp + ST)
b2 = Ku * ST / (Tp + ST)

// a transformer translation factor depends from a lineer position 
// of secondary coil
// the linear position {p} is an intergral from the motor speed {w}
// p(n+1) = p(n) + kw * w(n + 1)
a1 = G * ST
b1 = 0

// A state-space matrix

A = [1, a1; 0, a2]
// For a control we use the "voltage" which is changing a motor speed
B = [0; b2]

// measurement matrix
kv = (Umax - Umin)/(Xmax - Xmin)
U0 = Umax - kv*Xmax
// U = kv*P + U0
C = [kv, 0]

// feedforward
D = 0

// check a controlability
// n = 2, R = [B, A*B]
R = [B, A*B]
disp('The controlability test')
disp(rank(R) == 2)

// check observability
// n = 2 
O = [C; C*A]
disp('The observability test')
disp(rank(O) == 2)

p = [-pr + pi*%i, -pr - pi*%i]

Ks = ppol(A,B, p)
Nc = (C^-1 .* Ks)(1)

Ac = A - B*Ks

spec(Ac)

Bc = B * Nc

Sc = syslin(ST, Ac, Bc, C, D)

// simulation

// use short impulse signal for 
//ui = zeros(1, 200); ui(1) = 1;
//us = ones(1, 200); us(1) = 0;
t=0:ST:4*Tmr
n = length(t)
ui = zeros(1, n); ui(1) = 1;
us = ones(1, n); us(1) = 0;
uf = abs(sin(0.02*t))
ug = us + 0.1 * sin(0.05*t)

// assign signal
r = us * 50.0

// initial state
x0 = [Xmin;0]

//////////////////////////////////////////

// closed-loop system (ideal)
//x2 = ltitr(Ac, Bc, u, x0)
//y2 = C*x2 + D*u;
x2 = zeros(2, n); x2(:,1) = x0
y2 = zeros(1, n)
u2 = r
for i=1:n-1 do   
    ux = Nc * r(i) - Ks *x2(:, i) // error signal
    // include a limiter of manipulation signal
    if ux > Qlimp then ux = Qlimp; end
    if ux < Qlimn then ux = Qlimn; end
    u2(i) = ux
    // plant simulation
    y2(:, i) = C * x2(:, i);
    x2(:, i + 1) = A * x2(:, i) + B * u2(i);
end
y2(:, n) = C * x2(:, n);

// closed-loop system with the observer (Kalman filter)

x3 = zeros(2, n); x3(:,1) = x0
y3 = zeros(1, n)
x3hat = zeros(2, n); x3hat(:,1) = x0
y3hat = y3
u3 = r
P = zeros(A)

S = Usd^2
ns3 = grand(1, n, "nor", 0, Usd * Fsd)
y3n = y3
ye = y3

// prepare a Kalman filter gain matrix 
for j=1:600 do
    Pminus = A * P * A' + Q
    Kp = (C * Pminus * C' + S)^-1
    Lk = Pminus * C' * Kp // [n x 1]
    P = (eye(P) - Lk*C) * Pminus
end

// KF stady state test

x3(1,1) = 600
x3(2,1) = 0

for i=1:n-1 do   
    // plant simulation
    y3(:, i) = C * x3(:, i);    
    x3(:, i + 1) = A * x3(:, i);
    // simulate a noised measurement
    y3n(:,i) = y3(:,i) + ns3(:,i)
    
    // controller simulation
    y3hat(:, i) = C * x3hat(:, i);
    // Kalman filter
    x3hatminus = A * x3hat(:, i)
    // estimation error
    ey = y3n(:,i) - C * x3hatminus
    ye(:,i) = ey
    // X expectation
    x3hat(:, i + 1) = x3hatminus + Lk * ey
end
x4 = x3
x4hat = x3hat
y4 = y3
y4n = y3n
y4hat = y3hat

// Controller with KF as observer
x3(1,1) = 0
x3(2,1) = 0

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
    // estimation error
    ey = y3n(:,i) - C * x3hatminus
    ye(:,i) = ey
    // X expectation
    x3hat(:, i + 1) = x3hatminus + Lk * ey

    // calculate manipulation 
    ux = Nc * r(i) - Ks *x3hat(:, i)
    //ux = r(i) // open the system
    // include a limiter of manipulation signal
    if ux > Qlimp then ux = Qlimp; end
    if ux < Qlimn then ux = Qlimn; end
    u3(i+1) = ux
end
y3(:, n) = C * x3(:, n);
y3hat(:, n) = C * x3hat(:, n);

xe3 = x3 - x3hat
ye3 = y3 - y3hat

// Show calculated parameters
disp('The dynamic object model parameters')
disp(A, "A =")
disp(B, "B =")
disp(C, "C =")

disp('Kalman filter matrix')
disp(Lk, "Lk =")

disp('Controller feeadback gains')
disp(Nc, "Nc =")
disp(Ks, "Ks =")

// plotting
if show_trends > 0 then
    f0 = figure();
    title("Wings position: target and actual (by step function)", "fontsize",3)
    plot([t', t', t'], [r', y2', y3'])
    
    if show_trends > 1 then        
        f4 = figure();
        title("Wings positioning error", "fontsize",3)
        subplot(211)
        plot(t', xe3')
        subplot(212)
        plot(t', ye3')
    end
    
    if show_trends > 2 then
        f1 = figure();
        subplot(211)
        title("Wings positioning control: ideal")
        plot([t', t'], [r', y2'])
        subplot(212)
        title("Wings positioning control: noised")
        plot([t', t'], [r', y3'])
        
        f3 = figure();
        subplot(211)
        title("Manipulation: ideal measurements")
        plot(t', u2')
        subplot(212)
        title("Manipulation: noised measurements")
        plot(t', u3')
    end
    
    if show_trends > 3 then
        f2 = figure();
        subplot(221)
        title("Model internals: position with ideal measurement")
        plot(t', x2')
        subplot(222)
        title("Model internals: position with noised measurements")
        plot(t', x3')
        subplot(223)
        title("Model internals: filtered measurement")
        plot([t', t'], [y4', y4hat'])
        subplot(224)
        title("Model internals: direct measurement")
        plot([t', t'], [y4n', y4hat'])
     end
end
