## The script calculates the parameters for
## the state-space controller and Kalman filter.
## The reason why the SS-controller with KF has been used is
## demand to very accurate voltage regulation
## without windings position encoders
## and with very noised voltage measurement.

## measurement parameters
Umin = 37.6 ## [V] minimal voltage
Umax = 107 ## [V] maximal voltage

## windings drive parameters
Tp = 0.5 ## [s] ]the drive lag time
nm = 1350 ## [1/min] nominal motor speed
G = 0.523 ## [mm/360 deg] translation factor
OutMax = 100 ## [%] maximal manipulation signal
Xmin = 0 ## [mm] minimal position
Xmax = 884 ## [mm] maximal position

## Kalman filter parameters
Q = [0.5, 0.0; 0.0, 0.0] ## Kalman filter weighting matrix
Usd = 12 ## standart deviation of the measured voltage signal
Fsd = 1.0 ## noise magnitude factor to check sensivity of model

## Options
show_trends = 9

## PLC cycle/sampling time
ST = 0.1 ## [s] sampling time

################################################################################
Qlimp = OutMax
Qlimn = -OutMax

## a changing of position is doing by a motor
## which can be represented as First-Order-Unit
## y(n+1) = a*y(n) + b*u(n+1)
## a = Tp / (Tp + T)
## b = Ku * T / (Tp + T)

Ku = nm / (60 * Qlimp)
a2 = Tp / (Tp + ST)
b2 = Ku * ST / (Tp + ST)

## a transformer translation factor depends from a lineer position
## of secondary coil
## the linear position {p} is an intergral from the motor speed {w}
## p(n+1) = p(n) + kw * w(n + 1)
a1 = G * ST
b1 = 0

## A state-space matrix

A = [1, a1; 0, a2]
## For a control we use the "voltage" which is changing a motor speed
B = [0; b2]

## measurement matrix
kv = (Umax - Umin)/(Xmax - Xmin)
U0 = Umax - kv*Xmax
## U = kv*P + U0
C = [kv, 0]

## feedforward
D = 0

## check a controlability
## n = 2, R = [B, A*B]
R = [B, A*B]
disp('The controlability test')
disp(rank(R) == 2)

## check observability
## n = 2
O = [C; C*A]
disp('The observability test')
disp(rank(O) == 2)

## simulation

Tmr = 60

t=0:ST:4*Tmr;
n = length(t);
ui = zeros(1, n); ui(1) = 1;
us = ones(1, n); us(1) = 0;
uf = abs(sin(0.02*t));
ug = us + 0.1 * sin(0.05*t);

## initial state
x0 = [Xmin;0];

##########################################

x1 = zeros(2, n);
x1(:,1) = x0;
y1 = zeros(1, n);
r1 = zeros(1, n);
k1 = int16(n / 2);
r1(2:2+k1) = ones(1, k1 +  1)  .* 100;
u1 = zeros(size(r1));

for i=1:n-1
    ux = r1(i);
    ## include a limiter of manipulation signal
    if ux > Qlimp
      ux = Qlimp;
    endif
    if ux < Qlimn
      ux = Qlimn;
    endif
    if ((x1(1, i) >= Xmax) && (ux > 0)) || ((x1(1, i) <= Xmin) && (ux < 0))
      ux  = 0;
    endif
    u1(i) = ux;
    ## plant simulation
    x1(:, i + 1) = A * x1(:, i) + B * u1(i);
    y1(:, i) = C * x1(:, i);
endfor


############################################

## closed-loop system (ideal)

disp('*** Closed-loop system ***')
pkg load control

## SS-Controller parameters
## complex poles
pr = -0.95 ## real part sets the reaction time and smoothness
## Note: PR must be in range (-1.0 .. 0.0)
## Small value means faster reaction
## Big value means slower and smoother reaction
## Practical range is (-0.98 .. -0.5)
pi = 0.0 ## image part sets the oscillation compensation
## For this type of drive there is no oscillation part
p = [-pr + pi*i, -pr - pi*i]
Ks = place(A, B, p)
Nc = ((C.^-1)'*Ks)(1)

Ac = A - B * Ks
Bc = B * Nc

## check a controlability
## n = 2, R = [B, A*B]
R2 = [Bc, Ac*Bc];
disp('The controlability test of closed-loop system')
disp(rank(R2) == 2)

## check observability
## n = 2
O2 = [C; C*Ac];
disp('The observability test of closed-loop system')
disp(rank(O2) == 2)

## assign signal
r = us * 50.0;
x2 = zeros(2, n);
x2(:,1) = x0;
y2 = zeros(1, n);
u2 = zeros(1, n);
for i=1:n-1
    ux = Nc * r(i) - Ks * x2(:, i); ## error signal
    ## include a limiter of manipulation signal
    if ux > Qlimp
      ux = Qlimp;
    endif
    if ux < Qlimn
      ux = Qlimn;
    endif
    if ((x2(1, i) >= Xmax) && (ux > 0)) || ((x2(1, i) <= Xmin) && (ux < 0))
      ux  = 0;
    endif
    u2(i) = ux;
    ## plant simulation
    y2(:, i) = C * x2(:, i);
    x2(:, i + 1) = A * x2(:, i) + B * u2(i);
endfor

y2(:, n) = C * x2(:, n);

## closed-loop system with the observer (Kalman filter)
disp('*** Closed-loop system with the observer (Kalman filter) ***')

x3 = zeros(2, n);
x3(:,1) = x0;
y3 = zeros(1, n);
x3hat = zeros(2, n);
 x3hat(:,1) = x0;
y3hat = y3;
u3 = r;
P = zeros(size(A));

S = Usd^2
ns3 = randn(1, n) .* (Usd * Fsd);
y3n = y3;
ye = y3;

## prepare a Kalman filter gain matrix
for j=1:600
    Pminus = A * P * A' + Q;
    Kp = (C * Pminus * C' + S)^-1;
    Lk = Pminus * C' * Kp; ## [n x 1]
    P = (eye(size(P)) - Lk*C) * Pminus;
endfor

## KF stady state test

x3(1,1) = 600;
x3(2,1) = 0;

for i=1:n-1
    ## plant simulation
    y3(:, i) = C * x3(:, i);
    x3(:, i + 1) = A * x3(:, i);
    ## simulate a noised measurement
    y3n(:,i) = y3(:,i) + ns3(:,i);

    ## controller simulation
    y3hat(:, i) = C * x3hat(:, i);
    ## Kalman filter
    x3hatminus = A * x3hat(:, i);
    ## estimation error
    ey = y3n(:,i) - C * x3hatminus;
    ye(:,i) = ey;
    ## X expectation
    x3hat(:, i + 1) = x3hatminus + Lk * ey;
endfor
x4 = x3;
x4hat = x3hat;
y4 = y3;
y4n = y3n;
y4hat = y3hat;

## Controller with KF as observer
x3(1,1) = 0;
x3(2,1) = 0;

for i=1:n-1
    ## plant simulation
    y3(:, i) = C * x3(:, i);
    x3(:, i + 1) = A * x3(:, i) + B * u3(i);
    ## simulate a noised measurement
    y3n(:,i) = y3(:,i) + ns3(:,i);

    ## controller simulation
    y3hat(:, i) = C * x3hat(:, i);
    ## Kalman filter
    x3hatminus = A * x3hat(:, i) + B * u3(:, i);
    ## estimation error
    ey = y3n(:,i) - C * x3hatminus;
    ye(:,i) = ey;
    ## X expectation
    x3hat(:, i + 1) = x3hatminus + Lk * ey;

    ## calculate manipulation
    ux = Nc * r(i) - Ks *x3hat(:, i);
    ##ux = r(i) ## open the system
    ## include a limiter of manipulation signal
    if ux > Qlimp
      ux = Qlimp;
    endif
    if ux < Qlimn
        ux = Qlimn;
    endif
    u3(i+1) = ux;
endfor
y3(:, n) = C * x3(:, n);
y3hat(:, n) = C * x3hat(:, n);

xe3 = x3 - x3hat;
ye3 = y3 - y3hat;

## Show calculated parameters
disp('The dynamic object model parameters')
disp("A =")
disp(A)
disp("B =")
disp(B)
disp("C =")
disp(C)

disp('Kalman filter matrix')
disp("Lk =")
disp(Lk)

disp('Controller feeadback gains')
disp("Nc =")
disp(Nc)
disp("Ks =")
disp(Ks)

## plotting
if show_trends > 0
  f0 = figure();

  subplot(311)
  plot(t', u1', "-;MNL;")
  xlabel ("time (sec)");
  ylabel ("Control Signal");

  subplot(312)
  plot(t', x1(1,:)')
  xlabel ("time (sec)");
  ylabel ("Wing Position, mm (X[1])");

  subplot(313)
  plot(t', x1(2,:)' .* 60)
  xlabel ("time (sec)");
  ylabel ("Motor speed, rpm (X[2])");

  if show_trends >= 1
    f1 = figure();
    subplot(311)

    plot(t', u2', "-;MNL;")
    xlabel ("time (sec)");
    ylabel ("Control Signal");

    subplot(312)
    plot(t', x2(1,:)')
    xlabel ("time (sec)");
    ylabel ("Wing Position, mm (X[1])");

    subplot(313)
    plot(t', r', "-r;Vset;", t', y2', "-b;Vreal;")
    xlabel ("time (sec)");
    ylabel ("Voltage, V (Y)");
  endif

  if show_trends >= 2
    f2 = figure();
    plot(t', r', "-;Set;", t', y2', "-b;Ideal Ctrl;",  t', y3', "-r;SS/OKF Ctrl;")
    xlabel ("time (sec)");
    ylabel ("Position (Y)");
    title("Wings position: target and actual (by step function)", "fontsize", 13)
  endif

  if show_trends >= 3
    f3 = figure();
    subplot(211)
    plot(t', xe3')
    title("Wings positioning error (state)")

    subplot(212)
    plot(t', ye3')
    title("Wings positioning error (meas)")
  endif

  if show_trends >= 4
      f4 = figure();
      subplot(211)
      plot([t', t'], [r', y2'])
      title("Wings positioning control: ideal")
      subplot(212)
      plot([t', t'], [r', y3'])
      title("Wings positioning control: noised")

      f4b = figure();
      subplot(211)
      plot(t', u2')
      title("Manipulation: ideal measurements")

      subplot(212)
      plot(t', u3')
      title("Manipulation: noised measurements")

  endif

  if show_trends >= 5
      f5 = figure();
      subplot(221)
      plot(t', x2')
      title("Model internals: position with ideal measurement")

      subplot(222)
      plot(t', x3')
      title("Model internals: position with noised measurements")

      subplot(223)
      plot([t', t'], [y4', y4hat'])
      title("Model internals: filtered measurement")

      subplot(224)
      plot([t', t'], [y4n', y4hat'])
      title("Model internals: filtered measurement")
   endif
endif

