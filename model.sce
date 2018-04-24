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

// measurement 
x1min = 200 // [mm] minimal position
x1max = 1000 // [mm] maximal position
Umin = 30 // [V] minimal voltage
Umax = 95 // [V] maximal voltage

kv = (Umax - Umin)/(x1max - x1min)

C = [kv, 0]
// feedworward
D = [0]

// check a controlability
// n = 2, R = [B, A*B]
R = [B, A*B]
rank(R) == 2

// check observability
// n = 2 
O = [C; C*A]
rank(O) == 2

// colsed-loop model
K = 0.1

Ac = A - B*K*C

So = syslin('d', A, B, C, D)

Sc = syslin('d', Ac, B, C, D)


