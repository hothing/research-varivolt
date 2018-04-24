ny=2;
nu=3;
nx=4;

Sp = ssrand(ny,nu,nx);
[A,B,C,D]= abcd(Sp);

Kc = -ppol(A,B,[-1,-1,-1,-1]);  //Controller gain
Kf = -ppol(A',C',[-2,-2,-2,-2]);Kf=Kf';    //Observer gain

Sob = obscont(Sp, Kc, Kf)
Scl = Sp /. (-Sob);
spa = spec(Scl('A'))   //closed loop system

[J,r] = obscont(Sp, Kc, Kf);
Q = ssrand(nu,ny,3);
Q('A') = Q('A')-(max(real(spec(Q('A'))))+0.5)*eye(Q('A')) 
//Q is a stable parameter

Sct = lft(J,r,Q); // full system

spct = spec(h_cl(Sp,K))  // closed-loop A matrix (should be stable);

disp(spa)
disp(spct)


// step simulation
t=0:0.05:5;

//deff('u=timefun(t)','u=abs(sin(t))')
deff('u=timefun(t)','u=[1 + 0*t; 1 + 0*t; 1 + 0*t]')
[y, x] =  csim(timefun, t, Sct)
plot2d([t',t'],[y', 0*t'])
