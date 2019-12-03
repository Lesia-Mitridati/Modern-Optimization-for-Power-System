
PARAMETERS
M /100/
;

FREE VARIABLES
x
y
Obj_leader
;

POSITIVE VARIABLES
mu_1
mu_2
mu_3
mu_4
;

BINARY VARIABLES
u1
u2
u3
u4
;

EQUATIONS
OF,constr_UL,constr_LL1,constr_LL2,constr_LL3,constr_LL4
constr_lagrange,KKT_11,KKT_12,KKT_21,KKT_22,KKT_31,KKT_32,KKT_41,KKT_42;

*primal constraints
OF.. Obj_leader =e= -y;
constr_UL.. x =g= 2;
constr_LL1..  x+y =g= 7;
constr_LL2..  -y =g= -9;
constr_LL3.. -x+2*y =g= -1;
constr_LL4..  -6*x - y =g= -45;

*KKT conditions
constr_lagrange.. 1 - mu_1 + mu_2 - 2*mu_3 + mu_4 =e= 0;

KKT_11..  x+y-7 =l= M*(1-u1);
KKT_12..  mu_1 =l= M*u1;

KKT_21..  -y+9 =l= m*(1-u2);
KKT_22..  mu_2 =l= M*u2;

KKT_31..  -x+2*y+1 =l= M*(1-u3);
KKT_32..  mu_3 =l= M*u3;

KKT_41..  -6*x - y+45 =l= M*(1-u4);
KKT_42..  mu_4 =l= M*u4;


MODEL MPEC
/OF,constr_UL,constr_LL1,constr_LL2,constr_LL3,constr_LL4
constr_lagrange,KKT_11,KKT_12,KKT_21,KKT_22,KKT_31,KKT_32,KKT_41,KKT_42/
;

SOLVE MPEC using MIP Minimizing Obj_leader;
Display x.l,y.l;