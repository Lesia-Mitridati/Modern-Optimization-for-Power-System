OPTIONS mpec = KNITRO;

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

EQUATIONS
OF,constr_UL,constr_LL1,constr_LL2,constr_LL3,constr_LL4,constr_lagrange;

*primal constraints
OF.. Obj_leader =e= -y;
constr_UL.. x =g= 2;
constr_LL1..  x+y =g= 7;
constr_LL2..  -y =g= -9;
constr_LL3.. -x+2*y =g= -1;
constr_LL4..  -6*x - y =g= -45;

*KKT conditions
constr_lagrange.. 1 - mu_1 + mu_2 - 2*mu_3 + mu_4 =e= 0;

MODEL MPEC
/OF
constr_UL
constr_LL1.mu_1
constr_LL2.mu_2
constr_LL3.mu_3
constr_LL4.mu_4
constr_lagrange/
;

SOLVE MPEC using mpec min Obj_leader;
Display x.l,y.l;