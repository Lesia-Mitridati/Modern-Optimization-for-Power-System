option lp=cplex;

sets
IP iterations /1*100/
s scenario /s1*s4/
alias (IP,K)
;

parameters
SetP1_DA(K)         Set of values obtained for variable P1 over iterations
SetP2_DA(K)         Set of values obtained for variable P2 DA over iterations
SetPW_DA(K)         Set of values obtained for variable PW DA over iterations
SetZ_MP(K)          Set of values obtained for obj fun in DA over iterations
SetP2_RT(s,K)         Set of values obtained for variable P2 RT over iterations
SetP_spill(s,K)       Set of values obtained for variable P spill over iterations
SetP_shed(s,K)        Set of values obtained for variable P shed over iterations
SetZ_sub(s,K)         Set of values obtained for subproblem's objective function over iterations
SetDual1(s,K)     Set of values obtained for dual variable \rho P1 over iterations
SetDual2(s,K)     Set of values obtained for dual variable \rho P2 over iterations
SetDualW(s,K)     Set of values obtained for dual variable \rho PW over iterations
SetA(K)         Set of values obtained for auxiliary variable \alpha over iterations
SetUB(K)        Set of values of UB over iterations
W_scenario(s)   wind scenarios /s1 30,s2 60,s3 70,s4 10/
;

scalars
d   convergence gap
m   order of current iteration
W_sub   Real Time wind realization parameter in subproblems
;

* Initial convergence gap
d=1000;

************ Master problem **************
variable
z_MP
a
;
positive variables
P1_DA
P2_DA
PW_DA
;
equation  OF,const_master1,const_master2,const_master3,const_master4,cut,a_down;

OF..            z_MP =e= 10*P1_DA + 20*P2_DA + a;
const_master1..      P1_DA =l= 100;
const_master2..      P2_DA =l= 30;
const_master3..      PW_DA =l= 70;
const_master4..      P1_DA + P2_DA + PW_DA =e= 120;
cut(K)$(m>1 and ord(K)<m).. a =g= 0.25*sum(s,20*SetP2_RT(s,K)+80*SetP_shed(s,K))
                                    + sum(s,SetDual1(s,K)*(P1_DA-SetP1_DA(K)))
                                    + sum(s,SetDual2(s,K)*(P2_DA-SetP2_DA(K)))
                                    + sum(s,SetDualW(s,K)*(PW_DA-SetPW_DA(K)));
a_down.. a =g= -1000000;

model master /OF,const_master1,const_master2,const_master3,const_master4,cut,a_down/ ;

************ Subproblems**************
variable
z_sub
P1_sub
P2_sub
PW_sub
P2_RT
;
positive variables
P_spill
P_shed
;
equation Sub_OF,cons_sub1,cons_sub2,cons_sub3,cons_sub4,cons_sub5,cons_sub6,cons_sub7,cons_sub8,cons_sub9,cons_sub10;

Sub_OF..         z_sub  =e= 0.25*(20*P2_RT + 80*P_shed);
cons_sub1..      P_spill =l= W_sub;
cons_sub2..      P_shed =l= 120;
cons_sub3..      P2_sub + P2_RT =l= 30;
cons_sub4..      P2_sub + P2_RT =g= 0;
cons_sub5..      P2_RT =l= 30;
cons_sub6..      P2_RT =g= -30;
cons_sub7..      P2_RT + (W_sub - PW_sub - P_spill) + P_shed =e= 0;
cons_sub8..   P1_sub =e= P1_DA.l;
cons_sub9..   P2_sub =e= P2_DA.l;
cons_sub10..  PW_sub =e= PW_DA.l;

model sub /Sub_OF,cons_sub1,cons_sub2,cons_sub3,cons_sub4,cons_sub5,cons_sub6,cons_sub7,cons_sub8,cons_sub9,cons_sub10/ ;

************* Iteration **********

SetP1_DA('1')=0;
SetP2_DA('1')=0;
SetPW_DA('1')=0;
SetZ_MP('1')=0;
loop(s,
         SetZ_sub(s,'1')=0 ;
         SetP2_RT(s,'1')=0;
         SetP_spill(s,'1')=0;
         SetP_shed(s,'1')=0;
         SetDual1(s,'1')= 0;
         SetDual2(s,'1')=0 ;
         SetDualW(s,'1')=0;
);

loop(IP$(d gt 0.1),
       m=ord(IP);

       solve master using lp minimizing z_MP;
       SetP1_DA(IP)=P1_DA.l;
       SetP2_DA(IP)=P2_DA.l ;
       SetPW_DA(IP)=PW_DA.l ;
       SetZ_MP(IP)=z_MP.l ;
       SetA(IP)=a.l;

       loop(s,
         W_sub =  W_scenario(s)

         solve sub using lp minimizing z_sub;

         SetZ_sub(s,IP)=z_sub.l ;
         SetP2_RT(s,IP)=P2_RT.l;
         SetP_spill(s,IP)=P_spill.l;
         SetP_shed(s,IP)=P_shed.l ;
         SetDual1(s,IP)= cons_sub8.m  ;
         SetDual2(s,IP)=cons_sub9.m  ;
         SetDualW(s,IP)=cons_sub10.m  ;
       );
       SetUB(IP)=sum(s,SetZ_sub(s,IP)) ;
       d=SetUB(IP)-SetA(IP);
);

display P1_DA.l,P2_DA.l,PW_DA.l,z_MP.l,SetZ_MP,SetA,SetUB;

