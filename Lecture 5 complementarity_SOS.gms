
sets
u /u1 * u2 /

variables
x
y
z

positive variables
m1
m2
m3
m4

sos1 variables
c1(u)
c2(u)
c3(u)
c4(u)



equations
con1
con2
con3
con4
con5
lagr
comp1
comp11
comp2
comp22
comp3
comp33
comp4
comp44
obje
;

con1.. -x+2=l=0 ;
con2.. -x-y+7=l=0;
con3.. y-9=l=0    ;
con4.. x-2*y-1=l=0 ;
con5.. 6*x+y-45=l=0 ;
lagr.. 1-m1+m2-2*m3+m4=e=0;
comp1.. (x+y-7) + m1 =e= c1('u1') + c1('u2');
comp11.. (x+y-7) - m1 =e= c1('u1') - c1('u2');
comp2.. (9-y) + m2 =e= c2('u1') + c2('u2')     ;
comp22.. (9-y) - m2 =e= c2('u1') - c2('u2')     ;
comp3.. (-x+2*y+1) + m3 =e= c3('u1') + c3('u2')   ;
comp33.. (-x+2*y+1) - m3 =e= c3('u1') - c3('u2')   ;
comp4.. (-6*x-y+45) + m4 =e= c4('u1') + c4('u2')    ;
comp44.. (-6*x-y+45) - m4 =e= c4('u1') - c4('u2')    ;
obje.. z=e=-y ;

model monday /all/;
solve monday using mip minimizing z;


