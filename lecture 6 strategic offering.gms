
  sets
  i   generators /i1*i3/
  d   inelastic loads /d1*d3/
  n   buses /n1*n3/
  l   lines /l1*l3/
;

  parameters

  M0 / 10000 /

  P_max(i)   installed capacity /
  i1  20
  i2  10
  i3  25/

  c(i)   marginal cost /
  i1  16
  i2  19
  i3  15/

  Load(d) Load level /
  d1 5
  d2 20
  d3 15/

  B(l)   Transmission lines susceptance /
  l1 100
  l2 125
  l3 150/
;

  parameters

  Fmax(l)   Transmission lines capacity /
  l1 5
  l2 10
  l3 10/
  ;


  Free variables
  cost           Total expected system cost
  P(i)           DA dispatch of generators
  theta(n)       Voltage angles
  lambda(n)      LMPs
  gamma          node reference N3 dual variable
;

  Positive variables
  offer         price offer
  mu_max(i)     max production dual variables
  mu_min(i)     min production dual variables
  rho_max(l)    max flow dual variables
  rho_min(l)    min flow dual variables
;

  Integer variables

  u_mu_max(i)     max production binary variables
  u_mu_min(i)     min production binary variables
  u_rho_max(l)    max flow binary variables
  u_rho_min(l)    min flow binary variables
;

  equations
 costfn
 offer_max
 node_balance_1,node_balance_2,node_balance_3
 Prod_max,Prod_min
 flow_max_1,flow_max_2,flow_max_3,flow_min_1,flow_min_2,flow_min_3
 slack_bus
 stat_g1,stat_g2,stat_g3,stat_theta1,stat_theta2,stat_theta3
 comp_gmax_1,comp_gmax_2
 comp_gmin_1,comp_gmin_2
 comp_fmax_1,comp_fmax_21,comp_fmax_22,comp_fmax_23
 comp_fmin_1,comp_fmin_21,comp_fmin_22,comp_fmin_23
;

  costfn.. cost =e= c('i1')*P('i1')
                         - (-c('i2')*P('i2')-c('i3')*P('i3')-P_max('i2')*mu_max('i2')-P_max('i3')*mu_max('i3')
                                 +(Load('d1')*lambda('n1')+Load('d2')*lambda('n2')+Load('d3')*lambda('n3'))
                                 -sum(l,Fmax(l)*rho_max(l))-sum(l,Fmax(l)*rho_min(l)));
* UL constraint
  offer_max.. offer =l= 20;

* DA constraints
  node_balance_1.. P('i1') + B('l1')*(theta('n2')-theta('n1')) + B('l2')*(theta('n3')-theta('n1')) =e= Load('d1');
  node_balance_2.. P('i2') + B('l1')*(theta('n1')-theta('n2')) + B('l3')*(theta('n3')-theta('n2')) =e= Load('d2');
  node_balance_3.. P('i3') + B('l2')*(theta('n1')-theta('n3')) + B('l3')*(theta('n2')-theta('n3')) =e= Load('d3');
  slack_bus.. theta('n3')=e=0;
  Prod_max(i).. P(i)-P_max(i)=l=0;
  Prod_min(i).. -P(i)=l=0;
  flow_max_1.. B('l1')*(theta('n1')-theta('n2'))-Fmax('l1')=l=0;
  flow_min_1.. -B('l1')*(theta('n1')-theta('n2'))-Fmax('l1')=l=0;
  flow_max_2.. B('l2')*(theta('n3')-theta('n1'))-Fmax('l2')=l=0;
  flow_min_2.. -B('l2')*(theta('n3')-theta('n1'))-Fmax('l2')=l=0;
  flow_max_3.. B('l3')*(theta('n2')-theta('n3'))-Fmax('l3')=l=0;
  flow_min_3.. -B('l3')*(theta('n2')-theta('n3'))-Fmax('l3')=l=0;

* KKT conditions

  stat_g1..  offer - lambda('n1') + mu_max('i1') - mu_min('i1') =e=0;
  stat_g2..  c('i2') - lambda('n2') + mu_max('i2') - mu_min('i2') =e=0;
  stat_g3..  c('i3') - lambda('n3') + mu_max('i3') - mu_min('i3')  =e=0;

  stat_theta1.. B('l1')*(lambda('n1')-lambda('n2')+rho_max('l1')-rho_min('l1'))
                 + B('l2')*(lambda('n1')-lambda('n3')-rho_max('l2')+rho_min('l2')) =e= 0;
  stat_theta2.. B('l1')*(lambda('n2')-lambda('n1')-rho_max('l1')+rho_min('l1'))
                 + B('l3')*(lambda('n2')-lambda('n3')+rho_max('l3')-rho_min('l3')) =e= 0;
  stat_theta3.. B('l2')*(lambda('n3')-lambda('n1')+rho_max('l2')-rho_min('l2'))
                 + B('l3')*(lambda('n3')-lambda('n2')-rho_max('l3')+rho_min('l3')) + gamma =e= 0;

  comp_gmax_1(i).. mu_max(i) =l= M0*u_mu_max(i);
  comp_gmax_2(i).. P_max(i)-P(i)=l= M0*(1-u_mu_max(i));

  comp_gmin_1(i).. mu_min(i) =l= M0*u_mu_min(i);
  comp_gmin_2(i).. P(i)=l= M0*(1-u_mu_min(i));

  comp_fmax_1(l).. rho_max(l) =l= M0*u_rho_max(l);
  comp_fmax_21.. Fmax('l1') - B('l1')*(theta('n1')-theta('n2')) =l= M0*(1-u_rho_max('l1'));
  comp_fmax_22.. Fmax('l2') - B('l2')*(theta('n3')-theta('n1')) =l= M0*(1-u_rho_max('l2'));
  comp_fmax_23.. Fmax('l3') - B('l3')*(theta('n2')-theta('n3')) =l= M0*(1-u_rho_max('l3'));

  comp_fmin_1(l).. rho_min(l) =l= M0*u_rho_min(l);
  comp_fmin_21.. Fmax('l1') + B('l1')*(theta('n1')-theta('n2')) =l= M0*(1-u_rho_min('l1'));
  comp_fmin_22.. Fmax('l2') + B('l2')*(theta('n3')-theta('n1')) =l= M0*(1-u_rho_min('l2'));
  comp_fmin_23.. Fmax('l3') + B('l3')*(theta('n2')-theta('n3')) =l= M0*(1-u_rho_min('l3'));

  model market /  all / ;
  solve market using mip minimizing cost;

  display
  cost.l, P.l , offer.l , lambda.l;

