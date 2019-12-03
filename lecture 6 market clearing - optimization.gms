
  sets
  i   generators /i1*i3/
  d   inelastic loads /d1*d3/
  n   buses /n1*n3/

  Mapi(i,n) /i1.n1 , i2.n2 , i3.n3/
  Mapd(d,n) /d1.n1 , d2.n2 , d3.n3/
  Mapnm(n,n) /n1.n2 , n2.n3 , n3.n1 , n2.n1 , n3.n2 , n1.n3/;
  alias (n,m);

  parameters
  P_max(i)   installed capacity /
  i1  20
  i2  10
  i3  25/

  c(i)   marginal cost /
  i1  16
  i2  19
  i3  15/

  L(d) Load level /
  d1 5
  d2 20
  d3 15/
  ;

  Table B(n,n)   Transmission lines susceptance
         n1         n2   n3
  n1     0          100  125
  n2     100        0    150
  n3     125        150  0;

  Table Fmax(n,n)   Transmission lines capacity
         n1         n2   n3
  n1     0          5    10
  n2     5          0    10
  n3     10         10   0;


  Free variables
  cost           Total expected system cost
  P(i)          DA dispatch of generators
  theta(n)    Voltage angles
;

  equations
 costfn,node_balance,Prod_max,Prod_min,flow_max,flow_min,slack_bus;

  costfn.. cost=e= sum(i,c(i)*P(i));

* DA constraints
  node_balance(n).. sum(i$Mapi(i,n),P(i))-sum(d$Mapd(d,n),L(d))-sum(m$Mapnm(n,m),B(n,m)*(theta(n)-theta(m)))=e=0;

  Prod_max(i).. P(i)=l=P_max(i);
  Prod_min(i).. P(i)=g=0;
  flow_max(n,m)$Mapnm(n,m).. B(n,m)*(theta(n)-theta(m))=l=Fmax(n,m);
  flow_min(n,m)$Mapnm(n,m).. B(n,m)*(theta(n)-theta(m))=g=-Fmax(n,m);
  slack_bus.. theta('n3')=e=0;


* KKT conditions


  model market /  all / ;
  solve market using lp minimizing cost;

  display
  cost.l, P.l;

