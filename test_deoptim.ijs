NB. test math/deoptim addon

Note 'To run all tests:'
  load 'math/deoptim'
  load 'math/deoptim/test_deoptim'
)

NB. =========================================================
NB. create nouns and verbs for testing

NB. Rosenbrock function
Rosenbrock=: 3 : 0
  'x1 x2'=. 2{. y
  (100 * *:(x2 - *:x1)) + *:(1 - x1)
)

NB. Fitting Chebychev polynomial problem
T8=: 1 0 _32 0 160 0 _256 0 128
lowerlimit=: T8 p. 1.2        NB. calc lower limit at 1.2
xVect=: ((],-)1.2) ,(%~i:) 30 NB. points to sample at

objfn=: 3 : 0
  res=. ([: +/ [: *: -. * 1 < |) 2}.y  NB. between _1 & 1
  res + ([:+/[:(*: * 0&>) lowerlimit -~ ]) 2{. y NB. _1.2 & 1.2
)

NB. ChebchevPoly v Verb for evaluating set of possible coefficients
ChebchevPoly=: objfn@:p.&xVect

Note 'test'
  2.1e_8>72.6606669-T8coeff p. 1.2     NB. t8(1.2)

)

Control=: 1e_6;1000;10;0.8;0.9;'';1;10;4

Note ''
  NB.vtr genmax npop f cr popln strategy refresh digits
 cntrl=: 1e_6;1000;10;0.8;0.9;'';1;10;4
 cntrl=: 1e_6;13000;10;0.7;0.5;'';2;10;4
 cntrl=: 1e_6;1000;10;0.8;0.9;'';3;10;4
 cntrl=: 1e_6;1000;10;0.5;0.9;'';4;10;4
 cntrl=: 0;5000;10;0.8;0.9;'';3;10;4
 res=: cntrl deoptim 'ChebchevPoly_base_';(|:((#T8),2)$_1000 1000)
 plot _1.01 1.01;'(0{::res) p. y'
)

test=: 3 : 0
  tmp=. deoptim 'Rosenbrock_base_';|:2 2$_10 10
  tmp=. Control deoptim 'Rosenbrock_base_';|:2 2$_10 10
  tmp=. getDEoptim 'Rosenbrock_base_';|:2 2$_10 10
  tmp=. (('strategy';2),,:'genmax';200)getDEoptim 'Rosenbrock_base_';|:2 2$_10 10

  assert. 0 = ChebchevPoly T8
  tmp=. Control deoptim 'ChebchevPoly_base_';|:9 2$_1000 1000
  tmp=. Control getDEoptim 'ChebchevPoly_base_';|:9 2$_1000 1000
  'test_deoptim passed'
)

smoutput test''