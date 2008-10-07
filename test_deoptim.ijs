NB. test math/deoptim addon

Note 'To run all tests:'
  load 'math/deoptim'
  load 'math/deoptim/test_deoptim'
)

NB. =========================================================
NB. create nouns and verbs for testing

Rosenbrock=: 3 : 0
  'x1 x2'=. 2{. y
  (100 * *:(x2 - *:x1)) + *:(1 - x1)
)

res=: deoptim_pdeoptim_ 'Rosenbrock_base_';(|:2 2$_10 10)

NB. polynomial test ref paper
NB. Chebychev polynomials
T4=: 1 _8 8
T8=: 1 _32 160 _256 128
T16=: 1 _128 2688 _21504 84480 _180224 212992 _131072 32768
NB.T4coeff=: 1 0 _8 0 8
NB.T8coeff=: 1 0 _32 0 160 0 _256 0 128
NB.T16coeff=: 1 0 _128 0 2688 0 _21504 0 84480 0 _180224 0 212992 0 _131072 0 32768

coeff2k=: [: }:@:, 0 ,.~ ] NB. }:,x,.0
Tcoeff=: coeff2k T8
Chchevpoly=: (coeff2k~ p. ]) f.  NB. eval chebychev polynomial

NB. eg: lowlimit Tcoeff
lowlimit=: p.&1.2  NB. calc lower limit at 1.2

NB. C code equivalent of  (|.coeff) p. val
p_C=: 4 : 0  
 coeff=. x
 val=. y
 px=. 0{coeff
 for_i. >:i. <:#coeff do.
   px=. (val * px) + i{coeff
 end.
 px
)

v1 =: 54.735 31.069 _79.616 _22.279 _76.028 _42.021 _34.925 _44.066 _53.706

xVect=: ((],-)1.2) ,(%~i:) (9<#Tcoeff){30 50

objfn=: 3 : 0
  res=. ([: +/ [: *: -. * 1 < |) 2}.y  NB. between _1 & 1
  res + ([:+/[:(*: * 0&>) (lowlimit Tcoeff) -~ ]) 2{. y NB. _1.2 & 1.2
)

Note 'test'
  2.1e_8>72.6606669-T8coeff p. 1.2     NB. t8(1.2)
  2.1e_8>10558.1450229-T16coeff p. 1.2 NB. t16(1.2)
  2.1e_8>72.6606669-t8 Chchevpoly 1.2     NB. t8(1.2)
  2.1e_8>10558.1450229-T16 Chchevpoly 1.2 NB. t16(1.2)
  0 = T16 objfn@:Chchevpoly xVect
  0 = T8  objfn@:Chchevpoly xVect
)

testfunc=: objfn@:p.&xVect
testfuncC=: objfn@:p_C&xVect

Note ''
  NB.vtr itermax npop f cr popln strategy refresh digits
 cntrl=: 1e_6;1000;10;0.8;0.9;'';1;10;4
 cntrl=: 1e_6;13000;10;0.7;0.5;'';2;10;4
 cntrl=: 1e_6;1000;10;0.8;0.9;'';3;10;4
 cntrl=: 1e_6;1000;10;0.5;0.9;'';4;10;4
 cntrl=: 0;5000;10;0.8;0.9;'';3;10;4
 res=: cntrl deoptim_pdeoptim_ 'testfunc_base_';(|:((#Tcoeff),2)$_1000 1000)
 resC=: cntrl deoptim_pdeoptim_ 'testfuncC_base_';(|:((#Tcoeff),2)$_1000 1000)
 plot _1.01 1.01;'((0;0){::res) p. y'
 plot _1.01 1.01;'((0;0){::resC) p_C y'
)
