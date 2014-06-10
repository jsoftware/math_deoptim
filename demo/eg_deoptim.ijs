NB. =========================================================
NB. GUI for DE example
NB. http://www.icsi.berkeley.edu/~storn/code.html

require 'plot strings'
require 'math/deoptim'
coclass 'pegde'
coinsert 'pdeoptim jgl2'

steps=: {. + (1&{ - {.) * (i.@>: % ])@{:

NB. =========================================================
NB. Problem definition - Chebyshev polynomial fitting problem

NB. Coefficients for Chebychev polynomials
T4=: 1 _8 8
T8=: 1 _32 160 _256 128
T16=: 1 _128 2688 _21504 84480 _180224 212992 _131072 32768
FuncNames=: 'T4 T8 T16'

updateFunc=: 3 : 0
  Tcoeff=: }:,y,.0
  xVect=: ((],-)1.2) ,(%~i:) (9<#Tcoeff){30 50 NB. number of samples
  LowLimit=: Tcoeff p. 1.2
)

NB.*testfunc v Evaluates set of vars and calculates value
testfunc=: 3 : 0
  res=. y p. xVect
  objfn res
)

reportProgress=: 3 : 0
  'BestVars BestVal Generations Digits'=. 4{.y
  result=: pack 'BestVars BestVal Generations nFEval Digits'
  updateOutput''
)

NB.*objfn v Objective function for calcluting value of a set of vars
NB. value is sum of squared errors from 63 or 103 samples
NB. errors when value outside _1 & 1 where x is _1 to 1
NB. errors when value is less than lowlimit Tcoeff at _1.2 & 1.2
objfn=: 3 : 0
  value=. ([: +/ [: *: -. * 1 < |) 2}.y  NB. between _1 & 1
  value + ([:+/[:(*: * 0&>) LowLimit -~ ]) 2{. y NB. _1.2 & 1.2
)

NB.*optimize v Applies DE to parameters given & returns DE result
optimize=: 3 : 0
  vtr=. 1e_12
  genmax=. 2000
  refresh=. 30
  args=. (pack 'vtr genmax refresh') pset y
  cntrl=. 'vtr genmax npop f cr strategy refresh' psel args
  bounds=. |:((#Tcoeff),2)$(- , ]) 'bounds' pget args
  cntrl getDEoptim 'testfunc';bounds
)

NB. =========================================================
NB. Form Definition
EGDE=: 0 : 0
pc egde;pn "Differential Evolution Demo";
bin vh;
cc cblFunc combolist;
cc cblStrat combolist;
set cblStrat items DEBest1 DERand1 DERandToBest1 DEBest2;
set cblStrat select 2;
bin s;
cc bStart button;cn "Start";
cc bClear button;cn "Clear";
bin zv;

groupbox Parameters;
bin hv;
cc  lblBounds static center;cn "Bounds:";
cc  edBounds edit readonly center;
set edBounds text 100;
cc  slBounds slider 1 1 3 20 491 100;
bin zv;
cc  lblF static center;cn "F:";
cc  edF edit readonly center;
set edF text 0.8;
cc  slF slider 1 0 2 10 201 81;
bin zv;
cc  lblCR static center;cn "CR:";
cc  edCR edit readonly center;
set edCR text 0.9;
cc  slCR slider 1 0 2 10 101 91;
bin zv;
cc  lblNP static center;cn "NP:";
cc  edNP edit readonly center;
set edNP text 5;
cc  slNP slider 1 5 1 1 10 1;
bin z;
groupboxend;
bin s1z;

groupbox Output;
bin vhv;
cc lblIter static center;cn "Generations:";
cc valIter static center;
bin zv;
cc lblnEval static center;cn "Evaluations:";
cc valnEval static center;
bin zv;
cc lblVal static center;cn "Cost-Value:";
cc valVal static center;
bin zzs1v;
cc lblMem static;cn "Best member:";
cc valMem static;
bin z;
groupboxend;
bin z;

minwh 247 256;cc deplot isigraph;
set deplot stretch 1;
pas 4 4;
)

NB. ---------------------------------------------------------
create=: 3 : 0
  wd EGDE
  initControls''
  deplot=: ('egde';'deplot') conew 'jzplot'
  updateOutput''
  wd 'pshow;'
)

destroy=: 3 : 0
  wd'pclose'
  destroy__deplot''
  codestroy''
)

NB. ---------------------------------------------------------
egde_close=: destroy

egde_bStart_button=: 3 : 0
  'bounds npop f cr strategy'=. 0". > edBounds;edNP;edF;edCR;cblStrat_select
  func=. cblFunc
  strategy=. >:strategy
  result=: optimize pack 'bounds npop f cr strategy func'
  updateOutput''
)

egde_bClear_button=: 3 : 0
  erase 'result'
  updateOutput''
)

egde_cblFunc_select=: 3 : 0
  updateFunc cblFunc~
  egde_bClear_button''
)

egde_slBounds_changed=: 3 : 0
  wd 'set edBounds text *',slBounds
)

egde_slF_changed=: 3 : 0
  wd 'set edF text *',": 0 >. 100 %~ <: 0".slF
)

egde_slCR_changed=: 3 : 0
  wd 'set edCR text *',": 0 >. 100 %~ <: 0".slCR
)

egde_slNP_changed=: 3 : 0
  wd 'set edNP text *',": 0".slNP
)

egde_deplot_update=: 3 : 0
  pd__deplot 'reset'
  pd__deplot 'xrange _1.5 1.5;yrange _1 10;grids 0 0;'
  pd__deplot 'color red;'
  pd__deplot _2 _1.2 _1.2 1.2 1.2 2; 2|.(4#LowLimit), _1 _1
  pd__deplot _1 _1 1 1 ; 1000 1 1 1000
  if. 0 ~: 'Generations' pget result do.
    pd__deplot 'color blue;'
    vls=. steps _1.4 1.4 100
    pd__deplot vls;('BestVars' pget result) p. vls
  end.
  glpaint__deplot ''
)

NB. ---------------------------------------------------------
updateOutput=: 3 : 0
  if. 0~: 4!:0 <'result' do.
    result=: (;:'Generations nFEval BestVal BestVars'),.<0
  end.
  wd 'set valIter  text *', ": 'Generations' pget result
  wd 'set valnEval text *', ":      'nFEval' pget result
  wd 'set valVal   text *', ":     'BestVal' pget result
  wd 'set valMem   text *', ' 'joinstring '0.3' 8!:0 'BestVars' pget result
  wd 'msgs'
  egde_deplot_update''
)

initControls=: 3 : 0
  wd 'set cblFunc items ',FuncNames
  wd 'set cblFunc select ',": (;:FuncNames) i. <'T8'
  updateFunc T8
  ''
)

''conew 'pegde'
