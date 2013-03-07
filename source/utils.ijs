NB. =========================================================
NB. utils for Differential Evolution addon

NB.*getSampleNR v Sample from array without replacement
NB.-y: array to sample
NB.-x: optional number and size of samples
NB.-      defaults to 1 sample of same size as y
NB.-note: getSampleNR=: ((1 , #) $: ]) : (] {~ {:@[ ? {.@[ $ #@])
getSampleNR=: 3 : 0
  (1 , #y) getSampleNR y
:
  'num sz'=. x
  (sz ? num $ #y){y
)

NB.*getSampleR v Sample from array with replacement
NB.-y: array to sample
NB.-x: optional number and size of samples
NB.-     defaults to 1 sample of same size as y
NB.-note: getSampleR=: ((1 , #) $: ]) : (] {~ [ ?@$ #@])
getSampleR=: 3 : 0
  (1 , #y) getSampleR y
:
  (x ?@$ #y){y
)
