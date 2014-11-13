;;;
; NAME: script13_0409
;
; NOTES:
;  1) plot s12 and k11 lcdm constraints
;;;


PRO plot_k11
restore, '/home/kstory/lps12/end2end/run_09/combined_spectrum_ApJ2012.sav'
sk = dl_all
erk = diag_nobeam

restore, '/home/rkeisler/ps09/combined_spectrum_20110110_185553_2008kmask.sav'
sr = dl_all
err = diag_nobeam

vec=indgen(47)+9

ratio = (sk/sr)
ratio /= mean(ratio[vec])
plot, l[vec], ratio[vec],yr=[0.8,1.2]
errplot, l[vec], (ratio-erk/sk)[vec], (ratio+erk/sk)[vec]

er2use = erk/sk

chisq = total( (((1-ratio)/er2use)[vec])^2.)

stop
END

