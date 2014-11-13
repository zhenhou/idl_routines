;;;
; NAME: kplot_hist
; PURPOSE:
;   Easily plot histogram of runlist cuts
;
; MODIFICATION HISTORY:
;  08/20/2011: (KTS) Created
;;;

;...................................................................
; Main function
;
pro kplot_hist, data, ktitle, nbins=nbins, xmax=xmax
compile_opt IDL2, HIDDEN

if (~keyword_set(nbins)) then nbins = 50

hmin=min(data) & hmax=max(data)
hh=  histogram(data, MIN=hmin, MAX=hmax, nbins=nbins)
hh_norm = DOUBLE(hh) / max(hh)
xvec = hmin + findgen(nbins) * (hmax-hmin) / double(nbins)
low_medwt = median(data) / 1e6 & high_medwt = median(data) * 2.

if ~keyword_set(xmax) then xmax = max(xvec)
;xmax = 20;0.002
;if (max(xvec) lt 0.002) then xmax = max(xvec)

plot, xvec, hh_norm, title=ktitle, xrange=[0,xmax]
;oplot, [low_medwt, low_medwt], [0,1], color=!red
;oplot, [high_medwt, high_medwt], [0,1], color=!red
end

