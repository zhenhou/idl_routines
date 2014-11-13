;;;
; NAME: save_mode_coupling_lps12.pro
; PURPOSE:
;   Make a sav file with mode coupling information
;
; INPUTS:
;   stub,        must match that used in from try_end2end.pro
;   
; OUTPUTS:
;
; NOTES:
;
; MODIFICATION HISTORY:
;  04/25/2012: (KTS) Created from /home/rkeisler/ps09/save_mode_coupling.pro
;;;
pro save_mode_coupling_lps12, stub, idx_list=idx_list

if ~keyword_set(idx_list) then $
  idx_list=[-1, 0] ;,1,3,4,5] ; Temporary
  ;idx_list=[0,1] ; Temporary

field_list = indgen(20)
nfields = n_elements(field_list)

f = lps12_fieldstruct()
nidx = n_elements(idx_list)
colors = kscolor_array()

end2end_dir = '/home/kstory/lps12/end2end/'
out_dir     = '/home/kstory/lps12/twod_tfs/'
savfile = out_dir+'save_mode_coupling'+stub+'.sav'

; add to the sav file
if (file_test(savfile) eq 1) then restore, savfile

for idx=0,nfields-1 do begin

    ; skip fields not in idx_list
    if intersect(idx, idx_list) eq -1 then continue
    f = lps12_fieldstruct() ; temporary
    fst = f[idx]
    print, 'analyze field ', idx, ' '+fst.name

    ;restore,end2end_dir+'end_'+fst.name+'_'+stub+'_kweight.sav'
    restore,end2end_dir+'end_'+fst.name+'_'+stub+'.sav'
    mode_coupling = smooth(transfer_iter[*,0,4]/transfer_iter[*,0,0],3)
    l_mode_coupling = ellkern
    if idx eq 0 then begin
        nl = n_elements(l_mode_coupling)
        mode_couplings = fltarr(nfields,nl)
        plot,l_mode_coupling,mode_coupling, $
          yr=[-1,1]*0.05+1,/yst,xr=[0,3e3], /xst, $
          xtitle='!12l!X!N',ytitle='(TF final)/(TF initial)', $
          chars=1.8,/nodata, $
          title='Effect of Mode-Coupling Correction'
        oplot,[-1,1]*999999.,[1,1],lines=2
          
    endif
    mode_couplings[idx,*] = mode_coupling
    oplot,l_mode_coupling,mode_couplings[idx,*],color=colors[idx]
endfor
legend,field_vec(),/top,/right,textc=colors
save,l_mode_coupling,mode_couplings,filename=savfile
end


