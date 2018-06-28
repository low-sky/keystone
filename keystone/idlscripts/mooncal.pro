;+
;  Compute gain correction factors based on OnMoon, OffMoon
;  Observations
; usage:
;    moonCal, <onMoonScan>, <offMoonScan>, <factors>, <tCmb>, <doOnOff>
;      <onMoonScan>   On center of the moon observation
;      <offMoonScan>  Off center of the moon observation
;      <factors>      Scale factors to multiply data for 
;                     acurate calibration
;      <tCmb>         optionally enter CMB temp (default 2.725 K)
;      <doOnOff>      optionally (>0) calibrate using On-Off, not sig-ref/ref
;Example:
;filein,'/home/scratch/glangsto/moon/AGBT08C_062_03.raw.acs.fits'
; pair of Moon scans
;onMoonScan  = 220
;offMoonScan = 221
; now get moon parameters
;-
;IDL procedure for KFPA computing Moon Calibration factors
;HISTORY 
;12JAN06 GIL optionally use last int to solve problem of moving scan
;11APR06 GIL create a single formated string for plotting
;11MAR25 GIL modularize
;11MAR24 GIL use standard moon model
;11FEB21 GIL estimate errors in opacity, tau
;11FEB10 GIL add moon phase code
;11JAN21 GIL work on Moon data
;11JAN02 GIL out of band W51
;10NOV30 GIL out of band W51
;10SEP01 GIL 800 MHz test
;10AUG19 GIL initial KFPA re-commissioning session
;10MAY02 GIL initial version
;10MAR07 GIL check K band spectra
;08AUG13 GIL plot the test NOD scans

pro moonCal, onMoonScan, offMoonScan, factors, tCmb, doOnOff=doOnOff, doLastInt=doLastInt, ifnum=ifnum

  if n_elements(ifnum) eq 0 then ifnum = 0 

  if (not keyword_set( offMoonScan)) then begin
      usage,'moonCal',/verbose
      return
  endif

  ; if no output factor array, create
  if (not keyword_set( factors)) then   factors = dindgen(16)
  if (not keyword_set( tCmb)) then   tCmb = 2.725
  if (not keyword_set( doOnOff)) then   doOnOff = 0
  if (not keyword_set( doLastInt)) then   doLastInt = 1

 
  unfreeze
  gettp,onMoonScan,plnum=1, ifnum=ifnum
  elMoon = !g.s[0].elevation
  endfix & show
  freeze
  gettp,offMoonScan,plnum=1, ifnum=ifnum
  endfix & oshow
  elRef = !g.s[0].elevation
  unfreeze
  onInfo = find_scan_info(onMoonScan)
  offInfo = find_scan_info(offMoonScan)

  dc = !g.s[0]
  mjd = dc.MJD
  projId = strtrim(string(dc.projid),2)
  print,dc.timestamp," -> ", mjd
  tAmb = dc.tambient - 273.15 ; convert to 'C'
  freqMHz= dc.reference_frequency/1.E6
  freqGHz= dc.reference_frequency/1.E9

; get opacity and estimate error, by differencing tau at +/- 1 hour
  hour = 1./12.
  tauM1 = gettau(mjd-hour, freqMHz)
  tau = gettau(mjd, freqMHz)
  tauP1 = gettau(mjd+hour, freqMHz)
  dTau = tau - tauM1
  if (dTau < 0.) then dTau = -1.*dTau
  dTau2 = tau - tauP1
  if (dTau2 < 0.) then dTau2 = -1. * dTau2
;use larger of two values as error estimate
  if (dTau2 > dTau) then dTau = dTau2
  print, tau, dTau, format='(%"zenithOpacity    : %7.3f +/- %7.3f")'
  opacityFactor=1.0D0 & opacityFactor2 = 1.0D0
  tempAir = 1.0D0 & tempAir2 = 1.0D0
  etaA = 1.D0 & etaB = 1.D0
  etaGBT, freqMHz, etaA, etaB

; compute factors for on moon
  tsystau, tau, opacityFactor, tAmb, elMoon, freqGHz, tempAir
  tsystau, tau+dTau, opacityFactor2, tAmb, elMoon, freqGHz, tempAir2
  dOpacity = opacityFactor2-opacityFactor
  print,opacityFactor, dOpacity, format='(%"Opacity Factor(M): %7.3f+/-%7.3f")'
  moonFactor = opacityFactor & dMoonFactor = dOpacity
; compute factors for on Ref
  tsystau, tau,      opacityFactor, tAmb, elRef, freqGHz, tempAir
  tau2 = tau+dTau
  tsystau, tau2, opacityFactor2, tAmb, elRef, freqGHz, tempAir2
  dOpacity = opacityFactor2-opacityFactor
  if dOpacity < 0. then dOpacity = - dOpacity
  print,opacityFactor, dOpacity, format='(%"Opacity Factor(R): %7.3f+/-%7.3f")'
  refFactor = opacityFactor & dRefFactor = dOpacity

  natm, elMoon, na
  print, 'Freq, etaA, etaB : ',freqMHz, etaA, etaB
  moonFactor = 1.0D0 & refFactor = 1.0
  gainel, elMoon, moonFactor
  gainel, elRef, refFactor
;  print, 'Gain Factor (m,r): ',moonFactor, refFactor

;get the moon model temperature
  moonModel, mjd, freqGHz, tMoon, tErrMoon, /doprint

  etal = 0.99

  natm, elMoon, na 
  opacityMoon = 1.0D0 & tempAirMoon = 1.0D0
  opacityRef = 1.0D0 & tempAirRef = 1.0D0
  tsystau, tau, opacityMoon, tAmb, elMoon, freqGHz, tempAirMoon
  tsystau, tau, opacityRef, tAmb, elRef, freqGHz, tempAirRef
  opacityTotal = na*tau
  print,'Opacity factor :  ', 1./opacityMoon
  allFactors = (etal*etaB*moonFactor)/opacityMoon
  print,'All Corrections:  ', allFactors

  tOne = tMoon*allFactors
  tTwo = (tCmb*etal*etaB/opacityRef)
  tThree = (tempAirMoon*etal*((1./opacityMoon)-(1./opacityRef)))

  tMoonObs = (tOne-tTwo)+tThree

  print,'Predicted On-Off Moon Temp K: ',tMoonObs
  print,'Predicted On Moon Temp(-Trx): ',tOne
  print,'Predicted Off Moon CMB Temp : ',tTwo
  print,'Predicted Elevation Change K: ',tThree

  si = scan_info(onMoonScan, /quiet)
  nfeed = si.n_feeds
  npol = si.n_polarizations
  nif = si.n_ifs
  nsamp = si.n_samplers

  onMoonTsys = dblarr(nif, nfeed, npol)
  offMoonTsys = dblarr(nif, nfeed, npol)
  unfreeze
  
  dTs = dindgen(nif, nfeed, npol)
  unfreeze
  showsigref,onMoonScan, offMoonScan, dTs, freqObs=freqObs

  ; if doing On source - Off source calibration ont (sig-ref)/ref
  if (doOnOff) then begin   
     showtsys,onMoonScan,onMoonTsys
     showtsys,offMoonScan,offMoonTsys
     for iii = 0,13 do begin 
        dTs[iii]=onMoonTsys[iii]-offMoonTsys[iii] 
        if (dTs[iii] lt 1) then dTs[iii] = 1.
     endfor
  endif  ; end if On - Off Cal

  factors = tMoonObs[0] / dTs
  if total(finite(factors,/infinity)) gt 0 then begin
     factors[where(finite(factors,/infinity))]=!values.f_nan
  endif

; EWR: Why is this here??
;  freexy
;  if (doOnOff) then freey
;  getsigref,onMoonScan,OffMoonScan,plnum=0
;  firstPol = !g.s[0].polarization
;  endfix & show
;  freeze
;  getsigref,onMoonScan,OffMoonScan,plnum=1
;  secondPol = !g.s[0].polarization
;  endfix & oshow
;  unfreeze

; keep values in range for the printout
  for iii = 0,13 do begin 
    if (factors[iii] lt .01) then factors[iii] = .01
    if (factors[iii] gt 9.999) then factors[iii] = 9.999
  endfor

  ;; print,"Calibration Factors"
  ;; print,format='($,"--gain-factors-left  ")'
  ;; for iii = 0, 11,2 do begin print,factors[iii],format='($,f5.3,",")'& endfor
  ;; for iii = 12,13,2 do begin print,factors[iii],format='($,f5.3," ")'& endfor
  ;; print,""
  ;; print,format='($,"--gain-factors-right ")'
  ;; for iii = 1, 11,2 do begin print,factors[iii],format='($,f5.3,",")'& endfor
  ;; for iii = 13,13,2 do begin print,factors[iii],format='($,f5.3," ")'& endfor
  ;; print,""
  moonTempLogFile = 'MoonTemp.log'
  print,'Summary in :',moonTempLogFile
  openw,alun,moonTempLogFile,/get_lun,/append
  printf,alun,mjd,freqGHz,tMoon,tMoonObs, format='($,f10.3,f7.3,f6.1,f6.1," ")'
  printf,alun,'----------------------'
  for k=0,(npol-1) do begin
     for i=0,(nif-1) do begin
        for j=0,(nfeed-1) do begin
           printf, alun, 'IF:'+strcompress(string(i),/rem)+$
                   " FD:"+strcompress(string(j),/rem)+$
                   " PL:"+strcompress(string(k),/rem)+$
                   " Freq:"+strcompress(string(freqObs[i,j,k]),/rem)+$
                   " Factor:"+strcompress(string(factors[i,j,k]),/rem)
        endfor
     endfor
  endfor
  Printf,alun,'----------------------'
  printf,alun,opacityTotal,elMoon,projId,onMoonScan,offMoonScan,$\
        format='($," ",f5.3,f5.1," ",a,i4,i4)'
  printf,alun,""
  close,alun
  free_lun,alun
  return
end

