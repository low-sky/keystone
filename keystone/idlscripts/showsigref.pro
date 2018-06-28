;+
; Show all of the Tsource values associated with a scan,refscan pair.
;
; showsigref prints (signal-ref)/ref calibrated values for all samplers'
; usage:
;   showsigref, <sigScan>,<refScan>,<TAValues>[,/bysampler]
;        where  <sigScan>   signal scan of a target source
;               <refScan>   referenc scan for position switched cal.
;               <TAValues>  Output Antenna Temps
;               <bySampler> Optionally list values by sampler not beam
;      Produces a table of Ta calibrated source temperature values
; <p>
; This uses getsigref to retrieve the system temperature for each sampler
; in a scan.  If the bysampler keyword is set then they are displayed
; under each sampler name with at most 8 samplers per row.  If that
; keyword is not set then they are displayed by feed and polarization
; with one row for each spectral window (ifnum).  Any triad of
; spectral window, feed, and polarization that doesn't have any data
; associated with is left blank in the output.
;
; <p>
; Because this must retrieve the raw data to calculated Tsys it often
; takes time to finish.
;
; @param scan {in}{required}{type=integer} The scan number of
; interest.
; @keyword bysampler {in}{optional}{type=boolean} When set the system
; temperature values will be identified by sampler name.
;-
; HISTORY
; 2011Feb14 GIL Initial version based on Bob Garwood's showtsys.pro
;
; @version $Id$
;
pro showsigref, scan, refscan, tAValues, bysampler=bysampler, doOneInt=doOneInt
    compile_opt idl2

    if n_elements(refscan) eq 0 then begin
        usage,'showsigref',/verbose
        return
    endif
    if (not keyword_set(doOneInt)) then doOneInt = 0
    si = scan_info(scan,/quiet)

    if size(si,/type) ne 8 then begin
        print,'Scan not found'
        return
    endif

    if n_elements(tAValues) ne 16 then tAValues = dindgen(16)

    if si.n_cal_states ne 2 then begin
        print,'Requires 2 cal states to determine Tsys'
        return
    endif

    nfeed = si.n_feeds
    npol = si.n_polarizations
    nif = si.n_ifs
    nsamp = si.n_samplers
    someave = 1.0d0
    somerms = 1.0d0

    isFrozen = !g.frozen
    freeze

    nOut = 0                    ; prepare to fill the output array
    if keyword_set(bysampler) then begin
                                ; awkward, may be badly sorted due to
                                ; things like A9 sorting after A10
        ; sort on ports only
        ports = fix(strmid(si.samplers,1))
        sortedSamplers = si.samplers[sort(ports)]
        lines = ceil(nsamp/8.0)
        for i=0,lines-1 do begin
            first = i*8
            last = min([first+7,nsamp-1])
            for k=first,last do begin
                print,sortedSamplers[k],format='($,"  ",A3," ")'
            endfor
            print
            for k=first,last do begin
                getsigref,scan,refscan,sampler=sortedSamplers[k],/q
                nChan = n_elements(*!g.s[0].data_ptr)
                somedata = (*!g.s[0].data_ptr)[(nChan/5):(nChan*4/5)]
                linejrect, somedata, someave, somerms
;                someave  = avg(somedata)
;                somerms  = stddev(somedata)
                !g.s[0].tsys=someave
                tAValues[nOut] = !g.s[0].tsys
                nOut = nOut+1
                print,!g.s[0].tsys,format='($,f5.1," ")'
            endfor
            print
            print
        endfor
    endif else begin
        ; label line first
        print,format='($,"#IF")
        for j=0,(nfeed-1) do begin
            for k=0,(npol-1) do begin
                print,si.feeds[j],strmid(si.polarizations[k],0,1),format='($,i2,a1,"   ")'
            endfor
        endfor
        ; print label once
        print,"     MHz     Az(d)  El(d) Scan"

        for i=0,(nif-1) do begin
            print,i,format='($,i2," ")'
            for j=0,(nfeed-1) do begin
                for k=0,(npol-1) do begin
                    select,count,scan=scan,ifnum=i,fdnum=j,plnum=k,int=iNum,/quiet
                    if count ne 0 then begin
                        if (doOneInt) then getsigref,scan,refscan,ifnum=i,fdnum=j,plnum=k,int=doOneInt,/quiet
                        if (not doOneInt) then getsigref,scan,refscan,ifnum=i,fdnum=j,plnum=k,/quiet
                        nChan = n_elements(*!g.s[0].data_ptr)

                        somedata = (*!g.s[0].data_ptr)[(nChan/5):(nChan*4/5)]
                        linereject, somedata, someave,somerms
;                        someave  = avg(somedata)
;                        somerms  = stddev(somedata)
                        !g.s[0].tsys=someave
                        tAValues[nOut] = !g.s[0].tsys
                        nOut = nOut+1
                        print,!g.s[0].tsys,format='($,f5.1," ")'
                    endif else begin
                        print,format='($,"      ")'
                    endelse
                endfor
            endfor
            print,!g.s[0].center_frequency*1.E-6,format='($,f11.3)'
            print,!g.s[0].azimuth, !g.s[0].elevation,scan,$\
                  format='($,f7.1," ",f6.1,i5)'
            print
        endfor
    endelse
    if not isFrozen then unfreeze
end
 
