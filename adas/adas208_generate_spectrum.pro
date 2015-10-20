FUNCTION compute_lambda_breaks,adf04,lmax=lmax,peclimit=peclimit,amin=amin,debug=debug
	IF NOT keyword_set(peclimit) THEN peclimit=48
	IF NOT keyword_set(lmax) THEN lmax=1000.0
	IF NOT keyword_set(amin) THEN amin=1.0e8
	read_adf04,file=adf04,fulldata=data
	ntrn=n(data.upper)+1
	lambda=dblarr(ntrn)
	FOR i=0,ntrn-1 DO lambda[i]=cm2ang(data.wa[data.upper[i]-1]-data.wa[data.lower[i]-1])
	tmp=where(lambda LT lmax AND data.aval GT amin)
	IF tmp[0] EQ -1 THEN RETURN,-1
	lambda=lambda[tmp]			;find breaks on for lambda < lmax
	order=sort(lambda)
	lambda=lambda[order]
	nlam=n(lambda)+1
	IF nlam LT peclimit THEN RETURN,[0,max(lambda)]	;if < # transition limit (why there is one is beyond me...)
        nbrk=fix(nlam/peclimit)
	lbrk=dblarr(2,nbrk+1)
	lbrk[0,0]=0.0
	lbrk[1,nbrk]=max(lambda)
	FOR i=0,nbrk-1 DO BEGIN
		lbrk[1,i]=lambda[(i+1)*peclimit]			;looks like run_adas208 does GE WMIN and LT wmax
		lbrk[0,i+1]=lambda[(i+1)*peclimit]
	ENDFOR
	FOR i=0,nbrk DO IF lbrk[0,i] EQ lbrk[1,i] THEN lbrk[1,i]*=1.00001	;will crash if grouping of degenerate wavelengths
	IF keyword_set(debug) THEN stop
	RETURN,lbrk
END

PRO generate_208_pecs,adf04,telimit=telimit,nfiles=nfiles,lmax=lmax,amin=amin,pass_dir=pass_dir,file=file,rec=rec
	nfiles=-1
	IF NOT keyword_set(telimit) THEN telimit=35 	;number of te limit (why there is one is beyond me...)
	IF NOT keyword_set(pass_dir) THEN pass_dir='/home/reinke/adas/pass/'
	file=last(strsplit(adf04,'/',/extract))
	file=strsplit(file,'.',/extract)
	dens=[1.0e12,1.0e13,1.0e14,1.0e15]
	te=10^make(2.0,4.0,35)	
	lbrk=compute_lambda_breaks(adf04,lmax=lmax,amin=amin)
	IF lbrk[0] EQ -1 THEN RETURN
	nfiles=n(lbrk[0,*])+1
	print, adf04
	print, amin,lmax
	FOR i=0,n(lbrk[0,*]) DO BEGIN
		print, 'running adas208 for wavelength band '+strtrim(i,1)+' of '+strtrim(n(lbrk[0,*]),1)+': '+string(lbrk[0,i])+' to '+string(lbrk[1,i])
		pop=0
		run_adas208,adf04=adf04,te=te,dens=dens,/pec,pass_dir=pass_dir,wmin=lbrk[0,i],wmax=lbrk[1,i],pop=pop,amin=amin,rec=rec
		wait,2.0
		spawn,'mv '+pass_dir+'pec.pass '+pass_dir+file[0]+'_wl'+strtrim(i,1)+'.pec'
		close,/all
        ENDFOR
END

PRO generate_arf_pecs,type=type,start=start	;generate from Ar ADAS adf04 files from Adam Foster
	IF NOT keyword_set(start) THEN start=0
	z=18
	IF NOT keyword_set(type) THEN type='ls'
	adfpath='/home/adas/adas/adf04/coparf#18/'
	FOR i=start,z-1 DO BEGIN
		iz=(z-1)-i
		adf04=adfpath+'arf40_'+type+'#ar'+strtrim(iz,1)+'.dat'
		generate_208_pecs,adf04,nfiles=nfiles,lmax=300.0,amin=1.0e8,pass_dir=pass_dir,file=file,/rec
		IF nfiles NE -1 THEN BEGIN
			FOR j=0,nfiles-1 DO BEGIN
				spawn, 'scp '+pass_dir+file[0]+'_wl'+strtrim(j,1)+'.pec mlreinke@cmodws61.psfc.mit.edu:~/atomic_physics/adas/adf15/pec_wl/'
			ENDFOR
		ENDIF
        ENDFOR
END

PRO generate_arg_pecs,type=type,start=start		;generate from Ar OPEN ADAS adf04 files from Alessandra Giunta
	IF NOT keyword_set(start) THEN start=0
	z=18
	IF NOT keyword_set(type) THEN type='ls'
	adfpath='/home/reinke/adas/adf04/ar/'
	FOR i=start,z-1 DO BEGIN
		iz=(z-1)-i
		adf04=adfpath+type+'#ar'+strtrim(iz,1)+'.dat'
		generate_208_pecs,adf04,nfiles=nfiles,lmax=300.0,amin=1.0e6,pass_dir=pass_dir,file=file
		IF nfiles NE -1 THEN BEGIN
			FOR j=0,nfiles-1 DO BEGIN
				spawn, 'scp '+pass_dir+file[0]+'_wl'+strtrim(j,1)+'.pec mlreinke@cmodws61.psfc.mit.edu:~/atomic_physics/adas/adf15/pec_wl/'
			ENDFOR
		ENDIF
        ENDFOR
END

PRO generate_arm_pecs,type=type,start=start		;generate from Ar OPEN ADAS adf04 files that are arf40 + RR in adf48
	IF NOT keyword_set(start) THEN start=0
	z=18
	IF NOT keyword_set(type) THEN type='ls'
	adfpath='/home/reinke/adas/adf04/ar/'
	FOR i=start,12 DO BEGIN
		iz=(z-1)-i
		adf04=adfpath+'mlr13_'+type+'#ar'+strtrim(iz,1)+'.dat'
		generate_208_pecs,adf04,nfiles=nfiles,lmax=300.0,amin=1.0e8,pass_dir=pass_dir,file=file,/rec
		IF nfiles NE -1 THEN BEGIN
			FOR j=0,nfiles-1 DO BEGIN
				spawn, 'scp '+pass_dir+file[0]+'_wl'+strtrim(j,1)+'.pec mlreinke@cmodws61.psfc.mit.edu:~/atomic_physics/adas/adf15/pec_wl/'
			ENDFOR
		ENDIF
        ENDFOR
END

PRO generate_mom_pecs,type=type,start=start		;generate from Mo OPEN ADAS adf04 files from Martin O'Mullane
	IF NOT keyword_set(start) THEN start=0
	z=42
	IF NOT keyword_set(type) THEN type='ls'
	adfpath='/home/adas/adas/adf04/copmm#42/'
	FOR i=start,z-1 DO BEGIN
		iz=(z-1)-i
		adf04=adfpath+type+'#mo'+strtrim(iz,1)+'.dat'
		generate_208_pecs,adf04,nfiles=nfiles,lmax=300.0,amin=1.0e6,pass_dir=pass_dir,file=file
		IF nfiles NE -1 THEN BEGIN
			FOR j=0,nfiles-1 DO BEGIN
				spawn, 'scp '+pass_dir+file[0]+'_wl'+strtrim(j,1)+'.pec mlreinke@cmodws61.psfc.mit.edu:~/atomic_physics/adas/adf15/pec_wl/'
			ENDFOR
		ENDIF
        ENDFOR
END
