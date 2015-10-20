;               type       I     integer type of conversion
;                                    = 1 => standard form out, standard form in
;                                      2 => Eissner  form out, standard form in
;                                      3 => standard form out, Eissner  form in
;                                      4 => Eissner  form out, Eissner  form in


FUNCTION convlev,lev,type=type
	IF strpos(strlowcase(lev[0]),'s') NE -1 THEN type='std' ELSE type='eis'
	output=lev
	IF type EQ 'eis' THEN BEGIN
		FOR i=0,n(lev) DO BEGIN
			xxcftr,in=lev[i],out=out,type=3
			output[i]=strupcase(out)
        	ENDFOR
	ENDIF
	FOR i=0,n(lev) DO output[i]=strtrim(strtrim(output[i],1))
	RETURN,output
END

FUNCTION matchlev,elev,es,el,rlev,rs,rl,error=error
	nlev=n(elev)+1
	e2r=intarr(nlev)
	error=intarr(nlev)
	FOR i=0,nlev-1 DO BEGIN
		tmp=where(rlev EQ elev[i] AND rs EQ es[i] AND rl EQ el[i])
		IF n(tmp) NE 0 THEN error[i]=1
		IF tmp[0] EQ -1 THEN error[i]=-1
		e2r[i]=tmp[0]
	ENDFOR
	RETURN,e2r
END

FUNCTION adf_add_recdata,data,level,parent,rec
	tags=tag_names(data)
	newdata=create_struct(tags[0],data.(0))
	FOR i=1,n(tags) DO BEGIN
		CASE tags[i] OF
			'LEVEL_REC' : BEGIN
				newdata=create_struct(newdata,tags[i],level)
                        END
			'PARENT_REC' : BEGIN
				newdata=create_struct(newdata,tags[i],parent)
                        END
			'REC' : BEGIN
				newdata=create_struct(newdata,tags[i],rec)
                        END
			ELSE : newdata=create_struct(newdata,tags[i],data.(i))
                ENDCASE
	ENDFOR	
	RETURN,newdata
END

FUNCTION  make_rec_comment,adf04,adf48,e2r
	
	nelev=n(e2r)+1
	comment=strarr(8+nelev)
	cbnd='C---------------------------------------------------'
	comment[0]=cbnd
	comment[n(comment)]=cbnd
	comment[1]='C     adf04 source: '+adf04
	comment[2]='C     adf48 source: '+adf48
	comment[3]='C'
	comment[4]='C     adf04->adf48 level'
	FOR i=0,nelev-1 DO comment[i+5]='C         '+strtrim(i+1,1)+' -> '+strtrim(e2r[i]+1,1)
	comment[n(comment)-2]='C'
	comment[n(comment)-1]='C  Made Using ADAS208_ADD_RLINES.PRO by: '+logname()+' '+systime(0)
	RETURN,comment
END

PRO add_rlines,adf04,adf48,path,debug=debug
	read_adf04,file=adf04,fulldata=edata
	xxdata_09,file=adf48,fulldata=rdata
	elev=convlev(edata.cstrga)
	es=edata.isa
	el=edata.ila
	rlev=convlev(rdata.cstrga)
	rs=rdata.isa
	rl=rdata.ila
	e2r=matchlev(elev,es,el,rlev,rs,rl,err=err)
	IF total(abs(err)) NE 0 THEN print, 'warning found - not all adf04 levels matched'
	tmp=where(err EQ 0)
	level=tmp+1
	ntmp=n(tmp)+1
	nte=n(edata.te)+1
	rec=dblarr(ntmp,nte)
	FOR i=0,ntmp-1 DO rec[i,*]=10^interpol(alog10(reform(rdata.diel_res[e2r[tmp[i]],0,*])),alog10(rdata.tea),alog10(edata.te))	;hard coded to parent # 0
	tmp=where(finite(rec[*,0]))		;might have interpolate rates = 0 if there are rec from excited states
	rec=rec[tmp,*]
	level=level[tmp]
	parent=level*0+rdata.iprti[0]
	data=adf_add_recdata(edata,level,parent,rec)
	comment= make_rec_comment(adf04,adf48,e2r)
	write_adf04,comments=comment,fulldata=data,outfile=path
	IF keyword_set(debug) THEN stop
END

PRO add_rec_test
	adf04='/home/adas/adas/adf04/adas#10/cop98#10_ls#ne9.dat'
	adf48='/home/adas/adas/adf48/nrb05##/nrb05##_ne10ls.dat'
	path='/home/reinke/adas/adf04/test.dat'
	add_rlines,adf04,adf48,path
	print, 'compare R lines in '+path+' to '+adf04
END

PRO add_rec2arf
	adf48_base='/home/adas/adas/adf48/'
	adf48_dir='nrb05#'+['#','h','he','li','be','b','c','n','o','f','ne','na','mg']
	adf04_base='/home/adas/adas/adf04/coparf#18/'
	outpath='/home/reinke/adas/adf04/ar/'
	FOR i=0,12 DO BEGIN
		qstr=strtrim(17-i,1)
		adf04=adf04_base+'arf40_ls#ar'+qstr+'.dat'
		adf48=adf48_base+adf48_dir[i]+'/'+adf48_dir[i]+'_ar'+strtrim(18-i,1)+'ls.dat'
		path=outpath+'mlr13_ls#ar'+qstr+'.dat'
		add_rlines,adf04,adf48,path
	ENDFOR
END
;adf04='/home/adas/adas/adf04/coparf#18/arf40_ls#ar16.dat'
;adf48='/home/adas/adas/adf48/nrb05#h/nrb05#h_ar17ls.dat'
