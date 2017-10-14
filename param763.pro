PRO params

  nmonte = 10000.
  restore,'~/Dropbox/MyStructures/othertargets.dat'
  k2 = othertargets[where(othertargets.name eq 'EPIC_210696763')]
  plx = 1000./90.3;136.2
  plx_err = (1000./136.2)*(3./136.2)
  v = 16.416
  v_err = 0.105
  b = 18.673
  b_err = 0.135
  r = 16.00
  r_err = 0.01
  g = 17.44
  g_err = 0.01
  i = 16.41
  i_err = 0.01
  z = 13.95
  z_err = 0.01
  j = 12.5390
  j_err = 0.021
  h = 11.9580
  h_err = 0.023
  k = 11.7150
  k_err = 0.018
  
  vj = v-j
  teff = (2.840-1.3453*vj+0.3906*vj^2.-0.0546*vj^3.+0.002913*vj^4.)*3500.
  teff_err = sqrt(55^2.+60^2.) ;; not quite right, should mc this
  print,'T_vj = '+string(teff,format="(I4)")+'+/-'+string(teff_err,format="(I4)")
  rj = r-j
  teff = (2.445-1.2578*rj+0.4340*rj^2.-0.0720*rj^3.+0.004502*rj^4.)*3500.
  teff_err = sqrt(58^2.+60^2.) ;; not quite right, should mc this
  print,'T_rj = '+string(teff,format="(I4)")+'+/-'+string(teff_err,format="(I4)")
  print,'T_sp = '+string(k2.uh22model_Values[6],format="(I4)")
  
  feh = 0.03
  feh_err = 0.05
  plx = abs(plx+plx_err*randomn(seed,nmonte))
  mk = k-5.0*(alog10(1000./plx)-1.)
  mh = h-5.0*(alog10(1000./plx)-1.)
  mj = j-5.0*(alog10(1000./plx)-1.)
  
  readcol,'~/Desktop/BHAC15_iso.2mass.txt',bmass,bteff,blum,blogg,brad,blith,bmj,bmh,bmk,format='d,d,d,d,d,d,d,d,d'
  plot,bteff,bmk,psym=8,/xstyle,/ystyle
  oploterror,[teff],[median(mk)],[60],[stdev(mk)],psym=8,color=cgcolor('red'),errcolor=cgcolor('red')

  
  r2 = (1.9305-0.3466*mk+0.01647*mk^2.)*(1.+feh*0.04458)
  r2 += r2*0.0289*randomn(seed,nmonte)
  print,'Parallax rad: ',median(R2),stdev(r2)
  mass = 0.5858+0.3872*mk-0.1217*mk^2.+0.0106*mk^3.-2.7262d-4*mk^4.
  mass += mass*0.018*randomn(seed,nmonte)
  print,'Parallax Mass: ', median(mass),stdev(mass)
  den = mass/r2^3.
  tmp = den[sort(den)]
  print,'Density: ',median(den),tmp[n_Elements(tmp)*0.851]-median(den),median(den)-tmp[n_Elements(tmp)*0.159]


  t = k2.uh22model_Values[6]
  t += randomn(seed,nmonte)*60.
  t /= 3500.
  r3 = (16.7700-54.3210*t+57.6627*t^2.-19.6994*t^3.)*(1+0.4565*feh)
  print,'T_eff-Rad: ',median(r3),stdev(r3)
  mass = 0.0078368671+0.48699518*r3+1.6777663*r3^2.-1.2430190*r3^3.
  mass += 0.035*mass*randomn(seed,nmonte)
  print,'T_eff-mass: ',median(mass),stdev(mass)

  mk_inf = 13.130919 - 40.535438*r3+103.52254*r3^2.-141.62716*r3^3.+73.910274*r3^4.
  mk_inf+=randomn(seed,nmonte)*0.08
  dist_inf = 10.0^(((k-mk_inf)/5.)+1.)

  print,'inferred distance: ',median(dist_inf),stdev(dist_inf)

  stop

  bcj = 0.9672+0.4291*rj-0.05677*rj^2.+2.528d-3*rj^3.-0.05249*feh
  bcj+=0.012*randomn(seed,nmonte)
  ;;bcj = 0.8879+0.3563*vj-0.02791*vj^2.-0.04857*feh
  ;;bcj+=0.012*randomn(seed,nmonte)
                                ;bch = 1.915 + 0.2135*rj - 0.01582*rj^2.+0.09088*feh
                                ;bch+=0.021*randomn(seed,nmonte)
                                ;bck = 1.572 + 0.6529*rj - 0.1260*rj^2 + 9.746d-3*rj^3.+0.08987*feh
                                ;bck+=0.03*randomn(seed,nmonte)
  mbol = j+bcj
                                ;mbol = h+bch
                                ;mbol = k+bck
  mbol_sun =-26.8167
  c1 = (4.*!pi * 1d4) / 3.865d33 ;;3.939d33
  c2 = 3.08567758d18
  c3 = 4.84835548d-6
  c = c1*c2^2.
  fbol = (10.0^((mbol_sun - mbol)/2.5))/(c*c3^2.)
  
  teff = k2.uh22model_Values[6]
  teff += 60.*randomn(seed,nmonte)
  theta = ((fbol*1d12)^.5/(teff/2341.)^2.)                                                  ;; formula we always use, but reversed
  dist = 1000./plx
  theta_rads = (theta/1000./(60*60))*(!pi/180)                                       ;; convert to radians
  coeff = (3.0856d13/695500)/2.                                             ;; km/pc, Rsun/km, 2 for radius instead of diameter
  radius = (dist*theta_rads)*coeff                                          ;; atan(theta) ~ theta, been tested, it's fine.

  print,median(radius),stdev(radius)
  stop
  c1 = (4.*!pi * 1d4) / 3.865d33 ;
  c2 = 3.08567758d18
  c = c1*c2^2.0
  l = c*(fbol*(1000./plx)^2.0)
  logl = alog10(l)
  dist = 1000./plx
  theta = (2.*r2*6.96d10)/(dist*3.086d18)
  theta*=1000.*60*60*180./!pi
  teff = 2341.*((fbol/1d-12)/theta^2.)^(1./4.)
  print,median(teff),stdev(teff)


  

  stop



END


PRO rpm

  nmonte = 5000
  hyades = mrdfits('~/Dropbox/Structures/Hyades.fit',1,header)

  pmra = hyades.PMRA
  pmdec = hyades.pmde
  pmtot = sqrt(hyades.pmra^2.+hyades.pmde^2.)/1000.
  
  rpm_k = hyades.rmag +5.0*alog10(pmtot)+5.
  thick=4
  plot,hyades.rmag-hyades.jmag,rpm_k,psym=8,/xstyle,/ystyle,xthick=thick,ythick=thick,xtitle='r-J color',ytitle='H!LK!N',charthick=thick,charsize=1.5,yrange=[20,5],xrange=[0,6]
  ;;cghistoplot,rpm_k,/outline
  r = 15.219+0.031*randomn(seed,nmonte)
  j = 11.303+0.021*randomn(seed,nmonte)
  k = 10.444+0.019*randomn(seed,nmonte)
  pmra = 120.55+randomn(seed,nmonte)*3.3
  pmdec = -21.09+randomn(seed,nmonte)*3.2
  pmtot = sqrt(pmra^2.+pmdec^2.)/1000.
  rpm_k_pl = r+5.0*alog10(pmtot)+5.
  oploterror,[median(r-j)],[median(rpm_k_pl)],[stdev(r-j)],[stdev(rpm_k_pl)],psym=8,errthick=thick,color=cgcolor('red'),errcolor=cgcolor('red')
  oplot,[median(r-j)],[median(rpm_k_pl)],psym=8,color=cgcolor('red')
  
  stop

  bins = 100
  restore,'~/Dropbox/Structures/targets6.dat'
  t = targets[where(Targets.apass.r gt 0 and targets.apass.v gt 10.)]
  pmtot = sqrt(t.pmra^2.+t.pmdec^2.)
  x = t.apass.r-t.j
  y = t.apass.r+5.0*alog10(pmtot)+5
  l = where(x gt 1.8 and y gt 9)
  x = x[l]
  y = y[l]
  xarr = generatearray(min(x),max(x),bins)
  yarr = generatearray(min(y),max(y),bins)
  xbin = xarr[2]-xarr[1]
  ybin = yarr[2]-yarr[1]
  data = hist_2d(x,y,bin1=xbin,bin2=ybin,min1=min(xarr)-xbin/10.,max1=max(xarr)+xbin/10.,min2=min(yarr)-ybin/10.,max2=max(yarr)+ybin/10.)
  per = dblarr(max(data))
  for j = 0,max(data)-1 do per[j] = total(data[where(data gt j)])/total(data)
  levels = [(closest(0.9973,per))[0],(closest(0.9545,per))[0],(closest(0.6827,per))[0]]+2
  ;;contour,data,xarr,yarr,levels=levels,C_COLORS = [100,150,200],/FILL,/overplot


  stop


END


PRO prob_member

  pmra = 11.7;18.4;
  pmra_err = 2.;7.
  pmdec = -45.9;-56.6;
  pmdec_err = 2.9
  distance = 120.0   ;; photometric distance
  distance_err = 20. ;; error on photometric distance
  ;;vrad = 39.28      ;; average of Mace values
  ;;rv_err = 0.1      ;; stdev(rv)/sqrt(n_elements(rv))
  ;;rv = [39.33603,39.09584,39.13028,39.56312]
  ;;rv_err = [0.094348156,0.10122231,0.091012937,0.11221297]
  ;;vrad = mean(rv)
  ;;vrad_err = 0.1;stdev(rv)/sqrt(n_elements(rv))
  vrad = 22.58; 18.
  vrad_err = 0.12;7.
  ra = 52.834602
  dec = 18.316229

  xyz_errors, ra, dec, distance, 1/3600., 1/3600., distance_err, X, Y, Z, EX, EY, EZ
  uvw_errors, ra, dec, pmra, pmdec, vrad, distance, .1/3600., .1/3600., pmra_err, pmdec_err, vrad_err, distance_err, $
              U, V, W, E_U, E_V, E_W
  print,u,v,w
  print,e_u,e_v,e_w
  uvwxyz_star = [u,v,w,x,y,z]
  uvwxyz_star_err = [e_u,e_v,e_w,ex,ey,ez]

  uvwxyz_p = [-6.7,-25.0,-12.8,-107.0,25.9,-48.3]
  uvwxyz_p_err = [2.0,2.0,2.0,9,8,5] ;[0.9,0.5,0.5]

  uvwxyz_field = [-10.92,-13.35,-6.79,-0.18,2.10,3.27]
  uvwxyz_field_err = [23.22,13.44,8.97,53.29,51.29,50.70]

  ;print,'$U$ \kms & $'+string(u,format="(D5.1)")+'\pm'+string(e_u,format="(D5.1)")+'$ & This paper \\'
  ;print,'$V$ \kms & $'+string(v,format="(D5.1)")+'\pm'+string(e_v,format="(D5.1)")+'$ & This paper \\'
  ;print,'$W$ \kms & $'+string(w,format="(D5.1)")+'\pm'+string(e_w,format="(D5.1)")+'$ & This paper \\'
  ;print,'$X$ pc & $'+string(x,format="(D5.1)")+'\pm'+string(ex,format="(D5.1)")+'$ & This paper \\'
  ;print,'$X$ pc & $'+string(y,format="(D5.1)")+'\pm'+string(ey,format="(D5.1)")+'$ & This paper \\'
  ;print,'$X$ pc & $'+string(z,format="(D5.1)")+'\pm'+string(ez,format="(D5.1)")+'$ & This paper \\'
  ;; print,uvwxyz_h
  ;; print,uvwxyz_h_err
  ;; print,uvwxyz_star
  ;; print,uvwxyz_star_err
  ;; print,uvwxyz_h-uvwxyz_star
  ;; stop

  uvwxyz_tot_err = sqrt(uvwxyz_star_err^2.+uvwxyz_p_err^2.)
  P_th_H = 1./((sqrt(2*!pi)*uvwxyz_tot_err)) * exp((-1./2.)*((uvwxyz_star-uvwxyz_p)/uvwxyz_tot_err)^2.)

  uvwxyz_err_tot = sqrt(uvwxyz_star_err^2.+uvwxyz_field_err^2.)
  P_th_f = 1./((sqrt(2*!pi))* uvwxyz_err_tot) * exp((-1./2.)*((uvwxyz_star-uvwxyz_field)/uvwxyz_err_tot)^2.)

  a = 300d;519.                      ;307.
  b = 204223d;721485.                   ;204223.
  tot_h = (p_th_h[0]*p_th_h[1]*p_th_h[2]*p_th_h[3]*p_th_h[4]*p_th_h[5])*(a/b)
  tot_f = (p_th_f[0]*p_th_f[1]*p_th_f[2]*p_th_f[3]*p_th_f[4]*p_th_f[5])*((b-a)/b)
  print,tot_h/(tot_f+tot_h)

  ;; v = 15.881
  ;; b = 17.76
  ;; r = 15.219
  ;; g = 16.562
  ;; i = 13.571
  ;; j = 11.303
  ;; h = 10.732
  ;; k = 10.444
  ;; restore,'~/Desktop/apass/apassN.dat'
  ;; rad = 10.
  ;; hy = hyades[where(abs(hyades.raj2000-ra) lt rad and abs(hyades.dej2000-dec) lt rad*3)]
  ;; print,n_elements(hy)
  ;; apass = apassn[where(abs(apassn.ra-ra) lt rad and abs(apassn.dec-dec) lt rad*3 and apassN.v gt 10 and apassN.r lt 20)]; and abs(apassn.r-r) lt 2 and abs(apassn.v-v) lt 2)]
  ;; print,n_elements(apass)
  ;; stop
  
  ;; inhyades = 0
  ;; nothyades = 0
  ;; hyades = mrdfits('~/Dropbox/Structures/Hyades.fit',1,header)
  ;; readcol,'~/Desktop/K2Campaign4targets.csv',epic,ra,dec,kp,go,format='d,d,d,d,a',delimiter=',',/silent
  ;; l = where(ra gt 0 and dec gt 0)
  ;; ra = ra[l]
  ;; dec = dec[l]
  ;; epic = epic[l]
  ;; kp = kp[l]
  ;; go = go[l]
  ;; for i = 0,n_elements(epic)-1 do begin
  ;;    gcirc,2,ra[i],dec[i],hyades.raj2000,hyades.dej2000,dist
  ;;    l = where(dist lt 10)
  ;;    if l[0] ne -1 then begin
  ;;       inhyades++
  ;;    endif else begin
  ;;       nothyades++
  ;;    endelse
  ;; endfor
  ;; print,inhyades,nothyades
  ;; stop
  
  
END



PRO uvwxyz_hyades

  thick = 4
  hyades = mrdfits('~/Dropbox/Structures/Hyades.fit',1,header)

  hip_dists = dblarr(n_elementS(hyades))
  hip_dists_err = hip_dists
  for i = 0,n_Elements(hyades)-1 do begin
     info = queryvizier('I/311/HIP2',[hyades[i].RAJ2000,hyades[i].DEJ2000],0.1,/silent)
     if isarray(info) then begin
        hip_dists[i] = info.plx
        hip_dists_err[i] = info.e_plx
        print,info.plx,info.e_plx,hyades[i].plx
     endif
  endfor
  stop
  
  dist = 1000./hyades.plx
  dist_err = dist*(hyades.e_plx/hyades.plx)
  xyz_errors,hyades.RAJ2000,hyades.DEJ2000,dist, .1/3600., .1/3600.,dist_err,xh,yh,zh,exh,eyh,ezh

  dist = 1000./hyades.plx
  dist_err = dist*(hyades.e_plx/hyades.plx)
  uvw_errors,hyades.RAJ2000,hyades.DEJ2000,hyades.PMRA,hyades.pmde,hyades.RV,dist,.1/3600.,.1/3600.,hyades.e_pmra,hyades.e_pmde,1,dist_err,ut,vt,wt,e_ut,e_vt,e_wt

  l = where(hyades.seq eq 169)
  !p.multi=[0,3,2]
  contouring,xh,yh
  oploterror,[xh[l]],[yh[l]],[exh[l]],[eyh[l]],errthick=thick,psym=8
  contouring,xh,zh
  oploterror,[xh[l]],[zh[l]],[exh[l]],[ezh[l]],errthick=thick,psym=8
  contouring,yh,zh
  oploterror,[yh[l]],[zh[l]],[eyh[l]],[ezh[l]],errthick=thick,psym=8
  contouring,ut,vt
  oploterror,[ut[l]],[vt[l]],[e_ut[l]],[e_vt[l]],errthick=thick,psym=8
  contouring,ut,wt
  oploterror,[ut[l]],[wt[l]],[e_ut[l]],[e_wt[l]],errthick=thick,psym=8
  contouring,vt,wt
  oploterror,[vt[l]],[wt[l]],[e_vt[l]],[e_wt[l]],errthick=thick,psym=8
  stop
  
END

PRO contouring,x,y

  charsize=1.5
  bins = 10
  xarr = generatearray(min(x),max(x),bins)
  yarr = generatearray(min(y),max(y),bins)
  xbin = xarr[2]-xarr[1]
  ybin = yarr[2]-yarr[1]
  data = hist_2d(x,y,bin1=xbin,bin2=ybin,min1=min(xarr)-xbin/10.,max1=max(xarr)+xbin/10.,min2=min(yarr)-ybin/10.,max2=max(yarr)+ybin/10.)
  per = dblarr(max(data))
  for j = 0,max(data)-1 do per[j] = total(data[where(data gt j)])/total(data)
  levels = [(closest(0.9973,per))[0],(closest(0.9545,per))[0],(closest(0.6827,per))[0]]+2
  contour,data,xarr,yarr,levels=levels,C_COLORS = [100,150,200],/FILL,charsize=charsize ;,/overplot

END


PRO phot_dist

  nmonte = 10000.
  restore,'~/Dropbox/MyStructures/othertargets.dat'
  k2 = othertargets[where(othertargets.name eq 'EPIC_210696763')]

  v = 16.416
  v_err = 0.105
  b = 18.673
  b_err = 0.135
  r = 16.00
  r_err = 0.05
  g = 17.44
  g_err = 0.1
  imag = 16.41
  i_err = 0.1
  z = 13.95
  z_err = 0.05
  j = 12.5390
  j_err = 0.021
  h = 11.9580
  h_err = 0.023
  k = 11.7150
  k_err = 0.018

  restore,'~/Dropbox/Radii/feiden5.dat'
  ;;feiden = feiden[where(feiden.teff gt 300 and feiden.teff lt 3600)]
  ;feiden = feiden[where(feiden.halpha gt 0)]; 'GJ 544 B')]

  ;;mk = k-5.0*(alog10(dist)-1.)

  errs = [b_err,v_err,r_err,i_err,j_err,h_err,k_err]
  chi = dblarr(n_elements(feiden))
  dist = chi
  dist_err = dist
  nmonte = 5000
  for i = 0,n_elements(feiden)-1 do begin
     diffs = [b-feiden[i].synthetic.b,v-feiden[i].synthetic.v,r-feiden[i].synthetic.r,imag-feiden[i].synthetic.i,j-feiden[i].synthetic.j,h-feiden[i].synthetic.h,k-feiden[i].synthetic.k]
     offset = wmean(diffs,errs,error=offset_err)
     f_errs = [feiden[i].synthetic.b_err,feiden[i].synthetic.v_err,feiden[i].synthetic.r_err,feiden[i].synthetic.i_err,feiden[i].synthetic.j_err,feiden[i].synthetic.h_err,feiden[i].synthetic.k_err]
     tot_err = sqrt(errs^2.+f_errs^2.)
     diffs-=offset
     chi[i] = total(diffs^2./tot_err^2.)/n_elements(errs)
     dist[i] = 10.0^(offset/5.)*feiden[i].distance

     tmp1 = feiden[i].distance+randomn(seed,nmonte)*feiden[i].distance_err
     tmp2 = offset+randomn(seed,nmonte)*offset_err
     dist_tmp = 10.0^(tmp2/5.)*tmp1
     dist_err[i] = stdev(dist_tmp)
     if chi[i] lt 10 then begin  ;and (dist[i] gt 50 or dist[i] lt 35) then begin
        print,diffs,diffs/errs
        print,dist[i]
        stop
     endif
  endfor
  set_plot,'x'
  print,feiden[wherE(chi eq min(chi))].teff
  plot,dist,chi,psym=8,/ylog,/xstyle,/ystyle,xtitle='Distance',ytitle='Chi^2'

  ll = where(chi lt min(chi)*3.0)
  print,median(dist[ll]),stdev(dist[ll]),robust_sigma(dist[ll])

  set_plot,'PS'
  device,filename='Chi_dist.eps',/encapsul,/color
  chis = greek('chi')
  nu = greek('nu')
  plot,chi,dist,psym=8,/xlog,/xstyle,/ystyle,xrange=[0.8,400],xtitle=chis+'!L'+nu+'!U2!N',ytitle='Distance (pc)',charsize=1.5,charthick=4,xthick=4,ythick=4,yrange=[30,200]
  x1 = min(chi[ll])-0.05
  x2 = max(chi[ll])+0.05
  y1 = min(dist[ll])-1
  y2 = max(dist[ll])+1
  oplot,[x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],thick=4,color=cgcolor('red')
  oplot,chi,dist,psym=8
  ;;polyfill,[0.7,5,5,0.7,0.7],[25,25,70,70,25],thick=4,color=cgcolor('red'),fill=0
  device,/close
  stop
  
END


PRO rv

  set_plot,'x'
                                ;epochs = [2457295.8693865747, 2457293.8594444450, 2457289.8955324078, 2457288.9227083339]-2454833d0
                                ;rv = [39.33603, 39.09584, 39.13028, 39.56312]*1000.
                                ;rv_err = [0.094348156, 0.10122231, 0.091012937, 0.11221297]*1000.

  readcol,'JDs.txt',dummy,object,spt,name,dash,poop,epochs,format='d,a,a,a,a,a,d'
  
  readcol,'hyades_velocities.txt',index,rv,rv_err,nv,nt,format='d,d,d,d,d'
  readcol,'hyades_mann_absolute.txt',jd,rv,rv_err,format='d,d,d'
  
  rv*=1000d0
  rv_err*=1000d0
  thick = 4
  
  for i = 0,n_elements(rv)-1 do begin
     del = ' & '
     print,string(epochs[i]-2454833d0,format="(D12.5)")+del+string(rv[i]-mean(rv),format="(I5)")+del+string(rv_err[i],format="(I5)")+'\\'
  endfor
  epochs-=(2454833d0+2229.57935d0) ;; +3.48454d0*64d0


  set_plot,'PS'
  yrange = [-550,550]
  device,filename='RVs.eps',/encapsul,/color
  !p.multi=[0,1,2]
  !y.margin=[1,1]
  period = 3.48454d0            ;*2d0
  ep = epochs mod period
  rv_err*=0.0
  rv_err+=153d0
  x = [ep/period,(ep/period)+1d0]
  y = [rv-mean(rv),rv-mean(rv)]
  x_err = [ep*0.,ep*0.]
  y_err = [rv_err,rv_err]
  ploterror,x,y,x_err,y_err,psym=8,errthick=thick,xthick=4,ythick=4,xrange=[0,2],yrange=yrange,/ystyle,/xstyle,ytitle='Relative RV (m/s)                  ',errcolor=cgcolor('black'),background=cgcolor('white'),COLOR=cgcolor('black'),charsize=1.5,charthick=thick,xtickname=strarr(10)+' '
  l = where(x lt 0.5 or x gt 1.5)
  oploterror,x[l],y[l],x_err[l],y_err[l],psym=8,color=cgcolor('grey'),errthick=thick,errcolor=cgcolor('grey')

  radius = 0.295227
  mass = 0.29363864
  e = 0
  omega = !pi/2
  t0 = 2229.5795898438

  fit = finalfit(epochs,rv-mean(Rv),rv_err)

  tmp1 = generatearray(0,period,125)
  deviates = mporbit(fit.params,date=tmp1,rv=tmp1,err=tmp1,ovel1=newk)
  ;;oplot,[tmp1/period,(tmp1/period)+1d0],[newk,newk],color=cgcolor('red'),thick=thick

  mstar = 0.294
  mstar_err = 0.021
  ecc = sqrt(fit.params[2]^2 + fit.params[3]^2)
  mass = fit.params[1] * (mstar * 1.98892e30)^(.66666)*(fit.period*86400d0)^(.3333333)*sqrt(1-ecc^2)
  mass /= (2*3.14159*6.673e-11)^(.333333)
  mass/= 1.8987e27
  nmonte = 5000
  mstar = mstar+mstar_err*randomn(seed,nmonte)
  ecosw = fit.params[2]+fit.perror[2]*randomn(seed,nmonte)
  esinw = fit.params[3]+fit.perror[3]*randomn(seed,nmonte)
  k = fit.params[1]+fit.perror[1]*randomn(seed,nmonte)
  period = fit.period+dblarr(nmonte)
  ecc = sqrt(ecosw^2d0+esinw^2d0)
  mass = k * (mstar * 1.98892e30)^(.66666)*(period*86400d0)^(.3333333)*sqrt(1-ecc^2)
  mass /= (2*3.14159*6.673e-11)^(.333333)
  mass/= 1.8987e27
  ;; errors
  l = where(finite(mass) eq 1)
  print,median(mass[l]),stdev(mass[l])


  period = 3.48454d0                                                                          ;*2d0
  p1 = {value:0.1d, fixed:0b, limited:[1, 1], limits:[0.0d, 10d]}                             ;; mass
  p2 = {value:period/2d0, fixed:0b, limited:[0, 0], limits:[period/2d0-0.05d,period/2d0+0.05d]} ;; T0
  p3 = {value:mean(rv)*1d0, fixed:1b, limited:[0, 0], limits:[min(rv),max(rv)]}
  p4 = {value:1d-1, fixed:0b, limited:[1,1],limits:[0d,0.4d0]}      ;; eccentricity
  p5 = {value:!pi/2d0, fixed:0b, limited:[1,1],limits:[0d,2d0*!pi]} ;; omega
  parinfo = [p1,p2,p4,p5]
  params = mpfitfun('K',epochs,rv-mean(Rv),rv_err,parinfo=parinfo,PERROR=PERROR)

  ;; see what jitter is needed to explain data:
  newk = k(epochs,params)
  rchisq = total((newk-(rv-mean(rv)))^2./rv_err^2.)/(n_elements(rv)-n_elements(parinfo))
  syserr = 1d0
  while rchisq gt 1 do begin
     syserr*=1.5
     rchisq = total((newk-(rv-mean(rv)))^2./(syserr^2.+rv_err^2.))/(n_elements(rv)-n_elements(parinfo))
  endwhile
  print,syserr
  
  
  print,'Mass = ',params[0],'+/-',perror[0]
  ;newe = generatearray(min(epochs),max(epochs),100)
  ;newk = k(newe,params)
  ;phase = [(newe mod period)/period,(newe mod period)/period+1d0]
  ;newk = [newk,newk]
  ;newk = newk[sort(phase)]
  ;phase = phase[sort(phase)]
  ;oplot,phase,newk,color=cgcolor('red'),thick=thick
  
  params[0] = 0.05
  newe = generatearray(min(epochs),max(epochs),100)
  newk = k(newe,params)
  phase = [(newe mod period)/period,(newe mod period)/period+1d0]
  newk = [newk,newk]
  newk = newk[sort(phase)]
  phase = phase[sort(phase)]
  oplot,phase,newk,color=cgcolor('teal'),thick=thick

  params[0] = 1.0
  newe = generatearray(min(epochs),max(epochs),100)
  newk = k(newe,params)
  phase = [(newe mod period)/period,(newe mod period)/period+1d0]
  newk = [newk,newk]
  newk = newk[sort(phase)]
  phase = phase[sort(phase)]
  oplot,phase,newk,color=cgcolor('blue'),thick=thick

  params[0] = 3.
  newe = generatearray(min(epochs),max(epochs),100)
  newk = k(newe,params)
  phase = [(newe mod period)/period,(newe mod period)/period+1d0]
  newk = [newk,newk]
  newk = newk[sort(phase)]
  phase = phase[sort(phase)]
  oplot,phase,newk,color=cgcolor('red'),thick=thick
  sharpcorners,thick=thick


  legend,['3xJupiter','Jupiter','Neptune'],color=[cgcolor('red'),cgcolor('blue'),cgcolor('teal')],psym=[8,8,8],/bottom,/right,charsize=1.3,charthick=4,box=0,textcolor=[cgcolor('red'),cgcolor('blue'),cgcolor('teal')],position=[1.95,-500]
  ;;oploterror,ep,rv-mean(rv),ep*0.,rv_err,psym=8,errthick=4,errcolor=cgcolor('black')

  !y.margin=[3,-1]
  rotp = 1.88 ;;/2d0
  ep = (epochs-0.2) mod rotp
  x = [ep/rotp,ep/rotp+1d0]
  y = [rv-mean(rv),rv-mean(rv)]
  x_err = x*0d0
  y_err = [rv_err,rv_err]
  ploterror,x,y,x_err,y_err,psym=8,errthick=4,xthick=4,ythick=4,xrange=[0,2],yrange=yrange,/ystyle,/xstyle,xtitle='Phase',errcolor=cgcolor('black'),background=cgcolor('white'),COLOR=cgcolor('black'),charsize=1.5,charthick=thick
  l = where(x lt 0.5 or x gt 1.5)
  oploterror,x[l],y[l],x_err[l],y_err[l],psym=8,color=cgcolor('grey'),errthick=thick,errcolor=cgcolor('grey')


  x = generatearray(0,2,100)
  y=0.015*8*1000d0*sin((2d*!pi) * x)
  oplot,x,y,thick=thick,color=cgcolor('red')

  legend,['Rotational Jitter'],color=[cgcolor('red')],box=0,/bottom,/right,psym=[8],charsize=1.3,charthick=4,textcolor=[cgcolor('red')],position=[1.95,-500]

  sharpcorners,thick=thick
  device,/close

  stop
  !p.multi=[0,1,1]
  !except=0
  set_plot,'x'
  device,retain=0
  ep = epochs mod period
  es = 0.;generatearray(0.,0.5,10)
  ws = 0.;generatearray(0,!pi,10)
  mp = generatearray(0.0,20.0,150)
  t0 = generatearray(2229.,2230,100);generatearray(2456.0742,2456.1,100) ;generatearray(2229.5795898438,2229.5795898438+3.48454,100);; 2457062.57910d0+65*3.48454-2454833d0
  bestchisq = 1d4
  maxmass = dblarr(1)
  maxe = maxmass
  for i = 0,n_elements(es)-1 do begin
     for j = 0,n_elementS(ws)-1 do begin
        for l = 0,n_elements(mp)-1 do begin
           for m = 0,n_elements(t0)-1 do begin
              v = k(epochs,[mp[l],t0[m],es[i],ws[j]])
              rchisq = total((v-(rv-mean(rv)))^2./(100d^2.+rv_err^2.))/n_elements(v)
              if rchisq lt 9. then begin
                 maxmass = [maxmass,mp[l]]
                 maxe = [maxe,es[i]]
                 ;;print,mp[l],rchisq
              endif
              if rchisq lt bestchisq then begin
                                ;if mp[l] gt 2 then begin; lt bestchisq then begin
                 bestchisq = rchisq
                 ploterror,ep,rv-mean(rv),ep*0.,rv_err,psym=8,errthick=4,xthick=4,ythick=4,xrange=[0,period],yrange=[-600,600],errcolor=cgcolor('orange')
                 oplot,ep,v,psym=8,color=cgcolor('red')
                 print,mp[l],rchisq
              endif
           endfor
        endfor
     endfor
     print,max(maxmass);,bestchisq
  endfor
  print,max(maxmass)
  print,maxe[where(maxmass eq max(maxmass))]
  
  stop


END

FUNCTION K,x,p,e=e,w=w
  period = 3.48454
                                ;if n_elements(e) eq 0 then e=0
                                ;if n_elements(w) eq 0 then omega = !pi/2 else omega = w
  mp = p[0]
  t0 = p[1]
  ;;rv_sys = p[2]
  e = p[2]
  omega = p[3]
  theta = (x-t0)*(2.*!pi/period) mod (2.*!pi)
  mass = 0.29363864
  k = 203.255*(1./period)^(1./3.)*(mp)*((1./mass)^(2./3.))*(cos(theta-omega)+e*cos(omega))
                                ;k-=rv_sys
  return,k
END


;04130560+1514520 04 13 05.61  15 14 52.00 99    15.372     0.150 C    10.444     0.017 44.8 R   0.20 -9999.0000 N n/a -9999.000 -9999.000 -9999.    -9999.    


function finalfit,date,rv,rv_err,circ=circ,forcek=forcek

  nparam = 6
  period = 3.484551d0
;Setup parameters for fit
  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                       limits:[0.D,0]}, 6)
;T0 should be after minimum date and within a period.
  Tmid = min(date) + period/2.0D
  parinfo[0:4].value = [period,172d0,0.01D,0.01D,Tmid]

;ecosw and esinw go between -1 and 1 and K should be greater than 0
;T0 should be after minimum date and within a period.
  parinfo[0].value = 3.484551
  parinfo[0].fixed = 1
  parinfo[1].limited[0] = 1
  parinfo[1].limits[0] = 0d0    ;2*amp/3
  if n_elements(forcek) gt 0 then begin
     parinfo[1].fixed = 1
     parinfo[1].value = forcek
  endif
  parinfo[2].limited = [1,1]
  parinfo[3].limited = [1,1] 
  parinfo[4].limited = [1,1]
  parinfo[2].limits = [-0.5D,0.5D];from transit fit[-1.0D,1.0D]
  parinfo[3].limits = [-0.5d,0.5d];from transit fit[-1.0D,1.0D]
  parinfo[4].limits = [min(date),min(date)+ period]

  ;period = params[0]
  ;k = params[1]
  ;ecosw = params[2]
  ;esinw = params[3]
  ;ep_peri = params[4]

;If circular set then ecosw and esinw are both 0 and fixed
  if(keyword_set(circ)) then begin
     parinfo[2].fixed = 1
     parinfo[3].fixed = 1
     parinfo[2].value = 0.0
     parinfo[3].value = 0.0
  endif

;Setup structure to take x,y,and err

  functargs = {date:date, rv:rv, err:rv_err} 

;Setup structure to store results
  fit = {params:dblarr(nparam),chi:0.D,cov:dblarr(nparam,nparam),dof:0.D,period:0.D,tperi:0.D,perror:dblarr(nparam)}

  fit.params = MPFIT('mporbit', PARINFO=parinfo,FUNCTARGS=functargs,BESTNORM=chi,DOF=dof,COV=cov,PERROR=perror,MAXITER=200,QUIET=soft,ERRMSG=errmsg)

;Store parameter errors
  if(keyword_set(perror)) then begin
     fit.perror = perror
     fit.cov = cov
  endif

;Store chi2 info
  fit.chi = chi
  fit.dof = dof

;Store Period
  fit.period = (2.0D * !DPI / fit.params[0])

;Store T_peri
;Get eccentricity and argp
ecc = sqrt(fit.params[2]^2 + fit.params[3]^2)

if(ecc eq 0.0) then $
   argp = 0.0 $
else $
   argp = atan(fit.params[3],fit.params[2])

fit.tperi =  convert_t2p(fit.params[4],fit.period,argp,ecc)

;Send fit values back
  return,fit
end


function mporbit,params,date=date,rv=rv,err=err,ovel1=ovel1
  common dataset, dset, dstar

;;Separate out the parinfo structure
  period = params[0]
  k = params[1]
  ecosw = params[2]
  esinw = params[3]
  ep_peri = params[4]
  v0 = params[5]
  slope = 0.0                   ;params[6]
  k2 = 0.0                      ;params[7]

;Apply the offset to the rv
  for i=8L,n_elements(params)-1 do begin
     num = i-8
     rv[where(dset eq num)] = rv[where(dset eq num)] + params[i]
  endfor 

  date1 = date;[where(dstar eq 0)]
  rv1 = rv;[where(dstar eq 0)]
  err1 = err;[where(dstar eq 0)]

;Get eccentricity and argp
  ecc = sqrt(ecosw^2 + esinw^2)

  if(ecc eq 0.0) then $
     argp = 0.0 $
  else $
     argp = atan(esinw,ecosw)

;I need E the eccentric anomaly.  This will be returned from function eanomaly
  enew = eanomaly(date1,ep_peri,period,ecc) 

;Get the true anomaly in both sine and cosine.  These are identical equations
;just used sin^2 + cos^2 = 1 to transform
  cosf = (cos(enew) - ecc)/(1.D - ecc*cos(enew))
  sinf = sqrt(1.D - ecc^2)*sin(enew)/(1.D - ecc*cos(enew))

;Find the radial velocity
;     vel = v0 + k*(cos(argp + f) + ecc * cos(argp))
;Since I don't have f must use the cosine sum formula
  ovel1 = v0 + k*(cos(argp)*cosf -sin(argp)*sinf +ecc*cos(argp)) + slope*date1

;Determine the standard deviant
  deviates1 = (rv1 - ovel1) / err1

  return,deviates1
end


;Inputs: T transit
;Outputs: T periastron
function convert_t2p,date,period,arg_p,ecc

;Using equations from Kane & von Braun (2009)
;E = 2*arctan ( sqrt(1-e)/sqrt(1+e) *tan((pi/2 - argp) /2)
;At epoch of transit arg_p + f = pi/2

  eanomoly = sqrt(1-ecc)/sqrt(1+ecc) * tan((!DPI/2 - arg_p)/2)
  eanomoly = 2 * atan(eanomoly)

  manomoly = eanomoly - ecc*sin(eanomoly)

  correct = period * manomoly / (2 * !DPI)

  date = date - correct

  return,date
end

;Inputs: JD, Epoch of Periastron, Period, Eccentricity
;Ouputs: eccentric anomaly

function eanomaly,data,ep_peri,period,ecc

  ;;ep_peri=2457289.0746898437
;The eccentric anomaly is related to time through the mean anomaly.  These
;relations are M = 2*pi/period *(t-t_0) and M = E - e * sin E.
;This cannot be solved analytically, so the Newton-Raphson iteration method
;must be used.  The formulation for this equation is then:
;E_(n+1) = M - e(E_n*cos E_n - sin E_n) / (1 - e*cos E_n)
;Charles & Tatum 1998 suggest that E_0 should equal pi

;Convergence criteria (same as for RVSIM)
  conv = 1.0d-4
  
;m must go between 0 and 2pi
  m = (data - ep_peri)/period
  m = 2D*!DPI * (m - floor(m))      
  
  eold = dblarr(n_elements(m))+!DPI
  enew = dblarr(n_elements(m))
;Run the Newton-Raphson method
  for i=0L,n_elements(m)-1 do begin
     n=0
     while(1) do begin
        enew[i] = (m[i] - ecc * (eold[i] * cos(eold[i]) - sin(eold[i]))) $
                  /(1-ecc*cos(eold[i]))
;Added in the $eold < $conv to avoid $eold of 0 issues.
        if(eold[i] lt conv) then break 
        if(abs((enew[i] - eold[i])/eold[i]) lt conv) then $
           break
        if(n > 50) then begin
           print, 'Too many iterations!'
           break
        endif

        eold[i] = enew[i]            
        n++
     endwhile
  endfor
  return,enew             
end


PRO Hyades_map

  hyades = mrdfits('~/Dropbox/Structures/Hyades.fit',1,header)
  hyades2 = mrdfits('~/Dropbox/Structures/Hyades2.fit',1,header)
  hyades2 = hyades2[where(hyades2.src eq 1)]

  distance = 45.7
  distance_err = 3.4

  dists = 1000./hyades.PLX
  dists_err = (hyades.e_plx/hyades.plx)*dists
  
  xyz_errors, hyades.raj2000, hyades.dej2000, dists, 0.1/3600., 0.1/3600., dists_err, X, Y, Z, EX, EY, EZ
  uvw_errors,  hyades.raj2000, hyades.dej2000, hyades.PMRA, hyades.pmde, hyades.rv, dists, 0.1/3600., 0.1/3600., hyades.e_pmra, hyades.e_pmde, 1.1, dists_err, $
               U, V, W, E_U, E_V, E_W

   

  dists = 1000./hyades2.PLX
  dists_err = 0.1*dists
  
  xyz_errors, hyades2.raj2000, hyades2.dej2000, dists, 0.1/3600., 0.1/3600., dists_err, X2, Y2, Z2, EX2, EY2, EZ2
  uvw_errors,  hyades2.raj2000, hyades2.dej2000, hyades.PMRA, hyades.pmde, hyades.rv, dists, 0.1/3600., 0.1/3600., 2, 2, 1.1, dists_err, $
              U2, V2, W2, E_U2, E_V2, E_W2

  u = [u,u2]
  v = [v,v2]
  w = [w,w2]
  x = [x,x2]
  y = [y,y2]
  z = [z,z2]

  dists = [1000./hyades.plx,1000./hyades2.plx]
  k = [hyades.ksmag,hyades2.kmag]
  mk = k-5.0*(alog10(dists)-1.0)

  ;; set_plot,'x'
  ;; p = plot3d(x,y,z,'o',/sym_filled,xrange=[min(x),max(x)],yrange=[min(y),max(y)],zrange=[min(z),max(z)],axis_style=2,depth_cue=[0,2],/perspective,rgb_table=33,vert_colors=(BYTSCL(mk)),shadow_color='black',xy_shadow=1,yz_shadow=1,xz_shadow=1,xtitle='X',ytitle='Y',ztitle='Z',SYM_OBJECT=ORB())
  ;; ax = p.AXES
  ;; ax[2].HIDE = 1
  ;; ax[6].HIDE = 1
  ;; ax[7].HIDE = 1
  ;; stop
  

  pmra = 120.55                 ;ppmxl
  pmra_err = 3.3
  pmdec = -21.09 ;;ppmxl
  pmdec_err = 3.2 
  distance = 45.7
  distance_err = 3.4
  ra = 63.273378
  dec = 15.247777
  vrad = 38637.730/1d3
  vrad_err = 153./1d3
  k = 10.444

  vrad = 39.76d0
  vrad_err = 0.1d0
  ra = am_racnv('04 29 38.993')
  dec = am_deccnv('+22 52 57.80')
  pmra = 83.0;;85.8;81.8
  pmra_err = 0.9;1.2;1.0
  pmdec = -35.7;-34.0;-35.2
  pmdec_err = 0.9               ;0.9
  distance = 59.4;58.12
  distance_err = 2.8            ;2.
  distance = 63
  distance_err = 10
  k = 8.368
  
  !p.font=0
  uvw_errors, ra, dec, pmra, pmdec, vrad, distance, 0.1/3600., 0.1/3600., pmra_err, pmdec_err, vrad_err, distance_err, $
              Uh, Vh, Wh, E_U, E_V, E_W
  xyz_errors, ra,dec, distance, 0.1/3600., 0.1/3600., distance_err, X2, Y2, Z2, EX2, EY2, EZ2
  mkh = k-5.0*(alog10(distance)-1.)
  print,X2,y2,z2
  print,ex2,ey2,ez2
  readcol,'BHAC15_800.txt',bmass,bteff,blum,blogg,brad,blith,bmj,bmh,bmk,format='d,d,d,d,d,d,d,d,d',/silent
  mkh = interpol(bmass,bmk,mkh)
  mk = interpol(bmass,bmk,mk)
  masses = [0.2,0.5,0.8,1.1,1.4]

  size1 = 2.2
  size2 = 1.75
  ocolor1 = 'white'              ;'black';
  ocolor2 = 'blue';'black'
  set_plot,'PS'
  device,filename='XYZ.eps',/encapsul,/color,xsize=25,ysize=10
  position = dblarr(4)
  !x.margin=[8,2]
  !y.margin=[4,5]
  charsize = 2.2
  thick=7
  symsize=0.6
  !p.multi=[0,3,1]
  ;plot,u,v,psym=8,/xstyle,/ystyle
  ;plot,u,w,psym=8,/xstyle,/ystyle
                                ;plot,v,w,psym=8,/xstyle,/ystyle
  position[0] = !x.window[0]
  cgloadct,0
  plotsym,0,/fill
  plot,x,y,psym=8,/xstyle,/ystyle,xthick=4,ythick=4,charsize=charsize,charthick=thick,xtitle='X (pc)',ytitle='Y (pc)',symsize=symsize
  ct = 3
  clip = [25,225]
  cgloadct,ct,clip=clip
  colors = generatearray(min(masses),max(masses),100)
  vmin = 0
  vmax = 255
  colors2 = generatearray(vmin,vmax,100)
  for i = 0,n_elements(x)-1 do begin
     c=round(colors2[closest(mk[i],colors)])
     oplot,[x[i]],[y[i]],psym=8,color=c,symsize=symsize
  endfor
  plotsym,3,/fill
  c=round(colors2[closest(mkh,colors)])
  oplot,[x2],[y2],psym=8,color=cgcolor(ocolor1),symsize=size1
  oplot,[x2],[y2],psym=8,color=cgcolor(ocolor2),symsize=size2
  oplot,[x2],[y2],psym=8,color=c,symsize=1.0
    sharpcorners,thick=thick

  plotsym,0,/fill
  cgloadct,0
  plot,x,z,psym=8,/xstyle,/ystyle,xthick=4,ythick=4,charsize=charsize,charthick=thick,xtitle='X (pc)',ytitle='Z (pc)',symsize=symsize
  ct = 3
  clip = [25,225]
  cgloadct,ct,clip=clip
  colors = generatearray(min(masses),max(masses),100)
  vmin = 0
  vmax = 255
  colors2 = generatearray(vmin,vmax,100)
  for i = 0,n_elements(x)-1 do begin
     c=round(colors2[closest(mk[i],colors)])
     oplot,[x[i]],[z[i]],psym=8,color=c,symsize=symsize
  endfor
  plotsym,3,/fill
  c=round(colors2[closest(mkh,colors)])
  oplot,[x2],[z2],psym=8,color=cgcolor(ocolor1),symsize=size1
  oplot,[x2],[z2],psym=8,color=cgcolor(ocolor2),symsize=size2
  oplot,[x2],[z2],psym=8,color=c,symsize=1.0
  xyouts,-45,23,'Mass (Solar)',charsize=charsize-0.7,charthick=thick,alignment=0.5
  sharpcorners,thick=thick
  
  plotsym,0,/fill
  cgloadct,0
  plot,y,z,psym=8,/xstyle,/ystyle,xthick=4,ythick=4,charsize=charsize,charthick=thick,xtitle='Y (pc)',ytitle='Z (pc)',symsize=symsize
  ct = 3
  clip = [25,225]
  cgloadct,ct,clip=clip
  colors = generatearray(min(masses),max(masses),100)
  vmin = 0
  vmax = 255
  colors2 = generatearray(vmin,vmax,100)
  for i = 0,n_elements(x)-1 do begin
     c=round(colors2[closest(mk[i],colors)])
     oplot,[y[i]],[z[i]],psym=8,color=c,symsize=symsize
  endfor
  plotsym,3,/fill
  c=round(colors2[closest(mkh,colors)])
  oplot,[y2],[z2],psym=8,color=cgcolor(ocolor1),symsize=size1
  oplot,[y2],[z2],psym=8,color=cgcolor(ocolor2),symsize=size2
  oplot,[y2],[z2],psym=8,color=c,symsize=1.0
  sharpcorners,thick=thick


  cgloadct,0
  cgloadct,ct,clip=clip
  position[0] = 0.07
  position[2] = !x.window[1]
  position[1] = !Y.window[1]+0.01
  position[3] = !y.window[1]+0.04
  locs = masses
  cgcolorbar,ticknames=strtrim(string(masses,format="(D3.1)"),2),divisions=4,xtickv=round([interpol(colors2,colors,locs[0]),interpol(colors2,colors,locs[1]),interpol(colors2,colors,locs[2]),interpol(colors2,colors,locs[3]),interpol(colors2,colors,locs[4])]),/top,position=position,title=' ',ticklen=-0.25,charsize=2.0,AnnotateColor=cgcolor('black',-100),charthick=thick
  sharpcorners,thick=thick

  device,/close
  stop
  
END

; $Id: //depot/idl/IDL_70/idldir/examples/doc/objects/orb__define.pro#1 $
;
; Copyright (c) 1997-2007, ITT Visual Information Solutions. All
;       rights reserved.
;+
; NAME:
;	ORB
;
; PURPOSE:
;	This object serves as a graphical representation of an orb,
;	(or sphere), which subclasses from the IDLgrModel class.
;
; CATEGORY:
;	Object graphics.
;
; CALLING SEQUENCE:
;	To initially create:
;	       	oOrb = OBJ_NEW('orb') 
;
;	To retrieve a property value:
;		oOrb->GetProperty
;
;	To set a property value:
;		oOrb->SetProperty
;
;	To print to the standard output stream the current properties of 
;	the orb:
;		oOrb->Print
;
;	To destroy:
;		OBJ_DESTROY, oOrb
;
; KEYWORD PARAMETERS:
;   ORB::INIT:
;	<Note that keywords accepted by IDLgrModel::Init and/or
;	 IDLgrPolygon::Init are also accepted here.>
;	POS:	A three-element vector, [x,y,z], specifying the position
;               of the center of the orb, measured in data units . 
;		Defaults to [0,0,0].
;	RADIUS: A floating point number representing the radius of the
;               orb (measured in data units).  The default is 1.0.
;	DENSITY: A floating point number representing the density at which
;               the vertices should be generated along the surface of the
;               orb.  The default is 1.0.
;	TEX_COORDS: Set this keyword to a nonzero value if texture map
;               coordinates are to be generated for the orb.
;
;   ORB::GETPROPERTY:
;	POS:	Set this keyword to a named variable that upon return will
;		contain a three-element vector, [x,y,z], specifying the 
;		position of the center of the orb, measured in data units . 
;	RADIUS: Set this keyword to a named variable that upon return will
;		contain a floating point number representing the radius of the
;               orb (measured in data units).
;	DENSITY: Set this keyword to a named variable that upon return will
;		contain a floating point number representing the density at 
;		which the vertices are generated along the surface of the
;               orb.
;
;   ORB::SETPROPERTY:
;	<Note that keywords accepted by IDLgrModel::SetProperty and/or
;	 IDLgrPolygon::SetProperty are also accepted here.>
;	POS:	A three-element vector, [x,y,z], specifying the position
;               of the center of the orb. Defaults to [0,0,0].
;	RADIUS: A floating point number representing the radius of the
;               orb (measured in data units).  The default is 1.0.
;	DENSITY: A floating point number representing the density at which
;               the vertices should be generated along the surface of the
;               orb.  The default is 1.0.
;
; EXAMPLE:
;	Create an orb centered at the origin with a radius of 0.5:
;		oOrb = OBJ_NEW('Orb', POS=[0,0,0], RADIUS=0.5) 
;
; MODIFICATION HISTORY:
; 	Written by:	RF, September 1996.
;-

;----------------------------------------------------------------------------
; ORB::INIT
;
; Purpose:
;  Initializes an orb object.
;
;  This function returns a 1 if initialization is successful, or 0 otherwise.
;
FUNCTION Orb::Init, POS=pos, RADIUS=radius, DENSITY=density, $
                       TEX_COORDS=tex_coords, _EXTRA=e

    IF (self->IDLgrModel::Init(_EXTRA=e) NE 1) THEN RETURN, 0

    self.pos = [0.0,0.0,0.0]
    self.radius = 1.0
    self.density = 1.0

    IF (N_ELEMENTS(pos) EQ 3) THEN $
        self.pos = pos

    IF (N_ELEMENTS(radius) EQ 1) THEN $
        self.radius = radius

    IF (N_ELEMENTS(density) EQ 1) THEN $
        self.density = density

    IF (N_ELEMENTS(tex_coords) EQ 1) THEN $
        self.texture = tex_coords

    ; Initialize the polygon object that will be used to represent
    ; the orb.
    self.oPoly = OBJ_NEW('IDLgrPolygon', SHADING=1, /REJECT, _EXTRA=e)

    self->Add,self.oPoly

    ; Build the polygon vertices and connectivity based on property settings.
    self->BuildPoly

    RETURN, 1
END

;----------------------------------------------------------------------------
; ORB::CLEANUP
;
; Purpose:
;  Cleans up all memory associated with the orb.
;
PRO orb::Cleanup

    ; Cleanup the polygon object used to represent the orb.
    OBJ_DESTROY, self.oPoly

    ; Cleanup the superclass.
    self->IDLgrModel::Cleanup
END

;----------------------------------------------------------------------------
; ORB::SETPROPERTY
;
; Purpose:
;  Sets the value of properties associated with the orb object.
;
PRO orb::SetProperty, $
    DENSITY=density, $
    PARENT=parent, $ ; Pass along to IDLgrModel only.
    POS=pos, $
    RADIUS=radius, $
    _EXTRA=e

    ; Pass along extraneous keywords to the superclass and/or to the
    ; polygon used to represent the orb.
    self->IDLgrModel::SetProperty, _EXTRA=e
    self.oPoly->SetProperty, _EXTRA=e

    IF (N_ELEMENTS(pos) EQ 3) THEN $
        self.pos = pos

    IF (N_ELEMENTS(radius) EQ 1) THEN $
        self.radius = radius

    IF (N_ELEMENTS(density) EQ 1) THEN $
        self.density = density

    ; Rebuild the polygon according to keyword settings.
    self->BuildPoly
END

;----------------------------------------------------------------------------
; ORB::GETPROPERTY
;
; Purpose:
;  Retrieves the value of properties associated with the orb object.
;
PRO orb::GetProperty, POS=pos, RADIUS=radius, DENSITY=density,$
                         POBJ=pobj, _REF_EXTRA=re

    ; Retrieve extra properties from polygon first, then model
    ; so that the model settings (for common keywords) will prevail.
    self.oPoly->GetProperty, _EXTRA=re
    self->IDLgrModel::GetProperty, _EXTRA=re

    pos = self.pos
    radius = self.radius 
    density = self.density 
    pobj = self.oPoly
END

PRO orb::Print
    PRINT, self.pos
    PRINT, self.radius
    PRINT, self.density
END

;----------------------------------------------------------------------------
; ORB::BUILDPOLY
;
; Purpose:
;  Sets the vertex and connectivity arrays for the polygon used to
;  represent the orb.
;
PRO orb::BuildPoly
    ; Build the orb.

    ; Number of rows and columns of vertices is based upon the density
    ; property.
    nrows = LONG(20.0*self.density)
    ncols = LONG(20.0*self.density)
    IF (nrows LT 2) THEN nrows = 2
    IF (ncols LT 2) THEN ncols = 2

    ; Create the vertex list and the connectivity array.
    nverts = nrows*ncols + 2
    nconn = (ncols*(nrows-1)*5) + (2*ncols*4)
    conn = LONARR(ncols*(nrows-1)*5 + 2*ncols*4)
    verts = FLTARR(3, nverts)
    IF (self.texture NE 0) THEN $
        tex = FLTARR(2,nverts)

    ; Fill in the vertices.
    i = 0L
    j = 0L
    k = 0L
    tzinc = (!PI)/FLOAT(nrows+1)
    tz = (!PI/2.0) - tzinc 
    FOR k=0,(nrows-1) DO BEGIN
        z = SIN(tz)*self.radius
        r = COS(tz)*self.radius
        t = 0
        IF (self.texture NE 0) THEN BEGIN
            tinc = (2.*!PI)/FLOAT(ncols-1)
        ENDIF ELSE BEGIN
            tinc = (2.*!PI)/FLOAT(ncols)
        ENDELSE
        FOR j=0,(ncols-1) DO BEGIN
            verts[0,i] = r*COS(t) + self.pos[0]
            verts[1,i] = r*SIN(t) + self.pos[1]
            verts[2,i] = z + self.pos[2]

            IF (self.texture NE 0) THEN BEGIN
                tex[0,i] = t/(2.*!PI)
                tex[1,i] = (tz+(!PI/2.0))/(!PI)
            ENDIF

            t = t + tinc
            i = i + 1L
        ENDFOR
        tz = tz - tzinc
    ENDFOR
    top = i
    verts[0,i] = self.pos[0]
    verts[1,i] = self.pos[1]
    verts[2,i] = self.radius*1.0 + self.pos[2]
    i = i + 1L
    bot = i
    verts[0,i] = self.pos[0]
    verts[1,i] = self.pos[1]
    verts[2,i] = self.radius*(-1.0) + self.pos[2]

    IF (self.texture NE 0) THEN BEGIN
        tex[0,i] = 0.5
        tex[1,i] = 0.0
        tex[0,i-1] = 0.5
        tex[1,i-1] = 1.0
    ENDIF

    ; Fill in the connectivity array.
    i = 0
    FOR k=0,(nrows-2) DO BEGIN
        FOR j=0,(ncols-1) DO BEGIN
            conn[i] = 4

            conn[i+4] = k*ncols + j

            w = k*ncols + j + 1L
            IF (j EQ (ncols-1)) THEN w = k*ncols
            conn[i+3] = w

            w = k*ncols + j + 1L + ncols
            IF (j EQ (ncols-1)) THEN w = k*ncols + ncols
            conn[i+2] = w

            conn[i+1] = k*ncols + j + ncols

            i = i + 5L
            IF ((self.texture NE 0) AND (j EQ (ncols-1))) THEN $
                i = i - 5L
        ENDFOR
    ENDFOR
    FOR j=0,(ncols-1) DO BEGIN
        conn[i] = 3
        conn[i+3] = top
        conn[i+2] = j+1L
        IF (j EQ (ncols-1)) THEN conn[i+2] = 0
        conn[i+1] = j
        i = i + 4L
        IF ((self.texture NE 0) AND (j EQ (ncols-1))) THEN $
            i = i - 4L
    ENDFOR
    FOR j=0,(ncols-1) DO BEGIN
        conn[i] = 3
        conn[i+3] = bot
        conn[i+2] = j+(nrows-1L)*ncols
        conn[i+1] = j+(nrows-1L)*ncols+1L
        IF (j EQ (ncols-1)) THEN conn[i+1] = (nrows-1L)*ncols
        i = i + 4L
        IF ((self.texture NE 0) AND (j EQ (ncols-1))) THEN $
            i = i - 4L
    ENDFOR

    self.oPoly->SetProperty, DATA=verts, POLYGONS=conn

    IF (self.texture NE 0) THEN $
        self.oPoly->SetProperty, TEXTURE_COORD=tex
END

;----------------------------------------------------------------------------
; ORB__DEFINE
;
; Purpose:
;  Defines the object structure for an orb object.
;
PRO orb__define
    struct = { orb, $
               INHERITS IDLgrModel, $
               pos: [0.0,0.0,0.0], $
               radius: 1.0, $
               density: 1.0, $
               texture: 0, $
               oPoly: OBJ_NEW() $
             }
END








;+
; NAME:
;       SCATTER_SURFACE
;
; PURPOSE:
;
;       The purpose of this program is to demonstrate how to
;       create a simple scatter surface plot with axes and rotational
;       capability in object graphics.
;
; AUTHOR:
;
;       FANNING SOFTWARE CONSULTING
;       David Fanning, Ph.D.
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: davidf@dfanning.com
;       Coyote's Guide to IDL Programming: http://www.dfanning.com
;
; CATEGORY:
;
;       Widgets, Object Graphics.
;
; CALLING SEQUENCE:
;
;       Scatter_Surface, x, y, z
;
; REQUIRED INPUTS:
;
;       None. Fake data will be used if no data is supplied in call. Otherwise,
;       pass three vectors of the same length, representing the X, Y, and Z values
;       of the data.
;
; OPTIONAL KEYWORD PARAMETERS:
;
;       EXACT:  Set this keyword to get exact axis scaling.
;
;       _EXTRA: This keyword collects otherwise undefined keywords that are
;        passed to the surface initialization routine.
;
;       GROUP_LEADER: The group leader for this program. When the group leader
;       is destroyed, this program will be destroyed.
;
;       LANDSCAPE: Set this keyword if you are printing in landscape mode. The
;       default is Portrait mode. The Landscape keyword on the PRINTER object
;       is set, but not all printers will honor this keyword setting. If yours
;       does not, set Landscape mode in the Printer Setup dialog.
;
;       VECTOR: Set this keyword if you want vector printing (as opposed to
;       the default bitmap printing).
;
;       XTITLE: A string used as the X title of the plot.
;
;       YTITLE: A string used as the Y title of the plot.
;
;       ZTITLE: A string used as the Z title of the plot.
;
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;       To use this program with your 2D data, type:
;
;        IDL> Scatter_Surface, vertices, polygons
;-


FUNCTION Normalize, range, Position=position

    ; This is a utility routine to calculate the scaling vector
    ; required to position a vector of specified range at a
    ; specific position given in normalized coordinates. The
    ; scaling vector is given as a two-element array like this:
    ;
    ;   scalingVector = [translationFactor, scalingFactor]
    ;
    ; The scaling vector should be used with the [XYZ]COORD_CONV
    ; keywords of a graphics object or model. For example, if you
    ; wanted to scale an X axis into the data range of -0.5 to 0.5,
    ; you might type something like this:
    ;
    ;   xAxis->GetProperty, Range=xRange
    ;   xScale = Normalize(xRange, Position=[-0.5, 0.5])
    ;   xAxis, XCoord_Conv=xScale

On_Error, 1
IF N_Params() EQ 0 THEN Message, 'Please pass range vector as argument.'

IF (N_Elements(position) EQ 0) THEN position = [0.0, 1.0] ELSE $
    position=Float(position)
range = Float(range)

scale = [((position[0]*range[1])-(position[1]*range[0])) / $
    (range[1]-range[0]), (position[1]-position[0])/(range[1]-range[0])]

RETURN, scale
END
;-------------------------------------------------------------------------



Pro Scatter_Surface_Cleanup, tlb

    ; Come here when program dies. Free all created objects.

Widget_Control, tlb, Get_UValue=info
IF N_Elements(info) NE 0 THEN Obj_Destroy, info.thisContainer
Obj_Destroy, info.thisPrinter, info.thisWindow, info.thisPolyline
Obj_Destroy, info.thisTrackball, info.thisModel
Obj_Destroy, info.xaxis, info.yaxis, info.zaxis

END
;-------------------------------------------------------------------



PRO Scatter_Surface_Draw_Events, event

     ; Draw widget events handled here: expose events and trackball
     ; events. The trackball uses RSI-supplied TRACKBALL oject.

Widget_Control, event.top, Get_UValue=info, /No_Copy

drawTypes = ['PRESS', 'RELEASE', 'MOTION', 'SCROLL', 'EXPOSE']
thisEvent = drawTypes(event.type)

CASE thisEvent OF

   'EXPOSE':  ; Nothing required except to draw the view.
   'PRESS': BEGIN
       Widget_Control, event.id, Draw_Motion_Events=1 ; Motion events ON.
       needUpdate = info.thisTrackball->Update(event, Transform=thisTransform)
       IF needUpdate THEN BEGIN
          info.thisModel->GetProperty, Transform=modelTransform
          info.thisModel->SetProperty, Transform=modelTransform # thisTransform
       ENDIF
       END
   'RELEASE': BEGIN
       Widget_Control, event.id, Draw_Motion_Events=0 ; Motion events OFF.
       needUpdate = info.thisTrackball->Update(event, Transform=thisTransform)
       IF needUpdate THEN BEGIN
          info.thisModel->GetProperty, Transform=modelTransform
          info.thisModel->SetProperty, Transform=modelTransform # thisTransform
       ENDIF
       END
   'MOTION': BEGIN ; Trackball events
       needUpdate = info.thisTrackball->Update(event, Transform=thisTransform)
       IF needUpdate THEN BEGIN
          info.thisModel->GetProperty, Transform=modelTransform
          info.thisModel->SetProperty, Transform=modelTransform # thisTransform
       ENDIF
       END
   'SCROLL': ; Nothing required except to draw the view.
ENDCASE

    ; Draw the view.

info.thisWindow->Draw, info.thisView

    ;Put the info structure back.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END
;-------------------------------------------------------------------



PRO Scatter_Surface_Style, event

     ; Event handler to select surface style.

Widget_Control, event.top, Get_UValue=info, /No_Copy

    ; What style is wanted?

Widget_Control, event.id, Get_UValue=newStyle
CASE newStyle OF

   'DOTS': info.thisPolyline->SetProperty, Style=0
   'MESH': info.thisPolyline->SetProperty, Style=1
   'SOLID': info.thisPolyline->SetProperty, Style=2, Shading=1
   'HIDDEN': BEGIN
       Widget_Control, event.id, Get_Value=buttonValue
       IF buttonValue EQ 'Hidden Lines OFF' THEN BEGIN
          setting = 0
          hlvalue = 'Hidden Lines ON'
       ENDIF ELSE BEGIN
          setting = 1
          hlvalue = 'Hidden Lines OFF'
       ENDELSE
       Widget_Control, event.id, Set_Value=hlvalue
       info.thisPolyline->SetProperty, Hidden_Lines=setting
       ENDCASE

ENDCASE

    ; Redraw the graphic.

info.thisWindow->Draw, info.thisView

    ;Put the info structure back.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END
;-------------------------------------------------------------------



PRO Scatter_Surface_Output, event

   ; This event handler creates GIF and JPEG files.

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; Get a snapshop of window contents. (TVRD equivalent.)

info.thisWindow->GetProperty, Image_Data=snapshot

   ; JPEG or GIF file wanted?

Widget_Control, event.id, GET_UValue=whichFileType
CASE whichFileType OF

   'GIF': BEGIN

         ; Because we are using RGB color for this model, we have
         ; a 3-m-n array. Use Color_Quan to create a 2D image and
         ; appropriate color tables for the GIF file.

      image2D = Color_Quan(snapshot, 1, r, g, b)
      filename = Dialog_Pickfile(/Write, File='idl.gif')
      IF filename NE '' THEN Write_GIF, filename, image2d, r, g, b
      END

   'JPEG': BEGIN

      filename = Dialog_Pickfile(/Write, File='idl.jpg')
      IF filename NE '' THEN Write_JPEG, filename, snapshot, True=1
      END

ENDCASE

    ;Put the info structure back.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END
;-------------------------------------------------------------------


PRO Scatter_Surface_Exit, event

   ; Exit the program. This will cause the CLEANUP
   ; routine to be called automatically.

Widget_Control, event.top, /Destroy
END
;-------------------------------------------------------------------



PRO Scatter_Surface_Printing, event

   ; PostScript printing and printer setup handled here.

Widget_Control, event.top, Get_UValue=info, /No_Copy

   ; Configure printer and print if user OKs.

result = Dialog_PrinterSetup(info.thisPrinter)
IF result EQ 1 THEN BEGIN

      ; Background colors can use a lot of toner. Change background
      ; color to white and axes to black before printing.

   info.xaxis->GetProperty, Color=axisColor
   info.thisView->GetProperty, Color=backgroundColor
   info.thisPolyline->GetProperty, Color=surfaceColor

   info.xaxis->SetProperty, Color=[0,0,0]
   info.yaxis->SetProperty, Color=[0,0,0]
   info.zaxis->SetProperty, Color=[0,0,0]
   info.thisView->SetProperty, Color=[255, 255, 255]
   info.thisPolyline->SetProperty, Color=[70,70,70]

      ; I want the output on the page to have the same aspect ratio
      ; as I see in the display window. I use the ASPECT function
      ; from the Coyote library to modify the view appropriately.
      ; Note that the "position" is returned in Normalized units (Units=3).

   info.thisWindow->GetProperty, Dimensions=wdims
   plotAspect = Float(wdims[1]) / wdims[0]
   info.thisPrinter->GetProperty, Dimensions=pdims
   windowAspect = Float(pdims[1]) / pdims[0]
   position = cgAspect(plotAspect, WindowAspect=windowAspect, Margin=0)
   info.thisView->SetProperty, Dimensions=[position[2]-position[0], position[3]-position[1]], $
      Location=[position[0], position[1]], Units=3

      ; Print the document.

   Widget_Control, Hourglass=1
   info.thisPrinter->Draw, info.thisView, Vector=info.vector
   info.thisPrinter->NewDocument
    Widget_Control, Hourglass=1

      ; Set colors and the view back to original values.

   info.xaxis->SetProperty, Color=axisColor
   info.yaxis->SetProperty, Color=axisColor
   info.zaxis->SetProperty, Color=axisColor
   info.thisView->SetProperty, Color=backgroundColor, Location=[0,0], Dimensions=[0,0]
   info.thisPolyline->SetProperty, Color=surfaceColor

ENDIF

   ; Put the info structure back.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END
;-------------------------------------------------------------------



PRO Scatter_Surface_Resize, event

     ; The only events generated by this simple program are resize
     ; events, which are handled here.

     ; Get the info structure.

Widget_Control, event.top, Get_UValue=info, /No_Copy

    ; Resize the draw widget.

info.thisWindow->SetProperty, Dimension=[event.x, event.y]

    ; Redisplay the graphic.

info.thisWindow->Draw, info.thisView

    ; Update the trackball objects location in the center of the
    ; window.

info.thisTrackball->Reset, [event.x/2, event.y/2], $
    (event.y/2) < (event.x/2)

    ;Put the info structure back.

Widget_Control, event.top, Set_UValue=info, /No_Copy
END
;-------------------------------------------------------------------



PRO Scatter_Surface, x, y, z, zcolors, _Extra=extra, XTitle=xtitle, $
   YTitle=ytitle, ZTitle=ztitle, Group_Leader=groupLeader, $
   Hidden_Lines=hidden_lines, Vector=vector, Exact=exact, $
   Landscape=landscape

    ; Check for input parameters.

IF N_Elements(x) EQ 0 OR N_Elements(y) EQ 0 OR N_Elements(z) EQ 0 THEN BEGIN

   ; Create the random data. Set the seed so you see what I see.

   ; Create the random data. Set the seed so you see what I see.

    seed = 1L
    x = RANDOMU(seed, 32)
    y = RANDOMU(seed, 32)
    z = EXP(-3 * ((x - 0.5)^2 + (y - 0.5)^2))

ENDIF
IF N_Elements(zcolors) EQ 0 THEN zcolors = BYTSCL(z)

    ; Check for keywords.

IF N_Elements(xtitle) EQ 0 THEN xtitle='X Axis'
IF N_Elements(ytitle) EQ 0 THEN ytitle='Y Axis'
IF N_Elements(ztitle) EQ 0 THEN ztitle='Z Axis'
hidden_lines = Keyword_Set(hidden_lines)
landscape = Keyword_Set(landscape)
vector = Keyword_Set(vector)


    ; Create a view. Use RGB color. Charcoal background.
    ; The coodinate system is chosen so that (0,0,0) is in the
    ; center of the window. This will make rotations easier.

thisView = OBJ_NEW('IDLgrView', Color=[255,255,255], $
   Viewplane_Rect=[-1.2,-1.1,2.3,2.1])

    ; Create a model for the surface and axes and add it to the view.
    ; This model will rotate under the direction of the trackball object.

thisModel = OBJ_NEW('IDLgrModel')
thisView->Add, thisModel

    ; Create helper objects. First, create title objects
    ; for the axes and plot. Color them green.

xTitleObj = Obj_New('IDLgrText', xtitle, Color=[0,255,0])
yTitleObj = Obj_New('IDLgrText', ytitle, Color=[0,255,0])
zTitleObj = Obj_New('IDLgrText', ztitle, Color=[0,255,0])

    ; Create font objects.

helvetica10pt = Obj_New('IDLgrFont', 'Helvetica', Size=10)
helvetica14pt = Obj_New('IDLgrFont', 'Helvetica', Size=14)

    ; Create a trackball for surface rotations. Center it in
    ; the 400-by-400 window. Give it a 200 pixel diameter.

thisTrackball = OBJ_NEW('Trackball', [200, 200], 200)

   ; Create a color palette for coloring the symbols.

thisPalette = Obj_New('IDLgrPalette')
thisPalette->LoadCT, 33
thisPalette->GetProperty, Red=r, Green=g, Blue=b
Obj_Destroy, thisPalette

   ; Create the symbols for each point.

npts = N_Elements(x)
theseSymbols=ObjArr(npts)
theselines = ObjArr(npts)
FOR j=0,npts-1 DO BEGIN
   orb = Obj_New('Orb', radius=0.1, $
        COLOR=[r[zcolors[j]], g[zcolors[j]], b[zcolors[j]]])
   theseSymbols[j] = Obj_New('IDLgrSymbol', orb, size=0.25)
   theseLines[j] = Obj_New('IDLgrPolyLine', [x[j], x[j]], [y[j], y[j]], [z[j], Min(z)], $
        COLOR=[r[zcolors[j]], g[zcolors[j]], b[zcolors[j]]])
ENDFOR

    ; Create Polyline object..

thisPolyline = OBJ_NEW('IDLgrPolyline', x, y, z, $
   LineStyle=6, SYMBOL=theseSymbols)

    ; Get the data ranges of the surface.

thisPolyline->GetProperty, XRange=xrange, YRange=yrange, ZRange=zrange

    ; Create axes objects for the surface. Color them green.
    ; Axes are created after the surface so the range can be
    ; set correctly. Note how I set the font to 10 point Helvetica
    ; by creating the axis with the title object, then getting the
    ; actual axis text from the axis object itself, and switching it.
    ; Set the RECOMPUTE_DIMENSIONS keyword on the axis text objects
    ; so the text doesn't go crazy when we change the data range.

xAxis = Obj_New("IDLgrAxis", 0, Color=[0,0,0], Ticklen=0.1, $
   Minor=4, Title=xtitleObj, Range=xrange, Exact=Keyword_Set(exact))
xAxis->GetProperty, Ticktext=xAxisText
xAxisText->SetProperty, Font=helvetica10pt, Recompute_Dimensions=2

yAxis = Obj_New("IDLgrAxis", 1, Color=[0,0,0], Ticklen=0.1, $
   Minor=4, Title=ytitleObj, Range=yrange, Exact=Keyword_Set(exact))
yAxis->GetProperty, Ticktext=yAxisText
yAxisText->SetProperty, Font=helvetica10pt, Recompute_Dimensions=2

zAxis = Obj_New("IDLgrAxis", 2, Color=[0,0,0], Ticklen=0.1, $
   Minor=4, Title=ztitleObj, Range=zrange, Exact=Keyword_Set(exact))
zAxis->GetProperty, Ticktext=zAxisText
zAxisText->SetProperty, Font=helvetica10pt, Recompute_Dimensions=2

    ; The axes may not use exact axis scaling, so the ranges may
    ; have changed from what they were originally set to. Get
    ; and update the range variables.

xAxis->GetProperty, CRange=xrange
yAxis->GetProperty, CRange=yrange
zAxis->GetProperty, CRange=zrange

    ; Set scaling parameters for the surface and axes so that everything
    ; is scaled into the range -0.5 to 0.5. We do this so that when the
    ; surface is rotated we don't have to worry about translations. In
    ; other words, the rotations occur about the point (0,0,0).

xs = Normalize(xrange, Position=[-0.5,0.5])
ys = Normalize(yrange, Position=[-0.5,0.5])
zs = Normalize(zrange, Position=[-0.5,0.5])

    ; Scale the axes and place them in the coordinate space.
    ; Note that not all values in the Location keyword are
    ; used. (I've put really large values into the positions
    ; that are not being used to demonstate this.) For
    ; example, with the X axis only the Y and Z locations are used.

xAxis->SetProperty, Location=[9999.0, -0.5, -0.5], XCoord_Conv=xs
yAxis->SetProperty, Location=[-0.5, 9999.0, -0.5], YCoord_Conv=ys
zAxis->SetProperty, Location=[-0.5,  0.5, 9999.0], ZCoord_Conv=zs

    ; Scale the surface. Notice the surface is scaled *AFTER* the
    ; actual data range is obtained from the axes (CRANGE, above).
    ; Failure to do this can result in inaccurate results.

thisPolyline->SetProperty, XCoord_Conv=xs, YCoord_Conv=ys, ZCoord_Conv=zs
FOR j=0,npts-1 DO theseLines[j] -> SetProperty, XCoord_Conv=xs, YCoord_Conv=ys, ZCoord_Conv=zs

    ; Add the surface and axes objects to the model.

thisModel->Add, thisPolyline
thisModel -> Add, theseLines
thisModel->Add, xAxis
thisModel->Add, yAxis
thisModel->Add, zAxis

    ; Rotate the surface model to the standard surface view.

thisModel->Rotate,[1,0,0], -90  ; To get the Z-axis vertical.
thisModel->Rotate,[0,1,0],  30  ; Rotate it slightly to the right.
thisModel->Rotate,[1,0,0],  30  ; Rotate it down slightly.

    ; Create the widgets to view the surface. Set expose events
    ; on the draw widget so that it refreshes itself whenever necessary.
    ; Button events are on to enable trackball movement.

tlb = Widget_Base(Title='Resizeable Window Surface Example', Column=1, $
   TLB_Size_Events=1, MBar=menubase)
drawID = Widget_Draw(tlb, XSize=400, YSize=400, Graphics_Level=2, Retain=0, $
   Expose_Events=1, Event_Pro='Scatter_Surface_Draw_Events', Button_Events=1)

    ; Create FILE menu buttons for printing and exiting.

fileID = Widget_Button(menubase, Value='File')

outputID = Widget_Button(fileID, Value='Save As...', /Menu)
dummy = Widget_Button(outputID, Value='GIF File', $
   UValue='GIF', Event_Pro='Scatter_Surface_Output')
dummy = Widget_Button(outputID, Value='JPEG File', $
   UValue='JPEG', Event_Pro='Scatter_Surface_Output')

dummy = Widget_Button(fileID, Value='Print', Event_Pro='Scatter_Surface_Printing')

dummy = Widget_Button(fileID, /Separator, Value='Exit', $
   Event_Pro='Scatter_Surface_Exit')

    ; Shaded surfaces will not look shaded unless there is a light
    ; source. Create a positional light source for this surface.

thisLight = Obj_New('IDLgrLight', Type=1, $
    Location=[xrange[1], yrange[1], 4*zrange[1]], $
    Direction=[xrange[0], yrange[0], zrange[0]])
thisModel->Add, thisLight

    ; Create a fill light source so you can see the underside
    ; of the surface. Otherwise, just the top surface will be visible.

fillLight = Obj_New('IDLgrLight', Type=1, Intensity=0.5, $
   Location=[(xrange[1]-xrange[0])/2.0, (yrange[1]-yrange[0])/2.0, -2*Abs(zrange[0])], $
   Direction=[(xrange[1]-xrange[0])/2.0, (yrange[1]-yrange[0])/2.0, zrange[1]])
thisModel->Add, fillLight

    ; Scale the light sources.

thisLight->SetProperty, XCoord_Conv=xs, YCoord_Conv=ys, ZCoord_Conv=zs
fillLight->SetProperty, XCoord_Conv=xs, YCoord_Conv=ys, ZCoord_Conv=zs

   ; Realize the widgets.

Widget_Control, tlb, /Realize

    ; Get the window destination object, which is the value of
    ; an object draw widget. The view will be drawn in the window
    ; when the window is exposed.

Widget_Control, drawID, Get_Value=thisWindow

   ; Get a printer object for this graphic.

thisPrinter = Obj_New('IDLgrPrinter', Print_Quality=2, Landscape=landscape)

   ; Create a container object to hold all the other
   ; objects. This will make it easy to free all the
   ; objects when we are finished with the program.

thisContainer = Obj_New('IDL_Container')

   ; Add created objects to the container. No need to add objects
   ; that have been added to the model, since a model object is
   ; a subclass of a container object. But helper objects that
   ; are NOT added to the model directly MUST be destroyed properly.

thisContainer->Add, thisView
thisContainer->Add, thisPrinter
thisContainer->Add, thisTrackball
thisContainer->Add, xTitleObj
thisContainer->Add, yTitleObj
thisContainer->Add, zTitleObj
thisContainer->Add, thisModel
thisContainer->Add, helvetica10pt
thisContainer->Add, helvetica14pt
thisContainer->Add, theseSymbols
thisContainer->Add, theseLines

   ; Create an INFO structure to hold needed program information.

info = { thisContainer:thisContainer, $ ; The object container.
         thisWindow:thisWindow, $       ; The window object.
         thisPrinter:thisPrinter, $     ; The printer object.
         thisPolyline:thisPolyline, $     ; The surface object.
         thisTrackball:thisTrackball, $ ; The trackball object.
         thisModel:thisModel, $         ; The model object.
         thisView:thisView, $           ; The view object.
         xaxis:xaxis, $                 ; The X axis object.
         yaxis:yaxis, $                 ; The Y axis object.
         zaxis:zaxis, $                 ; The Z axis object.
         landscape:landscape, $         ; A flag for landscape output.
         vector:vector }                ; A flag for vector output.

   ; Store the info structure in the UValue of the TLB.

Widget_Control, tlb, Set_UValue=info, /No_Copy

   ; Call XManager. Set a cleanup routine so the objects
   ; can be freed upon exit from this program.

XManager, 'Scatter_Surface', tlb, Cleanup='Scatter_Surface_Cleanup', /No_Block, $
   Event_Handler='Scatter_Surface_Resize', Group_Leader=groupLeader
END
;-------------------------------------------------------------------
