PRO param_c13

;; kinematic distance:
  dist = 59.4;58.12
  dist_err = 2.8;2.6
;;       17.206060      0.73500748
  nmonte = 1d5
  dist = dist+dist_err*randomn(seed,nmonte)
  k = 8.368+0.019*randomn(seed,nmonte)
  h = 8.496+0.020*randomn(seed,nmonte)
  j = 9.096+0.022*randomn(seed,nmonte)

  mk = k-5.0*(alog10(dist)-1.)
  mh = h-5.0*(alog10(dist)-1.)
  mj = j-5.0*(alog10(dist)-1.)
  
  readcol,'BHAC15_800.txt',bmass,bteff,blum,blogg,brad,blith,bmj,bmh,bmk,format='d,d,d,d,d,d,d,d,d',/silent
  ;;readcol,'BHAC15_1000.txt',bmass,bteff,blum,blogg,brad,blith,bmj,bmh,bmk,format='d,d,d,d,d,d,d,d,d',/silent
  ;;readcol,'BHAC15_625.txt',bmass,bteff,blum,blogg,brad,blith,bmj,bmh,bmk,format='d,d,d,d,d,d,d,d,d',/silent
  ;;
  ;; restore,'mist.dat'
  ;; mist = mist[where(abs(mist.age-800) lt 50)]
  ;; bteff = mist.teff
  ;; bmass = mist.m
  ;; blum = mist.l
  ;; brad = mist.r
  ;; blogg = mist.g
  ;; bmj = mist.mj
  ;; bmh = mist.mh
  ;; bmk  = mist.mk
  ;mk = mk
  ;bmk = bmk
  
  blum = 10.0^blum
  mass = interpol(bmass,bmk,mk)
  radius = interpol(brad,bmk,mk)
  teff = interpol(bteff,bmk,mk)
  lum = interpol(blum,bmk,mk)
  den = mass/radius^3.
  print,'Mass='+string(median(mass),format="(D4.2)")+'+/-'+string(stdev(mass),format="(D4.2)")
  print,'Radius='+string(median(radius),format="(D4.2)")+'+/-'+string(stdev(radius),format="(D4.2)")
  ;;format_str,den,med,pl,mi
  ;;print,'Density='+string(med,format="(D4.2)")+'$^{'+string(pl,format="(D4.2)")+'}_{'+string(mi,format="(D4.2)")+'}$'
  print,'Teff='+string(median(teff),format="(I4)")+'+/-'+string(stdev(teff),format="(I3)")
  print,'Lum='+string(median(lum),format="(D5.3)")+'+/-'+string(stdev(lum),format="(D5.3)")
  rsun = 6.955d10
  msun = 1.989d33
  bigg = 6.6743d-8
  g = bigG*mass*msun/(radius*rsun)^2.
  logg = alog10(g)
  print,'Log(g)='+string(median(logg),format="(D4.2)")+'+/-'+string(stdev(logg),format="(D4.2)")

  fbol = 0.149d-12+0.005d-12*randomn(seed,nmonte)
  c1 = (4.*!pi * 1d4) / 3.865d33 ;
  c2 = 3.08567758d18
  c = c1*c2^2.0
  l = c*(fbol*dist^2.0)
  logl = alog10(l)
  print,'Lum from Fbol: '+string(median(l),format="(D5.3)")+'+/-'+string(stdev(l),format="(D5.3)")


  ;; print,'-----'
  ;; mass = interpol(bmass,blum,l)
  ;; radius = interpol(brad,blum,l)
  ;; teff = interpol(bteff,blum,l)
  ;; lum = interpol(blum,blum,l)
  ;; den = mass/radius^3.
  ;; print,'Mass='+string(median(mass),format="(D4.2)")+'+/-'+string(stdev(mass),format="(D4.2)")
  ;; print,'Radius='+string(median(radius),format="(D4.2)")+'+/-'+string(stdev(radius),format="(D4.2)")
  ;; format_str,den,med,pl,mi
  ;; print,'Density='+string(med,format="(D4.2)")+'$^{'+string(pl,format="(D4.2)")+'}_{'+string(mi,format="(D4.2)")+'}$'
  ;; print,'Teff='+string(median(teff),format="(I4)")+'+/-'+string(stdev(teff),format="(I3)")
  ;; print,'Lum='+string(median(lum),format="(D5.3)")+'+/-'+string(stdev(lum),format="(D5.3)")
  ;; rsun = 6.955d10
  ;; msun = 1.989d33
  ;; bigg = 6.6743d-8
  ;; g = bigG*mass*msun/(radius*rsun)^2.
  ;; logg = alog10(g)
  ;; print,'Log(g)='+string(median(logg),format="(D4.2)")+'+/-'+string(stdev(logg),format="(D4.2)")

  print,'------'
  ;; Teq
  a = 0.3
  a_r = [39.0,51.6,23.2] ;; 2.86, 1.44, 0.98
  a_r_err_pl = [7.2,5.2,2.4]
  a_r_err_mi = [1.5,4.5,1.2]
  for i = 0,2 do begin
     mult = randomn(seed,nmonte)
     mult[where(mult lt 0)]*=a_r_err_mi[i]
     mult[where(mult gt 0)]*=a_r_err_pl[i]
     ;;$39.0^{+7.2}_{-1.5}$ & $50.7^{+5.2}_{-4.5}$ & $23.1^{+2.4}_{-1.2}
     ar = a_r[i]+mult
     teq = teff*(1-a)^0.5*sqrt(0.5/ar)
     format_str,teq,med,pl,mi
     print,'Teq='+string(med,format="(I3)")+'$^{'+string(pl,format="(I3)")+'}_{'+string(mi,format="(I3)")+'}$'
  endfor
  stop

END



PRO format_str,arr,med,pl,mi

  tmp1 = arr[where(finite(arr) eq 1)]
  tmp2 = tmp1[sort(tmp1)]
  l = where(finite(tmp2) eq 1)
  tmp2 = tmp2[l]
  med = median(tmp2)
  pl = tmp2[n_elements(tmp2)*(0.5+0.68/2.)]-med
  mi = med-tmp2[n_elements(tmp2)*(0.5-0.68/2.)]

END
