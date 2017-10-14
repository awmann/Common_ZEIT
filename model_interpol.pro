PRO model_interpol,model

  nmonte = 1d5
  plx = 32.98+1.76*randomn(seed,nmonte)
  k =  7.193 + 0.024*randomn(seed,nmonte)
  mk  = k -5.0*(alog10(1000d0/plx)-1.0)
  feh = -0.2;-0.5
  age = 5000
  ;model = 'mist'
  ;model = 'dsep'
  ;print,model
  case model of
     'bhac':readcol,'BHAC15_1000.txt',bmass,bteff,blum,blogg,brad,blith,bmj,bmh,bmk,format='d,d,d,d,d,d,d,d,d',/silent
     'mist':begin
        restore,'~/Dropbox/Structures/mist.dat'
        l = where(abs(mist.age-age) lt 100 and abs(mist.feh-feh) lt 0.1)
        blum = mist[l].l
        brad = mist[l].r
        bmass = mist[l].m
        bmk = mist[l].mk
        bteff = mist[l].teff
     end
     'dsep':begin
        restore,'~/Desktop/dsep_solar.dat'
        r = (sqrt(6.6743d-8*dsep.m*1.989d33/(10d0^dsep.logg)))/6.955d10
        l = where(dsep.age eq 1.0 and abs(dsep.feh-feh) lt 0.21 and dsep.logg gt 4 and dsep.afe eq 0.0)
        blum = dsep[l].l
        brad = r[l]
        bmass = dsep[l].m
        bteff = dsep[l].teff
        bmk = dsep[l].k
        dsep = dsep[l]
     end
  endcase
  
  blum = 10.0^blum
  mass = interpol(bmass,bmk,mk)
  radius = interpol(brad,bmk,mk)
  teff = interpol(bteff,bmk,mk)
  lum = interpol(blum,bmk,mk)
  den = mass/radius^3.
  print,'Mass='+string(median(mass),format="(D4.2)")+'+/-'+string(stdev(mass),format="(D4.2)")
  print,'Radius='+string(median(radius),format="(D4.2)")+'+/-'+string(stdev(radius),format="(D4.2)")
  print,'Teff='+string(median(teff),format="(I4)")+'+/-'+string(stdev(teff),format="(I3)")
  print,'Lum='+string(median(lum),format="(D5.3)")+'+/-'+string(stdev(lum),format="(D5.3)")
  rsun = 6.955d10
  msun = 1.989d33
  bigg = 6.6743d-8
  g = bigG*mass*msun/(radius*rsun)^2.
  logg = alog10(g)
  print,'Log(g)='+string(median(logg),format="(D4.2)")+'+/-'+string(stdev(logg),format="(D4.2)")

  stop

END
