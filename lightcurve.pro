PRO wrapper

  epics = [211067768,204364515,204506777,210696763,205117205,210894022,211075914,$
           210974364,212002525,211972086,211082420,211093684,211990866,211923502,211913977,$
           210490365,210363145,211143504,210822691,210673842,210532473,210463428,$
           211970147,205046529,210651828,210659688,210931967,211102816,$
           211822797,211901114,203791768]

  ;;epics = [205117205]
  ;;epics = [211970147,211913977,211990866]
  ;;epics = [210994935,210832378,210803812,211041474,211047980]
  ;;epics = [212009427,211799258,212002525,211946007,211972086]
  ;;epics = [210363145,205117205,205046529,211901114,211822797,211969807,211916756,211970147,211913977,211990866]
                                ;epics = [211901114,211822797,211969807,211916756,211970147,211913977,211990866]
                                ;epics = [211038138,211989185] ;;,211989185,203758182]
                                ;epics = [211038138,210941195,211989185,203758182,204262368]
                                ;epics = [205046529]
                                ;epics = [210941195]

  epics = [247306074,247822311,247002634,247589423, 248045685,247556609,210743724,210923016,247556609,247639356,247002634,247306074,247950452,247738463]
  epics = epics[sort(epics)]
  epics = epics[uniq(epics)]
  ;;removed = 210709514,210718930,246816662,247623969
  ;;missing = 247239752,248131102
  epics = [247556609,203849738,202873945,204364515]
  for i = 0,n_elements(epics)-1 do begin
     ;;lightcurve,epics[i],type='A'
     lightcurve,epics[i],type='V'
     ;;lightcurve,epics[i],type='E'
     ;;lightcurve,epics[i],type='K'
     ;;f epics[i] gt 203000000d0 and epics[i] lt 206000000d0 then lightcurve,epics[i],type='C'
  endfor
  

END

PRO cand

  readcol,'planets.txt',epic,cluster,p,t0,format='d,a,d,d',/silent

  cluster = strtrim(cluster,2)
  for i = 0,n_elements(epic)-1 do begin
     smooth = 10
     if cluster[i] eq 'Pl' then smooth = 10
     if cluster[i] eq 'UsA' then smooth = 5
     
     ;;lightcurve,epic[i],period=p[i],t0=t0[i],smooth=smooth,type='V'
     ;;lightcurve,epic[i],period=p[i],t0=t0[i],smooth=smooth,type='A'
     ;;lightcurve,epic[i],period=p[i],t0=t0[i],smooth=smooth,type='E'
     stop
  endfor

END



PRO lightcurve,epic,print=print,$
               type=type,$
               datemin=datemin,$
               datemax=datemax,$
               period=period,$
               t0=t0,$
               smooth=smooth,$
               everyother=everyother  ;; V=vanderburg corrected, G=Gaidos corrected

  if n_elements(everyother) eq 0 then everyother = 0
  if n_elements(datemax) eq 0 then datemax = 1d10
  if n_elements(datemin) eq 0 then datemin = 0
  baset = 2454833d0             ;+2070d0
  plotsym,0,/fill
  dur = 0.17
  !y.margin=[4,1]
  fitpoly = 1
  debug = 0
  smooth = 30
  bins = 600
  fullrange = 0
  if n_elements(epic) eq 0 then epic = 205117205
  print,'EPIC: '+string(epic,format="(I10)")
  charsize = 3.0
  charthick = 5.0
  thick = 10.0
  ;;t0 = -1
  if n_elements(period) ne 0 then per = period
  if n_elements(per) eq 0 then begin
     case epic of
        246865365: begin
           per = 3.385841
           t0 = 2997.00104
           dur = 0.15
           ra = 78.33914
           dec = 16.4142
           smooth = 30
        end
        247002634: begin
           per = 1.5632
           t0 = 2987.21336649
           smooth = 10
           dur = 0.1
        end
        247589423: begin
           per = 17.3053994
           t0 = 2979.78904
           per = [17.3071156085d0,25.575132578d0,7.9753003744d0]
           per_err = [0.000250016500129d0,0.00236123372544d0,0.00077589788687d0]
           t0 = [2979.717720410d0, 2947.81148275d0,2984.75628769d0]
           t0_err = [0.0009587611721140d,0.00637929869208d0,0.0046197241677d0]     
           smooth = 20
           dur = [0.25,0.25,0.15]
           bins = 300
           ll = sort(per)
           per = per[ll]
           t0 = t0[ll]
           dur = dur[ll]
        end
        247611242: begin
           per = 24.1884697
           t0 = 2457853.299394d0-2454833d0
           dur = 0.1
           smooth = 10
        end
        248045685: begin
           per = 37.3381479
           t0 = 2457828.386980-2454833d0
        end
        247556609: begin ;;EB Hyades
           per = 13.8709
           t0 = 2983.4119
           dur = 0.2
        end
        210743724: begin
           per = 5.9462;;*2.0
           t0 = 2985.349657
           dur = 0.1
           smooth = 10
        end
        248131102: begin ;; Dipper?
           per = 1.83947514
           t0 = 2987.37458606
           dur = 0.1
           smooth = 15
        END
        247822311: begin
           per = 21.4792
           t0 = 3006.22
           dur = 0.2
           smooth = 8
        end
        247306074: begin
           per = 22.923
           t0 = 2284.68465
           dur = 0.2
           smooth = 15
        end
        247239752: begin ;; missing?
           per = 1.96057929
           t0 = 2987.55416
           dur = 0.1
           smooth = 15
        end
        ;; more questionable ones
        247639356: begin
           per = 1.6374
           t0 = 2987.32154745
           smooth = 15
           dur = 0.1
        end
        210709514: begin
           per = 30.631
           t0 = 2957.969
           dur = 0.2
           smooth = 20
        end
        210718930: begin
           per = 5.1617
           t0 = 2983.88399074
           dur = 0.1
           smooth = 20
        end
        210743724: begin
           per = 5.9462
           t0 = 2985.349657
           dur = 0.1
           smooth = 20
        end
        210923016: begin
           per = 4.637;4.6375
           t0 = 2984.0766
           dur = 0.12
           smooth = 5
        end

        247623969: begin
           per = 4.2353
           t0 = 2986.9685
           dur = 0.1
           smooth = 15
        end
        247738463: begin
           per = 32.9852
           t0 = 2909
           dur = 0.3
           smooth = 15
        end
        247950452: begin
           per = 9.8432
           t0 = 2987.0938
           dur = 0.15
           smooth = 15
        end
       


        ;; end new candidates
     
        211001695: begin
           per = 9.9940 
           t0 = 2221.304
           dur = 0.1
           smooth = 10
        end
        210736056: begin
           per = 29.053
           t0 = 57044.194729d0+2400000.5d0;;2211.195
           dur = 0.2
           smooth = 10
        end
        211795569: begin
           per = 6.5053
           t0 = 2306.228686
           dur = 0.1
           smooth = 20
        end
        211411955: begin
           per = 13.5355
           t0 = 2304.262574
           dur = 0.1
           smooth = 20
        end
        219800881: begin
           per = 13.85
           t0 = 2487.000002
           dur = 0.1
           smooth = 30
        end
        210827030: begin
           per = 1.0221/2.0
           t0 = 2231.05 ;2230.976451
           smooth = 10
           dur = 0.16
           print,t0
        end
        211810283: begin
           per = 29.1545
           t0 = 2301.185218
           smooth = 10
           dur = 0.5
        END
        203710077: begin
           per = 17.6554
           t0 = 2057.77422
           smooth = 5
        end
        203520293: begin
           per = 17.5871
           t0 = 2056.557162
           smooth = 5
        end
        203542463: begin
           per = 30.012
           t0 = 2047.2
           smooth = 5
        end
        210473042: begin
           per =12.258
           t0 = 2222.896
        end
        201500729: begin
           per = 12.4906
           t0 = 2221.99199
        end
        210829407: begin
           per = 19.6078
           t0 = 2222.9106
        end
        210861773: begin
           smooth = 20
           per = 9.992
           t0 = 2221.75505
        end
        210958195: begin
           smooth = 10
           per = 4.5434
           t0 = 2227.48123
        end
        211017469: begin
           per = 9.7333
           t0 = 2223.06644
        end
        211093117: begin
           per = 16.6223
           t0 = 2226.3104
        end
        212679181: begin
           per = 1.054486
           t0 = 2385.60731638
           smooth = 30
           dur = 0.1
        end
        204080947: begin
           per = 0.735717
           t0 = 2061.42877821
           smooth = 30
           dur = 0.1
        end
        205339436: begin
           per = 1.226188
           t0 = 2061.67393785
           smooth = 30
           dur = 0.1
        end
        205374937: begin
           per = 0.566
           t0 = -1
           range = [0.1,0.4]
           smooth = 50
        end
        204281213: begin
           per = 5.0
           t0 = -1
           smooth = 20
           range = [0,4.9]
        end
        203753577: begin
           per = 3.4007758
           t0 = 2062.094
           ;;range = [0,3]
           dur = 0.1
           smooth = 20
        end
        202925169: begin
           per = 13.6277
           t0 =  2249.42        ;-1
           smooth = 20
           dur = 0.2
                                ;range = [0,12]
        end
        211989185: begin
           per = 12.9433
           t0 = 2249.42
           ;;range = [10.6,11.4]
           dur = 0.15
           smooth = 20
        end
        210941195: begin
           per = 12.6582
           ;;t0 = -1
           ;;range = [9,12.]
           smooth = 30
           ;;     range[0] = (t0 mod per) - dur*1.5d0
           ;;     range[1] = (t0 mod per) + dur*1.5d0
           t0 = 2251.14
           dur = 0.2
        end
        204346718: begin
           per = 0.64292573
           t0 = 2063.42544
           smooth = 10
           dur = 5.977/24.0
           range = [0,0.5]
        end
        203560204: begin
           per = 29.33873041
           t0 = 2091.7063
           smooth = 10
           dur = 5.083/24d0
        end
        203850058: begin
           per = 2.88
           t0 = -1
           smooth = 1000
           range = [0,2.5]
        end
        211038138: begin
           per = 9.92
           t0 = 2089.42
           smooth = 20
           dur = 0.1
           ;;range = [6.0,6.4]
        end
        205151387: begin
           per = 9.55
           t0 = -1
           smooth = 1000
           range = [0,9.4]
        end
        204638512: begin
           per =5.0
           t0 = -1
           smooth = 1000
           range = [0,4.8]
        end
        202613991: begin
           per = 0.720
           t0 = -1
           dur = 0.1
           smooth = 20
           range = [0,0.65]
        end
        203758182: begin
           per = 17.7
           t0 = -1
           dur = 0.1
           smooth = 20
           range = [0.2,1.5]
        end
        202533810: begin
           per = 5.0
           t0 = -1
           smooth = 15
           dur = 0.2
           range = [0,4]
        end
        203849738: begin ;;EB
           per = 0.619898
           t0 = -1
           smooth = 10
           dur = 0.2
           range = [0,0.6]
        end
        202873945: begin ;;EB
           per = 0.625744
           t0 = 2061.1
           smooth = 10
           dur = 0.2
           range = [0,0.6]
        end
        204364515: begin ;;EB
           per = 1.45627
           t0 = 2061.298093
           smooth = 10
           dur = 0.2
           range = [0,1.3]
        end
        211082420: begin
           per = 1.23045
           t0 = 2229.74009535 ;;2231.14761
           smooth = 30
           dur = 0.1
        end
        211076026: begin
           t0 = 2240.67127169
           per = 15.410479
           smooth = 10
           dur = 0.1
        end
        211041474: begin
           t0 = 2245.84058319
           per = 18.481621
           smooth = 10
           dur = 0.1
        end
        211047980: begin
           per = 17.481
                                ;range = [3.4,4.2]
           t0 = 2241.40686747
           smooth = 10
           dur = 0.2
        end
        210994935: begin
           per = 14.634
           ;;range = [2.8,4.0]
           smooth = 10
           dur = 0.2
           t0 = 2242.51006433
        end
        210832378: begin
           per = 11.413
           ;;range = [2.1,2.4]
           t0 = 2239.22046906
           smooth = 10
           dur = 0.1
        end
        210803812: begin
           per = 11.412
                                ;range = [2.3,2.6]
           t0 = 2239.19997574
           smooth = 10
           dur = 0.1
        end
        211804579: begin
           range = [0,1]
           per = 1.5235
           dur = 0.1
           bins = 300
           smooth = 5
           t0 = 2307.1
           print,t0
        end
        211037886: begin
           per = 4.339512
           t0 = 2231.53851146
           smooth = 10
           dur = 0.03
           bins = 500
        end
        203837143: begin
           per = 16.237942	         
           t0 = 2070.94946841
           smooth = 10
           dur = 0.125
        end
        203848625: begin
           per = 10.057708
           range = [1.5,3.5]
           t0 = 2074.25964627
           smooth = 30
           dur = 0.3
        end
        203791768: begin
           per = 15.851775
                                ;range = [4,6]
           t0 = 2065.76832900
           fullrange = 0
        end
        211067768: begin
           per = 16.53034
           range = [1,5]
           smooth = 10
           fullrange=0
           t0=-1
        end
        204364515: begin
           per = 1.4559
           ;;range = [0,1.3]
           smooth = 20
           fullrange = 1
           t0 = 2061.59117
        end
        204506777: begin
           per = 0.81
           range = [0,0.7]
           fullrange = 1
           smooth = 15
           t0 = -1
        end
        210696763: begin
           t0 = -1
           per = 7.9844
           range = [7,9]
           per = 3.6510101
        end
        205117205: begin
           fitpoly = 2
           t0 = 2065.69282551   ;2065.69320392
           t0_err = 0.00165147421467
           per = 5.42486851247d0 ;5.42485318276d0
           per_err = 2.59760937651d-5
           ;;t0+=per/2.0
           smooth = 10
           fullrange=0
           dur = 4.3/24d0
           bins = 400
        end
        210750726: begin
           per = 4.612
           t0 = -1
           range = [0.8,1.1]
           bins = 250
        end
        212009427: begin
           fullrange = 1
           per = 1.55678        ;0.7764            ;
           bins = 40
           t0 = 2308.03069+0.7764
           smooth = 12
        end
        211969807: begin
           per = 1.97434
           range = [1.2,1.5]
           t0=2307.37773
           smooth = 20
           bins = 300
        end
        211969807: begin
           fullrange=1
           per = 34.691367
           range = [30,31]
           t0=2337.1311
           smooth = 40
           bins = 300
        end 
        211938988: begin
           per = 1.088631
           ;;range = [0.8,1.09]
           t0 = 2306.756393
           smooth = 20
           bins=400
        end
        211916756: begin
           per = 10.13582
           range = [6,8]
           t0 = 2307.74318
           smooth = 15
           dur = 0.2
           fullrange=1
        end
        210744818: begin
           per = 12.693054
           t0 = -1              ;2299.63767112
           range = [1,12]
           smooth = 10
        end
        211799258: begin
           t0 = 2320.15668
           per = 19.53307
           range = [15,16]
           smooth = 30
           bins = 400
           fullrange=1
        end
        211946007: begin
           t0 = 2308.14121
           per = 1.98289
        end
        210696763: begin
           t0 = 2230.15668
           per = 7.9844
           range = [7d,7.5d]
           smooth = 50
           bins = 250
           per = 3.65101
           range = [1.1,1.6]
           bins = 300
        end
        211901114: begin
           per = 1.6490351
           range = [0.75,0.9]
           t0 = -1
           smooth = 15
           bins=400
        end
        211822797: begin
           per = 21.170939
           range = [3,4.5]
           smooth = 25
           bins = 300
           t0 = -1
        end
        211102816: begin
           per = 10.55001
           range = [3.5,4.5]
           smooth = 40
           t0 = 2234.77992
           bins = 450
        end
        210931967: begin
           per = 0.34658302
           range = [0.0,0.1]
           smooth = 20
           bins = 400
        end
        210838726: begin
           per = 1.096
           range = [0.35,0.55]
           smooth = 25
           bins = 400
        end
        210659688: begin
           per = 2.35774
           range = [1.5,1.7]
           smooth = 50
           bins = 300
        end
        210651828: begin
           per = 12.76788
           range = [3.0,5.0]
           smooth = 40
           bins = 350
        end
        210527821: begin
           per = 9.58702
           range = [4,6]
           smooth = 30
           bins = 350
        end
        205046529: begin
           per = 1.83549136621d0
           t0 = 2061.47740498d0
           smooth = 10
           fullrange=0
           dur = 0.19520789
        end
        211970147: begin
           per = 9.9188643
           smooth = 10
           ;;range = [5.2,5.8]
           t0 = 2316.553
           bins = 600
           dur = 0.1
        end
        210463428: begin
           fullrange=1
           per = 4.72104
           smooth = 15
           t0 = 2234.12841
           range = [1.0,1.2]
           bins = 300
        end
        210532473: begin
           per = 12.63194
           smooth = 20
           t0 = 2242.65763
           range = [6.5,7.2]
        end
        210673842: begin
           per = 17.03056
           t0 = 2234.70520
           range = [3.4,4.2]
           smooth = 15
        end
        210822691: begin
           per = 8.0957196/2.0
           smooth = 50
           t0 = 2231.0476
           fullrange = 1
           debug = 0
        end
        211143504: begin
           per = 0.44613        ;*4.0
           t0 = 2231.32900
           range = dblarr(2)
           range[0] = (t0 mod per) - 0.2
           range[1] = (t0 mod per) + 0.2
        end
        210363145: begin
           per = 8.2
           t0 = 2237.80344
           fullrange = 0
           smooth = 15
           fullrange = 0
        end
        202900527: begin
           per = 13.008503
           t0 = 2072.757762
           range = dblarr(2)
           smooth = 20
        end
        210490365: begin
           per = 3.484552
           t0 = 2229.57935      ;+2454833d0
           smooth = 10
           fullrange=0
        end
                                ;211913977: begin
                                ;   per = 14.6891;4.8313492
                                ;   t0 = 2319.66411133
                                ;   ;;range = [13.2,13.9]
                                ;   smooth = 10
                                ;   fullrange=0
                                ;   bins = 400
                                ;   dur = 0.12
                                ;end
        211913977: begin
           per = 28.011
           t0 = 2333.1306 ;;2319.66411133
           smooth = 10
           fullrange=0
           bins = 600
           dur = 0.22
           fullrange=2
        end
        211923502: begin ;; this is a false positive from an EB
           per = 9.497737
           t0 = -1
           range = [3,9]
           fullrange = 0
           smooth = 20
        end
        211990866: begin
           per = 1.673615       ;1.673763;
           t0 = 2307.72411133
           smooth = 10
           dur = 0.15
           bins = 400
        end
        211093684: begin
           per = 7.05157
           t0 = 2231.71183
           smooth = 30
        end
        211972086: begin
           per = 3.00752
           t0 = 2306.88141
           smooth = 30
           dur = 0.2
        end
        212002525: begin
           per = 5.80671
           t0 = 2308.80914
           smooth = 50
        end
        211635844: begin
           per = 10.7619
           t0 = 2302.3128
           smooth = 20
           dur = 0.2
        end
        210974364: begin
           per = 32.747
           t0 = 7068.748d0+2450000d0-2454833d0
           smooth = 20
        end
        211075914: begin
           per = 42.8
           t0 = 7099.2d0+2450000d0-2454833d0
           smooth = 5
           fullrange=1
        end
        211705654: begin
           per = 2.5331
           t0 = 2305.444062
           smooth = 10
           dur = 0.1
        end
        210894022: begin
           per = 5.349
           t0 = -1
           range = [4,5]
           smooth = 20
           fullrange=0
        end
        210776021: begin
           t0 = -1
           per = 17.840
           range = [12,14]
           smooth = 25
        end
        204262368: begin
           per = 22.69
           range = [15,18]
        end
        else: begin
           print,'No parameters, goint to default'
           smooth = 20
           t0 = -1
           per = 30.0
           range = [0,29.5]
        end
     endcase
  endif
  if n_elements(t0) eq 0 then t0 = -1
  if t0[0] ne -1 and n_elements(range) lt 2 then begin
     range = dblarr(2,n_elements(per))
     for jj = 0,n_elements(t0)-1 do begin
        range[0,jj] = (t0[jj] mod per[jj]) - dur[jj]*1.5d0
        range[1,jj] = (t0[jj] mod per[jj]) + dur[jj]*1.5d0
     endfor
  endif
  ;;print,range
  ;;print,per
  ;;if type eq 'A' then smooth = 100
  if n_elements(range) lt 2 then range = [0,per*0.95]
  if n_elementS(smooth) eq 0 then smooth = 25

  per = 1d0*per
  if n_elements(type) eq 0 then type = 'V'
  set_plot,'x'
  device,retain=2
  charsize=2
  thick=4
  if n_elements(print) eq 0 then print = 0
  lc_corr = []
  ;;readcol,'hlsp_k2sff_k2_lightcurve_'+strtrim(string(epic,format="(I10)"),2)+'-c04_kepler_v1_llc-default-aper.txt',time,lc
  if type eq 'V' then begin
     if epic eq 246865365 then begin
        x = mrdfits('~/Desktop/H2/246865365-c13_llc.fits',1,header)
        lc = x.PDCSAP_FLUX
        time = x.time
     endif else begin
        if epic eq 247589423 then readcol,'~/Desktop/H2/ep247589423simultfit.csv',time,lc,lc_corr,format='d,d,',delimiter=',',/silent else begin
           spawn,'ls ~/Desktop/H2/*/*'+strtrim(string(epic,format="(I10)"),2)+'*.csv',file
           if file[0] eq '' then begin
              spawn,'ls ~/Desktop/H2/*'+strtrim(string(epic,format="(I10)"),2)+'*.txt',file
           endif 
           file = strtrim(file[0],2)
           if file ne '' then begin
              readcol,file,time,lc,/silent,format='d,d',delimiter=','
           endif
        endelse
     endelse
     ;;if epic eq 205117205 then begin
     ;;   readcol,'~/Dropbox/Python_transit/ep205117205.csv',time,lc,lc_corr,format='d,d,d',/silent
     ;;   lc_corr = lc/lc_corr
     ;;endif else begin
     ;;spawn,'ls /Volumes/Loki/K2Lightcurves/*/hlsp_k2sff_k2_lightcurve_'+strtrim(string(epic,format="(I10)"),2)+'*.fits',list
     ;;list = list[0]
     ;; if file_test(list) then begin
     ;;    x = mrdfits(list,1,header,/silent)
     ;;    time = x.t
     ;;    lc = x.fcor
     ;; endif else begin
     ;;    fname = '/Volumes/Vali/K2phot/ep'+strtrim(string(epic,format="(I10)"),2)+'.csv'
     ;;    if file_test(fname) and 1 eq 2 then begin
     ;;       readcol,fname,time,lc,format='d,d',/silent
     ;;    endif else begin
     ;;       fname = '/Volumes/Vali/K2phot/'+strtrim(string(epic,format="(I10)"),2)+'_Corrected.txt'
     ;;       if file_test(fname) and 1 eq 2 then begin
     ;;          readcol,fname,time,lc,format='d,d',/silent
     ;;       endif else begin
     ;;          fname = '/Volumes/Vali/K2phot/ep'+strtrim(string(epic,format="(I10)"),2)+'.csv'
     ;;          if file_test(fname) then begin
     ;;             readcol,fname,time,lc,format='d,d',/silent
     ;;             lc/=median(lc)
     ;;          endif else begin
     ;;             spawn,'ls /Volumes/Vali/K2phot/hlsp_k2sff_k2_lightcurve_'+strtrim(string(epic,format="(I10)"),2)+'-c*_kepler_v1_llc-default-aper.txt',filename
     ;;             if file_test(filename) then begin
     ;;                readcol,filename,time,lc,format='d,d',/silent
     ;;             endif else stop;return
     ;;          endelse
     ;;       endelse
     ;;    endelse
     ;; endelse
     ;;endelse
  endif
  if type eq 'C' then begin
     spawn,'ls /Volumes/Work/CoveyLC/*'+strtrim(string(epic,format="(I10)"),2)+'*.lc',filename
     if strtrim(filename[0],2) ne '' then begin
        readcol,filename[0],time,lc,err,fl,fl_err,pos,format='d,d,d,d,d,d',/silent
        lc/=median(lc)
        ;;lc = fl
     endif else begin
        return
     endelse
  endif
  if type eq 'E' then begin
     spawn,'ls /Volumes/WORK/EVEREST/hlsp_everest_k2_llc_'+strtrim(string(epic,format="(I10)"),2)+'*.fits',filename
     filename = filename[0]
     if file_test(filename) then begin
        x = mrdfits(filename[0],1,header,/silent)
        lc = x.FLUX
        time = x.time
        lc/=median(lc)
     endif else return
  endif
  if type eq 'K' then begin
     spawn,'ls /Volumes/Vali/K2phot/ktwo'+strtrim(string(epic,format="(I10)"),2)+'*llc.fits',filename
     if file_test(filename) then begin
        x = mrdfits(filename[0],1,header,/silent)
        lc = x.PDCSAP_FLUX
        time = x.time
        lc/=median(lc)
     endif else return
  endif
  if type eq 'A' then begin
     spawn,'ls /Volumes/Vali/K2SC/hlsp_k2sc_k2_llc_'+strtrim(string(epic,format="(I10)"),2)+'*_kepler_v1_lc.fits',filename
     if file_test(filename) eq 1 then begin
        tmp = mrdfits(filename,1,header)
        time = tmp.time
        lc = tmp.flux
        err = tmp.error
        err/=median(lc)
        lc/=median(lc)
        lc_corr = lc
        lc*=tmp.trend_t
     endif else return
  endif
  l = where(time gt 0 and lc gt 0 and finite(time) eq 1 and finite(lc) eq 1 and time lt datemax and time gt datemin)
  time = time[l]
  lc = lc[l]
  if n_elements(lc_corr) gt 0 then lc_corr = lc_corr[l]
  lc*=1d0
  time*=1d0
  
                                ;if epic eq 205117205 and type eq 'V' then begin ;; fix the discontinuity
                                ;   after = where(time ge 2104.1714)
                                ;   lc[after]/=(lc[1772]/lc[1771])
                                ;endif

  ;;if type eq 'V' then lc += 0.0004*randomn(seed,n_elements(lc))     
  
  if epic eq 210363145 and type eq 'V' then begin
     readcol,'ep210363145simultfit.csv',time,lc,lc_corr,format='d,d,d',/silent
  endif
  if epic eq 210490365 and type eq 'V' then begin
     readcol,'~/Dropbox/epic_mdwarf/epic210490365simultfit.csv',time,lc,lc_corr,format='d,d,d',/silent
     tmp = robust_poly_fit(time,lc,8,yfit)
     yfit = am_runmed(lc,100)
     lc/=yfit
  endif


  checkper = 0
  checkper=1
  if checkper eq 1 then begin
     ;;tmp = lc
     ;;l = wherE(lc gt 1.0 and lc lt 2)
     l = where(lc gt 0)
     doperiodogram,epic,time[l],lc[l]
     ;;lc = tmp
  endif

  !x.margin=[10,2]
  if debug eq 0 then begin
     !p.multi=[0,1,3]
     set_plot,'PS'
     adder = ''
     if datemin gt 0 then adder += '_'+strtrim(datemin,2)
     if datemax lt 1d5 then adder+=+'_'+strtrim(datemax,2)
     device,filename='Vetfiles/'+strtrim(string(epic,format="(I10)"),2)+'_vet1_'+type+adder+'.eps',/encapsul,/color,xsize=20,ysize=10
  endif else begin
     !p.multi=[0,2,3]
     set_plot,'x'
  endelse
  
  if type eq 'A' then lc = lc_corr
  ;;if epic eq 210363145 or (epic eq 210490365 and type eq 'V') or (epic eq 205117205 and type eq 'V') then begin
  ;;   lc = lc_corr
  ;;endif
  if epic eq 205117205 and type eq 'C' then begin ;; cut out outliers
     ll = where(lc gt 0.99)
     lc = lc[ll]
     time = time[ll]
  endif
  if epic eq 205046529 then begin
     ll = where(lc lt 1.005 and lc gt 0.97)
     lc = lc[ll]
     time = time[ll]
  endif

  if epic ne 247589423 then lc_corr = []
  makelightcurveplots,time,lc,epic=epic,fullrange=fullrange,t0=t0,per=per,range=range,fitpoly=fitpoly,debug=debug,$
                      dur=dur,$
                      smooth=smooth,$
                      print=print,$
                      type=type,$
                      datemin=datemin,$
                      datemax=datemax,$
                      charsize=charsize,$
                      thick=thick,$
                      charthick=charthick,$
                      bins=bins,$
                      lc_corr = lc_corr,$
                      everyother=everyother
  
  ;; if debug eq 0 then device,/close

  ;; if print eq 1 then begin
  ;;    close,15
  ;;    openw,15,'~/Desktop/TAP_K2/'+strtrim(string(epic,format="(I10)"),2)+'_full_'+type+'.txt'
  ;;    printf,15,'#  TAP combination of setup parameters and light curves.  Adjust parameter setup matrix to'
  ;;    printf,15,'#  change the setup before loading this file as a "Transit File" in a new instance of TAP'
  ;;    printf,15,'#'
  ;;    printf,15,'#  transit  long_int    int_length_min  cadence_multiplier'
  ;;    printf,15,'#         1         1         29.4244              10'
  ;;    printf,15,'#'
  ;;    printf,15,'#  transit  set               param                        value lock                      low_lim                       hi_lim  limited   prior=1_penalty=2            val            sig'
  ;;    printf,15,'#         1     1              Period                 '+string(per,format="(D10.6)")+'    0                 3.0                10.0000000000        0                   0     '+string(per,format="(D10.4)")+'       0.000057'
  ;;    printf,15,'#         1     1                   T                 0.026   0                0.0                1.        1                   0        0.00000        0.00000'
  ;;    printf,15,'#         1     1                   b                 0.002355859    0                 -1.2               1.2        1                   0        0        0'
  ;;    printf,15,'#         1     1               Rp/R*                  0.1    0                 0.01                 0.3        1                   0        0.00000        0.00000'
  ;;    if type eq 'G' then $
  ;;       printf,15,'#         1     1         Mid_Transit               0.78     0               0.6              0.9        1                   0      0        0'
  ;;    if type eq 'V' then $
  ;;       printf,15,'#         1     1         Mid_Transit               2229.5783691406      0               2229.5              2229.65        1                   0      0        0'
  ;;    printf,15,'#         1    1           Linear_LD                  0.45    0                 0.0000000000                 1.0000000000        1                   2         0.45         0.1'
  ;;    printf,15,'#         1    1             Quadratic_LD                  0.35    0                 0.0000000000                 1.0000000000        1                   2         0.35         0.1'
  ;;    printf,15,'#         1     1        Eccentricity                 0.00000000001    0                 0.0000001                 1.0000000        1                   0        0.0        0.0'
  ;;    printf,15,'#         1     1               Omega                  0.0000    0                 0.0000000000                 6.2831853072        1                   0        0.00000        0.00000'
  ;;    printf,15,'#         1     1             OOT_t^0                 1.0000000000    1                 1.00000                 1.0001000000        0                   0        0.00000        0.00000'
  ;;    printf,15,'#         1     1             OOT_t^1                 0.0000000000    1                 0.0000                 0.0000000100        0                   0        0.00000        0.00000'
  ;;    printf,15,'#         1     1             OOT_t^2                 0.0000000000    1                0.00000                 0.0001000000        0                   0        0.00000        0.00000'
  ;;    printf,15,'#         1     1           Sigma_Red                 0.0000000000    0                 0.0000000000                 0.0015000000        0                   0        0.00000        0.00100'
  ;;    printf,15,'#         1     1         Sigma_White                 0.0010000000    0                 0.0000000000                 0.0100000000        0                   0        0.00000        0.00000'
  ;;    printf,15,'# transit      1'

  ;;    index = 0
  ;;    if type eq 'G' then l = where(time mod period gt range[0] and time mod period lt range[1])
  ;;    for i = 0,n_elements(l)-1 do begin
  ;;       printf,15,time[l[i]],lc[l[i]]
  ;;    endfor
  ;;    index++
  ;;    close,15
  ;; endif
  close,/all
  ;;device,/close
  ;;stop
  
END


PRO evenodd,per,time,lc,modn,bintime=bintime,binlc=binlc,err=binlc_err,bins=bins,range=range,epic=epic,type=type,yrange=yrange,xline=xline,nox=nox ;,fullrange=fullrange

  if n_elements(print) eq 0 then print = 0
  if n_elements(bins) eq 0 then bins = 100
  tt = min(time)
  foldtime = time
  counter = 0
  while max(foldtime) gt per do begin
     l = where(foldtime gt per)
     if l[0] ne -1 then foldtime[l]-=per
     ll = where(foldtime gt per and foldtime le per*2.0)
     if ll[0] ne -1 then begin
        lc[ll]/=median(lc[ll])
        if counter mod 2 eq modn then lc[ll] = 0
     endif
     counter++
  endwhile
  qq = where(lc gt 0)
  foldtime = foldtime[qq]
  lc = lc[qq]
  tmp = lc[sort(lc)]
  if n_elements(yrange) lt 2 then yrange = [tmp[n_elements(tmp)*0.000],tmp[n_elements(tmp)*0.99]]
  qq = sort(foldtime)
  foldtime = foldtime[qq]
  lc = lc[qq]
  if n_elements(nox) eq 1 then begin
     plot,foldtime,lc,psym=8,symsize=0.5,/xstyle,/ystyle,charsize=1.25,charthick=3,xthick=3,ythick=3,yrange=yrange,xrange=range,xtickn=strarr(10)+' '
  endif else plot,foldtime,lc,psym=8,symsize=0.5,/xstyle,/ystyle,xtitle='Phased time (days)',charsize=1.25,charthick=3,xthick=3,ythick=3,yrange=yrange,xrange=range
  bintime = dblarr(bins)
  binlc = bintime
  binlc_err = bintime
  space = max(foldtime)/bins
  for i = 0d0,bins-1d0 do begin
     qq = where(foldtime lt space*(i+1.) and foldtime ge space*(i))
     if qq[0] ne -1 then begin
        binlc[i] = median(lc[qq])
        bintime[i] = median(foldtime[qq])
        tmp = robust_sigma(lc[qq])/sqrt(n_elements(qq)-1)
        if abs(tmp) gt 0.1 then tmp = stdev(lc[qq])/sqrt(n_elements(qq)-1)
        binlc_err[i] = tmp
     endif else begin
        binlc[i] = 1d
        binlc_err[i] = 1d
     endelse
  endfor
  oploterror,bintime,binlc,dblarr(n_elements(binlc)),binlc_err,psym=8,color=cgcolor('red'),errcolor=cgcolor('red')
  if n_elements(xline) eq 0 then begin
     tmp = binlc[sort(binlc)]
     xline = median(tmp[0:3])
  endif
  xline = 1d0-0.00260
  oplot,range,[xline,xline],linestyle=2,color=cgcolor('blue'),thick=4
  ;;qq = where(foldtime gt 4.2 and foldtime lt 4.28)
  ;;print,'even/odd',1d0-median(lc[qq]),1d0-mean(lc[qq]),robust_sigma(lc[qq]),robust_sigma(lc[qq])/sqrt(n_elements(qq)-1)
  
END

PRO plotter,per,time,lc,bintime=bintime,binlc=binlc,err=binlc_err,bins=bins,print=print,range=range,epic=epic,type=type,fullrange=fullrange,yrange=yrange,charsize=charsize,charthick=charthick

  if n_elements(print) eq 0 then print = 0
  if n_elements(bins) eq 0 then bins = 100
  tt = min(time)
  foldtime = time
  counter = 0
  ;; while max(foldtime) gt per do begin
  ;;    l = where(foldtime gt per)
  ;;    if l[0] ne -1 then begin
  ;;       foldtime[l]-=per
  ;;       ll = where(foldtime gt per and foldtime le per*2.0)
  ;;       ;;if counter mod 2 eq 1 then begin
  ;;       if ll[0] ne -1 then lc[ll]/=median(lc[ll])
  ;;       ;;endif else begin
  ;;       ;;   if ll[0] ne -1 then lc[ll]=0
  ;;       ;;endelse
  ;;    endif
  ;;    counter++
  ;; endwhile

  !p.multi=[n_elements(per),n_elements(per),3]
  qq = where(lc ne 0)
  time_base = time[qq]
  lc_base = lc[qq]

  for jj = 0,n_elements(per)-1 do begin
     lc = lc_base
     time = time_base
     foldtime = time mod per[jj]

     tmp = lc[where(foldtime gt range[0,jj] and foldtime lt range[1,jj])]
     tmp = tmp[sort(tmp)]
     if n_elements(yrange) le n_elements(per)*2 then yrange = [tmp[n_elements(tmp)*0.005],tmp[n_elements(tmp)*0.995]]
     if fullrange eq 1 then yrange = [min(lc),max(lc)]
     if fullrange eq 2 then yrange = [tmp[n_elements(tmp)*0.001],tmp[n_elements(tmp)*0.995]]
     if type eq 'A' then yrange = [tmp[n_elements(tmp)*0.005],tmp[n_elements(tmp)*0.97]]

     ;;if jj eq 1 then range[*,jj] = [0,12]
     qq = sort(foldtime)
     foldtime = foldtime[qq]
     lc = lc[qq]
     if jj eq fix(n_elements(per)/2) then xtitle='Phased time (days)' else xtitle=''
     plot,foldtime,lc,psym=8,symsize=0.5,/xstyle,/ystyle,xtitle=xtitle,yrange=yrange,xrange=range[*,jj],charsize=charsize,charthick=charthick,xthick=charthick,ythick=charthick ;,yrange=[0.986,1.005];,xrange=[0.5,1.0];,ytitle='Normalized Flux'
     
     ;;space = 20.0/60./24.          ;
     ;;bins = max(foldtime)/space
     bintime = dblarr(bins)
     binlc = bintime
     binlc_err = bintime
     space = max(foldtime)/bins
     for i = 0d0,bins-1d0 do begin
        qq = where(foldtime lt space*(i+1.) and foldtime ge space*(i))
        if qq[0] ne -1 then begin
           binlc[i] = median(lc[qq])
           bintime[i] = median(foldtime[qq])
           binlc_err[i] = robust_sigma(lc[qq]) ;/sqrt(n_elements(qq)-1)
           ;;binlc_err[i] = stdev(lc[qq]);/sqrt(n_elements(qq)-1)
        endif else begin
           print,'Too many bins?'
           binlc[i] = 1d
           bintime[i] = 1d
           binlc_err[i] = 1d    ;
        endelse
     endfor
     ;;print,'Space: ',space*60*24.
     print,'depth:' ,(1d0-min(binlc))*1000d0
     ;;ll = where(bintime gt 2.4 and bintime lt 2.6 and binlc gt 0.999)
     ;;binlc[ll]*=0.9985
     ;;oploterror,bintime,binlc,dblarr(n_elements(binlc)),binlc_err,psym=8,color=cgcolor('red'),errcolor=cgcolor('red'),errthick=4
     ;;oplot,bintime,binlc,psym=8,color=cgcolor('red')
  endfor

END


PRO TAP_fit

  thick = 4
  charsize = 1.5
  readcol,'~/Desktop/TAP_K2/TAPmcmc_20151014_0915/ascii_phased_data.ascii',time,lc,dummy,orig_time,model,residual,/silent
  ;;readcol,'~/Desktop/TAP_K2/TAPmcmc_20151014_0750//ascii_phased_model.ascii',junk,time2,model_flux,/silent
                                ;newtime = generatearray(min(time),max(time),n_elements(time)*10)
                                ;newmodel = interpol(model,time,newtime)
                                ;newmodel = TS_SMOOTH(newmodel,10)
                                ;model = smooth(model,12)
                                ;model = ts_smooth(model,12)

  ;;fancymodel = 0.0
  ;;for i = 4,n_elements(newmodel)-6,10 do begin
  ;;   fancymodel = [fancymodel,wmean(newmodel[i-4:i+5],newmodel[i-4:i+5])]
  ;;endfor
  ;;shrink,facnymodel
  ;;newmodel = ((dblarr(10.)+10.^(-1d0))#reform(newmodel,10,n_elements(time)))[0,*]

  ;; set_plot,'PS'
  ;; device,filename='TAP_fit.eps',/encapsul,/color,xsize=7,ysize=5
  ;; !p.multi=[0,1,2]
  ;; !y.margin=[-2,1]
  ;; plot,time*24,lc,psym=8,xrange=[-2,2],/xstyle,/ystyle,xthick=thick,ythick=thick,charthick=thick,charsize=charsize,ytitle='Normalized Flux',xtickname=strarr(10)+' ',yrange=[0.986,1.002]
  ;; oplot,time*24,model,thick=thick,color=cgcolor('red') ;,psym=8
  ;; ;;oplot,time2*24,model_flux,thick=thick,color=cgcolor('blue')
  ;;                               ;oplot,newtime*24,newmodel,thick=thick,color=cgcolor('blue'),psym=2
  ;; ;;oplot,time*24,fancymodel,thick=thick,color=cgcolor('green'),psym=1
  ;; !y.margin=[4,2]
  ;; plot,time*24,1000*((lc-model)),/xstyle,/ystyle,xtitle='Time Since Mid-Transit (hours)',psym=8,xthick=thick,ythick=thick,charthick=thick,charsize=charsize,ytitle='Residual (10!U-4!N)',xrange=[-2,2],yrange=[-2.5,2.5]
  ;; oplot,time*24,time*0,thick=thick,color=cgcolor('red')

  ;; device,/close

  l = where(time*24. gt -0.9 and time*24. lt 0.9,comp=g)
  ;;print,total(abs((lc[l]-model[l])/lc[l]))

  ;;print,robust_sigma((lc[l]-model[l]))
  g = where(time*24 gt 0.9 or time*24 lt -0.9)
  ;;print,robust_sigma((lc[g]-model[g]))  

END


function am_runmed,arr,length

  tmp = findgen(n_elements(arr))
  smoothed = arr*0d0
  for i = 0d0,1d0*n_elements(arr)-1d0 do begin
     l = where(abs(tmp-i) lt length)
     smoothed[i] = robust_mean(arr[l],5)
  endfor
  
  return,smoothed
end


PRO doperiodogram,epic,time,lc

  ;; first filter out low-power
  res = robust_poly_fit(time,lc,3,yfit)
  newlc = lc/yfit
  ;;yfit = am_runmed(lc,200)
  ;;lc/=yfit
  !p.multi=[0,1,1]
  !x.margin=[8,1]
  !y.margin=[5,1]
  !p.font=0
  charsize=1.5
  thick=7
  set_plot,'PS'
  device,filename='rotfiles/'+strtrim(string(epic,format="(I10)"),2)+'_rotpower.eps',xsize=25,ysize=20
  !p.multi=[0,1,1]
  !y.margin=[5,1]
  !x.margin=[8,3]
  power=periodogram(time,newlc,per=[0.1,30.],/quiet)
  plot,power[0,*],power[1,*],/xstyle,/ystyle,/xlog,xtitle='Period (days)',ytitle='Power',xthick=thick,ythick=thick,thick=thick,charsize=charsize
  period = power[0,*]
  pow = power[1,*]
  tmp1 = period[where(period lt 50 and period gt 0.1)]
  tmp2 = pow[where(period lt 50 and period gt 0.1)]
  p = tmp1[where(tmp2 eq max(tmp2))] ;;50.334360
  p = p[0]

  tmp1 = period
  tmp2 = pow

  res = gaussfit(tmp1,tmp2,a,nterms=4)
  oplot,tmp1,res,color=cgcolor('red'),thick=5
  print,'Best period: '+string(tmp1[where(tmp2 eq max(tmp2))],format="(D10.7)") + ' +/- ',a[2]
  bestper = tmp1[where(tmp2 eq max(tmp2))]
  bestper = a[1] ;;bestper[0]
  device,/close

  !p.multi=[0,1,1]
  !x.margin=[9,2]
  !y.margin = [4,1]
  device,filename='rotfiles/'+strtrim(string(epic,format="(I10)"),2)+'_phase.eps',xsize=25,ysize=20 ;,xsize=25,ysize=10
  cgloadct,1,clip=[25,220]
  col = (max(time)-min(time))/bestper
  colorinc = round(255/col)
  inc = 0

  tmp = newlc[sort(newlc)]
  yrange = [tmp[n_elements(tmp)*0.02],tmp[n_elements(tmp)*0.98-1]]
  plot,[0],[0],xrange=[0,1],yrange=yrange,xtitle='Rotational Phase',ytitle='Relative Flux',charsize=charsize,charthick=thick,xthick=thick,ythick=thick
  sharpcorners,thick=thick
  for time0 = min(time),max(time),bestper do begin
     loc = where(time ge time0 and time lt time0+bestper)
     if n_elements(loc) gt 1 then oplot,(time[loc]-time0)/bestper,newlc[loc],psym=8,color=colorinc*inc
     inc++
  endfor
  
  ;;plot,time mod bestper,newlc,psym=8,/xstyle,/ystyle
  device,/close
  
END


PRO makelightcurveplots,time,lc,epic=epic,fullrange=fullrange,t0=t0,per=per,range=range,fitpoly=fitpoly,debug=debug,$
                        dur=dur,$
                        smooth=smooth,$
                        print=print,$
                        type=type,$
                        datemin=datemin,$
                        datemax=datemax,$
                        charsize=charsize,$
                        thick=thick,$
                        charthick=charthick,$
                        bins = bins,$
                        lc_corr=lc_corr,$
                        everyother=everyother

  lc/=median(lc)
  baset = 2454833d0             ;2450000d0

  tmp = lc[sort(lc)]
  yrange = [tmp[n_elements(tmp)*0.01],tmp[n_elements(tmp)*0.99]]
  if fullrange eq 1 then yrange = [min(lc),tmp[n_elements(tmp)*0.999]]
  if fullrange eq 2 then yrange = [tmp[n_elements(tmp)*0.001],tmp[n_elements(tmp)*0.995]]
  t = time+2454833d0-baset
  counter = findgen(n_elements(t))
  if everyother eq 1 then qq = where(counter mod 2 eq 1) else qq = where(counter gt -100)
  !y.margin=[1.5,1]
  plot,t[qq],lc[qq]/median(lc[qq]),/xstyle,/ystyle,psym=8,symsize=0.5,yrange=yrange,charsize=charsize,charthick=charthick,xthick=charthick,ythick=charthick,xtickn=strarr(20)+' '
  count = 0
  col = ['red','blue','teal']
  for jj = 0,n_elements(t0)-1 do begin  
     t0_c = t0[jj]+2454833d0-baset
     if t0_c gt 0 then begin
        count = 0
        place = t0_c-(per[jj]*100.) + (per[jj]*count)
        while place lt max(t) do begin
           place = t0_c-(per[jj]*100.) + (per[jj]*count)
           oplot,[place,place],[min(lc),max(lc)],thick=thick,linestyle=2,color=cgcolor(col[jj])
           count++
        endwhile
     endif
  endfor
  tmp = am_runmed(lc,smooth/2.)
  for jj = 0,n_elements(t0)-1 do begin
     transitplace = where((time mod per[jj]) gt ((t0[jj] mod per[jj]) - dur[jj]/1.75) and (time mod per[jj]) lt ((t0[jj] mod per[jj]) + dur[jj]/1.75),comp=notransit)
     ;;transitplace = where((time mod per) gt ((t0 mod per) - dur/1.25) and (time mod per) lt ((t0 mod per) + dur/1.25),comp=notransit)
     ;;transitplace = where((time mod per) gt ((t0 mod per) - dur/1.1) and (time mod per) lt ((t0 mod per) + dur/1.1),comp=notransit)
     num = findgen(n_elements(lc))
     match2,transitplace,num,suba,subb
     for i = 0,n_Elements(transitplace)-1 do begin
        replace = where(abs(transitplace[i]-num) lt 25 and subb eq -1)
        if replace[0] ne -1 then begin
           if fitpoly eq 1 then begin
              fit = linfit(time[replace],tmp[replace])
              newval = fit[0]+fit[1]*time[transitplace[i]]
           endif else begin
              newval = interpol(tmp[replace],time[replace],time[transitplace[i]])
           endelse
           tmp[transitplace[i]]=newval
        endif else begin
           print,transitplace[i],num
        endelse
     endfor
  endfor
  oplot,time,tmp,color=cgcolor('blue'),thick=3
  !y.margin=[4,-1.5]
  lc/=tmp

  ;; if epic eq 210894022 and type eq 'K' then begin
  ;;    fname = '/Volumes/Vali/K2phot/ep'+strtrim(string(epic,format="(I10)"),2)+'.csv'
  ;;    if file_test(fname) then begin
  ;;       readcol,fname,time,lc,format='d,d'
  ;;       lc/=median(lc)
  ;;    endif else begin
  ;;       spawn,'ls /Volumes/Vali/K2phot/hlsp_k2sff_k2_lightcurve_'+strtrim(string(epic,format="(I10)"),2)+'-c*_kepler_v1_llc-default-aper.txt',filename
  ;;       if file_test(filename) then begin
  ;;          readcol,filename,time,lc,format='d,d'
  ;;          lc/=median(lc)
  ;;       endif else return
  ;;    endelse
  ;;    tmp = am_runmed(lc,smooth/2.)
  ;;    transitplace = where((time mod per) gt ((t0 mod per) - dur/1.75) and (time mod per) lt ((t0 mod per) + dur/1.75))
  ;;    num = findgen(n_elements(lc))
  ;;    match2,transitplace,num,suba,subb
  ;;    for i = 0,n_Elements(transitplace)-1 do begin
  ;;       replace = where(abs(transitplace[i]-num) lt 15 and subb eq -1)
  ;;       fit = linfit(time[replace],tmp[replace])
  ;;       newval = fit[0]+fit[1]*time[transitplace[i]]
  ;;       tmp[transitplace[i]]=newval
  ;;    endfor
  ;; endif
  

  ;; calculate the depth
  for jj = 0,n_elements(t0)-1 do begin
     if t0[jj] ne -1 then begin
        transitplace_narrow = where((time mod per[jj]) gt ((t0[jj] mod per[jj]) - dur[jj]/2.5) and (time mod per[jj]) lt ((t0[jj] mod per[jj]) + dur[jj]/2.5),comp=notransit)
        d = median(lc[transitplace_narrow])
        print,'Depth_med: ',100d0*(1d0-d),100d0*robust_sigma(lc[transitplace_narrow])
        print,'Depth_mean: ',100d0*(1d0-robust_mean(lc[transitplace_narrow],4)),100d0*robust_sigma(lc[transitplace_narrow])
        transitplace_wide = where((time mod per[jj]) gt ((t0[jj] mod per[jj]) - dur[jj]) and (time mod per[jj]) lt ((t0[jj] mod per[jj]) + dur[jj]),comp=notransit)
        print,'Depth_wide: ',1d0-median(lc[transitplace_wide]),robust_sigma(lc[transitplace_wide])
     endif
  endfor
  if epic eq 205117205 then begin
     ll = where(lc lt 1.003)
     lc = lc[ll]
     time = time[ll]
  endif
  if epic eq 210923016 then begin
     ll = where(lc lt 1.01 and lc gt 0.98)
     lc = lc[ll]
     time = time[ll]
  endif

  if n_elements(lc_corr) gt 0 then lc = lc_corr
  
  ;;cut = where(abs(lc - median(lc)) lt robust_sigma(lc)*10d0)
  ;;time = time[cut]
  ;;lc = lc[cut]
  tmp = lc[sort(lc)]
  yrange = [tmp[n_elements(tmp)*0.005],tmp[n_elements(tmp)*0.995]]
  if fullrange eq 1 then yrange = [tmp[n_elements(tmp)*0.00],tmp[n_elements(tmp)*0.999]]
  if fullrange eq 2 then yrange = [tmp[n_elements(tmp)*0.001],tmp[n_elements(tmp)*0.995]]
  if epic eq 210776021 then yrange=[tmp[n_elements(tmp)*0.01],tmp[n_elements(tmp)*0.9]]
  if epic eq 203848625 then yrange=[tmp[n_elements(tmp)*0.01],tmp[n_elements(tmp)*0.92]]
  t = time+2454833d0-baset
  plot,t,lc,/xstyle,/ystyle,psym=8,symsize=0.5,yrange=yrange,xtitle='Time (BJD-'+strtrim(string(baset,format="(I10)"),2)+')',charsize=charsize,charthick=charthick,xthick=charthick,ythick=charthick,ytitle='Relative Brightness'
  for jj = 0,n_elements(t0)-1 do begin
     t0_c = t0[jj]+2454833d0-baset
     if t0_c gt 0 then begin
        count = 0
        place = t0_c-(per[jj]*100.) + (per[jj]*count)
        while place lt max(t) do begin
           place = t0_c-(per[jj]*100.) + (per[jj]*count)
           oplot,[place,place],[0,10],thick=thick,linestyle=2,color=cgcolor(col[jj])
           ;;if place gt min(t) and place lt max(t) then print,place,per[jj],t0[jj]
           count++
        endwhile
     endif
  endfor
  !y.margin=[4,0]
  
  l = where(time mod per gt range[0] and time mod per lt range[1],comp=g)
  mad=MEANABSDEV(lc[g])
  stdev = robust_sigma(lc[g])
  ;;print,stdev*1000d0
  if print eq 1 then begin
     forprint,string(time,format="(D20.10)")+string(9b)+string(lc,format="(D20.10)")+string(9b),lc*0+stdev,lc*0+1,textout=strtrim(epic,2)+'_full_'+type+'.txt',/nocomment
  endif     
  base_time = time
  base_lc = lc

  if epic eq 203848625 then plotter,per,time,lc,bintime=bintime,binlc=binlc,err=binlc_err,bins=bins,range=range,epic=epic,type=type,fullrange=fullrange,print=print,charsize=charsize,charthick=charthick,yrange=[0.996,1.005] else $
     plotter,per,time,lc,bintime=bintime,binlc=binlc,err=binlc_err,bins=bins,range=range,epic=epic,type=type,fullrange=fullrange,yrange=yrange,print=print,charsize=charsize,charthick=charthick
  if debug eq 0 then device,/close

  ;;if debug eq 0 then device,filename='Vetfiles2/'+strtrim(string(epic,format="(I10)"),2)+'_vet2_'+type+'.eps',/encapsul,/color,xsize=20,ysize=10
  ;;time = base_time
  ;;lc = base_lc
  ;;if debug eq 0 then !p.multi=[0,1,2]

  ;;if fullrange ne 1 then yrange = -1
  ;;!y.margin=[2,1]
  ;;evenodd,per,time,lc,1,bintime=bintime,binlc=binlc,err=binlc_err,bins=bins,range=range,epic=epic,type=type,yrange=;;yrange,nox=1
  ;;time = base_time
  ;;lc = base_lc
  ;;!y.margin=[5,-2]
  ;;evenodd,per,time,lc,0,bintime=bintime,binlc=binlc,err=binlc_err,bins=bins,range=range,epic=epic,type=type,yrange=yrange
  ;;lc = base_lc
  
END
