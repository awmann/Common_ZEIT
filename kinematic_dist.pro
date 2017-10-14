PRO kinematic_dist

  pmra = 120.55     ;; ppmxl
  pmra_err = 3.3    ;; ppmxl
  pmdec = -21.09    ;; ppmxl
  pmdec_err = 3.2   ;; ppmxl
  ;;distance = 44.0   ;; photometric distance
  ;;distance_err = 7. ;; error on photometric distance
  ;vrad = 39.28      ;; average of Mace values
  ;rv_err = 0.07      ;; stdev(rv)/sqrt(n_elements(rv))
  rv = [39.33603,39.09584,39.13028,39.56312];[39.336,39.326,39.316,39.306];
  rv_err = [0.094348156,0.10122231,0.091012937,0.11221297]
  vrad = mean(rv)
  vrad_err = stdev(rv)/sqrt(n_elements(rv))
  ra = 63.273378
  dec = 15.247777

  ;;39.76+/-0.1
  vrad = 39.76d0
  vrad_err = 0.1d0
  ra = am_racnv('04 29 38.993')
  dec = am_deccnv('+22 52 57.80')
  pmra = 83.0;;85.8;81.8
  pmra_err = 0.9;1.2;1.0
  pmdec = -35.7;-34.0;-35.2
  pmdec_err = 0.9;0.9

  uh = -43.1
  vh = -19.2
  wh = -1.4
  h_b = [uh,vh,wh]
  h_err = [0.9,0.2,0.4]
  h_err = sqrt(h_err^2.+1.2^2.);[1.1,1.1,1.1]
  nmonte = 10000.
  bigdists = dblarr(nmonte)
  for j = 0,nmonte-1 do begin
     h = h_b+h_err*randomn(seed,3)
     distance = generatearray(30.0,65.0,2000)
     chisqs = dblarr(n_elementS(distance))
     for i = 0,n_elements(distance)-1 do begin
        uvw_errors, ra, dec, pmra, pmdec, vrad, distance[i], 1./3600., 1./3600., pmra_err, pmdec_err, vrad_err, 0.0, U, V, W, E_U, E_V, E_W
        chisq = total(([u,v,w]-h)^2./(h_err^2+[e_u,e_v,e_w]^2.))
        chisqs[i] = chisq/3.
     endfor
     ;;plot,distance,chisqs,psym=8,/ylog
     ;;print,distance[where(chisqs eq min(chisqs))]
     ;;print,(sqrt(2./deriv(distance,deriv(distance,chisqs))))[where(chisqs eq min(chisqs))]
     bigdists[j] = distance[where(chisqs eq min(chisqs))]
     if j mod 100 eq 20 then print,median(bigdists[0:j-1]),stdev(bigdists[0:j-1]),min(chisqs)
  endfor

  cghistoplot,bigdists,/outline
  print,median(bigdists),stdev(bigdists)
  print,1000./median(bigdists),(1000./median(bigdists))*(stdev(bigdists)/median(bigdists))
  stop
     

END

 
function uvw,in,out

  
  uvw_errors, in[0], in[1], in[2], in[3], vrad, distance, 1./3600., 1./3600., pmra_err, pmdec_err, vrad_err, distance_err, U, V, W, E_U, E_V, E_W

END


PRO kinematic_dist2

  pmra = 23.5;18.4;
  pmra_err = 4.;7.
  pmdec = -37.9;-56.6;
  pmdec_err = 4.0
  base_dist = 133.0   
  vrad = 7.485014 ; 7.193193 , 7.4
  vrad_err = 0.175423565 ; 0.170971171, 0.2
  ra = 55.228426
  dec = 12.572612

  ;uh = -43.1
  ;vh = -19.2
  ;wh = -1.4
  ;h_b = [uh,vh,wh]
  ;h_err = [0.9,0.2,0.4]
  ;h_err = [1.1,1.1,1.1]

  h_b = [-6.7,-25.0,-12.8];,-107.0,25.9,-48.3]
  h_err = [2.0,2.0,2.0];[0.9,0.5,0.5];[1.0,1.0,1.0] ;

  nmonte = 1000.
  bigdists = dblarr(nmonte)
  for j = 0,nmonte-1 do begin
     h = h_b+h_err*randomn(seed,3)
     distance = generatearray(base_dist/1.5,base_dist*1.5,2500)
     chisqs = dblarr(n_elementS(distance))
     for i = 0,n_elements(distance)-1 do begin
        uvw_errors, ra, dec, pmra, pmdec, vrad, distance[i], 1./3600., 1./3600., pmra_err, pmdec_err, vrad_err, 0.0, U, V, W, E_U, E_V, E_W
        chisq = total(([u,v,w]-h)^2./(h_err^2+[e_u,e_v,e_w]^2.))
        chisqs[i] = chisq/3.
     endfor
     ;;plot,distance,chisqs,psym=8,/ylog
     ;;print,distance[where(chisqs eq min(chisqs))]
     ;;print,(sqrt(2./deriv(distance,deriv(distance,chisqs))))[where(chisqs eq min(chisqs))]
     bigdists[j] = distance[where(chisqs eq min(chisqs))]
     ;;print,min(chisqs)
  endfor

  cghistoplot,bigdists,/outline
  print,median(bigdists),stdev(bigdists)
  print,1000./median(bigdists),(1000./median(bigdists))*(stdev(bigdists)/median(bigdists))
  stop
     

END


function uvw,in,out

  
  uvw_errors, in[0], in[1], in[2], in[3], vrad, distance, 1./3600., 1./3600., pmra_err, pmdec_err, vrad_err, distance_err, U, V, W, E_U, E_V, E_W

END
