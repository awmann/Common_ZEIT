PRO wrapper
  
  epics = [211901114,211822797,211969807,211916756,211970147,211913977,211990866]
  for i = 0,n_elements(epics)-1 do prob_member_generic,epics[i],'Pr'
  
END


PRO prob_member_generic,epic,cluster

  if n_elements(cluster) eq 0 then cluster = 'all'
  case epic of
     210736056: begin
        distance = 100
        distance_err = 5
        vrad = 7.5
        vrad_err = 0.2
        ra = 54.397720
        dec = 18.896429
     end
     211990866: begin
        distance = 167.27075
        distance_err = 40.705452
        vrad = 33.60
        vrad_err = 0.30
        ra = 129.60125
        dec = 20.106063
     end
     211913977: begin
        distance = 171.30104
        distance_err = 22.223772
        vrad = 34.31
        vrad_err = 0.17
        ra = 130.344085
        dec =18.933875
     end
     211970147: begin
        distance = 168.34491
        distance_err = 22.900765
        vrad = 34.845
        vrad_err = 0.17
        ra = 130.056046
        dec = 19.77881
     end
     211822797: begin
        distance = 168.79697
        distance_err = 27.587216
        vrad = 34.847566
        vrad_err = 0.17
        ra = 130.4103394
        dec = 17.63999939
     end
     211969807: begin
        distance = 184.62053
        distance_err = 22.555972
        vrad = 34.807586
        vrad_err = 0.17
        ra = 129.636756
        dec = 19.773829
     end
     211901114: begin
        distance = 164.50404
        distance_err = 25.1
        vrad = 34.062246
        vrad_err = 0.17
        ra = 130.398723
        dec = 18.743057
     end
     211916756: begin
        distance = 210.22289
        distance_err = 44.952405
        vrad = 35.849206
        vrad_err = 0.166678682
        ra =129.362743
        dec = 18.976685
     end
     ;;212009427: begin
     ;;   distance = 184.62053
     ;;   distance_err = 22.555972
     ;;   vrad = -42.361596
     ;;   vrad_err = 0.17
     ;;   ra = 127.874459
     ;;   dec = 20.410422
     ;;end
     210363145: begin
        distance = -1
        distance_err = -1;;40.952405
        vrad = 7.48;7.4
        vrad_err = 0.17;0.2
        ra = 55.228426
        dec = 12.572612
     end
     210894022: begin
        distance = 44.0   ;; photometric distance
        distance_err = 7. ;; error on photometric distance
        ;vrad = 38637.730/1d3
        ;vrad_err = 153./1d3
        vrad = -16.3
        vrad_err = 1.0
        ra = 059.8897545
        dec = 21.2986850
     end
  endcase
  info = queryvizier('I/322A/out',[ra,dec],0.1,/all)
  if n_elements(info) gt 1 then begin
     gcirc,2,ra,dec,info.raj2000,info.dej2000,dist
     info = info[where(dist eq min(Dist))]
     info = info[0]
  endif
  if ~isarray(info) then begin
     info = queryvizier('V/139/sdss9',[ra,dec],0.1,/all)
     if n_elements(info) gt 1 then begin
        gcirc,2,ra,dec,info.raj2000,info.dej2000,dist
        info = info[where(dist eq min(Dist))]
        info = info[0]
     endif
  endif
  if ~isarray(info) then begin
     info = queryvizier('I/317/sample',[ra,dec],0.1,/all)
     pmra = info.pmra
     pmdec = info.pmde
     pmra_err = info.epRA
     pmdec_err = info.epDE
     stop
  endif else begin
     pmra = info.pmra
     pmdec = info.pmde
     pmra_err = info.e_pmra
     pmdec_err = info.e_pmde
  endelse
  xyz_errors, ra, dec, distance, 1/3600., 1/3600., distance_err, X, Y, Z, EX, EY, EZ
  uvw_errors, ra, dec, pmra, pmdec, vrad, distance, .1/3600., .1/3600., pmra_err, pmdec_err, vrad_err, distance_err, $
              U, V, W, E_U, E_V, E_W
  print,'UVW: ',u,v,w,' +/- ',e_u,e_v,e_w
  print,'XYZ: ',x,y,z,' +/- ',ex,ey,ez
  uvwxyz_star = [u,v,w,x,y,z]
  uvwxyz_star_err = [e_u,e_v,e_w,ex,ey,ez]

  uvwxyz_field = [-10.92,-13.35,-6.79,-0.18,2.10,3.27]
  uvwxyz_field_err = [23.22,13.44,8.97,53.29,51.29,50.70]

  if cluster eq 'Pl' or cluster eq 'all' then begin
     print,'Pleiades:'
     uvwxyz_p = [-6.7,-25.0,-12.8,-107.0,25.9,-48.3]
     uvwxyz_p_err = [2.0,2.0,2.0,9,8,5] 
     uvwxyz_tot_err = sqrt(uvwxyz_star_err^2.+uvwxyz_p_err^2.)
     P_th_H = 1./((sqrt(2*!pi)*uvwxyz_tot_err)) * exp((-1./2.)*((uvwxyz_star-uvwxyz_p)/uvwxyz_tot_err)^2.)
     uvwxyz_err_tot = sqrt(uvwxyz_star_err^2.+uvwxyz_field_err^2.)
     P_th_f = 1./((sqrt(2*!pi))* uvwxyz_err_tot) * exp((-1./2.)*((uvwxyz_star-uvwxyz_field)/uvwxyz_err_tot)^2.)
     a = 600.                      
     b = 10000.                
     tot_h = (p_th_h[0]*p_th_h[1]*p_th_h[2]*p_th_h[3]*p_th_h[4]*p_th_h[5])*(a/b)
     tot_f = (p_th_f[0]*p_th_f[1]*p_th_f[2]*p_th_f[3]*p_th_f[4]*p_th_f[5])*((b-a)/b)
     print,tot_h/(tot_f+tot_h)
  endif

  if cluster eq 'H' or cluster eq 'all' then begin
     print,'Hyades:'
     uvwxyz_p = [-41.1,-19.2,-1.4,-43.1,0.7,-17.3]
     uvwxyz_p_err = [2.0,2.0,2.0,9.3265732,8.5119028,5.8616913] 
     uvwxyz_tot_err = sqrt(uvwxyz_star_err^2.+uvwxyz_p_err^2.)
     P_th_H = 1./((sqrt(2*!pi)*uvwxyz_tot_err)) * exp((-1./2.)*((uvwxyz_star-uvwxyz_p)/uvwxyz_tot_err)^2.)
     uvwxyz_err_tot = sqrt(uvwxyz_star_err^2.+uvwxyz_field_err^2.)
     P_th_f = 1./((sqrt(2*!pi))* uvwxyz_err_tot) * exp((-1./2.)*((uvwxyz_star-uvwxyz_field)/uvwxyz_err_tot)^2.)
     a = 519.                      
     b = 21485.                
     tot_h = (p_th_h[0]*p_th_h[1]*p_th_h[2]*p_th_h[3]*p_th_h[4]*p_th_h[5])*(a/b)
     tot_f = (p_th_f[0]*p_th_f[1]*p_th_f[2]*p_th_f[3]*p_th_f[4]*p_th_f[5])*((b-a)/b)
     print,tot_h/(tot_f+tot_h)
  endif
  
  if cluster eq 'Pr' or cluster eq 'all' then begin
     print,'Praesepe:'
     uvwxyz_p = [-41.5,-19.8,-9.7,-138.4,-67.0,97.6]
     uvwxyz_p_err = [2.0,2.0,2.0,9,8,5] 
     uvwxyz_tot_err = sqrt(uvwxyz_star_err^2.+uvwxyz_p_err^2.)
     P_th_H = 1./((sqrt(2*!pi)*uvwxyz_tot_err)) * exp((-1./2.)*((uvwxyz_star-uvwxyz_p)/uvwxyz_tot_err)^2.)
     uvwxyz_err_tot = sqrt(uvwxyz_star_err^2.+uvwxyz_field_err^2.)
     P_th_f = 1./((sqrt(2*!pi))* uvwxyz_err_tot) * exp((-1./2.)*((uvwxyz_star-uvwxyz_field)/uvwxyz_err_tot)^2.)
     a = 519.                      
     b = 21485.                
     tot_h = (p_th_h[0]*p_th_h[1]*p_th_h[2]*p_th_h[3]*p_th_h[4]*p_th_h[5])*(a/b)
     tot_f = (p_th_f[0]*p_th_f[1]*p_th_f[2]*p_th_f[3]*p_th_f[4]*p_th_f[5])*((b-a)/b)
     print,tot_h/(tot_f+tot_h)
  endif
  
END
