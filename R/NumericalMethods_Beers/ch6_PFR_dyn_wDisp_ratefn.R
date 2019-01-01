
ratefn = function(t,c,parms){
  
  D_uL = parms$D_uL
  ngrid = parms$ngrid
  nc = parms$nc
  cinlet = parms$cinlet
  k1 = parms$k1
  k2 = parms$k2
  
  amp_pct = parms$amp_pct
  amp_freq = parms$amp_freq
  
  incFluct = parms$incFluct
  
  if(t >= tstart & incFluct == 1){
    cB_ss = cinlet[2]
    cinlet[2] = cB_ss + cB_ss * amp_pct * sin((t - tstart)*2*pi/amp_freq)
  }
  
  dz = 1/(ngrid + 1)
  
  alpha_lo = 1/dz + D_uL/dz^2
  alpha_mid = -1/dz - 2*D_uL/dz^2
  alpha_hi = D_uL/dz^2
  
  PeLoc = dz/D_uL
  
  conc = matrix(rep(0,ngrid * nc),ncol = nc)
  dconc = matrix(rep(0,ngrid * nc),ncol = nc)
  rxnrate = matrix(rep(0,ngrid * nc),ncol = nc)
  
  
  k = 0
  for(i in 1:ngrid){
    for(j in 1:nc){
      k = k + 1
      conc[i,j] = c[k]
    }
  }
  
  for(i in 1:ngrid){
    rxnrate[i,1] = -k1 * conc[i,1] * conc[i,2]
    rxnrate[i,2] = -k1 * conc[i,1] * conc[i,2] - k2 * conc[i,2] * conc[i,3]
    rxnrate[i,3] = k1 * conc[i,1] * conc[i,2] - k2 * conc[i,2] * conc[i,3]
    rxnrate[i,4] = k2 * conc[i,2] * conc[i,3]
  }
  
  for(i in 2:(ngrid-1)){
    for (j in 1:nc){
        dconc[i,j] = alpha_lo * conc[i-1,j] + alpha_mid * conc[i,j] + alpha_hi * conc[i+1,j] + rxnrate[i,j]     
    }
  }
  
  for(j in 1:nc){
    
    dconc[1,j] = (alpha_mid + alpha_lo/(1+PeLoc))*conc[1,j] + 
      (alpha_lo * PeLoc * cinlet[j]/(1+PeLoc)) + alpha_hi * conc[2,j] + rxnrate[1,j]
    
    dconc[ngrid,j] = alpha_lo * conc[ngrid-1,j] + (alpha_mid + alpha_hi) * conc[ngrid,j] + rxnrate[ngrid,j]
    
  }
  
  dc = rep(0,ngrid*nc)
  k = 0
  for(i in 1:ngrid){
    for (j in 1:nc){
      k = k + 1
      dc[k] = dconc[i,j]
    }
  }
  
  return(list(dc))
  
}