######## Functions for spectrum deisotoping #######################################

## Function for counting the number of C, N, O, H atoms based on the numbers of Hex, HexNAc, dHex and NeuAc in a permethylated glycan 
atom.count=function(hex, dhex, hexnac, neuac)
{
  nc=hex*9+dhex*8+hexnac*11+neuac*16+2
  nh=hex*16+dhex*14+hexnac*19+neuac*27+6
  no=hex*5+dhex*4+hexnac*5+neuac*8+1
  nn=hexnac+neuac
  res=data.frame(nc,nh,no,nn)
  return(res)
}






## Function for calculating isotopic peak probability. k=0 corresponds to the first peak
p.iso.peak=function(nc,nh,no,nn,PC.1,PH.1,PO.1,PO.2,PN.1, k)
{
  ## probability of isotopic peak PK
  Pk=0
  for (C.1 in 0:k)
  {
    for (H.1 in 0:k)
    {
      for (N.1 in 0:k)
      {
        for (O.1 in 0:k)
        {
          for (O.2 in 0:k)
          {
            if ((C.1+H.1+N.1+O.1+2*O.2)==k)
            {
              Pk=Pk+dbinom(x=C.1, size = nc, prob = PC.1) *
                dbinom(x=H.1, size = nh, prob = PH.1) *
                dbinom(x=N.1, size = nn, prob = PN.1) *
                dmultinom(x=c((no-O.1-O.2),O.1,O.2),prob = c(1-PO.1-PO.2,PO.1,PO.2))
            }   
          }
        }
      }
    }
  }
  return(Pk)
}





## Function for calculating monosacharride composition of a glycan based on the mz value 
glycan.composition.finder=function(input.data=glycan.database, mz=3776.841)
{
  ## set glycan composition on core structures
  glycan.composition=c(0,4,0,5)
  
  ## name the monosaccharides
  
  names(glycan.composition)=c("dHex.count","HexNAc.count",
                              "NeuAc.count","Hex.count")
  
  
  ## find terminal modifications from glycan database, add to core structures
  terminal = subset.data.frame(input.data, 
                               abs(mass-mz)<0.3)
  terminal=terminal[1,]
  
  ## calculate the number of total monosaccharides
  if (terminal$Glycan.type %in% c('Polyhexose'))
  {
    glycan.composition=c(terminal$dHex.count,
                         terminal$HexNAc.count,
                         terminal$NeuAc.count,
                         terminal$Hex.count)
  }
    

  else if (terminal$Glycan.type %in% c('Hexa'))
  {
    glycan.composition=glycan.composition+c(terminal$dHex.count,
                                            terminal$HexNAc.count+4,
                                            terminal$NeuAc.count,
                                            terminal$Hex.count+4)
  }
    
  
  else if (terminal$Glycan.type %in% c('Penta'))
  {
    glycan.composition=glycan.composition+c(terminal$dHex.count,
                                            terminal$HexNAc.count+3,
                                            terminal$NeuAc.count,
                                            terminal$Hex.count+3)
  }
    
    
  else if (terminal$Glycan.type %in% c('Tetra'))
  {
    glycan.composition=glycan.composition+c(terminal$dHex.count,
                                            terminal$HexNAc.count+2,
                                            terminal$NeuAc.count,
                                            terminal$Hex.count+2)
  }
  else if (terminal$Glycan.type=='Tri')
  {
    glycan.composition=glycan.composition+c(terminal$dHex.count,
                                            terminal$HexNAc.count+1,
                                            terminal$NeuAc.count,
                                            terminal$Hex.count+1)
  }
  else 
  {
    glycan.composition=glycan.composition+c(terminal$dHex.count,
                                            terminal$HexNAc.count,
                                            terminal$NeuAc.count,
                                            terminal$Hex.count)
  }
  return(glycan.composition)
}





## Function for calculating isotopic peak distribution of a permethylated glycan. The first peak normalized to 1.
dstr.iso.peak=function (mz=3776.841)
{
  
  #### calculate the number of monosaccharides in a glycan
  glycan.composition=glycan.composition.finder(input.data=glycan.database, mz=mz)
  
  #### calculate number of atoms in a permethylated glycan
  glycan=atom.count(hex=glycan.composition[4],
                    hexnac=glycan.composition[2],
                    dhex=glycan.composition[1],
                    neuac=glycan.composition[3])
  
  #### distribution all peaks
  Intensities=c()
  for (i in 0:7)
  {
    Intensities=c(Intensities,
                  p.iso.peak(nc=glycan$nc,nh=glycan$nh,no=glycan$no,nn=glycan$nn,
                             PC.1=0.0107,PH.1=0.000115,PO.1=0.00038,
                             PO.2=0.00205, PN.1=0.00368,k=i))
  }
  
  Intensities.model=Intensities
  
  #### make the first peak as 1
  Intensities.model=Intensities.model/Intensities.model[1] 
  
  return(Intensities.model)
}





## Function for grouping the masses based on the mass differences, store the groups in a list
mass.grouper=function(masses=c(3211.568, 3213.583, 3215.588))
{
  masses=sort(masses)
  # add a number larger than the max number in masses
  # this is to ensure that the last number of the input masses will be 
  # processes in all cases
  masses=c(masses, max(masses)+100) 
  ans=list()
  i=1
  while (i<length(masses))
  {
    temp=c(masses[i])
    for (j in (i+1):length(masses))
    {
      if (min(abs(temp-masses[j]))<4.3)
      {
        temp=c(temp, masses[j])
      }
    }
    i=i+length(temp)
    ans=c(ans, list(temp))
  }
  return(ans)
}






## Deisotope function for a given spectrum and a mass group
deisotope=function(spectrum.number, mass.group=c(3211.568, 3213.583, 3215.588), 
                   number.of.peaks)
{
  ### experimental peaks
  spectrum.data=spectra.aligned[[spectrum.number]]
  Intensities=c()
  for (i in 0:(number.of.peaks-1))
  {
    Intensities=c(Intensities, 
       max(spectrum.data[abs(spectrum.data@mass-(mass.group[1]+i))<0.5]@intensity)-
       min(spectrum.data[abs(spectrum.data@mass-(mass.group[1]+i))<2]@intensity))
  }
  actual=Intensities
  
  if (sum(actual)==-Inf)
  {actual=rep(0,number.of.peaks)}
  
  ### predicted isotopic peak clusters for each mass.group, build a peak matrix
  gap=round(mass.group[2:length(mass.group)]-mass.group[1:(length(mass.group)-1)])
  gap=c(0,gap)
  peak.matrix=matrix(nrow = number.of.peaks, ncol = length(mass.group))
  for (i in 1:length(mass.group))
  {
    peak=dstr.iso.peak(mass.group[i])
    cluster=c(rep(0,sum(gap[1:i])), peak, rep(0,10))
    cluster=cluster[1:number.of.peaks]
    peak.matrix[,i]=cluster
  }
  
  ### build and solve the equation
  A=peak.matrix
  b=matrix(actual, ncol=1)
  ans=qr.solve(A,b)
  
  ### Calculate the predicted peaks
  predicted=A %*% ans
  predicted=as.vector(predicted)
  ans=as.vector(ans)
  
  ### calculate Rsquared
  Rsquared=summary(lm(predicted ~ actual,
                      data=data.frame(predicted,actual)))$r.squared 
  
  ### get the final answer
  ans=c(ans,Rsquared)
  
  ### plot the result
  plot(actual,xlab='', ylab='Intensity', type='h',
       ylim= c(0,max(actual)*1.2))
  
  par(new=T)
  
  plot(predicted, xlab=paste(sample.name.list[spectrum.number],'-',
                             sample.collection.time.list[spectrum.number]),
       ylab='Intensity', col= 'red',
       ylim=c(0, max(actual)*1.2), 
       sub=paste('R2=',round(Rsquared,digits = 2)),
       main=round(mass.group,1))
 
  
  ## optional: plot predicted individual peak cluster  
  #for (i in 1:length(mass.group))
  #{
  #  x=A[,i]*ans[i]
  #  x[x<0.01]=NA
  #  par(new=T)
  #  plot(x, ylab='Intensity', col='blue',pch=4+i,
  #       ylim=c(0, max(actual)*1.2))
  #}
  
  return(ans)
}


