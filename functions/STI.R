### Indici G e Gstar ###
star = FALSE

# Pseudo-Kronecker function
pkron = function(A, B){
  a = A[1,]
  nA = length(a)
  nB = NCOL(B)
  RES = matrix(NA, NROW(B), nA*NCOL(B))
  counter = 0
  for(i in 1:nA){
    RES[,(i-1)*nB + (1:nB)] = a[i] * B
  }
  return(RES)
}

# Spatial weight matrix
GetW = function(coords,
                id.name,
                buffer = 15,
                include.self = TRUE){
  ncells = length(unique(coords$id))
  W = matrix(0, ncells, ncells)
  ids = sort(unique(coords$id))
  for(i in 1:ncells){
    idists = spDistsN1(pts = data.matrix(coords[,c("lon", "lat")]),
                       pt = as.numeric(coords[which(coords[,id.name] ==
                                                      ids[i]),c("lon", "lat")]),
                       longlat = TRUE)
    if(include.self == TRUE){
      idists = coords[which(idists < buffer), id.name]
    }else{
      idists = setdiff(coords[which(idists < buffer), id.name], ids[i])
    }
    W[i, match(idists, ids)] = 1
  }
  return(W)
}

# Space–time connectivity matrix
STCM = function(spatCon, nCells, nTimes, include.self=TRUE){
  IT = diag(nTimes)
  Is = diag(nCells)
  WT = 1 - diag(nTimes)
  if(include.self){
    VST = pkron(IT, spatCon) + pkron(WT, Is)
    # VST = kronecker(IT, spatCon) + kronecker(WT, Is)
  } else {
    VST = pkron(WT, spatCon) + pkron(WT, Is)
    # VST = kronecker(WT, spatCon) + kronecker(WT, Is)
  }
  return(VST)
}

G = function(stConn, values, nCells, include.self=TRUE){
  num = den = numeric(nCells)
  TT = length(values) / nCells
  for(i in 1:nCells){
    num[i] = sum(stConn[i,] * values)
  }
  if(include.self){
    den = sum(values)
  } else {
    Xj = matrix(values, nrow=TT, byrow=T)
    den = numeric(nCells)
    for(iCell in 1:nCells){
      den[iCell] = sum(Xj[,-iCell])
    }
  }
  res = num / den
  return(res)
}

Gstandard = function(stConn, values, nCells, include.self=TRUE){
  num = den = numeric(nCells)
  TT = length(values) / nCells
  N = nCells * TT
  if(!include.self){
    for(i in 1:nCells){
      xbari = sum(values[-i]) / (N-1)
      Vi = sum(stConn[i,-i])
      SN1i = sum(stConn[i,-i]^2)
      s2Ni = sum(values[-i]^2) / (N-1) - xbari^2
      num[i] = sum(stConn[i,] * values) - Vi * xbari
      den[i] = sqrt(s2Ni * (((N-1) * SN1i - Vi^2) / (N-2)))
    }
    res = num / den
  } else {
    xbaristar = sum(values) / N
    s2Nistar = sum(values^2) / N - xbaristar^2
    for(i in 1:nCells){
      Vistar = sum(stConn[i,])
      SN1istar = sum(stConn[i,]^2)
      num[i] = sum(stConn[i,] * values) - Vistar * xbaristar
      den[i] = sqrt(s2Nistar * ((N * SN1istar - Vistar^2) / (N-1)))
    }
    res = num / den
  }
  res = num / den
  return(res)
}
