##############################################################################
#####   TREE-STRUCTURED SCALE EFFECTS IN BINARY AND ORDINAL REGRESSION   #####
#####                     Electronic Supplement                          #####
##############################################################################
#####                    Author: Moritz Berger                           #####
##############################################################################

# The following contains auxiliary functions called in LocScaleTree.R to fit 
# the tree-structured location-scale model. 

# Note: these functions are not thought to be used inpendent of 
#       LocScaleTree(). 

# adjust node 
lu <- function(last=last, cd=cd, d=d, erg){
  
  s <- c("l","u")
  
  for(i in 1:length(s)){
    s_new <- paste0(last,s[i])
    if(cd==d){
      erg[length(erg)+1] <- s_new
    } else{
      erg <- lu(s_new,cd+1,d,erg)
    }
  }
  return(erg)
}

# modify factors
mod_factors <- function(y, x){
  tab <- table(x,y)
  nx  <- rowSums(tab)
  ptab <- tab/nx
  pp  <- colMeans(ptab)
  ptabc <- t(apply(ptab,1, function(x) x-pp))
  sig   <- matrix(0, nrow=ncol(tab), ncol=ncol(tab))
  for(j in 1:nrow(tab)){
    sig <- sig+nx[j]*ptabc[j,]%*%t(ptabc[j,])
  }
  sig <- sig/(length(y)-1)
  v <- eigen(sig)$vectors[,1]
  sa <- apply(ptabc, 1, function(x) v%*%x) 
  xnew <- as.numeric(x)
  xnew <- sapply(1:length(xnew), function(j) which(order(sa)==xnew[j]))
  return(xnew)
}

# compute ordered values 
ord_values <- function(x){
  if(!all((x - round(x)) == 0) || length(unique(x))>50){
    ret <- quantile(x,seq(0.05,1,by=0.05))
  } else{
    ret <- unique(sort(x))
  }
  return(ret)
}

thresh <- function(ordered_values){
  ret <- ordered_values[-length(ordered_values)]
  return(ret)
}

# functions to build design 
design_one  <- function(x,threshold,upper){
  if(upper){
    ret <- ifelse(x > threshold,1,0)
  } else{
    ret <- ifelse(x > threshold,0,1)
  }
  return(ret)
}

design <- function(x,thresholds,upper){
  ret <- sapply(thresholds, function(j) design_one(x,j,upper))
  return(ret)
}

designlist <- function(X,vars,label,thresholds,var_names,upper=TRUE){
  ret <- lapply(vars, function(j) {
    ret <- design(X[,j],thresholds[[j]],upper)
    colnames(ret) <- paste0("s",which(var_names==j),1:length(thresholds[[j]]),ifelse(upper,"_u","_l"),label)
    return(ret)})
  names(ret) <- vars 
  return(ret)
}

designlists <- function(DM_kov,nvar,n_s,n_levels,ordered_values){
  
  v <- lapply(1:nvar,function(j) 1:(n_levels[j]-1)) 
  w <- lapply(1:nvar, function(j) rep(paste0("s",j),n_s[j]))
  
  design_upper <- lapply(1:nvar, function(j){
    design_matrix <- sapply(1:(n_levels[j]-1),function(k) { ifelse(DM_kov[,j] > ordered_values[[j]][k],1,0)})
    colnames(design_matrix) <- paste0(w[[j]],v[[j]],"_u")
    design_matrix
  })
  design_lower <- lapply(1:nvar, function(j){
    design_matrix <- abs(design_upper[[j]]-1)
    colnames(design_matrix) <- paste0(w[[j]],v[[j]],"_l")
    design_matrix
  })
  return(list(design_lower,design_upper))
}

# functions to fit models 
one_model <- function(var,co,kn,count,j,design_lower,design_upper,params,dat0){
  
  dat   <- data.frame(dat0,do.call(cbind,design_lower),do.call(cbind,design_upper))

  help1 <- params[[count]][[co]]
  help2 <- params[[count]][[-co]]
  help22<- paste(help2, collapse="+")
  help3 <- help1[-kn]
  help4 <- paste(help1[kn],c(colnames(design_lower[[var]])[j],colnames(design_upper[[var]])[j]),sep=":")
  help5 <- paste(c(help3,help4), collapse="+")
  if(co==1){
    help6 <- formula(paste("y~", help5))
    help7 <- formula(paste("~", help22))
  } 
  if(co==2){
    help6 <- formula(paste("y~", help22))
    help7 <- formula(paste("~", help5))
  }
  mod <- suppressWarnings(clm(help6, data=dat, scale=help7))
  return(mod)
}

allmodels <- function(var,co,kn,count,design_lower,design_upper,splits_evtl,params,dat0,mod0,n_s){
  
  deviances <- rep(0,n_s[var])
  splits_aktuell <- splits_evtl[[count]][[co]][[var]][kn,]
  splits_aktuell <- splits_aktuell[!is.na(splits_aktuell)]
  
  if(length(splits_aktuell)>0){
    for(j in splits_aktuell){
      mod <- one_model(var,co,kn,count,j,design_lower,design_upper,params,dat0)
      deviances[j] <- (-2)*logLik(mod0) - (-2)*logLik(mod)
    }
  }
  return(deviances)
}

one_permutation <- function(var,co,kn,count,nvar,n_levels,ordered_values,
                            DM_kov,which_obs,splits_evtl,params,dat0,mod0,n_s){
  
  obs_aktuell <- which_obs[[count]][[co]][kn,]
  obs_aktuell <- obs_aktuell[!is.na(obs_aktuell)]
  DM_kov_perm <- DM_kov
  DM_kov_perm[obs_aktuell,var] <- sample(DM_kov_perm[obs_aktuell,var],length(obs_aktuell))
  
  design_upper_perm      <- designlists(DM_kov_perm,nvar,n_s,n_levels,ordered_values)[[1]]
  design_lower_perm      <- designlists(DM_kov_perm,nvar,n_s,n_levels,ordered_values)[[2]]
  
  dv_perm <- allmodels(var,co,kn,count,design_lower_perm,design_upper_perm,splits_evtl,params,dat0,mod0,n_s)
  
  return(max(dv_perm))
  
}