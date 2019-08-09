##############################################################################
#####   TREE-STRUCTURED SCALE EFFECTS IN BINARY AND ORDINAL REGRESSION   #####
#####                     Electronic Supplement                          #####
##############################################################################
#####                    Author: Moritz Berger                           #####
##############################################################################

# Function to fit a tree-structured location-scale model 

### Input: 
# y: Response vector 
# DM_kov: Data.frame of covariates 
# alpha: Global significance level for the permutation tests 
# nperm: Number of permutations used for the permutation tests 
# scale_comp: If true (default), the scaling component is included in the model; 
#             otherwise, a simple proportional odds model is fitted 
# trace: If true (default), information about the estimation progress is printed 

### Description: 
# The function requires the R add-on package 'ordinal' and the self-implemented 
# functions in 'functions.R' 

### Output: 
# model:      Final model, fitted after the last iteration of tree-building 
# theta_hat:  Estimated threshold coefficients 
# beta_hat:   Estimated coefficients in the location term 
# gamma_hat:  Estimated coefficients in the variance term 
# splits:     Matrix with detailed information about all executed splits during 
#             the estimation process 
# pvalues:    P-values of each permutation test during tree-building
# devs:       Maximal value statistics T_j of the selected variables in each 
#             iteration of tree-building 
# crits:      Critical values of each permutation test during tree-building
# X:          Model matrix 
# y:          Response vector 


LocScaleTree <- function(y, 
                           DM_kov,
                           alpha, 
                           nperm, 
                           scale_comp=TRUE,
                           trace=TRUE){
  
  require(ordinal)
  
  n    <- length(y)
  nvar <- ncol(DM_kov)
  
  if(!is.null(names(DM_kov))){
    var_names <- names(DM_kov)
  } else{
    var_names <- paste0("x",1:nvar)
  }
  
  ordered_values <- lapply(1:nvar, function(j) ord_values(DM_kov[,j]))
  n_levels       <- sapply(ordered_values,length)
  thresholds     <- lapply(ordered_values,thresh)
  n_s            <- n_levels-1

  mod_potential <- list()
  devs          <- c()
  crits         <- c()
  pvalues       <- c() 
  splits        <- list() 
  splits[[1]] <- splits[[2]]  <- data.frame("variable"=numeric(),
                                            "split"=numeric(),
                                            "level"=numeric(),
                                            "node"=numeric(),
                                            "threshold"=numeric(),
                                            "number"=numeric(),
                                            "left"=numeric(),
                                            "right"=numeric(),stringsAsFactors=FALSE)
  params        <- list()
  vars_evtl     <- list()
  splits_evtl   <- list()
  which_obs     <- list()
  numbers       <- list()
  count         <- 1
  if(scale_comp){
    coI <- 1:2
  } else{
    coI <- 1 
  }
  
  params[[1]]      <- list("int","int")
  which_obs[[1]]   <- replicate(2, matrix(1:n,nrow=1),simplify=FALSE)
  vars_evtl[[1]]   <- replicate(2, nvar, simplify=FALSE)
  splits_evtl[[1]] <- replicate(ifelse(scale_comp,2,1), lapply(1:nvar, function(var) matrix(1:n_s[var],nrow=1)),simplify=FALSE)
  numbers[[1]]     <- replicate(2, 1, simplify=FALSE)
  
  dat0   <- data.frame("y"=as.factor(y),"int"=rep(1,n),DM_kov)
  mod0   <- mod_potential[[1]] <- clm(y~1, data=dat0, scale=~1)
  
  design_lower <- designlists(DM_kov,nvar,n_s,n_levels,ordered_values)[[1]]
  design_upper <- designlists(DM_kov,nvar,n_s,n_levels,ordered_values)[[2]]
  sig      <- TRUE
  anysplit <- TRUE
  
  while(sig & anysplit){
    
    # fit all models 
    dv <- lapply(1:nvar,function(var){
            lapply(coI, function(co){
              n_knots   <- length(params[[count]][[co]])
              deviances <- matrix(rep(0,n_s[var]*n_knots),ncol=n_knots)
              for(kn in 1:n_knots){
                deviances[,kn] <- allmodels(var,co,kn,count,design_lower,design_upper,splits_evtl,params,dat0,mod0,n_s)
              }
              return(deviances)
            })
          })

    
    # select optimum 
    variable <- which.max(lapply(1:nvar,function(j) max(unlist(dv[[j]]))))
    comp     <- which.max(lapply(coI, function(j) max(dv[[variable]][[j]])))
    split    <- as.numeric(which(dv[[variable]][[comp]]==max(dv[[variable]][[comp]]),arr.ind=TRUE)[,1])
    knoten   <- as.numeric(which(dv[[variable]][[comp]]==max(dv[[variable]][[comp]]),arr.ind=TRUE)[,2])
    if(length(split)>1){
      split  <- split[1]
      knoten <- knoten[1]
      warning(paste("Maximum in iteration ",count," not uniquely defined"))
    }
    param_old   <- params[[count]][[comp]][knoten]
    level       <- length(strsplit(param_old,":")[[1]])
    number      <- numbers[[count]][[comp]][knoten]
    left        <- max(numbers[[count]][[comp]])+1
    right       <- max(numbers[[count]][[comp]])+2
    
    param_new   <- paste(param_old,c(colnames(design_lower[[variable]])[split],colnames(design_upper[[variable]])[split]),sep=":")
    
    
    # compute permutation test 
    dev <- rep(NA,nperm)
    
    for(perm in 1:nperm){
      dev[perm] <- one_permutation(variable,comp,knoten,count,nvar,n_levels,ordered_values,
                                   DM_kov,which_obs,splits_evtl,params,dat0,mod0,n_s)
      if(trace){
        cat(".")
      }
    }
    
    # test decision 
    adaption <- vars_evtl[[count]][[comp]][knoten]
    crit_val <- quantile(dev,1-(alpha/adaption))
    Tj       <- max(dv[[variable]][[comp]])
    proof    <- Tj > crit_val
    devs[count]    <- Tj
    crits[count]    <- crit_val
    pvalues[count] <- sum(dev>Tj)/nperm
    
    if(proof){
      
      # fit new model 
      mod0  <- mod_potential[[count+1]] <- one_model(variable,comp,knoten,count,split,design_lower,design_upper,params,dat0)

      # adjust node 
      if(level>1){
        help_kn4 <- lu(c(),1,level-1,c())
        help_kn5 <- unlist(strsplit(param_old,""))
        help_kn6 <- paste0(help_kn5[which(help_kn5=="_")+1],collapse="")
        knoten2  <- which(help_kn4==help_kn6)
      } else{
        knoten2 <- knoten
      }
      
      splits[[comp]][count,"variable"] <- variable
      splits[[comp]][count,"split"] <- split
      splits[[comp]][count,"level"] <- level
      splits[[comp]][count,"node"]  <- knoten2
      splits[[comp]][count,"threshold"] <- thresholds[[variable]][[split]]
      splits[[comp]][count,c(6,7,8)] <- c(number,left,right)
      
      params[[count+1]]                             <- params[[count]]
      params[[count+1]][[comp]]                     <- rep("",length(params[[count]][[comp]])+1)
      params[[count+1]][[comp]][c(knoten,knoten+1)] <- param_new
      params[[count+1]][[comp]][-c(knoten,knoten+1)]<- params[[count]][[comp]][-knoten]
      
      n_knots                                                       <- length(params[[count+1]][[comp]])
      splits_evtl[[count+1]]                                        <- splits_evtl[[count]]
      for(var in 1:nvar){
        splits_evtl[[count+1]][[comp]][[var]]                       <- matrix(0,nrow=n_knots,ncol=n_s[var])
        splits_evtl[[count+1]][[comp]][[var]][c(knoten,knoten+1),]  <- matrix(rep(splits_evtl[[count]][[comp]][[var]][knoten,],2),nrow=2,byrow=T)
        splits_evtl[[count+1]][[comp]][[var]][-c(knoten,knoten+1),] <- splits_evtl[[count]][[comp]][[var]][-knoten,]
      }
      splits_evtl[[count+1]][[comp]][[variable]][knoten,splits_evtl[[count+1]][[comp]][[variable]][knoten,]>=split] <- NA 
      splits_evtl[[count+1]][[comp]][[variable]][(knoten+1),splits_evtl[[count+1]][[comp]][[variable]][(knoten+1),]<=split] <- NA
      
      # any split? 
      anysplit <- !all(is.na(unlist(splits_evtl[[count+1]])))
      
      vars_evtl[[count+1]]                             <- vars_evtl[[count]]
      vars_evtl[[count+1]][[comp]]                     <- rep(0,n_knots)
      vars_evtl[[count+1]][[comp]][c(knoten,knoten+1)] <- rep(vars_evtl[[count]][[comp]][knoten],2)
      vars_evtl[[count+1]][[comp]][-c(knoten,knoten+1)]<- vars_evtl[[count]][[comp]][-knoten]
      
      if(length(which(!is.na(splits_evtl[[count+1]][[comp]][[variable]][knoten,])))==0){ 
        vars_evtl[[count+1]][[comp]][knoten] <- vars_evtl[[count+1]][[comp]][knoten]-1 
      }
      if(length(which(!is.na(splits_evtl[[count+1]][[comp]][[variable]][knoten+1,])))==0){ 
        vars_evtl[[count+1]][[comp]][knoten+1] <- vars_evtl[[count+1]][[comp]][knoten+1]-1 
      }
      
      which_obs[[count+1]]                               <- which_obs[[count]]
      which_obs[[count+1]][[comp]]                       <- matrix(0,nrow=n_knots,ncol=n)
      which_obs[[count+1]][[comp]][c(knoten,knoten+1),]  <- matrix(rep(which_obs[[count]][[comp]][knoten,],2),nrow=2,byrow=T)
      which_obs[[count+1]][[comp]][-c(knoten,knoten+1),] <- which_obs[[count]][[comp]][-knoten,]
      threshh <- ordered_values[[variable]][1:n_s[variable]][split]
      which_obs[[count+1]][[comp]][knoten,DM_kov[,variable]>threshh] <- NA
      which_obs[[count+1]][[comp]][(knoten+1),DM_kov[,variable]<=threshh] <- NA
      
      numbers[[count+1]]                                  <- numbers[[count]]
      numbers[[count+1]][[comp]]                      <- numeric(length=n_knots)
      numbers[[count+1]][[comp]][c(knoten,knoten+1)]  <- c(left,right)
      numbers[[count+1]][[comp]][-c(knoten,knoten+1)] <- numbers[[count]][[comp]][-knoten] 
      
      # trace
      if(trace){
        cat(paste0("\n Split"," ",count,";"," ","Component"," ",comp,"\n"))
      }
      
      count <- count+1 
      
    } else{
      sig <- FALSE
    }
  }
  
  ################################################################################### 
  
  mod_opt     <- mod_potential[[count]]
  params_opt  <- params[[count]]
  
  theta_hat   <- mod_opt$alpha
  
  if(count>1){
    if(length(mod_opt$beta)>0){
      beta_hat <- mod_opt$beta
      beta_hat[is.na(beta_hat)] <- 0
      beta_hat <- beta_hat[params_opt[[1]]]
    } else{
      beta_hat <- NULL
    }
    if(length(mod_opt$zeta)>0){
      gamma_hat <- mod_opt$zeta
      gamma_hat[is.na(gamma_hat)] <- 0 
      gamma_hat <- gamma_hat[params_opt[[2]]]
    } else{
      gamma_hat <- NULL
    }
  } else{
    beta_hat  <- NULL
    gamma_hat <- NULL
  }
  
  names(splits) <- c("beta","gamma")
  for(i in 1:2){
    splits[[i]] <- na.omit(splits[[i]])
  }
  splits$beta$variable  <- names(DM_kov)[splits$beta$variable]
  splits$gamma$variable <- names(DM_kov)[splits$gamma$variable]
  
  to_return <- (list("model"=mod_opt,
                     "theta_hat"=theta_hat,
                     "beta_hat"=beta_hat,
                     "gamma_hat"=gamma_hat,
                     "splits"=splits,
                     "pvalues"=pvalues,
                     "devs"=devs,
                     "crits"=crits,
                     "X"=DM_kov,
                     "y"=y,
                     "call"=match.call()))
  
  class(to_return) <- "LocScaleTree"
  return(to_return)
  
}