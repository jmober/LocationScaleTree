##############################################################################
#####   TREE-STRUCTURED SCALE EFFECTS IN BINARY AND ORDINAL REGRESSION   #####
#####                     Electronic Supplement                          #####
##############################################################################
#####                    Author: Moritz Berger                           #####
##############################################################################

# Function to plot the trees of a fitted tree-structured location-scale model
# by LocScaleTree()

### Input: 
# x:            Fitted tree; object of class 'LocScaleTree'
# component:    Component of the model to be plotted; 
#               can be 'beta' (location term) or 'gamma' (variance term)
# ellipse_a,
# ellipse_b:    Values passed to draw_ellipse(); radii of the ellypses 
# ellipse_x,
# ellipse_y:    Values passed to draw_ellipse(); adjustment of the position 
#               of the ellypses in x- and y-direction 
# branch_adj:   (Optional) value to adjust the postion of the labels 
#               in y-direction 
# cex.lines:    Width of branches of the tree
# cex.branches: Size of labels of the tree
# cex.coefs:    Size of coefficients in the terminal nodes of the tree
# cex.main:     Size of the title of the tree
# cex.numbers:  Size of the numbers of the nodes 
# draw_numbers: If true (default), the numbers of the nodes are printed
# title:        (Optional) title, which is added to the tree

plot.LocScaleTree <- function(x, 
                              component=c("beta","gamma"), 
                              ellipse_a=0.8, 
                              ellipse_b=0.2, 
                              ellipse_x=0,
                              ellipse_y=0.2,
                              branch_adj=0,
                              cex.lines=2,
                              cex.branches=1,
                              cex.coefs=1,
                              cex.main=1,
                              cex.numbers=1, 
                              draw_numbers=TRUE,
                              title=NULL){
  
  if(is.null(x$splits)){
    cat("There is no plot available")
  } else{
    
    X        <- x$X 
    info     <- x$splits[[component]]
    
    if(nrow(info)==0){
      cat("There is no tree to plot.\n")
    } else{
      
      if(component=="beta"){
        coefs_hat <- x$beta_hat
      } 
      if(component=="gamma"){
        coefs_hat <- x$gamma_hat
      }
      
      require(plotrix)
      
      endnodes      <- list()
      endnodes[[1]] <- 1
      for(j in 1:nrow(info)){
        endnodes[[j+1]] <- numeric(length=(j+1))
        what <- c(info[j,"left"],info[j,"right"])
        delete <- info[j,"number"]
        where  <- which(endnodes[[j]]==delete)
        endnodes[[j+1]][c(where,where+1)] <- what
        endnodes[[j+1]][-c(where,where+1)] <- endnodes[[j]][-where]
      }
      endnodes <- endnodes[[length(endnodes)]]
      dir      <- sapply(1:length(endnodes), function(j){ifelse(endnodes[j] %in% info[,7], "l","r")})
      
      n_levels <- length(unique(info[,"level"]))
      
      hilfspunkte <- list()
      hilfspunkte[[1]] <- matrix(NA,nrow=2^n_levels,ncol=2)
      hilfspunkte[[1]][,1] <- 2^n_levels
      hilfspunkte[[1]][,2] <- rep(n_levels+1,2^n_levels)
      
      steps <- 2^((n_levels:1-1))
      
      for(i in 1:n_levels){
        
        hilfspunkte[[i+1]] <- hilfspunkte[[i]]
        hilfspunkte[[i+1]][,2] <- rep(n_levels+1-i,2^n_levels)
        
        help  <- c(-steps[i],steps[i])
        help1 <- rep(help,each=steps[i])
        help2 <- rep(help1,length=2^n_levels) 
        hilfspunkte[[i+1]][,1] <- hilfspunkte[[i]][,1]+help2 
        
        which_knots <- info[info[,"level"]==i,"node"]
        help3 <- seq(1,2^n_levels)
        help4 <- split(help3,rep(1:2^(i-1),each=2^n_levels/2^(i-1)))
        help5 <- unlist(lapply(which_knots, function(j) help4[[j]]))
        hilfspunkte[[i+1]][-help5,] <- hilfspunkte[[i]][-help5,]
        
      }
      
      
      plot.new()
      plot.window(ylim=c(0.5,n_levels+1),xlim=c(0,2^(n_levels+1)))
      
      for(j in length(hilfspunkte):2){
        for(i in 1:(2^n_levels)){
          lines(c(hilfspunkte[[j-1]][i,1],hilfspunkte[[j]][i,1]),c(hilfspunkte[[j-1]][i,2],hilfspunkte[[j]][i,2]),
                lwd=cex.lines)
        }
      }
      if(is.null(title)){
        title <- component
      }
      title(title,cex.main=cex.main)
      
      betas_hat <- format(round(coefs_hat,3),nsmall=3)
      points_betas <- unique(hilfspunkte[[n_levels+1]])
      for(i in 1:length(betas_hat)){
        if(dir[i]=="l"){
          fac <- -1
        } else{
          fac <- 1
        }
        draw.ellipse(x=points_betas[i,1]+fac*ellipse_x,y=points_betas[i,2]-ellipse_y,a=ellipse_a,b=ellipse_b,lwd=cex.lines,col=grey(0.8))
        text(points_betas[i,1]+fac*ellipse_x,points_betas[i,2]-ellipse_y,betas_hat[i],cex=cex.coefs)
      }
      
      if(draw_numbers){
        for(i in 1:length(betas_hat)){
          if(dir[i]=="l"){
            fac <- -1
          } else{
            fac <- 1
          }
          rect(points_betas[i,1]+fac*ellipse_x-max(0.2,ellipse_x/3),points_betas[i,2]+ellipse_b-max(0.05,ellipse_b/3),points_betas[i,1]+fac*ellipse_x+max(0.2,ellipse_x/3),points_betas[i,2]+ellipse_b+max(0.05,ellipse_b/3),col=grey(0.9),lwd=cex.numbers)
          text(points_betas[i,1]+fac*ellipse_x,points_betas[i,2]+ellipse_b,endnodes[i],cex=cex.numbers)
        }
      }

      for(i in 1:nrow(info)){
        help4 <- split(help3,rep(1:2^(info[i,"level"]-1),each=2^n_levels/2^(info[i,"level"]-1)))[[info[i,"node"]]]
        point_var <- unique(hilfspunkte[[info[i,"level"]]][help4,])
        points(point_var[1],point_var[2],cex=cex.lines-1,pch=19)
        point_left  <- c(point_var[1]-steps[info[i,"level"]]-branch_adj,point_var[2]-0.5)
        point_right <- c(point_var[1]+steps[info[i,"level"]]+branch_adj,point_var[2]-0.5)
        var   <- info[i,"variable"]
        thres <- info[i,"threshold"]
        sort_values <- unique(sort(X[,var]))
        if(thres==min(sort_values)){
          text(point_left[1],point_left[2],paste0(var,"=",round(thres,2)),cex=cex.branches,adj=c(1,0))
        } else{
          text(point_left[1],point_left[2],paste0(var,"<=",round(thres,2)),cex=cex.branches,adj=c(1,0))
        }
        if(thres==max(sort_values[-length(sort_values)])){
          text(point_right[1],point_right[2],paste0(var,"=",round(max(sort_values),2)),cex=cex.branches,adj=c(0,0))
        } else{
          text(point_right[1],point_right[2],paste0(var,">",round(thres,2)),cex=cex.branches,adj=c(0,0))
        }
      }
    }
  }
}