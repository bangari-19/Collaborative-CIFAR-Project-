#230222: AR=.6; Contemp = .35; lagged = .45, noise=.01
rm(list=ls())
seed=123456
set.seed(seed)
#######################################################
#### 1. load functions needed to simulate the data ####
#######################################################

#Group Relations for Set 1: AR: V2, V4, V6; Contemporary Connections: V5 -> V3, V6 -> V4; Lagged Connections: V1lagged -> V6, V4lagged -> V5

mat.generate.asw <- function(p.con,nvar,AR,dens,p.group,con.b,lag.b){
  repeat{
    A <- matrix(0, ncol = nvar, nrow = nvar, )        ###create null contemporaneous matrix
    Phi <- matrix(0, ncol = nvar, nrow = nvar)        ###create null lagged matrix (remember -> ARs in diagonal)
    pos.all <- 2*(nrow(A)*ncol(A))                    ###list all available (non-ar) spots (You could add - 2*nrow(A) to remove AR spots (see AW code))
    cnt.all <- dens*p.group*pos.all                   ###number of group paths given density and proportion of group paths
    indices <- which(Phi == 0, arr.ind = TRUE)        ###create data frame that lists all combinations of columns/rows in matrices. Run "indices <- indices[which(indices[,1] != indices[,2]), ]" to kick out diagonal
    row.col<-c(17, 24, 31, 28)                        ###select group-level associations
    diag.set<-c(8, 22, 36)                            ###select AR paths 
    print(round(p.con*length(row.col)))
    n.p.1<-row.col[1:round(p.con*length(row.col))]    ###identify contemporaneous group paths 
    n.p.2<-row.col[(length(n.p.1)+1):length(row.col)] ###identify lagged group paths 
    grp.con <- n.p.1                                  ###name contemporaneous group paths
    grp.lag <- n.p.2                                  ###name lagged group paths
    grp.diag <- diag.set                              ###name AR paths
    Phi[indices[grp.lag,]] <- lag.b                   ###set lagged paths (Identified below)
    A[indices[grp.con,]] <- con.b                     ###set contemporaneous paths (Identified below)
    Phi[indices[grp.diag,]] <- AR                     ###set AR paths (Identified below)
    break      
    
  }
  
  all <- cbind(Phi, A)                                          ###combine contemporaneous and lagged matrix 
  ind.pres <- which(all != 0, arr.ind = T)                      ###identify the group paths 
  level <- "grp"
  all.lvl <- matrix(NA, ncol = ncol(all), nrow = nrow(all))
  all.lvl[ind.pres] <- level
  
  all_sub1 <- all
  all_lvl1 <- all.lvl
  
  res <- list(sub1 = all_sub1, 
              lvl1 = all_lvl1)
  return(res)
}

#This function generates group-level relations, adds individual-level relations to the matrix, 
#adds noise, and simulates the time series (unnecessary for our purposes but helpful for checks!)
ts.generate.asw <- function (mat, lvl, t,dens,p.group,con.b,lag.b,p.con) {
  repeat {
    
    repeat{
      v <- ncol(mat)/2                                                                              ###calculate number of variables in matrix
      Phi <- mat[, 1:v]                                                                             ###pull out Phi
      A <- mat[, (v+1):(v*2)]                                                                       ###pull out A
      A_ind <- matrix(0, ncol=v, nrow=v)                                                            ###set A indices to zero
      indices.A <- which(A==0,arr.ind=T)                                                            ###finds indices of zero elements in matrix
      indices.Phi <- which(Phi==0,arr.ind=T)
      indices.A <- indices.A[which(indices.A[,1]!=indices.A[,2]),]                                  ###kick out diagonal (from A matrix only)
      pos.all<-(v*v-v)*2                                                                            ###determine the number of individual paths to add
      cnt.all <- dens*p.group*pos.all  
      rand<-sample(c(round(cnt.all/2),(round(cnt.all)-round(cnt.all/2))),size=2,replace=FALSE)      ###randomly select row/col combinations for individual-level paths
      row.col.A      <- sample(1:nrow(indices.A), rand[1], replace = F) 
      row.col.Phi      <- sample(1:nrow(indices.Phi), rand[2], replace = F) 
      Phi[indices.Phi[row.col.Phi,]] <- lag.b                                                       ###set betas for lagged and contemporaneous paths
      A[indices.A[row.col.A,]]     <- con.b
      noise.inds      <- which(A != 0, arr.ind = TRUE)                                              ####add noise to A betas, SD =.1
      A[noise.inds]   <- A[noise.inds] + rnorm(n = nrow(noise.inds), mean = 0, sd = .01)
      noise.inds      <- which(Phi != 0, arr.ind = TRUE)                                            ###add noise to Phi betas, SD =.1
      Phi[noise.inds] <- Phi[noise.inds] + rnorm(n = nrow(noise.inds), mean = 0, sd = .01)
      break
    }
    
    st <- (t+50)                                                                                    ###This chunk is from Alex's code and is used to simulate the time series 
    noise <- matrix(rnorm(v*st,0,1),v) #
    I     <- diag(v) # identity matrix
    time  <- matrix(0,nrow=v, ncol=(st+1))
    time1 <- matrix(0,nrow=v, ncol=st)
    
    for (i in 1:st){                                                                                ###simulate data points for each time step
      time1[,i]  <- solve(I-A)%*%(Phi%*%time[,i] + noise[,i])
      time[,i+1] <- time1[,i]
    }               
    time1  <- time1[,(51:(50+t))]                                                                   
    series <- t(time1)
    paths  <- cbind(Phi, A)
    if (abs(max(series, na.rm = TRUE)) < 20 & abs(min(series, na.rm = TRUE)) > .01 
        & abs(min(series, na.rm = TRUE)) < 20) break
  }
  
  lvl[is.na(lvl) & paths != 0] <- "ind"
  
  list   <- list("series"  = series,
                 "paths"   = paths,
                 "levels"  = lvl)
  return(list)
}


#######################################################
#### 2. Simulate Data# ################################
#######################################################

# enter simulation parameters
v             <- c(6) # Number of variables
n             <- c(150) # number of individuals
t             <- c(100) # Number of time points
rep           <- seq(1) # replications per condition 
ar            <-c(.6) # ar paths to try
conditions    <- expand.grid(t, n, v, ar,rep)
all           <- rbind(conditions)
colnames(all) <- c("t", "n", "v", "ar", "rep")
negcon <- c(.35, -.35)                                     #to create the negative numbers
neglag <- c(.45, -.45)                                     #to create the negative numbers

# This loop generates a folder name for each condition and replication
for (i in 1:nrow(all)){
  all$folder[i] <- paste("t",all$t[i],
                         "ar",all$ar[i],"rep",all$rep[i],sep="_")
}

rownames(all) <- NULL
colnames(all) <- c("t","n","v","ar","rep","folder")
all$t         <- as.numeric(as.character(all$t))
all$n         <- as.numeric(as.character(all$n))
all$ar         <- as.numeric(as.character(all$ar))


# name directories to place simulated data in (Change)
dir.create('/Users/reneehlozek/Dropbox/CIFAR_Sims/230222_ORIG_26Feb2025/')
data.path <- '/Users/reneehlozek/Dropbox/CIFAR_Sims/230222_ORIG_26Feb2025/data'
true.path <- '/Users/reneehlozek/Dropbox/CIFAR_Sims/230222_ORIG_26Feb2025/true'
level.path <- '/Users/reneehlozek/Dropbox/CIFAR_Sims/230222_ORIG_26Feb2025/levels'

dir.create(data.path)
dir.create(true.path)
dir.create(level.path)

# Creates actual folders for each iteration
folders <- all$folder
for (i in 1:nrow(all)){
  data <- file.path(data.path, folders[i])
  dir.create(data)
  true <- file.path(true.path, folders[i])
  dir.create(true)
  level <- file.path(level.path, folders[i])
  dir.create(level)
}

# Does the simulations
for (i in 1:nrow(all)){
  if (!length(list.files(file.path(data.path,all$folder[i]))) %in% c(50)) {
    
    # generate group matrix for each simulated data set
    res <- mat.generate.asw(p.con = .50, 
                            nvar = all$v[i], AR=all$ar[i],
                            p.group = .50,dens = .20,
                            con.b = sample(negcon, 1), lag.b = sample(neglag, 1))
    
    #for each individual,generate matrix and ts
    for (a in 1:all$n[i]){
      
      out <- ts.generate.asw(mat = res$sub1,
                             lvl = res$lvl1,
                             t   = all$t[i],
                             p.group = .50,dens = .20,
                             con.b = sample(negcon, 1), lag.b = sample(neglag, 1),
                             p.con = .50)
      
      out$series <- round(out$series,digits=5)
      
      
      write.csv(out$series,
                file.path(data.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
      
      write.csv(out$paths,
                file.path(true.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
      
      write.csv(out$levels,
                file.path(level.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
    }
    
  }
}

