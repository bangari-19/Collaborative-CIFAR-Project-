#230222: AR=.6; Contemp = .35; lagged = .45, noise=.01
rm(list=ls())
install.packages('SyncRNG')
library(SyncRNG)
s <- SyncRNG(seed=123456)
for (i in 1:10) {
  cat(s$randi(), '\n')
}

#######################################################
#### 1. load functions needed to simulate the data ####
#######################################################

null_groups <- function(indices,lvl,nrow) {

  result <- sapply(1:nrow(indices), function(i) (indices[i, 1] - 1) * nrow + indices[i, 2])
  c <- which(lvl=="grp")
  any_in_grp <- which(result %in% c)
  print(any_in_grp)
  print('any_in_grp')

  return(indices)

}

convert_to_arr_ind_format <- function(row, col, matrix_dim) {
  arr_ind_result <- which(apply(expand.grid(seq_len(matrix_dim[1]), seq_len(matrix_dim[2])), 1, function(x) all(x == c(row, col))))
  return(arr_ind_result)
}

check_reverse_offdiagonal <- function(mat){
regen <- FALSE
  for (i in 1:nrow(mat)) {
    for (j in 1:i) {
      if ((j!=i) && (mat[i, j] != 0.0) && (mat[j, i] != 0.0)) {
        #cat('i=', i, ', j=', j, ', matij=', mat[i, j], ', matji=', mat[j, i], ', we must regenerate \n')
        regen <- TRUE
        break
      }
    }
  }
  return(regen)

}
check_symmetric_indices <- function(mat) {
  # Get the non-zero indices of the matrix using arr.ind = TRUE
  non_zero_indices <- which(mat != 0, arr.ind = TRUE)
  new_nonz_indices <- non_zero_indices
  # Initialize a vector to store results
  results <- character()
  # Check for each non-zero index if its corresponding off-diagonal symmetric index is in the list
  countnonzero=0
  for (i in seq(1, nrow(non_zero_indices), by = 2)) {
    row_index <- non_zero_indices[i, "row"]
    col_index <- non_zero_indices[i, "col"]

    # Check if the symmetric counterpart is off-diagonal and in the list
    symmetric_index <- c(col_index, row_index)
    if (row_index != col_index && any(apply(non_zero_indices, 1, function(x) all(x == symmetric_index)))) {
      countnonzero <-countnonzero+1
      indexarr <- convert_to_arr_ind_format(col_index,row_index,c(6,12))
      mat[indexarr]<-0
    }
  }

  return(list(mat=mat,nonzeroels=countnonzero))
}
#Group Relations for Set 1: AR: V2, V4, V6; Contemporary Connections: V5 -> V3, V6 -> V4; Lagged Connections: V1lagged -> V6, V4lagged -> V5

mat.generate.asw <- function(p.con,nvar,AR,dens,p.group,con.b,lag.b){
  repeat{
    A <- matrix(0, ncol = nvar, nrow = nvar, )        ###create null contemporaneous matrix
    Phi <- matrix(0, ncol = nvar, nrow = nvar)        ###create null lagged matrix (remember -> ARs in diagonal)
    pos.all <- 2*(nrow(A)*ncol(A))                    ###list all available (non-ar) spots (You could add - 2*nrow(A) to remove AR spots (see AW code))
    cnt.all <- dens*p.group*pos.all                   ###number of group paths given density and proportion of group paths
    indices <- which(Phi == 0, arr.ind = TRUE)        ###create data frame that lists all combinations of columns/rows in matrices. Run "indices <- indices[which(indices[,1] != indices[,2]), ]" to kick out diagonal
    row.col<-c(17, 24, 31, 28)                        ###select group-level associations in lagged
    diag.set<-c(8, 22, 36)                            ###select AR paths 
    n.p.1<-row.col[1:round(p.con*length(row.col))]    ###identify contemporaneous group paths 
    n.p.2<-row.col[(length(n.p.1)+1):length(row.col)] ###identify lagged group paths 
    # print(n.p.1)s
    # print("group cont")
    # print(n.p.2)
    # print("group lagged")
    grp.con <- n.p.1                                  ###name contemporaneous group paths
    grp.lag <- n.p.2                                  ###name lagged group paths
    grp.diag <- diag.set 
    # print(grp.diag)
    # print("AR paths group")                             ###name AR paths
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

# Specify any integer

#This function generates group-level relations, adds individual-level relations to the matrix, 
#adds noise, and simulates the time series (unnecessary for our purposes but helpful for checks!)
ts.generate.asw <- function (mat, lvl, t,dens,p.group,con.b,lag.b,p.con) {
  repeat {
    set.seed(10)
    repeat{
      set.seed(10)
      v <- ncol(mat)/2                                                                             ###calculate number of variables in matrix
      Phi <- mat[, 1:v]                                                                             ###pull out Phi
      A <- mat[, (v+1):(v*2)]                                                                       ###pull out A
      A_ind <- matrix(0, ncol=v, nrow=v)                                                            ###set A indices to zero
      #A[3,5]<-A[5,3]
      group_A_reverse <-check_reverse_offdiagonal(A)
      if(group_A_reverse) {
      print('something is wrong with your group level contem paths')
      print(huh)}
      indices.A <- which(A==0,arr.ind=T) 

      ###finds indices of zero elements in matrix 
      #(ie checks for places where group-level paths not already there)

      # By definition, this already checks that this has no group level paths.
      # Now check if there are reverse paths at the group level

      group_phi_reverse <-check_reverse_offdiagonal(Phi)
      if(group_phi_reverse) {print('something is wrong with your group level lagged paths')}
      
      indices.Phi <- which(Phi==0,arr.ind=T)

      # Actually do the reverse of what we did before. 
      #Find the row, col of the group level paths, and remove them from indices.Phi before simulating indiv level paths

      # now check for reverse group paths and mask them
      indices.A <- indices.A[which(indices.A[,1]!=indices.A[,2]),]                                  ###kick out diagonal (from A matrix only)
      pos.all<-(v*v-v)*2                                                                            ###determine the number of individual paths to add
      cnt.all <- dens*p.group*pos.all  

      

      rand<-sample(c(round(cnt.all/2),(round(cnt.all)-round(cnt.all/2))),size=2,replace=TRUE)   
         ###randomly select row/col combinations for individual-level paths

         # The key here is that we want to remove paths that have group level 
         # now check for the reverse/off-diagonal symmetric ones, and regenerate if found
        # print('group level A')
        # print(A) # contains the group level A
        # print('group level Phi')
        # print(Phi) # contains the group level Phi
        
        row.col.A      <- sample(1:nrow(indices.A), rand[1], replace = F) 
        row.col.Phi      <- sample(1:nrow(indices.Phi), rand[2], replace = F)

        Phi[indices.Phi[row.col.Phi,]] <- lag.b                                                       ###set betas for lagged paths
        A[indices.A[row.col.A,]]     <- con.b                                                         ###set betas for contemporaneous paths
        # print('group + individ level A')
        # print(A) # contains the group level A
        # print('group + individ level Phi')
        # print(Phi) # contains the group level Phi
        # Hardcoding diagonal symmetry to test
        # A[3,5]=A[5,3]
        # Phi[5,4]=10*Phi[4,5]
        # print('before checks A')
        # print(A)
        # print('before checks Phi')
        # print(Phi)
        # print('-----')
        # break
        # break

        regen <-check_reverse_offdiagonal(A)

        while (regen==TRUE){
          A <- mat[, (v+1):(v*2)] # reset A to the group level
          row.col.A      <- sample(1:nrow(indices.A), rand[1], replace = F) 
          A[indices.A[row.col.A,]]     <- con.b
          # print('inside regen A')
          # print(A)
          regen <-check_reverse_offdiagonal(A)
        }
        
        regen <-check_reverse_offdiagonal(Phi)
        while (regen==TRUE){
          Phi <- mat[, 1:v]
          row.col.Phi      <- sample(1:nrow(indices.Phi), rand[2], replace = F)
          Phi[indices.Phi[row.col.Phi,]] <- lag.b 
          # print('inside regen Phi')
          # print(Phi)
          regen <-check_reverse_offdiagonal(Phi)
        }



      # 1. No random individual paths where groups exist (for AR or not)
      # 2. We want to avoid reverse direction paths in individuals/groups (reverse triangle null) 

      nonoisepaths  <- cbind(Phi, A) ### bind Phi and A before adding noise to them to check

      noise.inds      <- which(A != 0, arr.ind = TRUE)                                              ####add noise to A betas, SD =.1
      A[noise.inds]   <- A[noise.inds] + rnorm(n = nrow(noise.inds), mean = 0, sd = .01)
      noise.inds      <- which(Phi != 0, arr.ind = TRUE)                                            ###add noise to Phi betas, SD =.1
      Phi[noise.inds] <- Phi[noise.inds] + rnorm(n = nrow(noise.inds), mean = 0, sd = .01)
      break
    }
    
    st <- (t+100)    #(t+50)                                                            ###This chunk is from Alex's code and is used to simulate the time series 
    noise <- matrix(rnorm(v*st,0,5),v) #
    I     <- diag(v) # identity matrix
    time  <- matrix(0,nrow=v, ncol=(st+1))
    time1 <- matrix(0,nrow=v, ncol=st)
    
    for (i in 1:st){                                                                                ###simulate data points for each time step
      time1[,i]  <- solve(I-A)%*%(Phi%*%time[,i] + noise[,i])
      time[,i+1] <- time1[,i]
    }               
    time1  <- time1[,(101:(100+t))] #time1[,(51:(50+t))]                                                                   
    series <- t(time1)
    paths  <- cbind(Phi, A)
    if (abs(max(series, na.rm = TRUE)) < 20 & abs(min(series, na.rm = TRUE)) > .01 
        & abs(min(series, na.rm = TRUE)) < 20) break
  }
  
  lvl[is.na(lvl) & paths != 0] <- "ind"
  
  list   <- list("series"  = series,
                 "paths"   = paths,
                 "nonoisepaths" = nonoisepaths,
                 "levels"  = lvl)
  return(list)
}


#######################################################
#### 2. Simulate Data# ################################
#######################################################

# enter simulation parameters
v             <- c(6) # Number of variables
n             <- c(100) # number of individuals
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
dir.create('/Users/reneehlozek/Dropbox/CIFAR_Sims/240415')

data.path <- '/Users/reneehlozek/Dropbox/CIFAR_Sims/240415/data'
true.path <- '/Users/reneehlozek/Dropbox/CIFAR_Sims/240415/true'
nonoisetrue.path <- '/Users/reneehlozek/Dropbox/CIFAR_Sims/240415/true_nonoise'
level.path <- '/Users/reneehlozek/Dropbox/CIFAR_Sims/240415/levels'

dir.create(data.path)
dir.create(true.path)
dir.create(nonoisetrue.path)
dir.create(level.path)

# Creates actual folders for each iteration
folders <- all$folder
for (i in 1:nrow(all)){
  data <- file.path(data.path, folders[i])
  dir.create(data)
  true <- file.path(true.path, folders[i])
  dir.create(true)
  nonoisetrue <- file.path(nonoisetrue.path, folders[i])
  dir.create(nonoisetrue)
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

      write.csv(out$nonoisepaths,
                file.path(nonoisetrue.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
      
      write.csv(out$levels,
                file.path(level.path, all$folder[i], paste0("ind_", a, ".csv")),
                row.names = FALSE)
    } # ends for loop
  
  } # ends if length loop
} # ends simulation for loop