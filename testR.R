convert_to_arr_ind_format <- function(row, col, matrix_dim) {
  arr_ind_result <- which(apply(expand.grid(seq_len(matrix_dim[1]), seq_len(matrix_dim[2])), 1, function(x) all(x == c(row, col))))
  return(arr_ind_result)
}

check_reverse_offdiagonal <- function(mat){
regen <- FALSE
# print('checking mat')
# print(mat)
# print(mat[1,0])
  for (i in 1:nrow(mat)) {
    for (j in 1:i) {
      # cat('i,j ',i,j,'\n')
      # cat('matij ',mat[i,j],'\n')
      # cat('matji ', mat[j,i],'\n')
      if ((j!=i) && (mat[i, j] != 0.0) && (mat[j, i] != 0.0)) {
        cat('i=', i, ', j=', j, ', matij=', mat[i, j], ', matji=', mat[j, i], ', we must regenerate \n')
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
#   print(dim(mat)[2])
#   print('dimension of matrix')
#   print(mat)
#   print('mat')
#   print(new_nonz_indices)
#   print('as nonz')
  # Check for each non-zero index if its corresponding off-diagonal symmetric index is in the list
  for (i in seq(1, nrow(non_zero_indices), by = 2)) {
    row_index <- non_zero_indices[i, "row"]
    col_index <- non_zero_indices[i, "col"]

    # Check if the symmetric counterpart is off-diagonal and in the list
    symmetric_index <- c(col_index, row_index)
    if (row_index != col_index && any(apply(non_zero_indices, 1, function(x) all(x == symmetric_index)))) {
#     result <- paste("Index (", non_zero_indices[i, "row"], ",", non_zero_indices[i, "col"], ") has a symmetric 
 #     off-diagonal counterpart in the list.\n")
    #   new_nonz_indices <-new_nonz_indices[setdiff(names(new_nonz_indices)), non_zero_indices[i, "row"]]
    #   print(new_nonz_indices)
      indexarr <- convert_to_arr_ind_format(row_index,col_index,c(3,3))
      mat[indexarr]<-0

  #    results <- c(results, result)
    }
  }


  return(mat)
}

# # Example usage:
# mat <- matrix(c(1, 2, 0, 2, 0, 3, 0, 3, 4), nrow = 3, byrow = TRUE)
# print(mat)
# print('mat before')

# newmat <- check_symmetric_indices(mat)
# # cat("Results:\n")
# # cat(output, "\n")
# print(newmat)
# print('mat after')


# # Example usage:
# row <- 2
# col <- 1
# matrix_dim <- c(6, 12)

# arr_ind_result <- convert_to_arr_ind_format(row, col, matrix_dim)
# print('row')
# print(row)
# print('col')
# print(col)
# print('this is the index')
# print(arr_ind_result)


number<-3
mat <- diag(4)
row <- sample(1:4, number, replace = TRUE)
col <- sample(1:4, number, replace = TRUE)
for (i in 1:number) {
  mat[row[i], col[i]] <- runif(1)
}
mat[1,2]=2
mat[2,1]=1

print(mat)
regen<-FALSE
print('regen')
print(regen)
regen<-check_reverse_offdiagonal(mat)
print('regen')
print(regen)
