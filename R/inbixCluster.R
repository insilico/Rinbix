# ----------------------------------------------------------------------------
# inbixCluster.R - Bill White - 10/10/15
#
# Rinbix package cluster/modularity functions.

# ----------------------------------------------------------------------------
#' Load ripM derived modules from RData file with filename.
#' 
#' \code{loadRipmResultsFromRdata} 
#' 
#' @param rdata_filename String filename of saved module object.
#' @return ripM module object.
#' @export
loadRipmResultsFromRdata <- function(rdata_filename) {
  # save modules to separate files
  load(file=rdata_filename)
  TRUE  
}

# -----------------------------------------------------------------------------
#' Attempt merge by sum of powers of passed matrix.
#' 
#' \code{mergeSumPowers} 
#' 
#' @param Aadj Matrix adjacency matrix.
#' @param startOrder Numeric starting order for merging.
#' @param maxOrder Numeric maximum order for merging.
#' @param minModuleSize Numeric minimum allowed size of module.
#' @param maxModuleSize Numeric maximum allowed size of module.
#' @param verbose Flag send verbose messages to standard out.
#' @param indentLevel Numeric tab indent level (recursion level).
#' @return list with list of modules (lists of variable names) and number of iterations.
mergeSumPowers <- function(Aadj, startOrder=2, maxOrder=4, 
                           minModuleSize=30, maxModuleSize=200, 
                           verbose=FALSE, indentLevel=1) {
  indent <- paste(rep("  ", indentLevel), collapse="")
  if(verbose) cat(indent, "[MERGE FUNCTION] Optimization function starting\n")
  n <- ncol(Aadj)
  this_order <- startOrder
  opt_mod_list_all <- NULL
  opt_mod_list_all_q <- NULL
  best_mod_idx <- 1
  best_mod_list <- NULL
  best_mod_list_len <- n
  best_mod_list_len_idx <- 1
  best_mod_order <- startOrder
  best_mod_q <- 0
  merge_success <- FALSE
  process_best_mod_list <- TRUE
  keep_looping <- TRUE
  iteration <- 0
  all_q_zero <- TRUE
  while(keep_looping & (this_order <= maxOrder)) {
    iteration <- iteration + 1
    #indent <- paste(rep("  ", indentLevel), iteration, collapse="")
    if(verbose) cat(indent, "[MERGE FUNCTION] Calling sum of matrix powers order:", this_order, "\n")
    this_matrix <- sum_powers(Aadj, this_order)
    opt_mod_list <- modularity(this_matrix)
    opt_mod_list_all <- c(opt_mod_list_all, list(opt_mod_list))
    opt_mod_list_all_q <- c(opt_mod_list_all_q, opt_mod_list$q)
    opt_mod_list_q <- round(opt_mod_list$q, 6)
    if(verbose) cat(indent, "[MERGE FUNCTION] Q =", opt_mod_list_q, "\n")
    if(opt_mod_list_q == 0) { 
      if(verbose) cat(indent, "[MERGE FUNCTION] Q == 0 detected, cannot split module, trying next power\n")
    } else {
      all_q_zero <- FALSE
      opt_mod_list_len <- length(opt_mod_list$modules)
      if(verbose) cat(indent, "[MERGE FUNCTION] Number of modules:", opt_mod_list_len, "\n")
      if(opt_mod_list_len < best_mod_list_len) {
        best_mod_idx <- iteration
        best_mod_list_len <- opt_mod_list_len
      }
      opt_mod_sizes <- sapply(1:length(opt_mod_list$modules), 
                              function(x) length(opt_mod_list$modules[[x]]))
      if(verbose) cat(indent, "[MERGE FUNCTION] Module sizes:", opt_mod_sizes, "\n")
      opt_big_mod_mask <- opt_mod_sizes >= maxModuleSize  # for split.
      opt_small_mod_mask <- opt_mod_sizes <= minModuleSize  # for merge.
      if(verbose) cat(indent, "[MERGE FUNCTION] Modules >= ", maxModuleSize, "nodes: ", opt_big_mod_mask, "\n")
      if(verbose) cat(indent, "[MERGE FUNCTION] Modules <= ", minModuleSize, "nodes: ", opt_small_mod_mask, "\n")
      # if we still have small modules, go to next power and try again
      if(all(opt_small_mod_mask == FALSE)) {
        # found a good merge!
        merge_success <- TRUE
        keep_looping <- FALSE
        best_mod_idx <- iteration
        best_mod_list <- opt_mod_list
      }
    }
    this_order <- this_order + 1
  }

  # small module optimization is done, determine what to return to caller
  if(merge_success) {
    # nice parition found with no small modules
    if(verbose) cat(indent, "[MERGE FUNCTION] Found a power sum with no small modules, order:", best_mod_order, "\n")
      best_mod_list <- opt_mod_list_all[[best_mod_idx]]
  } else {
    # merge was not perfect, but if not all Q=0, pick the smallest module
    if(all_q_zero) {
      best_mod_list <- list(modules=list(1:n))
      best_mod_q <- 0
      best_mod_order <- 0
    } else {
      if(verbose) cat(indent, "[MERGE FUNCTION] Did not find a power sum with no small modules\n")
      if(verbose) cat(indent, "[MERGE FUNCTION] Picking minimum split size, index:", best_mod_idx, "\n")
      best_mod_list <- opt_mod_list_all[[best_mod_idx]]
      best_mod_q <- opt_mod_list_all_q[[best_mod_idx]]
      merge_success <- TRUE
    }
  }

  if(verbose) cat(indent, "[MERGE FUNCTION] Optimization function done\n")
  list(success=merge_success, modules.obj=best_mod_list, q=best_mod_q, order=best_mod_order)
}

# -----------------------------------------------------------------------------
#' Detect modular structure in a network.
#' 
#' \code{modularity}
#' 
#' @param G Matrix Adjacency matrix.
#' @param verbose Flag to output messages to stdout.
#' @return list of groups, Q and module assignments.
#' @export
modularity <- function(G, verbose=FALSE) {
  MODULARITY_THRESHOLD <- 0.000001
  n <- nrow(G);
  geneNames <- colnames(G)
  
  # create adjacency matrix by thresholding and/or conversion to binary
  # zero the diagonal
  diag(G) <- 0
  # create real symmetric modularity matrix B
  k <- colSums(G)
  m <- 0.5 * sum(k)
  B <- G - k %*% t(k) / (2.0 * m);
  
  # column indices
  firstModule <- seq(from=1, to=n)
  # stack for recursive subdivision of modules
  processStack <- list(firstModule)
  # list of vectors of column indices
  modules <- NULL
  # Q is a measure of "goodness" of the module decomposition
  Q <- 0
  
  # "recursive" loop to find submodules
  iteration <- 0
  while(length(processStack) > 0) {
    iteration <- iteration + 1
    # pop the stack for a list of column indices
    thisModule <- unlist(processStack[[length(processStack)]])
    #cat("this module:", thisModule, "\n")
    processStack[length(processStack)] <- NULL
    # create Bg, a submatrix of B based on current module indices
    newDim <- length(thisModule)
    Bg <- matrix(ncol=newDim, nrow=newDim, data=c(0))
    for(l1 in 1:newDim) {
      for(l2 in 1:newDim) {
        Bg[l1, l2] = B[thisModule[l1], thisModule[l2]];
      }
    }
    # adjust the diagonal
    rowsums <-rowSums(Bg)
    for(i in 1:length(rowsums)) {
      Bg[i, i] = Bg[i, i] - rowsums[i]
    }
    # get the best split of the modules based on eigenvector decomposition
    sub_modules <- modularityBestSplit(Bg, m)
    deltaQ <- sub_modules$Q
    s <- sub_modules$s_out
    # assign indices based on two groups
    s1 <- thisModule[s==-1]
    s2 <- thisModule[s==1]
    if((length(s1) == 0) || (length(s2) == 0)) {
      # stopping criteria for recursive splits
      modules[[length(modules)+1]] <- thisModule
      if(iteration == 1) {
        Q <- deltaQ
      }
    }
    else {
      if(deltaQ <= MODULARITY_THRESHOLD) {
        # stopping criteria for recursive splits
        modules[[length(modules)+1]] <- thisModule
      } else {
        # "recursive step": push the two groups onto the stack and repeat
        processStack[[length(processStack)+1]] <- s1
        processStack[[length(processStack)+1]] <- s2
        # update cummulative Q
        Q <- Q + deltaQ
      }
    }
  }
  
  # ----------------------------------------------------------------------------
  # output modules
  # ----------------------------------------------------------------------------
  if(verbose) {
    cat("Q:", Q, "\n")
    cat("Number of modules:", length(modules), "\n")
  }
  groupAssignments <- NULL
  for(i in 1:length(modules)) {
    modIdx <- modules[[i]]
    modNum <- length(modIdx)
    modGenes <- geneNames[modIdx]
    #cat("Module", i, ":", modGenes, "\n")
    modGroup <- rep(i, modNum)
    thisGroup <- cbind(modGenes, modGroup)
    groupAssignments <- rbind(groupAssignments, thisGroup)
  }
  colnames(groupAssignments) <- c("Gene", "Group")
  
  list(groups=groupAssignments, q=Q, modules=modules)
}

# -----------------------------------------------------------------------------
#' Find the best split of a modularity matrix using eigenvector decomposition.
#' 
#' \code{modularityBestSplit} 
#' 
#' @param B Modularity matrix.
#' @param m m (average degree?)
#' @return List with modularity value Q and best split vector.
#' @export
modularityBestSplit <- function(B, m) {
  # function to split columns of matrix into two groups
  # get the maximum eigenvector
  eigResult <- eigen(B);
  eigval <- eigResult$values;
  eigvec <- eigResult$vectors
  # in R, this is the first vector
  maxeig_val <- eigval[1];
  maxeig_vec <- eigvec[,1];
  # use the sign of the eigenvector values to assign group status +/-1
  s_out <- ifelse(maxeig_vec < 0, -1, 1)
  # calculate Q for this split
  Q_mat <- t(s_out) %*% B %*% s_out
  Q <- Q_mat[1,1]
  Q <- Q * (1.0 / (m * 4.0))
  
  # return Q and list assignments
  list(Q=Q, s_out=s_out)
}

# ----------------------------------------------------------------------------
#' Run Recursive Indirect Paths Modularity (rip-M).
#' 
#' \code{ripM} 
#' 
#' @param Acorr Matrix correlation matrix.
#' @param thresholdType String threshold type: "hard" or "soft".
#' @param thresholdValue Numeric hard threshold correlation value or soft threshold power.
#' @param startMergeOrder Numeric power n to raise adjacencyMatrix^n in first merge attempt.
#' @param maxMergeOrder Numeric power n to raise adjacencyMatrix^n in final merge attempt.
#' @param maxModuleSize Numeric maximum allowed size of module.
#' @param minModuleSize Numeric minimum allowed size of module.
#' @param useAbs Flag take absolute value of the correlation matrix.
#' @param useWeighted Flag use weighted adjacency matrix versus binary.
#' @param hubSelection String how to choose the hub for each module (weighted, posthresh).
#' @param simpleModularity Flag perform simple modularity without split/merge.
#' @param verbose Flag send verbose messages to standard out .
#' @return list with: list of modules (list of lists of variable names), list of hubs, 
#' list of node degrees.
#' @export
ripM <- function(Acorr, thresholdType="hard", thresholdValue=0.8, 
                 startMergeOrder=2, maxMergeOrder=4,
                 maxModuleSize=200, minModuleSize=30, 
                 useAbs=TRUE, useWeighted=FALSE, 
                 hubSelection="weighted", simpleModularity=FALSE, 
                 verbose=FALSE) {
  # ---------------------------------------------------------------------------
  # preprocessing
  Aorig <- Acorr
  if(useAbs) {
    if(verbose) cat("Taking absolute value of Acorr\n")
    Acorr <- abs(Acorr)
  }
  # thresholding
  if(thresholdType == "soft") {
    if(verbose) cat("Soft threshold Acorr^", thresholdValue, "\n", sep="")
    Acorr <- Acorr ^ thresholdValue
  } else {
    if(verbose) cat("Hard threshold abs(Acorr) >", thresholdValue, "\n")
    passThreshold <- Acorr > thresholdValue
    # create adjacency matrix from thresholding
    if(useWeighted) {
      if(verbose) cat("Keeping weighted values that pass threshold [", thresholdValue, "]\n")
      Acorr[!passThreshold] <- 0
    } else {
      if(verbose) cat("Keeping values that pass threshold [", thresholdValue, "] as binary 0/1\n")
      Acorr[passThreshold] <- 1
      Acorr[!passThreshold] <- 0
    }
  }

  # ---------------------------------------------------------------------------
  # ripM
  rip_modules <- list()
  if(simpleModularity) {
    # call plain old modularity
    if(verbose) cat("ripM SIMPLE MODULARITY algorithm starting...\n")
    simple_mod_list <- modularity(Acorr)
    mod_list_vars <- list()
    for(mod_idx in 1:length(simple_mod_list$modules)) {
      this_module <- unlist(simple_mod_list$modules[[mod_idx]])
      #cat("This module length:", length(this_module), "\n")
      this_var_names <- colnames(Acorr)[this_module]
      mod_list_vars <- c(mod_list_vars, list(this_var_names))
      #mod_list_vars <- c(mod_list_vars, list(this_module))
    }
    rip_modules$module_list <- mod_list_vars
    rip_modules$iterations <- 1
  } else {
    # call stack-based kernel
    if(verbose) cat("ripM MERGE-AND-SPLIT algorithm starting...\n")
    rip_modules <- ripMKernelStack(Acorr, startOrder=startMergeOrder, maxOrder=maxMergeOrder, 
                                   minModuleSize=minModuleSize, maxModuleSize=maxModuleSize, 
                                   verbose, indentLevel=1)
  }
  
  # ---------------------------------------------------------------------------
  # post-processing return values
  # get the highest degree node (hub) for each module
  # determine node degrees
  if(verbose) cat("ripM preparing return values: hubs, degrees, sizes. adjacency...\n")
  diag(Acorr) <- 0
  global_degrees <- rowSums(Acorr)
  names(global_degrees) <- colnames(Acorr)
  module_hubs <- NULL
  module_degrees <- NULL
  for(module_idx in 1:length(rip_modules$module_list)) {
    this_module <- rip_modules$module_list[[module_idx]]
#     this_global_degrees <- global_degrees[this_module]
#     names(this_global_degrees) <- this_module
    if(length(this_module) > 1) {
      this_module_matrix <- Acorr[this_module, this_module]
      this_module_matrix_degrees <- rowSums(this_module_matrix)
    } else {
      this_module_matrix <- Acorr[this_module, ]
      this_module_matrix_degrees <- c(sum(this_module_matrix))
    }
    names(this_module_matrix_degrees) <- this_module
    this_module_hub <- this_module_matrix_degrees[which.max(this_module_matrix_degrees)]
    this_module_hub_degree <- this_module_matrix_degrees[names(this_module_hub)]
    module_hubs <- c(module_hubs, this_module_hub)
    module_degrees <- c(module_degrees, this_module_hub_degree)
  }
  rip_modules$hubs <- module_hubs
  rip_modules$sizes <- sapply(rip_modules$module_list, FUN=length)
  rip_modules$degrees <- module_degrees
  rip_modules$adj <- Acorr
  
  # return
  rip_modules
}

# -----------------------------------------------------------------------------
#' Run Recursive Indirect Paths Modularity (ripM) stack-based kernel.
#' 
#' \code{ripMKernelStack} 
#' 
#' @param Aadj Matrix adjacency matrix.
#' @param startOrder Numeric starting order for merging.
#' @param maxOrder Numeric maximum order for merging.
#' @param minModuleSize Numeric minimum allowed size of module.
#' @param maxModuleSize Numeric maximum allowed size of module.
#' @param verbose Flag send verbose messages to standard out .
#' @param indentLevel Flag tab indent level (recursion level).
#' @return list with list of modules (lists of variable names) and number of iterations.
ripMKernelStack <- function(Aadj, startOrder=2, maxOrder=4, 
                            minModuleSize=30, maxModuleSize=200, 
                            verbose=FALSE, indentLevel=1) {
  if(ncol(Aadj) == 0) { return(list()) }
  var_names <- colnames(Aadj)
  module_list <- NULL
  indent <- paste(rep("  ", indentLevel), collapse="")
  if(verbose) cat(indent, "----------------------------------------------------------\n")
  if(verbose) cat(indent, "ripM STACK algorithm starting...\n")
  iteration <- 0
  max_iteration <- 100
  first_module <- seq(from=1, to=ncol(Aadj))
  process_stack <- list(first_module)
  module_list <- list()
  prev_module <- rep(-1, ncol(Aadj))
  while((length(process_stack) > 0) && (iteration < max_iteration)) {
    iteration <- iteration + 1
    if(verbose) cat(indent, "*****************************************************\n")
    if(verbose) cat(indent, "Iteration:", iteration, "\n")

    # pop the stack for a list of column indices (a module)
    if(verbose) cat(indent, "[POP] Preparing module matrix from module indices\n")
    this_module <- unlist(process_stack[[length(process_stack)]])
    if(verbose) cat(indent, "This module length:", length(this_module), "\n")
    process_stack[length(process_stack)] <- NULL
    sub_matrix <- Aadj[this_module, this_module]
    this_var_names <- colnames(sub_matrix)
    # keep from trying to process the same module over and over
    if((length(this_module) == length(prev_module)) && all(this_module == prev_module)) {
      cat("PREVIOUSLY SEEN ON THE STACK - SAVING!\n")
      this_mod_vars <- this_var_names[prev_module]
      module_list[[length(module_list)+1]] <- this_mod_vars
      prev_module <- NULL
      next
    }
    
    # run modularity on the first module on the stack  
    if(verbose) cat(indent, "Running modularity on this module matrix\n")
    mod_list <- modularity(sub_matrix)
    this_q <- round(mod_list$q, 6)
    if(verbose) cat(indent, "Q =", this_q, "\n")
    if((length(mod_list$modules) == 1) || (this_q == 0)) {
      if(verbose) cat(indent, "Zero Q detected. Cannot split module.\n")
      this_mod_idx <- mod_list$modules[[1]]
      this_mod_len <- length(this_mod_idx)
      this_mod_vars <- this_var_names[this_mod_idx]
      if(verbose) cat(indent, "Saving irreducible module, size:", this_mod_len,"\n")
      module_list[[length(module_list)+1]] <- this_mod_vars
      next
    }
    if(verbose) cat(indent, "Number of modules:", length(mod_list$modules), "\n")
    mod_sizes <- sapply(1:length(mod_list$modules), function(x) length(mod_list$modules[[x]]))
    if(verbose) cat(indent, "Module sizes:", mod_sizes, "\n")
    
    # first-level constraints check
    big_mod_mask <- mod_sizes >= maxModuleSize  # for split.
    small_mod_mask <- mod_sizes < maxModuleSize  # for merge.
    if(verbose) cat(indent, "Modules >= ", maxModuleSize, "nodes: ", big_mod_mask, "\n")
    if(verbose) cat(indent, "Modules < ", maxModuleSize, "nodes: ", small_mod_mask, "\n")
    
    # look at each module and split/merge
    small_mods_list <- list()
    small_mods_list$modules <- mod_list$modules[small_mod_mask]
    small_mod_idx <- NULL
    for(i in 1:length(mod_sizes)) {
      this_mod <- mod_list$modules[[i]]
      this_mod_length <- length(this_mod)
      this_mod_vars <- this_var_names[this_mod]
      if(big_mod_mask[i]) {
        # call split on the subnetwork of nodes from this module
        # !!!NOTE: MAP SUBMATRIX INDICES BACK TO ORIGINAL ADJ MATRIX INDICES!!!
        if(verbose) cat(indent, "[PUSH] Pushing large module onto recursion stack, size:", this_mod_length, "\n")
        orig_adj_idx <- which(var_names %in% this_mod_vars)
        process_stack[[length(process_stack)+1]] <- orig_adj_idx
      } else {
        if(small_mod_mask[i]) {
          # add to merge list indices
          if(verbose) cat(indent, "Collecting small module indices, size:", this_mod_length, "\n")
          small_mod_idx <- c(small_mod_idx, this_mod)
        }
      }
    }
    
    # see if we can make a better module split with higher powers of remaining modules  
    if(length(small_mods_list$modules) > 1) {
      if(verbose) cat(indent, "[MERGE] Attempting SMALL module optimization with higher-power matrix sums\n")
      opt_mod_matrix <- sub_matrix[small_mod_idx, small_mod_idx]
      if(verbose) cat(indent, "[MERGE] Collected small modules into matrix size:", 
          dim(opt_mod_matrix)[1], "\n")
      opt_var_names <- this_var_names[small_mod_idx]
      colnames(opt_mod_matrix) <- opt_var_names
      diag(opt_mod_matrix) <- 0
      
      # attempt merging small modules coalesced into a matrix
      merge_results <- mergeSumPowers(opt_mod_matrix, 
                                      startOrder=startOrder, 
                                      maxOrder=maxOrder, 
                                      minModuleSize=minModuleSize, 
                                      maxModuleSize=maxModuleSize, 
                                      verbose=verbose,
                                      indentLevel=2)
      if(merge_results$success) {
        final_mod_list <- merge_results$modules.obj
        if(verbose) cat(indent, "[MERGE] Success!\n")
        if(verbose) cat(indent, "[MERGE] Processing best module list\n")
        final_mod_sizes <- unlist(lapply(final_mod_list$modules, length))
        # recurse any big modules from the small module optimization
        # optimization-level constraints check
        final_big_mod_mask <- final_mod_sizes >= maxModuleSize  # for split.
        if(verbose) cat(indent, "[MERGE] modules >= ", maxModuleSize, "nodes: ", final_big_mod_mask, "\n")
        # save small modules if any got through
        final_small_mod_mask <- final_mod_sizes <= minModuleSize  
        if(verbose) cat(indent, "[MERGE] modules <= ", minModuleSize, "nodes: ", final_small_mod_mask, "\n")
        # save Goldilocks modules
        final_good_mod_mask <- (final_mod_sizes > minModuleSize) & (final_mod_sizes < maxModuleSize)
        if(verbose) cat(indent, "[MERGE] Goldilocks modules: ", final_good_mod_mask, "\n")
        # process the small module optimization results
        for(i in 1:length(final_mod_sizes)) {
          this_mod_idx <- final_mod_list$modules[[i]]
          this_mod_len <- length(this_mod_idx)
          this_mod_vars <- opt_var_names[this_mod_idx]
          if(any(is.na(this_mod_vars))) {
            browser()
            stop(indent, "NAs in module node names")
          }
          if(final_big_mod_mask[i]) {
            if(verbose) cat(indent, "[MERGE] [PUSH] large module onto recursion stack, size:", this_mod_len, "\n")
            # !!!NOTE: MAP MERGE INDICES BACK TO ORIGINAL ADJ MATRIX INDICES!!!
            orig_adj_idx <- which(var_names %in% this_mod_vars)
            #orig_adj_idx <- which(this_mod_vars %in% var_names)
            process_stack[[length(process_stack)+1]] <- orig_adj_idx
          } else {
            if(final_small_mod_mask[i]) {
              if(verbose) cat(indent, "[MERGE] Saving unredeemable small module, size:", this_mod_len, "\n")
              module_list[[length(module_list)+1]] <- this_mod_vars
            } else {
              # this module passed constraints, save it
              if(verbose) cat(indent, "[MERGE] Saving Goldilocks module, size:", this_mod_len,"\n")
              module_list[[length(module_list)+1]] <- this_mod_vars
            }
          }
        } # END: for all final merge modules
      } else {
        if(verbose) cat(indent, "[MERGE] Merge failed, saving small module list as single module\n")
        module_list[[length(module_list)+1]] <- opt_var_names
      } # END: if merge success
      # END: small module merge optimization
      if(verbose) cat(indent, "[MERGE] Ending SMALL module optimization successfully\n")
    } else {
      if(length(small_mods_list$modules) == 1) {
        this_mod_idx <- small_mods_list$modules[[1]]
        this_mod_len <- length(this_mod_idx)
        this_mod_vars <- this_var_names[this_mod_idx]
        if(verbose) cat(indent, "[MERGE] Saving orphan module, size:", this_mod_len,"\n")
        module_list[[length(module_list)+1]] <- this_mod_vars
        if(verbose) 
          cat(indent, "[MERGE] Ending SMALL module optimization successfully - single module\n")
      }
    }
    # keep from looping infinitely over same module
    prev_module <- this_module
  } # END: while stack processing
  if(verbose) cat(indent, "[DONE]\n")
  if(iteration == max_iteration) {
    stop(paste(indent, "FATAL ERROR: Max iteration", max_iteration, "reached"))
  }  
  if(verbose) cat(indent, "ripM exiting normally\n")
  
  list(module_list=module_list, iterations=iteration)
}

# -----------------------------------------------------------------------------
#' Assign all variables names in a rip-M result object module list index values.
#' 
#' \code{ripmVariableModuleAssignments} 
#' 
#' @param ripm_result Rip-M result object.
#' @return data frame with two columns: var(iable) and module (integer).
#' @export
ripmVariableModuleAssignments <- function(ripm_result) {
  variable_assignments <- NULL
  for(i in 1:length(ripm_result$module_list)) {
    this_mod_list <- ripm_result$module_list[[i]]
    for(j in 1:length(this_mod_list)) {
      #cat(i, j, this_mod_list[j], "\n")
      variable_assignments <- rbind(variable_assignments,
                                    data.frame(var=this_mod_list[j], cluster=i, stringsAsFactors=F))
    }
  }
  variable_assignments <- variable_assignments[order(variable_assignments$var), ]
  variable_assignments
}

# ----------------------------------------------------------------------------
#' Save ripM derived modules to GMT format file.
#' http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/
#' Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29
#' 
#' \code{saveModuleListGmt} 
#' 
#' @param module_obj Module object from ripM call.
#' @param output_prefix String filename prefix for output GMT file.
#' @export
saveModuleListGmt <- function(module_obj, output_prefix) {
  # save modules to gmt
  cat("Saving", length(module_obj$module_list), "modules to GMT and transposed GMT.\n")
  results <- NULL
  for(mod_idx in 1:length(module_obj$module_list)) {
    this_mod <- module_obj$module_list[[mod_idx]]
    results <- rbind(results, data.frame(paste(paste("hub_", module_obj$hubs[mod_idx], sep=""), 
                                               "\t", paste(this_mod, collapse="\t"), sep="")))
  }
  gmt_filename <- paste(output_prefix, ".gmt", sep="")
  cat("Saving GMT:", gmt_filename, "\n")
  #lapply(results, cat, "\n", file=gmt_filename, append=FALSE)
  write.table(results, gmt_filename, quote=F, row.names=F, col.names=F)
  # transpose for spreadsheets
  #gmt_filename <- paste(output_prefix, ".gmtt", sep="")
  #cat("Saving transposed GMT:", mod_filename, "\n")
  #write.table(t(results), gmt_filename, quote=F, row.names=F, col.names=F)

  TRUE  
}

# ----------------------------------------------------------------------------
#' Save ripM derived modules to RData file with filename prefix.
#' 
#' \code{saveRipmResultsToRdata} 
#' 
#' @param Module_obj module object from ripM call.
#' @param output_prefix String filename prefix for module output files.
#' @export
saveRipmResultsToRdata <- function(results_list, output_prefix) {
  # save modules to separate files
  save_filename <- paste(output_prefix, "_results.Rdata", sep="")
  cat("Saving", save_filename, "\n")
  save(results_list, file=save_filename)
  TRUE  
}

# ----------------------------------------------------------------------------
#' Save ripM derived modules to separate files with filename prefix.
#' 
#' \code{saveRipmModules} 
#' 
#' @param module_obj Module object from ripM call.
#' @param output_prefix String filename prefix for module output files.
#' @export
saveRipmModules <- function(module_obj, output_prefix) {
  # save modules to separate files
  cat("Saving", length(module_obj$module_list), "modules to separate files.\n")
  results <- NULL
  for(mod_idx in 1:length(module_obj$module_list)) {
    this_mod <- module_obj$module_list[[mod_idx]]
    mod_filename <- paste(output_prefix, "_module", mod_idx, ".tab", sep="")
    cat("Saving", mod_filename, "\n")
    write.table(this_mod, mod_filename, quote=F, row.names=F, col.names=F)
  }
  TRUE  
}

# -----------------------------------------------------------------------------
#' Assign all variables names in a WGCNA result object cluster list index values.
#' 
#' \code{wgcnaVariableModuleAssignments} 
#' 
#' @param wgcna_result WGCNA result object.
#' @return data frame with two columns: var(iable) and module (integer).
#' @export
wgcnaVariableModuleAssignments <- function(wgcna_result) {
  variable_assignments <- NULL
  for(i in 1:length(wgcna_result$mods.genes)) {
    this_mod_list <- wgcna_result$mods.genes[[i]]
    for(j in 1:length(this_mod_list)) {
      #cat(i, j, this_mod_list[j], "\n")
      variable_assignments <- rbind(variable_assignments,
                                    data.frame(var=this_mod_list[j], cluster=i, stringsAsFactors=F))
    }
  }
  variable_assignments <- variable_assignments[order(variable_assignments$var), ]
  variable_assignments
}

