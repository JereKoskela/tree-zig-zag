##############################################################################################
# Enter the relative location of the tree-zig-zag/Results/ folder from your working directory.
# Run getwd() to find your working directory.
working_dir <- ".tree-zig-zag/Results/"
##############################################################################################

library("mcmcse")
library("tictoc")
tic()
mutation_models <- c("finitesites", "infinitesites")
scenarios <- c("default", "large", "long")
# velocities used in the experiments, in the same order as the loops over mutation models and scenarios
v_theta_zigzags <- c(4, 4, 20, 8, 6, 40)
v_theta_hybrids <- c(4, 4, 20, 8, 6, 40)
# dummy variable for looping through velocity vectors
k <- 1
# number of batches from which to compute ESS estimates for the zigzag and hybrid models
batches <- 100

for (mut in mutation_models) {
  for (sce in scenarios) {
    zigzag <- as.matrix(read.table(paste0(working_dir, "zigzag-", mut, "-", sce, ".txt")))
    hybrid <- as.matrix(read.table(paste0(working_dir, "hybrid-", mut, "-", sce, ".txt")))
    metro <- as.matrix(read.table(paste0(working_dir, "metropolis-", mut, "-", sce, ".txt")))

    n_zz <- dim(zigzag)[1]
    n_hy <- dim(hybrid)[1]
    ncol <- dim(zigzag)[2]
    v_theta_zigzag <- v_theta_zigzags[k]
    v_theta_hybrid <- v_theta_hybrids[k]

    times_zigzag <- c(0,cumsum(abs(diff(zigzag[,1])))) / v_theta_zigzag
    times_hybrid <- c(0,abs(diff(hybrid[,1]))) / v_theta_hybrid
    times_hybrid[(1:n_hy)[hybrid[,ncol+1]==1]] <- 0
    times_hybrid <- cumsum(times_hybrid)

    par_vec <- c(5,5,4,2)
    cex_var <- 2

    ymin <- floor(min(metro[,1], hybrid[,1], zigzag[,1]))
    ymax <- ceiling(max(metro[,1], hybrid[,1], zigzag[,1]))
    png(paste0("zigzag-", mut, "-", sce, "-mutationrate.png"))
    par(mar=par_vec)
    plot(times_zigzag, zigzag[,1], type="l",ylim=c(ymin,ymax),ylab="Mutation rate",xlab="Time",main="Zig-zag",las=1,cex.lab=cex_var,cex.axis=cex_var,cex.main=cex_var)
    dev.off()
    png(paste0("hybrid-", mut, "-", sce, "-mutationrate.png"))
    par(mar=par_vec)
    plot(times_hybrid, hybrid[,1], type="l",ylim=c(ymin,ymax),ylab="Mutation rate",xlab="Time",main="Hybrid",las=1,cex.lab=cex_var,cex.axis=cex_var,cex.main=cex_var)
    dev.off()
    png(paste0("metropolis-", mut, "-", sce, "-mutationrate.png"))
    par(mar=par_vec)
    plot(1:length(metro[,1]), metro[,1], type="l",ylim=c(ymin,ymax),ylab="Mutation rate",xlab="Step",main="Metropolis",las=1,cex.lab=cex_var,cex.axis=cex_var,cex.main=cex_var)
    dev.off()

    ymin <- floor(min(metro[,2], hybrid[,2], zigzag[,2]))
    ymax <- ceiling(max(metro[,2], hybrid[,2], zigzag[,2]))
    png(paste0("zigzag-", mut, "-", sce, "-treeheight.png"))
    par(mar=par_vec)
    plot(times_zigzag, zigzag[,2], type="l",ylim=c(ymin,ymax),ylab="Tree height",xlab="Time",main="Zig-zag",las=1,cex.lab=cex_var,cex.axis=cex_var,cex.main=cex_var)
    dev.off()
    png(paste0("hybrid-", mut, "-", sce, "-treeheight.png"))
    par(mar=par_vec)
    plot(times_hybrid, hybrid[,2], type="l",ylim=c(ymin,ymax),ylab="Tree height",xlab="Time",main="Hybrid",las=1,cex.lab=cex_var,cex.axis=cex_var,cex.main=cex_var)
    dev.off()
    png(paste0("metropolis-", mut, "-", sce, "-treeheight.png"))
    par(mar=par_vec)
    plot(1:length(metro[,2]), metro[,2], type="l",ylim=c(ymin,ymax),ylab="Tree height",xlab="Step",main="Metropolis",las=1,cex.lab=cex_var,cex.axis=cex_var,cex.main=cex_var)
    dev.off()

    vels <- matrix(rep(NA, (n_zz - 1) * length(zigzag[1,])), ncol = length(zigzag[1,]))
    for (i in 1:(n_zz - 1)) {
      if (times_zigzag[i + 1] > times_zigzag[i]) {
        vels[i,] <- as.numeric((zigzag[i + 1,] - zigzag[i,]) / (times_zigzag[i + 1] - times_zigzag[i]))
      } else {
        vels[i,] <- rep(0, length(zigzag[1,]))
      }
    }
    means <- colSums(diff(times_zigzag) * (zigzag[1:(n_zz - 1),] + diff(times_zigzag) * vels / 2)) / times_zigzag[n_zz]
    vars <- colSums(diff(times_zigzag) * (zigzag[1:(n_zz - 1),]^2 + zigzag[1:(n_zz - 1),] * vels * diff(times_zigzag) + diff(times_zigzag)^2 * vels^2 / 3)) / times_zigzag[n_zz] - means * means
    batch_means <- matrix(rep(NA, batches * length(zigzag[1,])), ncol = length(zigzag[1,]))
    for (i in 1:batches) {
      inds <- (1:n_zz)[times_zigzag <= i * times_zigzag[n_zz] / batches & times_zigzag >= (i - 1) * times_zigzag[n_zz] / batches]
      batch_means[i,] <- sqrt(1 / (times_zigzag[inds[length(inds)]] - times_zigzag[inds[1]])) * colSums(diff(times_zigzag[inds]) * (zigzag[inds[1:(length(inds) - 1)],] + vels[inds[1:(length(inds) - 1)],] * diff(times_zigzag[inds]) / 2))
    }
    ess_zigzag <- times_zigzag[n_zz] * vars / apply(batch_means, 2, var)

    ncol <- dim(hybrid)[2] - 1
    vels <- matrix(rep(NA, (n_hy - 1) * ncol), ncol = ncol)
    for (i in 1:(n_hy - 1)) {
      if (times_hybrid[i + 1] > times_hybrid[i]) {
        vels[i,] <- as.numeric((hybrid[i + 1,1:ncol] - hybrid[i,1:ncol]) / (times_hybrid[i + 1] - times_hybrid[i]))
      } else {
        vels[i,] <- rep(0, ncol)
      }
    }
    means <- diff(times_hybrid) %*% (hybrid[2:n_hy,1:ncol] + hybrid[1:(n_hy-1),1:ncol]) / (2 * times_hybrid[n_hy])
    tmp <- hybrid[1:(n_hy-1),1:2] - sweep(vels, 1, times_hybrid[1:(n_hy-1)], "*")
    vars <- colSums(sweep(tmp^2, 1, diff(times_hybrid), "*") + sweep(vels * tmp, 1, times_hybrid[2:n_hy]^2 - times_hybrid[1:(n_hy-1)]^2, "*") + sweep(vels^2, 1, times_hybrid[2:n_hy]^3 - times_hybrid[1:(n_hy-1)]^3, "*") / 3) / times_hybrid[n_hy] - means * means
    batch_means <- matrix(rep(NA, batches * ncol), ncol = ncol)
    for (i in 1:batches) {
      inds <- (1:n_hy)[times_hybrid <= i * times_hybrid[n_hy] / batches & times_hybrid >= (i - 1) * times_hybrid[n_hy] / batches]
      if (inds[1] > 1) {
        inds <- c(inds[1] - 1, inds)
      }
      if (inds[length(inds)] < n_hy) {
        inds <- c(inds, inds[length(inds)] + 1)
      }
      tmp_times <- times_hybrid[inds]
      tmp_z <- hybrid[inds,1:ncol]
      if (i > 1) {
        tmp_times[1] <- (i - 1) * times_hybrid[n_hy] / batches
        tmp_z[1,] <- hybrid[inds[2],1:ncol] - vels[inds[1],] * (times_hybrid[inds[2]] - (i - 1) * times_hybrid[n_hy] / batches)
      }
      tmp_times[length(tmp_times)] <- i * times_hybrid[n_hy] / batches
      if (inds[length(inds)] < n_hy) {
        tmp_z[length(tmp_times),1:ncol] <- hybrid[inds[length(inds) - 1],1:ncol] + vels[inds[length(inds)],] * (i * times_hybrid[n_hy] / batches - times_hybrid[inds[length(inds)] - 1])
      }
      batch_means[i,] <- sqrt(batches / times_hybrid[n_hy]) * colSums(sweep(tmp_z[1:(length(tmp_times) - 1),], 1, diff(tmp_times), "*") - sweep(vels[inds[1:(length(inds) - 1)],], 1, tmp_times[1:(length(tmp_times) - 1)] * diff(tmp_times), "*") + sweep(vels[inds[1:(length(inds) - 1)],], 1, tmp_times[2:length(tmp_times)]^2 - tmp_times[1:(length(tmp_times) - 1)]^2, "*") / 2)
    }
    ess_hybrid <- times_hybrid[n_hy] * vars / apply(batch_means, 2, var)

    write(ess_zigzag, file=paste0("ess-", mut, "-", sce, ".txt"))
    write(ess_hybrid, file=paste0("ess-", mut, "-", sce, ".txt"), append=TRUE)
    write(as.numeric(ess(metro)), file=paste0("ess-", mut, "-", sce, ".txt"), append=TRUE)
    k <- k + 1
  }
}
toc()