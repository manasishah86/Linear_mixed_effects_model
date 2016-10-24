#We now consider if there is any change in beta-diversity over the course of prengnacy, 
#where beta-diversity is defined here as the average distance between subjects (at that gestational time).
divdfs <- list()
ANSIZE=2   # Plot visual parameter
PVALLEN=5

FMT_Types <- c("Fresh", "Frozen", "Lyophilized")

FMT_ps_list = lapply(FMT_Types, function(x, input_physeq = CDI_Recipient){
  prune_samples(get_variable(input_physeq, "FMT_Type") == x, input_physeq)
})

names(FMT_ps_list) <- FMT_Types

require("doParallel")
require("foreach")
registerDoParallel(makeCluster(3))
# Function to pull out SampleIDs by the Time
getSamByGW <- function(x, cc) {
  if(x %in% cc$Time) { 
    return(cc$SampleID[which(x==cc$Time)][[1]])
  } else { return(NA) }
}

# Function to calculate distance between two samples
getDist2Sam <- function(sam1, sam2, ps, ...) {
  phyloseq::distance(prune_samples(c(sam1, sam2), ps), ...)[[1]]
}

FMTSTART=0; FMTEND=30
time <- c(0, 7, 14, 30)


# Function to get lm coeffs for the beta-diversity vs. gw fit
get_mb <- function(ps, gwdists=NULL, permute=F) {
  samdf <- data.frame(sample_data(ps))
  subdfs <- split(samdf, samdf$SampleID)
  
  if(permute) {
    for(sub in names(subdfs)) {
      if(runif(1) > 0.5) {
        # Reversing in the GW10-40 range
        subdfs[[sub]]$Time <- (FMTSTART+FMTEND-subdfs[[sub]]$Time)
      }
    }
  }
  
  betas <- seq_along(time)
  names(betas) <- time
  for(gw in time) {
    sams <- sapply(subdfs, function(x) getSamByGW(gw, x))
    sams <- sams[!is.na(sams)]
    if(is.null(gwdists)) {
      if(meth=="weighted unifrac") {
        dists <- phyloseq::distance(prune_samples(sams, ps), method="unifrac", weighted=T, parallel=T)
      } else {
        dists <- phyloseq::distance(prune_samples(sams, ps), method=meth, parallel=T)
      }
    } else {
      dists <- gwdists[[as.character(gw)]]
      dists <- as.matrix(dists)
      #      print(gw)
      #      print(sams)
      dists <- dists[sams, sams]
      dists <- as.dist(dists)
    }
    
    betas[[as.character(gw)]] <- mean(dists)
  }
  lmo <- lm(beta ~ gw, data.frame(beta=betas, gw=time))
  m <- lmo$coefficients[["gw"]]
  b <- lmo$coefficients[["(Intercept)"]]
  return(list(m=m, b=b, lmo=lmo))
}

meths <- c("bray", "jsd")
par(mfrow=c(2,2))
NPERM=1000
betadivs <- data.frame()
mbs <- data.frame()
mbs_perm <- data.frame()
for(meth in meths) {
  for(FMT_Type in FMT_Types) {
    
    ps <- FMT_ps_list[[FMT_Type]]
    ps <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
    samdf <- data.frame(sample_data(ps))
    
    # Calculate the distance matrix for each pair of complementary time
    
######################################################################################################################################    
    getDist <- function(gw, meth) phyloseq::distance(prune_samples(samdf$Time %in% c(gw, FMTSTART+FMTEND-gw), ps), method=meth, parallel=T)
    getWUF <- function(gw) phyloseq::distance(prune_samples(samdf$Time %in% c(gw, FMTSTART+FMTEND-gw), ps), method="unifrac", weighted=T, parallel=T)
    if(meth == "wuf") gwdists <- lapply(time, getWUF)
    else gwdists <- lapply(time, getDist, meth=meth)
    names(gwdists) <- time
    
    fit <- get_mb(ps, gwdists=gwdists, permute=F)
    m <- fit$m;b <- fit$b
    
    # Do the time-reversal permutations
    pms <- rep(-999, NPERM); pbs <- rep(-999, NPERM)
    for(ii in seq(1,NPERM)) {
      pfit <- get_mb(ps, gwdists=gwdists, permute=T)
      pms[[ii]] <- pfit$m
      pbs[[ii]] <- pfit$b
    }
    
    # plot the slope (m) estimated from the data along with the histograms of ms from
    # the permutations. Show pval calculated from that.
    betas <- fit$lmo$model$beta
    plot(time, betas, ylab=meth, main=FMT_Type)
    abline(a=fit$b, b=fit$m, col="blue")
    pval <- mean(pms < -abs(m)) + mean(pms > abs(m))
    hist(pms, breaks=30, main=paste(FMT_Type, meth, "p:", pval))
    abline(v=m, col="red", lw=3)
    
    # Save the linear fit params and calculated beta-diversities
    mbs <- rbind(mbs, data.frame(FMT_Type=FMT_Type, meth=meth, m=m, b=b, pval=pval))
    mbs_perm <- rbind(mbs_perm, data.frame(FMT_Type=FMT_Type, meth=meth, m=pms, b=pbs))
    betadivs <- rbind(betadivs, data.frame(FMT_Type=FMT_Type, meth=meth, Time=time, beta=betas))
  }
}


meths <- c("bray", "jsd")
preds <- t(mapply(function(m,b) m*time+b, mbs_perm$m, mbs_perm$b))
colnames(preds) <- as.character(time)
df_perm <- cbind(mbs_perm, preds)
dfs <- split(df_perm, list(df_perm$FMT_Type, df_perm$meth))
betadivs$lo <- c(NA)
betadivs$hi <- c(NA)
split_betas <- split(betadivs, list(betadivs$FMT_Type, betadivs$meth))
split_mbs <- split(mbs, list(mbs$FMT_Type, mbs$meth))
for(FMT_Type in FMT_Types) {
  for(meth in meths) {
    key <- paste0(FMT_Type,".",meth)
    df <- dfs[[key]]
    predvar <- apply(df[,as.character(time)], 2, var)
    lo <- apply(df[,as.character(time)], 2, function(x) quantile(x, 0.025))
    hi <- apply(df[,as.character(time)], 2, function(x) quantile(x, 0.975))
    pred <- split_mbs[[key]][[1,"m"]]*time + split_mbs[[key]][[1,"b"]]
    split_betas[[key]]$lo <- lo
    split_betas[[key]]$hi <- hi
    split_betas[[key]]$pred <- pred
    split_betas[[key]]$SE <- sqrt(predvar)
  }
}
beta_divs <- unsplit(split_betas, list(betadivs$FMT_Type, betadivs$meth))


meths <- c("bray", "jsd")
pvs <- list(bray=list(ymin=0, ymax=1, ylab="Beta-diversity (Bray)", xan = 25, yan=1),
            jsd=list(ymin=0, ymax=1, ylab="Beta-diversity (JSD)", xan = 25, yan=1))
            #wuf=list(ymin=0, ymax=0.5, ylab="Beta-diversity (WUF)", xan = 30, yan=0.48))

b_fits <- list()
for(FMT_Type in FMT_Types) {
  df <- divdfs[[FMT_Type]]
  for(meth in meths) {
    dfkey <- paste0(FMT_Type,".",meth)
    key <- paste(FMT_Type,meth)
    df <- split_betas[[dfkey]]
    pv <- pvs[[meth]]
    pval <- mbs[[which(mbs$FMT_Type==FMT_Type & mbs$meth==meth), "pval"]]
    
    # LM modeling/plotting
    pb <- ggplot(data=df, aes(x=Time)) + geom_point(aes(y=beta), size=0.5)
    pb <- pb + geom_line(aes(y=pred), color="red")
    pb <- pb + geom_ribbon(aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.2, fill="green")
    pb <- pb + ylab(pv$ylab) + ylim(pv$ymin, pv$ymax)
    pb <- pb + xlab("FMT Timeline (Days)")
    pb <- pb + annotate("text", x=pv$xan, y=pv$yan, size=5, label = paste("p =", substr(as.character(pval),1,PVALLEN))) + ggtitle(FMT_Type)
    #pb <- pb + ylab(NULL) + xlab(NULL)
    b_fits[[key]] <- pb
  }
}
library("gridExtra")
listorder <- paste(FMT_Types, rep(meths, each=3))
b_fits <- b_fits[listorder]

do.call(grid.arrange, b_fits)

tiff("betadiv_timeseries.tiff", width = 40, height = 30, units = "cm", res = 300,
     compression = "lzw", colortype = "true")
theme_set(theme_bw(base_size = 18))
theme_update(axis.text = element_text(size = 10))
do.call(grid.arrange, b_fits)
invisible(dev.off())


#grid.arrange(b_fits,  ncol=2, top="Beta-diversity:", "Bray", "JSD")
##################################################################################################################################################################################################################################################################


FMTSTART=7; FMTEND=30
time <- c(7, 14, 30)


# Function to get lm coeffs for the beta-diversity vs. gw fit
get_mb <- function(ps, gwdists=NULL, permute=F) {
  samdf <- data.frame(sample_data(ps))
  subdfs <- split(samdf, samdf$SampleID)
  
  if(permute) {
    for(sub in names(subdfs)) {
      if(runif(1) > 0.5) {
        # Reversing in the GW10-40 range
        subdfs[[sub]]$Time <- (FMTSTART+FMTEND-subdfs[[sub]]$Time)
      }
    }
  }
  
  betas <- seq_along(time)
  names(betas) <- time
  for(gw in time) {
    sams <- sapply(subdfs, function(x) getSamByGW(gw, x))
    sams <- sams[!is.na(sams)]
    if(is.null(gwdists)) {
      if(meth=="weighted unifrac") {
        dists <- phyloseq::distance(prune_samples(sams, ps), method="unifrac", weighted=T, parallel=T)
      } else {
        dists <- phyloseq::distance(prune_samples(sams, ps), method=meth, parallel=T)
      }
    } else {
      dists <- gwdists[[as.character(gw)]]
      dists <- as.matrix(dists)
      #      print(gw)
      #      print(sams)
      dists <- dists[sams, sams]
      dists <- as.dist(dists)
    }
    
    betas[[as.character(gw)]] <- mean(dists)
  }
  lmo <- lm(beta ~ gw, data.frame(beta=betas, gw=time))
  m <- lmo$coefficients[["gw"]]
  b <- lmo$coefficients[["(Intercept)"]]
  return(list(m=m, b=b, lmo=lmo))
}

meths <- c("bray", "jsd")
par(mfrow=c(2,2))
NPERM=1000
betadivs <- data.frame()
mbs <- data.frame()
mbs_perm <- data.frame()
for(meth in meths) {
  for(FMT_Type in FMT_Types) {
    
    ps <- FMT_ps_list[[FMT_Type]]
    ps <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
    samdf <- data.frame(sample_data(ps))
    
    # Calculate the distance matrix for each pair of complementary time
    
    ######################################################################################################################################    
    getDist <- function(gw, meth) phyloseq::distance(prune_samples(samdf$Time %in% c(gw, FMTSTART+FMTEND-gw), ps), method=meth, parallel=T)
    getWUF <- function(gw) phyloseq::distance(prune_samples(samdf$Time %in% c(gw, FMTSTART+FMTEND-gw), ps), method="unifrac", weighted=T, parallel=T)
    if(meth == "wuf") gwdists <- lapply(time, getWUF)
    else gwdists <- lapply(time, getDist, meth=meth)
    names(gwdists) <- time
    
    fit <- get_mb(ps, gwdists=gwdists, permute=F)
    m <- fit$m;b <- fit$b
    
    # Do the time-reversal permutations
    pms <- rep(-999, NPERM); pbs <- rep(-999, NPERM)
    for(ii in seq(1,NPERM)) {
      pfit <- get_mb(ps, gwdists=gwdists, permute=T)
      pms[[ii]] <- pfit$m
      pbs[[ii]] <- pfit$b
    }
    
    # plot the slope (m) estimated from the data along with the histograms of ms from
    # the permutations. Show pval calculated from that.
    betas <- fit$lmo$model$beta
    plot(time, betas, ylab=meth, main=FMT_Type)
    abline(a=fit$b, b=fit$m, col="blue")
    pval <- mean(pms < -abs(m)) + mean(pms > abs(m))
    hist(pms, breaks=30, main=paste(FMT_Type, meth, "p:", pval))
    abline(v=m, col="red", lw=3)
    
    # Save the linear fit params and calculated beta-diversities
    mbs <- rbind(mbs, data.frame(FMT_Type=FMT_Type, meth=meth, m=m, b=b, pval=pval))
    mbs_perm <- rbind(mbs_perm, data.frame(FMT_Type=FMT_Type, meth=meth, m=pms, b=pbs))
    betadivs <- rbind(betadivs, data.frame(FMT_Type=FMT_Type, meth=meth, Time=time, beta=betas))
  }
}


meths <- c("bray", "jsd")
preds <- t(mapply(function(m,b) m*time+b, mbs_perm$m, mbs_perm$b))
colnames(preds) <- as.character(time)
df_perm <- cbind(mbs_perm, preds)
dfs <- split(df_perm, list(df_perm$FMT_Type, df_perm$meth))
betadivs$lo <- c(NA)
betadivs$hi <- c(NA)
split_betas <- split(betadivs, list(betadivs$FMT_Type, betadivs$meth))
split_mbs <- split(mbs, list(mbs$FMT_Type, mbs$meth))
for(FMT_Type in FMT_Types) {
  for(meth in meths) {
    key <- paste0(FMT_Type,".",meth)
    df <- dfs[[key]]
    predvar <- apply(df[,as.character(time)], 2, var)
    lo <- apply(df[,as.character(time)], 2, function(x) quantile(x, 0.025))
    hi <- apply(df[,as.character(time)], 2, function(x) quantile(x, 0.975))
    pred <- split_mbs[[key]][[1,"m"]]*time + split_mbs[[key]][[1,"b"]]
    split_betas[[key]]$lo <- lo
    split_betas[[key]]$hi <- hi
    split_betas[[key]]$pred <- pred
    split_betas[[key]]$SE <- sqrt(predvar)
  }
}
beta_divs <- unsplit(split_betas, list(betadivs$FMT_Type, betadivs$meth))


meths <- c("bray", "jsd")
pvs <- list(bray=list(ymin=0, ymax=1, ylab="Beta-diversity (Bray)", xan = 25, yan=1),
            jsd=list(ymin=0, ymax=1, ylab="Beta-diversity (JSD)", xan = 25, yan=1))
#wuf=list(ymin=0, ymax=0.5, ylab="Beta-diversity (WUF)", xan = 30, yan=0.48))

b_fits <- list()
for(FMT_Type in FMT_Types) {
  df <- divdfs[[FMT_Type]]
  for(meth in meths) {
    dfkey <- paste0(FMT_Type,".",meth)
    key <- paste(FMT_Type,meth)
    df <- split_betas[[dfkey]]
    pv <- pvs[[meth]]
    pval <- mbs[[which(mbs$FMT_Type==FMT_Type & mbs$meth==meth), "pval"]]
    
    # LM modeling/plotting
    pb <- ggplot(data=df, aes(x=Time)) + geom_point(aes(y=beta), size=0.5)
    pb <- pb + geom_line(aes(y=pred), color="red")
    pb <- pb + geom_ribbon(aes(ymin=pred-2*SE, ymax=pred+2*SE), alpha=0.2, fill="green")
    pb <- pb + ylab(pv$ylab) + ylim(pv$ymin, pv$ymax)
    pb <- pb + xlab("FMT Timeline (Days)")
    pb <- pb + annotate("text", x=pv$xan, y=pv$yan, size=5, label = paste("p =", substr(as.character(pval),1,PVALLEN))) + ggtitle(paste(FMT_Type, "Recipient at three time points"))
    #pb <- pb + ylab(NULL) + xlab(NULL)
    b_fits[[key]] <- pb
  }
}
library("gridExtra")
listorder <- paste(FMT_Types, rep(meths, each=3))
b_fits <- b_fits[listorder]

do.call(grid.arrange, b_fits)

tiff("betadiv_3time.tiff", width = 40, height = 30, units = "cm", res = 300,
     compression = "lzw", colortype = "true")
theme_set(theme_bw(base_size = 18))
theme_update(axis.text = element_text(size = 10))
do.call(grid.arrange, b_fits)
invisible(dev.off())




#pb <- pb + annotate("text", x=pv$xan, y=pv$yan, size=5, label = paste("p =", substr(as.character(pval),1,PVALLEN))) + ggtitle(paste(FMT_Type, "Recipient at three time point"))


#################################################################################################################################
#Another way to look for systematic trends during FMT is to ordinate, and indicate time via color
ords <- list()
for(FMT_Type in FMT_Types) {
  ps <- FMT_ps_list[[FMT_Type]]
  ords[[FMT_Type]] <- ordinate(ps, "NMDS", "bray")
}
require("RColorBrewer")
timepal <- colorRampPalette(brewer.pal(9,"RdBu"))(41)

contours <- list()
for(FMT_Type in FMT_Types) {
  ps <- FMT_ps_list[[FMT_Type]]
  pord <- plot_ordination(ps, ords[[FMT_Type]], type = "samples")
  #  print(pord + geom_point(aes(color=Time), size = 3, alpha=0.5) + ggtitle(paste(FMT_Type, "NMDS", "bray")))
  data <- pord$data
  data$Time_point <- c(NA)
  data$Time_point[data$Time == 0] <- "Baseline FMT"
  #data$Time_point[data$Time %in% c(7, 14, 30)] <- "Post FMT"
  #data$Time_point[data$Time %in% c(7, 14)] <- "Two weeks post FMT"
  data$Time_point[data$Time == 7] <- "One Week Post FMT"
  data$Time_point[data$Time == 14] <- "Two Weeks post FMT"
  data$Time_point[data$Time == 30] <- "One Month post FMT"
  cont <- ggplot(data=data, aes(x=NMDS1, y=NMDS2))
  cont <- cont + geom_point(size=0.5) + geom_density2d(data=data[!is.na(data$Time_point),], aes(color=Time_point)) + xlim(-2,2) + ylim(-2,2)
  contours[[FMT_Type]] <- cont + ggtitle(FMT_Type)
}
# SUPPLEMENTARY FIGURE
do.call(grid.arrange, c(contours, ncol=2))

tiff("contour_plot.tiff", width = 30, height = 20, units = "cm", res = 300,
     compression = "lzw", colortype = "true")
theme_set(theme_bw(base_size = 18))
theme_update(axis.text = element_text(size = 10))
do.call(grid.arrange, c(contours, ncol=2))
invisible(dev.off())



