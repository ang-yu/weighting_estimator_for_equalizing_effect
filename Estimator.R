

# R1, R2 must be binary variables in numeric format, and having 1 indicate membership, 0 indicate non-membership
# R1 is the target group (the advantaged group), R2 is the base group (the disadvantaged group)
# make sure sum(data$R1==1 & data$R2==1)=0
# Y must also be numeric
# there shouldn't be a variable in more than one of the sets Y,W,R1,R2,Q,L,C  

equalize <- function(Y, W, R1, R2, Q=NULL, L=NULL, C=NULL, data, metric="Risk Difference", K=1000, alpha=0.05) {
  
  data_inuse <- na.omit(data[,c(Y, W, R1, R2, Q, L, C)])
  # Even if R1 and R2 do not exhaust the data, the cases in neither R1 nor R2 are still retained.
  # Otherwise, when C is present, the original disparity between R1 and R2 plus the original disparity between R2 and R3 will not equal the OD between R1 and R3.
  
  Q <- paste(Q,collapse="+")
  C <- paste(C,collapse="+")
  L <- paste(L,collapse="+")
  
  if (Q!="" & L!="" & C!="") {
    nume_formula <- as.formula(paste(W, paste(R1,R2,Q,C,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(R1,R2,Q,C,L,sep="+"),sep="~"))
    adjustment_R1_formula <- as.formula(paste(R1,C,sep="~"))
    adjustment_R2_formula <- as.formula(paste(R2,C,sep="~"))
  } else if (Q!="" & L!="") {
    nume_formula <- as.formula(paste(W, paste(R1,R2,Q,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(R1,R2,Q,L,sep="+"),sep="~"))
  } else if (Q!="" & C!="") {
    nume_formula <- as.formula(paste(W, paste(R1,R2,Q,C,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(R1,R2,Q,C,sep="+"),sep="~"))
    adjustment_R1_formula <- as.formula(paste(R1,C,sep="~"))
    adjustment_R2_formula <- as.formula(paste(R2,C,sep="~"))
  } else if (L!="" & C!="") {
    nume_formula <- as.formula(paste(W, paste(R1,R2,C,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(R1,R2,C,L,sep="+"),sep="~"))
    adjustment_R1_formula <- as.formula(paste(R1,C,sep="~"))
    adjustment_R2_formula <- as.formula(paste(R2,C,sep="~"))
  } else if (Q!="") {
    nume_formula <- as.formula(paste(W, paste(R1,R2,Q,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(R1,R2,Q,sep="+"),sep="~"))
  } else if (L!="") {
    nume_formula <- as.formula(paste(W, paste(R1,R2,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(R1,R2,L,sep="+"),sep="~"))
  } else if (C!="") {
    nume_formula <- as.formula(paste(W, paste(R1,R2,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(R1,R2,C,sep="+"),sep="~"))
    adjustment_R1_formula <- as.formula(paste(R1,C,sep="~"))
    adjustment_R2_formula <- as.formula(paste(R2,C,sep="~"))
  } else {
    nume_formula <- as.formula(paste(W, paste(R1,R2,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(R1,R2,sep="+"),sep="~"))
  }
  
  nume_mol <- glm(nume_formula, family=binomial(link = "logit"), data=data_inuse)
  nume_newdata <- data_inuse
  nume_newdata[,R1] <- 1
  nume_newdata[,R2] <- 0
  suppressWarnings( nume_pred <- predict(nume_mol, newdata = nume_newdata, type = "response") )
  
  deno_mol <- glm(deno_formula, family=binomial(link = "logit"), data=data_inuse)
  deno_newdata <- data_inuse
  deno_newdata[,R1] <- 0
  deno_newdata[,R2] <- 1
  suppressWarnings( deno_pred <- predict(deno_mol, newdata = deno_newdata, type = "response") )
  
  if (C!="") {
    suppressWarnings( adjustment_R1_pred <- predict(glm(adjustment_R1_formula, family=binomial(link = "logit"), data=data_inuse), type="response") )
    suppressWarnings( adjustment_R2_pred <- predict(glm(adjustment_R2_formula, family=binomial(link = "logit"), data=data_inuse), type="response") )
    original_R1 <- mean(data_inuse[,Y][data_inuse[,R1]==1]*mean(data_inuse[,R1]==1)/adjustment_R1_pred[data_inuse[,R1]==1])
    original_R2 <- mean(data_inuse[,Y][data_inuse[,R2]==1]*mean(data_inuse[,R2]==1)/adjustment_R2_pred[data_inuse[,R2]==1])
    post_R2 <- mean(data_inuse[,Y][data_inuse[,R2]==1]*(nume_pred/deno_pred)[data_inuse[,R2]==1]*mean(data_inuse[,R2]==1)/adjustment_R2_pred[data_inuse[,R2]==1])
  } else {
    original_R1 <- mean(data_inuse[,Y][data_inuse[,R1]==1])
    original_R2 <- mean(data_inuse[,Y][data_inuse[,R2]==1])
    post_R2 <- mean(data_inuse[,Y][data_inuse[,R2]==1]*(nume_pred/deno_pred)[data_inuse[,R2]==1])
  }
  
  if (metric=="Risk Difference") {
    original <- original_R1-original_R2
    remaining <- original_R1-post_R2
    reduction <- post_R2-original_R2
  } else if (metric=="Risk Ratio") {
    original <- original_R1/original_R2
    remaining <- original_R1/post_R2
    reduction <- post_R2/original_R2
  }
  
  boot_original <- boot_remaining <- boot_reduction <- rep(NA, K)
  set.seed(36)
  
  for (i in 1:K) {
    
    indices <- sample(1:nrow(data_inuse), nrow(data_inuse), replace = TRUE)
    data_boot <- data_inuse[indices,]
    
    nume_mol <- glm(nume_formula, family=binomial(link = "logit"), data=data_boot)
    nume_newdata <- data_boot
    nume_newdata[,R1] <- 1
    nume_newdata[,R2] <- 0
    suppressWarnings( nume_pred <- predict(nume_mol, newdata = nume_newdata, type = "response") )
    
    deno_mol <- glm(deno_formula, family=binomial(link = "logit"), data=data_boot)
    deno_newdata <- data_boot
    deno_newdata[,R1] <- 0
    deno_newdata[,R2] <- 1
    suppressWarnings( deno_pred <- predict(deno_mol, newdata = deno_newdata, type = "response") )
    
    if (C!="") {
      suppressWarnings( adjustment_R1_pred <- predict(glm(adjustment_R1_formula, family=binomial(link = "logit"), data=data_boot), type="response") )
      suppressWarnings( adjustment_R2_pred <- predict(glm(adjustment_R2_formula, family=binomial(link = "logit"), data=data_boot), type="response") )
      original_R1 <- mean(data_boot[,Y][data_boot[,R1]==1]*mean(data_boot[,R1]==1)/adjustment_R1_pred[data_boot[,R1]==1])
      original_R2 <- mean(data_boot[,Y][data_boot[,R2]==1]*mean(data_boot[,R2]==1)/adjustment_R2_pred[data_boot[,R2]==1])
      post_R2 <- mean(data_boot[,Y][data_boot[,R2]==1]*(nume_pred/deno_pred)[data_boot[,R2]==1]*mean(data_boot[,R2]==1)/adjustment_R2_pred[data_boot[,R2]==1])
    } else {
      original_R1 <- mean(data_boot[,Y][data_boot[,R1]==1])
      original_R2 <- mean(data_boot[,Y][data_boot[,R2]==1])
      post_R2 <- mean(data_boot[,Y][data_boot[,R2]==1]*(nume_pred/deno_pred)[data_boot[,R2]==1])
    }
    
    if (metric=="Risk Difference") {
      original_boot <- original_R1-original_R2
      remaining_boot <- original_R1-post_R2
      reduction_boot <- post_R2-original_R2
    } else if (metric=="Risk Ratio") {
      original_boot <- original_R1/original_R2
      remaining_boot <- original_R1/post_R2
      reduction_boot <- post_R2/original_R2
    }
    
    boot_original[i] <- original_boot
    boot_remaining[i] <- remaining_boot
    boot_reduction[i] <- reduction_boot
    
  }
  
  boot_original_sd <- sqrt(mean((boot_original-mean(boot_original))^2)) # use this formula instead of sd() so that the denominator is K instead of K-1
  boot_remaining_sd <- sqrt(mean((boot_remaining-mean(boot_remaining))^2))
  boot_reduction_sd <- sqrt(mean((boot_reduction-mean(boot_reduction))^2))
  
  # Basic bootstrap confidence intervals
  boot_original_lower <- 2*original-boot_original[order(boot_original)][floor(K*(1-alpha/2))]
  boot_original_upper <- 2*original-boot_original[order(boot_original)][ceiling(K*alpha/2)]
  
  boot_remaining_lower <- 2*remaining-boot_remaining[order(boot_remaining)][floor(K*(1-alpha/2))]
  boot_remaining_upper <- 2*remaining-boot_remaining[order(boot_remaining)][ceiling(K*alpha/2)]
  
  boot_reduction_lower <- 2*reduction-boot_reduction[order(boot_reduction)][floor(K*(1-alpha/2))]
  boot_reduction_upper <- 2*reduction-boot_reduction[order(boot_reduction)][ceiling(K*alpha/2)]
  
  boot_original_ci <- paste("(",paste(sprintf("%.4f",boot_original_lower),sprintf("%.4f",boot_original_upper),sep=","),")",sep="")
  boot_remaining_ci <- paste("(",paste(sprintf("%.4f",boot_remaining_lower),sprintf("%.4f",boot_remaining_upper),sep=","),")",sep="")
  boot_reduction_ci <- paste("(",paste(sprintf("%.4f",boot_reduction_lower),sprintf("%.4f",boot_reduction_upper),sep=","),")",sep="")
  
  output <- matrix(nrow=3,ncol=3,dimnames=list(c("Original Disparity","Remaining Disparity","Reduction in Disparity"),c("Estimate","Std. Error","Conf. Int.")))
  output[1,] <- c(sprintf("%.4f",original),sprintf("%.4f",boot_original_sd),boot_original_ci)
  output[2,] <- c(sprintf("%.4f",remaining),sprintf("%.4f",boot_remaining_sd),boot_remaining_ci)
  output[3,] <- c(sprintf("%.4f",reduction),sprintf("%.4f",boot_reduction_sd),boot_reduction_ci)
  
  return(as.data.frame(output))
}


