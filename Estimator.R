

# R1, R2 must be binary variables in numeric format, and having 1 indicate membership, 0 indicate non-membership
# R1 is the target group (the advantaged group), R2 is the base group (the disadvantaged group)
# make sure sum(data$R1==1 & data$R2==1)=0
# Y must also be numeric
# there shouldn't be a variable in more than one of the sets Y,M,R1,R2,Q,L,C  

equalize <- function(Y, M, R1, R2, Q=NULL, L=NULL, C=NULL, data, metric="Risk Difference", K=1000) {
  
  data_inuse <- na.omit(data[,c(Y, M, R1, R2, Q, L, C)])
  # Even if R1 and R2 do not exhaust the data, the cases in neither R1 nor R2 are still retained.
  # Otherwise, when C is present, the original disparity between R1 and R2 plus the original disparity between R2 and R3 will not equal the OD between R1 and R3.
  
  Q <- paste(Q,collapse="+")
  C <- paste(C,collapse="+")
  L <- paste(L,collapse="+")
  
  if (Q!="" & L!="" & C!="") {
    nume_formula <- as.formula(paste(M, paste(R1,R2,Q,C,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(M, paste(R1,R2,Q,C,L,sep="+"),sep="~"))
    adjustment_R1_formula <- as.formula(paste(R1,C,sep="~"))
    adjustment_R2_formula <- as.formula(paste(R2,C,sep="~"))
  } else if (Q!="" & L!="") {
    nume_formula <- as.formula(paste(M, paste(R1,R2,Q,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(M, paste(R1,R2,Q,L,sep="+"),sep="~"))
  } else if (Q!="" & C!="") {
    nume_formula <- as.formula(paste(M, paste(R1,R2,Q,C,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(M, paste(R1,R2,Q,C,sep="+"),sep="~"))
    adjustment_R1_formula <- as.formula(paste(R1,C,sep="~"))
    adjustment_R2_formula <- as.formula(paste(R2,C,sep="~"))
  } else if (L!="" & C!="") {
    nume_formula <- as.formula(paste(M, paste(R1,R2,C,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(M, paste(R1,R2,C,L,sep="+"),sep="~"))
    adjustment_R1_formula <- as.formula(paste(R1,C,sep="~"))
    adjustment_R2_formula <- as.formula(paste(R2,C,sep="~"))
  } else if (Q!="") {
    nume_formula <- as.formula(paste(M, paste(R1,R2,Q,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(M, paste(R1,R2,Q,sep="+"),sep="~"))
  } else if (L!="") {
    nume_formula <- as.formula(paste(M, paste(R1,R2,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(M, paste(R1,R2,L,sep="+"),sep="~"))
  } else if (C!="") {
    nume_formula <- as.formula(paste(M, paste(R1,R2,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(M, paste(R1,R2,C,sep="+"),sep="~"))
    adjustment_R1_formula <- as.formula(paste(R1,C,sep="~"))
    adjustment_R2_formula <- as.formula(paste(R2,C,sep="~"))
  } else {
    nume_formula <- as.formula(paste(M, paste(R1,R2,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(M, paste(R1,R2,sep="+"),sep="~"))
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
    original <- sprintf("%.4f",original_R1-original_R2)
    remaining <- sprintf("%.4f",original_R1-post_R2)
    reduction <- sprintf("%.4f",post_R2-original_R2)
  } else if (metric=="Risk Ratio") {
    original <- sprintf("%.4f",original_R1/original_R2)
    remaining <- sprintf("%.4f",original_R1/post_R2)
    reduction <- sprintf("%.4f",post_R2/original_R2)
  }
  
  boot_original <- boot_remaining <- boot_reduction <- rep(NA, K)
  set.seed(36)
  
  for (i in 1:K) {
    
    indices <- sample(1:nrow(data_inuse), nrow(data_inuse), replace = T)
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
      original_boot <- sprintf("%.4f",original_R1-original_R2)
      remaining_boot <- sprintf("%.4f",original_R1-post_R2)
      reduction_boot <- sprintf("%.4f",post_R2-original_R2)
    } else if (metric=="Risk Ratio") {
      original_boot <- sprintf("%.4f",original_R1/original_R2)
      remaining_boot <- sprintf("%.4f",original_R1/post_R2)
      reduction_boot <- sprintf("%.4f",post_R2/original_R2)
    }
    
    boot_original[i] <- original_boot
    boot_remaining[i] <- remaining_boot
    boot_reduction[i] <- reduction_boot
    
  }
  
  boot_sd_original <- sprintf("%.4f",sd(boot_original))
  boot_sd_remaining <- sprintf("%.4f",sd(boot_remaining))
  boot_sd_reduction <- sprintf("%.4f",sd(boot_reduction))
  
  cat(paste(paste(" Original Disparity",original,sep=": "),paste("Standard Error",boot_sd_original,sep=": "),sep=", "),"\n",
      paste(paste("Remaining Disparity",remaining,sep=": "),paste("Standard Error",boot_sd_remaining,sep=": "),sep=", "),"\n",
      paste(paste("Reduction in Disparity",reduction,sep=": "),paste("Standard Error",boot_sd_reduction,sep=": "),sep=", "))
  
}

