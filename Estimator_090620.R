

equalize <- function(Y, W, R1, R2, Q=NULL, L=NULL, C=NULL, data, percent=100, metric="difference", K=1000, alpha=0.05, truncation_threshold=1, common_support=F) {
  
  options(error = NULL) # so that the error message when stopped doesn't invoke debugging interface. 
  options(scipen=999) # so that the weight output is not printed in scientific notation
  
  # check data
  if(!is.data.frame(data)) stop("data must be a data.frame.",call.=FALSE)
  if(missing(Y)) stop("Y must be provided.")
  if(missing(W)) stop("W must be provided.")
  if(missing(R1)) stop("R1 must be provided.")
  if(missing(R2)) stop("R2 must be provided.")
  
  data_inuse <- na.omit(data[,c(Y, W, R1, R2, Q, L, C)])
  # Even if R1 and R2 do not exhaust the data, the cases in neither R1 nor R2 are still retained.
  # Otherwise, when C is present, the original disparity between R1 and R2 plus the original disparity between R2 and R3 will not equal the original disparity between R1 and R3.
  # Get rid of rows with any infinite value
  data_inuse <- data_inuse[is.finite(rowSums(data_inuse[,unlist(lapply(data_inuse, is.numeric))])),]
  
  # check percent
  if(!is.numeric(percent)) stop("percent must be a numeric value.",call.=FALSE)
  if(!percent>0) stop("percent must be larger than 0.",call.=FALSE)
  
  # check metric
  if (!metric%in%c("difference","risk ratio","odds ratio")) stop("metric must be difference, risk ratio, or odds ratio.",call.=FALSE)
  
  # check truncation_threshold
  if (!is.numeric(truncation_threshold)) stop("truncation_threshold must be a numeric value.",call.=FALSE)
  if (truncation_threshold<=0 | truncation_threshold>1) stop("truncation_threshold must be in (0,1].",call.=FALSE)
  
  # check alpha
  if (!is.numeric(alpha)) stop("alpha must be a numeric value.",call.=FALSE)
  if (alpha<=0 | alpha>=1) stop("alpha must be in (0,1).",call.=FALSE)
  
  # check Y
  if(!is.character(Y)) stop("Y must be a character scalar vector.",call.=FALSE)
  if(length(Y)>1) stop("Y must be only one variable",call.=FALSE)
  if(!is.numeric(data_inuse[,Y])) stop("Y must be numeric",call.=FALSE) 
  if(any(!unique(data_inuse[,Y])%in%c(0,1)) & metric!="difference") stop("Y must only have values of 0 and 1 when the outcome metric is risk ratio or odds ratio.",call.=FALSE) 

  # check W
  if(!is.character(W)) stop("W must be a character scalar vector.",call.=FALSE)
  if(length(W)>1) stop("W must be only one variable.",call.=FALSE)
  if(length(unique(data_inuse[,W]))!=2) stop("W must be binary.",call.=FALSE)
  if(!is.numeric(data_inuse[,W])) stop("W must be numeric.",call.=FALSE)
  if(max(data_inuse[,W])!=1 | min(data_inuse[,W])!=0) stop("W must only have values of 0 and 1.",call.=FALSE)
  
  # check R1 and R2
  if(!is.character(R1)) stop("R1 must be a character scalar vector.",call.=FALSE)
  if(length(R1)>1) stop("R1 must be only one variable.",call.=FALSE)
  if(length(unique(data_inuse[,R1]))!=2 & length(unique(data_inuse[,R1]))!=1) stop("R1 must be either binary or constant.",call.=FALSE)
  if(!is.numeric(data_inuse[,R1])) stop("R1 must be numeric",call.=FALSE) 
  if((max(data_inuse[,R1])!=1 | min(data_inuse[,R1])!=0) & length(unique(data_inuse[,R1]))==2) stop("R1 must have values of either 0 and 1 or just 1 .",call.=FALSE)
  if(length(unique(data_inuse[,R1]))==1) if(unique(data_inuse[,R1]!=1)) stop("R1 must have values of either 0 and 1 or just 1 .",call.=FALSE)
  # It's fine for R1 to be constant 1 in the sample (when the target level of W is the sample mean).
  # When R1==1 for all cases, the model will work as intended. 
  
  if(!is.character(R2)) stop("R2 must be a character scalar vector.",call.=FALSE)
  if(length(R2)>1) stop("R2 must be only one variable.",call.=FALSE)
  if(length(unique(data_inuse[,R2]))!=2) stop("R2 must be binary.",call.=FALSE)
  if(!is.numeric(data_inuse[,R2])) stop("R2 must be numeric",call.=FALSE) 
  if(max(data_inuse[,R2])!=1 | min(data_inuse[,R2])!=0) stop("R2 must only have values of 0 and 1.",call.=FALSE)
  # note that R1 and R2 are allowed to be overlapped. For example, R1 may be intervened to have the average W across R1 and R2.

  # check Q, L, and C
  if(!is.null(Q) & !is.character(Q)) stop("Q must be a character vector.",call.=FALSE)
  if(!is.null(L) & !is.character(L)) stop("L must be a character vector.",call.=FALSE)
  if(!is.null(C) & !is.character(C)) stop("C must be a character vector.",call.=FALSE)
  
  # check if there is overlap between the covariates vectors
  if (!is.null(Q) & !is.null(L) & !is.null(C)) {
    if(any(
      any(Y==W, Y==R1, Y==R2, Y%in%Q, Y%in%L, Y%in%C),
      any(W==R1, W==R2, W%in%Q, W%in%L, W%in%C),
      any(R1%in%Q, R1%in%L, R1%in%C),
      any(Q%in%L, Q%in%C),
      any(L%in%C)
    )) stop("Each variable can only have one role.",call.=FALSE)
    # note that R1==R2 is allowed. It is at least possible to conceive of an intervention to make R2 members have R2 average of W, regardless of L.
  } else if (!is.null(Q) & !is.null(L)) {
    if(any(
      any(Y==W, Y==R1, Y==R2, Y%in%Q, Y%in%L),
      any(W==R1, W==R2, W%in%Q, W%in%L),
      any(R1%in%Q, R1%in%L),
      any(Q%in%L)
    )) stop("Each variable can only have one role.",call.=FALSE)
  } else if (!is.null(Q) & !is.null(C)) {
    if(any(
      any(Y==W, Y==R1, Y==R2, Y%in%Q, Y%in%C),
      any(W==R1, W==R2, W%in%Q, W%in%C),
      any(R1%in%Q, R1%in%C),
      any(Q%in%C)
    )) stop("Each variable can only have one role.",call.=FALSE)
  } else if (!is.null(L) & !is.null(C)) {
    if(any(
      any(Y==W, Y==R1, Y==R2, Y%in%L, Y%in%C),
      any(W==R1, W==R2, W%in%L, W%in%C),
      any(R1%in%L, R1%in%C),
      any(L%in%C)
    )) stop("Each variable can only have one role.",call.=FALSE)
  } else if (!is.null(Q)) {
    if(any(
      any(Y==W, Y==R1, Y==R2, Y%in%Q),
      any(W==R1, W==R2, W%in%Q),
      any(R1%in%Q)
    )) stop("Each variable can only have one role.",call.=FALSE)
  } else if (!is.null(L)) {
    if(any(
      any(Y==W, Y==R1, Y==R2, Y%in%L),
      any(W==R1, W==R2, W%in%L),
      any(R1%in%L)
    )) stop("Each variable can only have one role.",call.=FALSE)
  } else if (!is.null(C)) {
    if(any(
      any(Y==W, Y==R1, Y==R2, Y%in%C),
      any(W==R1, W==R2, W%in%C),
      any(R1%in%C)
    )) stop("Each variable can only have one role.",call.=FALSE)
  } else {
    if(any(
      any(Y==W, Y==R1, Y==R2),
      any(W==R1, W==R2)
    )) stop("Each variable can only have one role.",call.=FALSE)
  }
  

  # check common_support
  if (!is.logical(common_support)) stop("common_support must be logical.",call.=FALSE)
  if (is.null(Q) & is.null(C) & isTRUE(common_support)) stop("common support restriction should only be applied when Q or C is specified.",call.=FALSE)
  if (sum(data_inuse[,R1]==1 & data_inuse[,R2]==1)==sum(data_inuse[,R1]==1) & isTRUE(common_support)) stop("common support restriction should only be applied when R2 is not strictly a subset of R1.",call.=FALSE)
  
  R1_Q <- paste(sapply(Q, function(x) paste(R1,x,sep=":")), collapse="+")
  R2_Q <- paste(sapply(Q, function(x) paste(R2,x,sep=":")), collapse="+")
  
  R1_C <- paste(sapply(C, function(x) paste(R1,x,sep=":")), collapse="+")
  R2_C <- paste(sapply(C, function(x) paste(R2,x,sep=":")), collapse="+")
  
  R1_L <- paste(sapply(L, function(x) paste(R1,x,sep=":")), collapse="+")
  R2_L <- paste(sapply(L, function(x) paste(R2,x,sep=":")), collapse="+")
  
  Q <- paste(Q,collapse="+")
  C <- paste(C,collapse="+")
  L <- paste(L,collapse="+")
  
  # generate formulas needed according to the presence or absence of Q, L, and C.
  if (Q!="" & L!="" & C!="") {
    nume_formula <- as.formula(paste(W, paste(R1,R2,Q,C,R1_Q,R2_Q,R1_C,R2_C,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(R1,R2,Q,C,L,R1_Q,R2_Q,R1_C,R2_C,R1_L,R2_L,sep="+"),sep="~"))
    adjustment_R1_formula <- as.formula(paste(R1,C,sep="~"))
    adjustment_R2_formula <- as.formula(paste(R2,C,sep="~"))
  } else if (Q!="" & L!="") {
    nume_formula <- as.formula(paste(W, paste(R1,R2,Q,R1_Q,R2_Q,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(R1,R2,Q,L,R1_L,R2_L,R1_Q,R2_Q,sep="+"),sep="~"))
  } else if (Q!="" & C!="") {
    nume_formula <- as.formula(paste(W, paste(R1,R2,Q,C,R1_Q,R2_Q,R1_C,R2_C,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(R1,R2,Q,C,R1_Q,R2_Q,R1_C,R2_C,sep="+"),sep="~"))
    adjustment_R1_formula <- as.formula(paste(R1,C,sep="~"))
    adjustment_R2_formula <- as.formula(paste(R2,C,sep="~"))
  } else if (L!="" & C!="") {
    nume_formula <- as.formula(paste(W, paste(R1,R2,C,R1_C,R2_C,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(R1,R2,C,L,R1_C,R2_C,R1_L,R2_L,sep="+"),sep="~"))
    adjustment_R1_formula <- as.formula(paste(R1,C,sep="~"))
    adjustment_R2_formula <- as.formula(paste(R2,C,sep="~"))
  } else if (Q!="") {
    nume_formula <- as.formula(paste(W, paste(R1,R2,Q,R1_Q,R2_Q,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(R1,R2,Q,R1_Q,R2_Q,sep="+"),sep="~"))
  } else if (L!="") {
    nume_formula <- as.formula(paste(W, paste(R1,R2,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(R1,R2,L,R1_L,R2_L,sep="+"),sep="~"))
  } else if (C!="") {
    nume_formula <- as.formula(paste(W, paste(R1,R2,C,R1_C,R2_C,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(R1,R2,C,R1_C,R2_C,sep="+"),sep="~"))
    adjustment_R1_formula <- as.formula(paste(R1,C,sep="~"))
    adjustment_R2_formula <- as.formula(paste(R2,C,sep="~"))
  } else {
    nume_formula <- as.formula(paste(W, paste(R1,R2,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(R1,R2,sep="+"),sep="~"))
  }
  
############################ under construction ############################
  # when common_support is specified (only use cases in the common support), dropped base group (R2) members who don't have target group (R1) counterparts in terms of Q and C.
  # note that logical checks in the beginning have made sure when common_support=T, at least one of Q and C is supplied.
  # the logical checks also make sure that R1 is not a constant.
  if (common_support==T) {
    if (Q!="" & C!="") {
      common_support_formula <- as.formula(paste(R1, paste(Q,C,sep="+"),sep="~"))
    } else if (Q!="") {
      common_support_formula <- as.formula(paste(R1, paste(Q,sep="+"),sep="~"))
    } else {
      common_support_formula <- as.formula(paste(R1, paste(C,sep="+"),sep="~"))
    }
  }

  # common_support_mol <- glm(common_support_formula, family=binomial(link = "logit"), data=data_inuse)
  # common_support_pred <- predict(common_support_mol, newdata = data_inuse[data_inuse[,R2]==1,], type = "response") 
  # min(common_support_pred)
  # sum(common_support_pred<0.01)
  
############################ under construction ############################
  
  # get the numerator and denominator for the intervention weight
  nume_mol <- glm(nume_formula, family=binomial(link = "logit"), data=data_inuse)
  nume_newdata <- data_inuse
  nume_newdata[,R1] <- 1
  nume_newdata[,R2] <- 0
  suppressWarnings( nume_pred <- predict(nume_mol, newdata = nume_newdata, type = "response") )
  nume_pred[data_inuse[,W]==0] <- 1-nume_pred[data_inuse[,W]==0]  # the prediction for people with W=0 should be 1-(P(M=1|covariates))
  
  deno_mol <- glm(deno_formula, family=binomial(link = "logit"), data=data_inuse)
  deno_newdata <- data_inuse
  deno_newdata[,R1] <- 0
  deno_newdata[,R2] <- 1
  suppressWarnings( deno_pred <- predict(deno_mol, newdata = deno_newdata, type = "response") )
  # When R1==1 for every case in the sample, the two lines "nume_newdata[,R1] <- 1" and "deno_newdata[,R1] <- 0" will not affect the results at all. 
  # When R1==1 for every case in the sample, "glm(nume_formula, family=binomial(link = "logit"), data=data_inuse)" will be equivalent to logit with only R2, Q, C, R2:Q, R:C, which is as intended.
  deno_pred[data_inuse[,W]==0] <- 1-deno_pred[data_inuse[,W]==0] 
  
  ## about the two logit models above: it doesn't matter if I specify no intercept, the predicted values will be the exact same.
  # However, if the subsample where either R1 or R2 equal to 1 is selected, the results will be somewhat different. I choose not to use such subsample, as noted above. 
  # Technically, running the models with interactions between R1/R2 and all covariates is like doing subsample modeling within values of R1 and R2.
  # When R1 and R2 are exhaustive, the model only makes use of information from three cells. 
  # But when they are non-exhaustive, the model also uses information from the cell (R1=0, R2=0). 
  # In particular, think of the main effects of Q,L,C, which will be estimated using information from all four cells when R1 and R2 are non-exhaustive. 
  # I think essentially, by using the non-exhaustive sample, I'm assuming some conditional homogeneity of the Q,L,C effects across R1/R2 cells to improve precision.
  
  if (C!="") {

      suppressWarnings( adjustment_R1_pred <- predict(glm(adjustment_R1_formula, family=binomial(link = "logit"), data=data_inuse), type="response") )
      suppressWarnings( adjustment_R2_pred <- predict(glm(adjustment_R2_formula, family=binomial(link = "logit"), data=data_inuse), type="response") )
      # get the truncation threshold, which is the specified quantile of all weights (two adjustment weights and one intervention weight pooled together)
      threshold <- quantile(c(mean(data_inuse[,R1]==1)/adjustment_R1_pred[data_inuse[,R1]==1], 
                  mean(data_inuse[,R2]==1)/adjustment_R2_pred[data_inuse[,R2]==1], 
                  (nume_pred/deno_pred)[data_inuse[,R2]==1]*mean(data_inuse[,R2]==1)/adjustment_R2_pred[data_inuse[,R2]==1]),
                  probs=truncation_threshold)
      # then any cases with weight higher than the threshold are not used.
      # when truncation_threshold==1, all cases are used. 
      retaining_index_R1 <- mean(data_inuse[,R1]==1)/adjustment_R1_pred[data_inuse[,R1]==1]<=threshold
      retaining_index_R2 <- mean(data_inuse[,R2]==1)/adjustment_R2_pred[data_inuse[,R2]==1]<=threshold & (nume_pred/deno_pred)[data_inuse[,R2]==1]*mean(data_inuse[,R2]==1)/adjustment_R2_pred[data_inuse[,R2]==1]<=threshold
      original_R1 <- mean(data_inuse[,Y][data_inuse[,R1]==1][retaining_index_R1]*mean(data_inuse[,R1]==1)/adjustment_R1_pred[data_inuse[,R1]==1][retaining_index_R1])
      original_R2 <- mean(data_inuse[,Y][data_inuse[,R2]==1][retaining_index_R2]*mean(data_inuse[,R2]==1)/adjustment_R2_pred[data_inuse[,R2]==1][retaining_index_R2])
      post_R2 <- mean(data_inuse[,Y][data_inuse[,R2]==1][retaining_index_R2]*(percent/100)*(nume_pred/deno_pred)[data_inuse[,R2]==1][retaining_index_R2]*mean(data_inuse[,R2]==1)/adjustment_R2_pred[data_inuse[,R2]==1][retaining_index_R2])
      # highest and median used weight is stored.
      highest_weight <- max(c(mean(data_inuse[,R1]==1)/adjustment_R1_pred[data_inuse[,R1]==1][retaining_index_R1], 
                            mean(data_inuse[,R2]==1)/adjustment_R2_pred[data_inuse[,R2]==1][retaining_index_R2],
                            (nume_pred/deno_pred)[data_inuse[,R2]==1][retaining_index_R2]*mean(data_inuse[,R2]==1)/adjustment_R2_pred[data_inuse[,R2]==1][retaining_index_R2]))
      median_weight <-  median(c(mean(data_inuse[,R1]==1)/adjustment_R1_pred[data_inuse[,R1]==1][retaining_index_R1], 
                               mean(data_inuse[,R2]==1)/adjustment_R2_pred[data_inuse[,R2]==1][retaining_index_R2],
                               (nume_pred/deno_pred)[data_inuse[,R2]==1][retaining_index_R2]*mean(data_inuse[,R2]==1)/adjustment_R2_pred[data_inuse[,R2]==1][retaining_index_R2]))
  } 
  else {
    
    # When C is absent, original_R1 and original_R2 are not subject to weight truncation
    original_R1 <- mean(data_inuse[,Y][data_inuse[,R1]==1])
    original_R2 <- mean(data_inuse[,Y][data_inuse[,R2]==1])
    # When C is absent, the threshold is simply the specified quantile of the only weight in the model: the intervention weight
    threshold <- quantile((nume_pred/deno_pred)[data_inuse[,R2]==1], probs=truncation_threshold)
    retaining_index <- (nume_pred/deno_pred)[data_inuse[,R2]==1]<=threshold
    post_R2 <- mean(data_inuse[,Y][data_inuse[,R2]==1][retaining_index]*(percent/100)*(nume_pred/deno_pred)[data_inuse[,R2]==1][retaining_index])
    highest_weight <- max((nume_pred/deno_pred)[data_inuse[,R2]==1][retaining_index])
    median_weight <- median((nume_pred/deno_pred)[data_inuse[,R2]==1][retaining_index])
    
  }
  
  if (metric=="difference") {
    original <- original_R1-original_R2
    remaining <- original_R1-post_R2
    reduction <- post_R2-original_R2
  } else if (metric=="risk ratio") {
    original <- original_R1/original_R2
    remaining <- original_R1/post_R2
    reduction <- post_R2/original_R2
  } else if (metric=="odds ratio") {
    original <- (original_R1/(1-original_R1))/(original_R2/(1-original_R2))
    remaining <- (original_R1/(1-original_R1))/(post_R2/(1-post_R2))
    reduction <- (post_R2/(1-post_R2))/(original_R2/(1-original_R2))
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
    nume_pred[data_boot[,W]==0] <- 1-nume_pred[data_boot[,W]==0] 
    
    deno_mol <- glm(deno_formula, family=binomial(link = "logit"), data=data_boot)
    deno_newdata <- data_boot
    deno_newdata[,R1] <- 0
    deno_newdata[,R2] <- 1
    suppressWarnings( deno_pred <- predict(deno_mol, newdata = deno_newdata, type = "response") )
    deno_pred[data_boot[,W]==0] <- 1-deno_pred[data_boot[,W]==0] 
    
    if (C!="") {
      
      suppressWarnings( adjustment_R1_pred <- predict(glm(adjustment_R1_formula, family=binomial(link = "logit"), data=data_boot), type="response") )
      suppressWarnings( adjustment_R2_pred <- predict(glm(adjustment_R2_formula, family=binomial(link = "logit"), data=data_boot), type="response") )
      # get the truncation threshold, which is the specified quantile of all weights (two adjustment weights and one intervention weight pooled together)
      threshold <- quantile(c(mean(data_boot[,R1]==1)/adjustment_R1_pred[data_boot[,R1]==1], 
                              mean(data_boot[,R2]==1)/adjustment_R2_pred[data_boot[,R2]==1], 
                              (nume_pred/deno_pred)[data_boot[,R2]==1]*mean(data_boot[,R2]==1)/adjustment_R2_pred[data_boot[,R2]==1]),
                            probs=truncation_threshold)
      # then any cases with weight higher than the threshold are not used.
      # when truncation_threshold==1, all cases are used. 
      retaining_index_R1 <- mean(data_boot[,R1]==1)/adjustment_R1_pred[data_boot[,R1]==1]<=threshold
      retaining_index_R2 <- mean(data_boot[,R2]==1)/adjustment_R2_pred[data_boot[,R2]==1]<=threshold & (nume_pred/deno_pred)[data_boot[,R2]==1]*mean(data_boot[,R2]==1)/adjustment_R2_pred[data_boot[,R2]==1]<=threshold
      original_R1_boot <- mean(data_boot[,Y][data_boot[,R1]==1][retaining_index_R1]*mean(data_boot[,R1]==1)/adjustment_R1_pred[data_boot[,R1]==1][retaining_index_R1])
      original_R2_boot <- mean(data_boot[,Y][data_boot[,R2]==1][retaining_index_R2]*mean(data_boot[,R2]==1)/adjustment_R2_pred[data_boot[,R2]==1][retaining_index_R2])
      post_R2_boot <- mean(data_boot[,Y][data_boot[,R2]==1][retaining_index_R2]*(percent/100)*(nume_pred/deno_pred)[data_boot[,R2]==1][retaining_index_R2]*mean(data_boot[,R2]==1)/adjustment_R2_pred[data_boot[,R2]==1][retaining_index_R2])
      
    } 
    else {
      
      # When C is absent, original_R1_boot and original_R2_boot are not subject to weight truncation
      original_R1_boot <- mean(data_boot[,Y][data_boot[,R1]==1])
      original_R2_boot <- mean(data_boot[,Y][data_boot[,R2]==1])
      # When C is absent, the threshold is simply the specified quantile of the only weight in the model: the intervention weight
      threshold <- quantile((nume_pred/deno_pred)[data_boot[,R2]==1], probs=truncation_threshold)
      retaining_index <- (nume_pred/deno_pred)[data_boot[,R2]==1]<=threshold
      post_R2_boot <- mean(data_boot[,Y][data_boot[,R2]==1][retaining_index]*(percent/100)*(nume_pred/deno_pred)[data_boot[,R2]==1][retaining_index])
      
    }
    
    if (metric=="difference") {
      original_boot <- original_R1_boot-original_R2_boot
      remaining_boot <- original_R1_boot-post_R2_boot
      reduction_boot <- post_R2_boot-original_R2_boot
    } else if (metric=="risk ratio") {
      original_boot <- original_R1_boot/original_R2_boot
      remaining_boot <- original_R1_boot/post_R2_boot
      reduction_boot <- post_R2_boot/original_R2_boot
    } else if (metric=="odds ratio") {
      original_boot <- (original_R1_boot/(1-original_R1_boot))/(original_R2_boot/(1-original_R2_boot))
      remaining_boot <- (original_R1_boot/(1-original_R1_boot))/(post_R2_boot/(1-post_R2_boot))
      reduction_boot <- (post_R2_boot/(1-post_R2_boot))/(original_R2_boot/(1-original_R2_boot))
    }
    
    boot_original[i] <- original_boot
    boot_remaining[i] <- remaining_boot
    boot_reduction[i] <- reduction_boot
    
  }
  # note that for the bootstrap, instead of using cases not weight-truncated in the original sample, weight truncation is performed independently for each bootstrap sample.
  
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
  
  results <- matrix(nrow=3,ncol=3,dimnames=list(c("Original Disparity","Remaining Disparity","Reduction in Disparity"),c("Estimate","Std. Error",paste((1-alpha)*100,"Conf. Int.",sep="% "))))
  results[1,] <- c(sprintf("%.4f",original),sprintf("%.4f",boot_original_sd),boot_original_ci)
  results[2,] <- c(sprintf("%.4f",remaining),sprintf("%.4f",boot_remaining_sd),boot_remaining_ci)
  results[3,] <- c(sprintf("%.4f",reduction),sprintf("%.4f",boot_reduction_sd),boot_reduction_ci)
  results <- as.data.frame(results)
  
  used_weights_info <- matrix(nrow=1,ncol=2,dimnames=list(NULL,c("Median","Highest")))
  used_weights_info[1,] <- c(format(round(median_weight, 2), nsmall = 2), format(round(highest_weight, 2), nsmall = 2))
  used_weights_info <- as.data.frame(used_weights_info)
  row.names(used_weights_info) <- "Weight"
  
  averages <- matrix(nrow=1,ncol=3,dimnames=list(NULL,c("Original R1","Original R2","Post R2")))
  averages[1,] <- c(format(round(original_R1, 4), nsmall = 4), format(round(original_R2, 4), nsmall = 4), format(round(post_R2, 4), nsmall = 4))  # always reporting 4 digits to the right of the decimal point
  averages <- as.data.frame(averages)
  row.names(averages) <- "Average"
  
  output <- list()
  output[[1]] <- metric
  output[[2]] <- results
  output[[3]] <- averages
  output[[4]] <- used_weights_info

  names(output) <- c("metric","results","averages","used_weights_info")
  
  return(output)
}




