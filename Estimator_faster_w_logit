
equalize <- function(Y, W, R1, R2, Q=NULL, L=NULL, C=NULL, data, percent=100, metric="difference", K=1000, alpha=0.05, truncation_threshold=1, common_support=F, survey_weight=NULL, model="logit") {
  
  options(error = NULL) # so that the error message when stopped doesn't invoke debugging interface. 
  options(scipen=999) # so that the weight output is not printed in scientific notation
  
  # check data
  if(!is.data.frame(data)) stop("data must be a data.frame.",call.=FALSE)
  if(missing(Y)) stop("Y must be provided.")
  if(missing(W)) stop("W must be provided.")
  if(missing(R1)) stop("R1 must be provided.")
  if(missing(R2)) stop("R2 must be provided.")
  
  if (FALSE %in% (c(Y, W, R1, R2, Q, L, C, survey_weight) %in% colnames(data)))  stop("The variables must be in the data.",call.=FALSE)
  data_nom <- na.omit(data[,c(Y, W, R1, R2, Q, L, C, survey_weight)])
  # Even if R1 and R2 do not exhaust the data, the cases in neither R1 nor R2 are still retained.
  # Otherwise, when C is present, the original disparity between R1 and R2 plus the original disparity between R2 and R3 will not equal the original disparity between R1 and R3.
  # Get rid of rows with any infinite value
  data_nom <- data_nom[is.finite(rowSums(data_nom[,unlist(lapply(data_nom, is.numeric))])),]
  # check if there's any character variable
  if (any(sapply(data_nom, is.character))) stop("Character variable shouldn't be used.",call.=FALSE)
  
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
  if(!is.numeric(data_nom[,Y])) stop("Y must be numeric",call.=FALSE) 
  if(any(!unique(data_nom[,Y])%in%c(0,1)) & metric!="difference") stop("Y must only have values of 0 and 1 when the outcome metric is risk ratio or odds ratio.",call.=FALSE) 
  
  # check W
  if(!is.character(W)) stop("W must be a character scalar vector.",call.=FALSE)
  if(length(W)>1) stop("W must be only one variable.",call.=FALSE)
  if(length(unique(data_nom[,W]))!=2) stop("W must be binary.",call.=FALSE)
  if(!is.numeric(data_nom[,W])) stop("W must be numeric.",call.=FALSE)
  if(max(data_nom[,W])!=1 | min(data_nom[,W])!=0) stop("W must only have values of 0 and 1.",call.=FALSE)
  
  # check R1 and R2
  if(!is.character(R1)) stop("R1 must be a character scalar vector.",call.=FALSE)
  if(length(R1)>1) stop("R1 must be only one variable.",call.=FALSE)
  if(length(unique(data_nom[,R1]))!=2 & length(unique(data_nom[,R1]))!=1) stop("R1 must be either binary or constant.",call.=FALSE)
  if(!is.numeric(data_nom[,R1])) stop("R1 must be numeric",call.=FALSE) 
  if((max(data_nom[,R1])!=1 | min(data_nom[,R1])!=0) & length(unique(data_nom[,R1]))==2) stop("R1 must have values of either 0 and 1 or just 1 .",call.=FALSE)
  if(length(unique(data_nom[,R1]))==1) if(unique(data_nom[,R1]!=1)) stop("R1 must have values of either 0 and 1 or just 1 .",call.=FALSE)
  # It's fine for R1 to be constant 1 in the sample (when the target level of W is the sample mean).
  # When R1==1 for all cases, the model will work as intended. 
  
  if(!is.character(R2)) stop("R2 must be a character scalar vector.",call.=FALSE)
  if(length(R2)>1) stop("R2 must be only one variable.",call.=FALSE)
  if(length(unique(data_nom[,R2]))!=2) stop("R2 must be binary.",call.=FALSE)
  if(!is.numeric(data_nom[,R2])) stop("R2 must be numeric",call.=FALSE) 
  if(max(data_nom[,R2])!=1 | min(data_nom[,R2])!=0) stop("R2 must only have values of 0 and 1.",call.=FALSE)
  # note that R1 and R2 are allowed to be overlapped. For example, R1 may be intervened to have the average W across R1 and R2.
  
  # check Q, L, and C
  if(!is.null(Q) & !is.character(Q)) stop("Q must be a character vector.",call.=FALSE)
  if(!is.null(L) & !is.character(L)) stop("L must be a character vector.",call.=FALSE)
  if(!is.null(C) & !is.character(C)) stop("C must be a character vector.",call.=FALSE)
  
  # check survey_weight
  if(!is.null(survey_weight) & !is.character(survey_weight)) stop("survey_weight must be a character scalar vector.",call.=FALSE)
  if(length(survey_weight)>1) stop("survey_weight must be only one variable.",call.=FALSE)
  if(!is.null(survey_weight) & !is.numeric(data_nom[,survey_weight])) stop("survey_weight must be numeric",call.=FALSE) 
  if(!is.null(survey_weight) & isFALSE(require(survey))) stop("The survey package needs to be installed.",call.=FALSE)
  if (is.null(survey_weight)) {
    data_nom$survey_w <- 1
    survey_weight <- "survey_w"}   # if no survey weight is supplied, simply use 1 for everyone
  
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
  
  # check "model"
  if (model!="logit" & model!="rf") stop("Only Logit and Random Forests are supported.",call.=FALSE)
  if (model=="rf" & isFALSE(require(ranger))) stop("The ranger package needs to be installed.",call.=FALSE)
  
  # check common_support argument
  if (!is.logical(common_support)) stop("common_support must be logical.",call.=FALSE)
  if (is.null(Q) & is.null(C) & isTRUE(common_support)) stop("common support restriction should only be applied when Q or C is specified.",call.=FALSE)
  if (sum(data_nom[,R1]==1 & data_nom[,R2]==1)==sum(data_nom[,R1]==1) & isTRUE(common_support)) stop("common support restriction should only be applied when R2 is not strictly a subset of R1.",call.=FALSE)
  # Below I detect nonempty levels in the factor variables in Q/C in the base group that doesn't have corresponding nonempty levels in the target group.
  # If common_support=T, then there's no problem, as those who don't have corresponding levels will not be used in model estimation.
  # Otherwise, non-corresponding levels become a problem, because when glm model will give an error when the model fitted using target group data is used to predict treatment for the base group, as in the numerator of the interventional weight.
  if (isFALSE(common_support)) {
    corresponding_level_indicator <- c(TRUE)
    if (!is.null(Q) & !is.null(C)) {
      if (length( which(sapply(data_nom[, c(Q,C)], is.factor)) )>0) {
        for (q in 1:length( which(sapply(data_nom[, c(Q,C)], is.factor)) ) )  {# iterate over all factor variables
          q_corresonding_level <- 
            levels( droplevels( data_nom[data_nom[,R2]==1, c(Q,C)[ which(sapply(data_nom[, c(Q,C)], is.factor))[q] ] ] ) ) %in%
            levels( droplevels( data_nom[data_nom[,R1]==1, c(Q,C)[ which(sapply(data_nom[, c(Q,C)], is.factor))[q] ] ] ) ) 
          # whether each nonempty level in qth factor in the base group has a corresponding nonempty level in the target group
          corresponding_level_indicator <- c(corresponding_level_indicator,q_corresonding_level)
        }
      }
    } else if (!is.null(Q)) {
      if (length( which(sapply(data_nom[, c(Q)], is.factor)) )>0) {
        for (q in 1:length( which(sapply(data_nom[, c(Q)], is.factor)) ) )  {# iterate over all factor variables
          q_corresonding_level <- 
            levels( droplevels( data_nom[data_nom[,R2]==1, c(Q)[ which(sapply(data_nom[, c(Q)], is.factor))[q] ] ] ) ) %in%
            levels( droplevels( data_nom[data_nom[,R1]==1, c(Q)[ which(sapply(data_nom[, c(Q)], is.factor))[q] ] ] ) ) 
          corresponding_level_indicator <- c(corresponding_level_indicator,q_corresonding_level)
        }
      }
    } else if (!is.null(C)) {
      if (length( which(sapply(data_nom[, c(C)], is.factor)) )>0) {
        for (q in 1:length( which(sapply(data_nom[, c(C)], is.factor)) ) )  {# iterate over all factor variables
          q_corresonding_level <- 
            levels( droplevels( data_nom[data_nom[,R2]==1, c(C)[ which(sapply(data_nom[, c(C)], is.factor))[q] ] ] ) ) %in%
            levels( droplevels( data_nom[data_nom[,R1]==1, c(C)[ which(sapply(data_nom[, c(C)], is.factor))[q] ] ] ) ) 
          corresponding_level_indicator <- c(corresponding_level_indicator,q_corresonding_level)
        }
      }
    }
    if (FALSE %in% corresponding_level_indicator) stop("There is at least one factor variable in Q or C that has a level where there are only base group members but not target group members. You may either change the factor variables into numeric dummies, or specify the common_support option to be TRUE.",call.=FALSE)
  }
  
  
  data_R1 <- data_nom[data_nom[,R1]==1,]
  data_R2 <- data_nom[data_nom[,R2]==1,]
  
  # When common_support is specified (only use cases in the common support), dropped base group (R2) members who don't have target group (R1) counterparts in terms of Q and C.
  # And R1 counterparts are defined as in the convex hull of Q, C vectors of R1 members.
  # Note that logical checks in the beginning have made sure when common_support=T, at least one of Q and C is supplied.
  # The logical checks also make sure that R1 is not a constant.
  # Also, the logical checks make sure that common_support==T only if R2 is not strictly a subset of R1, because in that case the convex hull test will pass all R2 members anyway.
  # Note that the convex hull test only rules out extrapolation, but not interpolation. 
  # Key reference is p.155 in King and Zeng (2006)
  if (common_support==T) {
    
    if (!is.null(Q) & !is.null(C)) {
      left_hand <- data_R1[,c(Q,C)]   # left-hand-side matrix for linear programming
      right_hand <- data_R2[,c(Q,C)]  # right-hand-side matrix
      common_support_formula <- as.formula( paste("~", paste("0",paste(Q,collapse="+"),paste(C,collapse="+"),sep="+"), sep="" ) )
    } else if (!is.null(Q)) {
      left_hand <- data_R1[,c(Q), drop = FALSE]
      right_hand <- data_R2[,c(Q), drop = FALSE]
      common_support_formula <- as.formula( paste("~", paste("0",paste(Q,collapse="+"),sep="+"), sep="" ) )
    } else {
      left_hand <- data_R1[,c(C), drop = FALSE]
      right_hand <- data_R2[,c(C), drop = FALSE]
      common_support_formula <- as.formula( paste("~", paste("0",paste(C,collapse="+"),sep="+"), sep="" ) )
    }
    
    left_hand <- model.matrix(common_support_formula, data = left_hand) # expand factor variables if there is any.
    right_hand <- model.matrix(common_support_formula, data = right_hand) # expand factor variables if there is any.
    
    left_hand_transpose <- rbind(t(left_hand), rep(1, nrow(left_hand))) # transpose and add a row of 1s
    obje <- c(rep(0, nrow(left_hand)))
    dire <- c(rep("=", ncol(left_hand) + 1))
    
    hull.test <- function (i) {
      right_hand_vector <- c(right_hand[i,], 1)  # take the ith row in the matrix and concatenate a 1
      lp.result <- lp(objective.in = obje, const.mat = left_hand_transpose, const.dir = dire, const.rhs = right_hand_vector)  
      # checking whether there is some nonzero coefficients that add to 1 and can combine all persons in A to get person B
      # note that the lp implicitly assumes every variable to be >= 0, which is needed for convex combination.
      return(lp.result$status)  # if 0, then in the hull; if 2, then not. 
    } 
    
    hull_result <- sapply(1:nrow(right_hand), hull.test)
    common_support_indicator <- rep(1, nrow(data_R2))
    common_support_indicator[hull_result==2] <- 0
    
  } else {
    common_support_indicator <- rep(1, nrow(data_R2))
    # if common_support=F, all R2 members are assumed to be in the common support
  }
  n_not_in_common_support <- sum(common_support_indicator==0)
  if (n_not_in_common_support==nrow(data_R2)) stop("Nobody is in the common support",call.=FALSE)
  
  # these with the underline are scalar vector, no longer vectors of multiple variable names
  Q_ <- paste(Q,collapse="+")
  C_ <- paste(C,collapse="+")
  L_ <- paste(L,collapse="+")
  
  # generate formulas needed according to the presence or absence of Q, L, and C.
  if (Q_!="" & L_!="" & C_!="") {
    nume_formula <- as.formula(paste(W, paste(Q_,C_,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(Q_,C_,L_,sep="+"),sep="~"))
    adjustment_R1_formula <- as.formula(paste(R1,C_,sep="~"))
    adjustment_R2_formula <- as.formula(paste(R2,C_,sep="~"))
  } else if (Q_!="" & L_!="") {
    nume_formula <- as.formula(paste(W, paste(Q_,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(Q_,L_,sep="+"),sep="~"))
  } else if (Q_!="" & C_!="") {
    nume_formula <- as.formula(paste(W, paste(Q_,C_,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(Q_,C_,sep="+"),sep="~"))
    adjustment_R1_formula <- as.formula(paste(R1,C_,sep="~"))
    adjustment_R2_formula <- as.formula(paste(R2,C_,sep="~"))
  } else if (L_!="" & C_!="") {
    nume_formula <- as.formula(paste(W, paste(C_,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(C_,L_,sep="+"),sep="~"))
    adjustment_R1_formula <- as.formula(paste(R1,C_,sep="~"))
    adjustment_R2_formula <- as.formula(paste(R2,C_,sep="~"))
  } else if (Q_!="") {
    nume_formula <- as.formula(paste(W, paste(Q_,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(Q_,sep="+"),sep="~"))
  } else if (L_!="") {
    nume_formula <- as.formula(paste(W, 1, sep="~")) # intercept only
    deno_formula <- as.formula(paste(W, paste(L_,sep="+"),sep="~"))
  } else if (C_!="") {
    nume_formula <- as.formula(paste(W, paste(C_,sep="+"),sep="~"))
    deno_formula <- as.formula(paste(W, paste(C_,sep="+"),sep="~"))
    adjustment_R1_formula <- as.formula(paste(R1,C_,sep="~"))
    adjustment_R2_formula <- as.formula(paste(R2,C_,sep="~"))
  } else {
    nume_formula <- as.formula(paste(W, 1 ,sep="~"))
    deno_formula <- as.formula(paste(W, 1 ,sep="~"))
  }
  
  # get the numerator and denominator for the intervention weight
  
  suppressWarnings( nume_mol <- svyglm(nume_formula, family=binomial(link = "logit"), design=svydesign(id=rownames(data_R1), data=data_R1, weights = data_R1[,survey_weight]), control = list(maxit = 15, epsilon = 1e-5)) )
  suppressWarnings( nume_pred <- predict(nume_mol, newdata = data_R2[common_support_indicator==1,], type = "response") )  # only those who are in the common support are used in prediction
  nume_pred[data_R2[common_support_indicator==1,W]==0] <- 1-nume_pred[data_R2[common_support_indicator==1,W]==0]  # the prediction for people with W=0 should be 1-(P(M=1|covariates))
  
  if (model=="logit") {
    suppressWarnings( deno_mol <- svyglm(deno_formula, family=binomial(link = "logit"), design=svydesign(id=rownames(data_R2), data=data_R2, weights = data_R2[,survey_weight]), control = list(maxit = 15, epsilon = 1e-5)) ) # even those not in the common support are used in the model (for efficiency) but they are selected out below.
    suppressWarnings( deno_pred <- predict(deno_mol, newdata = data_R2[common_support_indicator==1,], type = "response") ) 
    deno_pred[data_R2[common_support_indicator==1,W]==0] <- 1-deno_pred[data_R2[common_support_indicator==1,W]==0] 
    
    ## about the two logit models above: it doesn't matter if I specify no intercept, the predicted values will be the exact same.
    # However, if the subsample where either R1 or R2 equal to 1 is retained, the results will be somewhat different. I choose not to use such subsample, as noted above. 
    # Also, all cases in data_nom are used to fit the models, and will have predicted weights for convenience. But only those in data_nom_cs will actually have their outcomes multiplied by the interventional weight.
    
  } else {
    # if model=="rf", random forest is used for the denominator model, but not the numerator model and not the models for the adjustment weight. The latter two are supposed to invoke low-dimensional covariate Q and C. 
    deno_mol <- ranger(deno_formula, data=data_R2, case.weights = data_R2[,survey_weight])
    suppressWarnings( deno_pred <- predict(deno_mol, data = data_R2[common_support_indicator==1,], type = "response")$predictions ) 
    deno_pred[data_R2[common_support_indicator==1,W]==0] <- 1-deno_pred[data_R2[common_support_indicator==1,W]==0] 
  }
  
  # the intervention weight
  intervene_w <- rep(1, nrow(data_R2))
  intervene_w[common_support_indicator==1] <- (percent/100)*(nume_pred/deno_pred)
  # so, intervention weight is estimated for those R2 members who are in the interventional common support. Other R2 members just have interventional weight=1, because their outcomes are kept as are.
  # note that percentage intervention is already incorporated.
  
  if (C_!="") {
    
    suppressWarnings( adjustment_R1_pred <- predict( svyglm(adjustment_R1_formula, family=binomial(link = "logit"), design=svydesign(id=rownames(data_nom), data=data_nom, weights = data_nom[,survey_weight]), control = list(maxit = 15, epsilon = 1e-5)), type = "response") )
    suppressWarnings( adjustment_R2_pred <- predict( svyglm(adjustment_R2_formula, family=binomial(link = "logit"), design=svydesign(id=rownames(data_nom), data=data_nom, weights = data_nom[,survey_weight]), control = list(maxit = 15, epsilon = 1e-5)), type = "response") )
    
    # the R1 adjustment weight
    adjust_R1_w <- mean(data_nom[,R1]==1)/adjustment_R1_pred[data_nom[,R1]==1]
    
    # the R2 adjustment weight
    adjust_R2_w <- mean(data_nom[,R2]==1)/adjustment_R2_pred[data_nom[,R2]==1]
    
    # get the truncation threshold, which is the specified quantile of all weights but not the survey weight (two adjustment weights and one intervention weight pooled together)
    threshold <- quantile(c(adjust_R1_w, 
                            adjust_R2_w,
                            intervene_w*adjust_R2_w),
                          probs=truncation_threshold)
    # then any cases with weight higher than the threshold are not used.
    # when truncation_threshold==1, all cases are used. 
    retaining_indicator_R1 <- adjust_R1_w <= threshold
    retaining_indicator_R2 <- adjust_R2_w <= threshold & intervene_w*adjust_R2_w <= threshold
    
    original_R1 <- weighted.mean(x=data_R1[retaining_indicator_R1,Y], w=(adjust_R1_w*data_R1[,survey_weight])[retaining_indicator_R1])
    original_R2 <- weighted.mean(x=data_R2[retaining_indicator_R2,Y], w=(adjust_R2_w*data_R2[,survey_weight])[retaining_indicator_R2])
    post_R2 <- weighted.mean(data_R2[retaining_indicator_R2,Y], w=(intervene_w*adjust_R2_w*data_R2[,survey_weight])[retaining_indicator_R2])
    
    # highest and median used weight (not accounting for the survey weight) is stored.
    highest_weight <- max(c(adjust_R1_w[retaining_indicator_R1], 
                            adjust_R2_w[retaining_indicator_R2],
                            (intervene_w*adjust_R2_w)[retaining_indicator_R2]))
    median_weight <-  median(c(adjust_R1_w[retaining_indicator_R1], 
                               adjust_R2_w[retaining_indicator_R2],
                               (intervene_w*adjust_R2_w)[retaining_indicator_R2]))
  } else {
    
    # When C is absent, original_R1 and original_R2 are not subject to weight truncation
    original_R1 <- weighted.mean(x=data_R1[,Y], w=data_R1[,survey_weight])
    original_R2 <- weighted.mean(x=data_R2[,Y], w=data_R2[,survey_weight])
    
    # When C is absent, the threshold is simply the specified quantile of the only weight in the model: the intervention weight
    threshold <- quantile(intervene_w, probs=truncation_threshold)
    retaining_indicator <- intervene_w <= threshold
    post_R2 <- weighted.mean(data_R2[retaining_indicator,Y], w=(intervene_w*data_R2[,survey_weight])[retaining_indicator])
    
    highest_weight <- max(intervene_w[retaining_indicator])
    median_weight <- median(intervene_w[retaining_indicator])
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
  
  registerDoParallel()
  
  boot_results <- foreach(i=1:K,.options.RNG=36, .combine='rbind', .packages=c("tidyverse","survey","ranger")) %dorng%  {
    
    indices <- sample(1:nrow(data_nom), nrow(data_nom), replace = TRUE)
    data_boot <- data_nom[indices,]
    
    data_R1 <- data_boot[data_boot[,R1]==1,]
    data_R2 <- data_boot[data_boot[,R2]==1,]
    
    data_R2 <- data_R2 %>% group_by(s_id) %>% mutate(count = n())
    data_R2 <- data_R2[data_R2$count>1, ]
    data_R2 <- as.data.frame(data_R2)
    
    if (common_support==T) {
      
      if (!is.null(Q) & !is.null(C)) {
        left_hand <- data_R1[,c(Q,C)]   # left-hand-side matrix for linear programming
        right_hand <- data_R2[,c(Q,C)]  # right-hand-side matrix
        common_support_formula <- as.formula( paste("~", paste("0",paste(Q,collapse="+"),paste(C,collapse="+"),sep="+"), sep="" ) )
      } else if (!is.null(Q)) {
        left_hand <- data_R1[,c(Q), drop = FALSE]
        right_hand <- data_R2[,c(Q), drop = FALSE]
        common_support_formula <- as.formula( paste("~", paste("0",paste(Q,collapse="+"),sep="+"), sep="" ) )
      } else {
        left_hand <- data_R1[,c(C), drop = FALSE]
        right_hand <- data_R2[,c(C), drop = FALSE]
        common_support_formula <- as.formula( paste("~", paste("0",paste(C,collapse="+"),sep="+"), sep="" ) )
      }
      
      left_hand <- model.matrix(common_support_formula, data = left_hand) 
      right_hand <- model.matrix(common_support_formula, data = right_hand) 
      
      left_hand_transpose <- rbind(t(left_hand), rep(1, nrow(left_hand)))
      obje <- c(rep(0, nrow(left_hand)))
      dire <- c(rep("=", ncol(left_hand) + 1))
      
      hull_result <- sapply(1:nrow(right_hand), hull.test)
      common_support_indicator <- rep(1, nrow(data_R2))
      common_support_indicator[hull_result==2] <- 0
      
    } else {
      common_support_indicator <- rep(1, nrow(data_R2))
    }
    
    suppressWarnings( nume_mol <- svyglm(nume_formula, family=binomial(link = "logit"), design=svydesign(id=rownames(data_R1), data=data_R1, weights = data_R1[,survey_weight]), control = list(maxit = 15, epsilon = 1e-5)) )
    suppressWarnings( nume_pred <- predict(nume_mol, newdata = data_R2[common_support_indicator==1,], type = "response") )  
    nume_pred[data_R2[common_support_indicator==1,W]==0] <- 1-nume_pred[data_R2[common_support_indicator==1,W]==0]  
    
    if (model=="logit") {
      suppressWarnings( deno_mol <- svyglm(deno_formula, family=binomial(link = "logit"), design=svydesign(id=rownames(data_R2), data=data_R2, weights = data_R2[,survey_weight]), control = list(maxit = 15, epsilon = 1e-5)) ) # even those not in the common support are used in the model (for efficiency) but they are selected out below.
      suppressWarnings( deno_pred <- predict(deno_mol, newdata = data_R2[common_support_indicator==1,], type = "response") ) 
      deno_pred[data_R2[common_support_indicator==1,W]==0] <- 1-deno_pred[data_R2[common_support_indicator==1,W]==0] 
    } else {
      deno_mol <- ranger(deno_formula, data=data_R2, case.weights = data_R2[,survey_weight])
      suppressWarnings( deno_pred <- predict(deno_mol, data = data_R2[common_support_indicator==1,], type = "response")$predictions ) 
      deno_pred[data_R2[common_support_indicator==1,W]==0] <- 1-deno_pred[data_R2[common_support_indicator==1,W]==0] 
    }
    
    intervene_w <- rep(1, nrow(data_R2))
    intervene_w[common_support_indicator==1] <- (percent/100)*(nume_pred/deno_pred)
    
    if (!is.null(C)) {
      
      suppressWarnings( adjustment_R1_pred <- predict( svyglm(adjustment_R1_formula, family=binomial(link = "logit"), design=svydesign(id=rownames(data_boot), data=data_boot, weights = data_boot[,survey_weight]), control = list(maxit = 15, epsilon = 1e-5)), type = "response") )
      suppressWarnings( adjustment_R2_pred <- predict( svyglm(adjustment_R2_formula, family=binomial(link = "logit"), design=svydesign(id=rownames(data_boot), data=data_boot, weights = data_boot[,survey_weight]), control = list(maxit = 15, epsilon = 1e-5)), type = "response") )
      
      adjust_R1_w <- mean(data_boot[,R1]==1)/adjustment_R1_pred[data_boot[,R1]==1]
      
      adjust_R2_w <- mean(data_boot[,R2]==1)/adjustment_R2_pred[data_boot[,R2]==1]
      
      threshold <- quantile(c(adjust_R1_w, 
                              adjust_R2_w,
                              intervene_w*adjust_R2_w),
                            probs=truncation_threshold)
      
      retaining_indicator_R1 <- adjust_R1_w <= threshold
      retaining_indicator_R2 <- adjust_R2_w <= threshold & intervene_w*adjust_R2_w <= threshold
      
      original_R1_boot <- weighted.mean(x=data_R1[retaining_indicator_R1,Y], w=(adjust_R1_w*data_R1[,survey_weight])[retaining_indicator_R1])
      original_R2_boot <- weighted.mean(x=data_R2[retaining_indicator_R2,Y], w=(adjust_R2_w*data_R2[,survey_weight])[retaining_indicator_R2])
      post_R2_boot <- weighted.mean(data_R2[retaining_indicator_R2,Y], w=(intervene_w*adjust_R2_w*data_R2[,survey_weight])[retaining_indicator_R2])
    } else {
      
      original_R1_boot <- weighted.mean(x=data_R1[,Y], w=data_R1[,survey_weight])
      original_R2_boot <- weighted.mean(x=data_R2[,Y], w=data_R2[,survey_weight])
      
      threshold <- quantile(intervene_w, probs=truncation_threshold)
      retaining_indicator <- intervene_w <= threshold
      post_R2_boot <- weighted.mean(data_R2[retaining_indicator,Y], w=(intervene_w*data_R2[,survey_weight])[retaining_indicator])
      
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
    
    boot_output <- c(unname(original_boot),unname(remaining_boot),unname(reduction_boot))
  
    return(boot_output)
  }
  # Note that for the bootstrap, instead of using cases not weight-truncated in the original sample, weight truncation is performed independently for each bootstrap sample.
  # The same is applied to common support restriction. 
  
  boot_original <- boot_results[,1]
  boot_remaining <- boot_results[,2]
  boot_reduction <- boot_results[,3]
  
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
  output[[5]] <- nrow(data)-nrow(data_nom)
  output[[6]] <- n_not_in_common_support
  
  names(output) <- c("metric","results","averages","info on weights used (not including survey weights)",
                     "number of rows dropped due to missing values",
                     "number of R2 members who don't have comparable R1 members in Q or C")
  
  return(output)
  
}
