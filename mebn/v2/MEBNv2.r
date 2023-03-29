#
# Functions for constructing
# Mixed Effects Bayesian Network
#
# Jari Turkia
#

mebn.get_mode <- function(v) {
  d <- density(v)
  m <- d$x[which.max(d$y)]
  
  return(m)
}

##################################################

mebn.get_personal_target_guidelines <- function(personal_info,patient_in_dialysis)
{
    # http://webohjekirja.mylabservices.fi/ISLAB/index.php?test=1999
    # P-K, kaikki 3.4 - 4.7
    
    lower_limits <- c(3.4)
    upper_limits <- c(4.7)
    
    if (patient_in_dialysis == FALSE) {
      
      # http://webohjekirja.mylabservices.fi/ISLAB/index.php?test=1431
      
      # 0-5 pv:                   1.55 - 2.65 mmol/l
      # 12-15 V:                  0.95 - 1.75 mmol/l
      # 16 - 17 V:                0.9 - 1.5 mmol/l
      # 4-11 V:                   1.2 - 1.8 mmol/l
      # 6 PV-3 V:                 1.25 - 2.1 mmol/l
      # M 18-49 V:                0.71 - 1.53 mmol/l
      # M yli 50 V:               0.71 - 1.23 mmol/l
      # Naiset yli 18 v:          0.76 - 1.41 mmol/l
      
      # Restricted to adults
      # summary(dialysis$ika)
      #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      # 26.00   54.00   63.00   61.84   70.00   81.00 
      
      # # gender: female = 1, male = 0
      if (personal_info$sukupuoli == 1)
      {
        lower_limits <- cbind(lower_limits, 0.76)
        upper_limits <- cbind(upper_limits, 1.41)
      } else {
        if (personal_info$ika < 50)
        {
          lower_limits <- cbind(lower_limits, 0.71)
          upper_limits <- cbind(upper_limits, 1.53)
        } else {
          lower_limits <- cbind(lower_limits, 0.71)
          upper_limits <- cbind(upper_limits, 1.23)
        }
      }
    } else {

      # Dialysis (phase 5)
      # According to https://www.muma.fi/files/512/munuaisten_vajaatoiminta_ja_kliininen_ravitsemushoito.pdf
      
      lower_limits <- cbind(lower_limits, 0.9)
      upper_limits <- cbind(upper_limits, 1.78)
      
    }
    
    # P-Alb
    # http://webohjekirja.mylabservices.fi/ISLAB/index.php?test=4586
    
    # alle 5 vrk:               28 - 44 g/l
    # 5 vrk-39 V:               36 - 48 g/l
    # 40 - 69 V:                36 - 45 g/l
    # yli 70 V:                 34 - 45 g/l
    
    if (personal_info$ika < 40)
    {
      lower_limits <- cbind(lower_limits, 36)
      upper_limits <- cbind(upper_limits, 48)
    } else if (personal_info$ika < 70) {
      lower_limits <- cbind(lower_limits, 36)
      upper_limits <- cbind(upper_limits, 45)
    } else {
      # yli 70 V
      lower_limits <- cbind(lower_limits, 34)
      upper_limits <- cbind(upper_limits, 45)
    }
    
    personal_limits <- within(list(), {
      lower_limits <- lower_limits
      upper_limits <- upper_limits
    })
    
  return(personal_limits)
}
  
##################################################

mebn.markov_blanket <- function(g, target)
{
  from.parents.e <- E(g)[to(target)]
  to.children.e <- E(g)[from(target)]
  children.v <- V(g)[get.edges(g, to.children.e)[, 2]]
  total.e <- c(from.parents.e, to.children.e, from.spouses.to.children.e)
  delete.e <- E(g)[setdiff(E(g), total.e)]
  g.simple <- delete.edges(g, delete.e)
  g.simple
}

##################################################

# Normalized root mean squared error
mebn.NRMSE <- function(pred_value, true_value, mean_value)
{
  sqrt((pred_value-true_value)^2)/mean_value
}  
  
##################################################

mebn.new_graph_with_randomvariables <- function(datadesc)
{
  library(igraph)  

  # Initialize graph structure
  reaction_graph <- make_empty_graph() 
  
  predictor_columns <- datadesc[datadesc$Order==100,]
  assumed_targets <- datadesc[datadesc$Order==200,]

  # Add nodes to reaction graph for all the random variables  
  reaction_graph <- reaction_graph + vertices(as.vector(assumed_targets$Name), 
                                              label=as.vector(assumed_targets$Name), 
                                              type=200,
                                              color = "#74aaf2", 
                                              size = 1,
                                              shape = "circle")
                                              
                                              # shape = "distbox",
                                              # mean=50,
                                              # l95CI=30,
                                              # u95CI=70,
                                              # scalemin=0,
                                              # scalemax=100)
  
  reaction_graph <- reaction_graph + vertices(as.vector(predictor_columns$Name), 
                                              label=as.vector(predictor_columns$Name), 
                                              type=100,
                                              color = "#3cd164", 
                                              size = 1, 
                                              shape = "circle")
  
                                              # shape = "distbox",
                                              # mean=50,
                                              # l95CI=30,
                                              # u95CI=70,
                                              # scalemin=0,
                                              # scalemax=100)
  

  return(reaction_graph)
}

##################################################

mebn.fully_connected_bipartite_graph <- function(datadesc)
{
  library(igraph)  
  
  # Initialize graph structure
  reaction_graph <- make_empty_graph() 
  
  predictor_columns <- datadesc[datadesc$Order==100,]
  assumed_targets <- datadesc[datadesc$Order==200,]
  
  # Add nodes to reaction graph for all the random variables  
  reaction_graph <- reaction_graph + vertices(as.vector(assumed_targets$Name), 
                                              label=as.vector(assumed_targets$Name), 
                                              type=200,
                                              color = "#74aaf2", 
                                              size = 1,
                                              shape = "circle")
  
  reaction_graph <- reaction_graph + vertices(as.vector(predictor_columns$Name), 
                                              label=as.vector(predictor_columns$Name), 
                                              type=100,
                                              color = "#3cd164", 
                                              size = 1, 
                                              shape = "circle")
  
    
  for (t in as.vector(assumed_targets$Name))
    for (p in as.vector(predictor_columns$Name))
      reaction_graph <- reaction_graph + edge(c(p, t))
  
  return(reaction_graph)
}

##################################################

mebn.posterior <- function(igraph_node)
{

}

##################################################
  
mebn.add_priornodes <- function(datadesc, reaction_graph)
{
  require(igraph)  
  
  # Add prior information from datadesc for random variables
  # - filter prior information for the nutrition predictors
  
  predictor_columns <- datadesc[datadesc$Order==100,]
  
  predictor_names <- as.vector(predictor_columns$Name)
  predprior <- datadesc[(!is.na(datadesc$Lowerbound) | !is.na(datadesc$Upperbound)) & datadesc$Order == 100,]
  predprior$PriorName <- paste0(predprior$Name, "_prior")
  
  # - add prior nodes
  reaction_graph <- reaction_graph + vertices(as.vector(predprior$PriorName), 
                                              label=as.vector(predprior$PriorName), 
                                              color = "#92d3a5", 
                                              size = 1, 
                                              shape = "distbox",
                                              mean=as.vector(predprior$Lowerbound+(predprior$Upperbound-predprior$Lowerbound)/2),
                                              l95CI=as.vector(predprior$Lowerbound),
                                              u95CI=as.vector(predprior$Upperbound),
                                              scalemin=as.vector(predprior$ScaleMin),
                                              scalemax=as.vector(predprior$ScaleMax))
  
  # Connect priors to predictors
  for (p in 1:nrow(predprior))
  {
    print(paste0(predprior[p,]$PriorName, " -> ", predprior[p,]$Name))
    reaction_graph <- reaction_graph + edges(as.vector(predprior[p,]$PriorName), as.vector(predprior[p,]$Name), shape = "arrow", weight = 1)
  }
  
  return(reaction_graph)
}

##################################################

mebn.normalize <- function(x, org_min, org_max) { return ((x - org_min) / (org_max - org_min)) }

##################################################

mebn.renormalize <- function(x, org_min, org_max) { return ((x + org_min) * (org_max - org_min)) }

##################################################

mebn.standardize <- function(x, org_mean, org_sd) 
{ 
  n_scaled <- (x - org_mean) / org_sd
  #n_scaled <- x / org_sd
  return(n_scaled)
}

##################################################

mebn.rescale <- function(x, org_sd) 
{ 
  #n_rescaled <- ((x * org_sd) + org_mean)
  n_rescaled <- x * org_sd
  return(n_rescaled)
}

##################################################

mebn.scale_gaussians <- function(r, data, datadesc, log_transform_ln = FALSE) 
{ 
  # data parameters contains whole dataset 
  # - subset only field(s) in datadesc so that index r matches
  
  data <- data[as.vector(datadesc$Name)]
  s <- data[,r] # data from column r
  
  if (datadesc[r,]$Distribution == "Gaussian")
  {
    if (sum(s) != 0)
    {
      # TODO: Scale also returns the scaling factor. Store it to restore the original scale.
      s <- mebn.standardize(s, mean(s), sd(s))
      #s <- mebn.normalize(s, min(s), max(s))
    }
  } 
  else if (datadesc[r,]$Distribution == "LogNormal") 
  {
    if (sum(s) != 0)
    {
      # Log-transformation before standardization
      if (log_transform_ln == TRUE)
      {
        s <- log(s)
      }
      
      s <- mebn.standardize(s, mean(s), sd(s))
      #s <- mebn.normalize(s, min(s), max(s))
    }
  }
  
  # If column is a factor this will give numeric value instead of factor level
  if (is.factor(s))
  {
    s <- as.numeric(levels(s))[s]
  }
  
  return (s) 
}

##################################################

mebn.set_mv_model_parameters <- function(predictor_columns, target_columns, group_column, inputdata, targetdata = NULL, normalize_values, reg_params = NULL)
{
  ident <- function(x) { return (x) }
  predictors <- inputdata[as.vector(predictor_columns$Name)]
  targets <- inputdata[as.vector(target_columns$Name)]
  
  # Scale if the predictor is Gaussian
  N <- nrow(inputdata)
  Y <- t(targets)
  
  if (normalize_values == TRUE)
  {
    X <- sapply(1:nrow(predictor_columns), mebn.scale_gaussians, data = inputdata, datadesc = predictor_columns)
    
    # append intercept 
    X <- cbind(rep(1,N), X)
  }
  else
  {
    X <- cbind(rep(1,N), apply(predictors, 2, ident))
  }
  
  NH <- 0 
  if (!is.null(targetdata))
  {
    NH <- sum(targetdata[targetdata == 1])
  }
  
  params <- within(list(),
                   {
                     N <- N
                     NH <- NH
                     X <- X
                     v <- ncol(targets)
                     p <- k <- ncol(X)               # all predictors may have random effects
                     vk <- ncol(targets)*ncol(X);
                     Y <- Y
                     Z <- X     
                     J <- length(levels(inputdata[[group_column]]))
                     group <- as.integer(inputdata[[group_column]])
                     holdout <- targetdata
                     offset <- 100
                   })
  
  params <- c(params, reg_params)
  
  return(params)
}

##################################################

mebn.set_mv_two_level_model_parameters <- function(predictor_columns, target_columns, group_column, subject_column, inputdata, targetdata = NULL, normalize_values, reg_params = NULL)
{
  ident <- function(x) { return (x) }
  predictors <- inputdata[as.vector(predictor_columns$Name)]
  targets <- inputdata[as.vector(target_columns$Name)]
  
  # Scale if the predictor is Gaussian
  N <- nrow(inputdata)
  Y <- t(targets)
  
  if (normalize_values == TRUE)
  {
    X <- sapply(1:nrow(predictor_columns), mebn.scale_gaussians, data = inputdata, datadesc = predictor_columns)
    
    # append intercept 
    X <- cbind(rep(1,N), X)
  }
  else
  {
    X <- cbind(rep(1,N), apply(predictors, 2, ident))
  }
  
  NH <- 0 
  if (!is.null(targetdata))
  {
    NH <- sum(targetdata[targetdata == 1])
  }
  
  # lookup table for hierarchical grouping
  
  df.grouping <- inputdata[c(group_column,subject_column)]
  df.grouping$level1 <- as.integer(inputdata[[group_column]])
  df.grouping$level2 <- as.integer(inputdata[[subject_column]])
  df.grouping <- df.grouping[order(df.grouping$level2),]
  df.grouping <- unique(df.grouping)
  
  group_lookup <- df.grouping$level1 

  #n_g <- length(levels(inputdata[[group_column]]))
  #n_s <- length(levels(inputdata[[subject_column]]))
  
  params <- within(list(),
                   {
                     N <- N
                     NH <- NH
                     X <- X
                     v <- ncol(targets)
                     p <- k <- ncol(X)               # all predictors may have random effects
                     Y <- Y
                     Z <- X
                     n_g <- length(unique(inputdata[[group_column]]))
                     group <- as.integer(inputdata[[group_column]])
                     n_s <- length(unique(inputdata[[subject_column]]))
                     subject <- as.integer(inputdata[[subject_column]])
                     group_for_subject <- group_lookup 
                     holdout <- targetdata
                     offset <- 100
                   })

  params <- c(params, reg_params)
  
  return(params)
}

##################################################

mebn.set_mv_csr_parameters <- function(predictor_columns, target_columns, group_column, inputdata, targetdata = NULL, normalize_values, reg_params = NULL)
{
  ident <- function(x) { return (x) }
  predictors <- inputdata[as.vector(predictor_columns$Name)]
  targets <- inputdata[as.vector(target_columns$Name)]
  
  N <- nrow(inputdata)
  Y <- t(targets)
  v <- ncol(targets)
  J <- length(levels(inputdata[[group_column]]))
  group <- as.integer(inputdata[[group_column]])
  
  # Normalize Gaussians to unit scale
  if (normalize_values == TRUE)
  {
    X <- sapply(1:nrow(predictor_columns), mebn.scale_gaussians, data = inputdata, datadesc = predictor_columns)
  } else
  {
    X <- sapply(1:nrow(predictor_columns), mebn.ident, data = inputdata, datadesc = predictor_columns)
  }
  
  # append intercept 
  Z <- cbind(rep(1,N), X)
  
  p <- ncol(X)
  k <- ncol(Z)
  
  # block diagonal X 
  XM <- kronecker(diag(1, v), X)
  
  # extract parameters for CSR sparse matrix
  #XM_csr <- rstan::extract_sparse_parts(XM)
  
  # personal input for one response, block-diagonal regarding to patients
  ZJ <- matrix(0,N,J*k)
  
  for (n in 1:N)
  {   
    a <- (group[n]-1)*k+1
    b <- group[n]*k
    ZJ[n,a:b] = Z[n,1:k];
  }
  
  # block diagonal Z, regarding to responses 
  ZM <- kronecker(diag(1, v), ZJ)
  
  # extract parameters for CSR sparse matrix
  ZM_csr <- rstan::extract_sparse_parts(ZM)
  
  # stacked Y
  Y <- c(targets[,1], targets[,2], targets[,3])
  
  params <- within(list(),
                   {
                     N <- N
                     v <- ncol(targets)
                     p <- k <- ncol(X)               # all predictors may have random effects
                     J <- J
                     Y <- Y
                     X <- XM
                     Z_w_size <- length(ZM_csr$w)
                     Z_u_size <- length(ZM_csr$u)
                     Z_v_size <- length(ZM_csr$v)
                     Z_w <- ZM_csr$w     
                     Z_v <- ZM_csr$v
                     Z_u <- ZM_csr$u
                     offset <- 100
                   })
  params <- c(params, reg_params)
  
  return(params)
}

##################################################

mebn.set_model_parameters <- function(predictor_columns, target_column, group_column, inputdata, targetdata = NULL, normalize_values, reg_params = NULL)
{
  ident <- function(x) { return (x) }
  predictors <- inputdata[as.vector(predictor_columns$Name)]
  
  target_name <- as.vector(target_column$Name)
  
  # Scale if the predictor is Gaussian
  N <- nrow(inputdata)
  
  if (normalize_values == TRUE)
  {
    X <- sapply(1:nrow(predictor_columns), mebn.scale_gaussians, data = inputdata, datadesc = predictor_columns)
    
    # append intercept 
    X <- cbind(rep(1,N), X)
    
    #Y <- scale(inputdata[target_name][,], center = FALSE, scale = TRUE)[,1]
    #Y <- mebn.scale(inputdata[target_name], sd(inputdata[target_name][,]))[,]  
    #Y <- mebn.normalize(inputdata[target_name], min(inputdata[target_name][,]), max(inputdata[target_name][,]))[,]  
    Y <- inputdata[target_name][,]
  }
  else
  {
    X <- cbind(rep(1,N), apply(predictors, 2, ident))
    Y <- inputdata[target_name][,]
  }
  
  NH <- 0 
  if (!is.null(targetdata))
  {
    NH <- sum(holdout_index[holdout_index == 1])
  }
  
  params <- within(list(),
                   {
                     N <- N
                     NH <- NH
                     X <- X
                     p <- k <- ncol(X)               # all predictors may have random effects
                     Y <- Y
                     Z <- X     
                     J <- length(levels(inputdata[[group_column]]))
                     group <- as.integer(inputdata[[group_column]])
                     holdout <- targetdata
                     offset <- 100
                   })
  
  params <- c(params, reg_params)
  
  return(params)
}

##################################################

mebn.set_prediction_parameters <- function(predictor_columns, target_column, inputdata, normalize_values, model_params = NULL)
{
  ident <- function(x) { return (x) }
  
  predictors <- inputdata[as.vector(predictor_columns$Name)]
  target_name <- as.vector(target_column$Name)
  
  # Scale if the predictor is Gaussian
  N <- nrow(inputdata)
  
  if (normalize_values == TRUE)
  {
    X <- sapply(1:nrow(assumedpredictors), mebn.scale_gaussians, data = inputdata, datadesc = assumedpredictors)
    Y <- inputdata[target_name][,]
  }
  else
  {
    X <- cbind(rep(1,N), apply(predictors, 2, ident))
    Y <- inputdata[target_name][,]
  }
  
  # append intercept 
  Z <- cbind(rep(1,N), X)
  
  params <- within(list(),
                   {
                     N <- N
                     X <- X
                     p <- k <- ncol(X)              # all predictors may have random effects, intercept not included
                     Y <- Y
                     Z <- Z     
                   })
  
  params <- c(params, model_params)
  
  return(params)
}

##################################################

mebn.set_model_parameters2 <- function(predictor_columns, target_column, group_column, inputdata, targetdata = NULL, normalize_values, reg_params = NULL)
{
  ident <- function(x) { return (x) }
  predictors <- inputdata[as.vector(predictor_columns$Name)]
  
  target_name <- as.vector(target_column$Name)
  
  prior_sigma <- rep(-1, nrow(predictor_columns))
  dim(prior_sigma) <- nrow(predictor_columns)
  prior_mean <- rep(0, nrow(predictor_columns))
  dim(prior_mean) <- nrow(predictor_columns)
  
  # Set informative priors 
  if (!is.na(predictor_columns$Lowerbound) && !is.na(predictor_columns$Upperbound))
  {  
    prior_sigma <- c(predictor_columns$Upperbound - predictor_columns$Lowerbound)
    prior_mean <- c(predictor_columns$Lowerbound + prior_sigma / 2)
    
    dim(prior_sigma) <- length(prior_sigma)
    dim(prior_mean) <- length(prior_mean)
  }
  
  use_holdout_data <- 0
  N_new <- 0
  X_new <- matrix(NA, nrow=0, ncol=ncol(assumedpredictors))
  Y_new <- c()
  J_new <- 0
  Z_new <- 0
  
  if (!is.null(targetdata))
  {
    use_holdout_data <- 1
    N_new <- nrow(targetdata)
    J_new <- length(levels(targetdata[[group_column]]))
  }
  
  # Scale if the predictor is Gaussian
  N <- nrow(inputdata)
  
  if (normalize_values == TRUE)
  {
    X <- sapply(1:nrow(assumedpredictors), mebn.scale_gaussians, data = inputdata, datadesc = assumedpredictors)
    
    # append intercept 
    X <- cbind(rep(1,N), X)
    
    #Y <- scale(inputdata[target_name][,], center = FALSE, scale = TRUE)[,1]
    #Y <- mebn.scale(inputdata[target_name], sd(inputdata[target_name][,]))[,]  
    #Y <- mebn.normalize(inputdata[target_name], min(inputdata[target_name][,]), max(inputdata[target_name][,]))[,]  
    Y <- inputdata[target_name][,]
    
    # Prepare training data also
    if (use_holdout_data == 1)
    {
      # note: X_new does not need intercept column
      X_new <- sapply(1:nrow(assumedpredictors), mebn.scale_gaussians, data = targetdata, datadesc = assumedpredictors)
      Y_new <- targetdata[target_name][,]
      
      # but Z_new does
      Z_new <- cbind(rep(1,N_new), X_new)
    }
  }
  else
  {
    X <- cbind(rep(1,N), apply(predictors, 2, ident))
    Y <- inputdata[target_name][,]
    
    if (use_holdout_data == 1)
    {
      # note: X_new does not need intercept column
      X_new <- apply(predictors, 2, ident)
      Y_new <- targetdata[target_name][,]
      
      # but Z_new does
      Z_new <- cbind(rep(1,N_new), X_new)
    }    
  }

  params <- within(list(),
                   {
                     N <- N
                     X <- X
                     p <- k <- ncol(X)               # all predictors may have random effects
                     # Mean and variance of Gaussian prior predictors
                     X_prior_sigma <- prior_sigma    # no prior for the intercept 
                     X_prior_mean <- prior_mean      # sigma < 0 means noninformative prior
                     Y <- Y
                     Z <- X     
                     N_new <- N_new
                     X_new <- X_new
                     Z_new <- Z_new
                     Y_new <- Y_new
                     predict_with_holdout <- use_holdout_data
                     J <- length(levels(inputdata[[group_column]]))
                     J_new <- J_new
                     group <- as.integer(inputdata[[group_column]])
                     group_new <- as.integer(targetdata[[group_column]])
                     offset <- 15
                   })
  
  params <- c(params, reg_params)
  
  return(params)
}

##################################################

mebn.localsummary <- function(fit)
{
  #draws <- extract(fit)
  
  #  mean      se_mean         sd          10%         90%     n_eff      Rhat
  ms <- summary(fit, pars=c("beta_Intercept", "beta", "sigma_b", "sigma_e"), probs=c(0.10, 0.90), na.rm = TRUE)

  ModelSummary <- within(list(),
    {
      intmean       <- round(ms$summary[rownames(ms$summary) %in% "beta_Intercept",],5)[1]
      intmean_lCI   <- round(ms$summary[rownames(ms$summary) %in% "beta_Intercept",],5)[4]
      intmean_uCI   <- round(ms$summary[rownames(ms$summary) %in% "beta_Intercept",],5)[5]
      fixef         <- round(ms$summary[startsWith(rownames(ms$summary), "beta["),],5)[,1]
      fixef_lCI     <- round(ms$summary[startsWith(rownames(ms$summary), "beta["),],5)[,4]
      fixef_uCI     <- round(ms$summary[startsWith(rownames(ms$summary), "beta["),],5)[,5]
      ranef_sd      <- round(ms$summary[startsWith(rownames(ms$summary), "sigma_b["),],5)[,1]
      ranef_sd_lCI  <- round(ms$summary[startsWith(rownames(ms$summary), "sigma_b["),],5)[,4]
      ranef_sd_uCI  <- round(ms$summary[startsWith(rownames(ms$summary), "sigma_b["),],5)[,5]
      std_error     <- round(ms$summary[rownames(ms$summary) %in% "sigma_e",],5)[1]
      std_error_lCI <- round(ms$summary[rownames(ms$summary) %in% "sigma_e",],5)[4]
      std_error_uCI <- round(ms$summary[rownames(ms$summary) %in% "sigma_e",],5)[5]
    })
  
  # Create matrix D
  #sdM <- diag(ModelSummary$ranef_sd)
  #ModelSummary$D <- sdM %*% ModelSummary$C %*% t(sdM)
  
  return(ModelSummary)
}

##################################################

mebn.localsummary_from_multivariate <- function(fit, targetindex)
{
  #draws <- extract(fit)
  
  #  mean      se_mean         sd          10%         90%     n_eff      Rhat
  ms <- summary(fit, pars=c(paste0("beta_Intercept"), paste0("beta"), paste0("sigma_b")), probs=c(0.10, 0.90), na.rm = TRUE)
  #ms <- summary(fit, pars=c(paste0("beta_Intercept"), paste0("beta"), paste0("sigma_b"), paste0("sigma_e")), probs=c(0.10, 0.90), na.rm = TRUE)
  
  ModelSummary <- within(list(),
                         {
                           intmean       <- round(ms$summary[rownames(ms$summary) %in% paste0("beta_Intercept[",targetindex,"]"),],5)[1]
                           intmean_lCI   <- round(ms$summary[rownames(ms$summary) %in% paste0("beta_Intercept[",targetindex,"]"),],5)[4]
                           intmean_uCI   <- round(ms$summary[rownames(ms$summary) %in% paste0("beta_Intercept[",targetindex,"]"),],5)[5]
                           fixef         <- round(ms$summary[startsWith(rownames(ms$summary), paste0("beta[",targetindex,",")),],5)[,1]
                           fixef_lCI     <- round(ms$summary[startsWith(rownames(ms$summary), paste0("beta[",targetindex,",")),],5)[,4]
                           fixef_uCI     <- round(ms$summary[startsWith(rownames(ms$summary), paste0("beta[",targetindex,",")),],5)[,5]
                           ranef_sd      <- round(ms$summary[startsWith(rownames(ms$summary), paste0("sigma_b[",targetindex,",")),],5)[,1]
                           ranef_sd_lCI  <- round(ms$summary[startsWith(rownames(ms$summary), paste0("sigma_b[",targetindex,",")),],5)[,4]
                           ranef_sd_uCI  <- round(ms$summary[startsWith(rownames(ms$summary), paste0("sigma_b[",targetindex,",")),],5)[,5]
                         })

  # This is missing from the cross-correlation model
  
  #std_error     <- round(ms$summary[rownames(ms$summary) %in% paste0("sigma_e[",targetindex,"]"),],5)[1]
  #std_error_lCI <- round(ms$summary[rownames(ms$summary) %in% paste0("sigma_e[",targetindex,"]"),],5)[4]
  #std_error_uCI <- round(ms$summary[rownames(ms$summary) %in% paste0("sigma_e[",targetindex,"]"),],5)[5]
  
  # Create matrix D
  #sdM <- diag(ModelSummary$ranef_sd)
  #ModelSummary$D <- sdM %*% ModelSummary$C %*% t(sdM)
  
  return(ModelSummary)
}

##################################################

mebn.localsummary_from_two_level_multivariate <- function(fit, targetindex)
{
  #  mean      se_mean         sd          10%         90%     n_eff      Rhat
  ms <- summary(fit, pars=c(paste0("beta_Intercept"), paste0("beta"), paste0("sigma_b_g"), paste0("sigma_b_s")), probs=c(0.5, 0.95), na.rm = TRUE)
  
  # for mode
  int_parname <- paste0("beta_Intercept[",targetindex,"]")
  int_draws <- rstan::extract(fit, pars = int_parname)
  fixef_params <- grep(paste0("beta\\[",targetindex,","), rownames(ms$summary), value = TRUE)
  fixef_draws <- rstan::extract(fit, pars = fixef_params)
  
  ModelSummary <- within(list(),
                         {
                           intmean       <- round(ms$summary[rownames(ms$summary) %in% paste0("beta_Intercept[",targetindex,"]"),],5)[1]
                           intmode       <- round(mebn.get_mode(int_draws[[int_parname]]), 5)
                           intmean_lCI   <- round(ms$summary[rownames(ms$summary) %in% paste0("beta_Intercept[",targetindex,"]"),],5)[4]
                           intmean_uCI   <- round(ms$summary[rownames(ms$summary) %in% paste0("beta_Intercept[",targetindex,"]"),],5)[5]
                           fixef         <- round(ms$summary[startsWith(rownames(ms$summary), paste0("beta[",targetindex,",")),],5)[,1]
                           fixef_mode    <- round(unlist(lapply(fixef_draws, mebn.get_mode)), 5)
                           fixef_lCI     <- round(ms$summary[startsWith(rownames(ms$summary), paste0("beta[",targetindex,",")),],5)[,4]
                           fixef_uCI     <- round(ms$summary[startsWith(rownames(ms$summary), paste0("beta[",targetindex,",")),],5)[,5]
                           ranef1_sd      <- round(ms$summary[startsWith(rownames(ms$summary), paste0("sigma_b_g[",targetindex,",")),],5)[,1]
                           ranef1_sd_lCI  <- round(ms$summary[startsWith(rownames(ms$summary), paste0("sigma_b_g[",targetindex,",")),],5)[,4]
                           ranef1_sd_uCI  <- round(ms$summary[startsWith(rownames(ms$summary), paste0("sigma_b_g[",targetindex,",")),],5)[,5]
                           ranef2_sd      <- round(ms$summary[startsWith(rownames(ms$summary), paste0("sigma_b_s[",targetindex,",")),],5)[,1]
                           ranef2_sd_lCI  <- round(ms$summary[startsWith(rownames(ms$summary), paste0("sigma_b_s[",targetindex,",")),],5)[,4]
                           ranef2_sd_uCI  <- round(ms$summary[startsWith(rownames(ms$summary), paste0("sigma_b_s[",targetindex,",")),],5)[,5]
                         })
  
  # This is missing from the cross-correlation model
  
  #std_error     <- round(ms$summary[rownames(ms$summary) %in% paste0("sigma_e[",targetindex,"]"),],5)[1]
  #std_error_lCI <- round(ms$summary[rownames(ms$summary) %in% paste0("sigma_e[",targetindex,"]"),],5)[4]
  #std_error_uCI <- round(ms$summary[rownames(ms$summary) %in% paste0("sigma_e[",targetindex,"]"),],5)[5]
  
  # Create matrix D
  #sdM <- diag(ModelSummary$ranef_sd)
  #ModelSummary$D <- sdM %*% ModelSummary$C %*% t(sdM)
  
  return(ModelSummary)
}

##################################################

mebn.personal_effects <- function(fit, person_id)
{
  #  mean      se_mean         sd          10%         90%     n_eff      Rhat
  ms <- rstan::summary(fit, pars=c("b", "personal_effect"), probs=c(0.10, 0.90), na.rm = TRUE)
  
  ModelSummary <- within(list(),
                         {
                           b       <- round(ms$summary[startsWith(rownames(ms$summary),paste0("b[", person_id, ",")),1], 5)
                           b_lCI   <- round(ms$summary[startsWith(rownames(ms$summary),paste0("b[", person_id, ",")),4], 5)
                           b_uCI   <- round(ms$summary[startsWith(rownames(ms$summary),paste0("b[", person_id, ",")),5], 5)
                           personal_effect         <- round(ms$summary[startsWith(rownames(ms$summary), paste0("personal_effect[", person_id, ",")),1], 5)
                           personal_effect_lCI     <- round(ms$summary[startsWith(rownames(ms$summary), paste0("personal_effect[", person_id, ",")),4], 5)
                           personal_effect_uCI     <- round(ms$summary[startsWith(rownames(ms$summary), paste0("personal_effect[", person_id, ",")),5], 5)
                         })
  
  return(ModelSummary)
}

##################################################

mebn.personal_effects_from_multivariate <- function(fit, targetindex, person_id)
{
  #  mean      se_mean         sd          10%         90%     n_eff      Rhat
  ms <- rstan::summary(fit, pars=c(paste0("b"), paste0("personal_effect")), probs=c(0.10, 0.90), na.rm = TRUE)
  
  ModelSummary <- within(list(),
                         {
                           b       <- round(ms$summary[startsWith(rownames(ms$summary),paste0("b[",person_id,",", targetindex, ",")),1], 5)
                           b_lCI   <- round(ms$summary[startsWith(rownames(ms$summary),paste0("b[",person_id,",", targetindex, ",")),4], 5)
                           b_uCI   <- round(ms$summary[startsWith(rownames(ms$summary),paste0("b[",person_id,",", targetindex, ",")),5], 5)
                           personal_effect         <- round(ms$summary[startsWith(rownames(ms$summary), paste0("personal_effect[",person_id,",", targetindex, ",")),1], 5)
                           personal_effect_lCI     <- round(ms$summary[startsWith(rownames(ms$summary), paste0("personal_effect[",person_id,",", targetindex, ",")),4], 5)
                           personal_effect_uCI     <- round(ms$summary[startsWith(rownames(ms$summary), paste0("personal_effect[",person_id,",", targetindex, ",")),5], 5)
                         })
  
  return(ModelSummary)
}

##################################################

mebn.adjusted_effects_from_multivariate <- function(fit, targetindex, person_id, group_id)
{
  # Get full hierarchy of effects until the given level
  
  #  mean      se_mean         sd          10%         90%     n_eff      Rhat
  ms <- rstan::summary(fit, pars=c(paste0("g"), paste0("group_effect")), probs=c(0.05, 0.95), na.rm = TRUE)

  # for mode
  g_params <- grep(paste0("g\\[",group_id,",", targetindex, ","), rownames(ms$summary), value = TRUE)
  g_draws <- rstan::extract(fit, pars = g_params)
  group_effect_params <- grep(paste0("group_effect\\[",group_id,",", targetindex, ","), rownames(ms$summary), value = TRUE)
  group_effect_draws <- rstan::extract(fit, pars = group_effect_params)
  
  ModelSummary1 <- within(list(),
                          {
                            g       <- round(ms$summary[startsWith(rownames(ms$summary),paste0("g[",group_id,",", targetindex, ",")),1], 5)
                            g_mode  <- round(unlist(lapply(g_draws, mebn.get_mode)),5)
                            g_lCI   <- round(ms$summary[startsWith(rownames(ms$summary),paste0("g[",group_id,",", targetindex, ",")),4], 5)
                            g_uCI   <- round(ms$summary[startsWith(rownames(ms$summary),paste0("g[",group_id,",", targetindex, ",")),5], 5)
                            group_effect         <- round(ms$summary[startsWith(rownames(ms$summary), paste0("group_effect[",group_id,",", targetindex, ",")),1], 5)
                            group_effect_mode    <- round(unlist(lapply(group_effect_draws, mebn.get_mode)), 5)
                            group_effect_lCI     <- round(ms$summary[startsWith(rownames(ms$summary), paste0("group_effect[",group_id,",", targetindex, ",")),4], 5)
                            group_effect_uCI     <- round(ms$summary[startsWith(rownames(ms$summary), paste0("group_effect[",group_id,",", targetindex, ",")),5], 5)
                          })
  
  #  mean      se_mean         sd          10%         90%     n_eff      Rhat
  ms <- rstan::summary(fit, pars=c(paste0("b"), paste0("personal_effect")), probs=c(0.05, 0.95), na.rm = TRUE)

  # for mode
  b_params <- grep(paste0("b\\[",person_id,",", targetindex, ","), rownames(ms$summary), value = TRUE)
  b_draws <- rstan::extract(fit, pars = b_params)
  personal_effect_params <- grep(paste0("personal_effect\\[",person_id,",", targetindex, ","), rownames(ms$summary), value = TRUE)
  personal_effect_draws <- rstan::extract(fit, pars = personal_effect_params)
  
  ModelSummary2 <- within(list(),
                          {
                            b       <- round(ms$summary[startsWith(rownames(ms$summary),paste0("b[",person_id,",", targetindex, ",")),1], 5)
                            b_mode  <- round(unlist(lapply(b_draws, mebn.get_mode)), 5)
                            b_lCI   <- round(ms$summary[startsWith(rownames(ms$summary),paste0("b[",person_id,",", targetindex, ",")),4], 5)
                            b_uCI   <- round(ms$summary[startsWith(rownames(ms$summary),paste0("b[",person_id,",", targetindex, ",")),5], 5)
                            personal_effect         <- round(ms$summary[startsWith(rownames(ms$summary), paste0("personal_effect[",person_id,",", targetindex, ",")),1], 5)
                            personal_effect_mode    <- round(unlist(lapply(personal_effect_draws, mebn.get_mode)), 5)
                            personal_effect_lCI     <- round(ms$summary[startsWith(rownames(ms$summary), paste0("personal_effect[",person_id,",", targetindex, ",")),4], 5)
                            personal_effect_uCI     <- round(ms$summary[startsWith(rownames(ms$summary), paste0("personal_effect[",person_id,",", targetindex, ",")),5], 5)
                          })
  
  ModelSummary <- append(ModelSummary1,ModelSummary2)  
  
  return(ModelSummary)
}

##################################################

mebn.personal_effects_from_multivariate_namedtargets <- function(fit, targetname, person_id)
{
  #  mean      se_mean         sd          10%         90%     n_eff      Rhat
  ms <- summary(fit, pars=c(paste0("b_",targetname), paste0("personal_effect_",targetname)), probs=c(0.10, 0.90), na.rm = TRUE)
  
  ModelSummary <- within(list(),
                         {
                           b       <- round(ms$summary[startsWith(rownames(ms$summary),paste0("b_",targetname,"[", person_id, ",")),1], 5)
                           b_lCI   <- round(ms$summary[startsWith(rownames(ms$summary),paste0("b_",targetname,"[", person_id, ",")),4], 5)
                           b_uCI   <- round(ms$summary[startsWith(rownames(ms$summary),paste0("b_",targetname,"[", person_id, ",")),5], 5)
                           personal_effect         <- round(ms$summary[startsWith(rownames(ms$summary), paste0("personal_effect_",targetname,"[", person_id, ",")),1], 5)
                           personal_effect_lCI     <- round(ms$summary[startsWith(rownames(ms$summary), paste0("personal_effect_",targetname,"[", person_id, ",")),4], 5)
                           personal_effect_uCI     <- round(ms$summary[startsWith(rownames(ms$summary), paste0("personal_effect_",targetname,"[", person_id, ",")),5], 5)
                         })
  
  return(ModelSummary)
}

##################################################

mebn.get_localfit <- function(target_name, local_model_cache = "models", mem_cache = FALSE)
{
  modelcache <- paste0(local_model_cache, "/", target_name, "_blmm", ".rds")
  localfit <- NULL
  
  memcache <- gsub("/", "_", modelcache)  
  memcache <- gsub(".rds", "", memcache)  
  
  if (mem_cache == TRUE & exists(memcache))
  {
    localfit <- get(memcache)
    #print(paste0("Using cached ", memcache))
  }
  else if (file.exists(modelcache))
  {
    localfit <- readRDS(modelcache)
    #print(paste0("Loading file ", modelcache))
    
    # Set memcache to Global Environment
    if (mem_cache == TRUE) assign(memcache, localfit, 1)

    #print(paste0("Setting memory cache: ", memcache))
  }
  
  return(localfit)
}

##################################################

mebn.get_parents <- function(g, nodename)
{
  vindex <- as.numeric(V(g)[nodename])
  parents <- neighbors(g, vindex, mode = c("in"))
  
  return(parents)
}

##################################################

mebn.get_parents_with_type <- function(g, nodename, type)
{
  vindex <- as.numeric(V(g)[nodename])
  parents <- neighbors(g, vindex, mode = c("in"))
  
  # filter by type
  parents <- parents[parents$type==type]
  
  return(parents)
}

##################################################

mebn.LOO_comparison <- function(target_variables, graphdir1, graphdir2)
{
  library(loo)
  
  comparison<-data.frame(matrix(nrow=nrow(target_variables), ncol=3))
  colnames(comparison) <- c("distribution", graphdir1, graphdir2)
  
  n <- 1
  for (targetname in target_variables$Name)
  {
    # Get models to compare
    m1 <- mebn.get_localfit(paste0(graphdir1, "/", targetname))
    m2 <- mebn.get_localfit(paste0(graphdir2, "/", targetname))
    
    # Statistics for model 1
    m1_loglik <- extract_log_lik(m1, merge_chains = FALSE)
    
    if (exists("m1_rel_n_eff")) remove(m1_rel_n_eff)
    if (exists("m1_loo")) remove(m1_loo)
    
    if (!any(is.na(exp(m1_loglik))))
    {
      m1_rel_n_eff <- relative_eff(exp(m1_loglik))
      suppressWarnings(m1_loo <- loo(m1_loglik, r_eff = m1_rel_n_eff, cores = 4))
    }
    
    # Statistics for model 2
    m2_loglik <- extract_log_lik(m2, merge_chains = FALSE)
    
    if (exists("m2_rel_n_eff")) remove(m2_rel_n_eff)
    if (exists("m2_loo")) remove(m2_loo)
    
    if (!any(is.na(exp(m2_loglik))))
    {
      m2_rel_n_eff <- relative_eff(exp(m2_loglik))
      suppressWarnings(m2_loo <- loo(m2_loglik, r_eff = m2_rel_n_eff, cores = 4))
    }
    
    c1 <- targetname
    
    if (exists("m1_loo"))
      c2 <- paste0(round(m1_loo$estimates[1,1], 3), " (", round(m1_loo$estimates[1,2], 3), ")")
    else
      c2 <- "NA"
    
    if (exists("m2_loo"))
      c3 <- paste0(round(m2_loo$estimates[1,1], 3), " (", round(m2_loo$estimates[1,2], 3), ")")
    else
      c3 <- "NA"
    
    comparison[n,1] <- c1
    comparison[n,2] <- c2
    comparison[n,3] <- c3
    n <- n + 1
  }
  
  return(comparison)
}

##################################################

mebn.LOO_comparison3 <- function(target_variables, graphdir1, graphdir2, graphdir3, N)
{
  library(loo)
  
  comparison<-data.frame(matrix(nrow=nrow(target_variables), ncol=7))
  colnames(comparison) <- c("distribution", graphdir1, "bad_k", graphdir2, "bad_k", graphdir3, "bad_k")
  
  n <- 1
  for (targetname in target_variables$Name)
  {
    # Get models to compare
    m1 <- mebn.get_localfit(paste0(graphdir1, "/", targetname))
    m2 <- mebn.get_localfit(paste0(graphdir2, "/", targetname))
    m3 <- mebn.get_localfit(paste0(graphdir3, "/", targetname))
    
    # Statistics for model 1
    m1_loglik <- extract_log_lik(m1, merge_chains = FALSE)
    
    if (exists("m1_rel_n_eff")) remove(m1_rel_n_eff)
    if (exists("m1_loo")) remove(m1_loo)
    
    if (!any(is.na(exp(m1_loglik))))
    {
      m1_rel_n_eff <- relative_eff(exp(m1_loglik))
      suppressWarnings(m1_loo <- loo(m1_loglik, r_eff = m1_rel_n_eff, cores = 4))
    }
    
    # Statistics for model 2
    m2_loglik <- extract_log_lik(m2, merge_chains = FALSE)
    
    if (exists("m2_rel_n_eff")) remove(m2_rel_n_eff)
    if (exists("m2_loo")) remove(m2_loo)
    
    if (!any(is.na(exp(m2_loglik))))
    {
      m2_rel_n_eff <- relative_eff(exp(m2_loglik))
      suppressWarnings(m2_loo <- loo(m2_loglik, r_eff = m2_rel_n_eff, cores = 4))
    }

    # Statistics for model 3
    m3_loglik <- extract_log_lik(m3, merge_chains = FALSE)
    
    if (exists("m3_rel_n_eff")) remove(m3_rel_n_eff)
    if (exists("m3_loo")) remove(m3_loo)
    
    if (!any(is.na(exp(m3_loglik))))
    {
      m3_rel_n_eff <- relative_eff(exp(m3_loglik))
      suppressWarnings(m3_loo <- loo(m3_loglik, r_eff = m3_rel_n_eff, cores = 4))
    }
    
    comparison[n,1] <- targetname
    
    if (exists("m1_loo")) {
      comparison[n,2] <- m1_loo$estimates[1,2]
      
      high_pareto_k <- length(pareto_k_ids(m1_loo, threshold = 0.7))
      high_k_percent <- high_pareto_k/N*100
      
      comparison[n,3] <- high_k_percent
    }
    else
    {
      comparison[n,2] <- "NA"
      comparison[n,3] <- "NA"
    }
    
    if (exists("m2_loo"))
    {
      comparison[n,4] <- m2_loo$estimates[1,2]
      
      high_pareto_k <- length(pareto_k_ids(m2_loo, threshold = 0.7))
      high_k_percent <- high_pareto_k/N*100
      
      comparison[n,5] <- high_k_percent
    }
    else
    {
      comparison[n,4] <- "NA"
      comparison[n,5] <- "NA"
    }

    if (exists("m3_loo"))
    {
      comparison[n,6] <- m3_loo$estimates[1,2]
      
      high_pareto_k <- length(pareto_k_ids(m3_loo, threshold = 0.7))
      high_k_percent <- high_pareto_k/N*100
      
      comparison[n,7] <- high_k_percent
    }
    else
    {
      comparison[n,6] <- "NA"
      comparison[n,7] <- "NA"
    }

    n <- n + 1
  }
  
  return(comparison)
}

##################################################
mebn.linpred <- function(X, beta, Z, b, g_alpha)
{
  
  mu <- beta %*% X + b %*% Z
  
  # ?? 
  return (dgamma(g_alpha, g_alpha(mu)))
}

##################################################

mebn.AR_comparison <- function(target_variables, graphdir)
{
  library(rstan)
  
  ar_table<-data.frame(matrix(nrow=nrow(target_variables), ncol=4))
  colnames(ar_table) <- c("distribution", "AR(1)", "CI-10%", "CI-90%")
  
  n <- 1
  for (targetname in target_variables$Name)
  {
    m1 <- mebn.get_localfit(paste0(graphdir, "/", targetname))
    
    c2 <- "NA"
    c3 <- "NA"
    c4 <- "NA"
    
    s1 <- summary(m1, pars=c("ar1"), probs=c(0.10, 0.90), na.rm = TRUE)
    #m1_extract <- extract(m1, pars = c("ar1"))
  
    ar_table[n,1] <- targetname
    ar_table[n,2] <- round(s1$summary[1],5)
    ar_table[n,3] <- round(s1$summary[4],5)
    ar_table[n,4] <- round(s1$summary[5],5)
    n <- n + 1
  }
  
  return(ar_table)
}

##################################################

mebn.target_dens_overlays <- function(localfit_directory, target_variables, dataset)
{
  library(rstan)
  library(bayesplot)
  library(ggplot2)
  
  # TODO: Check if dir exists (localfit_directory)
  
  color_scheme_set("purple")
  bayesplot_theme_set(theme_default())
  
  dens_plots <- list()
  i <- 1
  
  for (targetname in target_variables$Name)
  {
    target_blmm <- mebn.get_localfit(paste0(localfit_directory,targetname))
    true_value <- as.vector(dataset[,targetname])
    
    posterior <- rstan::extract(target_blmm, pars = c("Y_rep"))
    posterior_y_50 <- posterior$Y_rep[1:100,]
    
    scalemin <- as.numeric(as.character(target_variables[target_variables$Name == targetname,]$ScaleMin))
    scalemax <- as.numeric(as.character(target_variables[target_variables$Name == targetname,]$ScaleMax))
    
    dens_plots[[i]] <- ppc_dens_overlay(true_value, posterior_y_50) + 
      coord_cartesian(xlim = c(scalemin,scalemax)) +
      ggtitle(targetname)
    
    i <- i + 1
  }
  
  bayesplot_grid(plots = dens_plots, legends = FALSE)
}

##################################################

mebn.multivariate_dens_overlays <- function(modelfile_path, target_variables, dataset)
{
  library(rstan)
  library(bayesplot)
  library(ggplot2)
  
  # TODO: Check if dir exists (localfit_directory)
  
  color_scheme_set("purple")
  bayesplot_theme_set(theme_default())
  
  multivariate_model_name <- paste0(target_variables$Name, collapse = "_")
  
  target_blmm <- mebn.get_localfit(multivariate_model_name, modelfile_path)
  
  dens_plots <- list()
  i <- 1

  for (targetname in target_variables$Name)
  {
    print(targetname)
    true_value <- as.vector(dataset[,targetname])
    
    posterior <- rstan::extract(target_blmm, pars = c("Y_rep"))
    posterior_y_50 <- posterior$Y_rep[3800:3900,i,]
    #posterior_y_50 <- posterior$Y_rep[,i,]
    
    scalemin <- as.numeric(as.character(target_variables[target_variables$Name == targetname,]$ScaleMin))
    scalemax <- as.numeric(as.character(target_variables[target_variables$Name == targetname,]$ScaleMax))
    
    dens_plots[[i]] <- ppc_dens_overlay(true_value, posterior_y_50) + 
      coord_cartesian(xlim = c(scalemin,scalemax)) +
      ggtitle(targetname)
    
    i <- i + 1
  }
  
  bayesplot_grid(plots = dens_plots, legends = FALSE)
}


##################################################

mebn.evaluate_predictions <- function(localfit_directory, target_variables, dataset)
{
  library(rstan)
  library(bayesplot)
  library(ggplot2)
  
  color_scheme_set("purple")

  rec_plots <- list()
  i <- 1
  
  for (targetname in target_variables$Name)
  {
    target_blmm <- mebn.get_localfit(paste0(localfit_directory,targetname))
    true_value <- as.vector(dataset[,targetname])
    
    draws <- as.matrix(target_blmm)
    
    obs <- length(true_value)
    
    p1 <- match("Y_pred[1]",colnames(draws))
    p2 <- match(paste0("Y_pred[", obs, "]"),colnames(draws))
    
    # https://mc-stan.org/bayesplot/reference/MCMC-recover.html
    rec_plots[[i]] <- mcmc_recover_hist(draws[, p1:p2], true_value[1:obs])
    i <- i + 1
  }
  
  bayesplot_grid(plots = rec_plots, legends = FALSE)
}

##################################################

mebn.sampling <- function(inputdata, targetdata, predictor_columns, target_column, group_column, local_model_cache = "models", stan_mode_file = "BLMM.stan", normalize_values = TRUE, reg_params = NULL)
{
  require(rstan)
  
  # Run Stan parallel on multiple cores
  rstan_options (auto_write=TRUE)
  options (mc.cores=parallel::detectCores ()) 
  
  target_name <- as.vector(target_column$Name)

  localfit <- mebn.get_localfit(target_name, local_model_cache)
  
  # Use cached model if it exists
  if (is.null(localfit))
  {
    stanDat <- mebn.set_model_parameters(predictor_columns, target_column, group_column, inputdata, targetdata, normalize_values, reg_params)

    localfit <- stan(file=stan_mode_file, data=stanDat, warmup = 1000, iter=2000, chains=1, init=0, control = list(adapt_delta = 0.80, max_treedepth = 12))
    
    print("Sampling done. Saving..")
    
    modelcache <- paste0(local_model_cache, "/", target_name, "_blmm", ".rds")
    
    saveRDS(localfit, file=modelcache)
    print("Done.")
  }
  
  return(localfit)
}  

##################################################

mebn.multivariate_sampling_fixed <- function(inputdata, targetdata, predictor_columns, target_columns, group_column, local_model_cache = "models", stan_mode_file = "BLMM.stan", normalize_values = TRUE, reg_params = NULL)
{
  require(rstan)
  
  # Run Stan parallel on multiple cores
  rstan_options (auto_write=TRUE)
  options (mc.cores=parallel::detectCores ()) 
  
  target_name <- paste0(target_columns$Name, collapse = "_")
  
  localfit <- mebn.get_localfit(target_name, local_model_cache)
  
  # Use cached model if it exists
  if (is.null(localfit))
  {
    # This functions is fixed to responses named "pk" and "fppi" - it will be generalized later
    stanDat <- mebn.set_mv_fixed_model_parameters(predictor_columns, target_columns, group_column, inputdata, targetdata, normalize_values, reg_params)
    
    localfit <- stan(file=stan_mode_file, data=stanDat, warmup = 1000, iter=2000, chains=4, init=0, control = list(adapt_delta = 0.80, max_treedepth = 12))
    
    print("Sampling done. Saving..")
    
    modelcache <- paste0(local_model_cache, "/", target_name, "_blmm", ".rds")
    
    saveRDS(localfit, file=modelcache)
    print("Done.")
  }
  
  return(localfit)
}  

##################################################

mebn.multivariate_sampling <- function(inputdata, targetdata, predictor_columns, target_columns, group_column, local_model_cache = "models", stan_mode_file = "BLMM.stan", normalize_values = TRUE, reg_params = NULL)
{
  require(rstan)
  
  # Run Stan parallel on multiple cores
  rstan_options (auto_write=TRUE)
  options (mc.cores=parallel::detectCores ()) 
  
  target_name <- paste0(target_columns$Name, collapse = "_")
  
  localfit <- mebn.get_localfit(target_name, local_model_cache)
  
  # Use cached model if it exists
  if (is.null(localfit))
  {
    stanDat <- mebn.set_mv_model_parameters(predictor_columns, target_columns, group_column, inputdata, targetdata, normalize_values, reg_params)
    
    localfit <- stan(file=stan_mode_file, data=stanDat, warmup = 1000, iter=2000, chains=4, init=0, control = list(adapt_delta = 0.80, max_treedepth = 12))
    
    print("Sampling done. Saving..")
    
    modelcache <- paste0(local_model_cache, "/", target_name, "_blmm", ".rds")
    
    saveRDS(localfit, file=modelcache)
    print("Done.")
  }
  
  return(localfit)
}  

##################################################

mebn.two_level_multivariate_sampling <- function(inputdata, targetdata, predictor_columns, target_columns, group_column, subject_column, local_model_cache = "models", stan_mode_file = "BLMM.stan", normalize_values = TRUE, reg_params = NULL)
{
  require(rstan)
  
  # Run Stan parallel on multiple cores
  rstan_options (auto_write=TRUE)
  options (mc.cores=parallel::detectCores ()) 
  
  target_name <- paste0(target_columns$Name, collapse = "_")
  
  localfit <- mebn.get_localfit(target_name, local_model_cache)
  
  # Use cached model if it exists
  if (is.null(localfit))
  {
    stanDat <- mebn.set_mv_two_level_model_parameters(predictor_columns, target_columns, group_column, subject_column, inputdata, targetdata, normalize_values, reg_params)
    
    #localfit <- stan(file=stan_mode_file, data=stanDat, warmup = 100, iter=200, chains=4, init=0)
    localfit <- stan(file=stan_mode_file, data=stanDat, warmup = 1000, iter=2000, chains=4, init=0, control = list(adapt_delta = 0.80, max_treedepth = 12))
    
    print("Sampling done. Saving..")
    
    modelcache <- paste0(local_model_cache, "/", target_name, "_blmm", ".rds")
    
    saveRDS(localfit, file=modelcache)
    print("Done.")
  }
  
  return(localfit)
}  

##################################################
  
mebn.multivariate_csr_sampling <- function(inputdata, targetdata, predictor_columns, target_columns, group_column, local_model_cache = "models", stan_mode_file = "BLMM.stan", normalize_values = TRUE, reg_params = NULL)
{
  require(rstan)
  
  # Run Stan parallel on multiple cores
  rstan_options (auto_write=TRUE)
  options (mc.cores=parallel::detectCores ()) 
  
  target_name <- paste0(target_columns$Name, collapse = "_")
  
  localfit <- mebn.get_localfit(target_name, local_model_cache)
  
  # Use cached model if it exists
  if (is.null(localfit))
  {
    stanDat <- mebn.set_mv_csr_parameters(predictor_columns, target_columns, group_column, inputdata, targetdata, normalize_values, reg_params)
    
    localfit <- stan(file=stan_mode_file, data=stanDat, warmup = 1000, iter=2000, chains=1, init=0, control = list(adapt_delta = 0.80, max_treedepth = 12))
    
    print("Sampling done. Saving..")
    
    modelcache <- paste0(local_model_cache, "/", target_name, "_blmm", ".rds")
    
    saveRDS(localfit, file=modelcache)
    print("Done.")
  }
  
  return(localfit)
}  
  
##################################################

mebn.predict <- function(inputdata, predictor_columns, target_column, local_model_cache = "models", stan_model_file = "BLMM.stan", normalize_values = TRUE, model_params = NULL)
{
  require(rstan)
  
  # Run Stan parallel on multiple cores
  rstan_options (auto_write=TRUE)
  options (mc.cores=parallel::detectCores ()) 
  
  target_name <- as.vector(target_column$Name)

  stanDat <- mebn.set_prediction_parameters(predictor_columns, target_column, inputdata, normalize_values, model_params)
    
  localfit <- stan(file=stan_model_file, 
                   data=stanDat, 
                   iter=1000, 
                   chains=4)

  return(localfit)
}  

##################################################

mebn.variational_inference <- function(inputdata, targetdata, predictor_columns, target_column, group_column, local_model_cache = "models", stan_mode_file = "BLMM.stan", normalize_values = TRUE, reg_params = NULL)
{
  require(rstan)
  
  # Run Stan parallel on multiple cores
  rstan_options (auto_write=TRUE)
  options (mc.cores=parallel::detectCores ()) 
  
  target_name <- as.vector(target_column$Name)
  
  localfit <- mebn.get_localfit(target_name, local_model_cache)
  
  # Use cached model if it exists
  if (is.null(localfit))
  {
    stanDat <- mebn.set_model_parameters(predictor_columns, target_column, group_column, inputdata, targetdata, normalize_values, reg_params)
    localmodel <- stan_model(file = stan_mode_file)
    localfit <- vb(localmodel, data=stanDat, output_samples=2500, iter=10000)
    
    modelcache <- paste0(local_model_cache, "/", target_name, "_blmm", ".rds")
    saveRDS(localfit, file=modelcache)
  }
  
  return(localfit)
}  

##################################################

mebn.get_rootnodes <- function(g)
{
  which(sapply(sapply(V(g), function(x) neighbors(g,x, mode="in")), length) == 0)
}

##################################################

mebn.write_gexf <- function(reaction_graph, gexf_path)
{
  require(rgexf)

  MakeRGBA <- function(RGBstring, alpha)
  {
    strtodec <- function(rgb, b, e) { strtoi(paste0("0x", substr(rgb, b, e))) }
    
    RGBA <- data.frame(strtodec(RGBstring, 2, 3), strtodec(RGBstring, 4, 5), strtodec(RGBstring, 6, 7), alpha)  
    colnames(RGBA) <- c("r", "g", "b", "alpha")  
    
    return(RGBA)
  }
  
  graphdata <- get.data.frame(reaction_graph)
  
  nodeviz <- list(color = MakeRGBA(V(reaction_graph)$color, 1.0), size = V(reaction_graph)$size, shape = V(reaction_graph)$shape)
  edgeviz <- list(shape = E(reaction_graph)$shape)
  
  edgesatt <- data.frame(E(reaction_graph)$value, E(reaction_graph)$value_lCI, E(reaction_graph)$value_uCI, E(reaction_graph)$b_sigma, E(reaction_graph)$b_sigma_lCI, E(reaction_graph)$b_sigma_uCI)
  
  if (length(edgesatt) == 6)
  {
    colnames(edgesatt) <- c("value", "value_lCI", "value_uCI", "b_sigma", "b_sigma_lCI", "b_sigma_uCI")
  }
  
  nodesatt <- data.frame(V(reaction_graph)$value, V(reaction_graph)$value_lCI, V(reaction_graph)$value_uCI,V(reaction_graph)$type)
  if (length(nodesatt) == 4)  
  {
    colnames(nodesatt) <- c("value", "value_lCI", "value_uCI", "type")
  }
  
  #edgelabels <- data.frame(paste0("beta = ", round(E(reaction_graph)$mean, 3), ", b_sigma = ", round(2*(E(reaction_graph)$uCI - E(reaction_graph)$mean), 3)))
  
  write.gexf(
    defaultedgetype = "directed",
    nodes = data.frame(V(reaction_graph)$name, V(reaction_graph)$label),
    edges = get.edgelist(reaction_graph),
    edgesWeight = graphdata[,3],
    #edgesLabel = edgelabels,
    nodesVizAtt = nodeviz,
    edgesVizAtt = edgeviz,
    edgesAtt = edgesatt,
    nodesAtt = nodesatt,
    output = gexf_path
  )
}

###################################

mebn.read_gexf <- function(gefx_path)
{
  gexf_graph <- read.gexf(gefx_path)
  ig <- rgexf::gexf.to.igraph(gexf_graph)
  
  return(ig)
}

###################################

mebn.BetaLevelTest <- function(LocalModelSummary, PredictorId)
{
  abs(LocalModelSummary$fixef[PredictorId]) > 0.001
}

###################################

mebn.RanefTest <- function(localsummary, PredictorId)
{
  abs(localsummary$fixef[PredictorId]) > 0.001 ||
    abs(localsummary$ranef_sd[PredictorId]) > 0.001
}

###################################

mebn.PersonalSignificanceTest <- function(personal_coef)
{
  abs(personal_coef) > 0.001
}

##################################################

mebn.bipartite_model <- function(reaction_graph, inputdata, targetdata = NULL, predictor_columns, assumed_targets, group_column, local_model_cache, stan_model_file, local_estimation, normalize_values = TRUE, reg_params = NULL)
{
  for (c in 1:dim(assumed_targets)[1])
  {
    target_column <- assumed_targets[c,]
    target_name <- as.vector(target_column$Name)
    
    print(target_name)
    
    localfit <- local_estimation(inputdata, targetdata, predictor_columns, target_column, group_column, local_model_cache, stan_model_file, normalize_values, reg_params)
    
    # Extract model summary
    localsummary <- mebn.localsummary(localfit)

    # - Loop through betas for current target
    predictor_names <- as.vector(predictor_columns$Name)
    
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      
      # Attach the random variable
      reaction_graph <- reaction_graph + edge(c(predictor_name, target_name), 
                                              weight = localsummary$fixef[p], 
                                              value = localsummary$fixef[p], 
                                              value_lCI = localsummary$fixef_lCI[p],
                                              value_uCI = localsummary$fixef_uCI[p],
                                              b_sigma = localsummary$ranef_sd[p],
                                              b_sigma_lCI = localsummary$ranef_sd_lCI[p],
                                              b_sigma_uCI = localsummary$ranef_sd_uCI[p],
                                              shape   = "confband")

      # Fixed-effect
      reaction_graph <- reaction_graph + vertex(paste0("beta_", predictor_name, "_", target_name), 
                                                label=paste0("beta_", predictor_name), 
                                                type="beta", 
                                                value = localsummary$fixef[p], 
                                                value_lCI = localsummary$fixef_lCI[p],
                                                value_uCI = localsummary$fixef_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("beta_", predictor_name, "_", target_name), paste0("beta_", predictor_name, "_", target_name), shape = "arrow", weight = 1, type = "beta") 
      
      # Add random-effect for significant predictors
      reaction_graph <- reaction_graph + vertex(paste0("b_", predictor_name, "_", target_name), 
                                                label=paste0("b_", predictor_name), 
                                                type="b", 
                                                value = 0, 
                                                value_lCI = 0,
                                                value_uCI = 0,
                                                size = 0.5, 
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + vertex(paste0("b_sigma_", predictor_name, "_", target_name), 
                                                label="b_sigma", 
                                                type="b_sigma", 
                                                value = localsummary$ranef_sd[p],
                                                value_lCI = localsummary$ranef_sd_lCI[p],
                                                value_uCI = localsummary$ranef_sd_uCI[p],
                                                size = localsummary$ranef_sd[p], 
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("b_sigma_", predictor_name, "_", target_name), 
                                              paste0("b_", predictor_name, "_", target_name), 
                                              shape = "arrow", 
                                              weight = localsummary$ranef_sd[p], 
                                              type = "b_sigma")
      
      reaction_graph <- reaction_graph + edge(paste0("b_", predictor_name, "_", target_name), 
                                              target_name, shape = "arrow", 
                                              weight=1, 
                                              type = "b")
    }
    
  } # loop targets
  
  return(reaction_graph)
}

##################################################

mebn.bipartite_multivariate_pk_fppi <- function(reaction_graph, inputdata, targetdata = NULL, predictor_columns, assumed_targets, group_column, local_model_cache, stan_model_file, local_estimation, normalize_values = TRUE, reg_params = NULL)
{
  # In this version, targets are not looped, but estimated together

  if (assumed_targets[1,]$Name == "pk" && assumed_targets[2,]$Name == "fppi" || assumed_targets[2,]$Name == "pk" && assumed_targets[1,]$Name == "fppi") {
    print("Multivariate estimation of targets 'pk' and 'fppi'")
  } else {
    print("ERROR: Current version of mebn.bipartite_multivariate is fixed to pk and fppi target columns. Aborting.")
    return(NULL);
  }
    
  localfit <- local_estimation(inputdata, targetdata, predictor_columns, assumed_targets, group_column, local_model_cache, stan_model_file, normalize_values, reg_params)

  #return(localfit);

  for (c in 1:dim(assumed_targets)[1])
  {
    target_column <- assumed_targets[c,]
    target_name <- as.vector(target_column$Name)
    
    print(target_name)

    # Extract local model summary from a multivariate model
    localsummary <- mebn.localsummary_from_multivariate(localfit, target_name)
    
    # - Loop through betas for current target
    predictor_names <- as.vector(predictor_columns$Name)
    
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      
      # Attach the random variable
      reaction_graph <- reaction_graph + edge(c(predictor_name, target_name), 
                                              weight = localsummary$fixef[p], 
                                              value = localsummary$fixef[p], 
                                              value_lCI = localsummary$fixef_lCI[p],
                                              value_uCI = localsummary$fixef_uCI[p],
                                              b_sigma = localsummary$ranef_sd[p],
                                              b_sigma_lCI = localsummary$ranef_sd_lCI[p],
                                              b_sigma_uCI = localsummary$ranef_sd_uCI[p],
                                              shape   = "confband")
      
      # Fixed-effect
      reaction_graph <- reaction_graph + vertex(paste0("beta_", predictor_name, "_", target_name), 
                                                label=paste0("beta_", predictor_name), 
                                                type="beta", 
                                                value = localsummary$fixef[p], 
                                                value_lCI = localsummary$fixef_lCI[p],
                                                value_uCI = localsummary$fixef_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("beta_", predictor_name, "_", target_name), paste0("beta_", predictor_name, "_", target_name), shape = "arrow", weight = 1, type = "beta") 
      
      # Add random-effect for significant predictors
      reaction_graph <- reaction_graph + vertex(paste0("b_", predictor_name, "_", target_name), 
                                                label=paste0("b_", predictor_name), 
                                                type="b", 
                                                value = 0, 
                                                value_lCI = 0,
                                                value_uCI = 0,
                                                size = 0.5, 
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + vertex(paste0("b_sigma_", predictor_name, "_", target_name), 
                                                label="b_sigma", 
                                                type="b_sigma", 
                                                value = localsummary$ranef_sd[p],
                                                value_lCI = localsummary$ranef_sd_lCI[p],
                                                value_uCI = localsummary$ranef_sd_uCI[p],
                                                size = localsummary$ranef_sd[p], 
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("b_sigma_", predictor_name, "_", target_name), 
                                              paste0("b_", predictor_name, "_", target_name), 
                                              shape = "arrow", 
                                              weight = localsummary$ranef_sd[p], 
                                              type = "b_sigma")
      
      reaction_graph <- reaction_graph + edge(paste0("b_", predictor_name, "_", target_name), 
                                              target_name, shape = "arrow", 
                                              weight=1, 
                                              type = "b")
    }
    
  } # loop targets
  
  return(reaction_graph)
}

##################################################

mebn.bipartite_multivariate <- function(reaction_graph, inputdata, targetdata = NULL, predictor_columns, assumed_targets, group_column, local_model_cache, stan_model_file, local_estimation, normalize_values = TRUE, reg_params = NULL)
{
  # In this version, targets are not looped, but estimated together
  
  localfit <- local_estimation(inputdata, targetdata, predictor_columns, assumed_targets, group_column, local_model_cache, stan_model_file, normalize_values, reg_params)
  
  #return(localfit);
  
  for (c in 1:dim(assumed_targets)[1])
  {
    target_column <- assumed_targets[c,]
    target_name <- as.vector(target_column$Name)
    
    print(paste0(target_name, " - ", c))
    
    # Extract local model summary from a multivariate model
    localsummary <- mebn.localsummary_from_multivariate(localfit, c)
    
    # - Loop through betas for current target
    predictor_names <- as.vector(predictor_columns$Name)
    
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      
      # Attach the random variable
      reaction_graph <- reaction_graph + edge(c(predictor_name, target_name), 
                                              weight = localsummary$fixef[p], 
                                              value = localsummary$fixef[p], 
                                              value_lCI = localsummary$fixef_lCI[p],
                                              value_uCI = localsummary$fixef_uCI[p],
                                              b_sigma = localsummary$ranef_sd[p],
                                              b_sigma_lCI = localsummary$ranef_sd_lCI[p],
                                              b_sigma_uCI = localsummary$ranef_sd_uCI[p],
                                              shape   = "confband")
      
      # Fixed-effect
      reaction_graph <- reaction_graph + vertex(paste0("beta_", predictor_name, "_", target_name), 
                                                label=paste0("beta_", predictor_name), 
                                                type="beta", 
                                                value = localsummary$fixef[p], 
                                                value_lCI = localsummary$fixef_lCI[p],
                                                value_uCI = localsummary$fixef_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("beta_", predictor_name, "_", target_name), paste0("beta_", predictor_name, "_", target_name), shape = "arrow", weight = 1, type = "beta") 
      
      # Add random-effect for significant predictors
      reaction_graph <- reaction_graph + vertex(paste0("b_", predictor_name, "_", target_name), 
                                                label=paste0("b_", predictor_name), 
                                                type="b", 
                                                value = 0, 
                                                value_lCI = 0,
                                                value_uCI = 0,
                                                size = 0.5, 
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + vertex(paste0("b_sigma_", predictor_name, "_", target_name), 
                                                label="b_sigma", 
                                                type="b_sigma", 
                                                value = localsummary$ranef_sd[p],
                                                value_lCI = localsummary$ranef_sd_lCI[p],
                                                value_uCI = localsummary$ranef_sd_uCI[p],
                                                size = localsummary$ranef_sd[p], 
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("b_sigma_", predictor_name, "_", target_name), 
                                              paste0("b_", predictor_name, "_", target_name), 
                                              shape = "arrow", 
                                              weight = localsummary$ranef_sd[p], 
                                              type = "b_sigma")
      
      reaction_graph <- reaction_graph + edge(paste0("b_", predictor_name, "_", target_name), 
                                              target_name, shape = "arrow", 
                                              weight=1, 
                                              type = "b")
    }
    
  } # loop targets
  
  return(reaction_graph)
}

##################################################

mebn.bipartite_two_level_multivariate <- function(reaction_graph, inputdata, targetdata = NULL, predictor_columns, assumed_targets, subject_column, group_column, local_model_cache, stan_model_file, local_estimation, normalize_values = TRUE, reg_params = NULL)
{
  # In this version, targets are not looped, but estimated together
  
  localfit <- local_estimation(inputdata, targetdata, predictor_columns, assumed_targets, group_column, subject_column, local_model_cache, stan_model_file, normalize_values, reg_params)
  
  #return(localfit);
  
  for (c in 1:dim(assumed_targets)[1])
  {
    target_column <- assumed_targets[c,]
    target_name <- as.vector(target_column$Name)
    
    print(paste0(target_name, " - ", c))
    
    # Extract local model summary from a multivariate model
    localsummary <- mebn.localsummary_from_two_level_multivariate(localfit, c)
    
    # - Loop through betas for current target
    predictor_names <- as.vector(predictor_columns$Name)
    
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      
      # Attach the random variable
      reaction_graph <- reaction_graph + edge(c(predictor_name, target_name), 
                                              weight = localsummary$fixef[p], 
                                              value = localsummary$fixef[p], 
                                              value_lCI = localsummary$fixef_lCI[p],
                                              value_uCI = localsummary$fixef_uCI[p],
                                              b1_sigma = localsummary$ranef1_sd[p],
                                              b1_sigma_lCI = localsummary$ranef1_sd_lCI[p],
                                              b1_sigma_uCI = localsummary$ranef1_sd_uCI[p],
                                              b2_sigma = localsummary$ranef2_sd[p],
                                              b2_sigma_lCI = localsummary$ranef2_sd_lCI[p],
                                              b2_sigma_uCI = localsummary$ranef2_sd_uCI[p],
                                              shape   = "confband")
      
      # Fixed-effect
      reaction_graph <- reaction_graph + vertex(paste0("beta_", predictor_name, "_", target_name), 
                                                label=paste0("beta_", predictor_name), 
                                                type="beta", 
                                                value = localsummary$fixef[p], 
                                                value_lCI = localsummary$fixef_lCI[p],
                                                value_uCI = localsummary$fixef_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("beta_", predictor_name, "_", target_name), paste0("beta_", predictor_name, "_", target_name), shape = "arrow", weight = 1, type = "beta") 
      
      # Add random-effect for significant predictors
      reaction_graph <- reaction_graph + vertex(paste0("b1_", predictor_name, "_", target_name), 
                                                label=paste0("b1_", predictor_name), 
                                                type="b1", 
                                                value = 0, 
                                                value_lCI = 0,
                                                value_uCI = 0,
                                                size = 0.5, 
                                                shape = "circle")

      reaction_graph <- reaction_graph + vertex(paste0("b2_", predictor_name, "_", target_name), 
                                                label=paste0("b2_", predictor_name), 
                                                type="b2", 
                                                value = 0, 
                                                value_lCI = 0,
                                                value_uCI = 0,
                                                size = 0.5, 
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + vertex(paste0("b1_sigma_", predictor_name, "_", target_name), 
                                                label="b1_sigma", 
                                                type="b1_sigma", 
                                                value = localsummary$ranef1_sd[p],
                                                value_lCI = localsummary$ranef1_sd_lCI[p],
                                                value_uCI = localsummary$ranef1_sd_uCI[p],
                                                size = localsummary$ranef1_sd[p], 
                                                shape = "circle")

      reaction_graph <- reaction_graph + vertex(paste0("b2_sigma_", predictor_name, "_", target_name), 
                                                label="b2_sigma", 
                                                type="b2_sigma", 
                                                value = localsummary$ranef2_sd[p],
                                                value_lCI = localsummary$ranef2_sd_lCI[p],
                                                value_uCI = localsummary$ranef2_sd_uCI[p],
                                                size = localsummary$ranef2_sd[p], 
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("b1_sigma_", predictor_name, "_", target_name), 
                                              paste0("b1_", predictor_name, "_", target_name), 
                                              shape = "arrow", 
                                              weight = localsummary$ranef1_sd[p], 
                                              type = "b1_sigma")

      reaction_graph <- reaction_graph + edge(paste0("b2_sigma_", predictor_name, "_", target_name), 
                                              paste0("b2_", predictor_name, "_", target_name), 
                                              shape = "arrow", 
                                              weight = localsummary$ranef2_sd[p], 
                                              type = "b2_sigma")
      
      reaction_graph <- reaction_graph + edge(paste0("b1_", predictor_name, "_", target_name), 
                                              target_name, shape = "arrow", 
                                              weight=1, 
                                              type = "b1")

      reaction_graph <- reaction_graph + edge(paste0("b2_", predictor_name, "_", target_name), 
                                              target_name, shape = "arrow", 
                                              weight=1, 
                                              type = "b2")
    }
    
  } # loop targets
  
  return(reaction_graph)
}


##################################################

mebn.bipartite_expfam_bn <- function(reaction_graph, predictor_columns, assumed_targets, target_models_dirs)
{
  for (c in 1:dim(assumed_targets)[1])
  {
    target_column <- assumed_targets[c,]
    target_name <- as.vector(target_column$Name)
    
    # Load previously sampled model
    local_model_cache <- target_models_dirs[target_models_dirs$Name==target_name,]$modelcache 
    
    print(target_name)
    print(local_model_cache)
    
    localfit <- mebn.get_localfit(target_name, local_model_cache)
    
    # Extract model summary
    localsummary <- mebn.localsummary(localfit)
    
    # - Loop through betas for current target
    predictor_names <- as.vector(predictor_columns$Name)
    
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      
      # Attach the random variable
      reaction_graph <- reaction_graph + edge(c(predictor_name, target_name), 
                                              weight = localsummary$fixef[p], 
                                              value = localsummary$fixef[p], 
                                              value_lCI = localsummary$fixef_lCI[p],
                                              value_uCI = localsummary$fixef_uCI[p],
                                              b_sigma = localsummary$ranef_sd[p],
                                              b_sigma_lCI = localsummary$ranef_sd_lCI[p],
                                              b_sigma_uCI = localsummary$ranef_sd_uCI[p],
                                              shape   = "confband")
      
      # Fixed-effect
      reaction_graph <- reaction_graph + vertex(paste0("beta_", predictor_name, "_", target_name), 
                                                label=paste0("beta_", predictor_name), 
                                                type="beta", 
                                                value = localsummary$fixef[p], 
                                                value_lCI = localsummary$fixef_lCI[p],
                                                value_uCI = localsummary$fixef_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("beta_", predictor_name, "_", target_name), paste0("beta_", predictor_name, "_", target_name), shape = "arrow", weight = 1, type = "beta") 
      
      # Add random-effect for significant predictors
      reaction_graph <- reaction_graph + vertex(paste0("b_", predictor_name, "_", target_name), 
                                                label=paste0("b_", predictor_name), 
                                                type="b", 
                                                value = 0, 
                                                value_lCI = 0,
                                                value_uCI = 0,
                                                size = 0.5, 
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + vertex(paste0("b_sigma_", predictor_name, "_", target_name), 
                                                label="b_sigma", 
                                                type="b_sigma", 
                                                value = localsummary$ranef_sd[p],
                                                value_lCI = localsummary$ranef_sd_lCI[p],
                                                value_uCI = localsummary$ranef_sd_uCI[p],
                                                size = localsummary$ranef_sd[p], 
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("b_sigma_", predictor_name, "_", target_name), 
                                              paste0("b_", predictor_name, "_", target_name), 
                                              shape = "arrow", 
                                              weight = localsummary$ranef_sd[p], 
                                              type = "b_sigma")
      
      reaction_graph <- reaction_graph + edge(paste0("b_", predictor_name, "_", target_name), 
                                              target_name, shape = "arrow", 
                                              weight=1, 
                                              type = "b")
    }
    
  } # loop targets
  
  return(reaction_graph)
}
##################################################

mebn.extract_personal_graph <- function(person_id, reaction_graph, personal_model_dir, predictor_columns, assumed_targets, local_distributions, normalized_personal_data, original_personal_data)
{
  library(igraph)
  library(rstan)
  
  # Personal graph and sampled distributions are stored in this directory
  dir.create(personal_model_dir, showWarnings = FALSE)
  
  for (c in 1:dim(assumed_targets)[1])
  {
    target_column <- assumed_targets[c,]
    target_name <- as.vector(target_column$Name)
    #print(target_name)
    
    localfit_directory <- local_distributions
    localfit <- mebn.get_localfit(target_name, localfit_directory)

    # Add parameters to target distribution
    target_vertex <- V(reaction_graph)[target_name]
    
    # - shape of gamma distribution, common for all the persons
    localfit_gamma_shape <- extract(localfit, pars = c("g_alpha"))
    g_alpha_file <- paste0(target_name, "_g_alpha.rds")
    saveRDS(localfit_gamma_shape, paste0(personal_model_dir, "/", g_alpha_file))
    target_vertex$g_alpha <- g_alpha_file
    
    reaction_graph <- set_vertex_attr(reaction_graph, "g_alpha", target_vertex, g_alpha_file)
    #V(g)["pk"]$g_alpha
    
    # - personal intercept
    localfit_p_intercept <- extract(localfit, pars = c(paste0("personal_intercept[",person_id,"]")))
    int_samples_file <- paste0(target_name, "_personal_intercept.rds")
    saveRDS(localfit_p_intercept, paste0(personal_model_dir, "/", int_samples_file))
    
    reaction_graph <- set_vertex_attr(reaction_graph, "personal_intercept", target_vertex, int_samples_file)
    
    # extract personal effects from the local distribution
    pe <- mebn.personal_effects(localfit, person_id)
    
    # - Loop through betas for current target
    predictor_names <- as.vector(predictor_columns$Name)
    
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      pred_column <- predictor_columns[p,]
      
      # Set personally estimated predictor distribution
      pred_vertex <- V(reaction_graph)[predictor_name]
      
      dist_string <- "error: distribution not implemented"
      norm_dist_string <- "error: distribution not implemented"
      
      if (pred_column$Distribution == "Bernoulli")
      {
        dist_string <- paste0("B(", toString(original_personal_data[1,p]), ")")
        norm_dist_string <- dist_string
      }
      else if (pred_column$Distribution == "Gaussian")
      {
        dist_string <- paste0("N(", toString(mean(original_personal_data[,p])),",",toString(sd(original_personal_data[,p])),")")
        norm_dist_string <- paste0("N(", toString(mean(normalized_personal_data[,p])),",",toString(sd(normalized_personal_data[,p])),")")
      }
      
      reaction_graph <- set_vertex_attr(reaction_graph, "distribution", pred_vertex, dist_string)
      reaction_graph <- set_vertex_attr(reaction_graph, "normdist", pred_vertex, norm_dist_string)
      
      # Attach the random variable
      # Data values are also stored in edge attributes for easier visualization
      # and where to find the full sampled distribution
      effect_property <- paste0("personal_effect[", person_id, ",", p ,"]")
      personal_beta <- extract(localfit, pars = effect_property)
      personal_beta <- unlist(personal_beta, use.names=FALSE) # convert from one item list to vector
      effect_name <- paste0(predictor_name,"_",target_name,"_beta")
      saveRDS(personal_beta, paste0(personal_model_dir, "/", effect_name, ".rds"))
      
      reaction_graph <- reaction_graph + edge(c(predictor_name, target_name),
                                              weight = pe$personal_effect[p], 
                                              value = pe$personal_effect[p], 
                                              value_lCI = pe$personal_effect_lCI[p],
                                              value_uCI = pe$personal_effect_uCI[p],
                                              distribution = paste0(effect_name,".rds"),
                                              b = pe$b[p], 
                                              b_lCI = pe$b_lCI[p],
                                              b_uCI = pe$b_uCI[p])
      
      # Personal effect (beta + b) vertex
      
      reaction_graph <- reaction_graph + vertex(paste0("personal_", predictor_name, "_", target_name), 
                                                label=paste0("personal_", predictor_name), 
                                                type="personal", color="#AAAAAA", 
                                                value = pe$personal_effect[p], 
                                                value_lCI = pe$personal_effect_lCI[p],
                                                value_uCI = pe$personal_effect_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("personal_", predictor_name, "_", target_name), target_name, shape = "arrow", weight = 1, type = "personal") 
      
      # Personal variation (b)
      reaction_graph <- reaction_graph + vertex(paste0("b_", predictor_name, "_", target_name), 
                                                label=paste0("b_", predictor_name), 
                                                type="b", color="#AAAAAA", 
                                                value = pe$b[p], 
                                                value_lCI = pe$b_lCI[p],
                                                value_uCI = pe$b_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("b_", predictor_name, "_", target_name), target_name, shape = "arrow", weight = 1, type = "b") 
    }
    
  } # loop targets
  
  write.graph(reaction_graph, paste0(personal_model_dir,"/personal_graph.graphml"), "graphml")
  
  return(reaction_graph)
}

##################################################

mebn.extract_personal_graph_from_mv <- function(person_id, reaction_graph, personal_model_dir, predictor_columns, target_variables, multivariate_modeldir, normalized_personal_data, original_personal_data, original_personal_concentrations, datadesc)
{
  library(igraph)
  library(rstan)
  
  # Personal graph and sampled distributions are stored in this directory
  dir.create(personal_model_dir, showWarnings = FALSE)

  multivariate_model_name <- paste0(target_variables$Name, collapse = "_")
  localfit <- mebn.get_localfit(multivariate_model_name, multivariate_modeldir)
  
  for (c in 1:dim(target_variables)[1])
  {
    target_column <- target_variables[c,]
    target_name <- as.vector(target_column$Name)

    # Add parameters to target distribution
    target_vertex <- V(reaction_graph)[target_name]
    
    # Add observed concentration distribution parameters
    reaction_graph <- set_vertex_attr(reaction_graph, "value_mean", target_vertex, mean(original_personal_concentrations[,c]))
    reaction_graph <- set_vertex_attr(reaction_graph, "value_sd", target_vertex, sd(original_personal_concentrations[,c]))
    
    # - shape of gamma distribution, common for all the persons
    localfit_gamma_shape <- extract(localfit, pars = c(paste0("g_alpha[",c,"]")))
    g_alpha_file <- paste0(target_name, "_g_alpha.rds")
    saveRDS(localfit_gamma_shape, paste0(personal_model_dir, "/", g_alpha_file))
    target_vertex$g_alpha <- g_alpha_file
    
    reaction_graph <- set_vertex_attr(reaction_graph, "g_alpha", target_vertex, g_alpha_file)

    # - personal intercept

    localfit_p_intercept <- extract(localfit, pars = c(paste0("personal_intercept[",person_id,",",c,"]")))
    int_samples_file <- paste0(target_name, "_personal_intercept.rds")
    saveRDS(localfit_p_intercept, paste0(personal_model_dir, "/", int_samples_file))
    
    reaction_graph <- set_vertex_attr(reaction_graph, "personal_intercept", target_vertex, int_samples_file)
    
    # extract personal effects from the local distribution
    pe <- mebn.personal_effects_from_multivariate(localfit, c, person_id)
    
    # - Loop through betas for current target
    predictor_names <- as.vector(predictor_columns$Name)
    
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      pred_column <- predictor_columns[p,]
      
      # Set personally estimated predictor distribution
      pred_vertex <- V(reaction_graph)[predictor_name]
      
      dist_string <- "error: distribution not implemented"
      norm_dist_string <- "error: distribution not implemented"
      
      if (pred_column$Distribution == "Bernoulli")
      {
        dist_string <- paste0("B(", toString(original_personal_data[1,p]), ")")
        norm_dist_string <- dist_string
      }
      else if (pred_column$Distribution == "Gaussian")
      {
        # TODO: Bayesian modeling
        
        org_mean <- mean(original_personal_data[,p])
        org_sd <- sd(original_personal_data[,p])
        norm_mean <- mean(normalized_personal_data[,p])
        norm_sd <- sd(normalized_personal_data[,p])
        
        # Standard deviation can be zero is all observations are same. This will break rest of the analysis, so we add very small sd that is there anyway.
        if (org_sd == 0)
        {
          org_sd <- 0.0001
        }
        if (norm_sd == 0)
        {
          norm_sd <- 0.0001
        }
        
        dist_string <- paste0("N(", toString(org_mean),",",toString(org_sd),")")
        norm_dist_string <- paste0("N(", toString(norm_mean),",",toString(norm_sd),")")

        reaction_graph <- set_vertex_attr(reaction_graph, "value_mean", pred_vertex, org_mean)
        
        # We cannot generate values from this distribution
        #if (is.na(org_sd))
        #{  
        #  reaction_graph$generable <- TRUE
        #}
      }
      
      reaction_graph <- set_vertex_attr(reaction_graph, "distribution", pred_vertex, dist_string)
      reaction_graph <- set_vertex_attr(reaction_graph, "normdist", pred_vertex, norm_dist_string)
      
      # Attach the random variable
      # Data values are also stored in edge attributes for easier visualization
      # and where to find the full sampled distribution
      effect_property <- paste0("personal_effect[",person_id,",",c,",",p,"]")
      personal_beta <- extract(localfit, pars = effect_property)
      personal_beta <- unlist(personal_beta, use.names=FALSE) # convert from one item list to vector
      effect_name <- paste0(predictor_name,"_",target_name,"_beta")
      saveRDS(personal_beta, paste0(personal_model_dir, "/", effect_name, ".rds"))
      
      reaction_graph <- reaction_graph + edge(c(predictor_name, target_name),
                                              weight = pe$personal_effect[p], 
                                              value = pe$personal_effect[p], 
                                              value_lCI = pe$personal_effect_lCI[p],
                                              value_uCI = pe$personal_effect_uCI[p],
                                              distribution = paste0(effect_name,".rds"),
                                              b = pe$b[p], 
                                              b_lCI = pe$b_lCI[p],
                                              b_uCI = pe$b_uCI[p])
      
      # Personal effect (beta + b) vertex
      
      reaction_graph <- reaction_graph + vertex(paste0("personal_", predictor_name, "_", target_name), 
                                                label=paste0("personal_", predictor_name), 
                                                type="personal", color="#AAAAAA", 
                                                distribution = paste0(effect_name,".rds"),
                                                value = pe$personal_effect[p], 
                                                value_lCI = pe$personal_effect_lCI[p],
                                                value_uCI = pe$personal_effect_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("personal_", predictor_name, "_", target_name), target_name, shape = "arrow", weight = 1, type = "personal") 
      
      # Personal variation (b)
      reaction_graph <- reaction_graph + vertex(paste0("b_", predictor_name, "_", target_name), 
                                                label=paste0("b_", predictor_name), 
                                                type="b", color="#AAAAAA", 
                                                value = pe$b[p], 
                                                value_lCI = pe$b_lCI[p],
                                                value_uCI = pe$b_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("b_", predictor_name, "_", target_name), target_name, shape = "arrow", weight = 1, type = "b") 
    }
    
  } # loop targets
  
  print(paste0("Writing personal graph '", personal_model_dir,"/personal_graph.graphml'"))
  
  write.graph(reaction_graph, paste0(personal_model_dir,"/personal_graph.graphml"), "graphml")
  
  return(reaction_graph)
}

##################################################

mebn.extract_multilevel_graph <- function(person_id, group_id, reaction_graph, adjusted_model_dir, predictor_columns, target_variables, multivariate_modeldir, normalized_personal_data, original_personal_data, original_personal_concentrations, datadesc)
{
  library(igraph)
  library(rstan)
  
  # Personal graph and sampled distributions are stored in this directory
  dir.create(adjusted_model_dir, showWarnings = FALSE)
  
  multivariate_model_name <- paste0(target_variables$Name, collapse = "_")
  localfit <- mebn.get_localfit(multivariate_model_name, multivariate_modeldir)
  
  for (c in 1:dim(target_variables)[1])
  {
    target_column <- target_variables[c,]
    target_name <- as.vector(target_column$Name)
    
    # Add parameters to target distribution
    target_vertex <- V(reaction_graph)[target_name]
    
    # Add observed concentration distribution parameters
    reaction_graph <- set_vertex_attr(reaction_graph, "value_mean", target_vertex, mean(original_personal_concentrations[,c]))
    reaction_graph <- set_vertex_attr(reaction_graph, "value_sd", target_vertex, sd(original_personal_concentrations[,c]))
    
    # - shape of gamma distribution, common for all the persons
    localfit_gamma_shape <- extract(localfit, pars = c(paste0("g_alpha[",c,"]")))
    g_alpha_file <- paste0(target_name, "_g_alpha.rds")
    saveRDS(localfit_gamma_shape, paste0(personal_model_dir, "/", g_alpha_file))
    target_vertex$g_alpha <- g_alpha_file
    
    reaction_graph <- set_vertex_attr(reaction_graph, "g_alpha", target_vertex, g_alpha_file)
    
    # - adjusted intercept
    
    localfit_p_intercept <- extract(localfit, pars = c(paste0("personal_intercept[",person_id,",",c,"]")))
    int_samples_file <- paste0(target_name, "_personal_intercept.rds")
    saveRDS(localfit_p_intercept, paste0(personal_model_dir, "/", int_samples_file))
    
    reaction_graph <- set_vertex_attr(reaction_graph, "personal_intercept", target_vertex, int_samples_file)

    # extract adjusted effects from the local distribution
    pe <- mebn.adjusted_effects_from_multivariate(localfit, c, person_id, group_id)
    
    # - Loop through betas for current target
    predictor_names <- as.vector(predictor_columns$Name)
    
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      pred_column <- predictor_columns[p,]

      # - amounts of variations at different levels. these are general effects.
      localsummary <- mebn.localsummary_from_two_level_multivariate(localfit, c)
      
      # add general effect
      reaction_graph <- reaction_graph + vertex(paste0(predictor_name, "_", target_name), 
                                                label=paste0(predictor_name, "_", target_name), 
                                                type="general_effect", color="#AAAAAA", 
                                                value = localsummary$fixef[p],
                                                mode = localsummary$fixef_mode[p],
                                                value_lCI = localsummary$fixef_lCI[p],
                                                value_uCI = localsummary$fixef_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0(predictor_name, "_", target_name), target_name, shape = "arrow", weight = 1, type = "general_effect") 
      #

      reaction_graph <- reaction_graph + vertex(paste0("sigma_g_", predictor_name, "_", target_name), 
                                                label=paste0("sigma_g_", predictor_name, "_", target_name), 
                                                type="sigma_g_sd", color="#AAAAAA", 
                                                value = localsummary$ranef1_sd[p], 
                                                value_lCI = localsummary$ranef1_sd_lCI[p],
                                                value_uCI = localsummary$ranef1_sd_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("sigma_g_", predictor_name, "_", target_name), target_name, shape = "arrow", weight = 1, type = "sigma_g") 
      
      reaction_graph <- reaction_graph + vertex(paste0("sigma_b_", predictor_name, "_", target_name),
                                                label=paste0("sigma_b_", predictor_name, "_", target_name),
                                                type="sigma_b", color="#AAAAAA",
                                                value = localsummary$ranef2_sd[p],
                                                value_lCI = localsummary$ranef2_sd_lCI[p],
                                                value_uCI = localsummary$ranef2_sd_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("sigma_b_", predictor_name, "_", target_name), target_name, shape = "arrow", weight = 1, type = "sigma_b")
      
      # Set personally estimated predictor distribution
      pred_vertex <- V(reaction_graph)[predictor_name]
      
      dist_string <- "error: distribution not implemented"
      norm_dist_string <- "error: distribution not implemented"
      
      if (pred_column$Distribution == "Bernoulli")
      {
        dist_string <- paste0("B(", toString(original_personal_data[1,p]), ")")
        norm_dist_string <- dist_string
      }
      else if (pred_column$Distribution == "Gaussian")
      {
        org_mean <- mean(original_personal_data[,p])
        org_sd <- sd(original_personal_data[,p])
        norm_mean <- mean(normalized_personal_data[,p])
        norm_sd <- sd(normalized_personal_data[,p])
        
        # Standard deviation can be zero is all observations are same. This will break rest of the analysis, so we add very small sd that is there anyway.
        if (org_sd == 0)
        {
          org_sd <- 0.0001
        }
        if (norm_sd == 0)
        {
          norm_sd <- 0.0001
        }
        
        dist_string <- paste0("N(", toString(org_mean),",",toString(org_sd),")")
        norm_dist_string <- paste0("N(", toString(norm_mean),",",toString(norm_sd),")")
        
        reaction_graph <- set_vertex_attr(reaction_graph, "value_mean", pred_vertex, org_mean)
      }
      else if (pred_column$Distribution == "LogNormal")
      {
        org_mean <- mean(original_personal_data[,p])
        org_sd <- sd(original_personal_data[,p])
        norm_mean <- mean(normalized_personal_data[,p])
        norm_sd <- sd(normalized_personal_data[,p])
        
        # Standard deviation can be zero is all observations are same. This will break rest of the analysis, so we add very small sd that is there anyway.
        if (org_sd == 0)
        {
          org_sd <- 0.0001
        }
        if (norm_sd == 0)
        {
          norm_sd <- 0.0001
        }
        
        dist_string <- paste0("LN(", toString(org_mean),",",toString(org_sd),")")
        norm_dist_string <- paste0("LN(", toString(norm_mean),",",toString(norm_sd),")")
        
        reaction_graph <- set_vertex_attr(reaction_graph, "value_mean", pred_vertex, mean(original_personal_data[,p]))
      }
      
      reaction_graph <- set_vertex_attr(reaction_graph, "distribution", pred_vertex, dist_string)
      reaction_graph <- set_vertex_attr(reaction_graph, "normdist", pred_vertex, norm_dist_string)
      
      # Attach the random variable
      # Data values are also stored in edge attributes for easier visualization
      # and where to find the full sampled distribution
      personal_effect_property <- paste0("personal_effect[",person_id,",",c,",",p,"]")
      personal_beta <- extract(localfit, pars = personal_effect_property)
      personal_beta <- unlist(personal_beta, use.names=FALSE) # convert from one item list to vector
      effect_name <- paste0(predictor_name,"_",target_name,"_beta")
      saveRDS(personal_beta, paste0(personal_model_dir, "/", effect_name, ".rds"))
      
      reaction_graph <- reaction_graph + edge(c(predictor_name, target_name),
                                              weight = pe$personal_effect[p], 
                                              value = pe$personal_effect[p],
                                              value_lCI = pe$personal_effect_lCI[p],
                                              value_uCI = pe$personal_effect_uCI[p],
                                              group = pe$group_effect[p],
                                              group_lCI = pe$group_effect_lCI[p],
                                              group_uCI = pe$group_effect_uCI[p],
                                              distribution = paste0(effect_name,".rds"),
                                              b = pe$b[p], 
                                              b_lCI = pe$b_lCI[p],
                                              b_uCI = pe$b_uCI[p])
      
      # Personal effect (beta + g + b) vertex
      
      reaction_graph <- reaction_graph + vertex(paste0("personal_", predictor_name, "_", target_name), 
                                                label=paste0("personal_", predictor_name), 
                                                type="personal", color="#AAAAAA", 
                                                distribution = paste0(effect_name,".rds"),
                                                value = pe$personal_effect[p], 
                                                mode = pe$personal_effect_mode[p],
                                                value_lCI = pe$personal_effect_lCI[p],
                                                value_uCI = pe$personal_effect_uCI[p],
                                                group = pe$group_effect[p],
                                                group_mode = pe$group_effect_mode[p],
                                                group_lCI = pe$group_effect_lCI[p],
                                                group_uCI = pe$group_effect_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("personal_", predictor_name, "_", target_name), target_name, shape = "arrow", weight = 1, type = "personal") 
      
      # Personal variation (b)
      reaction_graph <- reaction_graph + vertex(paste0("b_", predictor_name, "_", target_name), 
                                                label=paste0("b_", predictor_name), 
                                                type="b", color="#AAAAAA", 
                                                value = pe$b[p],
                                                mode = pe$b_mode[p],
                                                value_lCI = pe$b_lCI[p],
                                                value_uCI = pe$b_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("b_", predictor_name, "_", target_name), target_name, shape = "arrow", weight = 1, type = "b") 

    }
    
  } # loop targets
  
  print(paste0("Writing personal graph '", personal_model_dir,"/personal_graph.graphml'"))
  
  write.graph(reaction_graph, paste0(personal_model_dir,"/personal_graph.graphml"), "graphml")
  
  return(reaction_graph)
}

##################################################

mebn.extract_personal_graph_from_mv_fixedtargets <- function(person_id, reaction_graph, personal_model_dir, predictor_columns, target_variables, multivariate_modeldir, normalized_personal_data, original_personal_data, datadesc)
{
  library(igraph)
  library(rstan)
  
  # Personal graph and sampled distributions are stored in this directory
  dir.create(personal_model_dir, showWarnings = FALSE)
  
  multivariate_model_name <- paste0(target_variables$Name, collapse = "_")
  localfit <- mebn.get_localfit(multivariate_model_name, multivariate_modeldir)
  
  for (c in 1:dim(target_variables)[1])
  {
    target_column <- target_variables[c,]
    target_name <- as.vector(target_column$Name)
    
    # Add parameters to target distribution
    target_vertex <- V(reaction_graph)[target_name]
    
    # - shape of gamma distribution, common for all the persons
    localfit_gamma_shape <- extract(localfit, pars = c(paste0("g_alpha_",target_name)))
    g_alpha_file <- paste0(target_name, "_g_alpha.rds")
    saveRDS(localfit_gamma_shape, paste0(personal_model_dir, "/", g_alpha_file))
    target_vertex$g_alpha <- g_alpha_file
    
    reaction_graph <- set_vertex_attr(reaction_graph, "g_alpha", target_vertex, g_alpha_file)
    #V(g)["pk"]$g_alpha
    
    # - personal intercept
    localfit_p_intercept <- extract(localfit, pars = c(paste0("personal_intercept_",target_name,"[",person_id,"]")))
    int_samples_file <- paste0(target_name, "_personal_intercept.rds")
    saveRDS(localfit_p_intercept, paste0(personal_model_dir, "/", int_samples_file))
    
    reaction_graph <- set_vertex_attr(reaction_graph, "personal_intercept", target_vertex, int_samples_file)
    
    # extract personal effects from the local distribution
    pe <- mebn.personal_effects_from_multivariate(localfit, target_name, person_id)
    
    # - Loop through betas for current target
    predictor_names <- as.vector(predictor_columns$Name)
    
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      pred_column <- predictor_columns[p,]
      
      # Set personally estimated predictor distribution
      pred_vertex <- V(reaction_graph)[predictor_name]
      
      dist_string <- "error: distribution not implemented"
      norm_dist_string <- "error: distribution not implemented"
      
      if (pred_column$Distribution == "Bernoulli")
      {
        dist_string <- paste0("B(", toString(original_personal_data[1,p]), ")")
        norm_dist_string <- dist_string
      }
      else if (pred_column$Distribution == "Gaussian")
      {
        # TODO: Bayesian modeling
        
        org_mean <- mean(original_personal_data[,p])
        org_sd <- sd(original_personal_data[,p])
        norm_mean <- mean(normalized_personal_data[,p])
        norm_sd <- sd(normalized_personal_data[,p])
        
        dist_string <- paste0("N(", toString(org_mean),",",toString(org_sd),")")
        norm_dist_string <- paste0("N(", toString(norm_mean),",",toString(norm_sd),")")
        
        # We cannot generate values from this distribution
        #if (is.na(org_sd))
        #{  
        #  reaction_graph$generable <- TRUE
        #}
      }
      
      reaction_graph <- set_vertex_attr(reaction_graph, "distribution", pred_vertex, dist_string)
      reaction_graph <- set_vertex_attr(reaction_graph, "normdist", pred_vertex, norm_dist_string)
      
      # Attach the random variable
      # Data values are also stored in edge attributes for easier visualization
      # and where to find the full sampled distribution
      effect_property <- paste0("personal_effect_", target_name, "[", person_id, ",", p ,"]")
      personal_beta <- extract(localfit, pars = effect_property)
      personal_beta <- unlist(personal_beta, use.names=FALSE) # convert from one item list to vector
      effect_name <- paste0(predictor_name,"_",target_name,"_beta")
      saveRDS(personal_beta, paste0(personal_model_dir, "/", effect_name, ".rds"))
      
      reaction_graph <- reaction_graph + edge(c(predictor_name, target_name),
                                              weight = pe$personal_effect[p], 
                                              value = pe$personal_effect[p], 
                                              value_lCI = pe$personal_effect_lCI[p],
                                              value_uCI = pe$personal_effect_uCI[p],
                                              distribution = paste0(effect_name,".rds"),
                                              b = pe$b[p], 
                                              b_lCI = pe$b_lCI[p],
                                              b_uCI = pe$b_uCI[p])
      
      # Personal effect (beta + b) vertex
      
      reaction_graph <- reaction_graph + vertex(paste0("personal_", predictor_name, "_", target_name), 
                                                label=paste0("personal_", predictor_name), 
                                                type="personal", color="#AAAAAA", 
                                                value = pe$personal_effect[p], 
                                                value_lCI = pe$personal_effect_lCI[p],
                                                value_uCI = pe$personal_effect_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("personal_", predictor_name, "_", target_name), target_name, shape = "arrow", weight = 1, type = "personal") 
      
      # Personal variation (b)
      reaction_graph <- reaction_graph + vertex(paste0("b_", predictor_name, "_", target_name), 
                                                label=paste0("b_", predictor_name), 
                                                type="b", color="#AAAAAA", 
                                                value = pe$b[p], 
                                                value_lCI = pe$b_lCI[p],
                                                value_uCI = pe$b_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("b_", predictor_name, "_", target_name), target_name, shape = "arrow", weight = 1, type = "b") 
    }
    
  } # loop targets
  
  print(paste0("Writing personal graph '", personal_model_dir,"/personal_graph.graphml'"))
  
  write.graph(reaction_graph, paste0(personal_model_dir,"/personal_graph.graphml"), "graphml")
  
  return(reaction_graph)
}

##################################################

mebn.Generate <- function(reaction_graph, graph_dir, query, queried_nodes, proposal_distribution_params, stan_model_file, iterations = 10000, point_est = FALSE, posterior_samples = 20)
{
  library(rstan)
  library(igraph)
  library(stringr)
  
  predictor_nodes <- V(reaction_graph)[type == 100]
  target_nodes <- V(reaction_graph)[type == 200]

  predictors <- length(predictor_nodes)
  targets <- length(target_nodes)
  
  # Collect model parameters from the graph
  full_sample_size <- 2000 # warmup and sampling stored in RDS
  density_samples <- 1001:2000 # portion of samples to use 
  
  beta_samples <- array(-1, dim = c(targets, predictors, length(density_samples)))
  intercept_samples <- array(-1, dim = c(targets, length(density_samples)))
  alpha_samples <- array(-1, dim = c(targets, length(density_samples)))
  
  prior_stats <- array(-1, dim = c(predictors, 2))
  cond_index <- c()
  pk_beta_samples <- c()
  fppi_beta_samples <- c()
  
  for (t in 1:targets) {
    
    # intercept and alpha parameters are response specific
    target_node <- target_nodes[t]
    
    if (endsWith(target_node$personal_intercept, ".rds")) {
      target_intercept <- unlist(readRDS(paste0(graph_dir, "/", target_node$personal_intercept)), use.names=FALSE) # vector(full_sample_size)
      intercept_samples[t,] <- target_intercept[density_samples]
      #print("Found distribution for intercept parameter.")
    } else {
      print("Expected rds-file of density samples. Aborting.")
      return()
    }
    
    if (endsWith(target_node$g_alpha, ".rds")) {
      target_alpha <- unlist(readRDS(paste0(graph_dir, "/", target_node$g_alpha)), use.names=FALSE) # vector(full_sample_size)
      alpha_samples[t,] <- target_alpha[density_samples]
      #print("Found distribution for alpha parameter.")
    } else {
      print("Expected rds-file of density samples. Aborting.")
      return()
    }

    # beta is collected from target/predictor
    for (p in 1:predictors) {
      
      # get edge connection from predictor to target
      beta_edge_id <- reaction_graph[from = predictor_nodes[p], to = target_nodes[t], edges=TRUE]
      
      # get distribution samples related to this edge
      # - distribution-attribute stores the name of the RDS-file that has the sample matrix
      edge_beta_distribution <- readRDS(paste0(graph_dir, "/", E(reaction_graph)[beta_edge_id]$distribution))
      
      # combine betas in one matrix for Stan
      beta_samples[t,p,] <- edge_beta_distribution[density_samples]
    }
  }
  
  # evidence is collected from predictors 
  for (p in 1:predictors) {
    
    # Get soft evidence from predictor random variables
    if (startsWith(predictor_nodes[p]$distribution, "N"))
    {
      # Use normalized distribution values 
      prior_stats[p, 1] <- as.numeric(str_extract(predictor_nodes[p]$normdist,"(?<=N\\().+(?=,)")) # N(*,
      prior_stats[p, 2] <- as.numeric(str_extract(predictor_nodes[p]$normdist,"(?<=,).+(?=\\))")) # ,*)
    } 
    else if (startsWith(predictor_nodes[p]$distribution, "LN")) {
      prior_stats[p, 1] <- as.numeric(str_extract(predictor_nodes[p]$normdist,"(?<=LN\\().+(?=,)")) # LN(*,
      prior_stats[p, 2] <- as.numeric(str_extract(predictor_nodes[p]$normdist,"(?<=,).+(?=\\))")) # ,*)
      
      if (is.na(prior_stats[p, 1]) || is.na(prior_stats[p, 2]))
      {
        print("Nutrient evidence contains NA")
      }
    } 
    else if (startsWith(predictor_nodes[p]$distribution, "B")) {
      prior_stats[p, 1] <- as.numeric(str_extract(predictor_nodes[p]$distribution,"(?<=B\\().+(?=\\))")) # B(*)
      prior_stats[p, 2] <- -1
    }
    
    # Collect indexes of queried nodes (withing prior_stats)
    repeat_without_conditioning <- 0
    
    if (!is.null(queried_nodes)) {
      if (predictor_nodes[p]$name %in% queried_nodes) 
      {
        print(paste0(predictor_nodes[p]$name," is conditional_nutrients[",p,"]"));
        
        cond_index <- c(cond_index,p)
        
        # - and set prior for the queried nodes
        # - these priors override the previously estimated soft evidences
        
        # TODO: Finish this and remove the hack
        
        # prior_stats[i, 1] <- proposal_distribution_params
        
        # HACK: make dynamic
        if (predictor_nodes[p]$name == "kalium")
        {
          prior_stats[p, 1] <- proposal_distribution_params[1]
          prior_stats[p, 2] <- proposal_distribution_params[2]
          
          
          #print(paste0("found prior for kalium: U(",proposal_distribution_params[1],",",proposal_distribution_params[2],")"))
          
        } else if (predictor_nodes[p]$name == "fosfori")  
        {
          prior_stats[p, 1] <- proposal_distribution_params[3]
          prior_stats[p, 2] <- proposal_distribution_params[4]
          
          #print(paste0("found prior for phosphorous: U(",proposal_distribution_params[3],",",proposal_distribution_params[4],")"))
        }
      } 
    } else {
      repeat_without_conditioning <- 1
      cond_index <- c(predictors+100,predictors+100) # some dummy value for Stan
    }
  }
  
  # HACK
  #cond_index <- c(cond_index,-1)
  
  point_est_param <- 0
  if (point_est == TRUE) point_est_param <- 1

  params <- within(list(),
                   {
                     parameter_samples <- length(density_samples) # number of parameter samples that we have
                     p <- predictors
                     v <- targets
                     cond_nodes <- length(cond_index) # size of cond_index array
                     cond_index <- as.vector(cond_index)
                     evidence <- prior_stats
                     
                     intercept <- intercept_samples
                     alpha <- alpha_samples
                     beta <- beta_samples
                     
                     point_est_params <- point_est_param
                     posterior_samples <- posterior_samples
                     repeat_without_conditioning <- repeat_without_conditioning
                     
                     query <- query
                   })
  
  saveRDS(params, "last_params.rds")
  
  rstan_options (auto_write=TRUE)
  options (mc.cores=parallel::detectCores ()) 
  
  query_result <- stan(file=stan_model_file, warmup = 0, iter=iterations, chains=4, algorithm="Fixed_param", data=params, seed=483892929)
  
  return(query_result)
}

##################################################

mebn.personal_prediction <- function(reaction_graph, graph_dir, evidence, stan_model_file)
{
  library(rstan)
  library(igraph)

  predictor_nodes <- V(reaction_graph)[type == 100]
  target_nodes <- V(reaction_graph)[type == 200]
  
  predictors <- length(predictor_nodes)
  targets <- length(target_nodes)
  
  # Collect model parameters from the graph
  full_sample_size <- 4000 # warmup and sampling stored in RDS
  density_samples <- 2001:4000 # number of lastest samples to use 
  
  beta_samples <- array(-1, dim = c(targets, predictors, length(density_samples)))
  intercept_samples <- array(-1, dim = c(targets, length(density_samples)))
  alpha_samples <- array(-1, dim = c(targets, length(density_samples)))
  
  for (t in 1:targets) {
    
    # intercept and alpha parameters are responses specific
    target_node <- target_nodes[t]

    if (endsWith(target_node$personal_intercept, ".rds")) {
      target_intercept <- unlist(readRDS(paste0(graph_dir, "/", target_node$personal_intercept)), use.names=FALSE) # vector(full_sample_size)
      intercept_samples[t,] <- target_intercept[density_samples]
    } else {
      print("Expected rds-file of density samples. Aborting.")
      return()
    }

    if (endsWith(target_node$g_alpha, ".rds")) {
      target_alpha <- unlist(readRDS(paste0(graph_dir, "/", target_node$g_alpha)), use.names=FALSE) # vector(full_sample_size)
      alpha_samples[t,] <- target_alpha[density_samples]
    } else {
      print("Expected rds-file of density samples. Aborting.")
      return()
    }

    # beta parameters are predictor specific
    for (p in 1:predictors) {
      
      # get edge connection from predictor to target
      beta_edge_id <- reaction_graph[from = predictor_nodes[p], to = target_nodes[t], edges=TRUE]
      
      # get distribution samples related to this edge
      # - distribution-attribute stores the name of the RDS-file that has the sample matrix
      edge_beta_distribution <- readRDS(paste0(graph_dir, "/", E(reaction_graph)[beta_edge_id]$distribution))
      
      # combine betas in one matrix for Stan
      beta_samples[t,p,] <- edge_beta_distribution[density_samples]
    }
  }

  params <- within(list(),
                   {
                     parameter_samples <- length(density_samples) # number of parameter samples that we have
                     p <- predictors
                     v <- targets
                     X <- evidence
                     N <- nrow(evidence)
                     intercept <- intercept_samples
                     alpha <- alpha_samples
                     beta <- beta_samples
                   })
  
  #saveRDS(params, "last_params.rds")
  
  rstan_options (auto_write=TRUE)
  options (mc.cores=parallel::detectCores ()) 
  
  prediction <- stan(file=stan_model_file, warmup = 0, iter=2000, chains=4, algorithm="Fixed_param", data=params, seed=483892929)
  
  return(prediction)
}

##################################################

mebn.personal_graph <- function(person_id, reaction_graph, predictor_columns, assumed_targets, local_distributions)
{
  library(igraph)
  library(rstan)
  
  for (c in 1:dim(assumed_targets)[1])
  {
    target_column <- assumed_targets[c,]
    target_name <- as.vector(target_column$Name)
    
    localfit_directory <- local_distributions[local_distributions$Name==target_name,]$modelcache
    
    localfit <- mebn.get_localfit(target_name, localfit_directory)

    # extract personal effects from the local distribution
    pe <- mebn.personal_effects(localfit, person_id)
    
    # - Loop through betas for current target
    predictor_names <- as.vector(predictor_columns$Name)
    
    for (p in 1:length(predictor_names))
    {
      predictor_name <- predictor_names[p]
      
      # Attach the random variable
      # Data values are also stored in edge attributes for easier visualization
      reaction_graph <- reaction_graph + edge(c(predictor_name, target_name),
                                              weight = pe$personal_effect[p], 
                                              value = pe$personal_effect[p], 
                                              value_lCI = pe$personal_effect_lCI[p],
                                              value_uCI = pe$personal_effect_uCI[p],
                                              b = pe$b[p], 
                                              b_lCI = pe$b_lCI[p],
                                              b_uCI = pe$b_uCI[p])

      # Personal effect (beta + b)
      reaction_graph <- reaction_graph + vertex(paste0("personal_", predictor_name, "_", target_name), 
                                                label=paste0("personal_", predictor_name), 
                                                type="personal", color="#AAAAAA", 
                                                value = pe$personal_effect[p], 
                                                value_lCI = pe$personal_effect_lCI[p],
                                                value_uCI = pe$personal_effect_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("personal_", predictor_name, "_", target_name), target_name, shape = "arrow", weight = 1, type = "personal") 

      # Personal variation (b)
      reaction_graph <- reaction_graph + vertex(paste0("b_", predictor_name, "_", target_name), 
                                                label=paste0("b_", predictor_name), 
                                                type="b", color="#AAAAAA", 
                                                value = pe$b[p], 
                                                value_lCI = pe$b_lCI[p],
                                                value_uCI = pe$b_uCI[p],
                                                shape = "circle")
      
      reaction_graph <- reaction_graph + edge(paste0("b_", predictor_name, "_", target_name), target_name, shape = "arrow", weight = 1, type = "b") 
    }
    
  } # loop targets
  
  return(reaction_graph)
}

##################################################

mebn.plot_typical_effects <- function(reaction_graph, top_effects, graph_layout = NULL)
{
  library(igraph)
  
  # Parameter and hyperparameter nodes are removed and visualized otherwise
  visual_graph <- reaction_graph
  
  # Remove edges to latent variables
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="beta"))
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="b_sigma"))
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="b"))
  
  # Remove nodes of latent variable
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="beta"))
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="b_sigma"))
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="b"))
  
  # Filter only the most significant edges having large typical effect or large personal variance
  alledges <- E(visual_graph)
  top_neg_edges <- head(alledges[order(alledges$weight)], top_effects)
  top_pos_edges <- head(alledges[order(-alledges$weight)], top_effects)
  #top_pers_edges <- head(alledges[order(-alledges$b_sigma)], top_effects)
  
  # Comment out this row to see all the connections at the model
  #visual_graph <- delete.edges(visual_graph, alledges[-c(top_neg_edges, top_pos_edges, top_pers_edges)])
  visual_graph <- delete.edges(visual_graph, alledges[-c(top_neg_edges, top_pos_edges)])
  
  # Graph layout
  V(visual_graph)$size = 3 
  
  #E(visual_graph)$weight <- E(visual_graph)$weight * 2
  
  if (is.null(graph_layout))
    graph_layout <- mebn.layout_bipartite_horizontal(visual_graph, V(visual_graph)$type == "100") 
  
  # Align vertex labels according to graph level
  V(visual_graph)[V(visual_graph)$type == "100"]$label.degree = pi # left side
  V(visual_graph)[V(visual_graph)$type == "200"]$label.degree = 0 # right side
  
  # Color and size encoding for edges according to beta coefficient
  E(visual_graph)[E(visual_graph)$weight > 0]$color="#D01C1F"
  E(visual_graph)[E(visual_graph)$weight > 0]$lty=1
  E(visual_graph)[E(visual_graph)$weight < 0]$color="#4B878B"
  E(visual_graph)[E(visual_graph)$weight < 0]$lty=1

  # Black and white
  
  #E(visual_graph)[E(visual_graph)$weight > 0]$color="#444444"
  #E(visual_graph)[E(visual_graph)$weight < 0]$color="#AAAAAA"

  E(visual_graph)$width = abs(E(visual_graph)$weight) * 4

  plot(visual_graph, 
       layout=graph_layout, 
       rescale=TRUE,
       vertex.label.color="black",
       vertex.label.cex=0.7,
       vertex.label.dist=3.8,
       edge.arrow.size=0.5,
       edge.arrow.width=1,
       axes=FALSE,
       margin=0,
       layout.par = par(mar=c(3.0,0,0,0)))

  #axis(1, at = -1:1, labels=c("Nutrients", "Processes and organs", "Personal goals"), cex.axis=0.7)
  #axis(1, at = 1:4)
  
  return(graph_layout)
}

##################################################

mebn.layout_bipartite_horizontal <- function(layout_graph, rank_condition)
{
  require(igraph)
  
  layout <- layout_as_bipartite(layout_graph, rank_condition)

  # - flip layout sideways, from left to right
  gap <- 6
  layout <- cbind(layout[,2]*gap, layout[,1])
  
  return(layout)
}

##################################################
  
mebn.plot_personal_effects <- function(personal_graph, top_effects, graph_layout = NULL, plot_title="")
{
  library(igraph)
  
  # Parameter and hyperparameter nodes are removed and visualized otherwise
  visual_graph <- personal_graph
  
  # Remove edges to latent variables
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="beta"))
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="b_sigma"))
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="b"))
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="personal"))
  
  # Remove nodes of latent variable
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="beta"))
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="b_sigma"))
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="b"))
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="personal"))
  
  # Filter only the most significant edges having large typical effect
  alledges <- E(visual_graph)
  top_neg_edges <- head(alledges[order(alledges$weight)], top_effects)
  top_pos_edges <- head(alledges[order(-alledges$weight)], top_effects)

  # Comment out this row to see all the connections at the model
  visual_graph <- delete.edges(visual_graph, alledges[-c(top_neg_edges, top_pos_edges)])
  
  # Graph layout
  V(visual_graph)$size = 5

  # - put all blood test values in own rank

  if (is.null(graph_layout))
    graph_layout <- mebn.layout_bipartite_horizontal(visual_graph, V(visual_graph)$type == "100") 
  
  # Align vertex labels according to graph level
  V(visual_graph)[V(visual_graph)$type == "100"]$label.degree = pi # left side
  V(visual_graph)[V(visual_graph)$type == "200"]$label.degree = 0 # right side
  
  # Color and size encoding for edges according to beta + b coefficients
  E(visual_graph)[E(visual_graph)$weight > 0]$color="#444444"
#  E(visual_graph)[E(visual_graph)$weight > 0]$lty=1
  E(visual_graph)[E(visual_graph)$weight < 0]$color="#AAAAAA"
#  E(visual_graph)[E(visual_graph)$weight < 0]$lty=5
  E(visual_graph)$width = abs(E(visual_graph)$weight) * 6
  
  plot(visual_graph, 
       layout=graph_layout, 
       rescale=TRUE,
       vertex.label.color="black",
       vertex.label.cex=0.6,
       vertex.label.dist=3.5,
       edge.arrow.size=0.5,
       edge.arrow.width=1,
       curved = 0,
       margin=0,
       layout.par = par(mar=c(0.8,0,0.3,0)))       
    
  return(graph_layout)
}

##################################################

mebn.plot_personal_variations <- function(reaction_graph, top_effects)
{
  library(igraph)
  
  # Parameter and hyperparameter nodes are removed and visualized otherwise
  visual_graph <- reaction_graph
  
  # Remove edges to latent variables
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="beta"))
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="b_sigma"))
  visual_graph <- delete.edges(visual_graph, which(E(visual_graph)$type=="b"))
  
  # Remove nodes of latent variable
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="beta"))
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="b_sigma"))
  visual_graph <- delete.vertices(visual_graph, which(V(visual_graph)$type=="b"))
  
  # Filter only the most significant edges having large typical effect or large personal variance
  alledges <- E(visual_graph)
  top_neg_edges <- head(alledges[order(alledges$weight)], top_effects)
  top_pos_edges <- head(alledges[order(-alledges$weight)], top_effects)
  top_pers_edges <- head(alledges[order(-alledges$b_sigma)], top_effects)
  
  # Comment out this row to see all the connections at the model
  visual_graph <- delete.edges(visual_graph, alledges[-c(top_neg_edges, top_pos_edges, top_pers_edges)])
  
  # Graph layout
  V(visual_graph)$size = 5
  # - put all blood test values in own rank
  bipa_layout <- layout_as_bipartite(visual_graph, types = V(visual_graph)$type == "100")
  # - flip layout sideways, from left to right
  gap <- 4
  bipa_layout <- cbind(bipa_layout[,2]*gap, bipa_layout[,1])
  
  # Align vertex labels according to graph level
  V(visual_graph)[V(visual_graph)$type == "100"]$label.degree = pi # left side
  V(visual_graph)[V(visual_graph)$type == "200"]$label.degree = 0 # right side  
  
  # Color and size encoding for edges according to beta coefficient
  E(visual_graph)$color="gray"
  E(visual_graph)$width = abs(E(visual_graph)$b_sigma) * 10
  
  plot(visual_graph, 
       layout=bipa_layout, 
       rescale=TRUE,
       vertex.label.color="black",
       vertex.label.cex=0.7,
       vertex.label.dist=3.8,
       edge.arrow.size=0.5,
       edge.arrow.width=1,
       axes=FALSE,
       margin=0,
       layout.par = par(mar=c(0,0,0,0)))
}

##################################################

mebn.plot_clusters <- function(cluster_data, cluster_spread, clusters_index, assumedpredictors, assumedtargets, keep_only_effects, feature_index, sort_by_amount = FALSE)
{
  # Number of items in all clusters
  all_items <- sum(cluster_spread$Freq)
  
  # Build effect names
  cluster_data$predictor <- rep(assumedpredictors[feature_index,]$Description,nrow(assumedtargets))
  cl <- rep(1,length(feature_index)) # number of predictors
  t <- nrow(assumedtargets) # number of targets
  t_idx <- c()              # predictors x targets
  for (i in seq(0,t-1))
  {
    t_idx <- c(t_idx, cl+i)      
  }
  cluster_data$response <- assumedtargets$Description[t_idx]
  effect_levels <- paste0(cluster_data$predictor," -> ", cluster_data$response)
  cluster_data$effect <- factor(effect_levels, levels=effect_levels)
  
  # Plot only those effects that showed previously personal variance 
  cluster_data.filtered <- cluster_data
  
  if (!is.null(keep_only_effects)) cluster_data.filtered <- cluster_data[cluster_data$effect %in% keep_only_effects$effect,]
  
  # Prepare data from plotting 
  plot_data <- cluster_data.filtered[c("effect")]

  i <- cluster_index[1]
  plot_data$amount <- cluster_data.filtered[as.character(i)][,]
  #plot_data$amount <- cluster_data.filtered$'1'
  plot_data$cluster <- i
  
  cluster_labels <- c()
  
  # Items in cluster
  cluster_items <- cluster_spread[cluster_spread[1]== i,]$Freq
  cluster_items_perc <- round((cluster_items / all_items)*100, 1)
  
  cluster_labels <- c(cluster_labels, paste0("cluster ", i, " (", cluster_items, "/", all_items, ", ", cluster_items_perc, "%)"))
  
  for (i in cluster_index[-c(1)])
  {
    temp_data <- cluster_data.filtered[c("effect")]
    temp_data$amount <- cluster_data.filtered[as.character(i)][,]
    temp_data$cluster <- i
    
    cluster_items <- cluster_spread[cluster_spread[1]== i,]$Freq
    cluster_items_perc <- round((cluster_items / all_items)*100, 1)
    
    cluster_labels <- c(cluster_labels, paste0("cluster ", i, " (", cluster_items, "/", all_items, ", ", cluster_items_perc, "%)"))
    
    plot_data <- rbind(plot_data, temp_data)
  }
  
  names(cluster_labels) <- cluster_index

  plot_data$below_above <- ifelse(plot_data$amount < 0, "below", "above")
  
  if (sort_by_amount == TRUE) {
    ggplot(plot_data, aes(x=reorder(effect, amount), y=amount)) + 
      geom_bar(stat='identity', aes(fill=below_above), width=.5, show.legend = FALSE) +
      scale_fill_manual(values = c("#333333", "#999999")) +
      coord_flip() +
      facet_wrap(~cluster, labeller = labeller(cluster = cluster_labels)) +
      theme_bw() +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(), text=element_text(size=9)) 
  } else {  
    ggplot(plot_data, aes(x=reorder(effect, amount), y=amount)) + 
      geom_bar(stat='identity', aes(fill=below_above), width=.5, show.legend = FALSE) +
      scale_fill_manual(values = c("#333333", "#999999")) +
      coord_flip() +
      facet_wrap(~cluster, labeller = labeller(cluster = cluster_labels)) +
      theme_bw() +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(), text=element_text(size=9)) 
  }  
}

##################################################

mebn.compare_typicals <- function(bn1, bn2)
{
  library(igraph)
  
  # - find beta nodes of both normal and gamma distributions 
  model1_nodes <- V(bn1)
  m1_beta <- model1_nodes[model1_nodes$type=="beta"]
  
  model2_nodes <- V(bn2)
  m2_beta <- model2_nodes[model2_nodes$type=="beta"]
  
  # - construct a table for comparing the estimates
  typical_effects<-data.frame(matrix(NA, nrow=length(m1_beta), ncol=0))
  
  typical_effects$effect <- unlist(lapply(strsplit(gsub("beta_","", m1_beta$name), "_"), function(x) paste0(toString(datadesc[datadesc$Name==x[1],]$Description)," -> ", toString(datadesc[datadesc$Name==x[2],]$Description))))
  
  typical_effects$model1 <- round(m1_beta$value,6)
  typical_effects$model1_lCI <- round(m1_beta$value_lCI,6)
  typical_effects$model1_hCI <- round(m1_beta$value_uCI,6)
  typical_effects$model2 <- round(m2_beta$value,6)
  typical_effects$model2_lCI <- round(m2_beta$value_lCI,6)
  typical_effects$model2_hCI <- round(m2_beta$value_uCI,6)
  
  return(typical_effects)
}

##########################################################

mebn.GetEvidence <- function(reaction_graph, queried_nodes, point_est = "mean")
{
  library(igraph)
  library(stringr)
  
  nodes <- V(reaction_graph)[type == 100]
  predictors <- length(nodes)
  
  evidence <- array(-1, dim = c(predictors, 2))
  cond_index <- c()
  
  for (p in 1:predictors) {
    
    # Get soft evidence from predictor random variables
    if (startsWith(nodes[p]$distribution, "N"))
    {
      # Use normalized distribution values 
      evidence[p, 1] <- as.numeric(str_extract(nodes[p]$normdist,"(?<=N\\().+(?=,)")) # N(*,
      evidence[p, 2] <- as.numeric(str_extract(nodes[p]$normdist,"(?<=,).+(?=\\))")) # ,*)
    } 
    else if (startsWith(nodes[p]$distribution, "LN")) {
      evidence[p, 1] <- as.numeric(str_extract(nodes[p]$normdist,"(?<=LN\\().+(?=,)")) # LN(*,
      evidence[p, 2] <- as.numeric(str_extract(nodes[p]$normdist,"(?<=,).+(?=\\))")) # ,*)
      
      if (is.na(evidence[p, 1]) || is.na(evidence))
      {
        print("Nutrient evidence contains NA")
      }
    } 
    else if (startsWith(nodes[p]$distribution, "B")) {
      evidence[p, 1] <- as.numeric(str_extract(nodes[p]$distribution,"(?<=B\\().+(?=\\))")) # B(*)
      evidence[p, 2] <- -1
    }
    
    # Collect indexes of queried nodes (withing prior_stats)
    if (nodes[p]$name %in% queried_nodes) 
    {
      #print(paste0(nodes[p]$name," is node[",p,"]"));
      cond_index <- c(cond_index,p)
    }
    
    if (point_est == "CI5")
    {
      evidence_point <- apply(evidence, 1, function(x) if (x[2] == -1) {x[1]} else {qnorm(p = 0.05, mean = x[1], sd = x[2], lower.tail = TRUE, log.p = FALSE)})
    }
    else if (point_est == "CI95")
    {
      evidence_point <- apply(evidence, 1, function(x) if (x[2] == -1) {x[1]} else {qnorm(p = 0.95, mean = x[1], sd = x[2], lower.tail = TRUE, log.p = FALSE)})    
    }
    else
    {
      evidence_point <- evidence[,1]    
    }

  }
  
  result <- within(list(),
                   {
                     evidence <- evidence
                     evidence_point <- evidence_point
                     cond_index <- cond_index
                   })
  
  return(result)
}

mebn.GetIntercept <- function(reaction_graph, graph_dir, point_est = "mean")
{
  chains <- 4
  density_chain_samples <- 1000
  density_samples <- density_chain_samples * chains  # sampling stored in RDS for all the chains (4*1000)
  
  target_nodes <- V(reaction_graph)[type == 200]
  
  targets <- length(target_nodes)
  intercept_samples <- array(-1, dim = c(targets, density_samples))
  intercept_point <- array(-1, dim = c(targets))
  
  for (t in 1:targets) {
    
    # intercept parameters are response specific
    target_node <- target_nodes[t]
    rds_dir <- gsub("//", "/", paste0(graph_dir, "/", target_node$personal_intercept))
    
    if (endsWith(target_node$personal_intercept, ".rds")) {
      target_intercept <- unlist(readRDS(rds_dir), use.names=FALSE) # vector(full_sample_size)
      
      intercept_samples[t,] <- target_intercept
      #print("Found distribution for intercept parameter.")
      
      if (point_est == "median")
      {
        intercept_point[t] <- median(intercept_samples[t,])   
      } 
      else if (point_est == "mean")
      {
        intercept_point[t] <- mean(intercept_samples[t,])   
      }
      else if (point_est == "CI5")
      {
        intercept_point[t] <- quantile(intercept_samples[t,], probs=0.05)   
      } 
      else if (point_est == "CI95")
      {
        intercept_point[t] <- quantile(intercept_samples[t,], probs=0.95)   
      }
    
    } else {
      print("Expected rds-file of density samples. Aborting.")
      return()
    }
  }
  
  return(intercept_point)
}

mebn.GetBeta <- function(reaction_graph, graph_dir, point_est = "mean")
{
  chains <- 4
  density_chain_samples <- 1000
  density_samples <- density_chain_samples * chains  # sampling stored in RDS for all the chains (4*1000)
  
  predictor_nodes <- V(reaction_graph)[type == 100]
  target_nodes <- V(reaction_graph)[type == 200]
  
  predictors <- length(predictor_nodes)
  targets <- length(target_nodes)
  
  beta_samples <- array(-1, dim = c(targets, predictors, density_samples))
  beta_point <- array(-1, dim = c(targets, predictors))
  
  for (t in 1:targets) {
    
    # intercept parameters are response specific
    target_node <- target_nodes[t]
    
    # beta is collected from target/predictor
    for (p in 1:predictors) {
      
      # get edge connection from predictor to target
      beta_edge_id <- reaction_graph[from = predictor_nodes[p], to = target_nodes[t], edges=TRUE]
      
      # get distribution samples related to this edge
      # - distribution-attribute stores the name of the RDS-file that has the sample matrix
      edge_beta_distribution <- readRDS(paste0(graph_dir, "/", E(reaction_graph)[beta_edge_id]$distribution))
      
      # combine betas in one matrix
      beta_samples[t,p,] <- edge_beta_distribution
      
      # point estimates
      if (point_est == "median")
      {
        beta_point[t,p] <- median(beta_samples[t,p,])   
      } 
      else if (point_est == "mean")
      {
        beta_point[t,p] <- mean(beta_samples[t,p,])   
      }
      else if (point_est == "mode")
      {
        beta_point[t,p] <- mebn.getmode(beta_samples[t,p,])   
      }
      else if (point_est == "CI5")
      {
        beta_point[t,p] <- quantile(beta_samples[t,p,], probs=0.05)   
      } 
      else if (point_est == "CI95")
      {
        beta_point[t,p] <- quantile(beta_samples[t,p,], probs=0.95)   
      } 
    }      
  }  
  
  return(beta_point)
}


mebn.Query <- function(reaction_graph, graph_dir, query, queried_nodes, proposal_limits, conc_lower_limits, conc_upper_limits, stan_model_file, X_point_est = "mean", beta_point_est = "mean", param_point_est = "mean", posterior_samples = 100, X_sd_coef, repeat_only = 0, condition_in_repeat = 0)
{
  library(rstan)
  library(igraph)
  library(stringr)
  
  predictor_nodes <- V(reaction_graph)[type == 100]
  target_nodes <- V(reaction_graph)[type == 200]
  
  predictors <- length(predictor_nodes)
  targets <- length(target_nodes)
  
  # Collect model parameters from the graph
  chains <- 4
  density_chain_samples <- 1000
  density_samples <- density_chain_samples * chains  # sampling stored in RDS for all the chains (4*1000)
  
  beta_samples <- array(-1, dim = c(targets, predictors, density_samples))
  intercept_samples <- array(-1, dim = c(targets, density_samples))
  alpha_samples <- array(-1, dim = c(targets, density_samples))
  
  beta_point <- array(-1, dim = c(targets, predictors))
  intercept_point <- array(-1, dim = c(targets))
  alpha_point <- array(-1, dim = c(targets))
  
  predictor_evidence <- array(-1, dim = c(predictors, 2))
  cond_index <- c()
  
  for (t in 1:targets) {
    
    # intercept and alpha parameters are response specific
    target_node <- target_nodes[t]
    
    if (endsWith(target_node$personal_intercept, ".rds")) {
      target_intercept <- unlist(readRDS(paste0(graph_dir, "/", target_node$personal_intercept)), use.names=FALSE) # vector(full_sample_size)
      intercept_samples[t,] <- target_intercept
      #print("Found distribution for intercept parameter.")
    } else {
      print("Expected rds-file of density samples. Aborting.")
      return()
    }

    if (endsWith(target_node$g_alpha, ".rds")) {
      target_alpha <- unlist(readRDS(paste0(graph_dir, "/", target_node$g_alpha)), use.names=FALSE) # vector(full_sample_size)
      alpha_samples[t,] <- target_alpha
      #print("Found distribution for alpha parameter.")
    } else {
      print("Expected rds-file of density samples. Aborting.")
      return()
    }

    # beta is collected from target/predictor
    for (p in 1:predictors) {
      
      # get edge connection from predictor to target
      beta_edge_id <- reaction_graph[from = predictor_nodes[p], to = target_nodes[t], edges=TRUE]
      
      # get distribution samples related to this edge
      # - distribution-attribute stores the name of the RDS-file that has the sample matrix
      edge_beta_distribution <- readRDS(paste0(graph_dir, "/", E(reaction_graph)[beta_edge_id]$distribution))
      
      # combine betas in one matrix for Stan
      beta_samples[t,p,] <- edge_beta_distribution
      
      # point estimates
      if (beta_point_est == "median")
      {
        beta_point[t,p] <- median(beta_samples[t,p,])   
      } 
      else if (beta_point_est == "mean")
      {
        beta_point[t,p] <- mean(beta_samples[t,p,])   
      } 
      else if (beta_point_est == "mode")
      {
        beta_point[t,p] <- mebn.get_mode(beta_samples[t,p,])   
      } 
      else if (beta_point_est == "CI5")
      {
        beta_point[t,p] <- quantile(beta_samples[t,p,], probs=0.05)   
        print("beta_point_est: CI5")
      } 
      else if (beta_point_est == "CI95")
      {
        beta_point[t,p] <- quantile(beta_samples[t,p,], probs=0.95)   
        print("beta_point_est: CI95")
      } 
    }
    
    # point estimates
    if (param_point_est == "median")
    {
      intercept_point[t] <- median(intercept_samples[t,])   
      alpha_point[t] <- median(alpha_samples[t,]) 
    } else
    if (param_point_est == "mean")
    {
      intercept_point[t] <- mean(intercept_samples[t,])   
      alpha_point[t] <- mean(alpha_samples[t,]) 
    }
  }
  
  # evidence is collected from predictors 
  for (p in 1:predictors) {
    
    # Get soft evidence from predictor random variables
    if (startsWith(predictor_nodes[p]$distribution, "N"))
    {
      # Use normalized distribution values 
      predictor_evidence[p, 1] <- as.numeric(str_extract(predictor_nodes[p]$normdist,"(?<=N\\().+(?=,)")) # N(*,
      predictor_evidence[p, 2] <- as.numeric(str_extract(predictor_nodes[p]$normdist,"(?<=,).+(?=\\))")) # ,*)
    } 
    else if (startsWith(predictor_nodes[p]$distribution, "LN")) {
      predictor_evidence[p, 1] <- as.numeric(str_extract(predictor_nodes[p]$normdist,"(?<=LN\\().+(?=,)")) # LN(*,
      predictor_evidence[p, 2] <- as.numeric(str_extract(predictor_nodes[p]$normdist,"(?<=,).+(?=\\))")) # ,*)
      
      if (is.na(predictor_evidence[p, 1]) || is.na(predictor_evidence))
      {
        print("Nutrient evidence contains NA")
      }
    } 
    else if (startsWith(predictor_nodes[p]$distribution, "B")) {
      predictor_evidence[p, 1] <- as.numeric(str_extract(predictor_nodes[p]$distribution,"(?<=B\\().+(?=\\))")) # B(*)
      predictor_evidence[p, 2] <- -1
    }
    
    # Collect indexes of queried nodes (withing prior_stats)
    if (predictor_nodes[p]$name %in% queried_nodes) 
    {
      print(paste0(predictor_nodes[p]$name," is conditional_nutrients[",p,"]"));
      
      cond_index <- c(cond_index,p)
    }
  }
  
  offset_value <- 1000
  
  if (repeat_only != 1) {
    X_beta_point <- beta_point[,-cond_index]
    Q_beta_point <- beta_point[,cond_index]
    X_evidence <- predictor_evidence[-cond_index,]
    
    p <- predictors - length(cond_index)
  }
  else
  {
    X_beta_point <- beta_point
    Q_beta_point <- beta_point[,cond_index] # is not used
    X_evidence <- predictor_evidence
    
    p <- predictors 
  }
  
  if (X_point_est == "mean")
  {
    X_evidence_point <- X_evidence[,1]    
    print("X_point_est: mean")
  }
  else
  if (X_point_est == "CI5")
  {
    X_evidence_point <- apply(X_evidence, 1, function(x) if (x[2] == -1) {x[1]} else {qnorm(p = 0.05, mean = x[1], sd = x[2], lower.tail = TRUE, log.p = FALSE)})
    print("X_point_est: CI5")
  }
  else
  if (X_point_est == "CI95")
  {
    X_evidence_point <- apply(X_evidence, 1, function(x) if (x[2] == -1) {x[1]} else {qnorm(p = 0.95, mean = x[1], sd = x[2], lower.tail = TRUE, log.p = FALSE)})    
    print("X_point_est: CI95")
  }
  
  #print(X_evidence)
  
  params <- within(list(),
                   {
                     responses <- targets
                     p <- p
                     r <- length(cond_index)
                     
                     proposal_limits <- proposal_limits
                     Y_lower_limits <- conc_lower_limits
                     Y_upper_limits <- conc_upper_limits
                     X_evidence <- X_evidence
                     X_evidence_point <- X_evidence_point
                     
                     intercept_point <- intercept_point
                     alpha_point <- alpha_point
                     X_beta_point <- X_beta_point
                     Q_beta_point <- Q_beta_point
                     linear_transformation <- offset_value
                     
                     posterior_samples <- posterior_samples
                     X_sd_coef <- X_sd_coef
                     
                     repeat_only <- repeat_only
                     condition_in_repeat <- condition_in_repeat
                     Q_index <- cond_index                     
                   })
  
  rstan_options (auto_write=TRUE)
  options (mc.cores=parallel::detectCores ())
  
  query_result <- stan(file=stan_model_file, warmup = 1000, iter=4000, chains=4, chain_id=1L, algorithm="NUTS", control = list(adapt_delta = 0.99, max_treedepth = 15), data=params, seed=483892929)
  #query_result <- stan(file=stan_model_file, warmup = 1000, iter=5000, chains=4, chain_id=1L, algorithm="NUTS", control = list(adapt_delta = 0.99, max_treedepth = 15), data=params, seed=483892929)
  
  #m <- stan_model(file=stan_model_file)
  #query_result <- optimizing(m, data=params, seed=483892929, verbose=TRUE, init = 0)
  
  # use this if more parameters are needed in result
  # result <- within(list(),
  #                  {
  #                    alpha_point <- alpha_point
  #                    result <- query_result
  #                  })
                     
  return(query_result)
}


########################



