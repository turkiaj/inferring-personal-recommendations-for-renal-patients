
data { 
  int<lower=0> N;   // number of observations
  int<lower=0> NH;  // number of holdout observations
  int<lower=1> p;   // number of predictors
  int<lower=1> v;   // number of responses
  int<lower=1> n_g; // number of groups in data 
  int<lower=1> n_s; // number of subjects in data
  int<lower=1> k;   // number of personal predictors (for both group and subject level)
  int<lower=1,upper=n_g> group[N];   //group indicator
  int<lower=1,upper=n_s> subject[N]; //subject indicator
  int<lower=1,upper=n_g> group_for_subject[n_s]; // lookup-table for group/subject hiearchy
  matrix[N,p] X;    // fixed-effect design matrix
  matrix[N,k] Z;    // group/subject-level design matrix
  vector[N] Y[v];   // responses
  real offset;

  // indicates if observation is used for training (0) or prediction (1)
  int<lower=0,upper=1> holdout[N]; 
} 

transformed data { 
  
  matrix[v*(N-NH),v*(p-1)] XM_t; // common training input for all responses
  matrix[v*(N-NH),v*k] ZM_t;     // personal training input for all responses
  vector[v*(N-NH)] Y_t;          // stacked training responses from all responses
  int t=1;                       // index
  matrix[v*(N-NH),v] I;          // identity matrix, block-diagonal matching stacked data
  vector[v*(N-NH)] linear_transformation;  // transforms all data to positive region
  
  // QR reparameteratization of the model matrix
  // - reparamization matrices are stacked like X_t
  
  matrix[v*(N-NH), v*(p-1)] Q_ast;
  matrix[v*(p-1), v*(p-1)] R_ast;
  matrix[v*(p-1), v*(p-1)] R_ast_inverse;

  I = rep_matrix(0,v*(N-NH),v);

  // X_t and Z_t are stacked vectors from all responses
  // - fill X_t and Z_t initially with 0
  XM_t = rep_matrix(0,v*(N-NH),v*(p-1));
  ZM_t = rep_matrix(0,v*(N-NH),v*k);
  
  for (m in 1:v)
  {
    for (n in 1:N)
    {
      if (holdout[n] == 0)
      {     
        // two loops create block-diagonal 1-matrix (N-NH x v)
        I[t,m] = 1;

        // XM_t is block diagonal regarding to response inputs
        // - the intercept is removed from the model matrix 
        XM_t[t,(m-1)*(p-1)+1:m*(p-1)] = X[n,2:p];

        // Same block dialonal design matrix is used for group and subject levels 
        ZM_t[t,(m-1)*k+1:m*k] = Z[n,1:k];
        
        // linear transformation is applied for true responses also
        // this is transformed back in final results
        Y_t[t] = Y[m,n] + offset;

        t += 1;
      }
    }
  }
  
  linear_transformation =  I * rep_vector(offset, v);

  // thin and scale the QR decomposition
  Q_ast = qr_thin_Q(XM_t) * sqrt(v*(N-NH) - 1);
  R_ast = qr_thin_R(XM_t) / sqrt(v*(N-NH) - 1);
  R_ast_inverse = inverse(R_ast);  
}

parameters { 
  // all parameters are vectorized based on response, 1..v
  // except personal effects are stacked in same matrices for estimating cross-covariance
  
  vector[v] beta_Intercept;         // intercept 
  //vector[v*(p-1)] beta;           // poulation-level effects (fixed effects) are QR decomposed..
  vector[v*(p-1)] theta_q;          // coefficients of Q_ast
  
  // Group and subject level random effects have their of Choleskies with same size

  cholesky_factor_corr[v*k] L_g;   // Cholesky factor of group ranef corr matrix
  vector<lower=0>[v*k] stacked_sigma_b_g;  // group-level random-effect standard deviations
  vector[v*k] z_g[n_g];            // unscaled group-level effects
  
  cholesky_factor_corr[v*k] L_s;   // Cholesky factor of subject ranef corr matrix
  vector<lower=0>[v*k] stacked_sigma_b_s;  // subject-level random-effect standard deviations
  vector[v*k] z_s[n_s];            // unscaled subject-level effects

  vector[v] g_log_alpha;           // alpha (shape) parameter of the gamma distribution
}

transformed parameters {
  vector[v] g_alpha;                // alpha (shape) parameter of each v gamma distribution
  //real<lower=0> sigma_e[v];       // residual standard deviations for each distribution v 
  vector[v*k] b_stack[n_s];         // subject level effects stacked in one vector from each response
  vector[v*k] g_stack[n_g];         // group level effects stacked in one vector from each response

  // Premultiply diagonal matrix [sigma_b] with the Cholesky decomposition L of
  // the correlation matrix Sigma_b to get variance-covariance matrix of group-level effects

  // local scope for Lambdas
  {
    // diag(sigma_b) * L
    matrix[v*k, v*k] Lambda_g;       // Tau * Cholesky decomposition
    matrix[v*k, v*k] Lambda_s;       // for groups and subjects
    
    Lambda_g = diag_pre_multiply(stacked_sigma_b_g, L_g); 
    Lambda_s = diag_pre_multiply(stacked_sigma_b_s, L_s); 
  
    for(j in 1:n_s) 
      b_stack[j] =  Lambda_s * z_s[j];

    for(j in 1:n_g) 
      g_stack[j] =  Lambda_g * z_g[j];
   }
  
  // - log transform alpha parameter to keep it positive
  g_alpha = exp(g_log_alpha);

  // estimate of variance 
  // (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4024993/)
  //sigma_e = log(1 ./ g_alpha + 1);
}

model {
  
    beta_Intercept ~ cauchy(0,10); // prior for the intercept following Gelman 2008
    
    stacked_sigma_b_g ~ student_t(3, 0, 10);
    stacked_sigma_b_s ~ student_t(3, 0, 10);
    
    for (j in 1:n_s) {
      z_s[j] ~ normal(0,1);
    }

    for (j in 1:n_g) {
      z_g[j] ~ normal(0,1);
    }
    
    // - LJK prior for personal effect Cholesky factor
    L_s ~ lkj_corr_cholesky(1);
    L_g ~ lkj_corr_cholesky(1);

  // brackets introduce a new scope for local variables that are not published in model
  {
    vector[v*(N-NH)] mu;             // expected value 
    vector[v*(N-NH)] g_beta;         // beta (rate) of Gamma distribution
    int vn = 1;
    
    mu = I * beta_Intercept + linear_transformation + Q_ast * theta_q;
    
    for (m in 1:v)
    {
      for (n in 1:N-NH)
      {
        // - row n is picked from ZM_t and multiplied with column of coefs from a patient
         mu[vn] = mu[vn] + ZM_t[vn] * g_stack[group[n]] + ZM_t[vn] * b_stack[subject[n]];
         vn += 1;
      }
    }

    // identity link
    g_beta = I * g_alpha ./ mu;

    Y_t ~ gamma(I * g_alpha, g_beta);

    // target += gamma_lpdf(Y_t[i] | g_alpha[i], g_beta);
  }
  
  //print("g_alpha:");
  //print(g_alpha[1]);
  //print(g_alpha[2]);
  //print(g_alpha[3]);

}

generated quantities {

  vector[v*(p-1)] beta_stack;         // poulation-level effects (fixed effects)
  corr_matrix[v*k] C_g;               // correlation matrix of group level effects 
  corr_matrix[v*k] C_s;               // correlation matrix of subject level effects 
  vector[N-NH] Y_rep[v];              // repeated response
  vector[p-1] beta[v];                // poulation-level effects (fixed effects)
  vector[k] g[n_g,v];
  vector[k] b[n_s,v];
  vector[k] group_effects;
  vector[k-1] group_effect[n_g,v];
  real group_intercept[n_g,v];
  vector[k-1] personal_effect[n_s,v];
  real personal_intercept[n_s,v];
  vector[k] sigma_b_g[v];             // unstacked
  vector[k] sigma_b_s[v];             // unstacked
  
  // vector[N-NH] log_lik[v];       // log-likelihood for LOO

  // Correlation matrix of random-effects, C = L'L
  C_g = multiply_lower_tri_self_transpose(L_g);
  C_s = multiply_lower_tri_self_transpose(L_s);

  // Posterior predictive distribution for model checking

  // Unstack Y_rep to separate columns
  {
    real Y_rep_stack[v*(N-NH)];       // stacked training responses from all responses
    vector[v*(N-NH)] mu_hat;          // expected value 
    vector[v*(N-NH)] g_beta_hat;      // beta (rate) of Gamma distribution
    int vn = 1;

    beta_stack = R_ast_inverse * theta_q;     // coefficients on x

    mu_hat = I * beta_Intercept + linear_transformation + XM_t * beta_stack;

    for (m in 1:v)
    {
      for (n in 1:N-NH)
      {
        // - row n is picked from ZM_t and multiplied with column of coefs from a patient
         mu_hat[vn] = mu_hat[vn] + ZM_t[vn] * g_stack[group[n]] + ZM_t[vn] * b_stack[subject[n]];
         vn += 1;
      }
    }
       
    // identity link
    g_beta_hat = I * g_alpha ./ mu_hat;
  
    Y_rep_stack = gamma_rng(I * g_alpha, g_beta_hat);
    
    // transform repeated values back to original intercept
    for (m in 1:v)
    {
      beta[m] = beta_stack[(m-1)*(p-1)+1:m*(p-1)];

      for (n in 1:N-NH)
      {
          Y_rep[m,n] = Y_rep_stack[(m-1)*(N-NH)+n] - offset;
      }
      
      // unstack group and subject variations
      sigma_b_g[m] = stacked_sigma_b_g[(m-1)*k+1:m*k];
      sigma_b_s[m] = stacked_sigma_b_s[(m-1)*k+1:m*k];
      
      // Finally, sample personal effects for each nutrient
      
      // - unstack group effects
      for (i in 1:n_g) 
      {
        g[i,m] = g_stack[i, (m-1)*k+1:m*k];
        
        // group intercept
        group_intercept[i,m] = beta_Intercept[m] + g[i,m][1];
      
        // beta vector does not include intercept, b is also sliced not to include it
        group_effect[i,m] = beta[m] + g[i,m][2:k];
      }

      // loop all subjects
      for (j in 1:n_s)
      {
        // - unstack subject effects
        b[j,m] = b_stack[j, (m-1)*k+1:m*k];
        
        // - find group effects for this subject
        group_effects = g[group_for_subject[j],m];
        
        // personal intercept
        personal_intercept[j,m] = beta_Intercept[m] + group_effects[1] + b[j,m][1];
      
        // beta vector does not include intercept, b is also sliced not to include it
        personal_effect[j,m] = beta[m] + group_effects[2:k] + b[j,m][2:k];
      }

    }
  }
}


