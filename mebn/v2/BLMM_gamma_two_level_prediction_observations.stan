
data {
  int<lower=0> N;   // number of observations
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
  
  vector[v] beta_Intercept;       // common intercept 
  vector[v*(p-1)] beta_stack;     // poulation-level effects (fixed effects)

  vector[v*k] g_stack[n_g];       // group level effects stacked in one vector from each response
  vector[v] g_alpha;              // alpha (shape) parameter of each v gamma distribution
} 

transformed data { 
  
  matrix[v*(N),v*(p-1)] XM_t; // common training input for all responses
  matrix[v*(N),v*k] ZM_t;     // personal training input for all responses
  vector[v*(N)] Y_t;          // stacked training responses from all responses
  int t=1;                    // index
  matrix[v*(N),v] I;          // identity matrix, block-diagonal matching stacked data
  vector[v*(N)] linear_transformation;  // transforms all data to positive region
  
  I = rep_matrix(0,v*(N),v);

  // X_t and Z_t are stacked vectors from all responses
  // - fill X_t and Z_t initially with 0
  XM_t = rep_matrix(0,v*(N),v*(p-1));
  ZM_t = rep_matrix(0,v*(N),v*k);
  
  for (m in 1:v)
  {
    for (n in 1:N)
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
  
  linear_transformation =  I * rep_vector(offset, v);
}

parameters { 
  // Only subject-level effects are predicted
  vector[v*k] z_s[n_s];            // unscaled subject-level effects
  vector[v*k] b_stack[n_s];        // subject level effects stacked in one vector from each response
}

model {
  
  {
    vector[v*N] mu;             // expected value 
    vector[v*N] g_beta;         // beta (rate) of Gamma distribution
    int vn = 1;
    
    // b_stack is only parameter, other variables are fixed data
    
    mu = I * beta_Intercept + linear_transformation + XM_t * beta_stack;
    
    for (m in 1:v)
    {
      for (n in 1:N)
      {
        // - row n is picked from ZM_t and multiplied with column of coefs from a patient
         mu[vn] = mu[vn] + ZM_t[vn] * g_stack[group[n]] + ZM_t[vn] * b_stack[subject[n]];
         vn += 1;
      }
    }

    // identity link
    g_beta = I * g_alpha ./ mu;

    Y_t ~ gamma(I * g_alpha, g_beta);
  }
  
}

generated quantities {

  vector[N] Y_rep[v];                 // repeated response
  vector[k] b[n_s,v];
  vector[k-1] personal_effect[n_s,v];
  real personal_intercept[n_s,v];

  // Unstack Y_rep to separate columns
  {
    vector[p-1] beta[v];                // population-level effects (fixed effects)
    vector[k] g[n_g,v];
    vector[k] group_effects;
    vector[k-1] group_effect[n_g,v];
    real group_intercept[n_g,v];

    real Y_rep_stack[v*N];       // stacked training responses from all responses
    vector[v*N] mu_hat;          // expected value 
    vector[v*N] g_beta_hat;      // beta (rate) of Gamma distribution
    int vn = 1;

    mu_hat = I * beta_Intercept + linear_transformation + XM_t * beta_stack;

    for (m in 1:v)
    {
      for (n in 1:N)
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

      for (n in 1:N)
      {
          Y_rep[m,n] = Y_rep_stack[(m-1)*(N)+n] - offset;
      }
      
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


