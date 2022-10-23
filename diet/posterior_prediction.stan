
data {
    int<lower=1> N;      // number of observations
    int p;               // number of predictors
    int<lower=1> v;      // number of responses (columns in parameters)
    matrix[N,p] X;       // input
    
    // parameters
    int parameter_samples;   // number of parameter samples

    vector[parameter_samples] intercept[v];
    matrix[p, parameter_samples] beta[v];
    vector[parameter_samples] alpha[v];
}

model {
}

generated quantities {

    real posterior[v];

    {
      real offset = 1000;
      real mu;
      real a;

      // in every iteration draw all prior samples
      for (i in 1:v) {
        for (n in 1:N) {
          for (j in 1:p) {
            for (s in 1:parameter_samples) {
          
              // expected value
              mu = intercept[i,s] + X[n,j] * beta[i,j,s] + offset;
              a = alpha[i,s];
  
              posterior[i] = gamma_rng(a, a / mu) - offset;
            }
          }
        }
      }
    }
}
