library( dagitty )
library( dplyr )
library( rethinking )
library( tibble )
library( gt )
library( gto )
library( officer )
library( ggplot2 )
library( ggdist )
options( digits = 8 )

# Load data ----
load( "Data_Kinship_drives_social_networks_even_after_market_integration.RData")

#1 DAG ----
Dag <- 
  dagitty( 'dag {
bb="-0.5,-0.5,0.5,0.5"
Age [pos="0.000,-0.350"]
Density [outcome,pos="0.000,0.320"]
Gender [pos="0.160,-0.250"]
Market [pos="0.210,0.050"]
Residence [pos="-0.160,-0.250"]
U_Individual [latent,pos="-0.210,0.150"]
U_Village [latent,pos="-0.280,-0.050"]
Age -> Density
Age -> Market
Age -> Residence
Gender -> Density
Gender -> Market
Market -> Density
Residence -> Density
Residence -> Market
U_Individual -> Density
U_Village -> Density
}
')

drawdag( Dag , radius = 8.8 )

#2 Data preparation ----
# Standardise age 
# Natural logarithmic transformation and standardisation of market
Ego.data <- Ego.data %>% 
  mutate( Age_std = ( Age - mean( Age ) )/sd( Age ) ,
          Age_std_2 = Age_std^2 , 
          Market_log = log( Market + 1 ) , 
          Market_log_std = ( Market_log - mean( Market_log ) )/sd( Market_log ) )

#3 Models ----
##3.1 Age ----
# Other variable should be included into models
adjustmentSets( Dag , 
                exposure = "Age" , 
                outcome = "Density" ,
                effect = "total" )

###3.1.1 Models of all density ----
Model_GA_all_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_all = as.integer( Tie_all ) , 
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Age2 = as.numeric( sprintf( "%0.4f" , Age_std_2 ) ) , 
  Gender = as.integer( Gender ) , 
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

{set.seed(123)
  Model_GA_all <- ulam(
    alist(
      Tie_all ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + bGA2[Gender] * Age2 +  
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      bGA2[Gender] ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 )
    ) , data = Model_GA_all_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}

Pre_model_GA_all <- precis( Model_GA_all , depth = 3 , prob = 0.90 , 
                            pars = c( "bGA" , "bGA2" ) )
Pre_model_GA_all

###3.1.2 Models of biological kin density ----
Model_GA_genetic_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_genetic = as.integer( Tie_genetic ) ,
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Age2 = as.numeric( sprintf( "%0.4f" , Age_std_2 ) ) , 
  Gender = as.integer( Gender ) , 
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

{set.seed(123)
  Model_GA_genetic <- ulam(
    alist(
      Tie_genetic ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + bGA2[Gender] * Age2 +  
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      bGA2[Gender] ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 )
    ) , data = Model_GA_genetic_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GA_genetic <- precis( Model_GA_genetic , depth = 3 , prob = 0.90 , 
                                pars = c( "bGA" , "bGA2" ) )
Pre_model_GA_genetic

###3.1.3 Models of affinal kin density ----
Model_GA_affine_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_affine = as.integer( Tie_affine ) ,
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Age2 = as.numeric( sprintf( "%0.4f" , Age_std_2 ) ) , 
  Gender = as.integer( Gender ) , 
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

{set.seed(123)
  Model_GA_affine <- ulam(
    alist(
      Tie_affine ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + bGA2[Gender] * Age2 +  
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      bGA2[Gender] ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 )
    ) , data = Model_GA_affine_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GA_affine <- precis( Model_GA_affine , depth = 3 , prob = 0.90 , 
                               pars = c( "bGA" , "bGA2" ) )
Pre_model_GA_affine

###3.1.4 Models of friend density ----
Model_GA_fri_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_fri = as.integer( Tie_fri ) ,
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Age2 = as.numeric( sprintf( "%0.4f" , Age_std_2 ) ) , 
  Gender = as.integer( Gender ) , 
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

{set.seed(123)
  Model_GA_fri <- ulam(
    alist(
      Tie_fri ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + bGA2[Gender] * Age2 +  
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      bGA2[Gender] ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 )
    ) , data = Model_GA_fri_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GA_fri <- precis( Model_GA_fri , depth = 3 , prob = 0.90 , 
                            pars = c( "bGA" , "bGA2" ) )
Pre_model_GA_fri

###3.1.5 Figures ----
# Posterior samples of estimates 
Post_model_GA_all <- extract.samples( Model_GA_all )
Post_model_GA_genetic <- extract.samples( Model_GA_genetic )
Post_model_GA_affine <- extract.samples( Model_GA_affine )
Post_model_GA_fri <- extract.samples( Model_GA_fri )

# Function to predict p-matrix based on posterior samples of estimates and new data
predict_p_matrix_GA <- function( post_samples , new_data , n_samples = 6000 ) {
  n_obs <- nrow( new_data ) 
  p_matrix <- matrix( NA , nrow = n_samples, ncol = n_obs ) # Null p_matrix
  
  # Loop over posterior samples
  for (s in 1:n_samples) {
    
    # Extract posterior sample for current iteration
    bGA <- post_samples$bGA[s, ]
    bGA2 <- post_samples$bGA2[s, ]
    V_bar <- post_samples$V_bar[s, 1 ]
    
    # Compute logit(p) for each observation in new_data
    logit_p <- with( new_data,
                     bGA[Gender] * Age + bGA2[Gender] * Age2 + 
                       V_bar )
    
    # Convert logit(p) to p value
    p_matrix[s, ] <- rethinking::inv_logit(logit_p)
  }
  
  # Return the predicted p-matrix
  return(p_matrix)
}

# New data for women and men, respectively
New_data_GA_F <- data.frame( Age_real = as.integer( 20:86 ) , 
                             Age = ( as.integer( 20:86 ) - mean( Ego.data$Age ) )/ sd( Ego.data$Age ) , 
                             Age2 = ( ( as.integer( 20:86 ) - mean( Ego.data$Age ) )/ sd( Ego.data$Age ) )^2 , 
                             Gender = as.integer( rep( 1 , 67 ) ) ) 

New_data_GA_M <- data.frame( Age_real = as.integer( 20:86 ) , 
                             Age = ( as.integer( 20:86 ) - mean( Ego.data$Age ) )/ sd( Ego.data$Age ) , 
                             Age2 = ( ( as.integer( 20:86 ) - mean( Ego.data$Age ) )/ sd( Ego.data$Age ) )^2 , 
                             Gender = as.integer( rep( 2 , 67 ) ) )

# Predicted p values based on new data
p_matrix_GA_all_F <- predict_p_matrix_GA( post_samples = Post_model_GA_all , 
                                          new_data = New_data_GA_F, 
                                          n_samples = 6000 )
p_matrix_GA_all_M <- predict_p_matrix_GA( post_samples = Post_model_GA_all , 
                                          new_data = New_data_GA_M, 
                                          n_samples = 6000 )

p_matrix_GA_genetic_F <- predict_p_matrix_GA( post_samples = Post_model_GA_genetic , 
                                              new_data = New_data_GA_F, 
                                              n_samples = 6000 )
p_matrix_GA_genetic_M <- predict_p_matrix_GA( post_samples = Post_model_GA_genetic , 
                                              new_data = New_data_GA_M, 
                                              n_samples = 6000 )

p_matrix_GA_affine_F <- predict_p_matrix_GA( post_samples = Post_model_GA_affine , 
                                             new_data = New_data_GA_F, 
                                             n_samples = 6000 )
p_matrix_GA_affine_M <- predict_p_matrix_GA( post_samples = Post_model_GA_affine , 
                                             new_data = New_data_GA_M, 
                                             n_samples = 6000 )

p_matrix_GA_fri_F <- predict_p_matrix_GA( post_samples = Post_model_GA_fri , 
                                          new_data = New_data_GA_F, 
                                          n_samples = 6000 )
p_matrix_GA_fri_M <- predict_p_matrix_GA( post_samples = Post_model_GA_fri , 
                                          new_data = New_data_GA_M, 
                                          n_samples = 6000 )

####3.1.5.1 Gender-specific effects of age ----
jpeg( "Age_stratified_by_gender.jpg" , 
      width = 200 , height = 160 , units = "mm" , res = 300 )

{par(mfrow=c(2,2), mar = c(4.5, 4.5, 1.8, 1))
  # All
  plot( NULL,
        xlim = c(20,86),
        ylim = c(0,1),
        xlab = "" ,
        ylab = "Predicted density" ,
        main = "All" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 , 
        cex.main = 1.5 )
  # women
  lines( c(20:86),
         apply( p_matrix_GA_all_F , 2 , mean) ,
         lwd = 2 , col = "#CB2313" )
  shade( apply( p_matrix_GA_all_F , 2 , PI , prob = 0.90 ) ,
         c(20:86) , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GA_all_F , 2 , PI , prob = 0.60 ) , 
         c(20:86) , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GA_all_F , 2 , PI , prob = 0.30 ) , 
         c(20:86) , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  
  # men
  lines( c(20:86),
         apply( p_matrix_GA_all_M , 2 , mean) ,
         lwd = 2 , col = "#046C9A" )
  shade( apply( p_matrix_GA_all_M , 2 , PI , prob = 0.90 ) , 
         c(20:86) , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GA_all_M , 2 , PI , prob = 0.60 ) , 
         c(20:86) , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GA_all_M , 2 , PI , prob = 0.30 ) , 
         c(20:86) , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  
  # Biological kin
  plot( NULL,
        xlim = c(20,86),
        ylim = c(0,1),
        xlab = "" ,
        ylab = "" ,
        main = "Biological kin" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 , 
        cex.main = 1.5 )
  legend( x = 65 , y = 1.05 ,  
          box.col = "white",
          legend = c( "Female" , "Male" ) , 
          lty = c( 1 , 1 ) , 
          col = c( "#CB2313" , "#046C9A" ) , 
          lwd = 2 ,
          cex = 1.2 , 
          bty = "n" ,
          y.intersp = 1.2 ,
          x.intersp = 0.3 ,
          seg.len = 0.8  )
  
  # women
  lines( c(20:86),
         apply( p_matrix_GA_genetic_F , 2 , mean) ,
         lwd = 2 , col = "#CB2313" )
  shade( apply( p_matrix_GA_genetic_F , 2 , PI , prob = 0.90 ) , 
         c(20:86) , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GA_genetic_F , 2 , PI , prob = 0.60 ) , 
         c(20:86) , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GA_genetic_F , 2 , PI , prob = 0.30 ) , 
         c(20:86) , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  # men
  lines( c(20:86),
         apply( p_matrix_GA_genetic_M , 2 , mean) ,
         lwd = 2 , col = "#046C9A" )
  shade( apply( p_matrix_GA_genetic_M , 2 , PI , prob = 0.90 ) , 
         c(20:86) , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GA_genetic_M , 2 , PI , prob = 0.60 ) , 
         c(20:86) , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GA_genetic_M , 2 , PI , prob = 0.30 ) , 
         c(20:86) , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  
  # Affinal kin
  plot( NULL,
        xlim = c(20,86),
        ylim = c(0,1),
        xlab = "Age" ,
        ylab = "Predicted density" ,
        main = "Affinal kin" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 , 
        cex.main = 1.5 )
  # women
  lines( c(20:86),
         apply( p_matrix_GA_affine_F , 2 , mean) ,
         lwd = 2 , col = "#CB2313" )
  shade( apply( p_matrix_GA_affine_F , 2 , PI , prob = 0.90 ) , 
         c(20:86) , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GA_affine_F , 2 , PI , prob = 0.60 ) , 
         c(20:86) , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GA_affine_F , 2 , PI , prob = 0.30 ) , 
         c(20:86) , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  # men
  lines( c(20:86),
         apply( p_matrix_GA_affine_M , 2 , mean) ,
         lwd = 2 , col = "#046C9A" )
  shade( apply( p_matrix_GA_affine_M , 2 , PI , prob = 0.90 ) , 
         c(20:86) , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GA_affine_M , 2 , PI , prob = 0.60 ) , 
         c(20:86) , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GA_affine_M , 2 , PI , prob = 0.30 ) , 
         c(20:86) , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  
  # Friend kin
  plot( NULL,
        xlim = c(20,86),
        ylim = c(0,1),
        xlab = "Age" ,
        ylab = "" ,
        main = "Friend" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 , 
        cex.main = 1.5 )
  # women
  lines( c(20:86),
         apply( p_matrix_GA_fri_F , 2 , mean) ,
         lwd = 2 , col = "#CB2313" )
  shade( apply( p_matrix_GA_fri_F , 2 , PI , prob = 0.90 ) , 
         c(20:86) , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GA_fri_F , 2 , PI , prob = 0.60 ) , 
         c(20:86) , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GA_fri_F , 2 , PI , prob = 0.30 ) , 
         c(20:86) , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  # men
  lines( c(20:86),
         apply( p_matrix_GA_fri_M , 2 , mean) ,
         lwd = 2 , col = "#046C9A" )
  shade( apply( p_matrix_GA_fri_M , 2 , PI , prob = 0.90 ) , 
         c(20:86) , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GA_fri_M , 2 , PI , prob = 0.60 ) , 
         c(20:86) , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GA_fri_M , 2 , PI , prob = 0.30 ) , 
         c(20:86) , 
         col = col.alpha( "#046C9A" , 0.15 ) )
}

dev.off()

####3.1.5.2 Gender differences across age ----
jpeg( "Age_gender_diff.jpg" , 
      width = 200 , height = 160 , units = "mm" , res = 300 )

{par(mfrow=c(2,2), mar = c(4.5, 4.5, 1.8, 1))
  # All
  plot( NULL,
        xlim = c(20,86),
        ylim = c(-0.5,0.5),
        xlab = "" ,
        ylab = "Predicted density" ,
        main = "All" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 , 
        cex.main = 1.5 )
  lines( c(20:86),
         apply( p_matrix_GA_all_F - p_matrix_GA_all_M , 2 , mean ) ,
         lwd = 2 , col = "black" )
  shade( apply( p_matrix_GA_all_F - p_matrix_GA_all_M , 2 , PI , prob = 0.90 ) ,
         c(20:86) , 
         col = col.alpha( "dimgrey" , 0.15 ) )
  shade( apply( p_matrix_GA_all_F - p_matrix_GA_all_M , 2 , PI , prob = 0.60 ) ,
         c(20:86) , 
         col = col.alpha( "dimgrey" , 0.15 ) )
  shade( apply( p_matrix_GA_all_F - p_matrix_GA_all_M , 2 , PI , prob = 0.30 ) ,
         c(20:86) , 
         col = col.alpha( "dimgrey" , 0.15 ) )
  
  # Biological kin
  plot( NULL,
        xlim = c(20,86),
        ylim = c(-0.5,0.5),
        xlab = "" ,
        ylab = "" ,
        main = "Biological kin" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 , 
        cex.main = 1.5 )
  lines( c(20:86),
         apply( p_matrix_GA_genetic_F - p_matrix_GA_genetic_M , 2 , mean ) ,
         lwd = 2 , col = "black" )
  shade( apply( p_matrix_GA_genetic_F - p_matrix_GA_genetic_M , 2 , PI , prob = 0.90 ) ,
         c(20:86) , 
         col = col.alpha( "dimgrey" , 0.15 ) )
  shade( apply( p_matrix_GA_genetic_F - p_matrix_GA_genetic_M , 2 , PI , prob = 0.60 ) ,
         c(20:86) , 
         col = col.alpha( "dimgrey" , 0.15 ) )
  shade( apply( p_matrix_GA_genetic_F - p_matrix_GA_genetic_M , 2 , PI , prob = 0.30 ) ,
         c(20:86) , 
         col = col.alpha( "dimgrey" , 0.15 ) )
  
  # Affinal kin
  plot( NULL,
        xlim = c(20,86),
        ylim = c(-0.5,0.5),
        xlab = "Age" ,
        ylab = "Predicted density" ,
        main = "Affinal kin" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 , 
        cex.main = 1.5 )
  lines( c(20:86),
         apply( p_matrix_GA_affine_F - p_matrix_GA_affine_M , 2 , mean ) ,
         lwd = 2 , col = "black" )
  shade( apply( p_matrix_GA_affine_F - p_matrix_GA_affine_M , 2 , PI , prob = 0.90 ) , 
         c(20:86) , 
         col = col.alpha( "dimgrey" , 0.15 ) )
  shade( apply( p_matrix_GA_affine_F - p_matrix_GA_affine_M , 2 , PI , prob = 0.60 ) ,
         c(20:86) , 
         col = col.alpha( "dimgrey" , 0.15 ) )
  shade( apply( p_matrix_GA_affine_F - p_matrix_GA_affine_M , 2 , PI , prob = 0.30 ) ,
         c(20:86) , 
         col = col.alpha( "dimgrey" , 0.15 ) )
  
  # Friend kin
  plot( NULL,
        xlim = c(20,86),
        ylim = c(-0.5,0.5),
        xlab = "Age" ,
        ylab = "" ,
        main = "Friend" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 , 
        cex.main = 1.5 )
  lines( c(20:86),
         apply( p_matrix_GA_fri_F - p_matrix_GA_fri_M , 2 , mean ) ,
         lwd = 2 , col = "black" )
  shade( apply( p_matrix_GA_fri_F - p_matrix_GA_fri_M , 2 , PI , prob = 0.90  ) ,
         c(20:86) , 
         col = col.alpha( "dimgrey" , 0.15 ) )
  shade( apply( p_matrix_GA_fri_F - p_matrix_GA_fri_M , 2 , PI , prob = 0.60 ) ,
         c(20:86) , 
         col = col.alpha( "dimgrey" , 0.15 ) )
  shade( apply( p_matrix_GA_fri_F - p_matrix_GA_fri_M , 2 , PI , prob = 0.30 ) ,
         c(20:86) , 
         col = col.alpha( "dimgrey" , 0.15 ) )
}

dev.off()

###3.1.6 Outputs ----
Output_model_GA <- data.frame( Mean = c( Pre_model_GA_all$mean ,
                                         Pre_model_GA_genetic$mean ,
                                         Pre_model_GA_affine$mean ,
                                         Pre_model_GA_fri$mean ) ,
                               CI5 = c( Pre_model_GA_all$`5%` ,
                                        Pre_model_GA_genetic$`5%` ,
                                        Pre_model_GA_affine$`5%` ,
                                        Pre_model_GA_fri$`5%` ) ,
                               CI95 = c( Pre_model_GA_all$`95%` ,
                                         Pre_model_GA_genetic$`95%` ,
                                         Pre_model_GA_affine$`95%` ,
                                         Pre_model_GA_fri$`95%` ) ,
                               Type = rep( c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) , each = 4 ) ,
                               Variable = rep( c( "Age[female]" , "Age[male]" , 
                                                  "Age^2[female]" , "Age^2[male]" ) , 4 ) ) %>% 
  mutate( Type = factor( Type , 
                         levels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ,
                         labels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) ,
          Variable = factor( Variable , 
                             levels = c( "Age[female]" , "Age^2[female]" , 
                                         "Age[male]" , "Age^2[male]" ) ,
                             labels = c( "Age[female]" , "Age^2[female]" , 
                                         "Age[male]" , "Age^2[male]" ) ) ,
          Model = "Age" , 
          "Mean [CI]" = paste( format( round( Mean , 2 ), nsmall = 2 ) , 
                               " [" , 
                               format( round( CI5 , 2 ), nsmall = 2 ) , 
                               ", " , 
                               format( round( CI95 , 2 ), nsmall = 2 ) , 
                               "]", 
                               sep = "") )

Age.gender.output <- tibble( Variable = c( "Age[female]" , "Age[male]" , 
                                           "Age^2[female]" , "Age^2[male]" ) ,
                             Type = "All" ) %>% 
  left_join( Output_model_GA[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `All` = "Mean [CI]" ) %>% 
  mutate( Type = "Biological kin" ) %>% 
  left_join( Output_model_GA[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `Biological kin` = "Mean [CI]" ) %>% 
  mutate( Type = "Affinal kin" ) %>% 
  left_join( Output_model_GA[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `Affinal kin` = "Mean [CI]" ) %>% 
  mutate( Type = "Friend" ) %>% 
  left_join( Output_model_GA[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `Friend` = "Mean [CI]" ) %>% 
  select( - Type ) %>% 
  gt() %>% 
  fmt_number( ) %>% 
  sub_missing( columns = everything() , 
               missing_text = " " ) %>% 
  tab_style( style = list( cell_text( weight = "bold" ) ) ,
             locations = cells_column_labels( everything() ) ) %>% 
  tab_style( style = list( cell_text( font = "Times New Roman" ) ) ,
             locations = cells_body() ) %>% 
  cols_align( align = "left",
              columns = "Variable" ) %>%
  cols_align( align = "center",
              columns = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) %>% 
  opt_table_lines( extent = "default" ) %>%
  tab_options( column_labels.border.top.color = "black" ,
               column_labels.border.top.width = px( 3 ) ,
               column_labels.border.bottom.color = "black" ,
               table_body.hlines.color = "lightgrey" ,
               table.border.bottom.color = "black" ,
               table.border.bottom.width = px( 3 ) ) %>% 
  cols_width( "Variable" ~ px( 170 ) , 
              everything() ~ px( 150 ) )
Age.gender.output

##3.2 Residence ----
# Other variable should be included into models
adjustmentSets( Dag , 
                exposure = "Residence" , 
                outcome = "Density" ,
                effect = "total" )

###3.2.1 Models of all density ----
Model_GR_all_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_all = as.integer( Tie_all ) , 
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Age2 = as.numeric( sprintf( "%0.4f" , Age_std_2 ) ) , 
  Gender = as.integer( Gender ) , 
  PMR = as.integer( Residence ) , 
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

{set.seed(123)
  Model_GR_all <- ulam(
    alist(
      Tie_all ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + bGA2[Gender] * Age2 + 
        bGR[Gender,PMR] + 
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      bGA2[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 ) 
    ) , data = Model_GR_all_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}

Pre_model_GR_all <- precis( Model_GR_all , depth = 3 , prob = 0.90 , 
                            pars = c( "bGA" , "bGA2" , "bGR" ) )
Pre_model_GR_all

###3.2.2 Models of biological kin density ----
Model_GR_genetic_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_genetic = as.integer( Tie_genetic ) ,
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Age2 = as.numeric( sprintf( "%0.4f" , Age_std_2 ) ) , 
  Gender = as.integer( Gender ) , 
  PMR = as.integer( Residence ) , 
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

{set.seed(123)
  Model_GR_genetic <- ulam(
    alist(
      Tie_genetic ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + bGA2[Gender] * Age2 + 
        bGR[Gender,PMR] + 
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      bGA2[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 ) 
    ) , data = Model_GR_genetic_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GR_genetic <- precis( Model_GR_genetic , depth = 3 , prob = 0.90 , 
                                pars = c( "bGA" , "bGA2" , "bGR" ) )
Pre_model_GR_genetic

###3.2.3 Models of affinal kin density ----
Model_GR_affine_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_affine = as.integer( Tie_affine ) ,
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Age2 = as.numeric( sprintf( "%0.4f" , Age_std_2 ) ) , 
  Gender = as.integer( Gender ) , 
  PMR = as.integer( Residence ) , 
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

{set.seed(123)
  Model_GR_affine <- ulam(
    alist(
      Tie_affine ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + bGA2[Gender] * Age2 + 
        bGR[Gender,PMR] + 
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      bGA2[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 ) 
    ) , data = Model_GR_affine_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GR_affine <- precis( Model_GR_affine , depth = 3 , prob = 0.90 , 
                               pars = c( "bGA" , "bGA2" , "bGR" ) )
Pre_model_GR_affine

###3.2.4 Models of friend density ----
Model_GR_fri_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_fri = as.integer( Tie_fri ) ,
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Age2 = as.numeric( sprintf( "%0.4f" , Age_std_2 ) ) , 
  Gender = as.integer( Gender ) , 
  PMR = as.integer( Residence ) , 
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

{set.seed(123)
  Model_GR_fri <- ulam(
    alist(
      Tie_fri ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + bGA2[Gender] * Age2 + 
        bGR[Gender,PMR] + 
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      bGA2[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 ) 
    ) , data = Model_GR_fri_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GR_fri <- precis( Model_GR_fri , depth = 3 , prob = 0.90 , 
                            pars = c( "bGA" , "bGA2" , "bGR" ) )
Pre_model_GR_fri

###3.2.5 Figures ----
# Posterior samples of estimates 
Post_model_GR_all <- extract.samples( Model_GR_all )
Post_model_GR_genetic <- extract.samples( Model_GR_genetic )
Post_model_GR_affine <- extract.samples( Model_GR_affine )
Post_model_GR_fri <- extract.samples( Model_GR_fri )

####3.2.5.1 Gender-specific effects of residence ----
#####3.2.5.1.1 Effects of residence for women and men, respectively ----
# Women
tibble( Value = c( Post_model_GR_all$bGR[,1,1] , 
                   Post_model_GR_all$bGR[,1,2] , 
                   Post_model_GR_all$bGR[,1,3] ,
                   
                   Post_model_GR_genetic$bGR[,1,1] ,
                   Post_model_GR_genetic$bGR[,1,2] ,
                   Post_model_GR_genetic$bGR[,1,3] ,
                   
                   Post_model_GR_affine$bGR[,1,1] ,
                   Post_model_GR_affine$bGR[,1,2] ,
                   Post_model_GR_affine$bGR[,1,3] ,
                   
                   Post_model_GR_fri$bGR[,1,1] ,
                   Post_model_GR_fri$bGR[,1,2] ,
                   Post_model_GR_fri$bGR[,1,3] ) ,
        Gender = c( rep( "Female" , 72000 ) ) ,
        Residence = c(rep( c( rep( "Matrilocal" , 6000 ) , 
                              rep( "Bilocal" , 6000 ) , 
                              rep( "Patrilocal" , 6000 ) ) , 4 ) ) ,
        Type = c( rep( "All" , 18000 ) ,
                  rep( "Biological kin" , 18000 ) ,
                  rep( "Affinal kin" , 18000 ) ,
                  rep( "Friend" , 18000 ) ) ) %>% 
  mutate( Residence = factor( Residence , 
                              levels = c( "Matrilocal" , "Bilocal" , "Patrilocal" ) , 
                              labels = c( "Matrilocal" , "Bilocal" , "Patrilocal" ) ) , 
          Type = factor( Type , 
                         levels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) , 
                         labels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) ) %>% 
  ggplot( aes( x = Value , y = Residence , 
               fill = Residence , color = Residence) ) +
  facet_wrap( ~ Type ) +
  geom_vline( xintercept = 0, linetype = 2 , linewidth = 1 ,
              color = "dimgray" ) +
  stat_halfeye( .width = .90, height = 0.8 ) + 
  labs( title = "Female" ) +
  scale_fill_manual( values = alpha( c( "#c5272d" , "#037f77" , "#0001a1" ) , 0.4 ) ) +
  scale_color_manual( values = c( "#c5272d" , "#037f77" , "#0001a1" ) ) +
  scale_y_discrete( NULL, labels = ggplot2:::parse_safe ) +
  xlab( "Posterior estimates" ) +
  coord_cartesian( ylim = c( 1.5, 3.6 ) ) +
  scale_x_continuous( limits = c( -2 , 2 ) , breaks = c( -1 , 0 , 1 ) ) +
  theme(plot.margin = margin(10, 10, 10, 20),
        strip.background = element_rect(color = "black", fill = "white") ,
        strip.text.x = element_text(size = 16),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "Gainsboro"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA),
        panel.background = element_blank(),
        plot.title = element_text(size=16),
        legend.position = "none",
        axis.title.x = element_text(size = 16,
                                    margin = margin(t = 0.3, r = 0, b = 0, l = 0,unit = "cm")),
        axis.text.x = element_text(colour="black",size=14,
                                   margin = margin(t = 0.1, r = 0, b = 0, l = 0,unit = "cm")),
        axis.title.y = element_text(size = 16,
                                    margin = margin(t = 0, r = 0.5, b = 0, l = 0,unit = "cm")),
        axis.text.y = element_text(colour="black",size=14,
                                   margin = margin(t = 0, r = 0.2, b = 0, l = 0,unit = "cm")),
        legend.title=element_text(size=14,face = "bold"),
        legend.text=element_text(size=12,face = "bold"))

ggsave( filename = "Residence_women.jpeg" , 
        width = 160 , height = 150 , units = "mm" , dpi = 300 )

# Men
tibble( Value = c( Post_model_GR_all$bGR[,2,1] , 
                   Post_model_GR_all$bGR[,2,2] , 
                   Post_model_GR_all$bGR[,2,3] ,
                   
                   Post_model_GR_genetic$bGR[,2,1] ,
                   Post_model_GR_genetic$bGR[,2,2] ,
                   Post_model_GR_genetic$bGR[,2,3] ,
                   
                   Post_model_GR_affine$bGR[,2,1] ,
                   Post_model_GR_affine$bGR[,2,2] ,
                   Post_model_GR_affine$bGR[,2,3] ,
                   
                   Post_model_GR_fri$bGR[,2,1] ,
                   Post_model_GR_fri$bGR[,2,2] ,
                   Post_model_GR_fri$bGR[,2,3] ) ,
        Gender = c( rep( "Male" , 72000 ) ) ,
        Residence = c(rep( c( rep( "Matrilocal" , 6000 ) , 
                              rep( "Bilocal" , 6000 ) , 
                              rep( "Patrilocal" , 6000 ) ) , 4 ) ) ,
        Type = c( rep( "All" , 18000 ) ,
                  rep( "Biological kin" , 18000 ) ,
                  rep( "Affinal kin" , 18000 ) ,
                  rep( "Friend" , 18000 ) ) ) %>% 
  mutate( Residence = factor( Residence , 
                              levels = c( "Matrilocal" , "Bilocal" , "Patrilocal" ) , 
                              labels = c( "Matrilocal" , "Bilocal" , "Patrilocal" ) ) , 
          Type = factor( Type , 
                         levels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) , 
                         labels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) ) %>% 
  ggplot( aes( x = Value , y = Residence , 
               fill = Residence , color = Residence) ) +
  facet_wrap( ~ Type ) +
  geom_vline( xintercept = 0, linetype = 2 , linewidth = 1 ,
              color = "dimgray" ) +
  stat_halfeye( .width = .90, height = 0.8 ) + 
  labs( title = "Male" ) +
  scale_fill_manual( values = alpha( c( "#c5272d" , "#037f77" , "#0001a1" ) , 0.4 ) ) +
  scale_color_manual( values = c( "#c5272d" , "#037f77" , "#0001a1" ) ) +
  scale_y_discrete( NULL, labels = ggplot2:::parse_safe ) +
  xlab( "Posterior estimates" ) +
  coord_cartesian( ylim = c( 1.5, 3.6 ) ) +
  scale_x_continuous( limits = c( -2 , 2 ) , breaks = c( -1 , 0 , 1 ) ) +
  theme(plot.margin = margin(10, 10, 10, 20),
        strip.background = element_rect(color = "black", fill = "white") ,
        strip.text.x = element_text(size = 16),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "Gainsboro"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA),
        panel.background = element_blank(),
        plot.title = element_text(size=16),
        legend.position = "none",
        axis.title.x = element_text(size = 16,
                                    margin = margin(t = 0.3, r = 0, b = 0, l = 0,unit = "cm")),
        axis.text.x = element_text(colour="black",size=14,
                                   margin = margin(t = 0.1, r = 0, b = 0, l = 0,unit = "cm")),
        axis.title.y = element_text(size = 16,
                                    margin = margin(t = 0, r = 0.5, b = 0, l = 0,unit = "cm")),
        axis.text.y = element_text(colour="black",size=14,
                                   margin = margin(t = 0, r = 0.2, b = 0, l = 0,unit = "cm")),
        legend.title=element_text(size=14,face = "bold"),
        legend.text=element_text(size=12,face = "bold"))

ggsave( filename = "Residence_men.jpeg" , 
        width = 160 , height = 150 , units = "mm" , dpi = 300 )

#####3.2.5.1.2 Residence difference for women and men, respectively ----
# Women
tibble( Value = c( c( Post_model_GR_all$bGR[,1,1] - Post_model_GR_all$bGR[,1,2] ) , 
                   c( Post_model_GR_all$bGR[,1,1] - Post_model_GR_all$bGR[,1,3] ) , 
                   c( Post_model_GR_all$bGR[,1,2] - Post_model_GR_all$bGR[,1,3] ) , 
                   
                   c( Post_model_GR_genetic$bGR[,1,1] - Post_model_GR_genetic$bGR[,1,2] ) , 
                   c( Post_model_GR_genetic$bGR[,1,1] - Post_model_GR_genetic$bGR[,1,3] ) , 
                   c( Post_model_GR_genetic$bGR[,1,2] - Post_model_GR_genetic$bGR[,1,3] ) , 
                   
                   c( Post_model_GR_affine$bGR[,1,1] - Post_model_GR_affine$bGR[,1,2] ) , 
                   c( Post_model_GR_affine$bGR[,1,1] - Post_model_GR_affine$bGR[,1,3] ) , 
                   c( Post_model_GR_affine$bGR[,1,2] - Post_model_GR_affine$bGR[,1,3] ) , 
                   
                   c( Post_model_GR_fri$bGR[,1,1] - Post_model_GR_fri$bGR[,1,2] ) , 
                   c( Post_model_GR_fri$bGR[,1,1] - Post_model_GR_fri$bGR[,1,3] ) , 
                   c( Post_model_GR_fri$bGR[,1,2] - Post_model_GR_fri$bGR[,1,3] ) ) ,
        Gender = c( rep( "Female" , 72000 ) ) ,
        Residence = c(rep( c( rep( "Matrilocal - Bilocal" , 6000 ) , 
                              rep( "Matrilocal - Patrilocal" , 6000 ) , 
                              rep( "Bilocal - Patrilocal" , 6000 ) ) , 4 ) ) ,
        Type = c( rep( "All" , 18000 ) ,
                  rep( "Biological kin" , 18000 ) ,
                  rep( "Affinal kin" , 18000 ) ,
                  rep( "Friend" , 18000 ) ) ) %>% 
  mutate( Residence = factor( Residence , 
                              levels = c( "Bilocal - Patrilocal" , 
                                          "Matrilocal - Bilocal" , 
                                          "Matrilocal - Patrilocal" ) , 
                              labels = c( "Bilocal - Patrilocal" , 
                                          "Matrilocal - Bilocal" , 
                                          "Matrilocal - Patrilocal" ) ) , 
          Type = factor( Type , 
                         levels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) , 
                         labels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) ) %>% 
  ggplot( aes( x = Value , y = Residence , 
               fill = Residence , color = Residence) ) +
  facet_wrap( ~ Type ) +
  geom_vline( xintercept = 0, linetype = 2 , linewidth = 1 ,
              color = "dimgray" ) +
  stat_halfeye( .width = .90, height = 0.8 ) + 
  labs( title = "Female" ) +
  scale_fill_manual( values = alpha( c( "#f1b656" , "#397fc7" , "#040676" ) , 0.4 ) ) +
  scale_color_manual( values = c( "#f1b656" , "#397fc7" , "#040676" ) ) +
  scale_y_discrete( NULL, labels = ggplot2:::parse_safe ) +
  xlab( "Posterior estimates" ) +
  coord_cartesian( ylim = c( 1.5, 3.2 ) ) +
  scale_x_continuous( limits = c( -1 , 1 ) , breaks = c( -1 , 0 , 1 ) ) +
  theme(plot.margin = margin(10, 10, 10, 20),
        strip.background = element_rect(color = "black", fill = "white") ,
        strip.text.x = element_text(size = 16),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "Gainsboro"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA),
        panel.background = element_blank(),
        plot.title = element_text(size=16),
        legend.position = "none",
        axis.title.x = element_text(size = 16,
                                    margin = margin(t = 0.3, r = 0, b = 0, l = 0,unit = "cm")),
        axis.text.x = element_text(colour="black",size=14,
                                   margin = margin(t = 0.1, r = 0, b = 0, l = 0,unit = "cm")),
        axis.title.y = element_text(size = 16,
                                    margin = margin(t = 0, r = 0.5, b = 0, l = 0,unit = "cm")),
        axis.text.y = element_text(colour="black",size=14,
                                   margin = margin(t = 0, r = 0.2, b = 0, l = 0,unit = "cm")),
        legend.title=element_text(size=14,face = "bold"),
        legend.text=element_text(size=12,face = "bold"))

ggsave( filename = "Residence_diff_women.jpeg" , 
        width = 200 , height = 150 , units = "mm" , dpi = 300 )

# Men
tibble( Value = c( c( Post_model_GR_all$bGR[,2,1] - Post_model_GR_all$bGR[,2,2] ) , 
                   c( Post_model_GR_all$bGR[,2,1] - Post_model_GR_all$bGR[,2,3] ) , 
                   c( Post_model_GR_all$bGR[,2,2] - Post_model_GR_all$bGR[,2,3] ) , 
                   
                   c( Post_model_GR_genetic$bGR[,2,1] - Post_model_GR_genetic$bGR[,2,2] ) , 
                   c( Post_model_GR_genetic$bGR[,2,1] - Post_model_GR_genetic$bGR[,2,3] ) , 
                   c( Post_model_GR_genetic$bGR[,2,2] - Post_model_GR_genetic$bGR[,2,3] ) , 
                   
                   c( Post_model_GR_affine$bGR[,2,1] - Post_model_GR_affine$bGR[,2,2] ) , 
                   c( Post_model_GR_affine$bGR[,2,1] - Post_model_GR_affine$bGR[,2,3] ) , 
                   c( Post_model_GR_affine$bGR[,2,2] - Post_model_GR_affine$bGR[,2,3] ) , 
                   
                   c( Post_model_GR_fri$bGR[,2,1] - Post_model_GR_fri$bGR[,2,2] ) , 
                   c( Post_model_GR_fri$bGR[,2,1] - Post_model_GR_fri$bGR[,2,3] ) , 
                   c( Post_model_GR_fri$bGR[,2,2] - Post_model_GR_fri$bGR[,2,3] ) ) ,
        Gender = c( rep( "Male" , 72000 ) ) ,
        Residence = c(rep( c( rep( "Matrilocal - Bilocal" , 6000 ) , 
                              rep( "Matrilocal - Patrilocal" , 6000 ) , 
                              rep( "Bilocal - Patrilocal" , 6000 ) ) , 4 ) ) ,
        Type = c( rep( "All" , 18000 ) ,
                  rep( "Biological kin" , 18000 ) ,
                  rep( "Affinal kin" , 18000 ) ,
                  rep( "Friend" , 18000 ) ) ) %>% 
  mutate( Residence = factor( Residence , 
                              levels = c( "Bilocal - Patrilocal" , 
                                          "Matrilocal - Bilocal" , 
                                          "Matrilocal - Patrilocal" ) , 
                              labels = c( "Bilocal - Patrilocal" , 
                                          "Matrilocal - Bilocal" , 
                                          "Matrilocal - Patrilocal" ) ) , 
          Type = factor( Type , 
                         levels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) , 
                         labels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) ) %>% 
  ggplot( aes( x = Value , y = Residence , 
               fill = Residence , color = Residence) ) +
  facet_wrap( ~ Type ) +
  geom_vline( xintercept = 0, linetype = 2 , linewidth = 1 ,
              color = "dimgray" ) +
  stat_halfeye( .width = .90, height = 0.8 ) + 
  labs( title = "Male" ) +
  scale_fill_manual( values = alpha( c( "#f1b656" , "#397fc7" , "#040676" ) , 0.4 ) ) +
  scale_color_manual( values = c( "#f1b656" , "#397fc7" , "#040676" ) ) +
  scale_y_discrete( NULL, labels = ggplot2:::parse_safe ) +
  xlab( "Posterior estimates" ) +
  coord_cartesian( ylim = c( 1.5, 3.2 ) ) +
  scale_x_continuous( limits = c( -1 , 1 ) , breaks = c( -1 , 0 , 1 ) ) +
  theme(plot.margin = margin(10, 10, 10, 20),
        strip.background = element_rect(color = "black", fill = "white") ,
        strip.text.x = element_text(size = 16),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "Gainsboro"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA),
        panel.background = element_blank(),
        plot.title = element_text(size=16),
        legend.position = "none",
        axis.title.x = element_text(size = 16,
                                    margin = margin(t = 0.3, r = 0, b = 0, l = 0,unit = "cm")),
        axis.text.x = element_text(colour="black",size=14,
                                   margin = margin(t = 0.1, r = 0, b = 0, l = 0,unit = "cm")),
        axis.title.y = element_text(size = 16,
                                    margin = margin(t = 0, r = 0.5, b = 0, l = 0,unit = "cm")),
        axis.text.y = element_text(colour="black",size=14,
                                   margin = margin(t = 0, r = 0.2, b = 0, l = 0,unit = "cm")),
        legend.title=element_text(size=14,face = "bold"),
        legend.text=element_text(size=12,face = "bold"))

ggsave( filename = "Residence_diff_men.jpeg" , 
        width = 200 , height = 150 , units = "mm" , dpi = 300 )

# Proportions of mean difference in βR values comparing each pairwise of residences 
tibble( Type = c( rep( "All" , 3 ) ,
                  rep( "Biological kin" , 3 ) ,
                  rep( "Affinal kin" , 3 ) ,
                  rep( "Friend" , 3 ) ) ,
        Residence = c(rep( c( "Matrilocal - Patrilocal" , 
                              "Matrilocal - Bilocal" , 
                              "Bilocal - Patrilocal" ) , 4 ) ) ,
        Female = c( sprintf("%0.4f" , 
                            length( which( c( Post_model_GR_all$bGR[,1,1] - Post_model_GR_all$bGR[,1,3] ) > 0 ) ) /6000 ) , 
                    sprintf("%0.4f" , 
                            length( which( c( Post_model_GR_all$bGR[,1,1] - Post_model_GR_all$bGR[,1,2] ) > 0 ) ) /6000 ) , 
                    sprintf("%0.4f" , 
                            length( which( c( Post_model_GR_all$bGR[,1,2] - Post_model_GR_all$bGR[,1,3] ) > 0 ) ) /6000 ) , 
                    
                    sprintf("%0.4f" , 
                            length( which( c( Post_model_GR_genetic$bGR[,1,1] - Post_model_GR_genetic$bGR[,1,3] ) > 0 ) ) /6000 ) , 
                    sprintf("%0.4f" , 
                            length( which( c( Post_model_GR_genetic$bGR[,1,1] - Post_model_GR_genetic$bGR[,1,2] ) > 0 ) ) /6000 ) , 
                    sprintf("%0.4f" , 
                            length( which( c( Post_model_GR_genetic$bGR[,1,2] - Post_model_GR_genetic$bGR[,1,3] ) > 0 ) ) /6000 ) , 
                    
                    sprintf("%0.4f" , 
                            length( which( c( Post_model_GR_affine$bGR[,1,1] - Post_model_GR_affine$bGR[,1,3] ) > 0 ) ) /6000 ) , 
                    sprintf("%0.4f" , 
                            length( which( c( Post_model_GR_affine$bGR[,1,1] - Post_model_GR_affine$bGR[,1,2] ) > 0 ) ) /6000 ) , 
                    sprintf("%0.4f" , 
                            length( which( c( Post_model_GR_affine$bGR[,1,2] - Post_model_GR_affine$bGR[,1,3] ) > 0 ) ) /6000 ) , 
                    
                    sprintf("%0.4f" , 
                            length( which( c( Post_model_GR_fri$bGR[,1,1] - Post_model_GR_fri$bGR[,1,3] ) > 0 ) ) /6000 ) , 
                    sprintf("%0.4f" , 
                            length( which( c( Post_model_GR_fri$bGR[,1,1] - Post_model_GR_fri$bGR[,1,2] ) > 0 ) ) /6000 ) , 
                    sprintf("%0.4f" , 
                            length( which( c( Post_model_GR_fri$bGR[,1,2] - Post_model_GR_fri$bGR[,1,3] ) > 0 ) ) /6000 ) ) ) %>% 
  bind_cols (
    tibble( Male = c( sprintf("%0.4f" , 
                              length( which( c( Post_model_GR_all$bGR[,2,1] - Post_model_GR_all$bGR[,2,3] ) > 0 ) ) /6000 ) , 
                      sprintf("%0.4f" , 
                              length( which( c( Post_model_GR_all$bGR[,2,1] - Post_model_GR_all$bGR[,2,2] ) > 0 ) ) /6000 ) , 
                      sprintf("%0.4f" , 
                              length( which( c( Post_model_GR_all$bGR[,2,2] - Post_model_GR_all$bGR[,2,3] ) > 0 ) ) /6000 ) , 
                      
                      sprintf("%0.4f" , 
                              length( which( c( Post_model_GR_genetic$bGR[,2,1] - Post_model_GR_genetic$bGR[,2,3] ) > 0 ) ) /6000 ) , 
                      sprintf("%0.4f" , 
                              length( which( c( Post_model_GR_genetic$bGR[,2,1] - Post_model_GR_genetic$bGR[,2,2] ) > 0 ) ) /6000 ) , 
                      sprintf("%0.4f" , 
                              length( which( c( Post_model_GR_genetic$bGR[,2,2] - Post_model_GR_genetic$bGR[,2,3] ) > 0 ) ) /6000 ) , 
                      
                      sprintf("%0.4f" , 
                              length( which( c( Post_model_GR_affine$bGR[,2,1] - Post_model_GR_affine$bGR[,2,3] ) > 0 ) ) /6000 ) , 
                      sprintf("%0.4f" , 
                              length( which( c( Post_model_GR_affine$bGR[,2,1] - Post_model_GR_affine$bGR[,2,2] ) > 0 ) ) /6000 ) , 
                      sprintf("%0.4f" , 
                              length( which( c( Post_model_GR_affine$bGR[,2,2] - Post_model_GR_affine$bGR[,2,3] ) > 0 ) ) /6000 ) , 
                      
                      sprintf("%0.4f" , 
                              length( which( c( Post_model_GR_fri$bGR[,2,1] - Post_model_GR_fri$bGR[,2,3] ) > 0 ) ) /6000 ) , 
                      sprintf("%0.4f" , 
                              length( which( c( Post_model_GR_fri$bGR[,2,1] - Post_model_GR_fri$bGR[,2,2] ) > 0 ) ) /6000 ) , 
                      sprintf("%0.4f" , 
                              length( which( c( Post_model_GR_fri$bGR[,2,2] - Post_model_GR_fri$bGR[,2,3] ) > 0 ) ) /6000 ) ) )
  ) %>% 
  mutate( 
    Female = as.numeric( Female ) ,
    Male = as.numeric( Male ) ,
    Female = paste0( format( Female * 100 , nsmall = 2) , "%" ) ,
    Male = paste0( format( Male * 100 , nsmall = 2 ) , "%" ) ) %>% 
  gt() %>% 
  fmt_number( ) %>% 
  sub_missing( columns = everything() , 
               missing_text = " " ) %>% 
  tab_style( style = list( cell_text( weight = "bold" ) ) ,
             locations = cells_column_labels(everything() ) ) %>% 
  tab_style( style = list( cell_text(font = "Times New Roman") ) ,
             locations = cells_body() ) %>% 
  cols_align( align = "center",
              columns = c( "Type" , "Residence" , "Female" , "Male" ) ) %>% 
  opt_table_lines( extent = "default" ) %>%
  tab_options( column_labels.border.top.color = "black" ,
               column_labels.border.top.width = px( 3 ) ,
               column_labels.border.bottom.color = "black" ,
               table_body.hlines.color = "lightgrey" ,
               table.border.bottom.color = "black" ,
               table.border.bottom.width = px( 3 ) ) %>% 
  cols_width( "Residence" ~ px( 170 ) , 
              "Female" ~ px( 100 ) , 
              "Male" ~ px( 100 ) , 
              everything() ~ px( 150 ) )

####3.2.5.2 Gender difference across residence ----
# F - M
tibble( Value = c( c( Post_model_GR_all$bGR[,1,1] - Post_model_GR_all$bGR[,2,1] ) , 
                   c( Post_model_GR_all$bGR[,1,2] - Post_model_GR_all$bGR[,2,2] ) , 
                   c( Post_model_GR_all$bGR[,1,3] - Post_model_GR_all$bGR[,2,3] ) , 
                   
                   c( Post_model_GR_genetic$bGR[,1,1] - Post_model_GR_genetic$bGR[,2,1] ) , 
                   c( Post_model_GR_genetic$bGR[,1,2] - Post_model_GR_genetic$bGR[,2,2] ) , 
                   c( Post_model_GR_genetic$bGR[,1,3] - Post_model_GR_genetic$bGR[,2,3] ) , 
                   
                   c( Post_model_GR_affine$bGR[,1,1] - Post_model_GR_affine$bGR[,2,1] ) , 
                   c( Post_model_GR_affine$bGR[,1,2] - Post_model_GR_affine$bGR[,2,2] ) , 
                   c( Post_model_GR_affine$bGR[,1,3] - Post_model_GR_affine$bGR[,2,3] ) , 
                   
                   c( Post_model_GR_fri$bGR[,1,1] - Post_model_GR_fri$bGR[,2,1] ) , 
                   c( Post_model_GR_fri$bGR[,1,2] - Post_model_GR_fri$bGR[,2,2] ) , 
                   c( Post_model_GR_fri$bGR[,1,3] - Post_model_GR_fri$bGR[,2,3] ) ) ,
        Residence = c(rep( c( rep( "Matrilocal" , 6000 ) , 
                              rep( "Bilocal" , 6000 ) , 
                              rep( "Patrilocal" , 6000 ) ) , 4 ) ) ,
        Type = c( rep( "All" , 18000 ) ,
                  rep( "Biological kin" , 18000 ) ,
                  rep( "Affinal kin" , 18000 ) ,
                  rep( "Friend" , 18000 ) ) ) %>% 
  mutate( Residence = factor( Residence , 
                              levels = c( "Matrilocal" , "Bilocal",  "Patrilocal" ) , 
                              labels = c( "Matrilocal" , "Bilocal",  "Patrilocal" ) ) , 
          Type = factor( Type , 
                         levels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) , 
                         labels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) ) %>% 
  ggplot( aes( x = Value , y = Residence , 
               fill = Residence , color = Residence) ) +
  facet_wrap( ~ Type ) +
  geom_vline( xintercept = 0, linetype = 2 , linewidth = 1 ,
              color = "dimgray" ) +
  stat_halfeye( .width = .90, height = 1.0 ) + 
  labs( title = "Gender difference" ) +
  scale_fill_manual( values = alpha( c( "#c5272d" , "#037f77" , "#0001a1" ) , 0.3 ) ) +
  scale_color_manual( values = c( "#c5272d" , "#037f77" , "#0001a1" ) ) +
  scale_y_discrete( NULL, labels = ggplot2:::parse_safe ) +
  xlab( "Posterior estimates" ) +
  coord_cartesian( ylim = c( 1.5, 3.6 ) ) +
  scale_x_continuous( limits = c( -1 , 1.5 ) , breaks = c( -1 , 0 , 1 ) ) +
  theme(plot.margin = margin(10, 10, 10, 20),
        strip.background = element_rect(color = "black", fill = "white") ,
        strip.text.x = element_text(size = 16),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "Gainsboro"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black",fill=NA),
        panel.background = element_blank(),
        plot.title = element_text(size=16),
        legend.position = "none",
        axis.title.x = element_text(size = 16,
                                    margin = margin(t = 0.3, r = 0, b = 0, l = 0,unit = "cm")),
        axis.text.x = element_text(colour="black",size=14,
                                   margin = margin(t = 0.1, r = 0, b = 0, l = 0,unit = "cm")),
        axis.title.y = element_text(size = 16,
                                    margin = margin(t = 0, r = 0.5, b = 0, l = 0,unit = "cm")),
        axis.text.y = element_text(colour="black",size=14,
                                   margin = margin(t = 0, r = 0.2, b = 0, l = 0,unit = "cm")),
        legend.title=element_text(size=14,face = "bold"),
        legend.text=element_text(size=12,face = "bold"))

ggsave( filename = "Gender_diff_by_residence.jpeg" , 
        width = 150 , height = 120 , units = "mm" , dpi = 300 )

# Proportions of mean difference in βRvalues comparing women to men
tibble( Type = c( rep( "All" , 3 ) ,
                  rep( "Biological kin" , 3 ) ,
                  rep( "Affinal kin" , 3 ) ,
                  rep( "Friend" , 3 ) ) ,
        Residence = c(rep( c( "Matrilocal" , 
                              "Bilocal" , 
                              "Patrilocal" ) , 4 ) ) ,
        `Gender difference` = c( sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_all$bGR[,1,1] - Post_model_GR_all$bGR[,2,1] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_all$bGR[,1,2] - Post_model_GR_all$bGR[,2,2] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_all$bGR[,1,3] - Post_model_GR_all$bGR[,2,3] ) > 0 ) ) /6000 ) , 
                                 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_genetic$bGR[,1,1] - Post_model_GR_genetic$bGR[,2,1] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_genetic$bGR[,1,2] - Post_model_GR_genetic$bGR[,2,2] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_genetic$bGR[,1,3] - Post_model_GR_genetic$bGR[,2,3] ) > 0 ) ) /6000 ) , 
                                 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_affine$bGR[,1,1] - Post_model_GR_affine$bGR[,2,1] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_affine$bGR[,1,2] - Post_model_GR_affine$bGR[,2,2] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_affine$bGR[,1,3] - Post_model_GR_affine$bGR[,2,3] ) > 0 ) ) /6000 ) , 
                                 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_fri$bGR[,1,1] - Post_model_GR_fri$bGR[,2,1] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_fri$bGR[,1,2] - Post_model_GR_fri$bGR[,2,2] ) > 0 ) ) /6000 ) , 
                                 sprintf("%0.4f" , 
                                         length( which( c( Post_model_GR_fri$bGR[,1,3] - Post_model_GR_fri$bGR[,2,3] ) > 0 ) ) /6000 ) ) ) %>% 
  mutate( 
    `Gender difference` = as.numeric( `Gender difference` ) ,
    `Gender difference` = paste0( format( `Gender difference` * 100 , nsmall = 2) , "%" ) ) %>% 
  gt() %>% 
  fmt_number( ) %>% 
  sub_missing( columns = everything() , 
               missing_text = " " ) %>% 
  tab_style( style = list( cell_text( weight = "bold" ) ) ,
             locations = cells_column_labels(everything() ) ) %>% 
  tab_style( style = list( cell_text(font = "Times New Roman") ) ,
             locations = cells_body() ) %>% 
  cols_align( align = "center",
              columns = c( "Type" , "Residence" , `Gender difference` ) ) %>% 
  opt_table_lines( extent = "default" ) %>%
  tab_options( column_labels.border.top.color = "black" ,
               column_labels.border.top.width = px( 3 ) ,
               column_labels.border.bottom.color = "black" ,
               table_body.hlines.color = "lightgrey" ,
               table.border.bottom.color = "black" ,
               table.border.bottom.width = px( 3 ) ) %>% 
  cols_width( everything() ~ px( 150 ) )

###3.2.6 Outputs ----
Output_model_GR <- data.frame( Mean = c( Pre_model_GR_all$mean ,
                                         Pre_model_GR_genetic$mean ,
                                         Pre_model_GR_affine$mean ,
                                         Pre_model_GR_fri$mean ) ,
                               CI5 = c( Pre_model_GR_all$`5%` ,
                                        Pre_model_GR_genetic$`5%` ,
                                        Pre_model_GR_affine$`5%` ,
                                        Pre_model_GR_fri$`5%` ) ,
                               CI95 = c( Pre_model_GR_all$`95%` ,
                                         Pre_model_GR_genetic$`95%` ,
                                         Pre_model_GR_affine$`95%` ,
                                         Pre_model_GR_fri$`95%` ) ,
                               Type = rep( c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) , each = 10 ) ,
                               Variable = rep( c( "Age[female]" , "Age[male]" , 
                                                  "Age^2[female]" , "Age^2[male]" ,
                                                  "Female: matrilocal" , 
                                                  "Male: matrilocal" , 
                                                  "Female: bilocal" , 
                                                  "Male: bilocal" , 
                                                  "Female: patrilocal" , 
                                                  "Male: patrilocal" ) , 4 ) ) %>% 
  mutate( Type = factor( Type , 
                         levels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ,
                         labels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) ,
          Variable = factor( Variable , 
                             levels = c( "Age[female]" , "Age[male]" , 
                                         "Age^2[female]" , "Age^2[male]" ,
                                         "Female: matrilocal" , 
                                         "Female: bilocal" , 
                                         "Female: patrilocal" , 
                                         "Male: matrilocal" , 
                                         "Male: bilocal" , 
                                         "Male: patrilocal" ) ,
                             labels = c( "Age[female]" , "Age[male]" , 
                                         "Age^2[female]" , "Age^2[male]" ,
                                         "Female: matrilocal" , 
                                         "Female: bilocal" , 
                                         "Female: patrilocal" , 
                                         "Male: matrilocal" , 
                                         "Male: bilocal" , 
                                         "Male: patrilocal" ) ) ,
          Model = "Residence" , 
          "Mean [CI]" = paste( format( round( Mean , 2 ), nsmall = 2 ) , 
                               " [" , 
                               format( round( CI5 , 2 ), nsmall = 2 ) , 
                               ", " , 
                               format( round( CI95 , 2 ), nsmall = 2 ) , 
                               "]", 
                               sep = "") )

Res.gender.output <- tibble( Variable = c( "Age[female]" , "Age[male]" , 
                                           "Age^2[female]" , "Age^2[male]" ,
                                           "Female: matrilocal" , 
                                           "Female: bilocal" , 
                                           "Female: patrilocal" , 
                                           "Male: matrilocal" , 
                                           "Male: bilocal" , 
                                           "Male: patrilocal" ) ,
                             Type = "All" ) %>% 
  left_join( Output_model_GR[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `All` = "Mean [CI]" ) %>% 
  mutate( Type = "Biological kin" ) %>% 
  left_join( Output_model_GR[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `Biological kin` = "Mean [CI]" ) %>% 
  mutate( Type = "Affinal kin" ) %>% 
  left_join( Output_model_GR[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `Affinal kin` = "Mean [CI]" ) %>% 
  mutate( Type = "Friend" ) %>% 
  left_join( Output_model_GR[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `Friend` = "Mean [CI]" ) %>% 
  select( - Type ) %>% 
  gt() %>% 
  fmt_number( ) %>% 
  sub_missing( columns = everything() , 
               missing_text = " " ) %>% 
  tab_style( style = list( cell_text( weight = "bold" ) ) ,
             locations = cells_column_labels( everything() ) ) %>% 
  tab_style( style = list( cell_text( font = "Times New Roman" ) ) ,
             locations = cells_body() ) %>% 
  cols_align( align = "left",
              columns = "Variable" ) %>%
  cols_align( align = "center",
              columns = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) %>% 
  opt_table_lines( extent = "default" ) %>%
  tab_options( column_labels.border.top.color = "black" ,
               column_labels.border.top.width = px( 3 ) ,
               column_labels.border.bottom.color = "black" ,
               table_body.hlines.color = "lightgrey" ,
               table.border.bottom.color = "black" ,
               table.border.bottom.width = px( 3 ) ) %>% 
  cols_width( "Variable" ~ px( 170 ) , 
              everything() ~ px( 150 ) )
Res.gender.output

##3.3 Market integration----
# Other variable should be included into models
adjustmentSets( Dag , 
                exposure = "Market" , 
                outcome = "Density" ,
                effect = "total" )

###3.3.1 Models of all density ----
Model_GM_all_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_all = as.integer( Tie_all ) , 
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Age2 = as.numeric( sprintf( "%0.4f" , Age_std_2 ) ) , 
  Gender = as.integer( Gender ) , 
  PMR = as.integer( Residence ) , 
  Market = as.numeric( sprintf("%0.4f" , Market_log_std ) ) ,
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

# Stratify by gender
{set.seed(123)
  Model_GM_all <- ulam(
    alist(
      Tie_all ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + bGA2[Gender] * Age2 + 
        bGR[Gender,PMR] + 
        bM[Gender] * Market +
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      bGA2[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      bM[Gender] ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 )  
    ) , data = Model_GM_all_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}

Pre_model_GM_all <- precis( Model_GM_all , depth = 3 , prob = 0.90 , 
                            pars = c( "bGA" , "bGA2" , "bGR" , "bM" ) )
Pre_model_GM_all

# Stratify by gender and village
{set.seed(123)
  Model_GVM_all <- ulam(
    alist(
      Tie_all ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + bGA2[Gender] * Age + 
        bGR[Gender,PMR] + 
        bGVM[Gender,VID] * Market +
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      bGA2[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,VID]: bGVM ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 )  
    ) , data = Model_GM_all_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}

Pre_model_GVM_all <- precis( Model_GVM_all , depth = 3 , prob = 0.90 , 
                             pars = c( "bGA" , "bGA2" , "bGR" , "bGVM" , 
                                       "sigma_V" , "sigma_I" ) )
Pre_model_GVM_all

###3.3.2 Models of biological kin density ----
Model_GM_genetic_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_genetic = as.integer( Tie_genetic ) ,
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Age2 = as.numeric( sprintf( "%0.4f" , Age_std_2 ) ) , 
  Gender = as.integer( Gender ) , 
  PMR = as.integer( Residence ) , 
  Market = as.numeric( sprintf("%0.4f" , Market_log_std ) ) ,
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

# Stratify by gender
{set.seed(123)
  Model_GM_genetic <- ulam(
    alist(
      Tie_genetic ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + bGA2[Gender] * Age2 + 
        bGR[Gender,PMR] + 
        bM[Gender] * Market +
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      bGA2[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      bM[Gender] ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 ) 
    ) , data = Model_GM_genetic_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GM_genetic <- precis( Model_GM_genetic , depth = 3 , prob = 0.90 , 
                                pars = c( "bGA" , "bGA2" , "bGR" , "bM" ) )
Pre_model_GM_genetic

# Stratify by gender and village
{set.seed(123)
  Model_GVM_genetic <- ulam(
    alist(
      Tie_genetic ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + bGA2[Gender] * Age + 
        bGR[Gender,PMR] + 
        bGVM[Gender,VID] * Market +
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      bGA2[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,VID]: bGVM ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 )     
    ) , data = Model_GM_genetic_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GVM_genetic <- precis( Model_GVM_genetic , depth = 3 , prob = 0.90 , 
                                 pars = c( "bGA" , "bGA2" , "bGR" , "bGVM" , 
                                           "sigma_V" , "sigma_I" ) )
Pre_model_GVM_genetic

###3.3.3 Models of affinal kin density ----
Model_GM_affine_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_affine = as.integer( Tie_affine ) ,
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Age2 = as.numeric( sprintf( "%0.4f" , Age_std_2 ) ) , 
  Gender = as.integer( Gender ) , 
  PMR = as.integer( Residence ) , 
  Market = as.numeric( sprintf("%0.4f" , Market_log_std ) ) ,
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

# Stratify by gender
{set.seed(123)
  Model_GM_affine <- ulam(
    alist(
      Tie_affine ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + bGA2[Gender] * Age2 + 
        bGR[Gender,PMR] + 
        bM[Gender] * Market +
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      bGA2[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      bM[Gender] ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 )
    ) , data = Model_GM_affine_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GM_affine <- precis( Model_GM_affine , depth = 3 , prob = 0.90 , 
                               pars = c( "bGA" , "bGA2" , "bGR" , "bM" ) )
Pre_model_GM_affine

# Stratify by gender and village
{set.seed(123)
  Model_GVM_affine <- ulam(
    alist(
      Tie_affine ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + bGA2[Gender] * Age + 
        bGR[Gender,PMR] + 
        bGVM[Gender,VID] * Market +
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      bGA2[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,VID]: bGVM ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 ) 
    ) , data = Model_GM_affine_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GVM_affine <- precis( Model_GVM_affine , depth = 3 , prob = 0.90 , 
                                pars = c( "bGA" , "bGA2" , "bGR" , "bGVM" , 
                                          "sigma_V" , "sigma_I" ) )
Pre_model_GVM_affine

###3.3.4 Models of friend density ----
Model_GM_fri_list <- with( Ego.data , list(
  Possible_tie =  as.integer( Possible_tie ) , 
  Tie_fri = as.integer( Tie_fri ) ,
  Age = as.numeric( sprintf( "%0.4f" , Age_std ) ) ,
  Age2 = as.numeric( sprintf( "%0.4f" , Age_std_2 ) ) , 
  Gender = as.integer( Gender ) , 
  PMR = as.integer( Residence ) , 
  Market = as.numeric( sprintf("%0.4f" , Market_log_std ) ) ,
  ID = as.integer( ID ) ,
  VID = as.integer( VID ) ) )

# Stratify by gender
{set.seed(123)
  Model_GM_fri <- ulam(
    alist(
      Tie_fri ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + bGA2[Gender] * Age2 + 
        bGR[Gender,PMR] + 
        bM[Gender] * Market +
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      bGA2[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      bM[Gender] ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 )
    ) , data = Model_GM_fri_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GM_fri <- precis( Model_GM_fri , depth = 3 , prob = 0.90 , 
                            pars = c( "bGA" , "bGA2" , "bGR" , "bM" ) )
Pre_model_GM_fri

# Stratify by gender and village
{set.seed(123)
  Model_GVM_fri <- ulam(
    alist(
      Tie_fri ~ binomial( Possible_tie , p ) ,
      logit( p ) <- bGA[Gender] * Age + bGA2[Gender] * Age + 
        bGR[Gender,PMR] + 
        bGVM[Gender,VID] * Market +
        V[VID] + # village intercepts 
        I[ID] , # individual intercepts 
      
      # define effects using other parameters
      transpars> vector[VID]: V <<- V_bar + z_V*sigma_V,
      transpars> vector[ID]: I <<- z_I*sigma_I,
      
      bGA[Gender] ~ normal( 0 , 0.5 ) ,
      bGA2[Gender] ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,PMR]: bGR ~ normal( 0 , 0.5 ) ,
      
      matrix[Gender,VID]: bGVM ~ normal( 0 , 0.5 ) ,
      
      V_bar ~ normal( 0 , 0.5 ) ,
      z_V[VID] ~ normal( 0 , 1 ) ,
      sigma_V ~ exponential( 1 ) ,
      z_I[ID] ~ normal( 0 , 4 ) ,
      sigma_I ~ exponential( 1 ) 
    ) , data = Model_GM_fri_list , 
    iter = 1500 , warmup = 500 , chains = 6 , cores = 6 , 
    log_lik = TRUE , control = list(adapt_delta = 0.99)
  )}
Pre_model_GVM_fri <- precis( Model_GVM_fri , depth = 3 , prob = 0.90 , 
                             pars = c( "bGA" , "bGA2" , "bGR" , "bGVM" , 
                                       "sigma_V" , "sigma_I" ) )
Pre_model_GVM_fri

###3.3.5 Figures ----
# Posterior samples of estimates
Post_model_GM_all <- extract.samples( Model_GM_all )
Post_model_GM_genetic <- extract.samples( Model_GM_genetic )
Post_model_GM_affine <- extract.samples( Model_GM_affine )
Post_model_GM_fri <- extract.samples( Model_GM_fri )

Post_model_GVM_all <- extract.samples( Model_GVM_all )
Post_model_GVM_genetic <- extract.samples( Model_GVM_genetic )
Post_model_GVM_affine <- extract.samples( Model_GVM_affine )
Post_model_GVM_fri <- extract.samples( Model_GVM_fri )

####3.3.5.1 Effects of market integration by gender ----
# Function to predict p-matrix based on posterior samples of estimates and new data
predict_p_matrix_GM <- function( post_samples , new_data , n_samples = 6000 ) {
  n_market <- length( unique( new_data$N_market ) )  # Number of markets
  p_matrix <- matrix( NA , nrow = n_samples , ncol = n_market )  # Initialize matrix
  
  # Loop over posterior samples
  for (s in 1:n_samples) {
    
    # Extract posterior sample for current iteration
    bGA <- post_samples$bGA[s, ]
    bGA2 <- post_samples$bGA2[s, ]
    bGR <- post_samples$bGR[s, , ]
    bM <- post_samples$bM[s, ]
    V_bar <- post_samples$V_bar[s, 1]
    
    # Loop over each market
    for (m in 1:n_market) {
      market_subset <- new_data[new_data$N_market == m, ]  # Subset data for market
      logit_p_vals <- with( market_subset,
                            bGA[Gender] * Age + bGA2[Gender] * Age2 + 
                              bGR[cbind(Gender, PMR)] + 
                              bM[Gender] * Market + 
                              V_bar)
      
      # Convert to probabilities and average over Age and PMR
      p_matrix[s, m] <- mean( rethinking::inv_logit( logit_p_vals ) , na.rm = TRUE)
    }
  }
  
  return(p_matrix)
}

# New data for women and men, respectively 
New_data_GM_F <- data.frame( Age_act = as.integer( rep( 20:86 , 222 ) ) , 
                             Age = rep( ( as.integer( 20:86 ) - mean( Ego.data$Age ) )/ sd( Ego.data$Age ) , 222) , 
                             Age2 = rep( ( ( as.integer( 20:86 ) - mean( Ego.data$Age ) )/ sd( Ego.data$Age ) )^2 , 222 ) , 
                             Gender = as.integer( rep( 1 , 14874 ) ) ,
                             PMR = as.integer( rep( 1:3 , 4958 ) ) ,
                             N_market = rep( 1:74 , 201 ) ,
                             Market_act = rep( seq( 0 , 730 , 10 ) , 201 ) ,
                             Market = rep( c( ( log( c( seq( 0 , 730 , 10 ) + 1 ) ) - mean( Ego.data$Market_log ) ) / sd( Ego.data$Market_log ) ) ) , 201 ) 

New_data_GM_M <- data.frame( Age_act = as.integer( rep( 20:86 , 222 ) ), 
                             Age = rep( ( as.integer( 20:86 ) - mean( Ego.data$Age ) )/ sd( Ego.data$Age ) , 222) , 
                             Age2 = rep( ( ( as.integer( 20:86 ) - mean( Ego.data$Age ) )/ sd( Ego.data$Age ) )^2 , 222 ) , 
                             Gender = as.integer( rep( 2 , 14874 ) ) ,
                             PMR = as.integer( rep( 1:3 , 4958 ) ) ,
                             N_market = rep( 1:74 , 201 ) ,
                             Market_act = rep( seq( 0 , 730 , 10 ) , 201 ) ,
                             Market = rep( c( ( log( c( seq( 0 , 730 , 10 ) + 1 ) ) - mean( Ego.data$Market_log ) ) / sd( Ego.data$Market_log ) ) ) , 201 ) 

# Predicted data based on new data
p_matrix_GM_all_F <- predict_p_matrix_GM( post_samples = Post_model_GM_all , 
                                          new_data = New_data_GM_F, 
                                          n_samples = 6000 )
p_matrix_GM_all_M <- predict_p_matrix_GM( post_samples = Post_model_GM_all , 
                                          new_data = New_data_GM_M, 
                                          n_samples = 6000 )

p_matrix_GM_genetic_F <- predict_p_matrix_GM( post_samples = Post_model_GM_genetic , 
                                              new_data = New_data_GM_F, 
                                              n_samples = 6000 )
p_matrix_GM_genetic_M <- predict_p_matrix_GM( post_samples = Post_model_GM_genetic , 
                                              new_data = New_data_GM_M, 
                                              n_samples = 6000 )

p_matrix_GM_affine_F <- predict_p_matrix_GM( post_samples = Post_model_GM_affine , 
                                             new_data = New_data_GM_F, 
                                             n_samples = 6000 )
p_matrix_GM_affine_M <- predict_p_matrix_GM( post_samples = Post_model_GM_affine , 
                                             new_data = New_data_GM_M, 
                                             n_samples = 6000 )

p_matrix_GM_fri_F <- predict_p_matrix_GM( post_samples = Post_model_GM_fri , 
                                          new_data = New_data_GM_F, 
                                          n_samples = 6000 )
p_matrix_GM_fri_M <- predict_p_matrix_GM( post_samples = Post_model_GM_fri , 
                                          new_data = New_data_GM_M, 
                                          n_samples = 6000 )

#####3.3.5.1.1 Gender-specific effects of market integration ----
Market_seq <- c( ( log( c( seq( 0 , 730 , 10 ) + 1 ) ) - mean( Ego.data$Market_log ) ) / sd( Ego.data$Market_log ) )

jpeg( "Market_stratified_by_gender.jpg" , 
      width = 200 , height = 160 , units = "mm" , res = 300 )

{par(mfrow=c(2,2), mar = c(4.5, 4.5, 1.8, 1))
  # All
  plot( NULL,
        xlim = c(-2.5, 2.2),
        ylim = c(0,1),
        xlab = "" ,
        ylab = "Predicted density" ,
        main = "All" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 ,
        cex.main = 1.5 )
  # women
  lines( Market_seq , 
         apply( p_matrix_GM_all_F , 2 , mean ) ,
         lwd = 2 , col = "#CB2313" )
  shade( apply( p_matrix_GM_all_F , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GM_all_F , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GM_all_F , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  # men
  lines( Market_seq , 
         apply( p_matrix_GM_all_M , 2 , mean ) ,
         lwd = 2 , col = "#046C9A" )
  shade( apply( p_matrix_GM_all_M , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GM_all_M , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GM_all_M , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  
  # Biological kin
  plot( NULL,
        xlim = c(-2.5, 2.2),
        ylim = c(0,1),
        xlab = "" ,
        ylab = "" ,
        main = "Biological kin" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 ,
        cex.main = 1.5 )
  legend( x = 0.8 , y = 1.0 ,  
          box.col = "white",
          legend = c( "Female" , "Male" ) , 
          lty = c( 1 , 1 ) , 
          col = c( "#CB2313" , "#046C9A" ) , 
          lwd = 2 ,
          cex = 1.2 , 
          bty = "n" ,
          y.intersp = 1.2 ,
          x.intersp = 0.3 ,
          seg.len = 0.8  )
  
  # women
  lines( Market_seq , 
         apply( p_matrix_GM_genetic_F , 2 , mean ) ,
         lwd = 2 , col = "#CB2313" )
  shade( apply( p_matrix_GM_genetic_F , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GM_genetic_F , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GM_genetic_F , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  # men
  lines( Market_seq , 
         apply( p_matrix_GM_genetic_M , 2 , mean ) ,
         lwd = 2 , col = "#046C9A" )
  shade( apply( p_matrix_GM_genetic_M , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GM_genetic_M , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GM_genetic_M , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  
  # Affinal kin
  plot( NULL,
        xlim = c(-2.5, 2.2),
        ylim = c(0,1),
        xlab = "Market integration (std)" ,
        ylab = "Predicted density" ,
        main = "Affinal kin" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 ,
        cex.main = 1.5 )
  # women
  lines( Market_seq , 
         apply( p_matrix_GM_affine_F , 2 , mean ) ,
         lwd = 2 , col = "#CB2313" )
  shade( apply( p_matrix_GM_affine_F , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GM_affine_F , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GM_affine_F , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  # men
  lines( Market_seq , 
         apply( p_matrix_GM_affine_M , 2 , mean ) ,
         lwd = 2 , col = "#046C9A" )
  shade( apply( p_matrix_GM_affine_M , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GM_affine_M , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GM_affine_M , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  
  # Friend kin
  plot( NULL,
        xlim = c(-2.5,2.2),
        ylim = c(0,1),
        xlab = "Market integration (std)" ,
        ylab = "" ,
        main = "Friend" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 ,
        cex.main = 1.5 )
  # women
  lines( Market_seq , 
         apply( p_matrix_GM_fri_F , 2 , mean ) ,
         lwd = 2 , col = "#CB2313" )
  shade( apply( p_matrix_GM_fri_F , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GM_fri_F , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  shade( apply( p_matrix_GM_fri_F , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "#CB2313" , 0.15 ) )
  # men
  lines( Market_seq , 
         apply( p_matrix_GM_fri_M , 2 , mean ) ,
         lwd = 2 , col = "#046C9A" )
  shade( apply( p_matrix_GM_fri_M , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GM_fri_M , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
  shade( apply( p_matrix_GM_fri_M , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "#046C9A" , 0.15 ) )
}

dev.off()

#####3.3.5.1.2 Gender difference in effects of market integration ----
jpeg( "Market_gender_diff.jpg" , 
      width = 200 , height = 160 , units = "mm" , res = 300 )

{par(mfrow=c(2,2), mar = c(4.5, 4.5, 1.8, 1))
  # All
  plot( NULL,
        xlim = c(-2.5, 2.2),
        ylim = c( -0.4 , 0.4 ),
        xlab = "" ,
        ylab = "Predicted density" ,
        main = "All" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 ,
        cex.main = 1.5 )
  abline( h = 0 , lty = 2 , lwd = 2, col = 1 )
  lines( Market_seq , 
         apply( p_matrix_GM_all_F - p_matrix_GM_all_M , 2 , mean ) ,
         lwd = 2 , col = "black" )
  shade( apply( p_matrix_GM_all_F - p_matrix_GM_all_M , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  shade( apply( p_matrix_GM_all_F - p_matrix_GM_all_M , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  shade( apply( p_matrix_GM_all_F - p_matrix_GM_all_M , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  
  # Genetic kin
  plot( NULL,
        xlim = c(-2.5, 2.2),
        ylim = c( -0.4 , 0.4 ),
        xlab = "" ,
        ylab = "" ,
        main = "Biological kin" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 ,
        cex.main = 1.5 )
  abline( h = 0 , lty = 2 , lwd = 2, col = 1 )
  lines( Market_seq , 
         apply( p_matrix_GM_genetic_F - p_matrix_GM_genetic_M , 2 , mean ) ,
         lwd = 2 , col = "black" )
  shade( apply( p_matrix_GM_genetic_F - p_matrix_GM_genetic_M , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  shade( apply( p_matrix_GM_genetic_F - p_matrix_GM_genetic_M , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  shade( apply( p_matrix_GM_genetic_F - p_matrix_GM_genetic_M , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  
  # Affinal kin
  plot( NULL,
        xlim = c(-2.5, 2.2),
        ylim = c( -0.4 , 0.4 ),
        xlab = "Market integration (std)" ,
        ylab = "Predicted density" ,
        main = "Affinal kin" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 ,
        cex.main = 1.5 )
  abline( h = 0 , lty = 2 , lwd = 2, col = 1 )
  lines( Market_seq , 
         apply( p_matrix_GM_affine_F - p_matrix_GM_affine_M , 2 , mean ) ,
         lwd = 2 , col = "black" )
  shade( apply( p_matrix_GM_affine_F - p_matrix_GM_affine_M , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  shade( apply( p_matrix_GM_affine_F - p_matrix_GM_affine_M , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  shade( apply( p_matrix_GM_affine_F - p_matrix_GM_affine_M , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  
  # Friend kin
  plot( NULL,
        xlim = c(-2.5,2.2),
        ylim = c( -0.4 , 0.4 ),
        xlab = "Market integration (std)" ,
        ylab = "" ,
        main = "Friend" ,
        cex.axis = 1.2 ,
        cex.lab = 1.5 ,
        cex.main = 1.5 )
  abline( h = 0 , lty = 2 , lwd = 2, col = 1 )
  lines( Market_seq , 
         apply( p_matrix_GM_fri_F - p_matrix_GM_fri_M , 2 , mean ) ,
         lwd = 2 , col = "black" )
  shade( apply( p_matrix_GM_fri_F - p_matrix_GM_fri_M , 2 , PI , prob = 0.90 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  shade( apply( p_matrix_GM_fri_F - p_matrix_GM_fri_M , 2 , PI , prob = 0.60 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
  shade( apply( p_matrix_GM_fri_F - p_matrix_GM_fri_M , 2 , PI , prob = 0.30 ) , 
         Market_seq , 
         col = col.alpha( "dimgray" , 0.15 ) )
}

dev.off()

####3.3.5.2 Effects of market integration by gender and village ----
# Function to predict p-matrix based on posterior samples of estimates and new data
generate_and_predict_p_matrix_GVM <- function( post_samples , Gender_value = 1, VID_value = 1, n_samples = 6000 ) {
  
  # Generate new data
  new_data <- data.frame(
    Age = rep( ( as.integer( 20:86 ) - mean( Ego.data$Age ) )/ sd( Ego.data$Age ) , 222 ) , 
    Age2 = rep( ( ( as.integer( 20:86 ) - mean( Ego.data$Age ) )/ sd( Ego.data$Age ) )^2 , 222 ) , 
    Gender = as.integer( rep( Gender_value , 14874 ) ) ,
    PMR = as.integer( rep( 1:3 , 4958 ) ) , 
    N_market = rep( 1:74 , 201 ) , 
    Market = rep(( log( seq( 0 , 730 , 10 ) + 1 ) - mean( Ego.data$Market_log ) ) / sd( Ego.data$Market_log ) , 201 ) , 
    VID = as.integer( rep( VID_value, 14874 ) )
  )
  
  # Initialize p_matrix
  n_market <- length( unique( new_data$N_market ) )
  p_matrix <- matrix( NA , nrow = n_samples , ncol = n_market )
  
  # Loop over posterior samples
  for ( s in 1:n_samples ) {
    
    # Extract posterior sample for current iteration
    bGA <- post_samples$bGA[s, ]
    bGA2 <- post_samples$bGA2[s, ]
    bGR <- post_samples$bGR[s, , ]
    bGVM <- post_samples$bGVM[s, , ]
    V <- post_samples$V[s, ]
    
    # Loop over each market
    for (m in 1:n_market) {
      market_subset <- new_data[new_data$N_market == m, ]  # Subset data for market
      logit_p_vals <- with(market_subset,
                           bGA[Gender] * Age + bGA2[Gender] * Age2 + 
                             bGR[cbind(Gender, PMR)] + 
                             bGVM[cbind(Gender, VID)] * Market +
                             V[VID] )
      
      # Convert to probabilities and average over Age and PMR
      p_matrix[s, m] <- mean(rethinking::inv_logit(logit_p_vals), na.rm = TRUE)
    }
  }
  
  return(p_matrix)
}

# Predicted data for women based on new data
# All network density
p_matrix_GVM_all_F_V1 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 1, VID_value = 1 )
p_matrix_GVM_all_F_V2 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 1, VID_value = 2 )
p_matrix_GVM_all_F_V3 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 1, VID_value = 3 )
p_matrix_GVM_all_F_V4 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 1, VID_value = 4 )
p_matrix_GVM_all_F_V5 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 1, VID_value = 5 )
p_matrix_GVM_all_F_V6 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 1, VID_value = 6 )
p_matrix_GVM_all_F_V7 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 1, VID_value = 7 )
p_matrix_GVM_all_F_V8 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 1, VID_value = 8 )
p_matrix_GVM_all_F_V9 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 1, VID_value = 9 )
p_matrix_GVM_all_F_V10 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                             Gender_value = 1, VID_value = 10 )
p_matrix_GVM_all_F_V11 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                             Gender_value = 1, VID_value = 11 )
p_matrix_GVM_all_F_V12 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                             Gender_value = 1, VID_value = 12 )
p_matrix_GVM_all_F_V13 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                             Gender_value = 1, VID_value = 13 )
p_matrix_GVM_all_F_V14 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                             Gender_value = 1, VID_value = 14 )

# Biological kin density
p_matrix_GVM_genetic_F_V1 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 1, VID_value = 1 )
p_matrix_GVM_genetic_F_V2 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 1, VID_value = 2 )
p_matrix_GVM_genetic_F_V3 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 1, VID_value = 3 )
p_matrix_GVM_genetic_F_V4 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 1, VID_value = 4 )
p_matrix_GVM_genetic_F_V5 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 1, VID_value = 5 )
p_matrix_GVM_genetic_F_V6 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 1, VID_value = 6 )
p_matrix_GVM_genetic_F_V7 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 1, VID_value = 7 )
p_matrix_GVM_genetic_F_V8 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 1, VID_value = 8 )
p_matrix_GVM_genetic_F_V9 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 1, VID_value = 9 )
p_matrix_GVM_genetic_F_V10 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                 Gender_value = 1, VID_value = 10 )
p_matrix_GVM_genetic_F_V11 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                 Gender_value = 1, VID_value = 11 )
p_matrix_GVM_genetic_F_V12 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                 Gender_value = 1, VID_value = 12 )
p_matrix_GVM_genetic_F_V13 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                 Gender_value = 1, VID_value = 13 )
p_matrix_GVM_genetic_F_V14 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                 Gender_value = 1, VID_value = 14 )

# Affinal kin density
p_matrix_GVM_affine_F_V1 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 1, VID_value = 1 )
p_matrix_GVM_affine_F_V2 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 1, VID_value = 2 )
p_matrix_GVM_affine_F_V3 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 1, VID_value = 3 )
p_matrix_GVM_affine_F_V4 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 1, VID_value = 4 )
p_matrix_GVM_affine_F_V5 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 1, VID_value = 5 )
p_matrix_GVM_affine_F_V6 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 1, VID_value = 6 )
p_matrix_GVM_affine_F_V7 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 1, VID_value = 7 )
p_matrix_GVM_affine_F_V8 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 1, VID_value = 8 )
p_matrix_GVM_affine_F_V9 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 1, VID_value = 9 )
p_matrix_GVM_affine_F_V10 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                                Gender_value = 1, VID_value = 10 )
p_matrix_GVM_affine_F_V11 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                                Gender_value = 1, VID_value = 11 )
p_matrix_GVM_affine_F_V12 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                                Gender_value = 1, VID_value = 12 )
p_matrix_GVM_affine_F_V13 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                                Gender_value = 1, VID_value = 13 )
p_matrix_GVM_affine_F_V14 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                                Gender_value = 1, VID_value = 14 )

# Friend density
p_matrix_GVM_fri_F_V1 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 1, VID_value = 1 )
p_matrix_GVM_fri_F_V2 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 1, VID_value = 2 )
p_matrix_GVM_fri_F_V3 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 1, VID_value = 3 )
p_matrix_GVM_fri_F_V4 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 1, VID_value = 4 )
p_matrix_GVM_fri_F_V5 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 1, VID_value = 5 )
p_matrix_GVM_fri_F_V6 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 1, VID_value = 6 )
p_matrix_GVM_fri_F_V7 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 1, VID_value = 7 )
p_matrix_GVM_fri_F_V8 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 1, VID_value = 8 )
p_matrix_GVM_fri_F_V9 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 1, VID_value = 9 )
p_matrix_GVM_fri_F_V10 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                             Gender_value = 1, VID_value = 10 )
p_matrix_GVM_fri_F_V11 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                             Gender_value = 1, VID_value = 11 )
p_matrix_GVM_fri_F_V12 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                             Gender_value = 1, VID_value = 12 )
p_matrix_GVM_fri_F_V13 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                             Gender_value = 1, VID_value = 13 )
p_matrix_GVM_fri_F_V14 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                             Gender_value = 1, VID_value = 14 )

# Predicted data for men based on new data
# All network density
p_matrix_GVM_all_M_V1 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 2, VID_value = 1 )
p_matrix_GVM_all_M_V2 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 2, VID_value = 2 )
p_matrix_GVM_all_M_V3 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 2, VID_value = 3 )
p_matrix_GVM_all_M_V4 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 2, VID_value = 4 )
p_matrix_GVM_all_M_V5 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 2, VID_value = 5 )
p_matrix_GVM_all_M_V6 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 2, VID_value = 6 )
p_matrix_GVM_all_M_V7 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 2, VID_value = 7 )
p_matrix_GVM_all_M_V8 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 2, VID_value = 8 )
p_matrix_GVM_all_M_V9 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                            Gender_value = 2, VID_value = 9 )
p_matrix_GVM_all_M_V10 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                             Gender_value = 2, VID_value = 10 )
p_matrix_GVM_all_M_V11 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                             Gender_value = 2, VID_value = 11 )
p_matrix_GVM_all_M_V12 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                             Gender_value = 2, VID_value = 12 )
p_matrix_GVM_all_M_V13 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                             Gender_value = 2, VID_value = 13 )
p_matrix_GVM_all_M_V14 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_all , 
                                                             Gender_value = 2, VID_value = 14 )

# Biological kin density
p_matrix_GVM_genetic_M_V1 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 2, VID_value = 1 )
p_matrix_GVM_genetic_M_V2 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 2, VID_value = 2 )
p_matrix_GVM_genetic_M_V3 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 2, VID_value = 3 )
p_matrix_GVM_genetic_M_V4 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 2, VID_value = 4 )
p_matrix_GVM_genetic_M_V5 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 2, VID_value = 5 )
p_matrix_GVM_genetic_M_V6 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 2, VID_value = 6 )
p_matrix_GVM_genetic_M_V7 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 2, VID_value = 7 )
p_matrix_GVM_genetic_M_V8 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 2, VID_value = 8 )
p_matrix_GVM_genetic_M_V9 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                Gender_value = 2, VID_value = 9 )
p_matrix_GVM_genetic_M_V10 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                 Gender_value = 2, VID_value = 10 )
p_matrix_GVM_genetic_M_V11 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                 Gender_value = 2, VID_value = 11 )
p_matrix_GVM_genetic_M_V12 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                 Gender_value = 2, VID_value = 12 )
p_matrix_GVM_genetic_M_V13 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                 Gender_value = 2, VID_value = 13 )
p_matrix_GVM_genetic_M_V14 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_genetic , 
                                                                 Gender_value = 2, VID_value = 14 )

# Affinal kin density
p_matrix_GVM_affine_M_V1 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 2, VID_value = 1 )
p_matrix_GVM_affine_M_V2 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 2, VID_value = 2 )
p_matrix_GVM_affine_M_V3 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 2, VID_value = 3 )
p_matrix_GVM_affine_M_V4 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 2, VID_value = 4 )
p_matrix_GVM_affine_M_V5 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 2, VID_value = 5 )
p_matrix_GVM_affine_M_V6 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 2, VID_value = 6 )
p_matrix_GVM_affine_M_V7 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 2, VID_value = 7 )
p_matrix_GVM_affine_M_V8 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 2, VID_value = 8 )
p_matrix_GVM_affine_M_V9 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                               Gender_value = 2, VID_value = 9 )
p_matrix_GVM_affine_M_V10 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                                Gender_value = 2, VID_value = 10 )
p_matrix_GVM_affine_M_V11 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                                Gender_value = 2, VID_value = 11 )
p_matrix_GVM_affine_M_V12 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                                Gender_value = 2, VID_value = 12 )
p_matrix_GVM_affine_M_V13 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                                Gender_value = 2, VID_value = 13 )
p_matrix_GVM_affine_M_V14 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_affine , 
                                                                Gender_value = 2, VID_value = 14 )

# Friend density
p_matrix_GVM_fri_M_V1 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 2, VID_value = 1 )
p_matrix_GVM_fri_M_V2 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 2, VID_value = 2 )
p_matrix_GVM_fri_M_V3 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 2, VID_value = 3 )
p_matrix_GVM_fri_M_V4 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 2, VID_value = 4 )
p_matrix_GVM_fri_M_V5 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 2, VID_value = 5 )
p_matrix_GVM_fri_M_V6 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 2, VID_value = 6 )
p_matrix_GVM_fri_M_V7 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 2, VID_value = 7 )
p_matrix_GVM_fri_M_V8 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 2, VID_value = 8 )
p_matrix_GVM_fri_M_V9 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                            Gender_value = 2, VID_value = 9 )
p_matrix_GVM_fri_M_V10 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                             Gender_value = 2, VID_value = 10 )
p_matrix_GVM_fri_M_V11 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                             Gender_value = 2, VID_value = 11 )
p_matrix_GVM_fri_M_V12 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                             Gender_value = 2, VID_value = 12 )
p_matrix_GVM_fri_M_V13 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                             Gender_value = 2, VID_value = 13 )
p_matrix_GVM_fri_M_V14 <- generate_and_predict_p_matrix_GVM( post_samples = Post_model_GVM_fri , 
                                                             Gender_value = 2, VID_value = 14 )

#####3.3.5.2.1 Gender-specific effects of market integration stratified by villages ----
# Function to predict p-matrix and draw relevant figures
draw_plots_gender_specific <- function( plot_data_list, num_rows = 3, num_cols = 5) {
  # Set up the plotting space
  par(mfrow = c(num_rows, num_cols), mar = c(4.5, 4.5, 2.5, 1))
  
  for (plot_data in plot_data_list) {
    xtitle <- plot_data$xtitle
    ytitle <- plot_data$ytitle
    Text <- plot_data$Text
    x <- plot_data$x
    y <- plot_data$y
    Market_seq <- plot_data$Market_seq
    p_matrix_F <- plot_data$p_matrix_F
    p_matrix_M <- plot_data$p_matrix_M
    
    # Create the plot
    plot(NULL,
         xlim = c(-2.5, 2.2),
         ylim = c(0, 1),
         xlab = xtitle,
         ylab = ytitle,
         font.lab = 2,
         font.axis = 2,
         cex.axis = 1.5,
         cex.lab = 2)
    
    # Add lines and shaded areas for women
    lines(Market_seq,
          apply(p_matrix_F, 2, mean),
          lwd = 2, col = "#CB2313")
    shade(apply(p_matrix_F, 2, PI, prob = 0.90), Market_seq, col = col.alpha("#CB2313", 0.15))
    shade(apply(p_matrix_F, 2, PI, prob = 0.60), Market_seq, col = col.alpha("#CB2313", 0.15))
    shade(apply(p_matrix_F, 2, PI, prob = 0.30), Market_seq, col = col.alpha("#CB2313", 0.15))
    
    # Add lines and shaded areas for men
    lines(Market_seq,
          apply(p_matrix_M, 2, mean),
          lwd = 2, col = "#046C9A")
    shade(apply(p_matrix_M, 2, PI, prob = 0.90), Market_seq, col = col.alpha("#046C9A", 0.15))
    shade(apply(p_matrix_M, 2, PI, prob = 0.60), Market_seq, col = col.alpha("#046C9A", 0.15))
    shade(apply(p_matrix_M, 2, PI, prob = 0.30), Market_seq, col = col.alpha("#046C9A", 0.15))
    text( x , y , cex= 2 , font = 2 , Text )
  }
}

Plot_GVM_all_data_list <- list(
  # First row
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V1",
        x = 1.8 ,
        y = 0.05 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_all_F_V1 , 
        p_matrix_M = p_matrix_GVM_all_M_V1 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V2",
        x = 1.8 ,
        y = 0.05 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_all_F_V2 , 
        p_matrix_M = p_matrix_GVM_all_M_V2 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V3",
        x = 1.8 ,
        y = 0.05 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_all_F_V3 , 
        p_matrix_M = p_matrix_GVM_all_M_V3 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V4",
        x = 1.8 ,
        y = 0.05 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_all_F_V4 , 
        p_matrix_M = p_matrix_GVM_all_M_V4 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V5",
        x = 1.8 ,
        y = 0.05 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_all_F_V5 , 
        p_matrix_M = p_matrix_GVM_all_M_V5 ) ,
  
  # Second row
  list( xtitle = "" ,
        ytitle = "Predicted density" ,
        Text = "V6",
        x = 1.8 ,
        y = 0.05 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_all_F_V6 , 
        p_matrix_M = p_matrix_GVM_all_M_V6 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V7",
        x = 1.8 ,
        y = 0.05 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_all_F_V7 , 
        p_matrix_M = p_matrix_GVM_all_M_V7 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V8",
        x = 1.8 ,
        y = 0.05 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_all_F_V8 , 
        p_matrix_M = p_matrix_GVM_all_M_V8 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V9",
        x = 1.8 ,
        y = 0.05 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_all_F_V9 , 
        p_matrix_M = p_matrix_GVM_all_M_V9 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V10",
        x = 1.8 ,
        y = 0.05 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_all_F_V10 , 
        p_matrix_M = p_matrix_GVM_all_M_V10 ) ,
  
  # Third row
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V11",
        x = 1.8 ,
        y = 0.05 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_all_F_V11 , 
        p_matrix_M = p_matrix_GVM_all_M_V11 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V12",
        x = 1.8 ,
        y = 0.05 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_all_F_V12 , 
        p_matrix_M = p_matrix_GVM_all_M_V12 ) ,
  list( xtitle = "Market integration (std)" ,
        ytitle = "" ,
        Text = "V13",
        x = 1.8 ,
        y = 0.05 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_all_F_V13 , 
        p_matrix_M = p_matrix_GVM_all_M_V13 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V14",
        x = 1.8 ,
        y = 0.05 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_all_F_V14 , 
        p_matrix_M = p_matrix_GVM_all_M_V14 ) 
)

jpeg( "Market_gender_village_all.jpg" , 
      width = 365 , height = 230 , units = "mm" , res = 300 )
draw_plots_gender_specific( Plot_GVM_all_data_list )
mtext( "All" , side = 3 , line = -2 , outer = TRUE , cex = 1.5 , font = 2 )
dev.off( )

Plot_GVM_genetic_data_list <- list(
  # First row
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V1",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_genetic_F_V1 , 
        p_matrix_M = p_matrix_GVM_genetic_M_V1 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V2",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_genetic_F_V2 , 
        p_matrix_M = p_matrix_GVM_genetic_M_V2 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V3",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_genetic_F_V3 , 
        p_matrix_M = p_matrix_GVM_genetic_M_V3 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V4",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_genetic_F_V4 , 
        p_matrix_M = p_matrix_GVM_genetic_M_V4 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V5",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_genetic_F_V5 , 
        p_matrix_M = p_matrix_GVM_genetic_M_V5 ) ,
  
  # Second row
  list( xtitle = "" ,
        ytitle = "Predicted density" ,
        Text = "V6",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_genetic_F_V6 , 
        p_matrix_M = p_matrix_GVM_genetic_M_V6 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V7",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_genetic_F_V7 , 
        p_matrix_M = p_matrix_GVM_genetic_M_V7 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V8",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_genetic_F_V8 , 
        p_matrix_M = p_matrix_GVM_genetic_M_V8 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V9",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_genetic_F_V9 , 
        p_matrix_M = p_matrix_GVM_genetic_M_V9 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V10",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_genetic_F_V10 , 
        p_matrix_M = p_matrix_GVM_genetic_M_V10 ) ,
  
  # Third row
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V11",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_genetic_F_V11 , 
        p_matrix_M = p_matrix_GVM_genetic_M_V11 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V12",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_genetic_F_V12 , 
        p_matrix_M = p_matrix_GVM_genetic_M_V12 ) ,
  list( xtitle = "Market integration (std)" ,
        ytitle = "" ,
        Text = "V13",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_genetic_F_V13 , 
        p_matrix_M = p_matrix_GVM_genetic_M_V13 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V14",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_genetic_F_V14 , 
        p_matrix_M = p_matrix_GVM_genetic_M_V14 ) 
)

jpeg( "Market_gender_village_genetic.jpg" , 
      width = 365 , height = 230 , units = "mm" , res = 300 )
draw_plots_gender_specific( Plot_GVM_genetic_data_list )
mtext( "Biological kin" , side = 3 , line = -2 , outer = TRUE , cex = 1.5 , font = 2 )
dev.off( )

Plot_GVM_affine_data_list <- list(
  # First row
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V1",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_affine_F_V1 , 
        p_matrix_M = p_matrix_GVM_affine_M_V1 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V2",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_affine_F_V2 , 
        p_matrix_M = p_matrix_GVM_affine_M_V2 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V3",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_affine_F_V3 , 
        p_matrix_M = p_matrix_GVM_affine_M_V3 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V4",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_affine_F_V4 , 
        p_matrix_M = p_matrix_GVM_affine_M_V4 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V5",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_affine_F_V5 , 
        p_matrix_M = p_matrix_GVM_affine_M_V5 ) ,
  
  # Second row
  list( xtitle = "" ,
        ytitle = "Predicted density" ,
        Text = "V6",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_affine_F_V6 , 
        p_matrix_M = p_matrix_GVM_affine_M_V6 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V7",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_affine_F_V7 , 
        p_matrix_M = p_matrix_GVM_affine_M_V7 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V8",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_affine_F_V8 , 
        p_matrix_M = p_matrix_GVM_affine_M_V8 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V9",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_affine_F_V9 , 
        p_matrix_M = p_matrix_GVM_affine_M_V9 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V10",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_affine_F_V10 , 
        p_matrix_M = p_matrix_GVM_affine_M_V10 ) ,
  
  # Third row
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V11",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_affine_F_V11 , 
        p_matrix_M = p_matrix_GVM_affine_M_V11 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V12",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_affine_F_V12 , 
        p_matrix_M = p_matrix_GVM_affine_M_V12 ) ,
  list( xtitle = "Market integration (std)" ,
        ytitle = "" ,
        Text = "V13",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_affine_F_V13 , 
        p_matrix_M = p_matrix_GVM_affine_M_V13 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V14",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_affine_F_V14 , 
        p_matrix_M = p_matrix_GVM_affine_M_V14 ) 
)

jpeg( "Market_gender_village_affinal.jpg" , 
      width = 365 , height = 230 , units = "mm" , res = 300 )
draw_plots_gender_specific( Plot_GVM_affine_data_list )
mtext( "Affinal kin" , side = 3 , line = -2 , outer = TRUE , cex = 1.5 , font = 2 )
dev.off( )

Plot_GVM_fri_data_list <- list(
  # First row
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V1",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_fri_F_V1 , 
        p_matrix_M = p_matrix_GVM_fri_M_V1 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V2",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_fri_F_V2 , 
        p_matrix_M = p_matrix_GVM_fri_M_V2 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V3",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_fri_F_V3 , 
        p_matrix_M = p_matrix_GVM_fri_M_V3 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V4",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_fri_F_V4 , 
        p_matrix_M = p_matrix_GVM_fri_M_V4 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V5",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_fri_F_V5 , 
        p_matrix_M = p_matrix_GVM_fri_M_V5 ) ,
  
  # Second row
  list( xtitle = "" ,
        ytitle = "Predicted density" ,
        Text = "V6",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_fri_F_V6 , 
        p_matrix_M = p_matrix_GVM_fri_M_V6 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V7",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_fri_F_V7 , 
        p_matrix_M = p_matrix_GVM_fri_M_V7 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V8",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_fri_F_V8 , 
        p_matrix_M = p_matrix_GVM_fri_M_V8 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V9",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_fri_F_V9 , 
        p_matrix_M = p_matrix_GVM_fri_M_V9 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V10",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_fri_F_V10 , 
        p_matrix_M = p_matrix_GVM_fri_M_V10 ) ,
  
  # Third row
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V11",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_fri_F_V11 , 
        p_matrix_M = p_matrix_GVM_fri_M_V11 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V12",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_fri_F_V12 , 
        p_matrix_M = p_matrix_GVM_fri_M_V12 ) ,
  list( xtitle = "Market integration (std)" ,
        ytitle = "" ,
        Text = "V13",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_fri_F_V13 , 
        p_matrix_M = p_matrix_GVM_fri_M_V13 ) ,
  list( xtitle = "" ,
        ytitle = "" ,
        Text = "V14",
        x = 1.8 ,
        y = 0.95 ,
        Market_seq = Market_seq , 
        p_matrix_F = p_matrix_GVM_fri_F_V14 , 
        p_matrix_M = p_matrix_GVM_fri_M_V14 ) 
)

jpeg( "Market_gender_village_fri.jpg" , 
      width = 365 , height = 230 , units = "mm" , res = 300 )
draw_plots_gender_specific( Plot_GVM_fri_data_list )
mtext( "Friend" , side = 3 , line = -2 , outer = TRUE , cex = 1.5 , font = 2 )
dev.off( )

#####3.3.5.2.2 Gender difference in effects of market integration ----
# Function to predict p-matrix and draw relevant figures
draw_plots_gender_dif <- function( plot_data_list, num_rows = 3, num_cols = 5) {
  # Set up the plotting space
  par(mfrow = c(num_rows, num_cols), mar = c(4.5, 4.5, 2.5, 1))
  
  for (plot_data in plot_data_list) {
    xtitle <- plot_data$xtitle
    ytitle <- plot_data$ytitle
    Text <- plot_data$Text
    Market_seq <- plot_data$Market_seq
    p_matrix_F <- plot_data$p_matrix_F
    p_matrix_M <- plot_data$p_matrix_M
    
    # Create the plot
    plot(NULL,
         xlim = c(-2.5, 2.2),
         ylim = c( -1 , 1 ),
         xlab = xtitle,
         ylab = ytitle,
         font.lab = 2,
         font.axis = 2,
         cex.axis = 1.5,
         cex.lab = 2)
    
    # Add lines and shaded areas for women
    lines(Market_seq,
          apply(p_matrix_F - p_matrix_M, 2, mean),
          lwd = 2, col = "black")
    shade(apply(p_matrix_F - p_matrix_M, 2, PI, prob = 0.90), Market_seq, col = col.alpha("dimgray", 0.15))
    shade(apply(p_matrix_F - p_matrix_M, 2, PI, prob = 0.60), Market_seq, col = col.alpha("dimgray", 0.15))
    shade(apply(p_matrix_F - p_matrix_M, 2, PI, prob = 0.30), Market_seq, col = col.alpha("dimgray", 0.15))
    abline( h = 0 , lty = 2 , lwd = 2, col = 1 )
    text( 1.9 , -0.9 , cex= 2 , font = 2 , Text )
  }
}

jpeg( "Market_village_gender_diff_all.jpg" , 
      width = 365 , height = 230 , units = "mm" , res = 300 )
draw_plots_gender_dif( Plot_GVM_all_data_list )
mtext( "All" , side = 3 , line = -2 , outer = TRUE , cex = 1.5 , font = 2 )
dev.off( )

jpeg( "Market_village_gender_diff_genetic.jpg" , 
      width = 365 , height = 230 , units = "mm" , res = 300 )
draw_plots_gender_dif( Plot_GVM_genetic_data_list )
mtext( "Biological kin" , side = 3 , line = -2 , outer = TRUE , cex = 1.5 , font = 2 )
dev.off( )

jpeg( "Market_village_gender_diff_affinal.jpg" , 
      width = 365 , height = 230 , units = "mm" , res = 300 )
draw_plots_gender_dif( Plot_GVM_affine_data_list )
mtext( "Affinal kin" , side = 3 , line = -2 , outer = TRUE , cex = 1.5 , font = 2 )
dev.off( )

jpeg( "Market_village_gender_diff_fri.jpg" , 
      width = 365 , height = 230 , units = "mm" , res = 300 )
draw_plots_gender_dif( Plot_GVM_fri_data_list )
mtext( "Friend" , side = 3 , line = -2 , outer = TRUE , cex = 1.5 , font = 2 )
dev.off( )

###3.3.6 Outputs ----
Output_model_GM <- data.frame( Mean = c( Pre_model_GM_all$mean ,
                                         Pre_model_GM_genetic$mean ,
                                         Pre_model_GM_affine$mean ,
                                         Pre_model_GM_fri$mean ) ,
                               CI5 = c( Pre_model_GM_all$`5%` ,
                                        Pre_model_GM_genetic$`5%` ,
                                        Pre_model_GM_affine$`5%` ,
                                        Pre_model_GM_fri$`5%` ) ,
                               CI95 = c( Pre_model_GM_all$`95%` ,
                                         Pre_model_GM_genetic$`95%` ,
                                         Pre_model_GM_affine$`95%` ,
                                         Pre_model_GM_fri$`95%` ) ,
                               Type = rep( c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) , each = 12 ) ,
                               Variable = rep( c( "Age[female]" , "Age[male]" , 
                                                  "Age^2[female]" , "Age^2[male]" , 
                                                  "Female: matrilocal" , 
                                                  "Male: matrilocal" , 
                                                  "Female: bilocal" , 
                                                  "Male: bilocal" , 
                                                  "Female: patrilocal" , 
                                                  "Male: patrilocal" , 
                                                  "Market[female]" , "Market[male]" ) , 4 ) ) %>% 
  mutate( Type = factor( Type , 
                         levels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ,
                         labels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) ,
          Variable = factor( Variable , 
                             levels = c( "Age[female]" , "Age[male]" , 
                                         "Age^2[female]" , "Age^2[male]" ,
                                         "Female: matrilocal" , 
                                         "Female: bilocal" , 
                                         "Female: patrilocal" , 
                                         "Male: patrilocal" , 
                                         "Male: bilocal" , 
                                         "Male: matrilocal" , 
                                         "Market[female]" , "Market[male]" ) ,
                             labels = c( "Age[female]" , "Age[male]" , 
                                         "Age^2[female]" , "Age^2[male]" ,
                                         "Female: matrilocal" , 
                                         "Female: bilocal" , 
                                         "Female: patrilocal" , 
                                         "Male: patrilocal" , 
                                         "Male: bilocal" , 
                                         "Male: matrilocal" , 
                                         "Market[female]" , "Market[male]" ) ) ,
          Model = "Market integration" , 
          "Mean [CI]" = paste( format( round( Mean , 2 ), nsmall = 2 ) , 
                               " [" , 
                               format( round( CI5 , 2 ), nsmall = 2 ) , 
                               ", " , 
                               format( round( CI95 , 2 ), nsmall = 2 ) , 
                               "]", 
                               sep = "") )

Output_model_GVM <- data.frame( Mean = c( Pre_model_GVM_all$mean ,
                                          Pre_model_GVM_genetic$mean ,
                                          Pre_model_GVM_affine$mean ,
                                          Pre_model_GVM_fri$mean ) ,
                                CI5 = c( Pre_model_GVM_all$`5%` ,
                                         Pre_model_GVM_genetic$`5%` ,
                                         Pre_model_GVM_affine$`5%` ,
                                         Pre_model_GVM_fri$`5%` ) ,
                                CI95 = c( Pre_model_GVM_all$`95%` ,
                                          Pre_model_GVM_genetic$`95%` ,
                                          Pre_model_GVM_affine$`95%` ,
                                          Pre_model_GVM_fri$`95%` ) ,
                                Type = rep( c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) , each = 40 ) ,
                                Variable = rep( c( "Age[female]" , "Age[male]" , 
                                                   "Age^2[female]" , "Age^2[male]" , 
                                                   "Female: matrilocal" , 
                                                   "Male: matrilocal" , 
                                                   "Female: bilocal" , 
                                                   "Male: bilocal" , 
                                                   "Female: patrilocal" , 
                                                   "Male: patrilocal" , 
                                                   "V1_Market[female]" , "V1_Market[male]" , 
                                                   "V2_Market[female]" , "V2_Market[male]" , 
                                                   "V3_Market[female]" , "V3_Market[male]" , 
                                                   "V4_Market[female]" , "V4_Market[male]" , 
                                                   "V5_Market[female]" , "V5_Market[male]" , 
                                                   "V6_Market[female]" , "V6_Market[male]" , 
                                                   "V7_Market[female]" , "V7_Market[male]" , 
                                                   "V8_Market[female]" , "V8_Market[male]" , 
                                                   "V9_Market[female]" , "V9_Market[male]" , 
                                                   "V10_Market[female]" , "V10_Market[male]" , 
                                                   "V11_Market[female]" , "V11_Market[male]" , 
                                                   "V12_Market[female]" , "V12_Market[male]" , 
                                                   "V13_Market[female]" , "V13_Market[male]" , 
                                                   "V14_Market[female]" , "V14_Market[male]" , 
                                                   "Sigma_V" , "Sigma_I" ) , 4 ) ) %>% 
  mutate( Type = factor( Type , 
                         levels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ,
                         labels = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) ,
          Variable = factor( Variable , 
                             levels = c( "Age[female]" , "Age[male]" , 
                                         "Age^2[female]" , "Age^2[male]" ,
                                         "Female: matrilocal" , 
                                         "Female: bilocal" , 
                                         "Female: patrilocal" , 
                                         "Male: matrilocal" , 
                                         "Male: bilocal" , 
                                         "Male: patrilocal" , 
                                         "V1_Market[female]" , "V1_Market[male]" , 
                                         "V2_Market[female]" , "V2_Market[male]" , 
                                         "V3_Market[female]" , "V3_Market[male]" , 
                                         "V4_Market[female]" , "V4_Market[male]" , 
                                         "V5_Market[female]" , "V5_Market[male]" , 
                                         "V6_Market[female]" , "V6_Market[male]" , 
                                         "V7_Market[female]" , "V7_Market[male]" , 
                                         "V8_Market[female]" , "V8_Market[male]" , 
                                         "V9_Market[female]" , "V9_Market[male]" , 
                                         "V10_Market[female]" , "V10_Market[male]" , 
                                         "V11_Market[female]" , "V11_Market[male]" , 
                                         "V12_Market[female]" , "V12_Market[male]" , 
                                         "V13_Market[female]" , "V13_Market[male]" , 
                                         "V14_Market[female]" , "V14_Market[male]" , 
                                         "Sigma_V" , "Sigma_I"  ) ,
                             labels = c( "Age[female]" , "Age[male]" , 
                                         "Age^2[female]" , "Age^2[male]" ,
                                         "Female: matrilocal" , 
                                         "Female: bilocal" , 
                                         "Female: patrilocal" , 
                                         "Male: matrilocal" , 
                                         "Male: bilocal" , 
                                         "Male: patrilocal", 
                                         "V1_Market[female]" , "V1_Market[male]" , 
                                         "V2_Market[female]" , "V2_Market[male]" , 
                                         "V3_Market[female]" , "V3_Market[male]" , 
                                         "V4_Market[female]" , "V4_Market[male]" , 
                                         "V5_Market[female]" , "V5_Market[male]" , 
                                         "V6_Market[female]" , "V6_Market[male]" , 
                                         "V7_Market[female]" , "V7_Market[male]" , 
                                         "V8_Market[female]" , "V8_Market[male]" , 
                                         "V9_Market[female]" , "V9_Market[male]" , 
                                         "V10_Market[female]" , "V10_Market[male]" , 
                                         "V11_Market[female]" , "V11_Market[male]" , 
                                         "V12_Market[female]" , "V12_Market[male]" , 
                                         "V13_Market[female]" , "V13_Market[male]" , 
                                         "V14_Market[female]" , "V14_Market[male]" , 
                                         "Sigma_V" , "Sigma_I"  ) ) ,
          Model = "Market integration" , 
          "Mean [CI]" = paste( format( round( Mean , 2 ), nsmall = 2 ) , 
                               " [" , 
                               format( round( CI5 , 2 ), nsmall = 2 ) , 
                               ", " , 
                               format( round( CI95 , 2 ), nsmall = 2 ) , 
                               "]", 
                               sep = "") )

###3.3.6.1 Market integration stratified by gender ----
GM.output <- tibble( Variable =  c( "Age[female]" , "Age[male]" , 
                                    "Age^2[female]" , "Age^2[male]" ,
                                    "Female: matrilocal" , 
                                    "Female: bilocal" , 
                                    "Female: patrilocal" , 
                                    "Male: matrilocal" , 
                                    "Male: bilocal" , 
                                    "Male: patrilocal" ,
                                    "Market[female]" , "Market[male]" ) ,
                     Type = "All" ) %>% 
  left_join( Output_model_GM[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `All` = "Mean [CI]" ) %>% 
  mutate( Type = "Biological kin" ) %>% 
  left_join( Output_model_GM[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `Biological kin` = "Mean [CI]" ) %>% 
  mutate( Type = "Affinal kin" ) %>% 
  left_join( Output_model_GM[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `Affinal kin` = "Mean [CI]" ) %>% 
  mutate( Type = "Friend" ) %>% 
  left_join( Output_model_GM[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( `Friend` = "Mean [CI]" ) %>% 
  select( - Type ) %>% 
  gt() %>% 
  fmt_number( ) %>% 
  sub_missing( columns = everything() , 
               missing_text = " " ) %>% 
  tab_style( style = list( cell_text( weight = "bold" ) ) ,
             locations = cells_column_labels( everything() ) ) %>% 
  tab_style( style = list( cell_text( font = "Times New Roman" ) ) ,
             locations = cells_body() ) %>% 
  cols_align( align = "left",
              columns = "Variable" ) %>%
  cols_align( align = "center",
              columns = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) %>% 
  opt_table_lines( extent = "default" ) %>%
  tab_options( column_labels.border.top.color = "black" ,
               column_labels.border.top.width = px( 3 ) ,
               column_labels.border.bottom.color = "black" ,
               table_body.hlines.color = "lightgrey" ,
               table.border.bottom.color = "black" ,
               table.border.bottom.width = px( 3 ) ) %>% 
  cols_width( "Variable" ~ px( 170 ) , 
              everything() ~ px( 150 ) )
GM.output

###3.3.6.2 Market integration stratified by gender and village ----
GVM.output <- tibble( Variable =  c( "Age[female]" , "Age[male]" , 
                                     "Age^2[female]" , "Age^2[male]" ,
                                     "Female: matrilocal" , 
                                     "Female: bilocal" , 
                                     "Female: patrilocal" , 
                                     "Male: matrilocal" , 
                                     "Male: bilocal" , 
                                     "Male: patrilocal" ,
                                     "V1_Market[female]" , "V1_Market[male]" , 
                                     "V2_Market[female]" , "V2_Market[male]" , 
                                     "V3_Market[female]" , "V3_Market[male]" , 
                                     "V4_Market[female]" , "V4_Market[male]" , 
                                     "V5_Market[female]" , "V5_Market[male]" , 
                                     "V6_Market[female]" , "V6_Market[male]" , 
                                     "V7_Market[female]" , "V7_Market[male]" , 
                                     "V8_Market[female]" , "V8_Market[male]" , 
                                     "V9_Market[female]" , "V9_Market[male]" , 
                                     "V10_Market[female]" , "V10_Market[male]" , 
                                     "V11_Market[female]" , "V11_Market[male]" , 
                                     "V12_Market[female]" , "V12_Market[male]" , 
                                     "V13_Market[female]" , "V13_Market[male]" , 
                                     "V14_Market[female]" , "V14_Market[male]" ,
                                     "Sigma_V" , "Sigma_I" ) ,
                      Type = "All" ) %>% 
  left_join( Output_model_GVM[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( All = "Mean [CI]" ) %>% 
  mutate( Type = "Biological kin" ) %>% 
  left_join( Output_model_GVM[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( "Biological kin" = "Mean [CI]" ) %>% 
  mutate( Type = "Affinal kin" ) %>% 
  left_join( Output_model_GVM[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( "Affinal kin" = "Mean [CI]" ) %>% 
  mutate( Type = "Friend" ) %>% 
  left_join( Output_model_GVM[ , c( 4 , 5 , 7 ) ] ,
             by = c( "Variable" , "Type" ) ) %>% 
  dplyr::rename( Friend = "Mean [CI]" ) %>% 
  select( - Type ) %>% 
  gt() %>% 
  fmt_number( ) %>% 
  tab_style( style = list( cell_text( weight = "bold" ) ) ,
             locations = cells_column_labels( everything() ) ) %>% 
  tab_style( style = list( cell_text( font = "Times New Roman" ) ) ,
             locations = cells_body() ) %>% 
  cols_align( align = "left",
              columns = "Variable" ) %>%
  cols_align( align = "center",
              columns = c( "All" , "Biological kin" , "Affinal kin" , "Friend" ) ) %>% 
  opt_table_lines( extent = "default" ) %>%
  tab_options( column_labels.border.top.color = "black" ,
               column_labels.border.top.width = px( 3 ) ,
               column_labels.border.bottom.color = "black" ,
               table_body.hlines.color = "lightgrey" ,
               table.border.bottom.color = "black" ,
               table.border.bottom.width = px( 3 ) ) %>% 
  cols_width( "Variable" ~ px( 170 ) , 
              everything() ~ px( 150 ) )
GVM.output

