


# create matrix of DGP settings
settings <- as.data.frame(matrix(c(500,500,500,500,500,500,50,50,50,20,20,20), ncol=2))
settings$V3 <- c(F,F,F,F,T,T)
settings$V4 <- c("constant","con_lin","con_non","binary","con_non","binary")
settings$V5 <- c(0.5,NA,NA,NA,NA,NA)



S = 10
M <- 50
ntile <- 5

error_matrix <-  matrix(NA,S,3)
colnames(error_matrix) <- c("GATES_MAE","DO_GATES_MAE", "True_Treatment")


error_result <- list()
GATES_result <- list()


for(t in 1:nrow(settings)) {
  
  N <- settings[t,1]
  k <- settings[t,2] 
  ID <- c(1:N)
  
  prop <- matrix(NA,N,M)
  prop <- cbind(prop,ID)
  
  pred_tm_gates <- matrix(NA,N,M)
  pred_dr_gates <- matrix(NA,N,M)
  pred_doubleML_gates <- matrix(NA,N,M)
  
  pred_tm_gates <- cbind(pred_tm_gates,ID)
  pred_dr_gates <- cbind(pred_dr_gates,ID)
  pred_doubleML_gates <- cbind(pred_doubleML_gates,ID)
  
  
  
  list_res_TM <- vector("list", M)
  list_res_DR <- vector("list", M)
  list_res_DoubleML <- vector("list", M)
  
  
  for(j in 1:S){ 
    
    theta_set <- ifelse(settings[t,4]=="constant",settings[t,5],settings[t,4])
    
    dataset <- datagen(y="con", N=settings[t,1],k=settings[t,2],random_d=settings[t,3],theta=theta_set,var=1)
    dataset$ID <- c(1:N)
    
    k <- ncol(dataset)-4
    covariates <- c(paste0("V", 1:k)) 
    covariates
    covariates_d <- c(paste0("V", 1:k),"d")
    
    dataset$d <- as.factor(ifelse(dataset$d==1,1,0))
    
    
    
    for(i in 1:M){
      
      
      
      
      ##### Parameter and datasets ##### 
      
      trainIndex <- createDataPartition(dataset$d, p = .5, list = FALSE) 
      df_aux <- dataset[trainIndex,] 
      df_main  <- dataset[-trainIndex,]
      
      
      
      # On the auxiliary sample
      # -----------------------
      
      # Propensity score using regression forests
      
      
      
      rf_prop <- ranger(d~.,data=df_aux[covariates_d],probability = T, importance= "impurity")
      p_dr <- predict(rf_prop,data=df_aux[,covariates_d])$predictions[,2]
      p <- predict(rf_prop,data=df_main[,covariates_d])$predictions[,2]
      
      # Conditional mean proxy using regression forests
      
      aux_1 <- df_aux[which(df_aux$d==1),]
      aux_0 <- df_aux[which(df_aux$d==0),]
      
      
      form  <-  as.formula(paste("y", paste(covariates, collapse=" + "), sep=" ~ "))
      
      rf_1 <- ranger(form,data=aux_1)
      rf_0 <- ranger(form,data=aux_0)
      
      y1_dr <- predict(rf_1,df_aux)$predictions
      y0_dr <- predict(rf_0,df_aux)$predictions
      
      y1 <- predict(rf_1,df_main)$predictions
      y0 <- predict(rf_0,df_main)$predictions
      
      # On the main sample
      # -----------------------
      
      # Propensity score offset W - e(X)
      
      df_main$d <- as.numeric(as.character(df_main$d)) - p
      
      
      ind_1  <- (p_dr>0.02 & p_dr<0.98)
      ind <- (p>0.02 & p<0.98)
      
      
      y1_dr <- y1_dr[ind_1]
      y0_dr <- y0_dr[ind_1]
      p_dr <- p_dr[ind_1]
      df_aux <- df_aux[ind_1,]
      
      p <- p[ind]
      y1 <- y1[ind]
      y0 <- y0[ind]
      df_main <- df_main[ind,]
      
      
      
      # Score function for distribution
      df_main$S <- (y1 - y0) 
      prop[,i][df_main$ID] <- df_main$S
      
      
      # Divide observations into k-tiles
      
      S2 <- df_main$S +runif(length(df_main$S), 0, 0.00001) # Include white noise to guarantee that the score (S) differs from the baseline effect
     
      breaks    <- quantile(S2, seq(0,1, 0.2),  include.lowest =T)
      breaks[1] <- breaks[1] - 0.01 # Offset for lower tails 
      breaks[6] <- breaks[6] + 0.01 # Offset for upper tails
      
      SG <- cut(S2, breaks = breaks)
      
      
      
      SGX       <- model.matrix(~-1+SG) # -1 Ereases the Intercept. Possible is also to keep the Intercept. 
      DSG       <- data.frame(as.numeric(I(as.numeric(df_main[,"d"])))*SGX)
      
      colnames(DSG) <- c("G1", "G2", "G3", "G4", "G5")
      df_main[,c("G1", "G2", "G3", "G4", "G5", "weight")] <- cbind( DSG$G1, DSG$G2, DSG$G3, DSG$G4, DSG$G5, as.numeric((1/(p*(1-p)))))
      
      
      form1 <- as.formula(paste("y", "~", "G1+G2+G3+G4+G5 ", sep=""))
      df_main$y <- as.numeric(df_main$y)
      
      
      
      # Now regress on group membership variables
      
      model <- lm(form1,df_main, weights = df_main$weight)
      groups <- c(paste0("G",1:ntile))
      groups <- dput(as.character(groups))
      
      thetahat1 <- model%>% 
        .$coefficients %>%
        .[groups]
      
      ####
      gates_zero_help <- df_main[colnames(DSG)]
      
      gates_zero <- as.data.frame(which(gates_zero_help!=0,arr.ind = T))
      gates_zero[,c("ID")] <- rownames(gates_zero)
      gates_zero <- gates_zero[,-1]
      
      
      
      thetahat2 <- as.data.frame(thetahat1)
      rownames(thetahat2) <- c("1","2","3","4","5")
      thetahat2["col"] <- rownames(thetahat2)
      head(thetahat2)
      
      gates_y <- merge(thetahat2,gates_zero,"col")
      gates_y$ID <- as.integer(gates_y$ID)
      pred_tm_gates[,i][gates_y$ID] <- gates_y$thetahat1
      ####
      
      
      
      # Confidence intervals
      cihat <- confint(model,level=0.9)[groups,]
      list_res_TM[[i]] <- tibble(coefficient = dput(as.character(c(paste0("Group", 1:ntile)))),
                                 estimates = thetahat1,
                                 ci_lower_90 = cihat[,1],
                                 ci_upper_90 = cihat[,2])
      
      
      
      
      
      
      
      
      
      
      #### This part is Doubly-Robust ####################
      
      
      
      
      
      
      
      # Doubly Robust
      df_aux$d <- as.numeric(ifelse(df_aux$d==1,1,0))
      
      y_mo <- (y1_dr - y0_dr) + ((df_aux$d*(df_aux$y-y1_dr))/p_dr) - ((1-df_aux$d)*(df_aux$y-y0_dr)/(1-p_dr))
      
      
      rf_dr <- ranger(y_mo~.,data=df_aux[covariates], importance = "impurity")
      score_dr <- predict(rf_dr,data=df_main[covariates])$predictions
      
      
      
      
      # Divide observations into k-tiles
      df_main$S <- (score_dr) 
      S2 <- df_main$S +runif(length(df_main$S), 0, 0.00001) 
      
      SG <- cut(S2, breaks = ntile)
      
      
      ## Double Orthogonal Scores - using u_hat ###################
      # Predict conditional mean of Y without D
      
      form_mu  <-  as.formula(paste("y", paste(covariates, collapse=" + "), sep=" ~ "))
      rf_mu <- ranger(form_mu,data=df_aux)
      y_hat <- predict(rf_mu,data=df_main)$predictions
      
      
      df_main$u_hat <- df_main$y - y_hat
      
      
      SGX       <- model.matrix(~-1+SG) # -1 Ereases the Intercept. Possible is also to keep the Intercept. 
      DSG       <- data.frame(as.numeric(I(as.numeric(df_main[,"d"])))*SGX)
      
      colnames(DSG) <- c("G1", "G2", "G3", "G4", "G5")
      df_main[,c("S", "G1", "G2", "G3", "G4", "G5", "weight")] <- cbind(df_main$S, DSG$G1, DSG$G2, DSG$G3, DSG$G4, DSG$G5, as.numeric((1/(p*(1-p)))))
      
      
      form1 <- as.formula(paste("u_hat", "~", "G1+G2+G3+G4+G5 ", sep=""))
      
      
      
      
      # Now regress on group membership variables
  
      model <- lm(form1,df_main)
      groups <- c(paste0("G",1:ntile))
      groups <- dput(as.character(groups))
      
      thetahat1 <- model%>% 
        .$coefficients %>%
        .[groups]
      
      # Confidence intervals
      cihat <- confint(model,level=0.9)[groups,]
      list_res_DoubleML[[i]] <- tibble(coefficient = dput(as.character(c(paste0("Group", 1:ntile)))),
                                       estimates = thetahat1,
                                       ci_lower_90 = cihat[,1],
                                       ci_upper_90 = cihat[,2])
      
      
      
      
      ####
      gates_zero_help <- df_main[colnames(DSG)]
      
      gates_zero <- as.data.frame(which(gates_zero_help!=0,arr.ind = T))
      gates_zero[,c("ID")] <- rownames(gates_zero)
      gates_zero <- gates_zero[,-1]
      
      
      
      thetahat2 <- as.data.frame(thetahat1)
      rownames(thetahat2) <- c("1","2","3","4","5")
      thetahat2["col"] <- rownames(thetahat2)
      head(thetahat2)
      
      gates_y <- merge(thetahat2,gates_zero,"col")
      gates_y$ID <- as.integer(gates_y$ID)
      pred_doubleML_gates[,i][gates_y$ID] <- gates_y$thetahat1
      
  
      

      
      
    }
    
    
    GATES_TM <- list_res_TM[] %>% 
      bind_rows %>%
      na.omit() %>%
      group_by(coefficient) %>%
      summarize_all(median) 
    
    
    pred_tm_gates_median <- pred_tm_gates[,-ncol(pred_tm_gates)]
    apply(pred_tm_gates_median,1, median, na.rm = TRUE) # Calculate the row median which is then used to classify each obs. into a "group". 
    
    error_matrix[j,1] <- mean(abs(dataset$theta-pred_tm_gates_median),na.rm=T)
    
    
    GATES_DoubleML <- list_res_DoubleML[] %>% 
      bind_rows %>%
      group_by(coefficient) %>%
      summarize_all(median)
    
    pred_doubleML_gates_median <- pred_doubleML_gates[,-ncol(pred_doubleML_gates)]
    apply(pred_doubleML_gates_median,1, median, na.rm = TRUE) # Calculate the row median which is then used to classify each obs. into a "group". 
    
    error_matrix[j,2] <- mean(abs(dataset$theta-pred_doubleML_gates_median),na.rm=T)
    
    
    
    
   
    error_matrix[j,3] <- mean (dataset$theta)
    
    
    print(paste0("................................... ","The current iteration is: ", j, " out of " ,S))
  }
  
  error_result[[t]] <- error_matrix
  GATES_result[[t]] <- c(GATES_TM, GATES_DR)
  
  print(paste0("................................... ","This is DGP : ", t, " out of " ,nrow(settings)))
  
}


##########################



error_result


error_all <- matrix(NA,S*nrow(settings),3)
error_all


b=0

for(j in 1:nrow(settings)){
  for(i in 1:S){
    
    error_all[i+b,1] <- error_result[[j]][i]
    error_all[i+b,2] <- error_result[[j]][i,2]
    error_all[i+b,3] <- as.numeric(j)
  }
  b = b+S 
  
}

error_all <- as.data.frame(error_all)
colnames(error_all) <- c("GATES_MAE", "DO_GATES_MAE", "SETTING")


# wilcox-test for mean differences between groups
wilcox.test(error_matrix[1:10,1],error_matrix[1:10,2]) # Better use non-parametric test since the assumption that X and Y are ~ N(.) is not fulfilled. 

mean_error_1 <- c()
mean_error_2 <- c()
j = 1
for(i in 1:6){
  
  mean_error_1[i] <- mean(safe_error_all_newDGP_lowDim[j:j+9,1])
  mean_error_2[i] <- mean(safe_error_all_newDGP_lowDim[j:j+9,3])
  j = j +10 
}

round(mean_error_1,2)
round(mean_error_2,2)



ggplot(error_all, aes(x=GATES_MAE, y=DO_GATES_MAE)) + 
  xlim(0.0,1.0) +
  ylim(0.0,1.0) +
  geom_abline(mapping= aes(intercept=0.0,slope = 1.0, color="45 Degree line")) +
  scale_colour_manual(values="red") +
  labs(colour="") +
  geom_point() + 
  theme_cowplot() +
  facet_wrap( ~ SETTING, scales="free", ncol=3)  + # Facet wrap with common scales 
  guides(fill = FALSE, color = FALSE, linetype = FALSE, shape = FALSE) +
  labs(x = "TMA GATES", y = "DO GATES")



