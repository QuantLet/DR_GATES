


# create matrix of DGP settings
settings <- as.data.frame(matrix(c(1000,1000,2000,2000,2000,2000,10,20,10,20,20,20), ncol=2))
settings$V3 <- c(F,F,F,F,T,T)
settings$V4 <- c("constant","con","con","binary","con","binary")
settings$V5 <- c(0.5,NA,NA,NA,NA,NA)



S = 2
M <- 10
ntile <- 5

error_matrix <-  matrix(NA,S,3)
colnames(error_matrix) <- c("GATES_MAE", "DR_MAE", "True_Treatment")


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
pred_tm_gates <- cbind(pred_tm_gates,ID)
pred_dr_gates <- cbind(pred_dr_gates,ID)
  


list_res_TM <- vector("list", M)
list_res_DR <- vector("list", M)

for(j in 1:S){ 
  
  theta_set <- ifelse(settings[t,4]=="constant",settings[t,5],settings[t,4])
  
  dataset <- datagen(y="con", N=settings[t,1],k=settings[t,2],random_d=settings[t,3],theta=theta_set,var=1)
  dataset$ID <- c(1:N)
  
  k <- ncol(dataset)-3
  covariates <- c(paste0("V", 1:k)) 
  covariates
 
  
  dataset$d <- as.factor(ifelse(dataset$d==1,1,0))
  
  
  
  for(i in 1:M){
  
  
  
  
  ##### Parameter and datasets ##### 
  
  trainIndex <- createDataPartition(dataset$d, p = .5, list = FALSE) 
  df_aux <- dataset[trainIndex,] 
  df_main  <- dataset[-trainIndex,]
  
  
  
  # On the auxiliary sample
  # -----------------------
  
  # Propensity score using regression forests
  
  
  # set up the cross-validated hyper-parameter search
  xgb_grid_1 = expand.grid(
    nrounds = 200, 
    max_depth = c(2, 3), 
    eta = c(0.001), 
    gamma = c(1),
    colsample_bytree = c(1.0), 
    min_child_weight = c(0.5, 1, 1.5),
    subsample = 1
  )
  
  # pack the training control parameters
  xgb_trcontrol_1 = trainControl(
    method = "cv",
    number = 3
  )
  
  
  # train the model for each parameter combination in the grid,   using CV to evaluate
  xgb_prop = train(
   x = df_aux[,-c(1:3,ncol(dataset))],
  y = df_aux$d,
  trControl = xgb_trcontrol_1,
  tuneGrid = xgb_grid_1,
  method = "xgbTree"
  )
  
  
  
  
  # Conditional mean proxy using regression forests
  
  aux_1 <- df_aux[which(df_aux$d==1),]
  aux_0 <- df_aux[which(df_aux$d==0),]
  
  
  xgb_cate_1 = train(
    x = aux_1[,-c(1:3,ncol(dataset))],
    y = aux_1$y,
    trControl = xgb_trcontrol_1,
    tuneGrid = xgb_grid_1,
    method = "xgbTree"
  )
  xgb_cate_0 = train(
    x = aux_0[,-c(1:3,ncol(dataset))],
    y = aux_0$y,
    trControl = xgb_trcontrol_1,
    tuneGrid = xgb_grid_1,
    method = "xgbTree"
  )
  
  
  #varImp(xgb_cate_0)
  #varImp(xgb_cate_1)
  
  
  # On the main sample
  # -----------------------
  
  # Propensity score offset W - e(X)
  p <- predict(xgb_prop, data.matrix(df_main[,-c(1:3,ncol(dataset))]),type="prob")[,2]
  #p  <- rep((nrow(df_aux[df_aux[,"d"]==1,]) + nrow(df_main[df_main[,"d"]==1,]))/(nrow(df_aux) + nrow(df_main)), nrow(df_main))  # for random d
  
  df_main$d <- as.numeric(as.character(df_main$d)) - p
  
  
  # Group variable by predicted CATE quintile 
  
  
  y1 <- predict(xgb_cate_1, data.matrix(df_main[,-c(1:3,ncol(dataset))]))
  y0 <- predict(xgb_cate_0, data.matrix(df_main[,-c(1:3,ncol(dataset))]))
  
  
  # Score function for distribution
  df_main$G1 <- (y1 - y0)
  prop[,i][df_main$ID] <- df_main$G1
  
  
  # Divide observations into k-tiles
  df_main$S <- (y1 - y0) 
  S2 <- df_main$S +runif(length(df_main$S), 0, 0.00001) # Include white noise to guarantee that the score (S) differs from the baseline effect
  #breaks    <- quantile(S2, seq(0,1, 0.2),  include.lowest =T) # Quantiles create 5 groups by default - no parameter necessary 
  #breaks[1] <- breaks[1] - 0.001 # Offset for lower tails 
  #breaks[6] <- breaks[6] + 0.001 # Offset for upper tails
  SG <- cut(S2, breaks = ntile)
  
  
  
  SGX       <- model.matrix(~-1+SG) # -1 Ereases the Intercept. Possible is also to keep the Intercept. 
  DSG       <- data.frame(as.numeric(I(as.numeric(df_main[,"d"])))*SGX)
  
  colnames(DSG) <- c("G1", "G2", "G3", "G4", "G5")
  df_main[,c("B", "S", "G1", "G2", "G3", "G4", "G5", "weight")] <- cbind(y0, df_main$S, DSG$G1, DSG$G2, DSG$G3, DSG$G4, DSG$G5, as.numeric((1/(p*(1-p)))))
  
  
  form1 <- as.formula(paste("y", "~", "B+S+G1+G2+G3+G4+G5 ", sep=""))
  df_main$y <- as.numeric(df_main$y)
  
  
  
  # Now regress on group membership variables
  #XX <- str_c(covariates, collapse="+")
  #fmla <- "y ~" %>% str_c(XX) %>% str_c("+ d:G")
  #model <- glm(fmla, df_main, family = "binomial") Produces log-odds as output 
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




# Propensity score offset W - e(X)
p1 <- predict(xgb_prop, data.matrix(df_aux[,-c(1:3,ncol(dataset))]),type="prob")[,2]
p <- predict(xgb_prop, data.matrix(df_main[,-c(1:3,ncol(dataset))]),type="prob")[,2]

#p  <- rep((nrow(df_aux[df_aux[,"d"]==1,]) + nrow(df_main[df_main[,"d"]==1,]))/(nrow(df_aux) + nrow(df_main)), nrow(df_main))  # for random d



# Group variable by predicted CATE quintile 


y1 <- predict(xgb_cate_1, data.matrix(df_aux[,-c(1:3,ncol(dataset))]))
y0 <- predict(xgb_cate_0, data.matrix(df_aux[,-c(1:3,ncol(dataset))]))


# Doubly Robust
df_aux$d <- as.numeric(ifelse(df_aux$d==1,1,0))

y_mo <- (y1 - y0) + ((df_aux$d*(df_aux$y-y1))/p1) - ((1-df_aux$d)*(df_aux$y-y0)/(1-p1))

xgb_dr = train(
  x = df_aux[,-c(1:3,ncol(dataset))],
  y = y_mo,
  trControl = xgb_trcontrol_1,
  tuneGrid = xgb_grid_1,
  method = "xgbTree"
)

score_dr <- predict(xgb_dr, data.matrix(df_main[,-c(1:3,ncol(dataset))]))

# Score function for distribution
df_main$G1 <- (score_dr)
#prop[,i][df_main$ID] <- df_main$G1


# Divide observations into k-tiles
df_main$S <- (score_dr) 
S2 <- df_main$S +runif(length(df_main$S), 0, 0.00001) # Include white noise to guarantee that the score (S) differs from the baseline effect
#breaks    <- quantile(S2, seq(0,1, 0.2),  include.lowest =T) # Quantiles create 5 groups by default - no parameter necessary 
#breaks[1] <- breaks[1] - 0.001 # Offset for lower tails 
#breaks[6] <- breaks[6] + 0.001 # Offset for upper tails
SG <- cut(S2, breaks = ntile)

df_main$d <- as.numeric(as.character(df_main$d)) - p


SGX       <- model.matrix(~-1+SG) # -1 Ereases the Intercept. Possible is also to keep the Intercept. 
DSG       <- data.frame(as.numeric(I(as.numeric(df_main[,"d"])))*SGX)

colnames(DSG) <- c("G1", "G2", "G3", "G4", "G5")
df_main[,c("S", "G1", "G2", "G3", "G4", "G5", "weight")] <- cbind(df_main$S, DSG$G1, DSG$G2, DSG$G3, DSG$G4, DSG$G5, as.numeric((1/(p*(1-p)))))


form1 <- as.formula(paste("y", "~", "G1+G2+G3+G4+G5 ", sep=""))
df_main$y <- as.numeric(df_main$y)



# Now regress on group membership variables
#XX <- str_c(covariates, collapse="+")
#fmla <- "y ~" %>% str_c(XX) %>% str_c("+ d:G")
#model <- glm(fmla, df_main, family = "binomial") Produces log-odds as output 
model <- lm(form1,df_main)
groups <- c(paste0("G",1:ntile))
groups <- dput(as.character(groups))

thetahat1 <- model%>% 
  .$coefficients %>%
  .[groups]

# Confidence intervals
cihat <- confint(model,level=0.9)[groups,]
list_res_DR[[i]] <- tibble(coefficient = dput(as.character(c(paste0("Group", 1:ntile)))),
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
pred_dr_gates[,i][gates_y$ID] <- gates_y$thetahat1


  
  
  
  
  }
  
  
  GATES_TM <- list_res_TM[] %>% 
    bind_rows %>%
    na.omit() %>%
    group_by(coefficient) %>%
    summarize_all(median) 
  
  
  pred_tm_gates_median <- pred_tm_gates[,-ncol(pred_tm_gates)]
  apply(pred_tm_gates_median,1, median, na.rm = TRUE) # Calculate the row median which is then used to classify each obs. into a "group". 
  
  error_matrix[j,1] <- mean(abs(dataset$theta-pred_tm_gates_median),na.rm=T)
  
  

GATES_DR <- list_res[] %>% 
  bind_rows %>%
  group_by(coefficient) %>%
  summarize_all(median)

pred_dr_gates_median <- pred_dr_gates[,-ncol(pred_dr_gates)]
apply(pred_dr_gates_median,1, median, na.rm = TRUE) # Calculate the row median which is then used to classify each obs. into a "group". 

error_matrix[j,2] <- mean(abs(dataset$theta-pred_dr_gates_median),na.rm=T)
error_matrix[j,3] <- mean (dataset$theta)
 

print(paste0("................................... ","The current iteration is: ", j, " out of " ,S))
}

error_result[[t]] <- error_matrix
GATES_result[[t]] <- c(GATES_TM, GATES_DR)

print(paste0("................................... ","This is DGP : ", t, " out of " ,nrow(settings)))

}


##########################



error_result


error_all <- matrix(NA,S*nrow(settings),4)
error_all


b=0

for(j in 1:nrow(settings)){
  for(i in 1:S){

  error_all[i+b,1] <- error_result[[j]][i]
  error_all[i+b,2] <- error_result[[j]][i,2]
  error_all[i+b,3] <- error_result[[j]][i,3]
  error_all[i+b,4] <- as.numeric(j)
}
 b = b+S 
  
}

error_all <- as.data.frame(error_all)
colnames(error_all) <- c("GATES_MAE", "DR_MAE", "TRUE_TE", "SETTING")


# wilcox-test for mean differences between groups
#wilcox.test(error_matrix[1:10,1],error_matrix[1:10,2]) # Better use non-parametric test since the assumption that X and Y are ~ N(.) is not fulfilled. 

 

ggplot(error_all, aes(x=GATES_MAE, y=DR_MAE)) + 
  xlim(0.0,0.3) +
  ylim(0.0,0.3) +
  geom_abline(mapping= aes(intercept=0.0,slope = 1.0, color="45 Degree line")) +
  scale_colour_manual(values="red") +
  labs(colour="") +
  geom_point() + 
theme_cowplot() +
facet_wrap( ~ SETTING, scales="free", ncol=3)  + # Facet wrap with common scales 
  guides(fill = FALSE, color = FALSE, linetype = FALSE, shape = FALSE) +
  labs(x = "Two-Model Estimator", y = "Doubly-Robust Estimator")



