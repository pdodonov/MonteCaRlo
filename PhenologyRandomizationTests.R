library(CircStats)


#### Script for permutation test of phenology
### Input file: 3 or more columns
### Col 1: number of weeks (or days, or months) since the beginning of the experiment
### Col 2: if that individual is .1 (1) or .0 (0)
### Cols 3 forth: proportion of branches with each phenophase
# Output will be in the form of a list; each element of the list will be one phenological characteristic. Within each element:
# Rows: names of the parameters
# Column 1: .0 plants; column 2: .1 plants; column 3: .0 - .1; column 4: prob (.0 < .1; column 5: prob (.0 > 2); column 6: two-tailed significance.
#Function arguments:
  #Nvar: number of variables in the dataset
  #N: number of individuals in each group; "detect" for the program to detect it automatically from the second column.
  #Nperm: number of permutations; 1000 permutations means 999 permutations + the original data.
  #Perm: T (assess significance by permutation tests) or F
  #test.AI: T (calculate and test activity index) or F
  #test.II: T (calculate or teste intensity index) or F
  #test.max: differences in the maximum values
  #test.mean: differences in the mean values
  #test.mwp: differences in the mean values, excluding zeroes
  #test.median: differences in the median values
  #test.time: differences in the mean angle (mean date)
  #test.dispersion: differences in the temporal dispersion (circular vector length)
  
angle.rep=function(x) {
  angle=rep(x[,1], times=x[,2]*100)  
  }
angle.mean = function(x) {
  angle=rep(x[,1], times=x[,2]*100)  
  return(circ.mean(rad(angle)))
  }
angle.disp = function(x) {
  angle=rep(x[,1], times=x[,2]*100)
  return(est.rho(rad(angle)))
  }
activity.index <- function (x) { #a short function to calculate the activity index.
  a = sum(x>0)/length(x)
  return(a)
  }

### The function below compares the phenology of two groups.  
RTPD = function(x, Nvar = 1, N1 = 30, N2 = 30, Nperm = 5000, perm = TRUE, test.AI = TRUE, test.II = T, test.max = T, test.mean = T, test.mwp = T, test.median = T, test.time=T, test.dispersion=T) {
  Na = N1
  Nb = N2
  if (Nvar > 1) result = list()
  Ntot = Na + Nb
  dados=x
  Nobs = length(unique(dados[,1]))
  main = x[order(x[,1], x[,2]),] #order by week and separate between the two categories 
  main.list = list() #Organizes the variables into a list
  for (i in 1:Nvar) {
    main.list[[i]] = main[,2+i]
    }
  groups=main[,2]
  times=main[,1]
  names(main.list)=colnames(main)[3:(2+Nvar)]
  result = list()
  for (i in 1:Nvar) {
    main.use = main.list[[i]]
    if (test.II == TRUE) {
      II.0=aggregate(subset(main.use, groups==0), list(subset(times, groups==0)), mean)[,2]
      II.1=aggregate(subset(main.use, groups==1), list(subset(times, groups==1)), mean)[,2]
      if (test.max == TRUE) {
        IImaxU = max (II.0)
        IImaxB = max (II.1)
        IImaxdif = IImaxU - IImaxB
        if (perm == TRUE) {
          IImaxdif.perm = numeric(Nperm)
          IImaxdif.perm[1] = IImaxdif
          }
        }
      if (test.mean == TRUE) {
        IImeanU = mean (II.0) #Mean intensity and activity indices for .1 and .0 individuals
        IImeanB = mean(II.1)
        IImeandif = IImeanU - IImeanB
        if (perm == TRUE) {
          IImeandif.perm = numeric(Nperm)
          IImeandif.perm[1] = IImeandif
          }
        }
      if (test.mwp == TRUE) {
        IImwpU = mean (II.0[II.0>0])
        IImwpB = mean (II.1[II.1>0])
        IImwpdif = IImwpU - IImwpB
        if (perm == TRUE) {
          IImwpdif.perm = numeric(Nperm)
          IImwpdif.perm[1] = IImwpdif
          }
        }
      if (test.median == TRUE) {
        IImedU = median (II.0)
        IImedB = median (II.1)
        IImeddif = IImedU - IImedB
        if (perm == TRUE) {
          IImeddif.perm = numeric(Nperm)
          IImeddif.perm[1] = IImeddif
          }
        }
      if (test.time == TRUE) {
        IItimeU=angle.mean(data.frame(unique(times), II.0))
        IItimeB=angle.mean(data.frame(unique(times), II.1))
        IItimedif = as.numeric(circ.range(c(IItimeU, IItimeB)))
        if (perm == TRUE) {
          IItimedif.perm = numeric(Nperm)
          IItimedif.perm[1] = IItimedif
          }
        }
      if (test.dispersion == TRUE) {
        IIconcU=angle.disp(data.frame(unique(times), II.0))
        IIconcB=angle.disp(data.frame(unique(times), II.1))
        IIconcdif = IIconcU-IIconcB
        if (perm == TRUE) {
          IIconcdif.perm = numeric(Nperm)
          IIconcdif.perm[1] = IIconcdif
          }
        }
      }
    if (test.AI == TRUE) {
      AI.0=aggregate(subset(main.use, groups==0), list(subset(times, groups==0)), activity.index)[,2]
      AI.1=aggregate(subset(main.use, groups==1), list(subset(times, groups==1)), activity.index)[,2]
      if (test.max == TRUE) {
        AImaxU = max (AI.0)
        AImaxB = max (AI.1)
        AImaxdif = AImaxU - AImaxB
        if (perm == TRUE) {
          AImaxdif.perm = numeric(Nperm)
          AImaxdif.perm[1] = AImaxdif
          }
        }
      if (test.mean == TRUE) {
        AImeanU = mean (AI.0)
        AImeanB = mean(AI.1)
        AImeandif = AImeanU - AImeanB
        if (perm == TRUE) {
          AImeandif.perm = numeric(Nperm)
          AImeandif.perm[1] = AImeandif
          }
        }
      if (test.mwp == TRUE) {
        AImwpU = mean (AI.0[AI.0>0])
        AImwpB = mean (AI.1[AI.1>0])
        AImwpdif = AImwpU - AImwpB
        if (perm == TRUE) {
          AImwpdif.perm = numeric(Nperm)
          AImwpdif.perm[1] = AImwpdif
          }
        }
      if (test.median == TRUE) {
        AImedU = median (AI.0)
        AImedB = median (AI.1)
        AImeddif = AImedU - AImedB
        if (perm == TRUE) {
          AImeddif.perm = numeric(Nperm)
          AImeddif.perm[1] = AImeddif
          }
        }
      if (test.time == TRUE) {
        AItimeU=angle.mean(data.frame(unique(times), AI.0))
        AItimeB=angle.mean(data.frame(unique(times), AI.1))
        AItimedif = as.numeric(circ.range(c(AItimeU, AItimeB)))
        if (perm == TRUE) {
          AItimedif.perm = numeric(Nperm)
          AItimedif.perm[1] = AItimedif
          }
        }
      if (test.dispersion == TRUE) {
        AIconcU=angle.disp(data.frame(unique(times), AI.0))
        AIconcB=angle.disp(data.frame(unique(times), AI.1))
        AIconcdif = AIconcU-AIconcB
        if (perm == TRUE) {
          AIconcdif.perm = numeric(Nperm)
          AIconcdif.perm[1] = AIconcdif
          }
        }
      }
    if (perm == TRUE) { ##Permutation
      for (j in 2: Nperm) {
        groups.perm = rep(sample(groups[1:Ntot], replace=F), times=length(unique(times)))
          if(test.II == TRUE) {
            II.0.perm=aggregate(subset(main.use, groups.perm==0), list(subset(times, groups.perm==0)), mean)[,2]
            II.1.perm=aggregate(subset(main.use, groups.perm==1), list(subset(times, groups.perm==1)), mean)[,2]
            if(test.max == TRUE) {
              IImaxU.perm = max (II.0.perm)
              IImaxB.perm = max (II.1.perm)
              IImaxdif.perm[j] = IImaxU.perm - IImaxB.perm
              }
            if (test.mean == TRUE) {
              IImeanU.perm = mean (II.0.perm)
              IImeanB.perm = mean(II.1.perm)
              IImeandif.perm[j] = IImeanU.perm - IImeanB.perm
              }
            if (test.mwp == TRUE) {          
              IImwpU.perm = mean (II.0.perm[II.0.perm>0])
              IImwpB.perm = mean (II.1.perm[II.1.perm>0])
              IImwpdif.perm[j] = IImwpU.perm - IImwpB.perm
              }
            if (test.median == TRUE) {
              IImedU.perm = median (II.0.perm)
              IImedB.perm = median(II.1.perm)
              IImeddif.perm[j] = IImedU.perm - IImedB.perm
              }
            if (test.time == TRUE) {
              IItimeU.perm=angle.mean(data.frame(unique(times), II.0.perm))
              IItimeB.perm=angle.mean(data.frame(unique(times), II.1.perm))
              IItimedif.perm[j] = as.numeric(circ.range(c(IItimeU.perm, IItimeB.perm)))
              }
            if (test.dispersion == TRUE) {
              IIconcU.perm=angle.disp(data.frame(unique(times), II.0.perm))
              IIconcB.perm=angle.disp(data.frame(unique(times), II.1.perm))
              IIconcdif.perm[j] = IIconcU.perm-IIconcB.perm
              }
            }
        if (test.AI == TRUE) {
          AI.0.perm=aggregate(subset(main.use, groups.perm==0), list(subset(times, groups.perm==0)), activity.index)[,2]
          AI.1.perm=aggregate(subset(main.use, groups.perm==1), list(subset(times, groups.perm==1)), activity.index)[,2]
          if(test.max == TRUE) {
            AImaxU.perm = max (AI.0.perm)
            AImaxB.perm = max (AI.1.perm)
            AImaxdif.perm[j] = AImaxU.perm - AImaxB.perm
            }
          if (test.mean == TRUE) {
            AImeanU.perm = mean (AI.0.perm)
            AImeanB.perm = mean(AI.1.perm)
            AImeandif.perm[j] = AImeanU.perm - AImeanB.perm
            }  
          if (test.mwp == TRUE) {          
            AImwpU.perm = mean (AI.0.perm[AI.0.perm>0])
            AImwpB.perm = mean (AI.1.perm[AI.1.perm>0])
            AImwpdif.perm[j] = AImwpU.perm - AImwpB.perm
            }
          if (test.median == TRUE) {
            AImedU.perm = median (AI.0.perm)
            AImedB.perm = median(AI.1.perm)
            AImeddif.perm[j] = AImedU.perm - AImedB.perm
            }
          if (test.time == TRUE) {
            AItimeU.perm=angle.mean(data.frame(unique(times), AI.0.perm))
            AItimeB.perm=angle.mean(data.frame(unique(times), AI.1.perm))
            AItimedif.perm[j] = as.numeric(circ.range(c(AItimeU.perm, AItimeB.perm)))
            }
          if (test.dispersion == TRUE) {
            AIconcU.perm=angle.disp(data.frame(unique(times), AI.0.perm))
            AIconcB.perm=angle.disp(data.frame(unique(times), AI.1.perm))
            AIconcdif.perm[j] = AIconcU.perm-AIconcB.perm
            }
          }
        print(c(i,j))
        }
      }
    if (test.II == TRUE) {
      if (test.max == TRUE) {
        IImaxdif.lower = sum(IImaxdif.perm <= IImaxdif)/Nperm
        IImaxdif.higher = sum(IImaxdif.perm >= IImaxdif)/Nperm
        IImaxdif.signif = min(IImaxdif.lower, IImaxdif.higher)*2
        }
      if (test.mean == TRUE) {  
        IImeandif.lower = sum(IImeandif.perm <= IImeandif)/Nperm
        IImeandif.higher = sum(IImeandif.perm >= IImeandif)/Nperm
        IImeandif.signif = min(IImeandif.lower,IImeandif.higher)*2
        }
      if (test.mwp == TRUE) {
        IImwpdif.lower = sum(IImwpdif.perm <= IImwpdif)/Nperm
        IImwpdif.higher = sum(IImwpdif.perm >= IImwpdif)/Nperm
        IImwpdif.signif = min(IImwpdif.lower,IImwpdif.higher)*2
        }
      
      if (test.median == TRUE) {
        IImeddif.lower = sum(IImeddif.perm <= IImeddif)/Nperm
        IImeddif.higher = sum(IImeddif.perm >= IImeddif)/Nperm
        IImeddif.signif = min(IImeddif.lower,IImeddif.higher)*2
        }
      if (test.time == TRUE) {
        IItimedif.lower = NA
        IItimedif.higher = NA
        IItimedif.signif = sum(IItimedif.perm >= IItimedif)/Nperm
        }
      if (test.dispersion == TRUE) {
        IIconcdif.lower = sum(IIconcdif.perm <= IIconcdif)/Nperm
        IIconcdif.higher = sum(IIconcdif.perm >= IIconcdif)/Nperm
        IIconcdif.signif = min(IIconcdif.lower,IIconcdif.higher)*2
        }
      }
    if (test.AI == TRUE) {
      if (test.max == TRUE) {
        AImaxdif.lower = sum(AImaxdif.perm <= AImaxdif)/Nperm
        AImaxdif.higher = sum(AImaxdif.perm >= AImaxdif)/Nperm
        AImaxdif.signif = min(AImaxdif.lower, AImaxdif.higher)*2
        }
      if (test.mean == TRUE) {  
        AImeandif.lower = sum(AImeandif.perm <= AImeandif)/Nperm
        AImeandif.higher = sum(AImeandif.perm >= AImeandif)/Nperm
        AImeandif.signif = min(AImeandif.lower,AImeandif.higher)*2
        }
      if (test.mwp == TRUE) {
        AImwpdif.lower = sum(AImwpdif.perm <= AImwpdif)/Nperm
        AImwpdif.higher = sum(AImwpdif.perm >= AImwpdif)/Nperm
        AImwpdif.signif = min(AImwpdif.lower,AImwpdif.higher)*2
        }
      if (test.median == TRUE) {
        AImeddif.lower = sum(AImeddif.perm <= AImeddif)/Nperm
        AImeddif.higher = sum(AImeddif.perm >= AImeddif)/Nperm
        AImeddif.signif = min(AImeddif.lower,AImeddif.higher)*2
        }
      if (test.time == TRUE) {
        AItimedif.lower = NA
        AItimedif.higher = NA
        AItimedif.signif = sum(AItimedif.perm >= AItimedif)/Nperm
        }
      if (test.dispersion == TRUE) {
        AIconcdif.lower = sum(AIconcdif.perm <= AIconcdif)/Nperm
        AIconcdif.higher = sum(AIconcdif.perm >= AIconcdif)/Nperm
        AIconcdif.signif = min(AIconcdif.lower,AIconcdif.higher)*2
        }
      }
    answer = matrix(0, nrow=((test.max+test.mean+test.mwp+test.median+test.time+test.dispersion) * (test.II + test.AI) ) ,ncol=6)
    measures=c(c(test.max, test.mean, test.mwp, test.median, test.time, test.dispersion)*test.II, c(test.max, test.mean, test.mwp, test.median, test.time, test.dispersion)*test.AI)

    rownames(answer)=c("IImean","IImax","IImwp","IImeadian", "IItime", "IIconcentration","AImean","AImax","AImwp","AImedian","AItime","AIconcentration")[measures==1]
    colnames(answer)=c(".0",".1","Difference","Prob. lower", "Prob. higher", "Significance")
    line = 1
    if (test.AI == T) {
      if (test.max == T) {
        answer[line,] = c(AImaxU, AImaxB, AImaxdif, AImaxdif.lower, AImaxdif.higher, AImaxdif.signif)
        rownames(answer)[line] = "AImax"
        line = line + 1
        }
      if (test.mean == T) {
        answer[line,] = c(AImeanU, AImeanB, AImeandif, AImeandif.lower, AImeandif.higher, AImeandif.signif)
        rownames(answer)[line] = "AImean"
        line = line + 1
        }
      if (test.mwp == T) {
        answer[line,] = c(AImwpU, AImwpB, AImwpdif, AImwpdif.lower, AImwpdif.higher, AImwpdif.signif)
        rownames(answer)[line] = "AImwp"
        line = line + 1
        }
      if (test.median == T) {
        answer[line,] = c(AImedU, AImedB, AImeddif, AImeddif.lower, AImeddif.higher, AImeddif.signif)
        rownames(answer)[line] = "AImedian"
        line = line + 1
        }
      if (test.time == T) {
        answer[line,] = c(AItimeU, AItimeB, AItimedif, AItimedif.lower, AItimedif.higher, AItimedif.signif)
        rownames(answer)[line] = "AItimes"
        line = line + 1
        }
      if (test.dispersion == T) {
        answer[line,] = c(AIconcU, AIconcB, AIconcdif, AIconcdif.lower, AIconcdif.higher, AIconcdif.signif)
        rownames(answer)[line] = "AIconcentration"
        line = line + 1
        }
      }
    if (test.II == T) {
      if (test.max == T) {
        answer[line,] = c(IImaxU, IImaxB, IImaxdif, IImaxdif.lower, IImaxdif.higher, IImaxdif.signif)
        rownames(answer)[line] = "IImax"
        line = line + 1
        }
      if (test.mean == T) {
        answer[line,] = c(IImeanU, IImeanB, IImeandif, IImeandif.lower, IImeandif.higher, IImeandif.signif)
        rownames(answer)[line] = "IImean"
        line = line + 1
        }
      if (test.mwp == T) {
        answer[line,] = c(IImwpU, IImwpB, IImwpdif, IImwpdif.lower, IImwpdif.higher, IImwpdif.signif)
        rownames(answer)[line] = "IImwp"
        line = line + 1
        }
      if (test.median == T) {
        answer[line,] = c(IImedU, IImedB, IImeddif, IImeddif.lower, IImeddif.higher, IImeddif.signif)
        rownames(answer)[line] = "IImedian"
        line = line + 1
        }
      if (test.time == T) {
        answer[line,] = c(IItimeU, IItimeB, IItimedif, IItimedif.lower, IItimedif.higher, IItimedif.signif)
        rownames(answer)[line] = "IItimes"
        line = line + 1
        }
      if (test.dispersion == T) {
        answer[line,] = c(IIconcU, IIconcB, IIconcdif, IIconcdif.lower, IIconcdif.higher, IIconcdif.signif)
        rownames(answer)[line] = "IIconcentration"
        line = line + 1
        }
      }
    if (Nvar > 1) result[[i]]=answer
    }
  if (Nvar > 1) names(result) = names(main.list)
  if (Nvar > 1) result else answer
  }
  

### The function below assesses the significance of the mean angle, for each group separately.
RTMA <- function(x, Nvar = 1, N1 = 30, N2 = 30, Nperm = 5000, perm = TRUE) { #performs automatically for activity and intensity indices
  Na <- N1
  Nb <- N2
  result <- matrix(ncol = 8, nrow = Nvar)
  colnames(result) <- c("II_Vector_U", "II_Signif_U", "II_Vector_B", "II_Signif_B", "AI_vector_U", "AI_signif_U", "AI_vector_B", "AI_signif_B")
  rownames(result) <- colnames(x)[3:(2+Nvar)]
  Ntot <- Na + Nb
  dados <- x
  Nobs <- length(unique(dados[,1]))
  main <- x[order(x[,1], x[,2]),] #order by week and separate between the two categories 
  main.list <- list() #Organizes the variables into a list
  for (i in 1:Nvar) {
    main.list[[i]] <- main[,2+i]
    }
  groups<-main[,2]
  times<-main[,1]
  names(main.list)<-colnames(main)[3:(2+Nvar)]
  for (i in 1:Nvar) {
    main.use <- main.list[[i]]
    
    II.0 <- aggregate(subset(main.use, groups==0), list(subset(times, groups==0)), mean)[,2]
    II.1 <- aggregate(subset(main.use, groups==1), list(subset(times, groups==1)), mean)[,2]
      
    IIconcU <- angle.disp(data.frame(unique(times), II.0))
    IIconcB <- angle.disp(data.frame(unique(times), II.1))
      
    IIconcU.perm <- numeric(Nperm)
    IIconcB.perm <- numeric(Nperm)
    IIconcU.perm[1] <- IIconcU
    IIconcB.perm[1] <- IIconcB

    for(j in 2:Nperm) {
      II.0.perm <- sample(II.0)
      II.1.perm <- sample(II.1)

      IIconcU.perm[j] <- angle.disp(data.frame(unique(times),II.0.perm))
      IIconcB.perm[j] <- angle.disp(data.frame(unique(times),II.1.perm))
    }
    IIconcU.signif <- sum(IIconcU.perm >= IIconcU) / Nperm
    IIconcB.signif <- sum(IIconcB.perm >= IIconcB) / Nperm

    AI.0 <- aggregate(subset(main.use, groups==0), list(subset(times, groups==0)), activity.index)[,2]
    AI.1 <- aggregate(subset(main.use, groups==1), list(subset(times, groups==1)), activity.index)[,2]

    AIconcU <- angle.disp(data.frame(unique(times), AI.0))
    AIconcB <- angle.disp(data.frame(unique(times), AI.1))
      
    AIconcU.perm <- numeric(Nperm)
    AIconcB.perm <- numeric(Nperm)
    AIconcU.perm[1] <- AIconcU
    AIconcB.perm[1] <- AIconcB

    for(j in 2:Nperm) {
      AI.0.perm <- sample(AI.0)
      AI.1.perm <- sample(AI.1)

      AIconcU.perm[j] <- angle.disp(data.frame(unique(times),AI.0.perm))
      AIconcB.perm[j] <- angle.disp(data.frame(unique(times),AI.1.perm))
    }
    AIconcU.signif <- sum(AIconcU.perm >= AIconcU) / Nperm
    AIconcB.signif <- sum(AIconcB.perm >= AIconcB) / Nperm

    result[i,] <- c(IIconcU, IIconcU.signif, IIconcB, IIconcB.signif, AIconcU, AIconcU.signif, AIconcB, AIconcB.signif)
  }
  return(result)
}
