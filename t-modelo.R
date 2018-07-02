  library(mixtools)
  library(mvtnorm)
  library(abind)
  library(MASS)
  
  
  d <- c()
  
  d1 = mvrnorm(500,mu=c(10,10.5),Sigma = diag(2)/100)
  d2 = mvrnorm(1000,mu=c(8.0,8),Sigma = diag(2)*0.8)
  d3 = mvrnorm(800,mu=c(0,0),Sigma = diag(2)*0.5)
  d4 = mvrnorm(900,mu=c(2,1),Sigma = diag(2)/3)
  d = rbind(d1,d2,d3,d4)
  colnames(d)<-c("X","Y") 
  d <-data.frame(d)
  plot(d[,c(1,2)])
  
  
  
  data(iris)
  #Carregando os dados.
  dados <-  t(d[sample(1:dim(d)[1]),c(1,2)])
  #dados <-  t(iris[sample(1:dim(iris)[1]),c(1,3)])
  
  #Funções
  CalculoMuC <- function(n,v0,k0,d,cov,me,cov0,me0){
    
    fator1 = (k0*n)/(n+k0)
    fator2 = n+k0
    fator3 = v0+n-d+1  
    
    mu = (k0*me0+n*me)/(fator2)
    cv = (cov0+cov+fator1*((me-me0)%*%t(me-me0)))*(fator2+1)
    cv = cv/(fator3*(fator2+1))
    
    return(list(mu = mu,
                cv = cv,
                gl = fator3))
  }
  
  AtualizarMediaCov <- function(exemplo, Media, 
                                Cov, grupo, d, n,
                                inicial=TRUE){
    if(inicial == TRUE){
      Media         <- rbind(Media , exemplo)
      Cov           <- abind(Cov, matrix(0,nrow=d,ncol=d))
    }else{
      Media[grupo,] <- (Media[grupo,]*n + t(exemplo)) / (n + 1)
      Cov[,,grupo]  <- (Cov[,,grupo]*n + n*((exemplo - Media[grupo,] )%*%t((exemplo - Media[grupo,] ))/(n+1)))/(n+1)
    }
    
    return(list(Media = Media,
                Cova  = Cov))
  } 
  
  
  #------------------------------------------
  d     <- dim(dados)[1]
  
  #definindo parâmetros iniciais
  alpha   <- 1
  me0     <- rep(0,d)
  cov0    <- diag(d)*10
  v0      <- d + 2
  k0      <- 0.01
  
  #Vetores de atribuições dos grupos
  Atrib   <- c(1)
  NGrup   <- c(1)
  Grups   <- c(1)
  meGrup  <- matrix(NA, nrow=0,ncol=d)
  covGrup <- array(NA, dim=c(d, d, 0))

  
  
  #Atribuindo o primeiro exemplo no primeiro grupo
  result  <- AtualizarMediaCov(exemplo = dados[,1],
                               Media   = meGrup,
                               Cov     = covGrup,
                               d       = d,
                               inicial = TRUE)
  meGrup  <- result$Media
  covGrup <- result$Cov
    
  graf <- c()
  
  #Simulando Streaming
  for(i in 2:dim(dados)[2]){
  
    
    probs <- c()
    prob  <- 0
    for(k in 1:length(Grups)){
      
    #Gerando os parâmetros
    
    valores <- CalculoMuC(n    = NGrup[k],
                          v0   = v0,
                          k0   = k0,
                          d    = d,
                          cov  = covGrup[,,k],
                          me   = meGrup[k,],
                          me0  = me0,
                          cov0 = cov0)

      prob <- log(NGrup[k])
      prob <- prob + dmvt(dados[,i],sigma=valores$cv,delta=valores$mu,df=valores$gl)
      probs <- c(probs,prob)
    }
    
    #Calculando um novo grupo
      cv = (cov0*(k0+1))/(k0*(v0-d+1))
      prob <- log(alpha)
      prob <- prob + dmvt(dados[,i],sigma=cv,delta=me0,df=v0-d+1)
      probs <- c(probs,prob)
      
      #Normalizando  
      probs <- probs-max(probs)
      probs <- exp(probs)
      probs <- probs/sum(probs)
      
      escolha <- base::sample(1:(length(Grups)+1), 1, prob=probs)
      Atrib   <- c(Atrib, escolha)
      graf    <- rbind(graf,c(dados[,i],escolha))
      
      if(escolha == (length(Grups)+1)){
        
        Grups <- c(Grups,escolha)
        NGrup <- c(NGrup,1)
        
        #Novo grupo
        result  <- AtualizarMediaCov(exemplo = dados[,i],
                                     Media   = meGrup,
                                     Cov     = covGrup,
                                     d       = d,
                                     inicial = TRUE)
        meGrup  <- result$Media
        covGrup <- result$Cov
      }else{
        NGrup[escolha] <- NGrup[escolha] + 1
        
        #Atualiza grupo
        result  <- AtualizarMediaCov(exemplo = dados[,i],
                                     Media   = meGrup,
                                     Cov     = covGrup,
                                     d       = d,
                                     grupo   = escolha,
                                     n       = NGrup[escolha],
                                     inicial = FALSE)
        meGrup  <- result$Media
        covGrup <- result$Cov
      }


      
      #arq = paste("exemplo",i,".png",sep = "")
      #plot.new()
      #png(file=arq, width=400, height=400)

     #dev.off()      
  }
  
  
  tamanho = length(Grups)
  colnames(graf)<-c("X","Y","Z")
  graf <-data.frame(graf)
  #plot(graf[,c(1,2)],col=graf$Z)
  plot(graf[,c(1,2)],col=col=graf$Z))
  for(k in 1:length(Grups)){
    ellipse(mu=meGrup[k,], 
            sigma=(covGrup[,,k]),  
            alpha = .25, lwd=3, npoints = 250, 
            col=palette()[k+1])
  }


  #system("convert -delay 100 *.png example_1.gif")
  #file.remove(list.files(pattern="*.png"))
