require(ggplot2)

errorful.experiment1 <- function(tau,kappa, n, nmc){
   h.kappa <- (2/kappa)^3
   g0.kappa.tau <- 1 - pbeta(tau, 2, kappa, lower.tail = TRUE)
   g1.kappa.tau <- 1 - pbeta(tau, kappa, 2, lower.tail = TRUE)

   alpha1 <- c(200,10,10)
   alpha2 <- c(10,10,200)

   v <- rbinom(nmc, 1, 0.5)

   X1 <- rdirichlet(sum(v==1), alpha1)
   X2 <- rdirichlet(sum(v==0), alpha2)
   X <-  matrix(0, nrow = length(v), ncol = length(alpha1))
   X[which(v == 1),] <- X1
   X[which(v == 0),] <- X2

   labels.hat <- numeric(nmc)

   for(i in 1:nmc){
       v.i <- c(rep(1,n),rep(0,n))
       v.i <- v.i[sample(n,1:n,replace=FALSE)]
       Y.i1 <- rdirichlet(n, alpha1)
       Y.i2 <- rdirichlet(n, alpha2)
       Y.i <- matrix(0, nrow = 2*n, ncol = length(alpha1))
       Y.i[which(v.i == 1),] <- Y.i1
       Y.i[which(v.i == 0),] <- Y.i2
       P.i <- Y.i%*%X[i,]
       P.i.tilde <- h.kappa*(g1.kappa.tau*P.i + g0.kappa.tau*(1 - P.i))
       mm <- rbinom(2*n, 1, P.i.tilde)
       labels <- 2*v.i - 1
       labels.hat[i] <- sign(sum(mm*labels))
   }
   labels.hat[labels.hat == 0] <- 1
   return( sum(labels.hat == (2*v - 1))/nmc)
}

figure4 <- function(p = 0.9,q = 0.1,kappa.1 = 3.5, kappa.2 = 2){
    xseq <- seq(from = 0, to = 1, by = 0.01)

    diff.vec <- pbeta(xseq, kappa.2, kappa.1, lower.tail = TRUE) -
        pbeta(xseq, kappa.1 , kappa.2, lower.tail = TRUE)

    p.vec <- p*diff.vec + pbeta(xseq, kappa.2, kappa.1, lower.tail = FALSE)
    q.vec <- q*diff.vec + pbeta(xseq, kappa.2, kappa.1, lower.tail = FALSE)

    p.vec <- (2/kappa.1)^3*p.vec
    q.vec <- (2/kappa.1)^3*q.vec

    var.p <- p.vec*(1-p.vec)
    var.q <- q.vec*(1 - q.vec)

    aa <- q.vec*sqrt(var.p) - p.vec*sqrt(var.q)

    mean.seq <- p.vec - q.vec
    variance.seq <- var.p + var.q

    df <- data.frame(xseq, mean.seq, variance.seq)

    plot1 <- ggplot(df,aes(x = xseq, y = variance.seq)) + geom_line(color="blue") +
        geom_line(aes(x = xseq, y = mean.seq), color = "red") +
            xlab(expression(tau)) + ylab("")
    plot1
}

Lfn <- function(n,p0,p1){
    sum(dbinom(0:n,n,p0) * pbinom(0:n,n,p1)) -
        sum(dbinom(0:n,n,p0) * dbinom(0:n,n,p1))/2
    ## NB: equal priors (and condition on n0=n1) !!
}
Lfn2 <- function(n,p0,p1){
    sum(dbinom(1:n,n,p0) * pbinom(0:(n-1),n,p1)) +
        sum(dbinom(0:n,n,p0) * dbinom(0:n,n,p1))/2
    ## NB: equal priors (and condition on n0=n1) !!
}
h = function(kappa,hpower=3){
    (2/kappa)^hpower
}

ones = matrix( c(1,1) , nrow=2,ncol=1,byrow=T)
I <- diag(2)
J = matrix(1 , nrow=2,ncol=2)
G0 <- function(kappa,tau) 1-pbeta(tau,2,kappa)
G1 <- function(kappa,tau) 1-pbeta(tau,kappa,2)

B = matrix( c(.9,.1,.1,.9) , nrow=2,ncol=2,byrow=T)
##  doubly stochastic affinity model =>
## graph density rho(G) \approx 1/2 (if pi=[1/2,1/2]) !!

sbmpi = matrix( c(1/2,1/2) , nrow=2,ncol=1,byrow=T)


figure1 <- function(prr=0.05,pgg=0.05,prg=0.01)
{
    set.seed(1111)
    require(igraph)
    
    g <- sbm.game(100,matrix(c(prr,prg,prg,pgg),nrow=2),c(50,50))

    ## vertex color/size
    V(g)$color[1:50] <- "red"
    V(g)$color[51:100] <- "green"
    V(g)$color[40] <- "black" # pick an arbitrary one for black!
    V(g)$size <- 10
    V(g)$size[40] <- V(g)$size[3]*2

    par(mfrow=c(1,2))
    plot(g,layout=layout.fruchterman.reingold(g),vertex.label=NA)

    Xhat <- adjacency.spectral.embedding(g,2)$X
    plot(Xhat,col=V(g)$color,pch=19,cex=V(g)$size/10,
        xlab="Embedding Dimension #1", ylab="Embedding Dimension #2")
    # dev.print(pdf,"Figure1.pdf",width=10,height=5)
}

figure2 <- function(){
    par(mfrow=c(1,1))
    n=50
    kappavec = seq(2.0,8,by=0.1)
    nkappa = length(kappavec)
    tauvec = seq(0.30,0.70,by=0.001)
    ntau = length(tauvec)
    L = matrix(0,nrow=nkappa,ncol=ntau)
    for(kappaindex in 1:nkappa){
        for(tauindex in 1:ntau){
            kappa=kappavec[kappaindex]
            tau=tauvec[tauindex]
            Btilde = h(kappa) * ( (G1(kappa,tau)-G0(kappa,tau))*B + G0(kappa,tau)*J )
            L[kappaindex,tauindex] = Lfn(n/2,Btilde[1,2],Btilde[1,1])
        }
    }

    image(kappavec[3:41],tauvec,(L[3:41,]), xlab=expression(h(kappa)),
          ylab=expression(tau), cex.axis=1.25, cex.lab=1.5, col=heat.colors(120), axes = F)
    contour(kappavec[3:41],tauvec,(L[3:41,]),add=T,lwd=2,labcex=1.5,nlevels=19)
    taustarkappa=NULL ; for(i in 1:nkappa) taustarkappa[i]=tauvec[order(L[i,])[1]]
    points(kappavec[3:41],taustarkappa[3:41],pch=16)
    points(kappavec[3:41],taustarkappa[3:41],type="l",lwd=3)
    ooo = order(t(L))[1]
    kappastarindex = ceiling(ooo/ntau)
    taustarindex = ooo-floor(ooo/ntau)*ntau
    kappastar = kappavec[kappastarindex]   ##   3.5
    taustar = tauvec[taustarindex]         ##   0.6
    Lstar = L[kappastarindex,taustarindex] ##   0.1610279
    points(kappastar,taustar,pch=16,cex=3)
    Ltxt = bquote(L[paste(kappa^"*",",",tau^"*")] == .(round(Lstar,4)))
    axis(2)
    axis(1, seq(2,6,0.5), labels = round(h(seq(2,6,0.5)), 2))

    text(4.4,0.625,Ltxt,cex=3)
}

figure3 <- function(){

    mypar = par(mar=c(4,7,2,1)) 
    kappastar <- 3.5
    n <- 50
    pseq = seq(0,h(kappastar),by=0.001)
    np = length(pseq)
    Lbinomialpair = matrix(1/2,nrow=np,ncol=np)
    for(p0i in 1:(np-1)){
        for(p1i in p0i:np){
            p0 = pseq[p0i]
            p1 = pseq[p1i]
            Lbinomialpair[p0i,p1i] = Lfn(n/2,p0,p1)
        }
    }
    tauvec = seq(0,1,by=0.001)
    ntau = length(tauvec)
    x=y=NULL
    zzz=NULL 
    SVD = matrix(0,nrow=ntau,ncol=2)

    for(tauindex in 1:ntau){
        tau=tauvec[tauindex]
        Btilde = h(kappastar)*((G1(kappastar,tau)-G0(kappastar,tau))*B + G0(kappastar,tau)*J)
        SVD[tauindex,]=svd(Btilde,2,2)$d
        x[tauindex]=Btilde[1,1]
        y[tauindex]=Btilde[1,2]
        zzz[tauindex] = Lfn(n/2,y[tauindex],x[tauindex])
    }
    image(pseq,pseq, (Lbinomialpair), xlab=expression(tilde(B)[paste(1,",",1)]),
          ylab=expression(tilde(B)[paste(1,",",2)]),cex.axis=1.25,
          cex.lab=1.25,col=heat.colors(120))
    contour(pseq,pseq, (Lbinomialpair) ,add=T) # lwd=2,labcex=1.5,nlevels=19)
    abline(coef=c(0,1))
    points(y,x,type="l",lty=1,lwd=3)
    points(y[601],x[601],type="p",pch=16,cex=3) ##  tau = 0.6 = tau^*
    points(y[501],x[501],type="p",pch=16,cex=2) ##  tau = 0.5 = tau_Bayes
    par(mypar)
}

vertexassignment <- function( B = matrix( c(.9,.1,.1,.9) , nrow=2,ncol=2,byrow=T),
                             kappa=3.5, n=50,nmc=10, cvseq=seq(0.30,0.90,by=0.025),
                             tauprob1=0.5,conditionontau=T,r1=6,s1=2,r0=2,s0=6,
                             augmentdiag=T,mc1=0){

    Lhat = matrix(1,nrow=nmc,ncol=length(cvseq))
    J = matrix( 1 , nrow=2,ncol=2)
    for(mc in 1:nmc){
        set.seed(mc+mc1)
        ifelse(conditionontau, (tau = c(rep(1,n*tauprob1),rep(2,n*(1-tauprob1)))),
               (tau = (rbinom(n,1,1-tauprob1)+1)) )
        for(whichcv in 1:length(cvseq)){
            cv = cvseq[whichcv]
            p1cv = 1-pbeta(cv,r1,s1)
            p0cv = 1-pbeta(cv,r0,s0)
            Bcv = h(kappa)*(p1cv-p0cv)*B + p0cv*J
            Acv = matrix(0,nrow=n,ncol=n)
            Acv[tau==1,tau==1] = rbinom( sum(tau==1)^2 , 1, Bcv[1,1])
            Acv[tau==2,tau==2] = rbinom( sum(tau==2)^2 , 1, Bcv[2,2])
            Acv[tau==1,tau==2] = rbinom( sum(tau==1)*sum(tau==2) , 1, Bcv[1,2])
            Acv[tau==2,tau==1] = rbinom( sum(tau==2)*sum(tau==1) , 1, Bcv[2,1])
            Acv = Acv * upper.tri(Acv)
            Acv = Acv+t(Acv)
            if(augmentdiag) for(i in 1:n) Acv[i,i] = sum(Acv[i,])/(n-1)
            S = svd(Acv)
            ## xhat = Acv %*% S$u %*% diag(sqrt(S$d[1:2]))
            xhat = S$u[,1:2] %*% diag(sqrt(S$d[1:2]))
            Lhat[mc,whichcv] = sum(lda(xhat,tau,CV=T)$class != tau)/n
        }
    }
    return(Lhat)
}

error.bayes.mvn <- function(){
    rho <- c(0.5,0.5)
    B <- matrix(c(0.9,0.1,0.1,0.9), nrow = 2, ncol = 2)
    x <- eigen(B)
    xx <- x$vectors %*% diag(sqrt(x$values))

    kappa1 <- 3.5
    kappa2 <- 2
    tau.seq <- seq(from = 0.01, to = 0.99, by = 0.01)

    h.kappa <- (2/kappa1)^2

    n <- 100
    error.bayes.mvn <- numeric(length(tau.seq))

    for(i in 1:length(tau.seq)){
        G1.i <- pbeta(tau.seq[i], kappa1, kappa2, lower.tail = FALSE)
        G0.i <- pbeta(tau.seq[i], kappa2, kappa1, lower.tail = FALSE)
        B.i <- h.kappa*(G1.i*B + G0.i*(1 - B))
        x.i <- eigen(B.i)
        xx.i <- x.i$vectors %*% diag(sqrt(x.i$values))

        tmp <- block.var(xx.i, rho)
        Sigma1 <- tmp[[1]]
        Sigma2 <- tmp[[2]]

       error.bayes.mvn[i] <- bayes.mvn(n, xx.i[1,], xx.i[2,], Sigma1, Sigma2, rho[1], rho[2], 100000)
    }
    return(error.bayes.mvn)
}

block.var <- function(X,rho){
    n <-  nrow(X)
    var.list <- list()

    Delta <-  matrix(0, nrow = ncol(X), ncol = ncol(X))
    for( i in 1:n){
        Delta <- Delta + outer(X[i,],X[i,])*rho[i]
    }

    for(i in 1:n){
        tmp1 <- X[i,]%*%t(X)
        tmp2 <- tmp1 - tmp1^2
        B <- matrix(0, nrow = ncol(X), ncol = ncol(X))
        for(j in 1:n){
            B <-  B + outer(X[j,], X[j,])*tmp2[j]*rho[j]
        }
        var.list[[i]] <- solve(Delta)%*%B%*%solve(Delta)
    }
    return(var.list)
}

bayes.mvn <- function(n, mu1, mu2, Sigma1, Sigma2, pi1, pi2, nmc){

    tau <- sample(1:2, nmc, replace = TRUE, prob = c(pi1,pi2))
    m1 <- sum(tau == 1)
    m2 <-  sum(tau == 2)

    X1 <- mvrnorm(m1, sqrt(n)*mu1, Sigma1)
    X2 <-  mvrnorm(m2, sqrt(n)*mu2, Sigma2)
    labels <- c(rep(1,m1),rep(-1,m2))

    x <- rbind(X1,X2)

    Sigma1.inv <- solve(Sigma1)
    Sigma2.inv <- solve(Sigma2)

    c1 <- log(pi1) - 1/2*log(det(Sigma1))
    c2 <- log(pi2) - 1/2*log(det(Sigma2))

    xbar1 <- x - sqrt(n)*outer(rep(1,nrow(x)),mu1)
    xbar2 <- x - sqrt(n)*outer(rep(1,nrow(x)),mu2)

    density1 <- c1 - rowSums((1/2*xbar1%*%Sigma1.inv)*xbar1)
    density2 <- c2 - rowSums((1/2*xbar2%*%Sigma2.inv)*xbar2)

    labels.hat <- sign(density1 - density2)
    error <- sum(labels.hat != labels)/nmc
    return(error)
}

figure8 <- function(nmc=10000){
    print("NOTE: This will take a few minutes to run. To reproduce an approximate version of the figure reduce nmc")
    require(MASS)
    x=y=NULL
    zzz=NULL 
    n = 50
    kappastar <- 3.5
    tauvec = seq(0,1,by=0.001)
    ntau = length(tauvec)
    SVD = matrix(0,nrow=ntau,ncol=2)

    for(tauindex in 1:ntau){
        tau=tauvec[tauindex]
        Btilde = h(kappastar)*((G1(kappastar,tau)-G0(kappastar,tau))*B + G0(kappastar,tau)*J )
        SVD[tauindex,]=svd(Btilde,2,2)$d
        x[tauindex]=Btilde[1,1]
        y[tauindex]=Btilde[1,2]
        zzz[tauindex] = Lfn(n/2,y[tauindex],x[tauindex])
    }
    plot(tauvec,zzz,type="l",lwd=2,xlab=expression(tau),ylab=expression(L) , ylim = c(0,0.5),
        cex.axis=1.25,cex.lab=2.2)

    B = matrix( c(.9,.1,.1,.9) , nrow=2,ncol=2,byrow=T)
    c1 = B[1,1] + B[1,2]
    c2 = B[1,1]^2 + B[1,2]^2
    G0 = function(kappa,tau) 1-pbeta(tau,2,kappa)
    G1 = function(kappa,tau) 1-pbeta(tau,kappa,2)


    numerator <- function(kappa,tau){
        sqrt(h(kappa))*(B[1,1] - B[1,2])*(G1(kappa,tau) - G0(kappa,tau))
    }
    denominator1 <- function(kappa,tau){
        c1*(G1(kappa,tau) - G0(kappa,tau))*(1 - 2*h(kappa)*G0(kappa,tau))
    }
    denominator2 <- function(kappa,tau){
        2*G0(kappa,tau)*(1 - h(kappa)* G0(kappa,tau))
    }
    denominator3 <- function(kappa,tau){
        c2*h(kappa)*(G1(kappa,tau) - G0(kappa,tau))^2
    }
    denominator <- function(kappa,tau){
        denominator1(kappa,tau) + denominator2(kappa,tau) - denominator3(kappa,tau)
    }
    objective = function(kappa,tau){
        numerator <- sqrt(h(kappa))*(B[1,1] - B[1,2])*(G1(kappa,tau) - G0(kappa,tau))
        denominator <- c1*(G1(kappa,tau) - G0(kappa,tau))*(1 - 2*h(kappa)*G0(kappa,tau)) +
            2*G0(kappa,tau)*(1 - h(kappa)* G0(kappa,tau)) +
            c2*h(kappa)*(G1(kappa,tau) - G0(kappa,tau))^2
        numerator/sqrt(denominator)
    }

    kappavec = seq(2.0,8,by=0.1)
    nkappa = length(kappavec)
    tauvecOPT = seq(0.10,0.90,by=0.001)
    ntau = length(tauvecOPT)
    OBJ = matrix(0,nrow=nkappa,ncol=ntau)
    for(kappaindex in 1:nkappa){
        for(tauindex in 1:ntau){
            kappa=kappavec[kappaindex]
            tau=tauvecOPT[tauindex]
            OBJ[kappaindex,tauindex] = objective(kappa,tau)
        }
    }
    Lapproximation = pnorm( -sqrt(n/2) * OBJ)
    points(tauvecOPT,Lapproximation[which(kappavec==3.5),],
           type="l",lty=2,lwd=2,col="green")

    LSTFP = vertexassignment(nmc=nmc)

    points(seq(0.3,0.9,by=0.025),apply(LSTFP,2,mean),col="red",pch=16)
    points(seq(0.3,0.9,by=0.025),
           apply(LSTFP,2,mean) + apply(LSTFP,2,sd)/sqrt(nmc),col="red",pch=".")
    points(seq(0.3,0.9,by=0.025),
           apply(LSTFP,2,mean) - apply(LSTFP,2,sd)/sqrt(nmc),col="red",pch=".")

    error.bayes = error.bayes.mvn()
    points(seq(from = 0.01, to = 0.99, by = 0.01), error.bayes, col = "blue", pch = "o")
}

figure5 <- function(font=""){
    require(reshape2)
    load("../herm-elegans-connectomes/herm_Graph.Rd")
    Ag[Ag>0]=1


    A <- Ag
    y <- vcols
    p <- order(y)
    y = y[p]
    A = A[p,p]

    size <- as.numeric(table(y))
    sizeCum <- c(0,size[1],size[1]+size[2])
    divLines = c(size[1],size[1]+size[2])

    textdf = data.frame(
        label =  c("Sensory Neurons","Interneurons","Motor Neurons","Sensory Neurons","Interneurons","Motor Neurons"),
        x = c(-10,-10,-10,size/2+sizeCum),y=c(size/2+sizeCum,-10,-10,-10),angle=c(90,90,90,0,0,0))



    g<-ggplot(melt(A),aes(x=Var1,y=Var2,fill=value))+geom_tile(show_guide=FALSE)+
        scale_fill_gradient(low="white", high="black")+
        scale_y_reverse()+
        geom_vline(xintercept=divLines-0.5,alpha=.3,size=1/5)+
        geom_hline(yintercept=divLines-0.5,alpha=.3,size=1/5)+
        geom_text(data=textdf,aes(label=label,x=x,y=y,fill=NULL,angle=angle),family=font, size=4)+
        xlab("")+ylab("")+coord_equal()+
        theme(legend.position="none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            text=element_text(family=font, size=10))

    #ggsave(filename="CElAdj.pdf",plot=g,height=6,width=6)
    return(g)
}


curveLooErr <- function(s, p){

    load("../herm-elegans-connectomes/herm_Graph.Rd")
    Ag[Ag>0]=1


    A <- Ag
    y <- vcols

    curve <- data.frame(pObs=seq(from=0.01,to=1,by=.01))
    curve$pEdg <-  1-s*curve$pObs**p

    nmc <- 1000
    errordf <- data.frame(mc=numeric(),pObs=numeric(),pEdg=numeric(),looErr=numeric())
    for(mc in 1:nmc){
        cat(mc,"\r")
        for(i in 1:nrow(curve)){
            pEdg <- curve$pEdg[i]
            pObs <- curve$pObs[i]
            errordf <- rbind(errordf,data.frame(mc=mc,pObs=pObs,pEdg=pEdg,
                        looErr=looErr(errorful(A,pObs,pEdg),y)))
        }
        save(errordf, file=paste("curveErr_power",p,".Rd",sep=''))
    }
    return(errordf)
}


slope = .2
power = 3

getCurveErr <- function(error, xy, acc=.01){
    nPt = nrow(xy)
    subErr <- error[1,]
    for(i in 1:nPt){
        subErr <- rbind(subErr, 
            error[abs(error$pObs-xy$pObs[i])<acc & abs(error$pEdg-xy$pEdg[i])<acc,][1,])
    }
    return(subErr[2:(nPt+1),])

}

figure6 <- function() {
    print("This will take a really long time, maybe over 24 hours.")
    triangleLooErr()
    trianglePlot()
}

triangleLooErr <- function(){
    load("../herm-elegans-connectomes/herm_Graph.Rd")
    Ag[Ag>0]=1


    A <- Ag
    y <- vcols


    pObsSeq <-  seq(from=.01,to=1,by=.01)
    pEdgSeq <- seq(from=.5,to=1,by=.01)
    nmc <- 1000

    errordf <- data.frame(mc=numeric(),pObs=numeric(),pEdg=numeric(),looErr=numeric())
    cat("\n")
    save(errordf,file="looErr.Rd")
    for(mc in 1:nmc){
        cat(mc,"\r")
        for(pObs in pObsSeq){
            for(pEdg in pEdgSeq){
                if(pEdg >= 1-slope*pObs){
                    errordf <- rbind(errordf,data.frame(mc=mc,pObs=pObs,pEdg=pEdg,
                        looErr=looErr(errorful(A,pObs,pEdg),y)))
                }
            }
        }
        save(errordf,file="looErr.Rd")
    }
}

trianglePlot <- function(slope=2,pow=c(3,5,9)){
    load("./looErr.Rd")
    errorSummary <- ddply(errordf,.(pObs,pEdg),
        function(d)data.frame(m=mean(d$looErr),s=sd(d$looErr)))

    curve <- data.frame()
    curvep <- data.frame(pObs=seq(from=0.01,to=1,by=.01),exponent=0,pEdg=0)
    for(p in pow){
        curvep$pEdg <- 1-slope*curvep$pObs**p
        curvep$exponent = p
        curve <- rbind(curve,curvep)
    }

    curve$exponent <- as.factor(curve$exponent)

    g<-ggplot(errorSummary, aes(x=pObs,y=pEdg,fill=m))+
        geom_tile()+
        scale_fill_gradient(low="red",high="yellow")+
        geom_line(data=curve,aes(fill=NULL,color=exponent,group=exponent))+
        labs(x="edge assessment probability",y="edge classification accuracy",fill="mean error",
            title="")+
        stat_contour(data = errorSummary, aes(x=pObs,y=pEdg,z=m),alpha=.2,bins=20)+
        theme(text=element_text(family="cmr10", size=10)) 

    ggsave("pEdg_vs_pObs_vs_meanError_w_curve.pdf",plot=g,width=7,height=3.5)

    print(g)
    return(g)
}   

figure7 <- function(slope = .2) {
    print("This will take a long time.")
    curveLooErr(slope,3)
    curveLooErr(slope,5)
    curveLooErr(slope,9)
    curvePlot()
}

curvePlot <- function(){
    require(plyr)
    load("../herm-elegans-connectomes/herm_Graph.Rd")
    
    # Chance 
    size <- as.numeric(table(vcols))
    chance= 1-max(size/sum(size))

    load("./curveErr_power3.Rd")
    nmc <- max(errordf$mc)
    errorSummary <- ddply(errordf,.(pObs,pEdg),function(d)
        data.frame(looErr=mean(d$looErr),
                   sd=sd(d$looErr),
                   se=sd(d$looErr)/sqrt(nmc),
                   exponent=3))

    load("./curveErr_power5.Rd")
    nmc <- max(errordf$mc)
    errorSummary <- rbind(errorSummary,
        ddply(errordf,.(pObs,pEdg),function(d)
            data.frame(looErr=mean(d$looErr),
                       sd=sd(d$looErr),
                       se=sd(d$looErr)/sqrt(nmc),
                       exponent=5)))

    load("./curveErr_power9.Rd")
    nmc <- max(errordf$mc)
    errorSummary <- rbind(errorSummary,
        ddply(errordf,.(pObs,pEdg),function(d)
            data.frame(looErr=mean(d$looErr),
                       sd=sd(d$looErr),
                       se=sd(d$looErr)/sqrt(nmc),
                       exponent=9)))

    errorSummary$exponent <- as.factor(errorSummary$exponent)

    g<-ggplot(errorSummary,aes(x=pObs,y=looErr,ymax=looErr+2*se,ymin=looErr-2*se,fill=exponent,shape=exponent))+
        geom_line()+geom_point()+
        geom_ribbon(alpha=.3)+
        geom_abline(intercept=chance,slope=0)+
        labs(x="edge assessment probability",y="mean error",
            title="")+
        theme(text=element_text(family="cmr10", size=10)) 

    ggsave("pObs_vs_meanError_curve.pdf",plot=g,width=7,height=3.5)
    return(g)
}
