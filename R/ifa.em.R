"ifa.em" <-function(y,ni,it=15,eps=0.001,init=NULL,scaling=TRUE)
{
if (scaling) y<-scale(y)
ni<-matrix(ni)
L<-nrow(ni)
numobs<-nrow(y)
numvar<-ncol(y)
ybar<- apply(y, 2., mean)
y<-scale(y, ybar, scale=FALSE) 
lik<--100000000000
totni<-prod(ni)
maxni<-max(ni)

if (is.null(init)) output.init<-ifa.init.pca(y,L) else output.init<-init

psi<-output.init$psi
psi<-diag(diag(psi))
H<-output.init$H


w<-matrix(0,L,maxni)
mu<-matrix(0,L,maxni)
vu<-matrix(0,L,maxni)

for (i in 1:L) for (j in 1:ni[i]) {w[i,j]<-runif(1,0,1)
                                        mu[i,j]<-runif(1,-1,1)
                                        vu[i,j]<-runif(1,0,1)} 
w<-w/rowSums(w)


out<-ifa.em.alg(y,numobs,L,numvar,ni,totni,maxni,it,H,w,mu,vu,eps,psi,lik)

likelihood<-out$likelihood
sigma<-out$sigma
pqy<-out$pqy
H<-out$H
w<-out$w
mu<-out$mu
vu<-out$vu
psi<-out$psi
psi<-diag(diag(psi))
se<-out$std.err

output<-list(H=H,lik=likelihood,w=w,mu=mu,vu=vu,psi=psi,totni=totni,ni=ni,L=L,sigma=sigma,pqy=pqy,numvar=numvar,numobs=numobs,scaling=scaling,se=se)
invisible(output)
}
