ifa.em.alg<-function(y,numobs,L,numvar,ni,totni,maxni,it,H,w,mu,vu,eps,psi,lik)
{                 
likelihood<-NULL
hh<-0
ratio<-1000
H.dep<-NULL

while ((hh < it) & (ratio > eps )) {
 hh<-hh+1 
    
wq<-matrix(1,1,totni)
muq<-matrix(0,L,totni)
vuq<-array(0,c(L,L,totni))

pr<-totni
cont<-1
for (k in 1:L) {j<-1
                  pr<-pr/ni[k]
                  for (i in 1:totni) {muq[k,i]<-mu[k,j]
                                         vuq[k,k,i]<-vu[k,j]
                                         wq[,i]<-wq[,i]*w[k,j]
                                         cont<-cont+1
                                         if (cont>pr) {if (j==ni[k]) j<-1 else (j<-j+1) 
                                                        cont<-1}
                                        }}

sigma<-array(0,c(L,L,totni))
spsi<-solve(psi)
roqy<-array(0,c(L,numobs,totni))
            
sigma<-array(apply((array(t(H) %*% spsi %*% H,c(L,L,totni))+array(apply(vuq,3,solve),c(L,L,totni))),3,solve),c(L,L,totni))
            
            
for (i in 1:totni) {
roqy[,,i]<-sigma[,,i]%*%((t(H)%*%spsi%*%t(y)) + matrix(solve(vuq[,,i])%*%muq[,i],L,numobs))}
    
Exqy<-roqy
Exxqy<-array(0,c(L,L,numobs,totni))
pyq<-matrix(0,numobs,totni)

for (i in 1:numobs) {for (j in 1 :totni) {Exxqy[,,i,j]<-sigma[,,j]+(roqy[,i,j]%*%t(roqy[,i,j]))
                                            pyq[i,j]<-dmvnorm(t(y[i,]),t(H%*%muq[,j]),(H%*%vuq[,,j]%*%t(H)+psi)) }}

pqy<-matrix(0,numobs,totni)
den<-(wq%*%t(pyq))
temp<-t(matrix(wq,totni,numobs))
pqy<-pyq*temp/den[,]
pqy<-ifelse(is.na(pqy),mean(pqy, na.rm=TRUE),pqy)
pqiy<-array(0,c(numobs,L,maxni))
nummu<-array(0,c(numobs,L,maxni))
numvu<-array(0,c(numobs,L,maxni))

pr<-totni
cont<-1
for (k in 1:L) {j<-1
                  pr<-pr/ni[k]
                  for (i in 1:totni) {nummu[,k,j]<-nummu[,k,j]+(pqy[,i]*Exqy[k,,i])
                                         numvu[,k,j]<-numvu[,k,j]+(pqy[,i]*Exxqy[k,k,,i])
                                         pqiy[,k,j]<-pqiy[,k,j]+pqy[,i]
                                         cont<-cont+1
                                         if (cont>pr) {if (j==ni[k]) j<-1 else (j<-j+1) 
                                                        cont<-1}
                                        }}

Exy<-matrix(0,L,numobs)
Exxy<-array(0,c(L,L,numobs))

for (j in 1:totni) {Exy<-Exy+(t(matrix(pqy[,j],numobs,L))*Exqy[,,j])
                    Exxy<-Exxy+(array(rep(pqy[,j],each=L*L),c(L,L,numobs))*Exxqy[,,,j])}

Exy<-ifelse(is.na(Exy),mean(Exy, na.rm=TRUE),Exy)
Exxy<-ifelse(is.na(Exxy),mean(Exxy, na.rm=TRUE),Exxy)

EExxy<-rowMeans(Exxy, na.rm=TRUE, dims=2)

H<-matrix(0,numvar,L)
sEExxy<-solve(EExxy)
psi<-matrix(0,numvar,numvar)

for (i in 1:numobs){H <- H + (t(t(y[i,]))%*%t(Exy[,i])%*%sEExxy)}
H<-H/numobs

for (i in 1:numobs) {psi<-psi+ (t(t(y[i,]))%*%t(y[i,])-(t(t(y[i,]))%*%t(Exy[,i])%*%t(H)))}
psi<-psi/numobs

Enummu<-matrix(0,L,maxni)
Enumvu<-matrix(0,L,maxni)
Eden<-matrix(0,L,maxni)
mu<-matrix(0,L,maxni)
vu<-matrix(0,L,maxni)
w<-matrix(0,L,maxni)

for (i in 1:L) for (j in 1:ni[i]) 
        {Enummu[i,j]<-sum(nummu[,i,j]) 
        Eden[i,j]<-sum(pqiy[,i,j])
        Enumvu[i,j]<-sum(numvu[,i,j])}
        
for (i in 1:L) for (j in 1:ni[i]) 
        {mu[i,j]<-Enummu[i,j]/Eden[i,j]
         vu[i,j]<-Enumvu[i,j]/Eden[i,j]-(mu[i,j]*mu[i,j])
         w[i,j]<-mean(pqiy[,i,j])}


                
######## scaling ##########

sigmascale<-matrix(0,L)

for (i in 1:L) {
           cont1<-0
           cont2<-0
    for (j in 1:ni[i]) {
        cont1<-cont1+ ( w[i,j]*(mu[i,j]*mu[i,j]+vu[i,j]) )
        cont2<-cont2+(w[i,j]*mu[i,j])    }
                
    sigmascale[i]<-cont1-(cont2*cont2)
     }


for (i in 1:L) {
    for (j in 1:ni[i]) {
        mu[i,j]<-mu[i,j]/sqrt(sigmascale[i]) 
        vu[i,j]<-vu[i,j]/sigmascale[i]
                      }
    for (j in 1:numvar)  H[j,i]=H[j,i]*sqrt(sigmascale[i])
   }

temp<-sum(log(pyq%*%t(wq)))
likelihood<-c(likelihood,temp)
ratio<-abs((temp-lik)/lik)
if ((temp < lik) & (hh > 5)) ratio<-eps
lik<-temp

                                    }

####### asymptotic standard error

niter<-length(likelihood)
temp<-(likelihood-likelihood[niter])
r<-log(temp[-1]/temp[-niter])
r<-exp(mean(r[1:(niter-2)],na.rm=T))

se.H<-matrix(0,numvar,L)
for (i in 1:L) se.H[,i]<-sqrt(diag(psi)/(numobs*EExxy[i,i]*(1-r)))
se.psi<-diag(sqrt(2*diag(psi^2)/(numobs*(1-r))))
se.mu<- sqrt(vu/(w*numobs*(1-r)))
se.vu<- sqrt(vu^2/(w*numobs*(1-r)))
se.w<-sqrt(1/((1-w)*w*numobs*(1-r)))
std.err<-list(H=se.H,psi=se.psi,mu=se.mu,vu=se.vu,w=se.w)
                                    
out<-list(H=H,w=w,mu=mu,vu=vu,psi=psi,likelihood=likelihood,sigma=sigma,pqy=pqy,lik=lik,std.err=std.err)
return(out)
}
