"ifa.init.random" <-
function(y,L)
{
numvar<-ncol(y)
H<-matrix(runif(L*numvar,-1,1),numvar,L)
psi<-var(y)*0.1

output<-list(psi=psi,H=H)
}
