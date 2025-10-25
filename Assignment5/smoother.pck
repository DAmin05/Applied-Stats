smoother.pck <-
c("smoother.pck", "bin.mean", "gauss.mean", "gauss.reg", "gauss.mean.trunc", 
"gauss.reg.trunc", "my.hat.w")
bin.mean <-
function(x,y,nbin,xcol=2)
{
o1<-order(x)
x1<-x[o1]
y1<-y[o1]
r1<-range(x)
inc<-(r1[2]-r1[1])/nbin
yvec<-NULL
smat<-NULL
for(i in 1:nbin){
        bin.low<-r1[1]+(i-1)*inc
        bin.high<-r1[1]+i*inc
        
        I1<-x1>=bin.low
if(i<nbin){
        I2<-x1<bin.high
}else{
        I2<-x1<=(bin.high+200)
}
        I3<-as.logical(I1*I2)
        yval<-mean(y1[I3])
        n1<-sum(I3)
        matdum<-NULL
        for(i in 1:n1){
        matdum<-rbind(matdum,I3*1/n1)
        }
        smat<-rbind(smat,matdum)
        yvec<-c(yvec,rep(yval,n1))
}
n99<-length(x1)
dferror<-length(x1)-sum(diag(2*smat-smat%*%(t(smat))))
delta1<-sum(diag(t(diag(n99)-smat)%*%(diag(n99)-smat)))
R<-t(diag(n99)-smat)%*%(diag(n99)-smat)
delta2<-2*sum(diag(R%*%R))
lines(x1,yvec,col=xcol)
ypred<-y1
ypred<-smat%*%y1
resid<-y-ypred

list(smat=smat,df=sum(diag(smat)),dferror=dferror,delta1=delta1,delta2=delta2,resid=resid,pred=ypred,x=x)
                
}
gauss.mean <-
function(x,y,lambda,xcol=3,do.plot=T)
{
o1<-order(x)
x1<-x[o1]
y1<-y[o1]
r1<-range(x)
smat<-NULL
n1<-length(x1)
for(i in 1:n1){
        v1<-dnorm(x1,x1[i],lambda)
        v1<-v1/sum(v1)
        smat<-rbind(smat,v1)
}
yhat<-smat%*%y1
if(do.plot){
lines(x1,yhat,col=xcol)
}
n99<-length(x1)
dferror<-length(x1)-sum(diag(2*smat-smat%*%(t(smat))))
delta1<-sum(diag(t(diag(n99)-smat)%*%(diag(n99)-smat)))
R<-t(diag(n99)-smat)%*%(diag(n99)-smat)
delta2<-2*sum(diag(R%*%R))
resid<-y1-smat%*%y1
ypred<-y1
ypred[o1]<-smat%*%y1
PRESS<-sum((resid/(1-diag(smat)))^2)
list(smat=smat,df=sum(diag(smat)),dferror=dferror,delta1=delta1,delta2=delta2,resid=resid,pred=ypred,press=PRESS)
        
}
gauss.reg <-
function(x,y,lambda,xcol=4,do.plot=T)
{
o1<-order(x)
x1<-x[o1]
y1<-y[o1]
r1<-range(x)
smat<-NULL
n1<-length(x1)
for(i in 1:n1){
        v1<-dnorm(x1,x1[i],lambda)
        v1<-v1/sum(v1)
        H1<-my.hat.w(x1,v1)
        smat<-rbind(smat,H1[i,])
}
yhat<-smat%*%y1
if(do.plot){
lines(x1,yhat,col=xcol)
}
n99<-length(x1)
dferror<-length(x1)-sum(diag(2*smat-smat%*%(t(smat))))
delta1<-sum(diag(t(diag(n99)-smat)%*%(diag(n99)-smat)))
R<-t(diag(n99)-smat)%*%(diag(n99)-smat)
delta2<-2*sum(diag(R%*%R))
resid<-y1-smat%*%y1
ypred<-y1
ypred[o1]<-smat%*%y1
list(smat=smat,df=sum(diag(smat)),dferror=dferror,delta1=delta1,delta2=delta2,resid=resid,pred=ypred)
}
gauss.mean.trunc <-
function(x,y,lambda,nnn,xcol=5,do.plot=T)
{
o1<-order(x)
x1<-x[o1]
y1<-y[o1]
r1<-range(x)
smat<-NULL
n1<-length(x1)
trunc.val<-n1-nnn
for(i in 1:n1){
        v1<-dnorm(x1,x1[i],lambda)
        o2<-order(v1)
        thresh<-v1[o2[trunc.val]]
        v1<-v1*(v1>thresh)
        v1<-v1/sum(v1)
        smat<-rbind(smat,v1)
}
yhat<-smat%*%y1
if(do.plot){
lines(x1,yhat,col=xcol)
}
n99<-length(x1)
dferror<-length(x1)-sum(diag(2*smat-smat%*%(t(smat))))
delta1<-sum(diag(t(diag(n99)-smat)%*%(diag(n99)-smat)))
R<-t(diag(n99)-smat)%*%(diag(n99)-smat)
delta2<-2*sum(diag(R%*%R))
resid<-y1-smat%*%y1
ypred<-y1
ypred[o1]<-smat%*%y1
list(smat=smat,df=sum(diag(smat)),dferror=dferror,delta1=delta1,delta2=delta2,resid=resid,pred=ypred)
                
}
gauss.reg.trunc <-
function(x,y,lambda,nnn,xcol=6,do.plot=T)
{
o1<-order(x)
x1<-x[o1]
y1<-y[o1]
r1<-range(x)
smat<-NULL
n1<-length(x1)
trunc.val<-n1-nnn
for(i in 1:n1){
        v1<-dnorm(x1,x1[i],lambda)
        o1<-order(v1)
        thresh<-v1[o1[trunc.val]]
        v1<-v1*(v1>thresh)
        v1<-v1/sum(v1)
        H1<-my.hat.w(x1,v1)
        smat<-rbind(smat,H1[i,])
}
yhat<-smat%*%y1
if(do.plot){
lines(x1,yhat,col=xcol)
}
n99<-length(x1)
dferror<-length(x1)-sum(diag(2*smat-smat%*%(t(smat))))
delta1<-sum(diag(t(diag(n99)-smat)%*%(diag(n99)-smat)))
R<-t(diag(n99)-smat)%*%(diag(n99)-smat)
delta2<-2*sum(diag(R%*%R))
resid<-y1-smat%*%y1
ypred<-y1
ypred[o1]<-smat%*%y1
list(smat=smat,df=sum(diag(smat)),dferror=dferror,delta1=delta1,delta2=delta2,resid=resid,pred=ypred)

        
}
my.hat.w <-
function(x,wt){
x1<-cbind(1,x)
x1%*%solve(t(x1)%*%diag(wt)%*%x1)%*%t(x1)%*%(diag(wt))
}
