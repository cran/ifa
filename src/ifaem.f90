SUBROUTINE ifaem (y,numobs,L,numvar,ni,totni,maxni,it,H,w,mu,vu,eps,delta,likelihood,sigma,pqy,EExxy)
!DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:'ifaem_' :: ifaem

IMPLICIT NONE



INTEGER, INTENT(IN) :: numobs,numvar,L,totni,maxni,it
INTEGER, INTENT(IN), DIMENSION(L,1)::ni
DOUBLE PRECISION, INTENT(IN), DIMENSION(numvar,numobs) :: y
DOUBLE PRECISION, INTENT(IN):: eps
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(numvar,numvar) :: delta
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(L,maxni) :: vu,mu,w
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(numvar,L) :: H
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(it) :: likelihood
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(L,L,totni) :: sigma
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(numobs,totni) :: pqy
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(L,L) :: EExxy

DOUBLE PRECISION, ALLOCATABLE :: wq(:,:)
DOUBLE PRECISION, ALLOCATABLE:: pyq(:,:)
DOUBLE PRECISION, ALLOCATABLE:: pqiy(:,:,:)
DOUBLE PRECISION, ALLOCATABLE :: Exy(:,:)
DOUBLE PRECISION, ALLOCATABLE :: Exxy(:,:,:)
DOUBLE PRECISION, ALLOCATABLE::  nummu(:,:,:), numvu(:,:,:)
DOUBLE PRECISION, ALLOCATABLE:: muq(:,:)
DOUBLE PRECISION, ALLOCATABLE:: vuq(:,:,:)
DOUBLE PRECISION, ALLOCATABLE :: sdelta(:,:), var(:,:)


DOUBLE PRECISION, ALLOCATABLE :: svuq(:,:),ssigma(:,:)
DOUBLE PRECISION, ALLOCATABLE::roqy(:,:,:)
DOUBLE PRECISION, ALLOCATABLE::Exqy(:,:,:)
DOUBLE PRECISION, ALLOCATABLE::Exxqy(:,:,:,:)
DOUBLE PRECISION, ALLOCATABLE::matres(:,:)
DOUBLE PRECISION, ALLOCATABLE :: den(:,:)
DOUBLE PRECISION, ALLOCATABLE:: temp3(:,:)
DOUBLE PRECISION, ALLOCATABLE :: temp4(:,:,:)
DOUBLE PRECISION, ALLOCATABLE :: temp5(:,:,:)
DOUBLE PRECISION, ALLOCATABLE:: valore(:,:),media(:,:)
DOUBLE PRECISION, ALLOCATABLE::sigmascale(:)
DOUBLE PRECISION::prob,cont1,cont2,ratio,lik
DOUBLE PRECISION, ALLOCATABLE :: sEExxy(:,:)
DOUBLE PRECISION, ALLOCATABLE :: Enummu(:,:), Eden(:,:), Enumvu(:,:)
INTEGER::i,j,k,pr,cont,hh

hh=0
ratio=1000
lik=-1000000


do while ((hh .LT. it) .AND. (ratio .GT. eps ))
hh=hh+1
allocate (vuq(L,L,totni))
allocate (muq(L,totni))
allocate (wq(1,totni))

wq=1
muq=0
vuq=0
pr=totni
 cont=1


DO k=1,L,1
		j=1
		pr=pr/ni(k,1)
		DO i=1,totni,1
			muq(k,i)=mu(k,j)
			vuq(k,k,i)=vu(k,j)
			wq(:,i)=wq(:,i)*w(k,j)
			cont=cont+1
			if (cont .GT. pr) then
							if (j .EQ. ni(k,1)) then
							j=1
							else
							j=j+1
							end if
							cont=1
			end if
		END DO
END DO



allocate (sdelta(numvar,numvar))

 call inversa(delta,sdelta,numvar)


allocate (roqy(L,numobs,totni))
roqy=0

allocate (matres(L,numobs))
allocate (svuq(L,L))
allocate (ssigma(L,L))

svuq=0



DO i=1,totni,1
        do j=1,L,1
		svuq(j,j)=1/vuq(j,j,i)
		enddo
		call inversa(svuq+matmul(matmul(transpose(H),sdelta),H),ssigma,L)
		sigma(:,:,i)=ssigma
		matres=reshape(matmul(svuq,muq(:,i)),(/L,numobs/),matmul(svuq,muq(:,i)))
		roqy(:,:,i)=matmul(sigma(:,:,i),matres+matmul(matmul(transpose(H),sdelta),y))
END DO

deallocate(svuq)
deallocate(sdelta)
deallocate(ssigma)

allocate (Exqy(L,numobs,totni))
allocate (Exxqy(L,L,numobs,totni))
allocate (pyq(numobs,totni))
allocate (var(numvar,numvar))
allocate(valore(numvar,1))
allocate(media(numvar,1))

Exqy=roqy
Exxqy=0
pyq=0


DO i=1,numobs,1
		DO j=1,totni,1
			Exxqy(:,:,i,j)=sigma(:,:,j) + &
			matmul(reshape(roqy(:,i,j),(/L,1/)),transpose(reshape(roqy(:,i,j),(/L,1/))))
			var=matmul(matmul(H,vuq(:,:,j)),transpose(H))+delta
			valore=reshape(y(:,i),(/numvar,1/))
			media=matmul(H,reshape(muq(:,j),(/L,1/)))
			media=reshape(media,(/numvar,1/))
			call mnorm(valore,numvar,media,var,prob)
			if (prob .EQ. 0) then
							pyq(i,j)=0.0000000001
							else
							pyq(i,j)=prob
							end if
		END DO
END DO

deallocate(media)
deallocate(valore)
deallocate(var)
deallocate(muq)
deallocate(vuq)
deallocate (roqy)
deallocate (matres)
allocate (den(1,numobs))
pqy=0
den=matmul(wq,transpose(pyq))


DO j=1,totni,1
DO i=1,numobs,1
			pqy(i,j)=log(pyq(i,j))+log(wq(1,j))-log(den(1,i))
			pqy(i,j)=exp(pqy(i,j))
		END DO
END DO

deallocate(den)
allocate (pqiy(numobs,L,maxni))

pqiy=0
allocate (nummu(numobs,L,maxni))
allocate (numvu(numobs,L,maxni))
nummu=0
numvu=0

pr=totni
 cont=1
DO k=1,L,1
		j=1
		pr=pr/ni(k,1)
		DO i=1,totni,1
			nummu(:,k,j)=nummu(:,k,j)+(pqy(:,i)*Exqy(k,:,i))
			numvu(:,k,j)=numvu(:,k,j)+(pqy(:,i)*Exxqy(k,k,:,i))
			pqiy(:,k,j)=pqiy(:,k,j)+pqy(:,i)
			cont=cont+1
			if (cont .GT. pr) then
							if (j .EQ. ni(k,1)) then
							j=1
							else
							j=j+1
							end if
							cont=1
			end if
		END DO
END DO

allocate (Exy(L,numobs))
allocate (Exxy(L,L,numobs))
allocate (temp3(L,numobs))
allocate (temp4(numobs,L,L))
allocate (temp5(L,L,numobs))

Exy=0
Exxy=0
DO j=1,totni,1
	temp3=transpose(reshape(pqy(:,j),(/numobs,L/),pqy(:,j),ORDER=(/1,2/)))
	Exy=Exy + (temp3*Exqy(:,:,j))
	temp4= reshape(pqy(:,j),(/numobs,L,L/),pqy(:,j),ORDER=(/1,2,3/))
	temp5=reshape(temp4,(/L,L,numobs/),temp4,ORDER=(/3,2,1/))
	Exxy=Exxy + (temp5*Exxqy(:,:,:,j))
END DO

deallocate(temp4)
deallocate(temp3)
deallocate(temp5)
deallocate (Exqy)
deallocate (Exxqy)


EExxy=SUM(Exxy,DIM=3)
EExxy=EExxy/numobs

deallocate(Exxy)

H=0
delta=0

allocate(sEExxy(L,L))

  call inversa(EExxy,sEExxy,L)


H=MATMUL(MATMUL(y,TRANSPOSE(Exy)),sEExxy)
H=H/numobs

deallocate(sEExxy)

delta=(matmul(y,transpose(y))-matmul(matmul(y,transpose(Exy)),transpose(H)))
delta=delta/numobs

deallocate(Exy)
allocate(Enummu(L,maxni))
allocate(Eden(L,maxni))
allocate(Enumvu(L,maxni))

Enummu=0
Enumvu=0
Eden=0
mu=0
vu=0
w=0

DO i=1,L,1
	DO j=1,ni(i,1),1
		Enummu(i,j)=SUM(nummu(:,i,j))
		Eden(i,j)=SUM(pqiy(:,i,j))
		Enumvu(i,j)=SUM(numvu(:,i,j))
	END DO
END DO

deallocate(nummu)
deallocate(numvu)


DO i=1,L,1
	DO j=1,ni(i,1),1
		mu(i,j)=Enummu(i,j)/Eden(i,j)
		vu(i,j)=Enumvu(i,j)/Eden(i,j)-(mu(i,j)*mu(i,j))
		w(i,j)=SUM(pqiy(:,i,j))/numobs
	END DO
END DO

deallocate(Eden)
deallocate(Enummu)
deallocate(Enumvu)
deallocate(pqiy)

allocate(sigmascale(L))

DO i=1,L,1
    cont1=0
	cont2=0
	DO j=1,ni(i,1),1
			cont1=cont1+(w(i,j)*(mu(i,j)*mu(i,j)+vu(i,j)))
			cont2=cont2+(w(i,j)*mu(i,j))
	END DO
	sigmascale(i)=cont1-(cont2*cont2)

END DO

DO i=1,L,1
	DO j=1,ni(i,1),1
		mu(i,j)=mu(i,j)/SQRT(sigmascale(i))
		vu(i,j)=vu(i,j)/sigmascale(i)
	END DO
	DO j=1,numvar,1
	H(j,i)=H(j,i)*SQRT(sigmascale(i))
	END DO
END DO

deallocate(sigmascale)


likelihood(hh)=sum(log(matmul(pyq,transpose(wq))))
ratio=abs((likelihood(hh)-lik)/lik)
lik=likelihood(hh)

deallocate(wq)
deallocate(pyq)
end do

END SUBROUTINE ifaem

SUBROUTINE MNORM(x,k,m,v,prob)


IMPLICIT NONE

INTEGER, INTENT(IN) :: k
DOUBLE PRECISION, INTENT(IN), DIMENSION(k,1) :: x,m
DOUBLE PRECISION, INTENT(IN), DIMENSION(k,k) :: v

DOUBLE PRECISION, INTENT(OUT):: prob

double precision, DIMENSION(k,k):: s_v

double precision:: det_v
double precision::pi=3.14159265358979
double precision,dimension(1,1)::temp
double precision::temp1

 call determin(v,k,det_v)
 call inversa(v,s_v,k)

temp=-(0.5)*(matmul((matmul(transpose(x-m),s_v)),x-m))
temp1=temp(1,1)

prob=((sqrt((2*pi)**k)*sqrt(det_v))**(-1))*exp(temp1)

	
END SUBROUTINE



 Subroutine determin(Ain,dim,determ)

!mediante la scomposizione di Cholesky!

INTEGER, INTENT(IN):: dim
DOUBLE PRECISION,INTENT(IN), DIMENSION(dim,dim):: Ain
DOUBLE PRECISION, INTENT(OUT):: determ
DOUBLE PRECISION, DIMENSION(dim):: P
INTEGER:: i,j,k
DOUBLE PRECISION::temp,s
DOUBLE PRECISION, DIMENSION(dim,dim)::A

A=Ain
do i=1,dim,1
	do j=i,dim,1
		s=A(i,j)
		temp=0
		if (i .GT. 1) then
						do k=1,i-1,1
						temp=temp+(A(i,k)*A(j,k))
						end do
		end if
		s=s-temp
		if (i .EQ. j) then 
						
						P(i)=sqrt(s)
						
					  else
						A(j,i)=s/(P(i))
		end if

	end do
end do


determ=1
do i=1,dim,1
	determ=determ*P(i)*P(i)
end do
end subroutine



Subroutine INVERSA ( Ain,Aout, NP )

      !-------------------------------------------------------------------------
      !
      !	      Taken from "Numeric recipes".  The original program was
      !       GAUSSJ which solves linear equations by the Gauss_Jordon
      !       elimination method.  Only the parts required to invert
      !	      matrices have been retained.
      !
      !	      J.P. Griffith  6/88
      !
      !-------------------------------------------------------------------------


      PARAMETER (NMAX=50)

      
	  	DOUBLE PRECISION, DIMENSION(NP,NP):: Ain,Aout
		DIMENSION IPIV(NMAX), INDXR(NMAX), INDXC(NMAX)

      n = np

	  Aout=Ain

      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(Aout(J,K)).GE.BIG)THEN
                  BIG=ABS(Aout(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE

        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=Aout(IROW,L)
            Aout(IROW,L)=Aout(ICOL,L)
            Aout(ICOL,L)=DUM
14        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        PIVINV=1./Aout(ICOL,ICOL)
        Aout(ICOL,ICOL)=1.
        DO 16 L=1,N
          Aout(ICOL,L)=Aout(ICOL,L)*PIVINV
16      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=Aout(LL,ICOL)
            Aout(LL,ICOL)=0.
            DO 18 L=1,N
              Aout(LL,L)=Aout(LL,L)-Aout(ICOL,L)*DUM
18          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=Aout(K,INDXR(L))
            Aout(K,INDXR(L))=Aout(K,INDXC(L))
            Aout(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END


