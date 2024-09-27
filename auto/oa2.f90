!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)


      DOUBLE PRECISION Ra,Rb,theta,Delta,K,p,alpha
      INTEGER i,j

       Ra=U(1)
       Rb=U(2)
       theta=U(3)
        
       Delta=1.0
       K=PAR(1)
       p=PAR(2)
       alpha=PAR(3)

       F(1)= -Delta*Ra + 0.5*K*(1-Ra**2)*(p*Ra*cos(alpha)+(1-p)*Rb*cos(theta-alpha))
       F(2)= -Delta*Rb + 0.5*K*(1-Rb**2)*(p*Rb*cos(alpha)+(1-p)*Ra*cos(theta+alpha))
       F(3)=  0.5*K*(1+Ra**2)/(Ra)*(p*Ra*sin(alpha) -(1-p)*Rb*sin(theta-alpha)) &
- 0.5*K*(1+Rb**2)/(Rb)*(p*Rb*sin(alpha) +(1-p)*Ra*sin(theta+alpha))

      IF(IJAC.EQ.0)RETURN

        do i = 1, 3
                do j = 1, 3
                        DFDU(i,j)=0;
                end do
        end do

       DFDU(1,1) = -Delta-Ra*K*(p*Ra*cos(alpha)+(1-p)*Rb*cos(theta-alpha)) & 
+ 0.5*(1-Ra**2)*p*K*cos(alpha);
       DFDU(1,2) = 0.5*(1-p)*K*cos(theta-alpha)*(1-Ra**2);
       DFDU(1,3) = -0.5*(1-Ra**2)*K*(1-p)*Rb*sin(theta-alpha);

       DFDU(2,1) = 0.5*(1-p)*K*cos(theta+alpha)*(1-Rb**2);
       DFDU(2,2) = -Delta-Rb*K*(p*Rb*cos(alpha)+(1-p)*Ra*cos(theta+alpha)) &
+ 0.5*(1-Rb**2)*p*K*cos(alpha)
       DFDU(2,3) = -0.5*(1-Rb**2)*K*(1-p)*Ra*sin(theta+alpha)

       DFDU(3,1) = p*K*sin(alpha)*Ra - (1-p)*K*( (Ra**2-1)/(2*Ra**2)*Rb*sin(theta-alpha)+(1+Rb**2)/(2.0*Rb)*sin(theta+alpha) )
       DFDU(3,2) =-p*K*sin(alpha)*Rb - (1-p)*K*( (1+Ra**2)/(2.0*Ra)*sin(theta-alpha) + (Rb**2-1)/(2*Rb**2)*Ra*sin(theta+alpha) )
       DFDU(3,3) = -(1-p)*K*((1+Ra**2)/(2.0*Ra)*Rb*cos(theta-alpha) &
 + (1+Rb**2)/(2*Rb)*Ra*cos(theta+alpha))


    IF(IJAC.EQ.1)RETURN

      do i = 1, 3
              do j = 1, 3
                      DFDP(i,j)=0;
              end do
      end do

      ! Rows are equations, cols are parameters K,p,alpha
      DFDP(1,1)= 0.5*(1-Ra**2)*(p*Ra*cos(alpha)+(1-p)*Rb*cos(theta-alpha))
      DFDP(1,2)= 0.5*K*(1-Ra**2)*(Ra*cos(alpha)-Rb*cos(theta-alpha))
      DFDP(1,3)= 0.5*K*(1-Ra**2)*(-p*Ra*sin(alpha)+(1-p)*Rb*sin(theta-alpha))

      DFDP(2,1)= 0.5*(1-Rb**2)*(p*Rb*cos(alpha)+(1-p)*Ra*cos(theta+alpha))
      DFDP(2,2)= 0.5*K*(1-Rb**2)*(Rb*cos(alpha)-Ra*cos(theta+alpha))
      DFDP(2,3)= 0.5*K*(1-Rb**2)*(-p*Rb*sin(alpha)-(1-p)*Ra*sin(theta+alpha))

      DFDP(3,1)= 0.5*(1+Ra**2)/(Ra)*(p*Ra*sin(alpha) -(1-p)*Rb*sin(theta-alpha)) &
 -0.5*(1+Rb**2)/(Rb)*(p*Rb*sin(alpha) +(1-p)*Ra*sin(theta+alpha))
      DFDP(3,2)= K*0.5*(1+Ra**2)/(Ra)*(Ra*sin(alpha) +Rb*sin(theta-alpha)) &
 -K*0.5*(1+Rb**2)/(Rb)*(Rb*sin(alpha) -Ra*sin(theta+alpha))
      DFDP(3,3)= K*0.5*(1+Ra**2)/(Ra)*(p*Ra*cos(alpha) +(1-p)*Rb*cos(theta-alpha)) &
 -K*0.5*(1+Rb**2)/(Rb)*(p*Rb*cos(alpha) +(1-p)*Ra*cos(theta+alpha))


      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      ! COMMON block needed if IPS=9 (homoclinic bifurcations) :
 !     INTEGER ITWIST,ISTART,IEQUIB,NFIXED,NPSI,NUNSTAB,NSTAB,NREV
 !     COMMON /BLHOM/ ITWIST,ISTART,IEQUIB,NFIXED,NPSI,NUNSTAB,NSTAB,NREV


      !PAR(1)=6.0 ! K
      !PAR(2)=0.9 ! p
      !PAR(3)=1.2 ! alpha

      ! Antiphase
      PAR(1)=10.0 ! K
      PAR(2)=0.9 ! p
      PAR(3)=1.2 ! alpha

      !U(1)=0.58
      !U(2)=0.57
      !U(3)=3.14
      U(1)=0.2   
      U(2)=0.22
      U(3)=1.6336281798666925


!      ! If IEQUIB=1 then put the equilibrium in PAR(11+i), i=1,...,NDIM :
!
!        IF (IEQUIB.NE.0) THEN
!          PAR(12) = 0.72945103
!          PAR(13) = 0.71145556
!          PAR(14) = 0.0770218
!          PAR(15) = 0.71145556
!          PAR(16) = 0.72945103
!          PAR(17) = 0.0770218
!        ENDIF



      END SUBROUTINE STPNT

      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT


!-----------------------------------------------------------
      DOUBLE PRECISION FUNCTION GETUY_MIN2(U,NDX,NTST,NCOL)
      INTEGER, INTENT(IN) :: NDX,NCOL,NTST
      DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
      DOUBLE PRECISION  :: PCLU_Y(0:NTST*NCOL)

      PCLU_Y=U(2,:)

      GETUY_MIN2=MINVAL(PCLU_Y)

      END FUNCTION GETUY_MIN2
!-----------------------------------------------------------
!-----------------------------------------------------------
      DOUBLE PRECISION FUNCTION GETUY_MIN(U,NDX,NTST,NCOL)
      INTEGER, INTENT(IN) :: NDX,NCOL,NTST
      DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
      DOUBLE PRECISION  :: PCLU_Y(0:NTST*NCOL)

      PCLU_Y=U(1,:)

      GETUY_MIN=MINVAL(PCLU_Y)

      END FUNCTION GETUY_MIN
!-----------------------------------------------------------

      SUBROUTINE PVLS(NDIM,U,PAR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
      DOUBLE PRECISION, EXTERNAL :: GETP,GETUY_MIN2,GETUY_MIN
      INTEGER NDX,NCOL,NTST

      NDX=NINT(GETP('NDX',0,U))
      NTST=NINT(GETP('NTST',0,U))
      NCOL=NINT(GETP('NCOL',0,U))

      PAR(4)=GETUY_MIN(U,NDX,NTST,NCOL)
      PAR(5)=GETUY_MIN2(U,NDX,NTST,NCOL)

      END SUBROUTINE PVLS


