      SUBROUTINE DFT(S,FS,NSAMP,FMIN1,FMAX1,FMIN2,FMAX2,ALIM,AMP,
     ,F1,A1,P1,F2,A2,P2)
      IMPLICIT NONE
      REAL*8 S,FS,FMIN1,FMAX1,FMIN2,FMAX2,ALIM,AMP(100)
      REAL*8 F1,A1,P1,F2,A2,P2
      COMMON/DFTDATA/ SIG(500000), NS,FSUB,ASUB,PSUB,FSUP,ASUP,PSUP
      REAL*8  SIG, FSUB,ASUB,PSUB,FSUP,ASUP,PSUP
      DATA  FSUB/0.0/,ASUB/0.0/,PSUB/0.0/,FSUP/0.0/,ASUP/0.0/,PSUP/0.0/
      INTEGER NSAMP, NS, k, n
      REAL*8 DT, SRE, SIM, AM, FR, PH
      REAL*8 PI
      DATA NS/0/, PI/3.141592653589793D0/
      
      DT=1.0/FS
      NS=NS+1
      SIG(NS)=S
      IF (NS == NSAMP) THEN
         NS=0
         DO k=1,100
            AMP(k)=0.0
         ENDDO
         ASUB=0.0
         ASUP=0.0
         DO k=1,100
            SRE=0.0
            SIM=0.0
            DO n=0,NSAMP-1
               SRE=SRE+SIG(n+1)*COS(2*PI*n*(k-1)/NSAMP)
               SIM=SIM-SIG(n+1)*SIN(2*PI*n*(k-1)/NSAMP)
            ENDDO
            AM = 2*SQRT(SRE*SRE+SIM*SIM)/NSAMP
            FR = 1.0*(k-1)/DT/NSAMP
            PH = DATAN2(SIM,SRE)*180.0/PI
            IF (AM >=ALIM) THEN
              AMP(NINT(FR))=AM
              IF (FR>=FMIN1 .AND. FR<=FMAX1 .AND. AM>ASUB ) THEN
                 ASUB=AM
                 FSUB=FR
                 PSUB=PH
                 WRITE(1,*) 'DFT k, F, A, PH: ',k, FR, AM, PH
              ENDIF
              IF (FR>=FMIN2 .AND. FR<=FMAX2 .AND. AM>ASUP ) THEN
                 ASUP=AM
                 FSUP=FR
                 PSUP=PH
              ENDIF                 
            ENDIF
         ENDDO
      ENDIF
      A1=ASUB
      F1=FSUB
      P1=PSUB
      A2=ASUP
      F2=FSUP
      P2=PSUP
      RETURN
      END


      SUBROUTINE CORR3P(X0, F0, DELT, AM0, PH0, AMOUT, FROUT, PHOUT)
      IMPLICIT NONE
      REAL*8 X0, F0, DELT, AM0, PH0, AMOUT, FROUT, PHOUT
      COMMON/C3PDATA/ NCOR,SIG(500000),NSAMP,M, NSIG, F1, AMO, FRO, PHO
      REAL*8 F1, SIG, AMO, FRO, PHO
      INTEGER NCOR, NSAMP, M, NSIG
      
      COMMON/C3PCNST/ PI, i
      REAL*8 PI
      COMPLEX*16 i
      COMPLEX*16 XK1, XK2, XK3, CFFT
      COMPLEX*16 G_R, E_LPHA, C_R, X_R, XOUT
      REAL*8 PHA, DF
      DATA NSIG/0/, NCOR/0/, PI/3.141592653589793D0/, i/(0D0,1D0)/
      INTEGER N1,N2,N3
      
      IF (F0<1.0) RETURN
      
      IF (NCOR.EQ.0) THEN
          NSAMP=NINT(1D0/DELT/F0)
          M=CEILING(NSAMP/4D0)
          F1=1D0/DELT/NSAMP
          NCOR=1
      ENDIF
      
      NSIG=NSIG+1
      SIG(NSIG)=X0
      
      IF (NSIG.EQ.2*NSAMP) THEN
!       IF (NSIG.EQ.NSAMP+2*M) THEN
!       WRITE(3,*) NSIG
        N1=NSAMP+1
        N2=N1-M
        N3=N2-M
        XK1=CFFT(SIG(N1),NSAMP,1)
        XK2=CFFT(SIG(N2),NSAMP,1)*CDEXP(2*PI*i*M/NSAMP)
        XK3=CFFT(SIG(N3),NSAMP,1)*CDEXP(4*PI*i*M/NSAMP)
!       WRITE(3,*) NCOR,XK1, XK2, XK3
        
        G_R=(XK3+XK1*CDEXP(4*PI*M*I/NSAMP))/XK2/2.0D0
        E_LPHA=G_R+CDSQRT(G_R*G_R-CDEXP(4D0*PI*M*i/NSAMP))
        C_R=-(XK2-XK1*E_LPHA)/(E_LPHA-DCONJG(E_LPHA)*
     *     CDEXP(4D0*PI*i*M/NSAMP))
        PHA=-DATAN2(DIMAG(E_LPHA),DREAL(E_LPHA))/M
        X_R=(XK1-C_R)*(NSAMP*DSIN(PHA/2D0)/DSIN(PHA*NSAMP/2D0)) *
     *     CDEXP(-i*PHA*(NSAMP-1)/2D0)
        DF=(PHA*NSAMP*F1)/2D0/PI
!        WRITE(3,*) 'PHA ', PHA

!  Original code from Matlab:
!         g_r=(xk3_ysxl/ xk2_ysxl+exp(1i*4*pi*f_xiu/N_s)*xk1_ysxl/xk2_ysxl)/2;
!         
!         e_lpha=g_r+sqrt(g_r^2-exp(1i*4*pi*f_xiu/N_s));
!         c_r=-1*(xk2_ysxl-xk1_ysxl*e_lpha)/(e_lpha-conj(e_lpha)*exp(1i*4*pi*f_xiu/N_s));
!         
!         pha=angle(e_lpha)/(-1*f_xiu);
!         
!         x_r=(xk1_ysxl-c_r)*(N_s*sin(pha/2)/sin(pha*N_s/2))*exp(-1i*pha*(N_s-1)/2);
!         df=(pha*N_s*f1)/(2*pi);

!        WRITE(3,*) X_R, F1, DF

        AMO=CDABS(X_R)
        PHO=DATAN2(DIMAG(X_R),DREAL(X_R))*180D0/PI
        FRO=F1+DF
        WRITE(7,*) FRO, AMO, PHO
        NSIG=0
        NCOR=NCOR+1
      ENDIF            
      AMOUT=AMO
      FROUT=FRO
      PHOUT=PHO
      RETURN
      END
        
        COMPLEX*16 FUNCTION CFFT(X,NS,K)
        IMPLICIT NONE
        REAL*8 X(NS)
        INTEGER NS, K, N
        COMMON/C3PCNST/ PI, i
        REAL*8 PI
        COMPLEX*16 i
        COMPLEX*16 SUM
        
        SUM=(0D0,0D0)
        DO N=0,NS-1
           SUM=SUM+X(N+1)*CDEXP(-i*2D0*PI*N*K/NS)
        ENDDO
        CFFT=2D0*SUM/NS
        RETURN
        END
        
        
        SUBROUTINE FFTK(DELT,SIG,NS,K,AMP,PHASE,FREQ)
        IMPLICIT NONE
        REAL*8 DELT, SIG(NS), AMP, PHASE, FREQ
        INTEGER NS, K
        COMMON/C3PCNST/ PI, i
        REAL*8 PI
        COMPLEX*16 i
        COMPLEX*16 SUM
        INTEGER N
        SUM=(0D0,0D0)
        DO N=0,NS-1
           SUM=SUM+SIG(N+1)*CDEXP(-i*2D0*PI*N*(K-1)/NS)
        ENDDO
        AMP = 2D0*CDABS(SUM)/NS
        PHASE = DATAN2(DIMAG(SUM),DREAL(SUM))*180D0/PI
        FREQ = 1D0*(K-1)/DELT/NS
        !WRITE(4,*) FREQ,AMP,PHASE
        RETURN
        END
        
        
      SUBROUTINE DFT1(S,FS,NSAMP,FMIN1,FMAX1,FMIN2,FMAX2,ALIM,AMP,
     ,F1,A1,P1,F2,A2,P2)
      IMPLICIT NONE
      REAL*8 S,FS,FMIN1,FMAX1,FMIN2,FMAX2,ALIM,AMP(100)
      REAL*8 F1,A1,P1,F2,A2,P2
      COMMON/DFT1DATA/ SIG(500000), NS,FSUB,ASUB,PSUB,FSUP,ASUP,PSUP
      REAL*8  SIG, FSUB,ASUB,PSUB,FSUP,ASUP,PSUP
      DATA  FSUB/0.0/,ASUB/0.0/,PSUB/0.0/,FSUP/0.0/,ASUP/0.0/,PSUP/0.0/
      INTEGER NSAMP, NS, k, n
      REAL*8 DT, SRE, SIM, AM, FR, PH
      REAL*8 PI
      DATA NS/0/, PI/3.141592653589793D0/
      
      DT=1.0/FS
      NS=NS+1
      SIG(NS)=S
      IF (NS == NSAMP) THEN
         NS=0
         DO k=1,100
            AMP(k)=0.0
         ENDDO
         ASUB=0.0
         ASUP=0.0
         DO k=1,100
            SRE=0.0
            SIM=0.0
            DO n=0,NSAMP-1
               SRE=SRE+SIG(n+1)*COS(2*PI*n*(k-1)/NSAMP)
               SIM=SIM-SIG(n+1)*SIN(2*PI*n*(k-1)/NSAMP)
            ENDDO
            AM = 2*SQRT(SRE*SRE+SIM*SIM)/NSAMP
            FR = 1.0*(k-1)/DT/NSAMP
            PH = ATAN2(SIM,SRE)*180.0/PI
            IF (AM >=ALIM) THEN
              AMP(NINT(FR))=AM
              IF (FR>=FMIN1 .AND. FR<=FMAX1 .AND. AM>ASUB ) THEN
                 ASUB=AM
                 FSUB=FR
                 PSUB=PH
                 WRITE(2,*) 'DFT1 k, F, A, PH: ',k, FR, AM, PH
              ENDIF
              IF (FR>=FMIN2 .AND. FR<=FMAX2 .AND. AM>ASUP ) THEN
                 ASUP=AM
                 FSUP=FR
                 PSUP=PH
              ENDIF                 
            ENDIF
         ENDDO
      ENDIF
      A1=ASUB
      F1=FSUB
      P1=PSUB
      A2=ASUP
      F2=FSUP
      P2=PSUP
      RETURN
      END
        