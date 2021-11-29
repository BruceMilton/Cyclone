c      subroutine curvature(nturn,ip,theta,radius,pr,e,rc,psi)
      subroutine curvature(nturn,ip,theta,radius,pr,e)
      use cyclone_data
      use bfields
      use IO

      real*8 theta,radius,pr,e
      real*8 PSQ,GAM
      real*8 SINA,RHO,ALFA,RC
C********************************************************
C
C      CALCULATES INSTANTANEOUS CENTER OF CURVATURE
C
C*******************************************************
C
c      RADIUS=Y(1)
c      PR=Y(2)
      PSQ = E / E0
      GAM = 1.0D+0 + PSQ
      PSQ = PSQ * ( 2.0D+0 + PSQ ) * ASQ
      IF (.not. BFLDCAL(theta,radius,.false.,0) ) go to 99
!      call bfldcal(theta,radius,.false.,1,*99)
      PSQ=DSQRT(PSQ)
      RHO=PSQ/BZ
      SINA=PR/PSQ
      ALFA=ASIN(SINA)
      RC=SQRT(RADIUS*RADIUS+RHO*RHO-2.*RADIUS*RHO*COS(ALFA))
      IF(RC.EQ.0.)GOTO 90
      PSI=SINA*RHO/RC
      IF(ABS(PSI).gt. 1.0) then
        PRINT 5000,PSQ,RHO,SINA,ALFA,RC,PSI
5000    FORMAT(' P='F10.5'   RHO='F10.5'   SINA='F10.5'   ALFA='
     +          F10.5'  RC='F10.5'   SINPSI='F10.5)
      endif
      IF(PSI.GT.1.)PSI=1.
      IF(PSI.LT.-1.) PSI=-1.  
      PSI=ASIN(PSI)
      IF(RHO*COS(ALFA).GT.RADIUS .AND. PSI.GE.0. ) PSI=3.1415926-PSI
      IF(RHO*COS(ALFA).GT.RADIUS .AND. PSI.LT.0.) PSI=-3.1415926-PSI
      PSI=PSI+THeta
      PSI=PSI * TCON
      th=theta*tcon
90    continue
      if(io_control36.formatted) then
            write(36,4514)nturn,ip,th,radius,pr,e,rc,psi
4514        format(2i5,5f12.7,f12.2,2f12.7)
      else
            write(36)nturn,ip,th,radius,pr,e,rc,psi
      endif
!      ifile36_lines=ifile36_lines+1
      io_control36.lines=io_control36.lines+1
99    return
      END

