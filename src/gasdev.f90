!!!!!!!!/////// Gaussian random no. generator \\\\\\\\\\\\\\\\\\\\\\\\\\\
! Polar rejection method (Knop[1969]) related to Box-Muller transform
 Subroutine gasdev(g2,idum,mean,variance) 
  !use numz
  Implicit none
      INTEGER,INTENT(INOUT)::idum
      DOUBLE PRECISION,INTENT(IN)::variance,mean
      DOUBLE PRECISION,INTENT(OUT)::g2
      DOUBLE PRECISION::ran2,g1
      INTEGER:: iset
      DOUBLE PRECISION:: fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
    !  if (iset.eq.0) then
        DO
        v1=2.0d0*ran2(idum)-1.0d0
        v2=2.0d0*ran2(idum)-1.0d0
        rsq=v1**2+v2**2
        if((rsq<1.0d0).AND.(rsq/=0.0d0))EXIT
        ENDDO
        fac=variance*SQRT(-2.0d0*log(rsq)/(rsq))
        g1=v1*fac+mean
        g2=v2*fac+mean      
      END subroutine gasdev
