!!!!/////// Uniform Random number generators////////////////////////////////////

FUNCTION ran2(idum)
 ! USE numz
  IMPLICIT NONE
  DOUBLE PRECISION:: ran2
  !INTEGER,INTENT(inout),OPTIONAL::idum
  INTEGER,INTENT(inout)::idum
  !INTEGER :: idum
  INTEGER,PARAMETER::IM1=2147483563,IM2=2147483399,IMM1=IM1-1
  INTEGER,PARAMETER::IA1=40014,IA2=40692,IQ1=53668
  INTEGER,PARAMETER::IQ2=52774,IR1=12211,IR2=3791   
  INTEGER,PARAMETER::NTAB=32,NDIV=1+IMM1/NTAB
  DOUBLE PRECISION,PARAMETER::AM=1.0d0/IM1,EPS=1.2e-7,RNMX=1.0d0-EPS
  INTEGER::idum2,j,k,iv(NTAB),iy
  SAVE iv,iy,idum2
  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
  IF (idum<0) THEN
     idum=MAX(-idum,1)
     idum2=idum
      DO j=NTAB+8,1,-1
         k=idum/IQ1
         idum=IA1*(idum-k*IQ1)-k*IR1
         IF (idum<0) idum=idum+IM1
         IF (j.LE.NTAB) iv(j)=idum
      ENDDO
      iy=iv(1)
   ENDIF
   k=idum/IQ1
   idum=IA1*(idum-k*IQ1)-k*IR1
   IF (idum<0) idum=idum+IM1
   k=idum2/IQ2
   idum2=IA2*(idum2-k*IQ2)-k*IR2
   IF (idum2<0) idum2=idum2+IM2
   j=1+iy/NDIV
   iy=iv(j)-idum2
   iv(j)=idum
   IF(iy.LT.1)iy=iy+IMM1
   ran2=MIN(AM*iy,RNMX)
   RETURN
 END FUNCTION ran2
