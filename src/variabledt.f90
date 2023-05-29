!!$ This subroutine calculates common ratio or separation constant s for different timestep
!! calculating s from the geometric progression formula

Subroutine variabledt
Use precisions !,Only:real_p,int_p
Use constants !,only:wmaxlas,wmintas,wmaxtas,wminlas,vslas,vstas,clas,ctas,six,four,half
Use variables
Implicit None

REAL (real_p) :: t 


REAL (real_p) :: tol_nw, po, d, fo, fpo  !TOL,P0,D,F0,FP0
INTEGER (int_p) :: i, no!,I,N0,FLAG,OUP
!CHARACTER NAME1*30,AA*1
! LOGICAL ok

t = time_pulse/(dt_1) 
!The 2 is only for a sine pulse
!  DEFINE FUNCTIONS F AND F' (DENOTED FP)
!     f(sep)=sep**n - t*sep  + (t-1)
!     fp(sep)=n*(sep**(n-1)) - t
!!! tolerances, max iteration and initial approximation   
tol_nw = tol_outer
no = max_iter
po = two
DO i = 1, no
      fo = po**n - t*po  + (t-1)  !fo = f(po)    !F0=F(P0)
      !     STEP 3
      !     COMPUTE P(I)
      fpo=n*(po**(n-1)) - t    !fpo=fp(po)
      d=fo/fpo
      !    STEP 6
      po=po-d
      !      fo=f(po)
      !     STEP 4
      IF( ABS(d) .LT. tol_nw ) THEN
            !    PROCEDURE COMPLETED SUCCESSFULLY
            sep = po
            EXIT   
      END IF
ENDDO
!!!solve for sep using newton raphson method. code reference burden:numerical analysis 2001 edition
!...Han added only root prints this out
IF (rank .eq. 0) WRITE (*,'(a,1x,f20.15)') 'variabledt: sep=', sep
!IF (sep == 1.00)THEN
!  WRITE(*,*) "dont need this root!!"
!     STOP
!  ENDIF

!sep = 1 !!! will be an equation
!STOP
End Subroutine variabledt
  
 