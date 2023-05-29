!***********************************************************************
! Declaration of Precision of Variables

  MODULE precisions

  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: int_p  = SELECTED_INT_KIND(8)
  INTEGER, PARAMETER :: real_p = SELECTED_REAL_KIND(8)
  ! INTEGER, PARAMETER :: real_s = SELECTED_REAL_KIND(4)

  END MODULE precisions
!***********************************************************************
