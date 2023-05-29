  Subroutine Update
  Use precisions,only:int_p
  Use grid !,only: ncells
  Use Variables !,only:gcone,gc,tnotone,tnot,ndir,intensityone,intensity,nbands, &
                !     io_debug, debug_level, np_max, polar_type, highband, lowband, temp_bcone
  USE grid                  
  Implicit none
  integer(int_p)::ic,si,iband, p, pol, icount,i
  
  IF ((debug_level > 0).and.(rank == 0)) WRITE(io_debug,'(a)') "Starting update" 

  do ic = 1, ncells 
     do iband = 1,nbands
        DO p = 1, np_max
        !    pol = polar_type(iband, p)
         !   IF (pol .eq. 0) CYCLE     
            gcone(ic,p,iband) = gc(ic,p,iband)  
        ENDDO       
     enddo
     
     tnotone(ic) = tnot(ic)
       
  end do
  
  do iband = lowband, highband !1, (highband - lowband + 1)    !nbands
  !do iband = 1, nbands
    icount = iband - lowband + 1
    DO p = 1, np_max
     !  pol = polar_type(iband, p)
     !  IF (pol .eq. 0) CYCLE     
        do ic = 1,ncells
            do si = 1,ndir
                intensityone(ic,si, icount,p) = intensity(ic,si,icount,p)            
            enddo
        enddo
    ENDDO
  enddo
  
!O i = 1, nbcfaces
  temp_bcone(:) =temp_bc(:)
  ttopone(:) = ttop(:)
!NDDO
  
IF ((debug_level > 0).and.(rank == 0)) WRITE(io_debug,'(a)') "Finishing update"
  
End subroutine update
