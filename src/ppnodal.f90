  Subroutine ppnode
  
  Use PRECISIONS
  Use CONSTANTS,Only : zero,two
  Use VARIABLES,Only : tnot,tnotv, debug_level, io_debug
  Use GRID !,Only : ncells,nbcfaces,lcv,lfv,nfcell,bf_to_f,nnodes, &
           !       bnode,nvcell,nvface,xc,xv,xf,yc,yv,yf,zc,zf,zv, &
           !       temp_bc,bctype,lfc,wcv,wbfv,ADIA,ISOT,CONJ, SYMM, ALSI 
  Implicit None
  
  Integer(int_p) :: i,j,currv,currf,indx,icell,inode, ifc, ic1, ic2,ib
  INTEGER(int_p) :: currbf, jv_i, jv_ifc, inode_a, inode_s
  
  Real(real_p) :: disntemp,tcell, temp
  REAL(real_p), DIMENSION(:), ALLOCATABLE :: weight
  
  IF(debug_level > 0) WRITE(io_debug,'(a)') "Starting ppnodal_ppnode" 


 !OPEN(unit=66, file="wsei.dat",  status="unknown")
  
  ALLOCATE(weight(nnodes))
  tnotv(:) = zero
  weight(:) = zero
  
  !! interior nodes, it was there
  
  DO icell=1,ncells
        tcell = tnot(icell)		
        DO j=1,nvcell(icell)
            inode = lcv(icell,j)
            weight(inode) = weight(inode) + wcv(icell,j)
            tnotv(inode) = tnotv(inode) + tcell*wcv(icell,j)
        END DO
  END DO
  
!    DO inode = 1,nnodes
!       IF(bnode(inode) /= 1)CYCLE
!       tnotv(inode) = tnotv(inode)/weight(inode)
!     ENDDO
     
!! boundary nodes need to be overwritten. !!jul 29- assy

     DO inode = 1,nnodes
       IF(bnode(inode) == 1) THEN
            tnotv(inode) = 0
       END IF
     ENDDO
     weight(:) = 0

	DO i = 1,nbcfaces
	   IF (bctype(i) == CONJ) CYCLE
       ifc = bf_to_f(i)
       ic1 = lfc(ifc,1)
       ic2 = lfc(ifc,2)
!        condf = wcf(ifc)*cond(ic1) + (one - wcf(ifc))*cond(ic2)
!       	IF(ic1 /= ic2) condf = cond(ic1)*cond(ic2)/condf ! Conjugate
       temp = temp_bc(i)
       DO j=1, nvface(ifc)	!1,nfvertex(ifc)
           inode = lfv(ifc,j)
           weight(inode) = weight(inode) + wbfv(i,j)
           !IF (bctype(i) == CONJ) !CYCLE 		!*condf
            tnotv(inode) = tnotv(inode) + temp*wbfv(i,j)		!*condf
           !ENDIF
       END DO
!       IF (bctype(i) == ISOT) THEN
!           ib = bf_to_f(i)
!           Do j = 1,nvface(ib)
!              inode = lfv(ib,j)
!              tnotv(inode) = temp_bc(i)
!          End Do
!       ENDIF
     END DO
  !   DO inode = 1,nnodes
  !     WRITE(66,*) inode, tnotv(inode), weight(inode), bnode(inode)
 !    ENDDO
     
    DO inode = 1,nnodes
       IF(bnode(inode) /= 1)CYCLE
       tnotv(inode) = tnotv(inode)/weight(inode)
     ENDDO
     
  

  DO i = 1,nbcfaces
    IF (bctype(i) .ne. SYMM) CYCLE 	!! cycle any face except symmetry
       ifc = bf_to_f(i)
       ic1 = lfc(ifc,1)
       ic2 = lfc(ifc,2)
      ! print *, i, ifc, ic1,ic2 
	DO j = 1,nfcell(ic1)
          currf = lcf(ic1,j)
          IF (currf == ifc) CYCLE			!! cycle same face
          IF (bface(currf) == 0) CYCLE		!! cycle interior face
	   currbf = f_to_bf (currf) 
	   IF (bctype(currbf) == ALSI) THEN
                !WRITE (*,*) ic1, ic2, currbf,i
            DO jv_i = 1,nvface(ifc)
	      inode_s = lfv(ifc,jv_i)
	      DO jv_ifc = 1,nvface(currf)
	      	  inode_a = lfv(currf,jv_ifc)
	         IF (inode_a == inode_s) THEN
                    tnotv(inode_a) = temp_bc(currbf)
	            !WRITE (*,*) inode_s, inode_a, tnotv(inode_a)
                ENDIF
              ENDDO
            ENDDO   
	  ENDIF	
       ENDDO 
!   ENDIF    
  END DO
!WRITE (*,*) "node 10251", tnotv(10251), xv(10251), yv(10251)
Do i = 1,nbcfaces
    IF (bctype(i) == ISOT) THEN
    ! IF ((bctype(i) == ISOT) .or. (bctype(i)== SYMM)) THEN

        ifc = bf_to_f(i)
           Do j = 1,nvface(ifc)
              inode = lfv(ifc,j)
              tnotv(inode) = temp_bc(i)
          End Do
       ENDIF
  ENDDO
!        currf = bf_to_f(i)
!        ic1 = lfc(currf,1)
!       	ic2 = lfc(currf,2)
!        indx = bctype(i)
!        Select Case(indx)
!            Case(ISOT)
!        Do j = 1,nvface(currf)
!               currv = lfv(currf,j)
!!!$        disntemp = Sqrt( (xf(currf)-xv(currv))**two + (yf(currf)-yv(currv))**two + (zf(currf)-zv(currv))**two )
!               tnotv(currv) = temp_bc(i)
!        End Do
!     
!!        Case(ADIA)
!!         temp = temp_bc(i)
!!            Do j = 1,nvface(currf)
!!                currv = lfv(currf,j)
!!                weight(currv) = weight(currv) + wbfv(i,j) 		!*condf
!!                tnotv(currv) = tnotv(currv) + temp*wbfv(i,j)		!*condf
!!!                tnotv(currv) = tnot(icell)
!!            End Do
!!            DO inode= 1,nnodes
!!                IF(bnode(inode) /= 1)CYCLE
!!                    tnotv(inode) = tnotv(inode)/weight(inode)
!!            ENDDO
!            
!        End Select
!  End Do
    

  DEALLOCATE(weight)
  
  IF(debug_level > 0) WRITE(io_debug,'(a)') "Finishing ppnodal_ppnode" 
  
  End Subroutine ppnode
