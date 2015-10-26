!********************************************************************** 
      subroutine get_crack_tip(x,b) 
!     Calculate the crack tip position from atomistics
!     Atomistic crack tip position is computed using the Neighbor list
!     of the atoms
!     This is used rather than a coordination number to avoid any
!     other calculations
!     In this instance Surface atoms are identified using the number of
!     the near neighbors < 23
!     Perfect atoms = 26
!     Dislocations > 26
!     The crack tip position is max(x) of the surface atoms 

      use mod_global 
      use mod_dynamo
      implicit double precision (a-h, o-z) 
      double precision x(nxdm, *), b(ndf, *)
      double precision box_max(2), box_min(2) 
      double precision xtip_min(2) 
      xtip_min = -1.e30 
      do j=1,2 
         box_max(j)=-1e30 
         box_min(j)=1E30 
      enddo 
      npatoms=0 
      npad = 0 
      ntot = 0 
!     ********************************************************
!     ---- Calculate the Box size ignoring the interface atoms 
      do i = 1,numnp 
         if( IsRelaxed(i)==1 .or. IsRelaxed(i)==2 .or.                  &
     &        IsRelaxed(i)==-1) then                                    
            if (IsRelaxed(i) ==-1) then 
               if (npad == 0) then 
                  startpad = i 
               endif 
               npad = npad + 1 
            else 
               pe = pe + energy(i) 
               npatoms=npatoms+1 
            endif 
            ntot = ntot + 1 
            do j=1,2 
               if (IsRelaxed(i) .ne. -1) then 
                  if(  x(j,i) > box_max(j)) then 
                     box_max(j)= x(j,i) 
                  endif 
                  if(  x(j,i) < box_min(j)) then 
                     box_min(j)= x(j,i) 
                  endif 
               end if 
            enddo 
         endif 
      enddo 
      startpad = npatoms - npad + 1 
!     ********************************************************
!     Compute surface atoms and therefore crack tip atom pos
!     Surface atoms are classified as those with less than the 
!     26 nearest neighbors 
!     ********************************************************
      imin = -1000000
      do i = 1, numnp 
         if (IsRelaxed(i) == 1 .or. IsRelaxed(i) ==2) then 
            if (x(1,i) < box_max(1) .and. x(1,i) > box_min(1)) then 
               if (x(2,i) < box_max(2) .and. x(2,i) > box_min(2)) then 
                  if (NumNeighbors(i) .le. 23) then 
                     if (x(1,i) > xtip_min(1)) then 
                        xtip_min(1) = x(1,i) 
                        xtip_min(2) = x(2,i) 
                        imin = i 
                     end if 
                  end if 
               end if 
            end if 
         end if 
      end do 
      if (imin > 0) then 
         xtip(1) = xtip_min(1) 
         xtip(2) = xtip_min(2) 
         x1 = (x(1,imin)+box_min(1))/(box_max(1)-box_min(1)) 
         y1 = (x(2,imin)+box_min(2))/(box_max(2)-box_min(2)) 
         write(*, '(A20,2I8,5F16.11)')'Crack Tip Position = ', imin,       
     $        NumNeighbors(imin), xtip_min(1),x(1,imin),x(2,imin),x1, y1  
         write(*,'(A20,2F16.11)') 'New Crack Tip = ', x(1,imin)
     $        + b(1,imin), x(2,imin)+b(2,imin)
                                   
      else
         print *, 'Crack not in atomistic region **********'
         print *, 'Assuming positio for crack tip'
         xtip_min(1) = -1.e-30
         xtip_min(2) = 0.0d0
      end if
                                                                        
c$$$      xtip = xtip_min
      return
      end subroutine get_crack_tip 
!********************************************************************** 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc        
!     Move crack tip back to (x0,y0)                                      
!     This means moving the AtomDispl from current atom position        
!     to the new one                                                    
!     1) xtip(1) contains current x position of crack                   
!     2)                                                                
!                                                                       
      subroutine move_atomistic_crack(x, ix, b) 
      use mod_global 
      use mod_file 
      use mod_dynamo
      implicit none 
      include '../Disl/disl_parameters.par' 
                                                                        
      double precision b(ndf,*), x(nxdm,*) 
      integer ix(nen1,*), iq
                                                                        
      integer iAtom, j, i,iFem 
                                                                        
      integer logic 
      integer, allocatable ::  mmap (:)
      double precision :: xnew(nxdm), xold(nxdm), rr, box_min(2),       &
     &     box_max(2), umag, bb(nxdm)                                   
      double precision :: xtip_move
      integer ntot, npatoms 
      double precision xdef, ydef, zdef 
      character*80 input,filename,temp 
      character*4 key 
      character*6 cnt
      integer numatoms
      double precision, allocatable :: new_b(:,:)
      integer, allocatable :: new_numneigh(:), new_neighlist(:,:)
      double precision, allocatable :: utilde(:,:), new_utilde(:,:)

      allocate(new_b(ndf, numnp))
      allocate(new_numneigh(numnp))
      allocate(new_neighlist(MAXNEIGHBORS,numnp))
      umag = 1.0 
      iq = 0
      filename ="out/before_crack_atoms.cfg" 
      call iofile(filename,'formatted  ',logic,.false.) 
      call dump_atom(x, b, logic)
      close(logic)


      filename ="out/before_mm.vtk" 
      call iofile(filename,'formatted  ',logic,.false.) 
      call dump_mesh(x,b,ix, logic, iq)
      close(logic)

      allocate(mmap(numnp)) 
      mmap = 0 
      print *, 'Crack tip motion = ', xtip_actual*x_tip_dir
      xtip_move = int(xtip_actual(1)/x_move_mesh)*x_move_mesh
      print *, 'Actual crack tip motion = ', int(xtip_actual
     $     /x_move_mesh)*x_move_mesh
      
      do i = 1,2
         crack_motion(i) = crack_motion(i) + x_tip_dir(i)*xtip_move
      end do
!     Calculate the utilde vector everywhere
      allocate(utilde(3,numnp), new_utilde(3,numnp))
      utilde = 0.0d0
      new_utilde = 0.0d0
      do i = 1, numnp
         call disl_displ(x(1:3,i), utilde(1:3,i))
      end do
      
!     Create mapping from new coordinates to old coordinates            
      do iAtom = 1, numnp 
         if (isRelaxed(iAtom) .eq. 1 .or. isRelaxed(iAtom) .eq. 2 
     $       .or. isRelaxed(iAtom) .eq. -1) then 
            do j = 1,1
!!               xold(j) = x(j,iAtom) + x_tip_dir(j)*xtip_actual(j)
               xold(j) = x(j,iAtom) + xtip_move
            end do
            xold(2:nxdm) = x(2:nxdm, iAtom) 
            do i = 1, numnp 
               if (isRelaxed(i) .eq. 1 .or. isRelaxed(i) .eq. 2 .or.    &
     &              isRelaxed(i) .eq. -1) then                          
                  xnew(1:nxdm) = x(1:nxdm, i) 
!     Calculate distance between 2 atoms                                
                  rr = 0.d0 
                  do j = 1, nxdm 
                     rr = rr + (xnew(j)-xold(j))**2 
                  end do 
                  rr = dsqrt(rr) 
                  if (rr < 1.d-3) then 
                     mmap(iAtom) = i 
                     new_b(:,i) = b(:,i)
                     new_numneigh(i) = NumNeighbors(i)
                     new_neighlist(:,i) = NeighborList(:,i)
                     exit 
                  end if 
               end if 
            end do 
         end if 
   31    continue 
!$$$  write(*,*) 'MMap', iAtom, mmap(iAtom)                             
      end do 
      do iAtom = 1, numnp 
         if (isRelaxed(iAtom) .eq. 1 .or. isRelaxed(iAtom) .eq. 2 .or.  &
     &        isRelaxed(iAtom) .eq. -1) then                
c     Copy displacements from atom i to iAtom to move the crack tip 
c     Also copy velocities and neighbour lists ???
            i = mmap(iAtom) 
c$$$            write(*,'(A10,2I7,3(1X,E15.8),2I7)'),
c$$$     $           'Mapping = ',iAtom, i
c$$$     $           , x(1,iAtom),x(1:nxdm,i), b(1,iAtom), b(1,i)
            if (i .ne. 0) then 
c$$$               write(*,'(A13,2E15.8,A4,2E15.8, A8, 4(1X,E15.8)),2I7')
c$$$     $              'Mapping from ', x(1:2,i), ' to ', x(1:2,iAtom),
c$$$     $              ' Disp = ', b(1:2,iAtom), b(1:2, i),
c$$$     $              NumNeighbors(iAtom), NumNeighbors(i)
               b(1:ndf,iAtom) = new_b(1:ndf, i) 
c     *****************************************************************
c     Updates the neighbor list as well 
c$$$               NumNeighbors(iAtom) = new_numneigh(i)
c$$$               NeighborList(:,iAtom) = new_neighlist(:,i)
c     *****************************************************************
            end if
c     End copy displacements                                                                         
         end if 
      end do 


      call move_slip_planes(xtip_move)

      do i = 1, numnp
         call disl_displ(x(1:3,i), new_utilde(1:3,i))
      end do
      do i = 1, numnp
         if (isRelaxed(i) == 0) then 
            b(1:3,i) = b(1:3,i) - utilde(1:3,i)
     $           + new_utilde(1:3,i)
         end if
      end do

      deallocate(new_b)
      deallocate(mmap)
      deallocate(new_numneigh)
      deallocate(new_neighlist)


c$$$c     Just use the finite deformation quantities of the atomistics
c$$$c     b(:,i) = b(:,i) - xtip_actual
c$$$c     Since this will be rigid body translation
c$$$c     The deformation gradient within the atomistics will not change
c$$$c     Right now pad atoms are also moved 
c$$$c     They should be moved 
c$$$      do iAtom = 1, numnp
c$$$         if (isRelaxed(iAtom) .eq. 1 .or. isRelaxed(iAtom) .eq. 2
c$$$     $        .or. isRelaxed(iAtom) .eq. -1) then 
c$$$           do j = 1, 1
c$$$              b(j,iAtom) = b(j,iAtom) - x_tip_dir(j)*xtip_actual(j)
c$$$           end do
c$$$         end if
c$$$      end do


      filename ="out/after_crack_atoms.cfg" 
      call iofile(filename,'formatted  ',logic,.false.) 
      call dump_atom(x, b, logic)
      close(logic)

      filename ="out/after_mm.vtk" 
      call iofile(filename,'formatted  ',logic,.false.) 
      call dump_mesh(x,b,ix, logic)
      print *, 'Total Crack tip motion (deltaA) = ', crack_motion
!     Set newmesh parameter to ensure that the DB strains are recalculated
      newmesh = .true. 

      close(logic)

      end subroutine move_atomistic_crack 
!**********************************************************************

      subroutine move_slip_planes(xtip_move)
      use mod_global
      use mod_dd_slip
      integer islp
      double precision :: xtip_move
      do islp = 1, nslp
         xslp(islp) = xslp(islp) - xtip_move
      end do
      call gen_slip_ends

      call assign_new_global(xtip_move)
      end subroutine move_slip_planes
