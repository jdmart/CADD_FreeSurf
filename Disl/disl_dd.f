c     
c     $Id: disl_dd.f,v 1.1.1.1 2003/03/12 20:09:00 shastry Exp $
c     
      subroutine disl_setup()
      use mod_dd_slip
      implicit none
      include 'disl_parameters.par'
      character*80 error_message
c     
      data i_disl /0/
c     
      if(i_disl.eq.0) then
         i_disl = 1
         call gen_sources
         call gen_obstacles
         call plot_slip
c$$$  stop
c     call gen_slip_planes
      else
         error_message = 'Second call to disl_setup'
         call error_handler(error_message)
      endif
      ndisl = 0
      return
      end


      subroutine disl_accept(r0, bv, th_e, th_s)
      implicit none
      include 'disl_parameters.par'
      double precision r0(3), bv(3), th_e, th_s, rr, RTOL
      parameter(RTOL=1.d0)
      integer ishift
      logical moved
c     
      integer n_total
      common n_total
      data n_total /0/
c     
      character*80 error_message
      integer fe_locate, i, j
c     
      if(i_disl.ne.1) then
         error_message = 'disl_accept: call disl_setup first!'
         call error_handler(error_message)
      endif
c     
      n_total = n_total + 1
c     
c     forced shut-off of lumping:
c     
      n_total=1
c     
      if((n_total.eq.1).or.(n_total.eq.2)) then
         ndisl = ndisl+1
         print *, '***************  Holding ', ndisl,
     &        ' dislocations *****************'

         if(ndisl.gt.max_disl) then
            error_message = 'max_disl needs to be increased'
            call error_handler(error_message)
         endif
c     if(ndisl.ge.2) then
c     stop 'Amit nucleation stop'
c     endif
         elem_disl(ndisl)=fe_locate(r0,1)
c     

         do i=1,3
            r_disl(i,ndisl) = r0(i)
         enddo
         if(elem_disl(ndisl).gt.0) then
            moved=.false.
            ishift=-nint(cos(th_e))
 1          continue
            do i=1,ndisl-1
               rr=0.d0
               do j=1,2
                  rr=rr+(r_disl(j,i)-r0(j))**2
               enddo
               if(rr.lt.RTOL) then
                  r0(1)=r0(1)+ishift*bv(1)
                  r0(2)=r0(2)+ishift*bv(2)
                  moved=.true.
                  go to 1
               endif
            enddo
            if(moved) then
               elem_disl(ndisl)=fe_locate(r0,elem_disl(ndisl))
               write(*,*) ndisl,' overlap found.  moved from'
               write(*,*) r_disl(1:3,ndisl), 'to'
               write(*,*) r0(1:3)
            endif
         endif
c     
         do i=1,3
            burgers(i,ndisl) = bv(i)
            r_disl(i,ndisl) = r0(i)
         enddo
         burg_length(ndisl) = 
     &        dsqrt( burgers(1,ndisl)*burgers(1,ndisl) + 
     &        burgers(2,ndisl)*burgers(2,ndisl) ) 
         theta_e(ndisl)=th_e
         theta_s(ndisl)=th_s
         print *, 'In element: ', elem_disl(ndisl)
         if(elem_disl(ndisl).ne.0) call sliprange(bv,r0,disl_range(1
     $        ,ndisl),disl_index(ndisl))
c     
      else
         print *, 'Lumping with disl # ', ndisl + n_total  - 4
         do i = 1, 3
            burgers(i, ndisl + n_total - 4) = burgers(i, ndisl+n_total
     $           -4)+ bv(i)
         enddo
         burg_length(ndisl + n_total -4)= sqrt(burgers(1, ndisl+n_total
     $        -4)**2+ burgers(2, ndisl+n_total-4)**2)
         if(n_total.eq.4) then
            n_total = n_total - 4
         endif
      endif
      return
      end



      subroutine disl_pass( r0, rd, burg, th_e, th_s,
     &     x, b, is_relaxed, numnp ,subtract, store
     $     ,idis_slip, islp, s_dis)

      use mod_dd_slip
c     
c     given location x and xd and the burg/theta of a dislocation,
c     compute
c     and return b=b-utilde(x)+utilde(xd)
c     and add this new disl to the d.d. side.
c     
      implicit none
      double precision r0(3), rd(3), burg(3), b(3,*), x(3,*) 
      double precision th_e, th_s
      integer is_relaxed(*), numnp
      logical subtract, store, nucl
c     
      double precision u(3), ud(3)
      integer i, j
      integer, optional :: idis_slip
      integer, optional :: islp
      double precision, optional :: s_dis
c     
      if (present(islp) .and. present(s_dis)) then 
         nucl = .true.
      else
         nucl = .false.
      end if

      print *, 'Image Locations'
      print *, 'Subtract =', r0(1:2)
      print *, 'Image = ', rd(1:2)
c      if (store) call disl_accept(rd, burg, th_e, th_s)
      if (store) then 
         if (nucl) then 
            call nucleate_atomistic_disl(islp, s_dis, th_e, th_s, burg)
         else
            print *, 'Slip plane num and distance not given to pass'
         end if
      end if
      do i=1, numnp
         if(is_relaxed(i).ne.0) then
            if(subtract) then
               call disl_u(r0, burg, th_e, th_s, x(1,i), u)
            else
               do j=1,3
                  u(j)=0.d0
               enddo
            endif
            call disl_u(rd, burg, th_e, th_s, x(1,i), ud)
	    do 1 j=1,3
 1             b(j,i)=b(j,i)-u(j)+ud(j)
            endif
         enddo
!     if (.not. store) then 
         if (present(idis_slip)) then 
            if (idis_slip > 0) then 
               call remove_disl_global(idis_slip)
            end if
         end if
!     endif
         return
         end



      function fd_no_disl()
      implicit none
      include 'disl_parameters.par'
      logical fd_no_disl
c     
      fd_no_disl = (ndisl.eq.0)
      return
      end

      subroutine fd_peach_koeller_nuc(rhs)
      use mod_dd_slip
      implicit none
      include 'disl_parameters.par'
      integer islp, i,j,k,l,jslp, li,ii
      double precision rhs(*), sout(3), s(2)
      print *, 'Cos and Sin and Burger in fd_nuc', blen

      ii = 0
      do islp = 1, nslp
         li = locphi(islp)
         if (nnuc(islp) > 0) then 
            do i = 1, nnuc(islp)
               ii = ii + 1
               sout = 0.0d0
               s = 0.0d0
               if (elem_source(i,islp) > 0) then 
                  call fe_stress(elem_source(i,islp), rhs,
     $                 sout)
               endif
c$$$  do jslp = i, nslp
c$$$  if (nnuc(jslp) > 0) then 
c$$$  if (ndis(jslp) > 0) then 
c$$$  do j = 1, ndis(islp)
c$$$  call disl_s(r_disl(1,j), burgers(1,j)
c$$$  $                            r_nuc(3,i,islp), burgers
c$$$  end do
c$$$  endif
c$$$  endif
c$$$  end do
               s(2) = s(1) + sout(1)*blen*cosphi(li)+
     $              sout(3)*blen*sinphi(li)
               s(1) = s(2) + sout(3)*blen*cosphi(li)+
     $              sout(2)*blen*sinphi(li)
               taui(i,islp) = (s(1)*cosphi(li) +
     $              s(2)*sinphi(li))/blen
               write(*,fmt='(3I7,1X,8(1X,E15.7))')
     $              ii,i,islp,sout(1), sout(2), sout(3),
     $              taui(i,islp), t_FR(i,islp),tnlaps(i,islp)

            end do
         endif
      end do
      end subroutine fd_peach_koeller_nuc

      subroutine fd_peach_koeller(rhs)
      implicit none
      include 'disl_parameters.par'
      double precision rhs(*)
c     
      integer i, j, k
      double precision s_out(3)
      character*80 error_message
c     
      if(i_disl.ne.1) then
         error_message = 'fd_peach_koeller: call disl_setup first!'
         call error_handler(error_message)
      endif
c     
      do i=1, ndisl
         if (elem_disl(i).gt.0) then
	    call fe_stress(elem_disl(i), rhs, pk_stress(1,i))
	    do j = 1, ndisl
               if(j.ne.i) then
                  call disl_s(r_disl(1,j), burgers(1,j),
     &                 r_disl(1,i), s_out,theta_s(j))
                  do k = 1, 3
c$$$  write(*,*) ' pk_stress', pk_stress(k,i)
                     pk_stress(k,i) = pk_stress(k,i) + s_out(k)
                  enddo
               endif
	    enddo
         endif
      enddo

      do i=1, ndisl
         if (elem_disl(i).gt.0) then
	    pk_force(2,i) = -(pk_stress(1,i)*burgers(1,i)+
     &           pk_stress(3,i)*burgers(2,i))
	    pk_force(1,i) = pk_stress(3,i)*burgers(1,i)+
     &           pk_stress(2,i)*burgers(2,i)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     if(.not.NumericalPK) then
            pk_f(i) = pk_force(1,i)*burgers(1,i)
     &           + pk_force(2,i)*burgers(2,i)
            pk_f(i) = pk_f(i) / burg_length(i)
            print *, 'Pk-force =', i, pk_f(i)
c     if(abs(pk_f(i)).lt.PEIERLS) pk_f(i)=0.d0
!     endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         else
	    pk_force(1,i) = 0.0d0
	    pk_force(2,i) = 0.0d0
	    pk_f(i) = 0.0d0
         endif
      enddo
c     
      return
      end



