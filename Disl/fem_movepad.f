c
c $Id: fem_movepad.f,v 1.1.1.1 2003/03/12 20:09:00 shastry Exp $
c
	subroutine fem_move_pad(x,b,ix,bc,prop,cz,id,is_relaxed,e_out
	1    ,f_out,FullField, MoveDisl,strainE0,numel,avevirst,MDTemp
	1    ,iFem, Moved)
c
c123456789012345678901234567890123456789012345678901234567890123456789012
c         1         2         3         4         5         6         7
c
c if the node is on the interface u=b
c if the node is on the boundary u=prop*bc
c changes b 
c e_out for informational purposes only
c
c************************************************************************
c       Comment ---- here prop=time variable which is really 
c       the load increment and basically the total displacments/ forces 
c       are scaled with prop in all the fe_routines, but really time 
c       from the md/input files
c       Correllation with md parameters 
c       cz --- z_length ... length of md box in the z direction
c       prop = time 
c       f_out ---- dr (forces)
c       e_out ---- eadd
c       x, b, ix  are the same parameters 
c       f---- bc  boundary conditions (id=1, displacement, id=0, force)
c       db---- displacement increment/step
c************************************************************************
	use mod_dd_slip
	implicit none
	include 'fem_parameters.par'
	include 'disl_parameters.par'

	double precision x(3,*), b(3,*), bc(3,*), f_out(3, *)
	double precision avevirst(3,3,*), tincr1
	integer numel
	integer id(3,*), ix(4,*),is_relaxed(*)
	double precision prop, cz
CC--Jun Song: MDTemp=SysTemp defined in mod_global
	double precision e_out,strainE0,MDTemp
	logical FullField, MoveDisl
	logical Moved
c
	double precision rhs(maxeqs), forces(maxeqs), e0, eplus, eminus
	integer i,iFem
	integer fe_locate, elem_old, elem_plus, elem_minus
	character*80 error_message
	logical CentralDiff,ForwardDiff,BackDiff
	logical ChangeTime

!
!      print *, 'Ifem in fem_move_pad is', iFem
	if(i_flag.ne.1) then
	  error_message = 'fem_solve: call fem_setup first!'
	  call error_handler(error_message)
	endif
	if(i_disl.ne.1) then
	  error_message = 'fem_solve: call disl_setup first!'
	  call error_handler(error_message)
	endif
!
!
!
!       solve FEM and DD problem
	tincr = tincr_sav
	if (Moved) then 
	   tincr = 1.d-20
	end if
	call assign_disloc_global
	print *, 'Total dislocations', ndisl
	print *, 'Time entering fd_solve is', prop
	call fd_solve(b, bc, prop, cz, id, is_relaxed, rhs, forces, e0)
!
!  compute the P.-K. force on dislocations
	if(MoveDisl) then
	   if (ndisl > 0) then 
	      if (iFem .eq. 1) then 
		 print *,' --- Entering Move Disl --'
                 ! call fd_peach_koeller(rhs)
		 call move_disloc(rhs)
		 call nucleate_disl(rhs)
		 dtime_dd = dtime_dd + tincr
		 print *, 'Dtime_dd =', dtime_dd
	      else
		 tincr = tincr_sav/100.0d0
		 print *,' --- Entering Move Disl --'
                 ! call fd_peach_koeller(rhs)
		 call move_disloc(rhs)
		 call nucleate_disl(rhs)
		 dtime_dd = dtime_dd + tincr
		 print *, 'Dtime_dd =', dtime_dd
		 
	      endif
	   else
	      if (iFem .eq. 1) then 
		 call move_disloc(rhs)
		 call nucleate_disl(rhs)
		 dtime_dd = dtime_dd + tincr
		 print *, 'Dtime_dd =', dtime_dd	      
	      endif
	   endif
	endif
!       
	do i=1,nfixed
	  f_out(ifix_hold(1,i), imap(ifix_hold(2,i))) = 
     &      f_out(ifix_hold(1,i), imap(ifix_hold(2,i)))
     &      -forces(ifixed(i))*cz
	enddo
	e_out = e0*cz
!
!       move dislocations based on P.-K. force
!        if(MoveDisl) call move_dis(10.0,MDTemp)
!
!       move pad atoms according to tilda and hat fields
	call fd_movepad(x, rhs, b)
!
!       update b so b=u_hat+u_tilda
	if(FullField) then
	  call fd_full_field(rhs, b)
	endif
!
!       compute total strain energy of fem region
        call fe_stress_check(rhs,strainE0)

	return
	end



	subroutine fd_movepad(x, rhs, b)
	implicit none
	include 'fem_parameters.par'
	double precision x(3, *), rhs(*), b(3,*)
c
	double precision u_out(3)
	integer i, j, k, lmn, n
c
	do i = 1, npad
	  n = padmap(i)
	  lmn = padelmnt(i)
	  do j = 1, ndof
	    b(j,n) = 0.d0
	    do k = 1, knode
	      b(j,n) = b(j,n) + 
     &          rhs( (iconn(k,lmn)-1)*ndof+j ) * padtricoord(k,i)
	    enddo
	  enddo
	  call disl_displ(x(1, n), u_out)
	  do j = 1, ndof
	    b(j, n) = b(j, n) + u_out(j)
	  enddo
	enddo
	return
	end
	


	subroutine fe_tricoord(r1, r2, r3, p, coord)
c
c  Yet another wrapper for intri from feaplib
c
	implicit none
	double precision r1(3), r2(3), r3(3), p(3), coord(3)
c
	double precision x1(2), x2(2), x3(2), pt(2)
	logical intri, ontri
	integer i
	character*80 error_message
c
	do i = 1, 2
	  x1(i) = r1(i)
	  x2(i) = r2(i)
	  x3(i) = r3(i)
	  pt(i) = p(i)
	enddo
	if(.not.intri(x1,x2,x3,pt,coord,ontri)) then
	  error_message = 'fe_tricoord: the atom is not in the element'
	  call error_handler(error_message)
	endif
	return
	end

c
c $Log: fem_movepad.f,v $
c Revision 1.1.1.1  2003/03/12 20:09:00  shastry
c vijay-   Initial import.
c
c




!!!!!!!!dw added subroutine!!!!!!!!!!!!!!!!!!
        subroutine move_dis(alpha,temperature)
        implicit none
        include 'disl_parameters.par'
        double precision alpha, mobility, max_vel
        double precision sf_f, aa, bb
        double precision min_pos, temperature, time_step_con
c
        integer fe_locate, i, elem_old
        double precision b
        double precision ddis
        character*80 error_message
c
!!!!    hacked parameters
        min_pos=-10.0
        time_step_con=5.0e-4
CC--Jun Song: make sure temperature>0.0
	if(temperature .lt. 0.0d0) then
	  write(*,*)"Temperature less than Zero!!!"
	  stop
	endif

CC--JS: 6.242e-2 is unit conversion constant. Do not change
CC--Change the stacking fault E for different materials 
        mobility=time_step_con/(6.242e-2*5.0e-8*temperature) 
        ! 3rd parameter (5.0e-8) is damp coef from Olmstead paper
        ! "Atomistic simulations of dislocation mobility.."
        max_vel=time_step_con*2000.0
        sf_f=.089*6.242e-2 ! sf energy in J/m2,i.e., 0.089
!!!!    end of hacked parameters

        do i = 1, ndisl
          if(elem_disl(i).ne.0) then
            pk_f(i) = (pk_force(1,i)*burgers(1,i)
     &                + pk_force(2,i)*burgers(2,i))/burg_length(i)
!            write(*,*) i,' non sf disl force = ',pk_f(i)
            ddis=sqrt(r_disl(1,i)**2+r_disl(2,i)**2)
            pk_f(i) = pk_f(i) + sf_f
!!!!        upper limit on velocities
            if (abs(pk_f(i)*mobility).gt.max_vel) then
               pk_f(i)=max_vel/mobility*pk_f(i)/abs(pk_f(i))
            endif
!            write(*,*) 'old disl pos = ',r_disl(1, i), r_disl(2, i),ddis
!            write(*,*) i,' total disl force = ',pk_f(i)
!            if(r_disl(2,i).gt.-100.0.or.pk_f(i).gt.0.0) then
            if(r_disl(2,i).le.min_pos.or.pk_f(i).lt.0.0) then
             r_disl(1, i) = r_disl(1,i) + 
     &         mobility*pk_f(i)*burgers(1,i)/burg_length(i)
             r_disl(2, i) = r_disl(2,i) + 
     &         mobility*pk_f(i)*burgers(2,i)/burg_length(i)
            endif
!            endif
            write(*,*) 'new disl pos = ',r_disl(1, i), r_disl(2, i)
            write(*,*) ' '
            elem_old = elem_disl(i)
            elem_disl(i)=fe_locate(r_disl(1,i), elem_disl(i))
            if(elem_disl(i).eq.0) then
              if(elem_old.gt.0) then
                elem_disl(i) = -elem_old
              else
                elem_disl(i) = elem_old
              endif
            endif
          endif
        enddo
        return
        end
