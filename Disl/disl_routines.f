c
c $Id: disl_routines.f,v 1.1.1.1 2003/03/12 20:09:00 shastry Exp $
c
        subroutine fd_update_u_tilde_bc(u_bc)
        use mod_dd_slip
c
c123456789012345678901234567890123456789012345678901234567890123456789012
c         1         2         3         4         5         6         7
c
        implicit none
        include 'fem_parameters.par'
        double precision u_bc(*)
c
        double precision u_out(3), u11(3)
        character*80 error_message
        integer i
c
        if(i_flag.ne.1) then
          error_message = 'fd_update_u_tilde_bc: call fem_setup first!'
          call error_handler(error_message)
        endif
c
        do i=1,nfixed
          call disl_displ(x0(1, ifix_hold(2,i)), u_out)
!          call dislp(x0(1, ifix_hold(2,i)), u11)
!          write(*, fmt='(A6,I7,1X,4(1X,E15.8))') 'dislp ',i,
!     $         u11(1),u_out(1)
!     $         , u11(2), u_out(2)
          u_bc(i) = u_out(ifix_hold(1,i))
!          print *, 'Bc = ', i, u_bc(i)
        enddo
        return
        end



        subroutine fd_update_f_tilde_bc(f_segm)
        implicit none
        integer ninp
        parameter (ninp=2)
        include 'fem_parameters.par'
        double precision f_segm(ndof,*)
c
        double precision s(3), s_out(3), xnorm(ndof), x(ndof+1)
        double precision xinp(ninp), w(ninp)
        integer i, j, k
        character*80 error_message
        data xinp /-0.577350269189625731d0, 0.577350269189625731d0/
        data w /1.d0, 1.d0/
c
        if(i_flag.ne.1) then
          error_message = 'fd_update_f_tilde_bc: call fem_setup first!'
          call error_handler(error_message)
        endif
c
        do i =1, nsegm
          xnorm(1)=x0(2, isegm(2,i))-x0(2,isegm(1,i))
          xnorm(2)=-(x0(1, isegm(2,i))-x0(1,isegm(1,i)))
c
c  Integral over the segment using 2 Gaussian integration points
c
          do 1 j=1,3
 1          s(j)=0.d0
c
          do j=1, ninp
            do k=1,ndof
              x(k) = (1.d0+xinp(j))/2.0d0*x0(k, isegm(1,i))
     &             + (1.d0-xinp(j))/2.0d0*x0(k, isegm(2,i))
            enddo
            x(ndof+1) = (1.d0+xinp(j))/2.0d0*x0(ndof+1, isegm(1,i))
     &                + (1.d0-xinp(j))/2.0d0*x0(ndof+1, isegm(2,i))
c
c            print *, 'Nsegm', i, x(1),x(2)
            call disl_stress(x, s_out)
            do 2 k=1,3
 2            s(k)=s(k)+s_out(k)*w(j)/2.0d0
          enddo
c
          f_segm(1,i)=s(1)*xnorm(1)+s(3)*xnorm(2)
          f_segm(2,i)=s(2)*xnorm(2)+s(3)*xnorm(1)
        enddo
        return
        end



        subroutine fd_corrective_force_bc(rhs, f_segm)
        implicit none
        include 'fem_parameters.par'
        double precision rhs(*), f_segm(ndof, *)
c
        integer i, j, k
c
c Puts a force on constrained nodes. This is overwritten by
c fe_substitute()
c
        do 1 i=1, nsegm
        do 1 j=1,2
        do 1 k=1, ndof
 1      rhs((isegm(j,i)-1)*ndof+k)=
     &        rhs((isegm(j,i)-1)*ndof+k)-f_segm(k,i)/2.0d0
        return
        end



        subroutine fd_augment_force_energy(rhs, forces, f_segm, e0)
        implicit none
        include 'fem_parameters.par'
        double precision rhs(*), forces(*), f_segm(ndof,*), e0
c
        integer i, j, k
c
c  f~ * u^ term
c
        do i=1, nsegm
          do j=1, 2
            do k=1, ndof
              forces((isegm(j,i)-1)*ndof+k) =
     &          forces((isegm(j,i)-1)*ndof+k) + f_segm(k,i)/2.0d0
              e0=e0 +
     &          f_segm(k,i)/2.0d0*rhs((isegm(j,i)-1)*ndof+k)
            enddo
          enddo
        enddo
c
        return
        end



        subroutine fd_corrective_displ_bc(presv, u_bc)
        implicit none
        include 'fem_parameters.par'
        double precision presv(*), u_bc(*)
c
        integer i
c
        do 1 i=1,nfixed
 1      presv(i) = presv(i)-u_bc(i)
c
        return
        end



        subroutine fd_restore_displ_bc(rhs, u_bc)
        implicit none
        include 'fem_parameters.par'
        double precision rhs(*), u_bc(*)
c
        integer i, ic
c
        do i=1,nfixed
          ic = ifixed(i)
          rhs(ic) = rhs(ic) + u_bc(i)
        enddo
        return
        end



        subroutine fd_full_field(rhs, b)
        implicit none
        include 'fem_parameters.par'
        double precision rhs(*), b(3,*)
c
        double precision u(3)
        integer i, j
c
        do i=1,nnodes
          call disl_displ(x0(1, i), u)
          do j=1,ndof
	    rhs((i-1)*ndof+j)=rhs((i-1)*ndof+j)+u(j)
            b(j,imap(i))=rhs((i-1)*ndof+j)
          enddo
            b(ndof+1, imap(i))=u(3)
        enddo
        return
        end



        subroutine fd_rescale_peach_koeller(cz)
        implicit none
        include 'disl_parameters.par'
        double precision cz
        integer i, j
c
        do i=1, ndisl
          do j=1, 2
            pk_force(j, i) = pk_force(j, i) * cz
          enddo
            pk_f(i) = pk_f(i) * cz
        enddo
c
        return
        end


c
c $Log: disl_routines.f,v $
c Revision 1.1.1.1  2003/03/12 20:09:00  shastry
c vijay-   Initial import.
c
c Revision 1.3  2001/11/13 03:34:47  shilkrot
c Added the routines that compute the P.-K. force.
c
c Revision 1.2  2001/08/22 03:18:35  shilkrot
c Fixed the expression for the energy and polished fem_alan a little bit.
c This wersion works with dislocation passing.
c
c Revision 1.1  2001/07/12 06:30:20  shilkrot
c The routines used to apply the tilde field and to handle the array of
c dislocations.
c
c
