c
c $Id: fem_service.f,v 1.1.1.1 2003/03/12 20:09:00 shastry Exp $
c
        subroutine fe_plot(x,b,ix,rhs,numel,avevirst)
        implicit none
        include 'fem_parameters.par'
c
        double precision x(3,*),b(3,*),rhs(*),avevirst(3,3,*)
        integer ix(4,*),numel
        double precision xl(knode*ndof), xu(knode*ndof), s_out(3)
        integer i, j, k
        character*80 filename
        character*3 cnt
        integer, save::icount1
        data icount1/0/
        write(cnt,'(i3)')icount1
        if(icount1.lt.10)cnt(1:2)='00'
        if(icount1.lt.100)cnt(1:1)='0'
        filename(1:7)='Stress_'
        filename(8:10)=cnt
        filename(11:14)='.plt'
c
        open(903, file=filename)
        write (903,*)
     &      'ZONE N=',3*numel,' E=', numel, ' , F=FEPOINT, ET=TRIANGLE'
        do i=1, nelm
          call fe_find(i, xl)
          call fe_find_rhs(i, xu, rhs)
          call fe_stress(i, rhs, s_out)
          do j=1, knode
            write (903,'(7e16.8)') (xl(ndof*(j-1)+k), k=1,ndof),
     &                  (xl(ndof*(j-1)+k)+xu(ndof*(j-1)+k), k=1,ndof),
     &                     (s_out(k), k=1,3)
          enddo
        enddo
        do i=1,numel
          if(ix(4,i).ne.0)then
            do j=1,knode
               write(903,'(7e16.8)')(x(k,ix(j,i)),k=1,ndof),
     &                    (x(k,ix(j,i))+b(k,ix(j,i)),k=1,ndof),
     &                    avevirst(1,1,ix(j,i)),avevirst(2,2,ix(j,i)),
     &                    avevirst(1,2,ix(j,i))
            enddo
          endif
        enddo
        write (903,*)
        do i=1, numel
          write (903,*) (knode*(i-1)+k, k=1, knode)
        enddo
        close(903)
        icount1=icount1+1
c
        return
        end


        subroutine fe_pk_print()
        implicit none
        include 'disl_parameters.par'
        double precision pk
        integer i, j
c
        do i=1, ndisl
          if(elem_disl(i).gt.0) then
            print *, 'Disl: ', i
            print *, 'Element: ', i, elem_disl(i)
            print *, 'Stress: ', i, (pk_stress(j,i), j=1,3)
            pk =  pk_force(1,i)*burgers(1,i)
     &            + pk_force(2,i)*burgers(2,i)
            pk = pk / burg_length(i)
            print *, 'Force: ', i, (pk_force(j,i), j=1,2), pk
            print *, 'Gliding force: ', i, pk_f(i)
          endif
        enddo
c
        return
        end



        subroutine disl_print(iflag)
        implicit none
        include 'disl_parameters.par'
        integer i, j, iflag
c
        open(unit=401,file='dislocations.plt',status='unknown')
        write(401,*) 'zone'
        do i=1, ndisl
           write(401,'(5e15.6,2i5)') (r_disl(j,i), j=1,2),
     $                           (burgers(j,i),j=1,3),i,iflag
        enddo
        call flush(401)
c
        return
        end



        subroutine fd_move(id, dr)
        implicit none
        include 'disl_parameters.par'
        integer id, fe_locate
        double precision dr(2)
c
        integer elem_old
c
        r_disl(1,id) = r_disl(1,id)+dr(1)
        r_disl(2,id) = r_disl(2,id)+dr(2)
        elem_old = elem_disl(id)
        elem_disl(id)=fe_locate(r_disl(1,id), elem_disl(id))
        if(elem_disl(id).eq.0) then
          if(elem_old.gt.0) then
            elem_disl(id) = -elem_old
          else
            elem_disl(id) = elem_old
          endif
        endif
        return
        end


        subroutine fe_stress_check(rhs,e0)
        implicit none
        include 'fem_parameters.par'
        double precision rhs(*)
        double precision e(3), s(3), e0, area, fe_area
        integer i
c
        e0 = 0.0d0
        do i=1, nelm
          call fe_strain(i, rhs, e)
          call fe_stress(i, rhs, s)
          area = fe_area(i)
          e0 = e0 +(s(1)*e(1)+s(2)*e(2)+s(3)*e(3))*area
        enddo
        e0 = e0/2.d0
C        print *, 'stress_check: ', e0
        return
        end



        function fe_area(lmn)
        implicit none
        include 'fem_parameters.par'
        integer lmn
        double precision fe_area
c
        double precision xl(knode*ndof)
        double precision pn(knode), qn(knode)
        double precision xjac(ndof, ndof), det, area
        integer i, j
c
c Everything is hard-coded, since I don't know where I may need
c the derivatives of shape functions
c
c Shape functions are: N_1 = 1 - p - q; N_2 = p; N_3 = q;
c
        pn(1)=-1.0d0
        pn(2)=1.0d0
        pn(3)=0.0d0
        qn(1)=-1.0d0
        qn(2)=0.0d0
        qn(3)=1.0d0
c
        call fe_find(lmn, xl)
        do 2 i=1,ndof
        do 2 j=1,ndof
2         xjac(i,j)=0.0d0
c
        do i=1, knode
          xjac(1,1)=xjac(1,1)+pn(i)*xl((i-1)*ndof+1)
          xjac(1,2)=xjac(1,2)+qn(i)*xl((i-1)*ndof+1)
          xjac(2,1)=xjac(2,1)+pn(i)*xl((i-1)*ndof+2)
          xjac(2,2)=xjac(2,2)+qn(i)*xl((i-1)*ndof+2)
        enddo
c
        det = xjac(1,1)*xjac(2,2) - xjac(1,2)*xjac(2,1)
        fe_area = abs(det)/2.0d0
c
        return
        end



        function fe_locate_old(r,idummy)
c
c  Returns the element number corresponding to position r or zero.
c  Needs to be replaced by something smarter than this.
c
        implicit none
        include 'fem_parameters.par'
        integer fe_locate_old, idummy
        double precision r(3)
c
        integer i, i1, i2, i3
        logical fe_in_tri
        character*80 error_message
c
        if(i_flag.ne.1) then
          error_message = 'fe_locate_old: call fem_setup first!'
          call error_handler(error_message)
        endif
c
        fe_locate_old = 0
        do i=1, nelm
          i1=iconn(1,i)
          i2=iconn(2,i)
          i3=iconn(3,i)
          if(fe_in_tri(x0(1,i1), x0(1,i2), x0(1,i3), r)) then
            fe_locate_old = i
            return
          endif
        enddo
        return
        end


c
c $Log: fem_service.f,v $
c Revision 1.1.1.1  2003/03/12 20:09:00  shastry
c vijay-   Initial import.
c
c Revision 1.2  2001/12/13 07:31:24  shilkrot
c Implemented breadth first search to find the element number for
c a dislocation. Changed the interface of fe_locate to include the starting
c element for the search. Old fe_locate is in fem_services.
c Changed the interface of fem_setup. Now two arrays used as temp space are
c passed from outside as the last two parameters.
c
c Revision 1.1  2001/11/13 04:30:36  shilkrot
c Service routines like plotting, printing etc.
c
