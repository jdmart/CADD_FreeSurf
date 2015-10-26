c
c $Id: disl_cg.f,v 1.1.1.1 2003/03/12 20:09:00 shastry Exp $
c
        subroutine cg_setup
        implicit none
        include 'disl_parameters.par'
c
        double precision pk_b(max_disl)
        integer i_cg
        common /cg_common/ pk_b, i_cg
        data  i_cg /0/
c
        integer i
        character*80 error_message
c
        if(i_cg.ne.0) then
          error_message = 'cg_setup: improper exit from CG!'
          call error_handler(error_message)
        else
          i_cg = 1
        endif
        if(i_disl.ne.1) then
          error_message = 'cg_setup: call disl_setup first!'
          call error_handler(error_message)
        endif
c
        call cg_reset
c
        return
        end



        subroutine cg_exit
        implicit none
        include 'disl_parameters.par'
c
        double precision pk_b(max_disl)
        integer i_cg
        common /cg_common/ pk_b, i_cg
c
        character*80 error_message
c
        if(i_cg.ne.1) then
          error_message = 'cg_exit: exit from CG without entering!'
          call error_handler(error_message)
        else
          i_cg = 0
        endif
c
        return
        end



        subroutine cg_reset
        implicit none
        include 'disl_parameters.par'
        double precision cg_gdotg
c
        double precision pk_b(max_disl)
        integer i_cg
        common /cg_common/ pk_b, i_cg
c
        integer i
        character*80 error_message
c
        if(i_cg.ne.1) then
          error_message = 'cg_reset: CG without entering!'
          call error_handler(error_message)
        endif
c
        do i=1, ndisl
          pk_b(i) = 0.0d0
        enddo
        return
        end



        function cg_gdotg()
        implicit none
        include 'disl_parameters.par'
        double precision cg_gdotg
c
        double precision pk_b(max_disl)
        integer i_cg
        common /cg_common/ pk_b, i_cg
c
        integer i
        character*80 error_message
c
        if(i_cg.ne.1) then
          error_message = 'cg_gdotg: CG without entering!'
          call error_handler(error_message)
        endif
c
        cg_gdotg = 0.0d0
        do i=1, ndisl
          if (elem_disl(i).gt.0) then
            cg_gdotg = cg_gdotg + pk_f(i) * pk_f(i)
          endif
        enddo
        return
        end



        subroutine cg_z(z,cgsmax)
        implicit none
        include 'disl_parameters.par'
        double precision z, cgsmax
c
        double precision pk_b(max_disl)
        integer i_cg
        common /cg_common/ pk_b, i_cg
c
        integer i
        character*80 error_message
c
        if(i_cg.ne.1) then
          error_message = 'cg_z: CG without entering!'
          call error_handler(error_message)
        endif
c
        cgsmax = 0.0d0
        do i=1, ndisl
          if(elem_disl(i).gt.0) then
            pk_b(i) = pk_b(i)*burg_length(i) * z + pk_f(i)
            if(dabs(pk_b(i)).gt.cgsmax) then
              cgsmax=dabs(pk_b(i))
            endif
            pk_b(i) = pk_b(i) / burg_length(i)
          else
            pk_b(i) = 0.0d0
          endif
        enddo
        return
        end



        subroutine cg_saxpy(alpha)
        implicit none
        include 'disl_parameters.par'
        double precision alpha
c
        double precision pk_b(max_disl)
        integer i_cg
        common /cg_common/ pk_b, i_cg
c
        integer fe_locate, i, elem_old
        double precision b
        character*80 error_message
c
        if(i_cg.ne.1) then
          error_message = 'cg_saxpy: CG without entering!'
          call error_handler(error_message)
        endif
c
        do i = 1, ndisl
          if(elem_disl(i).ne.0) then
            r_disl(1, i) = r_disl(1,i) + alpha * pk_b(i) * burgers(1,i)
            r_disl(2, i) = r_disl(2,i) + alpha * pk_b(i) * burgers(2,i)
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
************************************************************************
        double precision function cg_gdots()
        implicit none
        include 'disl_parameters.par'
        double precision alpha
c
        double precision pk_b(max_disl)
        integer i_cg
        common /cg_common/ pk_b, i_cg
c
        integer i
        character*80 error_message
c
        if(i_cg.ne.1) then
          error_message = 'cg_gdots: CG without entering!'
          call error_handler(error_message)
        endif
c
        cg_gdots=0.d0
        do i=1, ndisl
          if(elem_disl(i).gt.0) then
            cg_gdots = cg_gdots - pk_b(i)*burg_length(i) * pk_f(i)
          endif
        enddo
        return
        end



        double precision function cg_sdots()
        implicit none
        include 'disl_parameters.par'
c
        double precision pk_b(max_disl)
        integer i_cg
        common /cg_common/ pk_b, i_cg
c
        integer i
        character*80 error_message
c
        if(i_cg.ne.1) then
          error_message = 'cg_sdots: CG without entering!'
          call error_handler(error_message)
        endif
c
        cg_sdots=0.d0
        do i=1, ndisl
          if(elem_disl(i).gt.0) then
            cg_sdots = cg_sdots +
     &               pk_b(i)*burg_length(i) * pk_b(i)*burg_length(i)
          endif
        enddo
        return
        end


c
c $Log: disl_cg.f,v $
c Revision 1.1.1.1  2003/03/12 20:09:00  shastry
c vijay-   Initial import.
c
c Revision 1.5  2002/03/05 03:00:45  shilkrot
c Moved pk_b out of dislocations.par into the cg common block.
c
c Revision 1.4  2002/02/28 22:43:08  shilkrot
c Removed the condition elem_disl != 0 for moving dislocations in cg_saxpy.
c
c Revision 1.3  2002/02/25 20:54:05  shilkrot
c Added a routine disl_reset.
c Changed the code to use burg_length.
c
c Revision 1.2  2001/12/13 07:31:24  shilkrot
c Implemented breadth first search to find the element number for
c a dislocation. Changed the interface of fe_locate to include the starting
c element for the search. Old fe_locate is in fem_services.
c Changed the interface of fem_setup. Now two arrays used as temp space are
c passed from outside as the last two parameters.
c
c Revision 1.1  2001/11/13 03:50:37  shilkrot
c Routines computing scalar products and moving dislocations for the
c conjugate gradient.
c
c
