c**---------------------------------------------------------------
c**   rest : read/write the restart file
c**
      subroutine rest(id,x,ix,itx,f,b,logic,key)
      use mod_global
      use mod_boundary
      use mod_grain
      implicit none
      include '../Disl/disl_parameters.par'
c
      double precision x(*),f(*),b(*)
      integer id(*),ix(*),logic,itx(*)
      character*4 key
      logical RepsUpToDate_Dummy
c
c---- read/write restart files
c
      character*4 xkey1,xkey2
      double precision xcrit
      common /ccrit/ xcrit(2)
      logical boundary
      integer icrit,isign,nodecrit,idfcrit
      double precision arcrate
      common /cstep/ boundary,icrit,arcrate,isign,nodecrit,idfcrit
      integer neqtot,i,j,k,nsdim,nxsjdim,niddim
     $     ,nixdim,ntdim,njddim,ng,nshpdim,nxdim,nfdim,nby

      ! get number of grains
c
      xkey1 = 'read'
      xkey2 = 'writ'
      if (key.eq.xkey1) then
      rewind logic
      read (logic) numnp,numel,neq
      if (numnp.gt.maxnp) then
      write(6,'(///''  **error** detected by subroutine rest''
     1 /'' number of nodes exceeds maximum'')')
      write(6,'(''  **required = '',i10)') numnp
      write(6,'(''  **maximum  = '',i10)') maxnp
      stop
      end if
      if (numel.gt.maxel) then
      write(6,'(///''  **error** detected by subroutine rest''
     1 /'' number of elements exceeds maximum'')')
      write(6,'(''  **required = '',i10)') numel
      write(6,'(''  **maximum  = '',i10)') maxel
      stop
      end if
      read (logic) time,dt,time,timeol
      read (logic) boundary,icrit,isign,nodecrit,idfcrit,
     1 xcrit(1),xcrit(2)
      neqtot = neq + nstad*numel
      read (logic) (b(i),i=1,neqtot)
      nsdim = nsdm*nquad*numel
      nshpdim = nshpdm*nquad*numel
      nxsjdim = nquad*numel
      niddim = ndf*numnp
      nxdim  = nxdm*numnp
      nixdim = nen1*numel
      nfdim  = ndf*numnp
      ntdim  = numnp
      njddim = ndf*numnp
      read (logic) (id(i),i=1,niddim)
      read (logic) (x(i),i=1,nxdim)
      read (logic) (ix(i),i=1,nixdim)
      read (logic) (itx(i),i=1,3*numel)
      read (logic) (f(i),i=1,nfdim)
      read (logic) nce,ncb
      if(nce.gt.NCEMAX) call IncreaseElist(nce-NCEMAX+100)
      read (logic) ((elist(i,j),i=1,2),j=1,nce)
c--Representative atom data
      read (logic) ng
      if (ng.ne.ngrains) then
         write(6,'(///''  **error** detected by subroutine rest''
     1    /'' number of grains in file different than current value'')')
         write(6,'(''  **in file  = '',i10)') ng
         write(6,'(''  **current  = '',i10)') ngrains
         stop
      endif
      read (logic) RepsUpToDate_Dummy
      read (logic) (energy(j),j=1,numnp)
      read (logic) (IsRelaxed(j),j=1,numnp)
c--   dislocations
      read(logic) ndisl
      if(ndisl.gt.max_disl) stop 'increase max_disl'
      read(logic) burgers(1:3,1:ndisl),r_disl(1:3,1:ndisl)
      read(logic) theta_e(1:ndisl),theta_s(1:ndisl)
      read(logic) burg_length(1:ndisl)
      read(logic) elem_disl(1:ndisl)
      do i=1,ndisl
         if(elem_disl(i).ne.0) then
            call sliprange(burgers(1:3,i),r_disl(1:3,i),disl_range(1:2,i
     $           ),disl_index(i))
         endif
      enddo
      end if

      if (key.eq.xkey2) then
      rewind logic
      write (logic) numnp,numel,neq
      write (logic) time,dt,time,timeol
      write (logic) boundary,icrit,isign,nodecrit,idfcrit,
     1 xcrit(1),xcrit(2)
      neqtot = neq + nstad*numel
      write (logic) (b(i),i=1,neqtot)
      nsdim = nsdm*nquad*numel
      nshpdim = nshpdm*nquad*numel
      nxsjdim = nquad*numel
      niddim = ndf*numnp
      nxdim  = nxdm*numnp
      nixdim = nen1*numel
      nfdim  = ndf*numnp
      ntdim  = numnp
      njddim = ndf*numnp
      write (logic) (id(i),i=1,niddim)
      write (logic) (x(i),i=1,nxdim)
      write (logic) (ix(i),i=1,nixdim)
      write (logic) (itx(i),i=1,3*numel)
      write (logic) (f(i),i=1,nfdim)
      write (logic) nce,ncb
      write (logic) ((elist(i,j),i=1,2),j=1,nce)
      write (logic) ngrains
c--Representative atom data
      write (logic) RepsUpToDate_Dummy
      write (logic) (energy(j),j=1,numnp)
      write (logic) (IsRelaxed(j),j=1,numnp)
c--   dislocations
      write(logic) ndisl
      write(logic) burgers(1:3,1:ndisl),r_disl(1:3,1:ndisl)
      write(logic) theta_e(1:ndisl),theta_s(1:ndisl)
      write(logic) burg_length(1:ndisl)
      write(logic) elem_disl(1:ndisl)
      end if
      return
      end


