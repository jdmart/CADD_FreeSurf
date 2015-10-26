c**---------------------------------------------------------------
c**     ma02  :  produces tecplot read-able files of the mesh and
c**              various contour variables.
c**
c**   Calling format:
c**   ma02,key,filename,index,scale,umag, nprint
c**
c**   key - type of plot - I only trust 'bcon', 'disp', 'disl' (see below)
c**   filename - name of output file (routine adds a numerical counter)
c**   index - not used
c**   scale - magnifies the output variable (scale=1.0 for no magnification)
c**   umag  - exagerates the displacements on the deformed mesh
c**           (umag=1.0 for true deformed mesh). Read from input
c**	      by the subroutine freein()
c**
c**   all above comments are pre 9/1/06
c**
c**   things have been updated 9/28/06 now key='stra' plots strain values
c**   and k='atom' plots an atom eye file
c**
c**   thing have been updated 10/4/06 so that key='stre' plots stress
c**   values 
c**
c**
c--
      subroutine ma02(id,x,ix,f,b,dr,db,input,flag)
      use mod_file
      use mod_dynamo
      use mod_global
      implicit none
      include '../Disl/disl_parameters.par'
c
      double precision x(nxdm,*),f(ndf,*),b(ndf,*),dr(*),db(*)
      integer ix(nen1,*),id(ndf,*)
c
c---- tecplot file
c
      character*80 input,filename,temp, fname_d1, fname_d2
      character*4 key
      character*6 cnt
      logical flag
cc--JS: update icount entries if adding new ikey      
      integer, save :: icount(11)=0, nprint
      integer i,lower,upper,next,ikey,j,logic,index,idum
      double precision dum,scale,umag
c
      key = input(1:4)
      lower = 4
      upper = next(lower,input)
      filename((lower-3):(upper-5)) = input((lower+1):(upper-1))
      do 10 i = upper-4,80
   10 filename(i:i) = ' '
c
      ikey=0
      if (key.eq.'noda') ikey=1
      if (key.eq.'disp') ikey=2 ! plot of displacements
      if (key.eq.'bcon') ikey=3 ! plot of boundary conditions
      if (key.eq.'ener') ikey=4
      if (key.eq.'disl') ikey=5 ! plot of discrete dislocations
      if (key.eq.'stre') ikey=6 ! plot of stresses is now working
      if (key.eq.'stra') ikey=7 ! plot of strains this is now working
      if (key.eq.'viri') ikey=8 ! plot of averaged virial stresses
      if (key.eq.'atom') ikey=9 ! plot of atoms to be read with atomeye
CC--Jun Song modifications
      if (key.eq.'pdbf') ikey=10 ! plot of atoms into pdb file
CC--Jun Song Mod ends
      if (key.eq.'virH') ikey=11 ! plot the viri stress for H atoms
      if (key.eq.'ovit') ikey=12 ! Write Cfg file for ovito
      if(ikey.eq.0) then
         write(*,*) 'ERROR: unknown key'
         return
      endif
c
      icount(ikey)=icount(ikey)+1
      write(cnt,'(i6)') icount(ikey)
      if (icount(ikey).lt.10) cnt(1:5)='00000'
      if (icount(ikey).lt.100)cnt(1:4)='0000'
      if (icount(ikey).lt.1000)cnt(1:3)='000'
      if (icount(ikey).lt.1000)cnt(1:2)='00'
      if (icount(ikey).lt.10000)cnt(1:1)='0'

      i=upper-4
      do j=1,upper-5
         if (filename(j:j).eq.'.') i=j
      enddo
      temp=filename(i:)
      filename(i:i+5)=cnt
      filename(i+6:)=temp

!      if (ndisl > 0 ) then 
c$$$       if (ikey .ne. 5) then 
         call iofile(filename,'formatted  ',logic,.true.)
c$$$      else
c$$$         if (ndisl > 0) then 
c$$$            call iofile(filename,'formatted  ',logic,.true.)
c$$$         endif
c$$$      endif

      lower = upper
      upper = next(lower,input)
      call freein(input,lower,upper,index,dum,1)
      lower = upper
      upper = next(lower,input)
      call freein(input,lower,upper,idum,scale,2)
      lower = upper
      upper = next(lower,input)
      call freein(input,lower,upper,idum,umag,2)
      lower = upper
c$$$      upper = next(lower,input)
c$$$      call freein(input,lower,upper,idum,nprint,2)

      

      print *, 'NPRINT = ', nprint, umag, scale
c
c
      if(ikey.eq.5) then
         print *, 'Total no. of dislocations = ', ndisl
c        if (ndisl > 0 .or. nprint .eq. 0) then 
c         if (ndisl > 0) then 
            call dumpdisl_vtk(scale,logic,umag,b,ndf)
c         else
c           print *, 'no dislocations to print in output file'
c        endif
      else
         if ( ikey < 6 .or. ikey > 8) then 
            call dumpit(x,ix,b,db,id,f,dr,scale,logic,key,index,umag)
         end if
         if (ikey >= 6 .or. ikey <=8) then 
            call dumpit_vtk(x,ix,b,db,id,f,dr,scale,logic,key,index
     $           ,umag)
         end if
      endif
      close(logic)
      return
      end

*******************************************************************
      subroutine GetStresses(strain,stress)
      implicit none
      double precision strain(3,3), stress(3,3)
      double precision vstrain(6), vstress(6)
      double precision t1,t2,theta,pi
      integer i,j

      character*80 error_message
      double precision cc(6,6)
      double precision xe, xnu, xlambda, xmu
      integer i_elas      
      common /elastic/ xe, xnu, xlambda, xmu, cc, i_elas
      double precision ev_convert
      ev_convert = 1.602176462

      if(i_elas.ne.1) then
          error_message = 'fe_dmat: call fe_elastic first!'
          call error_handler(error_message)
      endif



      vstrain(1)=strain(1,1)
      vstrain(2)=strain(2,2)
      vstrain(3)=strain(3,3) 
      vstrain(4)=strain(2,3)+strain(3,2)
      vstrain(5)=strain(1,3)+strain(3,1)
      vstrain(6)=strain(1,2)+strain(2,1)

      do i=1,6
        vstress(i)=0.0
        do j=1,6
!          vstress(i)=vstress(i)+cc(i,j)*vstrain(j)*100e9/1e6
          vstress(i)=(vstress(i)+cc(i,j)*vstrain(j))
        enddo
      enddo

      stress(1,1)=vstress(1)
      stress(2,2)=vstress(2)
      stress(3,3)=vstress(3)
      stress(2,3)=vstress(4)
      stress(1,3)=vstress(5)
      stress(1,2)=vstress(6)
      stress(3,2)=vstress(4)
      stress(3,1)=vstress(5)
      stress(2,1)=vstress(6)

!      section if we want to plot stress on slip system
!      pi=3.141592653
!      theta=-70.0*pi/180.0      
!      t1=cos(theta+pi/2)*stress(1,1)+sin(theta+pi/2)*stress(1,2);
!      t2=cos(theta+pi/2)*stress(1,2)+sin(theta+pi/2)*stress(2,2);
!     // get shear stress //
!      stress(2,3)=t1*cos(theta)+t2*sin(theta);
      
       

      end


      subroutine dump_atom(x, b, logic)
      use mod_global
      use mod_dynamo
      double precision x(nxdm,*),b(ndf,*)
      integer logic, ntot, numpnts
      integer n1, n2, n3, npatoms

      double precision dev(6),xdef,ydef,zdef,hyd,sigeff,y,epseff
      integer numtri,i,j,n,NDFMAX,ii,jj,nout, npad, startpad
      double precision dwx1, dwx2
      double precision box_max(2), box_min(2)
      double precision pe, ke
      integer isurf

      double precision :: umag

      umag = 1.0d0

      pe = 0.0d0
      do j=1,2
         box_max(j)=-1E100
         box_min(j)=1E100
      enddo
      npatoms=0
      npad = 0
      ntot = 0
      do i = 1,numnp
         if( IsRelaxed(i)==1 .or. IsRelaxed(i)==2
     $        .or.  IsRelaxed(i)==-1) then
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
               if(  x(j,i) + umag*b(j,i) > box_max(j)) then
                  box_max(j)= x(j,i) + umag*b(j,i)
               endif
               if(  x(j,i) + umag*b(j,i) < box_min(j)) then
                  box_min(j)= x(j,i) + umag*b(j,i)
               endif
            enddo
         endif
      enddo
      startpad = npatoms - npad + 1
      print *, 'Start pad atoms = ', startpad
      print *, 'Total Potential Energy = ', pe
      print *, 'Average Potential Energy = ', dble(pe/npatoms)      

      do i = 1, numnp
         if (i.eq.1) then

            write(logic,'(''Number of particles = '',i5)') ntot
            write(logic,'(''A = 1.0 Angstrom '')')
            write(logic,'(''H0(1,1) = '',f10.4,'' A'')') 
     $           box_max(1)-box_min(1)
c     write(logic,'(''H0(1,2) = '',f10.4,'' A'')') box_min(2)
            write(logic,'(''H0(1,2) = 0 A'')')
            write(logic,'(''H0(1,3) = 0 A'')')
            write(logic,'(''H0(2,1) = 0 A'')') 
c     write(logic,'(''H0(2,1) = '',f10.4,'' A'')') box_min(1)
            write(logic,'(''H0(2,2) = '',f10.4,'' A'')') 
     $           box_max(2)-box_min(2) 
            write(logic,'(''H0(2,3) = 0 A'')')
            write(logic,'(''H0(3,1) = 0 A'')')
            write(logic,'(''H0(3,2) = 0 A'')')
c     write(logic,'(''H0(3,1) = '',f10.4,'' A'')') box_min(1) 
c     write(logic,'(''H0(3,2) = '',f10.4,'' A'')') box_min(2) 
            write(logic,'(''H0(3,3) = '',f10.3,'' A'')') z_length
            write(logic,'(''.NO_VELOCITY. '')') 
            write(logic,'(''entry_count = 6 '')')
            write(logic,'(''auxiliary[0] = pote [eV] '')') 
            write(logic,'(''auxiliary[1] = coord [num] '')') 
            write(logic,'(''auxiliary[2] = Surface [num] '')') 
            write(logic,'(''1.000000 '')')
            write(logic,'(''Al '')') 
         end if
         
         xdef =(x(1,i)+umag*b(1,i)+box_min(1))/(box_max(1)-box_min(1))
         ydef =(x(2,i)+umag*b(2,i)+box_max(1))
     $        /(box_max(2)-box_min(2))
         zdef = (x(3,i) + umag*b(3,i))/z_length
C     C           zdef = b(3,i)
         if( IsRelaxed(i)==1 .or. IsRelaxed(i)==2 ) then
            if (numNeighbors(i) .eq.  26) then 
               isurf = 0
            end if
            if ((numNeighbors(i) >  22 .and.numNeighbors(i) <= 26)
     $           .or. numNeighbors(i) >  26) then 
               isurf = 1
            end if

            if (numNeighbors(i) <= 22) then 
               isurf = 2
            end if

            write(logic,'(4f16.11,1x,2I4)')xdef,ydef,zdef, 1.d0
     $           *energy(i), NumNeighbors(i), isurf
         endif
      end do
      do i = 1, numnp
         if (i .eq. 1) then 
            write(logic,'(''Ni '')') 
         endif

         xdef =(x(1,i)+umag*b(1,i)+box_min(1))
     $        /(box_max(1)-box_min(1))
         ydef =(x(2,i)+umag*b(2,i)+box_max(1))
     $        /(box_max(2)-box_min(2))
         zdef = (x(3,i) + umag*b(3,i))/z_length
C     C           zdef = b(3,i)
         if( IsRelaxed(i)==-1) then
            write(logic,'(4f16.11,1X,2I4)')xdef,ydef,zdef, energy(i)
     $           ,NumNeighbors(i),0
         endif
      end do


      end subroutine dump_atom

c***********************************************************************
      subroutine dumpit(x,ix,b,db,id,f,dr,scale,logic,key,index,umag)
      use mod_global
      use mod_dynamo
      implicit none
c
      double precision x(nxdm,*),b(ndf,*),db(ndf,*),f(ndf,*),dr(ndf,*)
     $     ,scale,umag
      integer id(ndf,*),ix(nen1,*),index,logic,numpnts,ntot
      character*4 key
c
      integer numtri,i,j,n,NDFMAX,ii,jj,nout, npad, startpad
      character out*100
      character TAtom*2
      double precision Mass
      parameter(NDFMAX=18)
      double precision dev(6),xdef,ydef,zdef,hyd,sigeff,y,epseff

      integer n1,n2,n3,npatoms
      double precision det, amat(3,2), epsloc(3,3), tarray(3,3)
      double precision, allocatable ::  nstrain(:,:,:)
      double precision, allocatable ::  avgnum(:)
      double precision dwx1, dwx2
      double precision box_max(2), box_min(2)
      double precision pe, ke
      integer isurf
      pe = 0.0
c
      numtri = numel
      if(key.ne.'atom'.and.key.ne.'virH') then
         if(key .ne. 'viri' .and. key .ne. 'stra' .and. key .ne. 'stre')
     $        then 
c$$$            write(logic,'('' TITLE = " '',a4,'' "'')') key
         end if
      endif


!     if we are doing strains or stresses then calculate them at each node
      if (key.eq.'stra'.or.key.eq.'stre'.or.key.eq.'viri') then

         allocate(nstrain(3,3,numnp))
         allocate(avgnum(numnp))
         do  i = 1,numnp
           avgnum(i)=0.0
           do j=1,3
             do ii=1,3
               nstrain(ii,j,i)=0.0
             enddo
           enddo  
         end do

         do i=1,numel
            n1=ix(1,i)
            n2=ix(2,i)
            n3=ix(3,i)
            det =(x(1,n1)-x(1,n3))*(x(2,n2)-x(2,n3)) -
     &           (x(1,n2)-x(1,n3))*(x(2,n1)-x(2,n3))
            amat(1,1)=(x(2,n2)-x(2,n3))/det
            amat(2,1)=(x(2,n3)-x(2,n1))/det
            amat(3,1)=(x(2,n1)-x(2,n2))/det
            amat(1,2)=(x(1,n3)-x(1,n2))/det
            amat(2,2)=(x(1,n1)-x(1,n3))/det
            amat(3,2)=(x(1,n2)-x(1,n1))/det 
c           ****compute the green strain time 2*****
            call GetElementStrain(b(1,n1),b(1,n2),b(1,n3),amat(1:3,1:2)
     $        ,epsloc(1:3,1:3))
            epsloc(1:3,1:3)=epsloc(1:3,1:3)*0.5

            if(key.eq.'stre'.or. key.eq.'viri') then
c           ***** compute stresses in units of eV/A^3 
             call GetStresses(epsloc,tarray)
            else
             tarray(1:3,1:3)=epsloc(1:3,1:3)
c***         convert to engineering strains
             tarray(2,3)=tarray(2,3)+tarray(3,2)
             tarray(1,3)=tarray(1,3)+tarray(3,1)
             tarray(1,2)=tarray(1,2)+tarray(2,1) 
            endif

            do j=1,3
              do jj=1,3
                nstrain(j,jj,n1)=nstrain(j,jj,n1)+tarray(j,jj)
                nstrain(j,jj,n2)=nstrain(j,jj,n2)+tarray(j,jj)
                nstrain(j,jj,n3)=nstrain(j,jj,n3)+tarray(j,jj)
              enddo
            enddo

            avgnum(n1)=avgnum(n1)+1.0
            avgnum(n2)=avgnum(n2)+1.0
            avgnum(n3)=avgnum(n3)+1.0

        enddo

        do i=1,numnp
         if (avgnum(i).gt.0) then 
          do j=1,3
           do jj=1,3
             nstrain(j,jj,i)=nstrain(j,jj,i)/avgnum(i)
            enddo 
          enddo
         endif
        enddo            

      end if


c     find max and min of coordintes so that a box can be made.
      if (key.eq.'atom' .or. key .eq. 'ovit') then
         pe = 0.0d0
        do j=1,2
          box_max(j)=-1E100
          box_min(j)=1E100
        enddo
        npatoms=0
        npad = 0
        ntot = 0
        do i = 1,numnp
          if( IsRelaxed(i)==1 .or. IsRelaxed(i)==2 .or.
     $          IsRelaxed(i)==-1) then
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
            if(  x(j,i) + umag*b(j,i) > box_max(j)) then
                box_max(j)= x(j,i) + umag*b(j,i)
            endif
            if(  x(j,i) + umag*b(j,i) < box_min(j)) then
                box_min(j)= x(j,i) + umag*b(j,i)
            endif
           enddo
          endif
        enddo
        startpad = npatoms - npad + 1
        print *, 'Start pad atoms = ', startpad
        print *, 'Total Potential Energy = ', pe
        print *, 'Average Potential Energy = ', dble(pe/npatoms)
       endif

c     begin loop over nodes for all cases 
c      numpnts=max(numnp,numnpp1)
      numpnts=max(numnp,numnpp1)
      do 10 i = 1,numpnts
         if (key.eq.'noda') then
            dev(1) = scale*b(1,i)
            dev(2) = scale*b(2,i)
            if (i.eq.1) then
               write(logic,'('' VARIABLES = X Y UX UY VX VY AX AY'')')
               write(logic,'('' ZONE T = "ZONE ONE", I = '',i5,'',J = ''
     $              ,i5,'', F = FEPOINT'')') numpnts,numtri
            end if
            xdef = x(1,i) + umag*b(1,i)
            ydef = x(2,i) + umag*b(2,i)
            write(logic,'(8e14.6)') xdef,ydef,(dev(j),j=1,6)
         end if
         if (key.eq.'disp') then
            dev(1) = scale*b(1,i)
            dev(2) = scale*b(2,i)
            if (ndf.ge.3) dev(3) = scale*b(3,i)
            if (ndf.eq.4) dev(4) = scale*b(4,i)
            if (i.eq.1) then
               if(ndf.gt.NDFMAX) then
                  write(*,*)
     $                 '***WARNING: tecplot file can only show first'
     $                 ,NDFMAX
                  write(*,*) '            degrees of freedom'
               endif
               if (ndf.eq.2) then
                  write(logic,'('' VARIABLES = X Y UX UY '')')
               elseif (ndf.eq.3) then
                  write(logic,'('' VARIABLES = X Y Z UX UY UZ'')')
               elseif (ndf.eq.4) then
                  write(logic,'('' VARIABLES = X Y Z UX UY UZ EZ'')')
               else
                  out='VARIABLES = X Y Z UX UY UZ'
                  nout=26
                  jj=1
                  do ii=4,min(ndf,NDFMAX)
                     if(jj.eq.1) then
                        out(nout+1:nout+3)=' SX'
                        jj=2
                     else if (jj.eq.2) then
                        out(nout+1:nout+3)=' SY'
                        jj=3
                     else if (jj.eq.3) then
                        out(nout+1:nout+3)=' SZ'
                        jj=1
                     endif
                     nout=nout+3
                  enddo
                  write(logic,*) out
               endif
               write(logic,'('' ZONE T = "ZONE ONE", I = '',i5,'',J = ''
     $              ,i5,'', F = FEPOINT'')') numpnts,numtri
            end if
            xdef = x(1,i) + umag*b(1,i)
            ydef = x(2,i) + umag*b(2,i)
CC            print*, 'spatial dimension', ndf
            if (ndf.ge.3) zdef = umag*b(3,i) + x(3,i)
            if (ndf.eq.2) then
               write(logic,'(4e14.6)') xdef,ydef,dev(1),dev(2)
            elseif (ndf.eq.3) then
               write(logic,'(6e14.6)') xdef,ydef,zdef,dev(1),dev(2)
     $              ,dev(3)
            elseif (ndf.eq.4) then
               write(logic,'(7e14.6)') xdef,ydef,zdef,dev(1),dev(2)
     $              ,dev(3),dev(4)
            else
               write(logic,'(21e14.6)') xdef,ydef,zdef,dev(1),dev(2)
     $              ,dev(3),(b(ii,i),ii=4,min(ndf,NDFMAX))
            endif
         end if

         if (key.eq.'bcon') then
c     dev(1) = scale*f(1,i)
c     dev(2) = scale*f(2,i)
c     if (ndf.ge.3) dev(3) = scale*f(3,i)
c            dev(1) = scale*dr(1,i)
c            dev(2) = scale*dr(2,i)
c            if (ndf.ge.3) dev(3) = scale*dr(3,i)
            dev(1) = dble(id(1,i))
            dev(2) = dble(id(2,i))
            if (ndf.ge.3) dev(3) = dble(id(3,i))
            if (i.eq.1) then
               if (ndf.eq.2) then
                  write(logic,'('' VARIABLES = X Y FX FY'')')
               else
                  write(logic,'('' VARIABLES = X Y FX FY FZ'')')
               endif
               write(logic,'('' ZONE T = "ZONE ONE", I = '',i5,'', J =''
     $              ,i5,'', F = FEPOINT'')') numpnts,numtri
            end if
            xdef = x(1,i) + umag*b(1,i)
            ydef = x(2,i) + umag*b(2,i)
            if (ndf.eq.2) then
               write(logic,'(5e14.6)') xdef,ydef,dev(1),dev(2)
            else
               write(logic,'(6e14.6)') xdef,ydef,dev(1),dev(2),dev(3)
            endif
         end if


         if (key.eq.'ener') then
            dev(1) = 1.d0*energy(i)
            if (i.eq.1) then
               write(logic,'('' VARIABLES = X Y Z E_atom'')')
               write(logic,'('' ZONE T = "ZONE ONE", I = '',i5,'',J = ''
     $              ,i5,'', F = FEPOINT'')') numpnts,numtri
            end if
            xdef = x(1,i) + umag*b(1,i)
            ydef = x(2,i) + umag*b(2,i)
c            zdef =  b(3,i)
            zdef = x(3,i) + umag*b(3,i)
            write(logic,'(6e14.6)') xdef,ydef,zdef,dev(1)
         end if


         if (key.eq.'atom')then

          if (i.eq.1) then

              write(logic,'(''Number of particles = '',i5)') ntot
              write(logic,'(''A = 1.0 Angstrom '')')
              write(logic,'(''H0(1,1) = '',f10.4,'' A'')') 
     $           box_max(1)-box_min(1)
c               write(logic,'(''H0(1,2) = '',f10.4,'' A'')') box_min(2)
              write(logic,'(''H0(1,2) = 0 A'')')
              write(logic,'(''H0(1,3) = 0 A'')')
              write(logic,'(''H0(2,1) = 0 A'')') 
c              write(logic,'(''H0(2,1) = '',f10.4,'' A'')') box_min(1)
              write(logic,'(''H0(2,2) = '',f10.4,'' A'')') 
     $           box_max(2)-box_min(2) 
              write(logic,'(''H0(2,3) = 0 A'')')
              write(logic,'(''H0(3,1) = 0 A'')')
              write(logic,'(''H0(3,2) = 0 A'')')
c              write(logic,'(''H0(3,1) = '',f10.4,'' A'')') box_min(1) 
c              write(logic,'(''H0(3,2) = '',f10.4,'' A'')') box_min(2) 
              write(logic,'(''H0(3,3) = '',f10.3,'' A'')') z_length
              write(logic,'(''.NO_VELOCITY. '')') 
              write(logic,'(''entry_count = 6 '')')
              write(logic,'(''auxiliary[0] = pote [eV] '')') 
              write(logic,'(''auxiliary[1] = coord [num] '')') 
              write(logic,'(''auxiliary[2] = Surface [num] '')') 
              write(logic,'(''1.000000 '')')
              write(logic,'(''Al '')') 
           end if
           
           xdef =(x(1,i)+umag*b(1,i)+box_min(1))/(box_max(1)-box_min(1))
           ydef =(x(2,i)+umag*b(2,i)+box_max(1))
     $          /(box_max(2)-box_min(2))
           zdef = (x(3,i) + umag*b(3,i))/z_length
CC           zdef = b(3,i)
           if( IsRelaxed(i)==1 .or. IsRelaxed(i)==2 ) then
              if (numNeighbors(i) .eq.  26) then 
                 isurf = 0
              end if
              if ((numNeighbors(i) >  22 .and.numNeighbors(i) <= 26)
     $             .or. numNeighbors(i) >  26) then 
                 isurf = 1
              end if

              if (numNeighbors(i) <= 22) then 
                 isurf = 2
              end if

              write(logic,'(4f16.11,1x,2I4)')xdef,ydef,zdef, 1.d0
     $             *energy(i), NumNeighbors(i), isurf
           endif

        endif

         if (key.eq.'ovit')then

          if (i.eq.1) then

              write(logic,'(''Number of particles = '',i5)') ntot
              write(logic,'(''A = 1.0 Angstrom '')')
              write(logic,'(''H0(1,1) = '',f10.4,'' A'')') 
     $           box_max(1)-box_min(1)
              write(logic,'(''H0(1,2) = 0 A'')')
              write(logic,'(''H0(1,3) = 0 A'')')
              write(logic,'(''H0(2,1) = 0 A'')') 
              write(logic,'(''H0(2,2) = '',f10.4,'' A'')') 
     $           box_max(2)-box_min(2) 
              write(logic,'(''H0(2,3) = 0 A'')')
              write(logic,'(''H0(3,1) = 0 A'')')
              write(logic,'(''H0(3,2) = 0 A'')')
              write(logic,'(''H0(3,3) = '',f10.3,'' A'')') z_length
           end if
           
           xdef =(x(1,i)+umag*b(1,i)+box_min(1))/(box_max(1)-box_min(1))
           ydef =(x(2,i)+umag*b(2,i)+box_max(1))
     $          /(box_max(2)-box_min(2))
           zdef = (x(3,i) + umag*b(3,i))/z_length
           if (IsRelaxed(i) == 1 .or. IsRelaxed(i) == 2) then 
              TAtom = "Al"
              Mass = 13.0
           end if
           if (IsRelaxed(i) == -1) then 
              TAtom = "Ni"
              Mass = 14.0
           end if

CC           zdef = b(3,i)
           if( IsRelaxed(i)==1 .or. IsRelaxed(i)==2 .or. IsRelaxed(i) ==
     $          -1) then
              write(logic,'(f4.0,1X,A2, 1X, 6f16.11)')Mass, TAtom,
     $          xdef,ydef,zdef,0.0,0.0,0.0
           endif

        endif



CC--Jun Song: Modifications: Output pdb file
         if (key.eq.'pdbf')then
CC--current position=origin+disp
         xdef = x(1,i)+umag*b(1,i)
         ydef = x(2,i)+umag*b(2,i)
         zdef = x(3,i)+umag*b(3,i)
CC--make sure points in simulation cell
         if((zdef .lt. 0.0d0).or.(zdef .gt. z_length)) then
           zdef = zdef-z_length*floor(zdef/z_length)
         endif
         
         if( IsRelaxed(i)==1 .or. IsRelaxed(i)==2 ) then
         if(atomSpecie(i)==1) then
         write(logic,'(A,i6,A,i12,f12.3,2f8.3)') "HETATM",i
     $   	,"Ni",1,xdef,ydef,zdef
         else
         write(logic,'(A,i6,A,i12,f12.3,2f8.3)') "HETATM",i
     $   	,"Al",2,xdef,ydef,zdef
         end if  
         endif
CC--JS output pad atoms also
         if(IsRelaxed(i)==-1) then
           write(logic,'(A,i6,A,i12,f12.3,2f8.3)') "HETATM",i
     $          ,"Fe",3,xdef,ydef,zdef
         endif
CC--output pad ends
         endif
CC--Jun Song Comment ends

         if (key.eq.'stra'.or.key.eq.'stre') then
            if (i.eq.1) then
              write(logic,
     $           '('' VARIABLES = X Y Z E11 E22 E33 E23 E13 E12 '')')
               write(logic,'('' ZONE T = "ZONE ONE", I = '',i5,'',J = ''
     $              ,i5,'', F = FEPOINT'')') numpnts,numtri
            end if
            xdef = x(1,i) + umag*b(1,i)
            ydef = x(2,i) + umag*b(2,i)
            zdef = x(3,i) + umag*b(3,i)
            write(logic,'(9e14.6)') xdef,ydef,zdef,nstrain(1,1,i)
     $          ,nstrain(2,2,i),nstrain(3,3,i),nstrain(2,3,i)
     $          ,nstrain(1,3,i),nstrain(1,2,i)
c            if(i.eq.numpnts) then
c               deallocate(nstrain,avgnum)
c            endif
         end if

         if (key.eq.'inte') then
            if (i.eq.1) then
               write(logic,'('' VARIABLES = X Y SED'')')
               write(logic,'('' ZONE T = "ZONE ONE", I = '',i5,'', J= ''
     $              ,i5,'', F = FEPOINT'')') numpnts,numtri
            end if
            xdef = x(1,i) + umag*b(1,i)
            ydef = x(2,i) + umag*b(2,i)
            write(logic,'(3e15.6)') xdef,ydef
         end if

         

         if (key.eq.'viri') then
            if (i.eq.1) then
              write(logic,
     $           '('' VARIABLES = X Y Z V11 V22 V33 V23 V13 V12 '')')
               write(logic,'('' ZONE T = "ZONE ONE", I = '',i5,'',J = ''
     $              ,i5,'', F = FEPOINT'')') numpnts,numtri
            end if
            xdef = x(1,i) + umag*b(1,i)
            ydef = x(2,i) + umag*b(2,i)
            zdef = x(3,i) + umag*b(3,i)
            if( IsRelaxed(i).eq.1 .or. IsRelaxed(i).eq.2) then
               write(logic,'(9e14.6)') xdef,ydef,zdef,avevirst(1,1,i)
     $          ,avevirst(2,2,i),avevirst(3,3,i),avevirst(2,3,i)
     $          ,avevirst(1,3,i),avevirst(1,2,i)
            else
               write(logic,'(6e14.6)') xdef,ydef,zdef,nstrain(1,1,i)
     $          ,nstrain(2,2,i),nstrain(3,3,i),nstrain(2,3,i)
     $          ,nstrain(1,3,i),nstrain(1,2,i)
     	    end if              
         end if

cc--Now output H position only!         
         if (key.eq.'virH') then
       	   if( IsRelaxed(i).eq.1 .and. atomSpecie(i).eq.2) then
             xdef = x(1,i) + umag*b(1,i)
             ydef = x(2,i) + umag*b(2,i)
             zdef = x(3,i) + umag*b(3,i)

             write(logic,'(6e14.6)') xdef,ydef,zdef
c     &            avevirst(1,1,i)
c     $          ,avevirst(2,2,i),avevirst(3,3,i),avevirst(2,3,i)
c     $          ,avevirst(1,3,i),avevirst(1,2,i)      
  	   end if
         end if 


 10   continue
      do i = 1, numnp
         if (key .eq. 'atom') then 
            if (i .eq. 1) then 
               write(logic,'(''Ni '')') 
            endif

            xdef =(x(1,i)+umag*b(1,i)+box_min(1))
     $           /(box_max(1)-box_min(1))
            ydef =(x(2,i)+umag*b(2,i)+box_max(1))
     $           /(box_max(2)-box_min(2))
            zdef = (x(3,i) + umag*b(3,i))/z_length
C     C           zdef = b(3,i)
            if( IsRelaxed(i)==-1) then
               write(logic,'(4f16.11,1X,2I4)')xdef,ydef,zdef, energy(i)
     $              ,NumNeighbors(i),0
            endif
         endif
      end do

      if((key.ne.'atom').and.(key.ne.'pdbf').and.(key.ne.'virH')) then
       do 20 n = 1,numel
         write(logic,'(4i6)') ix(1,n),ix(2,n),ix(3,n),ix(3,n)
 20    continue
      endif

      return
      end


************************************************************************
      subroutine dumpit_vtk(x,ix,b,db,id,f,dr,scale,logic,key,index
     $     ,umag)
      use mod_global
      implicit none
c
      double precision x(nxdm,*),b(ndf,*),db(ndf,*),f(ndf,*),dr(ndf,*)
     $     ,scale,umag
      integer id(ndf,*),ix(nen1,*),index,logic,numpnts,ntot
      character*4 key
c
      integer numtri,i,j,n,NDFMAX,ii,jj,nout, npad, startpad
      character out*100
      parameter(NDFMAX=18)
      double precision dev(6),xdef,ydef,zdef,hyd,sigeff,y,epseff

      integer n1,n2,n3,npatoms
      double precision det, amat(3,2), epsloc(3,3), tarray(3,3)
      double precision, allocatable ::  nstrain(:,:,:)
      double precision, allocatable ::  avgnum(:)
      double precision,allocatable :: virist(:,:), nstr1(:,:)
      double precision dwx1, dwx2
      double precision box_max(2), box_min(2)
      double precision :: rx1(3), ud1(3), b1(3), b2(3), b3(3), ud2(3),
     $     ud3(3), rx2(3), rx3(3)
      double precision :: sout1(3), sout2(3), sout3(3)
c     New Variable for vtk 
      integer ndf22
      double precision :: ev_convert, fact
      ev_convert = 1.602176462
c     Converts from ev/A^3 to MPa
      fact = 1.d0/(ev_convert)/1.d-5

c
      numtri = numel
c$$$      if(key.ne.'atom'.and.key.ne.'virH') then
c$$$c$$$        write(logic,'('' TITLE = " '',a4,'' "'')') key
c$$$      endif


!     if we are doing strains or stresses then calculate them at each node
      if (key.eq.'stra'.or.key.eq.'stre'.or.key.eq.'viri') then

         allocate(nstrain(3,3,numnp))
         allocate(nstr1(3,3))
         allocate(virist(3,3))
         allocate(avgnum(numnp))
         do  i = 1,numnp
            avgnum(i)=0.0
            do j=1,3
               do ii=1,3
                  nstrain(ii,j,i)=0.0
               end do
            end do  
         end do

         do i=1,numel
            n1=ix(1,i)
            n2=ix(2,i)
            n3=ix(3,i)
            b1(1:3)=b(1:3, n1)
            b2(1:3)=b(1:3, n2)
            b3(1:3)=b(1:3, n3)
c     ****************************************************
c     Adding back displacement contribution to satisfy superposition
            rx1(1:3)=x(1:3,n1)
            rx2(1:3)=x(1:3,n2)
            rx3(1:3)=x(1:3,n3)

            call disl_displ(rx1,ud1)
            call disl_displ(rx2,ud2)
            call disl_displ(rx3,ud3)

c     Can also call disl_stress to calculate stress directly at the
c     nodal points to make sure that everything is ok
c     Currently the plots show some anamalies due to interpolation and
c     slip planes passing through elements etc. 
c     Interpolation from FE displacments is not quite right
            b1 = b1 - ud1
            b2 = b2 - ud2 
            b3 = b3 - ud3
c     ****************************************************            
            det =(x(1,n1)-x(1,n3))*(x(2,n2)-x(2,n3)) -
     &           (x(1,n2)-x(1,n3))*(x(2,n1)-x(2,n3))
            amat(1,1)=(x(2,n2)-x(2,n3))/det
            amat(2,1)=(x(2,n3)-x(2,n1))/det
            amat(3,1)=(x(2,n1)-x(2,n2))/det
            amat(1,2)=(x(1,n3)-x(1,n2))/det
            amat(2,2)=(x(1,n1)-x(1,n3))/det
            amat(3,2)=(x(1,n2)-x(1,n1))/det 
c           ****compute the green strain time 2*****
            call GetElementStrain(b1,b2,b3,amat(1:3,1:2)
     $        ,epsloc(1:3,1:3))

c$$$            call GetElementStrain(b1,b2,b3,amat(1:3,1:2),epsloc(1:3
c$$$     $           ,1:3))

            epsloc(1:3,1:3)=epsloc(1:3,1:3)*0.5

            if(key.eq.'stre'.or. key.eq.'viri') then
c           ***** compute stresses in units of eV/A^3 
               call GetStresses(epsloc,tarray)
            else
               tarray(1:3,1:3)=epsloc(1:3,1:3)
c***  convert to engineering strains
               tarray(2,3)=tarray(2,3)+tarray(3,2)
               tarray(1,3)=tarray(1,3)+tarray(3,1)
               tarray(1,2)=tarray(1,2)+tarray(2,1) 
            end if


            do j=1,3
               do jj=1,3
                  nstrain(j,jj,n1)=nstrain(j,jj,n1)+tarray(j,jj)
                  nstrain(j,jj,n2)=nstrain(j,jj,n2)+tarray(j,jj)
                  nstrain(j,jj,n3)=nstrain(j,jj,n3)+tarray(j,jj)
               end do
            end do
            if (key .eq. 'stre' .or. key .eq. 'viri') then 
c           ***** Compute superposition stress due to dislocations
               sout1 = 0.0d0
               sout2 = 0.0d0
               sout3 = 0.0d0
               call disl_stress(rx1, sout1)
               call disl_stress(rx2, sout2)
               call disl_stress(rx3, sout3)

               nstrain(1,1,n1) = nstrain(1,1,n1) + sout1(1)
               nstrain(1,2,n1) = nstrain(1,2,n1) + sout1(3)
               nstrain(2,1,n1) = nstrain(1,2,n1) + sout1(3)
               nstrain(2,2,n1) = nstrain(2,2,n1) + sout1(2)

               nstrain(1,1,n2) = nstrain(1,1,n2) + sout2(1)
               nstrain(1,2,n2) = nstrain(1,2,n2) + sout2(3)
               nstrain(2,1,n2) = nstrain(1,2,n2) + sout2(3)
               nstrain(2,2,n2) = nstrain(2,2,n2) + sout2(2)

               nstrain(1,1,n3) = nstrain(1,1,n3) + sout3(1)
               nstrain(1,2,n3) = nstrain(1,2,n3) + sout3(3)
               nstrain(2,1,n3) = nstrain(1,2,n3) + sout3(3)
               nstrain(2,2,n3) = nstrain(2,2,n3) + sout3(2)
            end if

c$$$            call disl_stress(rx1, ud1)
c$$$            call disl_stress(rx2, ud2)
c$$$            call disl_stress(rx3, ud3)
c$$$
c$$$            do j = 1,3
c$$$               do jj = 1, 3
c$$$                  if (j==1 .and. jj==1) then 
c$$$                     nstrain(j,j,n1) = nstrain(j,j,n1) + ud1(1)
c$$$                     nstrain(j,j,n2) = nstrain(j,j,n2) + ud2(1)
c$$$                     nstrain(j,j,n3) = nstrain(j,j,n3) + ud3(1)
c$$$                  end if
c$$$                  if (j==1 .and. jj==2) then 
c$$$                     nstrain(j,j,n1) = nstrain(j,j,n1) + ud1(1)
c$$$                     nstrain(j,j,n2) = nstrain(j,j,n2) + ud2(1)
c$$$                     nstrain(j,j,n3) = nstrain(j,j,n3) + ud3(1)
c$$$                  end if
c$$$
c$$$               end do
c$$$            end do
            avgnum(n1)=avgnum(n1)+1.0
            avgnum(n2)=avgnum(n2)+1.0
            avgnum(n3)=avgnum(n3)+1.0
         end do

        do i=1,numnp
         if (avgnum(i).gt.0) then 
          do j=1,3
           do jj=1,3
             nstrain(j,jj,i)=nstrain(j,jj,i)/avgnum(i)
            enddo 
          enddo
         endif
        enddo            

      end if


c     find max and min of coordintes so that a box can be made.
      if (key.eq.'atom') then
         do j=1,2
            box_max(j)=-1E100
            box_min(j)=1E100
         end do
        npatoms=0
        npad = 0
        ntot = 0
        do i = 1,numnp
          if( IsRelaxed(i)==1 .or. IsRelaxed(i)==2 .or.
     $          IsRelaxed(i)==-1) then
             if (IsRelaxed(i) ==-1) then 
                if (npad == 0) then 
                   startpad = i
                end if
                npad = npad + 1
             else
                npatoms=npatoms+1
             end if
             ntot = ntot + 1
             do j=1,2
                if(  x(j,i) + umag*b(j,i) > box_max(j)) then
                   box_max(j)= x(j,i) + umag*b(j,i)
                end if
                if(  x(j,i) + umag*b(j,i) < box_min(j)) then
                   box_min(j)= x(j,i) + umag*b(j,i)
                end if
             end do
          end if
       end do
       startpad = npatoms - npad + 1
       print *, 'Start pad atoms = ', startpad
      end if

c     begin loop over nodes for all cases 
c      numpnts=max(numnp,numnpp1)
      numpnts=max(numnp,numnpp1)
      if (key .eq. 'stra' .or. key .eq. 'stre' .or. key .eq. 'viri')
     $     then
         write(logic,fmt='(A)') '# vtk DataFile Version 2.0'
         if (key .eq. 'stra') then 
            write(logic,fmt='(A)') 'Strains from CADD'
         end if
         if (key .eq. 'stre') then 
            write(logic,fmt='(A)') 'Stresses from CADD'
         end if
         if (key .eq. 'viri') then 
            write(logic,fmt='(A)') 'Virial Stresses from CADD'
         end if
         write(logic,fmt='(A)') 'ASCII'
         write(logic,fmt='(A)') 'DATASET UNSTRUCTURED_GRID'
         write(logic,fmt='(A6,1x,I7,1x,A5)') 'POINTS',numpnts,'float'
         do i = 1, numpnts   
            
            xdef = x(1,i) + umag*b(1,i)
            ydef = x(2,i) + umag*b(2,i)
            zdef = x(3,i) + umag*b(3,i)
            write(logic,'(3(1x,e14.6))') xdef,ydef, zdef
         end do
         write(logic,*)
         write(logic,'(A5,1X,I7,1X,I7)') 'CELLS', numel, 4*numtri
         
         do n = 1,numel
            write(logic,'(4i6)') 3, ix(1,n)-1,ix(2,n)-1,ix(3,n)-1
         end do
         write(logic,*)

         write(logic,'(A10,1X,I7)') 'CELL_TYPES', numel

         do n=1,numel
            write(logic,fmt='(5(1x,I7))') 5
         end do

         write(logic, *) 'POINT_DATA', numpnts

         if (key .eq. 'stra') then 
            write(logic, *) 'SCALARS strain float 3'
         else
            write(logic, *) 'TENSORS stress float'
         end if

c$$$         write(logic, *) 'LOOKUP_TABLE default'
         
         do i = 1, numpnts
            if (key .eq. 'stra') then 
               write(logic,'(6(1x,e14.6))') nstrain(1,1,i), nstrain(2
     $              ,2,i), nstrain(1,2,i)
c$$$               , nstrain(2,3,i), nstrain(1
c$$$     $              ,3,i), nstrain(1,2,i)
            end if
            virist(:,:) = avevirst(:,:,i)*fact
            nstr1(:,:) = nstrain(:,:,i)*fact
            if (key .eq. 'viri') then 
               if( IsRelaxed(i).eq. 1 .or. IsRelaxed(i).eq.2) then
                  write(logic,'(3e14.6)') virist(1,1),virist(1
     $                 ,2),virist(1,3)
                  write(logic,'(3e14.6)') virist(1,2),virist(2
     $                 ,2),virist(2,3)
                  write(logic,'(3e14.6)') virist(1,2),virist(2
     $                 ,2),virist(2,3)
                  write(logic,*)
               else
                  write(logic,'(3e14.6)') nstr1(1,1),nstr1(1,2
     $                 ),nstr1(1,3)
                  write(logic,'(3e14.6)') nstr1(1,2),nstr1(2,2
     $                 ),nstr1(2,3)
                  write(logic,'(3e14.6)') nstr1(1,3),nstr1(2,3
     $                 ),nstr1(3,3)
                  write(logic,*)
               end if
            end if
         end do
         write (logic,*) 'VECTORS F float'
         do i = 1, numpnts
            write(logic,'(3e14.6)') dble(id(1,i)), dble(id(2,i)), 0.0
         end do
         write(logic, *) 'CELL_DATA', numel
         write(logic, *) 'SCALARS DB integer'
         write(logic, *) 'LOOKUP_TABLE default'

         do n=1,numel
            write(logic, '(i5)') ix(nen1,n)
         end do

         
      end if
      return
      end
c***************************************************************
      subroutine dump_mesh(x,b,ix, logic, iq)
      use mod_global
      implicit none
c
      double precision x(nxdm,*),b(ndf,*),scale,umag
      integer ix(nen1,*),index,logic,numpnts,ntot
      integer, optional :: iq
      character*4 key
c
      integer numtri,i,j,n,NDFMAX,ii,jj,nout, npad, startpad
      character out*100
      parameter(NDFMAX=18)
      double precision dev(6),xdef,ydef,zdef,hyd,sigeff,y,epseff

      integer n1,n2,n3,npatoms
      double precision det, amat(3,2), epsloc(3,3), tarray(3,3)
      double precision, allocatable ::  nstrain(:,:,:)
      double precision, allocatable ::  avgnum(:)
      double precision dwx1, dwx2
      double precision box_max(2), box_min(2)
      double precision :: rx1(3), ud1(3), b1(3), b2(3), b3(3), ud2(3),
     $     ud3(3), rx2(3), rx3(3)

      umag = 1.0d0
      numpnts = numnp
      numtri = numel

      write(logic,fmt='(A)') '# vtk DataFile Version 2.0'
      write(logic,fmt='(A)') 'Deformed mesh from CADD'
      write(logic,fmt='(A)') 'ASCII'
      write(logic,fmt='(A)') 'DATASET UNSTRUCTURED_GRID'
      write(logic,fmt='(A6,1x,I7,1x,A5)') 'POINTS',numpnts,'float'
      do i = 1, numpnts   
         if (present(iq)) then 
            xdef = x(1,i)-xtip_actual(1) + umag*b(1,i)
            ydef = x(2,i)-xtip_actual(2) + umag*b(2,i)
            zdef = x(3,i) + umag*b(3,i)
         else
            xdef = x(1,i) + umag*b(1,i)
            ydef = x(2,i) + umag*b(2,i)
            zdef = x(3,i) + umag*b(3,i)
         end if
         write(logic,'(3(1x,e14.6))') xdef,ydef, zdef
      end do
      write(logic,*)
      write(logic,'(A5,1X,I7,1X,I7)') 'CELLS', numel, 4*numtri
      
      do n = 1,numel
         write(logic,'(4i6)') 3, ix(1,n)-1,ix(2,n)-1,ix(3,n)-1
      end do
      write(logic,*)
      
      write(logic,'(A10,1X,I7)') 'CELL_TYPES', numel
      
      do n=1,numel
         write(logic,fmt='(5(1x,I7))') 5
      end do
      end subroutine dump_mesh

c***************************************************************
