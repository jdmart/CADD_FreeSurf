c**---------------------------------------------------------------
c**   pdel : write a load-delta point to the pdel file
c**
      subroutine pdel(f,dr,id,logic,ndf,x,time)
      implicit none
      double precision f,dr,x,time
      integer id,logic,ndf
c
      double precision force
      double precision totener,strener
      common/enerdat/ totener,strener
c
      call pdelcalc(f,dr,id,force,x)
      write(logic,'(4e15.5)') time,force,totener,strener
      call flush(logic)
      return
      end
c**---------------------------------------------------------------------
c** InitialiseEnergy
c**
c--
      subroutine InitialiseEnergy(e_flag,f_flag,id,dr,f)
      use mod_global
      implicit none

c--passed variables
      logical e_flag,f_flag
      integer id(ndf,*)
      double precision dr(ndf,*),f(ndf,*)
c--Local Variables
      integer i,j

c     Initialize out-of-balance force vector
      if (f_flag) then
         do i=1,numnp
            do j=1,ndf
               if(id(j,i).eq.0) then
                  dr(j,i)=time*f(j,i)
               else
                  dr(j,i)=0.d0
               endif
            enddo
         enddo
         if(numnpp1.gt.numnp) dr(1:ndf,numnpp1) =0.d0
      endif

      if(e_flag) then
         do i = 1, numnp
            energy(i) = 0.0
         end do
      endif
      return
      end

!**
      subroutine StatusCalc(x,ix,silent)
      use mod_file
      use mod_global
      implicit none
!
      integer ix(nen1,numel)
      double precision x(nxdm,numnp),s(3)
      character*80 filename
      integer i,j,k,icheck,icheck2, logic
      logical first,silent,intri,ontri

      if (allocated(IsRelaxed)) deallocate(IsRelaxed)
      if (.not. allocated(IsRelaxed)) allocate(IsRelaxed(numnp))
!--- only atoms in the atomistic region (abs(ix(nen1,i)).eq.1) are "relaxed"
!--- in CG.
!-- IsRelaxed=0: continuum
!-- IsRelaxed=1: atomistic
!-- IsRelaxed=2: interface
!
      IsRelaxed(1:num2Dnode)=-1
!
      do i=1,num2Dnode
         first=.true.
         do j=1,numel
            do k=1,3
               if(ix(k,j).eq.i) then
                  if(first) then
                     icheck=ix(nen1,j)
                     if(icheck.lt.0) icheck=1
                     first=.false.
                     if(ix(nen1,j).eq.0) then
                        IsRelaxed(i)=0
                     else
                        IsRelaxed(i)=1
                     endif
                     go to 110
                  else
                     icheck2=ix(nen1,j)
                     if(icheck2.lt.0) icheck2=1
                     if(icheck.eq.icheck2) then
                        if(ix(nen1,j).eq.0) then
                           IsRelaxed(i)=0
                        else
                           IsRelaxed(i)=1
                        endif
                        go to 110
                     else
                        IsRelaxed(i)=2
                        go to 120
                     endif
                  endif
               endif
            enddo
 110        continue
         enddo
 120     continue
      enddo
      if(numnpp1.gt.numnp) IsRelaxed(numnpp1)=1

      if(silent) return
      nqc=0
      nspring=0
      do i=1,numnp
         if(IsRelaxed(i).eq.1) then
            nqc=nqc+1
         else if (IsRelaxed(i).eq.2) then
            nspring=nspring+1
         endif
      enddo
      write(*,*) 'Number of nodes in the atomistic region:',nqc
      write(*,*) 'Number of nodes in the elastic region:',numnp-nqc
     $     -nspring
      write(*,*) 'Number of nodes on the interface:',nspring
      if(nen1.ne.4) stop 'ERROR: nen must be 3!'
      filename='out/check.plt'
      open(unit=123,file=filename,status='unknown')
      write(123,*) 'VARIABLES = X Y Z ISRELAXED'
      write(123,*) 'zone, t=continuum'
      do i=1,numnp
         if(IsRelaxed(i).eq.0) write(123,1000) x(1:3,i),0
      enddo
      write(123,*) 'zone, t=atomistic'
      do i=1,numnp
         if(IsRelaxed(i).eq.1) write(123,1000) x(1:3,i),1
      enddo
      write(123,*) 'zone, t=interface'
      do i=1,numnp
         if(IsRelaxed(i).eq.2) write(123,1000) x(1:3,i),2
      enddo
      write(123,*) 'zone, t=Pad'
      do i=1,numnp
         if(IsRelaxed(i).eq.-1) write(123,1000) x(1:3,i),-1
      enddo
      close(123)
 1000 format(3e15.6,i10)
      filename='out/check.vtk'
      call iofile(filename,'formatted  ',logic,.true.)
      write(logic,fmt='(A)') '# vtk DataFile Version 2.0'
      write(logic,fmt='(A)') 'Atoms/Nodes colored by region'
      write(logic,fmt='(A)') 'ASCII'
      write(logic,fmt='(A)') 'DATASET POLYDATA'
      write(logic,fmt='(A6,1x,I7,1x,A5)') 'POINTS',numnp,'float'
      do i=1,numnp
         write(logic, '(3e13.5)') x(1,i), x(2,i), 0.0
      end do
      
      write(logic, *) 'POINT_DATA', numnp

      write(logic, *) 'SCALARS atoms int  1'
      write(logic, *) 'LOOKUP_TABLE default'
      do i = 1, numnp
         write(logic,*) isRelaxed(i)
      end do


      return
      end subroutine StatusCalc
************************************************************************
      subroutine latticecheck(x)
      use mod_global
      implicit none
      double precision x(nxdm,numnp)
      integer i,igrain(3)
      double precision pm(3)

! make sure nodes are on atomic sites
      do i = 1, num2Dnode
         call NearestBSite(x(1,i),1,.true.,pm,igrain(1))
         if((abs(pm(1)-x(1,i)).gt.0.0001).or.(abs(pm(2)-x(2,i)).gt.0
     $        .0001)) then
            write(*,*) '***WARNING: node is not on atomic site!'
            write(*,*) i,x(1,i),x(2,i),x(3,i)
            write(*,*) igrain(1),pm(1),pm(2),pm(3)
            write(*,*) pm(1)-x(1,i),pm(2)-x(2,i),pm(3)-x(3,i)
         endif
      end do
      end subroutine latticecheck
************************************************************************
c**---------------------------------------------------------------
c**     PlotEsi : makes the esi.tec file
c**
c--
      subroutine PlotEsi(x,ix)
      use mod_global
      use mod_file
      implicit none

      double precision x(nxdm,*)
      integer ix(nen1,*)
c
      character*80 vars,filename
      integer logic,n,i1,i2,i3,i,iregion
      double precision x1,x2,x3,y1,y2,y3

c     Write file header
      filename='out/esi.tec'
      call iofile(filename,'formatted  ',logic,.true.)
      write(logic,'('' VARIABLES = X Y IREGION'')')
      write(logic,'('' ZONE T = "ZONE ONE", I = '',i5,'', J = '',i5,
     1 '', F = FEPOINT'')') numel*3,numel

      do n = 1,numel
         iregion=ix(nen1,n)
         i1=ix(1,n)
         i2=ix(2,n)
         i3=ix(3,n)
          x1=x(1,i1)
          x2=x(1,i2)
          x3=x(1,i3)
          y1=x(2,i1)
          y2=x(2,i2)
          y3=x(2,i3)
          write(logic,'(2e13.5,i4)') x1,y1,iregion
          write(logic,'(2e13.5,i4)') x2,y2,iregion
          write(logic,'(2e13.5,i4)') x3,y3,iregion
      enddo
      do n=1,numel
         i=(n-1)*3
         write(logic,'(4i6)') i+1,i+2,i+3,i+3
      enddo
      close(logic)
      return
      end
************************************************************************




************************************************************************
c**---------------------------------------------------------------
c**     PlotEsi : makes the esi.tec file
c**
c--
      subroutine PlotEsi_vtk(x,ix)
      use mod_global
      use mod_file
      implicit none

      double precision x(nxdm,*)
      integer ix(nen1,*)
c
      character*80 vars,filename
      integer logic,n,i1,i2,i3,i,iregion
      double precision x1,x2,x3,y1,y2,y3
      integer scal(numel)

c     Write file header
      filename='out/esi.vtk'
      call iofile(filename,'formatted  ',logic,.true.)
c$$$      write(logic,'('' VARIABLES = X Y IREGION'')')
c$$$      write(logic,'('' ZONE T = "ZONE ONE", I = '',i5,'', J = '',i5,
c$$$     1 '', F = FEPOINT'')') numel*3,numel

      write(logic,fmt='(A)') '# vtk DataFile Version 2.0'
      write(logic,fmt='(A)') 'Detection Band Check from CADD'
      write(logic,fmt='(A)') 'ASCII'
      write(logic,fmt='(A)') 'DATASET UNSTRUCTURED_GRID'
      write(logic,fmt='(A6,1x,I7,1x,A5)') 'POINTS',numel*3,'float'
      do n = 1,numel
         iregion=ix(nen1,n)
         i1=ix(1,n)
         i2=ix(2,n)
         i3=ix(3,n)
          x1=x(1,i1)
          x2=x(1,i2)
          x3=x(1,i3)
          y1=x(2,i1)
          y2=x(2,i2)
          y3=x(2,i3)
          write(logic,'(3e13.5)') x1,y1,0.0
          write(logic,'(3e13.5)') x2,y2,0.0
          write(logic,'(3e13.5)') x3,y3,0.0
      enddo
      write(logic,*)
      write(logic,'(A5,1X,I7,1X,I7)') 'CELLS', numel, 4*numel

      do n=1,numel
         i=(n-1)*3-1
         write(logic,'(4i6)') 3,i+1,i+2,i+3
      enddo
      write(logic,'(A10,1X,I7)') 'CELL_TYPES', numel

      do n=1,numel
         write(logic,fmt='(5(1x,I7))') 5
      end do

      write(logic, *) 'CELL_DATA', numel
      write(logic, *) 'SCALARS DB integer'
      write(logic, *) 'LOOKUP_TABLE default'

      do n=1,numel
         write(logic, '(i5)') ix(nen1,n)
      end do

      close(logic)
      return
      end
************************************************************************
      subroutine checktime
      use mod_global
      implicit none
      write(*,*) '**** Current time is:',time
      write(*,*) '    and time step is:',dt
      end subroutine checktime
