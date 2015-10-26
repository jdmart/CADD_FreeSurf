c**---------------------------------------------------------------
c**    delaunay  :  generate the mesh
c**
c**   Algorithm :-
c**        allocate temp storage
c**        call contri
c**        write a mesh file for debugging purposes
c**        deallocate temp storage
c--
      subroutine delaunay(id,x,ix,f,b,itx)
      use mod_global
      use mod_boundary
      implicit none
c
      integer id(ndf,*),ix(nen1,*)
      integer itx(3,*)
      double precision x(nxdm,*),b(ndf,*),f(ndf,*)
c
      integer, pointer :: list(:),w(:),ixnew(:,:)
      integer i,j
      character*80 filename

      ! Allocate temporary storage
      allocate (list(numnp))
      allocate (w((numnp+3)*2))
      allocate (ixnew(3,(2*numnp+1)))

      ! create plot file of mesh
      filename='out/before.tec'
      call meshout(filename,nen1,nxdm,ndf,x,b,ix,numel,numnp)

      do i=1,numnp
         list(i)=i
      enddo
      call contri(numnp,numnp,nce,ncb,elist,x
     1     ,nxdm,list,w,ixnew,itx,numel)

      if (numel.gt.maxel) then
         write(6,'(///''  **error** detected by subroutine delaunay''
     1        /''  **number of elements exceeds maximum''
     2        /''  **elements requested = '',i5
     3        /''  **elements maximum   = '',i5)') numel,maxel
         stop
      end if
c
c     this is necessary because contri expects ix to be (3,numel) but we
c     have (4,numel).
      
      do i = 1,numel
         do j = 1,3
            ix(j,i) = ixnew(j,i)
         end do
      end do
      neq = ndf*numnp

      ! create plot file of mesh
      filename='out/after.tec'
      call meshout(filename,nen1,nxdm,ndf,x,b,ix,numel,numnp)

      ! Deallocate storage
      deallocate(w,ixnew,list)

      return
      end
c**---------------------------------------------------------------
c**     meshout  :  dump a mesh to a tecplot file
c**
c**   Non-Obvious Parameters :-
c**     filename   (in) : output file name
c**           ne   (in) : number of elements
c**           np   (in) : number of nodes
c--
      subroutine meshout(filename,nen1,nxdm,ndf,x,b,ix,ne,np)
      use mod_boundary
      use mod_file
      implicit none
c
      integer nen1,nxdm,ndf,ix(nen1,*),ne,np
      double precision x(nxdm,*),b(ndf,*)
      character*80 filename
c
      integer logic,i,j
      call iofile(filename,'formatted  ',logic,.true.)
      if (ne.gt.0.and.np.gt.0) then
         write(logic,*)'zone n=',np,', e=',ne,',f=fepoint,et=triangle'
      else if (np.gt.0) then
         write(logic,*)'zone'
      endif
      do i=1,np
         write(logic,'(4e14.5)') x(1,i),x(2,i),x(1,i)+b(1,i),x(2,i)+b(2
     $        ,i)
      enddo
      do i=1,ne
         write(logic,'(3i6)') (ix(j,i),j=1,3)
      enddo
      if (nce.eq.0) return
      write(logic,*)'zone n=',nce*2,', e=',nce,',f=fepoint,et=triangle'
      do i=1,nce
         j=elist(1,i)
         write(logic,'(4e14.5)') x(1,j),x(2,j),x(1,j)+b(1,j),x(2,j)+b(2
     $        ,j)
         j=elist(2,i)
         write(logic,'(4e14.5)') x(1,j),x(2,j),x(1,j)+b(1,j),x(2,j)+b(2
     $        ,j)
      enddo
      do i=1,nce
         write(logic,'(3i6)') 2*i,2*i-1,2*i
      enddo
      close(logic)
      end
