C Create free surface for indentation 
C Create atoms for rigid indenter 
C 
C 
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
c build a blunt *center, crack
      subroutine mp01(id,x,ix,f,b,itx)
      use mod_grain
      use mod_global
      use mod_file
      use mod_boundary
      use mod_crack
      use mod_material
      use mod_dd_slip
      implicit none
      include "./Disl/disl_parameters.par"
!     Input variables
      integer id, ix, itx(3,*)
      double precision x,f,b
      dimension id(ndf,1),x(nxdm,1),ix(nen1,1),f(ndf,1),b(ndf,1)

!     common declarations
      common/debugger/debug
      logical debug
      double precision tol

!     local variables
      logical n1,n2,n3,m1,m2,m3, useDetectionBand
      double precision, pointer:: nodeAngle(:)
      integer nXRegions,nYRegions, icell,ndxdy,numNodes,numx,numy,
     $     iSpace,ixrem,iyrem,iGrain,coincidentNodes,nodeStart,inode,
     $     nsort,np1,numnp0, countLast, i,j, k, node1, node2, node3
      integer nr1,nr2,nr3,nr4
      double precision XMax(0:20),YMax(0:20),nodeSite(3),
     $     dx,dy,dxdy,xxmin,xxmax,yymin,yymax,yyMinOrig,
     $     delx,dely,xx,yy, minDb, rcutmesh,Dist2, large
      data tol /1.d-6/
      logical placeNode,top,bot,left,right,mirror
      
      type Region
      	double precision xmin, xmax, ymin, ymax
      end type Region
      TYPE(Region) atomRegion, detectionBand, innerRegion, 
     $	mirrorAtomRegion, simulationCell
      logical insideRegion, inside      
      
      integer temp_slip, ii, jj, islp,iii
      double precision xslip_start, xslip_end, yslip_start, yslip_end
      double precision slip_angle(3), xxx1, yyy1, dxslip, dyslip
      double precision xslp1, xendslp1, ZBQLU01, rr1, lnuc, sn, smax
      integer*4 timeArray(3) 

C! VBS added this to read in tolerances
      common/elemconv/numelold
      integer lcrack,numelold,numnpc
cc--JS: Specially for crack asymmetry
      integer dwnumx, dwnumy, dwfactorx, dwfactory
      integer logic
      character*80 filename
cc--JS: dwfactor 2 for asym and 1 for sym      
      dwfactorx=1
      dwfactory = 1
cc
      write(*,*)
      write(*,*) 'Generating mesh with free surface on top indentation'
      write(*,*)


      x0crack = 0.01 
      y0crack = 0.01
!
      XMax(0)=-1.
      YMax(0)=-1.
      
      
! Read mesh data            
      read(5,*) nXRegions,nYRegions
      read(5,*) (XMax(i),i=1,nXRegions)
      read(5,*) (YMax(i),i=1,nYRegions)
      read(5,*) minDb
      read(5,*) rcutmesh
      read(5,*) x0crack, y0crack
      read(5,*) pad_width

      mirror=.false.      
     
!       if(rcutmesh.lt.0.d0) then
!          rcutmesh=abs(rcutmesh)
! 	write (6,*) 'Mirror not allowed'
! 	stop
!       endif

      
      do k=min(nXRegions,nYRegions)+1,max(nXRegions,nYRegions)
         if (nXRegions.lt.nYRegions) then
            XMax(k)=XMax(nXRegions)
         else
            YMax(k)=YMax(nYRegions)
         endif
      enddo

      
!     Extract lattice data
!     Normally, the second atom is dx away from the first due to the
!     way cell is sorted.
!     if there is only 1 atom on the lowest y-plane of cell, then dx is the
!     cell width in the x dirn.
      dx=grains(1)%cell(1,2)
      if (grains(1)%cell(2,2).gt.tol) dx=grains(1)%dcell(1)
c
      do icell=1,grains(1)%ncell
         if (grains(1)%cell(2,icell).gt.tol) then
            dy=grains(1)%cell(2,icell)
            dxdy=grains(1)%cell(1,icell)
            goto 9
         endif
      enddo
 9    continue
 
 
      if (dxdy.ne.0.) then
         ndxdy=nint(dx/dxdy)
      else
         ndxdy=1
      endif
      
      print*,'Generating nodes'
! Generate the nodes in each box, the mesh is coarsened with increasing
! distance from the center of the box.
      numx=int(XMax(nXRegions)/dx)+1
      numy=int(YMax(nYRegions)/dy)+1
      numNodes=0
      do 10 i=-numx,numx
         xx=i*dx
         if (abs(xx).gt.XMax(nXRegions)) go to 10
	 
         do 20 j=-numy,0
	 
            yy=j*dy
            if (abs(yy).gt.YMax(nYRegions)) go to 20
	    
!! Determine the region in which the node is being placed	    
            do k=max(nXRegions,nYRegions)-1,1,-1
cc--JS: Symmetric geometry-Commented if ther is crack
!               if ( (abs(xx).gt.XMax(k))
!     $              .or.(abs(yy).gt.YMax(k)) ) go to 21
cc--JS: Asymetry due to the presence of crack
                if ( (xx.gt.XMax(k))
     $              .or.(yy.gt.YMax(k)/dwfactory) ) go to 21
                if ( (xx.lt.-XMax(k)/dwfactorx)
     $              .or.(yy.gt.YMax(k)/dwfactory) ) go to 21
                if ( (xx.gt.XMax(k))
     $              .or.(yy.lt.-YMax(k)) ) go to 21
               if ( (xx.lt.-XMax(k)/dwfactorx)
     $              .or.(yy.lt.-YMax(k)) ) go to 21
            enddo
 21         iSpace=2**k
            if(iSpace.eq.0) iSpace=1

!! Decide if a node should be placed in this region
            ixrem=mod(abs(i),iSpace)
            iyrem=mod(abs(j),iSpace)
            placeNode=(ixrem+iyrem).eq.0
            if (.not.placeNode) go to 20

! Assign the node a position
            numNodes=numNodes+1
            x(1,numNodes)=xx+mod(j,ndxdy)*dxdy
            x(2,numNodes)=yy
! Qu modification to remove extra layers of atoms around crack faces
c$$$            if((j.eq.0.and.x(1,numNodes).lt.-0.01*dx).or.
c$$$     &           (j.eq.1.and.x(1,numNodes).lt.-0.01*dx) .or.
c$$$     &      (j.eq.-1.and.x(1,numNodes).lt.-dx) )then
c$$$c     ***************************************************************
c$$$c     Add this line for a sharp crack 
c$$$c     Comment this line and uncomment to the last 3 lines to 
c$$$c     for blunt crack with 3 layers missing 
c$$$c     ***************************************************************
c$$$c            if(j.eq.0.and.x(1,numNodes).lt.-0.1*dx) then 
c$$$              numNodes=numNodes-1
c$$$            endif

c$$$ Chakravarthy edits for indentation
            
c$$$
c$$$!     FOR NOW ONLY, assign nodes to the crack faces. Leo mentions
c$$$! this is a hack. But this appears legitimate.
c$$$!     Add back nodes according to desired crack shape
c$$$            if(j.eq.0.and.x(1,numNodes).lt.-0.01*dx) then
c$$$              numNodes = numNodes+2
c$$$              x(1,numNodes-1) = xx+mod(j,ndxdy)*dxdy
c$$$!              x(2,numNodes-1) = yy+2*dy
c$$$              x(2,numNodes-1) = yy+dy
c$$$              x(1,numNodes) = xx+mod(j,ndxdy)*dxdy
c$$$!              x(2,numNodes) = yy-2*dy 
c$$$              x(2,numNodes) = yy 
c$$$!              print *, i, x(1, numNodes)
c$$$              if(i.eq.-320) then
c$$$                countLast = numNodes-1
c$$$              endif
c$$$            endif
c$$$C Qu modification ends

 20      end do
 10    end do
        print *, 'Moving nodes to nearest atomic site'

 
c       Move nodes to the nearest atomic sites
        large = 1.e30
        xxmax=-large
	    xxmin=large
	    yymax=-large
	    yymin=large
	    do i=1,numNodes
	
		    nodeSite(1)=x(1,i)
		    nodeSite(2)=x(2,i)
		    if (i.ne.countLast) then
                       call NearestBsite(nodeSite,1,.false.,x(1,i)
     $                      ,iGrain)
		    end if

!!     find xxmax, xxmin, etc.
		    xxmax = max(x(1,i), xxmax)
		    xxmin = min(x(1,i), xxmin)
		    yymax = max(x(2,i), yymax)
		    yymin = min(x(2,i), yymin)
c           print *, i, x(1,i), x(2,i)
            if (x(1,i) .eq. 0.0d0) then
                if (x(2,i) .eq. 0.0d0) then
                    temp_slip = i
                end if
            end if
	    end do
	   xslip_start = x(1,temp_slip+1)
       yslip_start = (x(2,temp_slip)+x(2,temp_slip+1))/2.d0

       xxx1 = x(1,temp_slip+1)-x(1,temp_slip)
       yyy1 = x(2,temp_slip+1)-x(2,temp_slip)
       slip_angle(1) = atan2(yyy1,xxx1)
       dxslip = abs(xxx1)*2.d0
       dyslip = abs(yyy1)
       print *, xxx1, yyy1, slip_angle(1)
       print *, 'Slip plane start and end', temp_slip
       print *, 'x_start', xslip_start, x(1,temp_slip+1)
       print *, 'y_start', yslip_start


      xxx1 = x(1,temp_slip-1)-x(1,temp_slip)
      yyy1 = abs(x(2,temp_slip-1)-x(2,temp_slip))
      slip_angle(2) = atan2(yyy1,xxx1)
      print *, xxx1, yyy1, slip_angle(2)
      slip_angle(3) = 0.0d0
      print *, 'DX', dxslip
      print *, 'DY', dyslip
      slip_angle(3) = 0.0d0

      
      simulationCell%xmin = xxmin
      simulationCell%xmax = xxmax
      simulationCell%ymin = yymin
      simulationCell%ymax = yymax
c     hard coded for 1 grain
            
      if(mirror) then
         numnp=numNodes
         do i=1,numnp
	 
            if(x(2,i).lt.yymin+tol) then
               nodeSite(1)=x(1,i)
               nodeSite(2)=-136.85
               call NearestBsite(nodeSite,1,.false.,x(1,i),iGrain)
            endif
	    
            numNodes=numNodes+1
            nodeSite(1)=x(1,i)
            nodeSite(2)=2*yymin-x(2,i)
            call NearestBsite(nodeSite,1,.false.,x(1,numNodes),iGrain)
	    
         enddo
	 
         yyMinOrig=yymin
         yymin=2*yymin-yymax
	 
	     simulationCell%ymin = yymin
      endif


! Remove coincident nodes
      print*, 'Removing coincident nodes'
      coincidentNodes=0
      numnp=numNodes
      do 30 i=numnp,2,-1
         if (i.gt.numNodes) goto 30
         do j=1,i-1
            delx=abs(x(1,i)-x(1,j))
            dely=abs(x(2,i)-x(2,j))
            if (delx+dely.lt.2.*tol) then
               x(1,j)=x(1,numNodes)
               x(2,j)=x(2,numNodes)
               numNodes=numNodes-1
               coincidentNodes=coincidentNodes+1
            endif
         enddo
 30   enddo
 
      if (coincidentNodes.ne.0) then
	write(6,*) coincidentNodes,' coincident nodes removed'
      endif
      
      numnp=numNodes
      if (numnp.gt.maxnp) then
         write(6,*)'***ERROR: Insufficient storage for nodes'
         write(6,*) '         numnp = ',numnp
         stop
      endif
      write(6,*) 'Total nodes: numnp = ',numnp


!     Apply boundary conditions and define the boundary for the
!     triangulator

      nce=0
      nodeStart=0

c$$$      print*, 'Detecting the lower crack ledge'
c$$$!     FIND THE LOWER CRACK LEDGE
c$$$      do i =1, numnp
c$$$! Qu modification to remove extra layers of atoms around crack faces
c$$$         if(dabs(x(2,i)).lt.tol.and.x(1,i).lt.0.1*dx) then
c$$$!         if(dabs(x(2,i)+dy).lt.tol.and.x(1,i).lt.0.1*dx) then
c$$$!        if(dabs(x(2,i)+2*dy).lt.tol.and.x(1,i).lt.-dx) then
c$$$! Qu modification ends
c$$$!        if(dabs(x(2,i)).lt.tol.and.x(1,i).lt.XMax(1)-4*rcutmesh) then
c$$$!        if(dabs(x(2,i)).lt.tol.and.x(1,i).lt.-XMax(1)+2*dx) then
c$$$        nce=nce+1
c$$$        if (nce.gt.NCEMAX) then
c$$$           if(nce.gt.NCEMAX) call IncreaseElist(100)
c$$$        endif
c$$$        elist(1,nce)=i
c$$$        endif
c$$$      enddo
c$$$
c$$$! sort the lower ledge so that is goes CW (from right to left).
c$$$      allocate(nodeAngle(numnp))
c$$$      nodeAngle=0.
c$$$      do i=1,nce
c$$$         inode=elist(1,i)
c$$$         nodeAngle(inode)= datan2(dy,x(1,inode))
c$$$      enddo
c$$$      nsort=nce
c$$$      call qsortr(nsort,elist(1,1),nodeAngle,1,1,1.d0)
c$$$      deallocate(nodeAngle)
c$$$      
c$$$! remove the last node from the list
c$$$Cc!!! Bill's changes!!!!
c$$$      nce=nce-1
c$$$      nodeStart = nce
c$$$
      
!     find all external boundary nodes
      simulationCell%xmin = simulationCell%xmin + 10.d0
      simulationCell%xmax = simulationCell%xmax - 10.d0
      simulationCell%ymin = simulationCell%ymin + 10.d0
      simulationCell%ymax = simulationCell%ymax - 0.01d0

      print*, 'Detecting the outer cell boundary', numnp
      print *,  simulationCell%xmin,  simulationCell%xmax
      print *,  simulationCell%ymin,  simulationCell%ymax

      
      do 77 i=1,numnp

         top=(x(2,i) .gt. simulationCell%ymax)
         bot=(x(2,i) .lt. simulationCell%ymin)
         left=(x(1,i) .lt. simulationCell%xmin)
         right=(x(1,i) .gt. simulationCell%xmax)

!     store all boundary points, but put crack faces at the beginning of
!     elist.  While you are at it, apply the b.c.s

c$$$         if (top.or.bot.or.right.or.left) then
c$$$            nce=nce+1
c$$$            if (nce.gt.NCEMAX) then
c$$$               if(nce.gt.NCEMAX) call IncreaseElist(100)
c$$$            endif
c$$$            elist(1,nce)=i

! apply the b.c's
c$$$            if (top.or.bot.or.right.or.left) then
c$$$               id(1,i)=1
c$$$               id(2,i)=1
c$$$               print *, 'BCs on  node',i, x(1,i), x(2,i)
c$$$            endif
	    	    
c$$$         endif
 77   continue

c$$$
c$$$! sort the boundary so that is goes CW.
c$$$      allocate(nodeAngle(numnp))
c$$$      nodeAngle=0.
c$$$      do i=nodeStart+1,nce
c$$$         inode=elist(1,i)
c$$$! YET ANOTHER HACK
c$$$         nodeAngle(inode)=datan2(x(2,inode)-tol,x(1,inode))
c$$$      enddo
c$$$      nsort=nce-nodeStart
c$$$      call qsortr(nsort,elist(1,nodeStart+1),nodeAngle,1,1,1.d0)
c$$$      deallocate(nodeAngle)
c$$$
c$$$! remove the last node from the list
c$$$      nce=nce-1
c$$$      nodeStart = nce
c$$$
c$$$
c$$$!     FIND THE UPPER CRACK LEDGE
c$$$      print*, 'Detecting the upper crack ledge'
c$$$      do i =1, numnp
c$$$! Qu modification to remove extra layers of atoms around crack faces
c$$$        if(dabs(x(2,i)-1.0*dy).lt.tol.and.x(1,i).lt.0.1*dx) then
c$$$!        if(dabs(x(2,i)-2.0*dy).lt.tol.and.x(1,i).lt.0.1*dx) then
c$$$!        if(dabs(x(2,i)-dy).lt.tol.and.x(1,i).lt.XMax(1)-4*rcutmesh) then
c$$$!        if(dabs(x(2,i)-dy).lt.tol.and.x(1,i).lt.-XMax(1)+2*dx) then
c$$$        nce=nce+1
c$$$        if (nce.gt.NCEMAX) then
c$$$           if(nce.gt.NCEMAX) call IncreaseElist(100)
c$$$        endif
c$$$        elist(1,nce)=i
c$$$        endif
c$$$      enddo
c$$$
c$$$! sort the upper ledge so that is goes CW (from left to right).
c$$$      allocate(nodeAngle(numnp))
c$$$      nodeAngle=0.
c$$$      do i=nodeStart+1,nce
c$$$         inode=elist(1,i)
c$$$         nodeAngle(inode)= -datan2(x(2,inode),x(1,inode))
c$$$      enddo
c$$$      nsort=nce-nodeStart
c$$$      call qsortr(nsort,elist(1,nodeStart+1),nodeAngle,1,1,1.d0)
c$$$      deallocate(nodeAngle)
c$$$
c$$$! finish defining the boundary
c$$$      do i=1,nce-1
c$$$         elist(2,i)=elist(1,i+1)
c$$$      enddo
c$$$      elist(2,nce)=elist(1,1)
      ncb=nce

c$$$      do i=1,numnp
c$$$        if(dabs((x(2,i))-dy).lt.tol.and.x(1,i).gt.-dx
c$$$     &     .and.x(1,i).lt.0.5*dx)then
c$$$          nce=nce+1
c$$$          if(nce.gt.NCEMAX) then
c$$$           if(nce.gt.NCEMAX) call IncreaseElist(100)
c$$$          endif
c$$$
c$$$          elist(1,nce)=i
c$$$          ncb=nce
c$$$          elist(2,nce)=elist(1,1)
c$$$          elist(2,nce-1)=elist(1,nce)
c$$$        endif
c$$$      enddo
c$$$      do i=1,numnp
c$$$        if(dabs(x(2,i)).lt.tol.and.x(1,i).gt.-0.5*dx
c$$$     &     .and.x(1,i).lt.0.5*dx)then
c$$$          nce=nce+1
c$$$          if(nce.gt.NCEMAX) then
c$$$           if(nce.gt.NCEMAX) call IncreaseElist(100)
c$$$          endif
c$$$
c$$$          elist(1,nce)=i
c$$$          ncb=nce
c$$$          elist(2,nce)=elist(1,1)
c$$$          elist(2,nce-1)=elist(1,nce)
c$$$        endif
c$$$      enddo
c$$$
c$$$      do i=1,numnp
c$$$        if(dabs((x(2,i))+dy).lt.tol.and.x(1,i).gt.-dx
c$$$     &     .and.x(1,i).lt.tol)then
c$$$          nce=nce+1
c$$$          if(nce.gt.NCEMAX) then
c$$$           if(nce.gt.NCEMAX) call IncreaseElist(100)
c$$$          endif
c$$$          
c$$$          elist(1,nce)=i
c$$$          ncb=nce
c$$$          elist(2,nce)=elist(1,1)
c$$$          elist(2,nce-1)=elist(1,nce)
c$$$        endif
c$$$      enddo

c     Triangulate, sets all elements to material 1 for this mesh
      print*, 'Triangulating'
      numnpc = numnp
      call delaunay(id,x,ix,f,b,itx)
      print*, 'Done'

      write(*,*) 'BEFORE adding overlap'
      write(6,*) 'Number of nodes: numnp = ',numnp
      write(6,*) 'Number of elements: numel = ',numel
      if(numel.gt.maxel) stop 'too many elements'
      if(numnp.gt.maxnp) stop 'too many nodes'
      
      


!     Find max/min of atomistic region
	innerRegion%xmin = -Xmax(1)/dwfactorx
	innerRegion%xmax = Xmax(1)
	innerRegion%ymin = -Ymax(1)
	innerRegion%ymax = Ymax(1)/dwfactory	
	call findAtomRegionSize(x,dx, dy, tol, innerRegion, atomRegion)


       do 75 i = 1,numel
         ix(nen1,i) = 1
 75   continue
     
 
 
!     Find continuum region elements
      if (mirror) then
		mirrorAtomRegion%xmin = atomRegion%xmin
		mirrorAtomRegion%xmax = atomRegion%xmax
		mirrorAtomRegion%ymin = 2*yyMinOrig-atomRegion%ymin
		mirrorAtomRegion%ymax = 2*yyMinOrig-atomRegion%ymax
      endif
	
      
!	Check if each node of any continuum element is in the
!	atomistic region      
      do i=1,numel
	
       
      	node1 = ix(1,i)
	  node2 = ix(2,i)
	  node3 = ix(3,i)
	
!! Determine if any node is in the atomistic region      
         n1=.not.insideRegion(x(1:2,node1),atomRegion)
         n2=.not.insideRegion(x(1:2,node2),atomRegion)
         n3=.not.insideRegion(x(1:2,node3),atomRegion)

!! Determine if any node is in the mirrored atomistic region
         if(mirror) then	 	     
             m1=.not.insideRegion(x(1:2,node1),mirrorAtomRegion)
             m2=.not.insideRegion(x(1:2,node2),mirrorAtomRegion)
             m3=.not.insideRegion(x(1:2,node3),mirrorAtomRegion)	     
         else
            m1=.true.
            m2=.true.
            m3=.true.
         endif
	 
         if((n1.and.n2.and.n3).and.(m1.and.m2.and.m3)) ix(nen1,i)=0
      enddo

 
           
      
!     Add pad atoms in the interface region. Note that elements are 
!     not needed in this region, so just add atoms.
!      go to 1234
      print*, 'Adding pad atoms', pad_width
      numNodes=numnp
!      XMax(2)=XMax(1)+2.d0*rcutmesh
!      YMax(2)=YMax(1)+2.d0*rcutmesh
      XMax(2)=XMax(1) + pad_width
      YMax(2)=YMax(1) + pad_width
      numx=int(XMax(2)/dx)+1
      numy=int(YMax(2)/dy)+1
      dwnumx=int(XMax(2)/dwfactorx/dx)+1
      dwnumy=int(YMax(2)/dwfactory/dy)+1
      dwnumy = 0
      do i=-dwnumx,numx
         xx=i*dx
         do j=-numy,dwnumy
            yy=j*dy
	    
            numNodes=numNodes+1
	    
            nodeSite(1)=xx+mod(j,ndxdy)*dxdy
            nodeSite(2)=yy
	    
            call NearestBsite(nodeSite,1,.false.,x(1,numNodes),iGrain)

!!          Skip this node if it is the atomistic region
            if(abs(j).lt.0 .and. i.lt.0)then
              numNodes=numNodes-1
            elseif(insideRegion(x(1:2,numNodes),atomRegion)) then
               numNodes=numNodes-1
               endif
	    
            if (mirror) then
               nodeSite(1)=x(1,numNodes)
               nodeSite(2)=2*yyMinOrig-x(2,numNodes)
               numNodes=numNodes+1
               call NearestBsite(nodeSite,1,.false.,x(1,numNodes),
     $				iGrain)
            endif	    
         enddo
      enddo
      print*, 'Done adding pad atoms'      
            
      np1=numnp+1
      numnp0=numnp
      numnp=numNodes
      if(numel.gt.maxel) stop 'Too many elements'
      if(numnp.gt.maxnp) stop 'Too many nodes'

      
      
!     Determine nodal character -- continuum, interface, atomistic etc.
      num2Dnode=numnp
      call StatusCalc(x,ix,.true.)

      
!     remove the coincident nodes on the interface
      print*, 'Removing coincident nodes near the interface'
      do i=numnp,np1,-1
         do j=1,numnp0
            if(IsRelaxed(j).ne.0) then
               if(Dist2(x(1:2,i),x(1:2,j),2).lt.1.e-6) then
                  f(1:ndf,i)=f(1:ndf,numnp)
                  x(1:nxdm,i)=x(1:nxdm,numnp)
                  numnp=numnp-1
                  go to 300
               endif
            endif
         enddo
 300     continue
      enddo
      num2Dnode=numnp
      if(numperiodz.gt.1)then
         call IncAtoms(x)
      endif

      if(numnp.gt.maxnp) stop 'Too many nodes'

c      allocate(atomSpecie(numnp))
c      do i=1,numnp
c         atomSpecie(i)=1
c      enddo

C---   Qu added on 08/26/2005
C----  For multimaterils purpose

      numnp0=numnp
c#########################################################################
C------- insert H atoms in the crack tip region 
      print*, 'nmaterials', nmaterials 
      if(nmaterials.gt.1)	then
	    call insertHatom(x)
	    print*, ' Adding interstitial atom'
      endif
 
      if(numnp.gt.maxnp) stop 'Too many nodes'

      allocate(atomSpecie(numnp))
      do i=1,numnp0
         atomSpecie(i)=1
      enddo
      do i=numnp0+1,numnp
         print*, '# of H atom', i 
         atomSpecie(i)=2
      write(6,*), 'int', x(1,i), '',x(2,i), '', x(3,i) 
      end do


      iHnumber=i 
!     M modif for removing one atom
      print*, 'done'


 1234 write(*,*) 'Final mesh size'
      write(6,*) 'Total nodes: numnp = ',numnp
      write(6,*) 'Total elements: numel = ',numel
      write(6,*) 'rcutmesh = ', rcutmesh
      numnpp1=-1


!     Create a detection band
!     Identify a path, defined by a closed polygon, ccw around the
!     vertices, along which the detection band elements will be placed.
!     If the atomistic continuum interface is not a closed path, define
!     this polygon to extend outside the mesh and surround the atomistic
!     region.


	   useDetectionBand = .true.
	   if (.not. useDetectionBand) then
		 ndbpoly = 0
		    return
		end if
        print *, 'Using a detection band'

!     There is one detection band with 4 points
!     Qu modified detection band rings starts
	  detectionBand%xmin = atomRegion%xmin + rcutmesh
	  detectionBand%xmax = atomRegion%xmax - rcutmesh
	  detectionBand%ymin = atomRegion%ymin + rcutmesh
	  detectionBand%ymax = atomRegion%ymax - rcutmesh
        minDb=int(minDb/dy)*dy
	   nr1 = int((abs(detectionBand%xmin)-minDb)/dx)+1
	   nr2 = int((abs(detectionBand%xmAX)-minDb)/dx)+1
c$$$        nr2 = nr1
	   nr3 = int((abs(detectionBand%ymin)-minDb)/dy)+1
	   nr4 = int((abs(detectionBand%ymax)-minDb)/dy)+1
        nr4 = nr3

	   ndbpoly=max(nr1,nr2,nr3,nr4)+1
        ndbpoly = 2
        print *, 'Detection band', nr1, nr2, nr3, nr4, ndbpoly
        call allocate_db
 
      do i=1,ndbpoly
         detectionBand%xmin = atomRegion%xmin + rcutmesh + (i-1)*dx
c$$$         detectionBand%xmax = atomRegion%xmax - 2.0*rcutmesh
c$$$         detectionBand%ymin = atomRegion%ymin + 2.0*rcutmesh
c$$$         detectionBand%ymax = atomRegion%ymax - 2.0*rcutmesh
         
         detectionBand%xmax = atomRegion%xmax - rcutmesh - (i-1)*dx
         detectionBand%ymin = atomRegion%ymin + rcutmesh + (i-1)*dy
         detectionBand%ymax = atomRegion%ymax - rcutmesh 
c$$$         detectionBand%ymax = atomRegion%ymax - rcutmesh 

c$$$
c$$$         detectionBand%xmax = min( minDb+(i-1)*dx
c$$$     $                            ,atomRegion%xmax - rcutmesh)
c$$$         detectionBand%ymin = max(-minDb-(i-1)*dy
c$$$     $                            ,atomRegion%ymin + rcutmesh)
c$$$         detectionBand%ymax = min( minDb+(i-1)*dy
c$$$     $                            ,atomRegion%ymax - rcutmesh)
      
   
         dbpoly(1,1,i)=detectionBand%xmin
         dbpoly(2,1,i)=detectionBand%ymin
      
         dbpoly(1,2,i)=detectionBand%xmax
         dbpoly(2,2,i)=detectionBand%ymin
      
         dbpoly(1,3,i)=detectionBand%xmax
         dbpoly(2,3,i)=detectionBand%ymax
      
         dbpoly(1,4,i)=detectionBand%xmin
         dbpoly(2,4,i)=detectionBand%ymax

c$$$     For this particular detection band the first layer is closest
c$$$     to the boundary 
c$$$     Any other logic can be used ... 
         if (i == 1) then 
            dbboundnear(i) = .true. 
         end if


c$$$         print *, 'Detection band region', i,
c$$$     $        dbpoly(1,1,i), dbpoly(2,1,i)
c$$$         print *,
c$$$     $        'Detection band region', i,dbpoly(1,2,i), dbpoly(2,2,i)
c$$$         print
c$$$     $        *, 'Detection band region', i,dbpoly(1,3,i), dbpoly(2,3,i)
c$$$         print
c$$$     $        *, 'Detection band region', i,dbpoly(1,4,i), dbpoly(2,4,i)
         print *, 'Detection Band Region', i
         write(*, '(8f10.3)') dbpoly(1,1,i), dbpoly(2,1,i), dbpoly(1
     $        ,2,i),dbpoly(2,2,i), dbpoly(1,3,i), dbpoly(2,3,i),
     $        dbpoly(1,4,i), dbpoly(2,4,i)

      end do
! Qu modified detection band rings ends
      
      write(6,*) 'width of the detection band: rcutmesh = ', rcutmesh

      call gen_slip_planes(simulationCell%xmin,
     & simulationCell%xmax,simulationCell%ymin,
     & simulationCell%ymax,
     & atomRegion%xmin, atomRegion%xmax,
     & atomRegion%ymin, atomRegion%ymax,slip_angle,
     &     grains(1)%dcell(1), xslip_start,yslip_start,dxslip,dyslip
     $     ,pad_width)

      x_move_mesh = 5.d0*dxslip
      MoveMesh = .false.
      Moved = .false. 

      return
      end
      
      

	 
************************************************************************       
!	  Checks to see if a 2D point x is located inside thisRegion.
      logical function insideRegion(x, thisRegion)
      implicit none

      type Region
      	double precision xmin, xmax, ymin, ymax
      end type Region
      TYPE(Region) thisRegion
      double precision x(2)
      
      insideRegion = (x(1).gt. thisRegion%xmin .and.
     $			x(1) .lt. thisRegion%xmax .and.
     $			x(2) .gt. thisRegion%ymin .and.
     $			x(2) .lt. thisRegion%ymax)
      
      return
      end
      
************************************************************************
      subroutine qsortr(n,list,xkey,nxdm,ind,sign)
      implicit double precision (a-h,o-z)
c
      dimension xkey(nxdm,1)
      integer list(2,*),n,ll,lr,lm,nl,nr,ltemp,stktop,maxstk
c
      parameter (maxstk=32)
c
      integer lstack(maxstk), rstack(maxstk)
c
      ll = 1
      lr = n
      stktop = 0
   10 if (ll.lt.lr) then
         nl = ll
         nr = lr
         lm = (ll+lr) / 2
         guess = sign*xkey(ind,list(1,lm))
c
c     Find xkeys for exchange
c
   20    if (sign*xkey(ind,list(1,nl)).lt.guess) then
            nl = nl + 1
            go to 20
         endif
   30    if (guess.lt.sign*xkey(ind,list(1,nr))) then
            nr = nr - 1
            go to 30
         endif
         if (nl.lt.(nr-1)) then
            ltemp    = list(1,nl)
            list(1,nl) = list(1,nr)
            list(1,nr) = ltemp
            nl = nl + 1
            nr = nr - 1
            go to 20
         endif
c
c     Deal with crossing of pointers
c
         if (nl.le.nr) then
            if (nl.lt.nr) then
               ltemp    = list(1,nl)
               list(1,nl) = list(1,nr)
               list(1,nr) = ltemp
            endif
            nl = nl + 1
            nr = nr - 1
         endif
c
c     Select sub-list to be processed next
c
         stktop = stktop + 1
         if (nr.lt.lm) then
            lstack(stktop) = nl
            rstack(stktop) = lr
            lr = nr
         else
            lstack(stktop) = ll
            rstack(stktop) = nr
            ll = nl
         endif
         go to 10
      endif
c
c     Process any stacked sub-lists
c
      if (stktop.ne.0) then
         ll = lstack(stktop)
         lr = rstack(stktop)
         stktop = stktop - 1
         go to 10
      endif
c
      return
      end

      
      
!     Find max/min of atomistic region      
      subroutine findAtomRegionSize(atomCoord,dx, dy, tol, 
     $		innerRegion, atomRegion)
      use mod_global
      implicit none
      
      type Region
      	double precision xmin, xmax, ymin, ymax
      end type Region
      TYPE(Region) atomRegion, innerRegion      
      double precision atomCoord(ndf,*), dx, dy, tol

!	  local variables
      double precision xmin, xmax, ymin, ymax, x, y, large  
      integer iNode           
!	  functions
      logical insideRegion  
      
      large = 1.e30
            
      xmax=-large
      xmin= large
      
      ymax=-large
      ymin= large
      
      do iNode=1,numnp
      
         x=atomCoord(1,iNode)
         y=atomCoord(2,iNode)
	 
         if(insideRegion(atomCoord(1:2,iNode),innerRegion)) then	    
            xmax=max(xmax,x)
            ymax=max(ymax,y)
            xmin=min(xmin,x)
            ymin=min(ymin,y)	    
         endif
	 
      enddo
      
      xmax=xmax-dx+tol
      ymax=ymax+tol
      xmin=xmin+dx-tol
      ymin=ymin-tol
      
      xmax=xmax-dx
      xmin=xmin+dx
      ymax=ymax-dy
      ymin=ymin+dy
            
      atomRegion%xmin = xmin
      atomRegion%xmax = xmax
      atomRegion%ymin = ymin
      atomRegion%ymax = ymax
            
      return
      end
      
! Qu added begin
!
      subroutine IncAtoms(x)
      use mod_global
      use mod_grain
      implicit none

      double precision x
      dimension x(nxdm,1)

      integer i,j,numnp0

      numnp0=numnp
      do i=2,numperiodz
         do j=1,numnp0
            if(IsRelaxed(j).ne.0)then
               numnp=numnp+1
               IsRelaxed(numnp)=IsRelaxed(j)
               x(1:2,numnp)=x(1:2,j)
               x(3,numnp)=x(3,j)+(i-1)*grains(1)%dcell(3)
            endif
         enddo
       enddo

      return
      end     
!
! Qu added ends  
      
c	c*********************************************************************
      subroutine insertHatom(x)
      use mod_global
      use mod_grain
      implicit none

      double precision x
      dimension x(nxdm,1)

C----Local variables
      integer i,j,atomHnum,inHatoms
      character *80 input_file

      input_file = 'interstitial.inp'
      open(unit=10, file = input_file, status = 'old')
      read(10,*)inHatoms
cc--  JS
      NumTotH=inHatoms
      do i=1, inHatoms
         numnp=numnp+1
      print*, '# of H atom', numnp
         read(10,*)atomHnum,(x(j,numnp),j=1,3)
      enddo      
      return
      end    
 
************************************************************************
      
************************************************************************
      subroutine PerturbMesh(x,nxdm,numnp,perturb)
      end
      
      
! ************************************************************************
!       logical function inside(x,x1,x2,y1,y2)
!       implicit none
!       double precision x(2),x1,x2,y1,y2
!       inside=(x(1).gt.x1).and.(x(1).lt.x2).and.(x(2).gt.y1).and.(x(2).lt
!      $     .y2)
!       end
!             
	  double precision function getRandomNumber1()
	     implicit none
	     double precision :: random,  ZBQLU01
	
		random = -1.d0 + 2.d0 *ZBQLU01(0.0D0)
		
		if (random .gt. 0) then
			random = 1.d0
		else if (random .lt. 0) then
			random = -1.d0
		else if (dabs(random) .lt. 1.0d-6) then
			random = 0.d0
		endif
		getRandomNumber1 = random				
	     return
  	  end
