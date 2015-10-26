      module mod_dynamo 
!                                                                       
!     This include is for use with the "dynamo" solver only.            
!     It is the part of dyn84_2.inc that is needed by the dynamo solver,
!     but is not already included in poten_eam.inc                      
!                                                                       
!.....$Modified: Wed Jan 18 10:27:09 1995 by foiles $                   
      integer NEIMAX 
      parameter (NEIMAX=300) 
      double precision rneigh(NEIMAX),dneigh(3,NEIMAX) 
      integer jneigh(NEIMAX) 
      double precision, allocatable:: rold(:,:),dis(:,:)                &
     &     ,rdyn(:)                                                     
      integer, allocatable:: nnindx(:),nnlst(:,:),knbr(:),              &
     &     NeighborList(:,:), NumNeighbors(:)                           
      double precision perub(3),perlb(3),perlen(3),alat,xbound,ybound(2)&
     &     ,zbound(2)                                                   
      logical onlyoi,twoi,threei,fouri 
      integer ngtlst,nneimx,mxlstu,mxnnei,nforce 
      double precision rctsqn,dradn 
      integer nneigh,nneips,nmeth,newlst,nlstmx 
      contains 
                                                                        
!********************************************************************** 
      double precision function rv2(i,j,b,x) 
      use mod_global 
      integer i,j 
      double precision b(ndf,numnp),x(nxdm,numnp) 
      rv2=b(i,j)+x(i,j) 
      end function rv2 
!********************************************************************** 

!.....                                                                  
!.....this routine determines the neighbors of a given atom             
!.....                                                                  
      subroutine gneigh(i,b,x,NeedList) 
      use mod_global 
      implicit double precision (a-h,o-z) 
      dimension b(ndf,*),x(nxdm,*) 
      integer totNeighbors 
      logical NeedList, firstTime 
      data firstTime/.true./ 
                                                                        
!--   Vijay                                                             
      if (firstTime .eq. .true.) then 
         firstTime = .false. 
         allocate(NeighborList(MAXNEIGHBORS, numnp)) 
         allocate(NumNeighbors(numnp)) 
                                                                        
         do iAtom = 1, numnp 
!     print*, 'iAtom: ', iAtom                                          
            do j = 1, MAXNEIGHBORS 
               NeighborList(j, iAtom) = 0 
            enddo 
            NumNeighbors(iAtom) = 0 
         enddo 
      endif 
!--   Vijay	  	                                                         
                                                                        
                                                                        
!     print*, 'In gNeigh: iAtom: ', i                                   
      if(i.eq.numnpp1) then 
         rcutsq=CutoffR2(-1) 
      else 
         rcutsq=CutoffR2(1) 
      endif 
      rctsqn = (sqrt(rcutsq) + dradn)**2 
                                                                        
                                                                        
      totNeighbors = 0 
!.....                                                                  
!.....beginning of nmeth = 2                                            
!.....storage of neighbor indices                                       
!.....                                                                  
!.....newlst specifies whether the neighbor list is being updated       
!.....                                                                  
      if(SimStep .eq. 1) then 
         if (newlst.ne.0) goto 2500 
      endif 
                                                                        
!     C--Jun Song New: only update those atoms sepcified in chkdis      
      if(SimStep .gt. 1) then 
         if((newlst.ne.0).and.(UpdateNeigh(i).eq.1)) goto 2500 
      endif 
!     c-end                                                             
                                                                        
!.....                                                                  
!.....if here then use the neighbors found earlier                      
!.....                                                                  
 2001 continue 
      if (.not.NeedList) return 
!     c-JS         jend = nnindx(i) - nnindx(i-1)                       
      jend=nnindx(i) 
      do 2110 jtmp = 1,jend 
!     c-JS            knbr(jtmp) = nnlst(jtmp+nnindx(i-1))              
         knbr(jtmp) = nnlst(i,jtmp)	 
                                                                        
         dis(1,jtmp) = rv2(1,i,b,x) - rv2(1,knbr(jtmp),b,x) 
         dis(1,jtmp) = dis(1,jtmp)                                      &
     &        - perlen(1)*nint(dis(1,jtmp)/perlen(1))                   
         rdyn(jtmp) = dis(1,jtmp)**2 
         dis(2,jtmp) = rv2(2,i,b,x) - rv2(2,knbr(jtmp),b,x) 
         dis(2,jtmp) = dis(2,jtmp)                                      &
     &        - perlen(2)*nint(dis(2,jtmp)/perlen(2))                   
         rdyn(jtmp) = rdyn(jtmp) + dis(2,jtmp)**2 
         dis(3,jtmp) = rv2(3,i,b,x) - rv2(3,knbr(jtmp),b,x) 
         dis(3,jtmp) = dis(3,jtmp)                                      &
     &        - perlen(3)*nint(dis(3,jtmp)/perlen(3))                   
         rdyn(jtmp) = rdyn(jtmp) + dis(3,jtmp)**2 
 2110 continue 
      nneigh = 0 
      do 2100 jtmp = 1,jend 
!.....                                                                  
!     determine which pairs are separated by less than rcut             
!     and store the needed information about these pairs                
!.....                                                                  
!     if nearest periodic image is out of range, then all images        
!     will be                                                           
!.....                                                                  
         if(knbr(jtmp).eq.numnpp1) then 
            rcutsq2=INDRADSQ 
         else 
            rcutsq2=rcutsq 
         endif 
         if (rdyn(jtmp).gt.rcutsq2) go to 2100 
         j = knbr(jtmp) 
         nneigh = nneigh + 1 
         rneigh(nneigh) = rdyn(jtmp) 
         jneigh(nneigh) = j 
!--   Vijay	                                                            
         totNeighbors = totNeighbors + 1 
         NeighborList(totNeighbors, i) = j 
!--   Vijay                                                             
                                                                        
         do 2120 kcoord = 1,3 
 2120       dneigh(kcoord,nneigh) = dis(kcoord,jtmp) 
!     check periodic images in z direction                              
         if(onlyoi)go to 2100 
         disz = -sign(perlen(3),dis(3,jtmp)) 
         rim = rdyn(jtmp) + 2.*dis(3,jtmp)*disz + disz**2 
                                                                        
!     if next nearest image is out of range, subsequent ones will also b
         if (rim.gt.rcutsq2) go to 2100 
         nneigh = nneigh + 1 
         rneigh(nneigh) = rim 
         jneigh(nneigh) = j 
!--   Vijay	                                                            
         totNeighbors = totNeighbors + 1 
         NeighborList(totNeighbors, i) = j 
!--   Vijay                                                             
                                                                        
         dneigh(1,nneigh) = dis(1,jtmp) 
         dneigh(2,nneigh) = dis(2,jtmp) 
         dneigh(3,nneigh) = dis(3,jtmp) + disz 
         if(.not.threei)go to 2100 
         disz = -disz 
         rim = rdyn(jtmp) + 2.*dis(3,jtmp)*disz + disz**2 
         if (rim.gt.rcutsq2) go to 2100 
         nneigh = nneigh + 1 
         rneigh(nneigh) = rim 
         jneigh(nneigh) = j 
!--   Vijay	                                                            
         totNeighbors = totNeighbors + 1 
         NeighborList(totNeighbors, i) = j 
!--   Vijay	                                                            
         dneigh(1,nneigh) = dis(1,jtmp) 
         dneigh(2,nneigh) = dis(2,jtmp) 
         dneigh(3,nneigh) = dis(3,jtmp) + disz 
         if(.not.fouri)go to 2100 
         disz = -2.*sign(perlen(3),dis(3,jtmp)) 
         rim = rdyn(jtmp) + 2.*dis(3,jtmp)*disz + disz**2 
         if (rim.gt.rcutsq2) go to 2100 
         nneigh = nneigh + 1 
         rneigh(nneigh) = rim 
         jneigh(nneigh) = j 
!--   Vijay	                                                            
         totNeighbors = totNeighbors + 1 
         NeighborList(totNeighbors, i) = j 
!--   Vijay	                                                            
         dneigh(1,nneigh) = dis(1,jtmp) 
         dneigh(2,nneigh) = dis(2,jtmp) 
         dneigh(3,nneigh) = dis(3,jtmp) + disz 
 2100 continue 
!.....                                                                  
!.....now do diagonal (i=j) term for three or four images of self       
!.....both cases produce two images                                     
!.....                                                                  
      nneips = nneigh 
      if(threei.and.i.ne.numnpp1)then 
!.....first image of self                                               
         nneips = nneips + 1 
         disz = perlen(3) 
         rim = disz**2 
!     don't need to check against range here                            
         rneigh(nneips) = rim 
         jneigh(nneips) = i 
!--   Vijay	                                                            
         totNeighbors = totNeighbors + 1 
         NeighborList(totNeighbors, i) = i 
!--   Vijay                                                             
                                                                        
         dneigh(1,nneips) = 0.0 
         dneigh(2,nneips) = 0.0 
         dneigh(3,nneips) = disz 
!.....second image of self                                              
         nneips = nneips + 1 
         disz = -disz 
!.....rim = disz**2    done above                                       
         rneigh(nneips) = rim 
         jneigh(nneips) = i 
!--   Vijay	                                                            
         totNeighbors = totNeighbors + 1 
         NeighborList(totNeighbors, i) = i 
!--   Vijay	                                                            
         dneigh(1,nneips) = 0.0 
         dneigh(2,nneips) = 0.0 
         dneigh(3,nneips) = disz 
      endif 
                                                                        
      if(nneips.gt.NEIMAX)then 
         write(*,9212)nneips,NEIMAX 
 9212    format(' number of neighbors',i5,                              &
     &        ' exceeds array bound ',i5)                               
         stop 
      endif 
!.....                                                                  
!.....compute the maximum number of neighbors seen so far               
!.....                                                                  
      nneimx = max0(nneimx,nneips) 
!.....                                                                  
!.....compute the timing                                                
!.....                                                                  
!.....call timeused(icpu,io,isys,imem)                                  
!.....gnetim = (1.e-6)*float(icpu) - timin                              
!.....gnetmx = amax1(gnetmx,gnetim)                                     
!.....gnetmn = amin1(gnetmn,gnetim)                                     
                                                                        
!--   Vijay                                                             
!     print*, 'Atom: ', i, ' Neighbors: ', NumNeighbors(i)              
      NumNeighbors(i) = totNeighbors 
!--   Vijay		                                                           
      return 
!.....                                                                  
!.....end of neighbor finding when using old neighbor list              
!.....                                                                  
!-------------------------------------------------------------          
 2500 continue 
!.....                                                                  
!.....determine new neighbor list while getting the neighbors           
!.....                                                                  
!.....                                                                  
!     c-JS            nnindx(i) = nnindx(i-1)                           
      nnindx(i)=0 
                                                                        
      nneigh = 0 
      nneips = 0 
      if(.not.NeedList) go to 2800 
      do 2700 j = 1,numnp 
         if(i.eq.j.or.IsRelaxed(j).eq.0) go to 2700 
!.....                                                                  
!     compute the square of the distance to the closest periodic image  
!.....                                                                  
         dis(1,j) = rv2(1,i,b,x) - rv2(1,j,b,x) 
         dis(1,j) = dis(1,j) - perlen(1)*nint(dis(1,j)/perlen(1)) 
         rdyn(j) = dis(1,j)**2 
         if (rdyn(j).gt.rctsqn) goto 2700 
         dis(2,j) = rv2(2,i,b,x) - rv2(2,j,b,x) 
         dis(2,j) = dis(2,j) - perlen(2)*nint(dis(2,j)/perlen(2)) 
         rdyn(j) = rdyn(j) + dis(2,j)**2 
         if (rdyn(j).gt.rctsqn) goto 2700 
         dis(3,j) = rv2(3,i,b,x) - rv2(3,j,b,x) 
         dis(3,j) = dis(3,j) - perlen(3)*nint(dis(3,j)/perlen(3)) 
         rdyn(j) = rdyn(j) + dis(3,j)**2 
!.....                                                                  
!.....determine if these particles are within the storage distance      
!.....                                                                  
         if (rdyn(j).gt.rctsqn) goto 2700 
!.....                                                                  
!.....store the index of the particle                                   
!.....                                                                  
!     c--JS                                                             
         nnindx(i) = nnindx(i) + 1 
         nnlst(i,nnindx(i)) = j 
!.....                                                                  
!     determine which pairs are separated by less than rcut             
!     and store the needed information about these pairs                
!.....                                                                  
!     if nearest periodic image is out of range, then all images        
!     will be                                                           
!.....                                                                  
         if (rdyn(j).gt.rcutsq) go to 2700 
         nneigh = nneigh + 1 
         rneigh(nneigh) = rdyn(j) 
         jneigh(nneigh) = j 
!--   Vijay	                                                            
         totNeighbors = totNeighbors + 1 
         NeighborList(totNeighbors, i) = jtmp 
!--   Vijay	                                                            
         do 2750 kcoord = 1,3 
 2750       dneigh(kcoord,nneigh) = dis(kcoord,j) 
!     check periodic images in z direction                              
         if(onlyoi)go to 2700 
         disz = -sign(perlen(3),dis(3,j)) 
         rim = rdyn(j) + 2.*dis(3,j)*disz + disz**2 
!     if next nearest image is out of range, subsequent ones will also b
         if (rim.gt.rcutsq) go to 2700 
         nneigh = nneigh + 1 
         rneigh(nneigh) = rim 
         jneigh(nneigh) = j 
!--   Vijay	                                                            
         totNeighbors = totNeighbors + 1 
         NeighborList(totNeighbors, i) = jtmp 
!--   Vijay		                                                           
         dneigh(1,nneigh) = dis(1,j) 
         dneigh(2,nneigh) = dis(2,j) 
         dneigh(3,nneigh) = dis(3,j) + disz 
         if(.not.threei)go to 2700 
         disz = -disz 
         rim = rdyn(j) + 2.*dis(3,j)*disz + disz**2 
         if (rim.gt.rcutsq) go to 2700 
         nneigh = nneigh + 1 
         rneigh(nneigh) = rim 
         jneigh(nneigh) = j 
!--   Vijay	                                                            
         totNeighbors = totNeighbors + 1 
         NeighborList(totNeighbors, i) = jtmp 
!--   Vijay		                                                           
         dneigh(1,nneigh) = dis(1,j) 
         dneigh(2,nneigh) = dis(2,j) 
         dneigh(3,nneigh) = dis(3,j) + disz 
         if(.not.fouri)go to 2700 
         disz = -2.*sign(perlen(3),dis(3,j)) 
         rim = rdyn(j) + 2.*dis(3,j)*disz + disz**2 
         if (rim.gt.rcutsq) go to 2700 
         nneigh = nneigh + 1 
         rneigh(nneigh) = rim 
         jneigh(nneigh) = j 
!--   Vijay	                                                            
         totNeighbors = totNeighbors + 1 
         NeighborList(totNeighbors, i) = jtmp 
!--   Vijay		                                                           
         dneigh(1,nneigh) = dis(1,j) 
         dneigh(2,nneigh) = dis(2,j) 
         dneigh(3,nneigh) = dis(3,j) + disz 
 2700 continue 
!                                                                       
!     add indenter                                                      
!                                                                       
      if(numnpp1.gt.numnp.and.i.ne.numnpp1) then 
         nnindx(i) = nnindx(i) + 1 
         nnlst(i,nnindx(i)) = numnpp1 
         nneigh = nneigh + 1 
         jneigh(nneigh) = numnpp1 
         rneigh(nneigh)=0.d0 
         do kcoord=1,3 
            dneigh(kcoord,nneigh) = rv2(kcoord,i,b,x) -                 &
     &           rv2(kcoord,numnpp1,b,x)                                
            rneigh(nneigh)=rneigh(nneigh)+dneigh(kcoord,nneigh)         &
     &           **2                                                    
         enddo 
         if(rneigh(nneigh).gt.INDRADSQ) nneigh=nneigh-1 
      endif 
!.....                                                                  
!.....now do diagonal (i=j) term for three or four images of self       
!.....both cases produce two images                                     
!.....                                                                  
      nneips = nneigh 
      if(threei.and.i.ne.numnpp1)then 
!.....first image of self                                               
         nneips = nneips + 1 
         disz = perlen(3) 
         rim = disz**2 
!     don't need to check against range here                            
         rneigh(nneips) = rim 
         jneigh(nneips) = i 
!--   Vijay	                                                            
         totNeighbors = totNeighbors + 1 
         NeighborList(totNeighbors, i) = jtmp 
!--   Vijay		                                                           
         dneigh(1,nneips) = 0.0 
         dneigh(2,nneips) = 0.0 
         dneigh(3,nneips) = disz 
!.....second image of self                                              
         nneips = nneips + 1 
         disz = -disz 
!.....rim = disz**2    done above                                       
         rneigh(nneips) = rim 
         jneigh(nneips) = i 
!--   Vijay	                                                            
         totNeighbors = totNeighbors + 1 
         NeighborList(totNeighbors, i) = jtmp 
!--   Vijay		                                                           
         dneigh(1,nneips) = 0.0 
         dneigh(2,nneips) = 0.0 
         dneigh(3,nneips) = disz 
      endif 
      if(nneips.gt.NEIMAX)then 
         write(*,9222)nneips,NEIMAX 
 9222    format(' number of neighbors',i5,                              &
     &        ' exceeds array bound ',i5)                               
         stop 
      endif 
!.....                                                                  
!.....compute the maximum number of neighbors seen so far               
!.....                                                                  
      nneimx = max0(nneimx,nneips) 
      mxnnei = max0(mxnnei,nnindx(i)) 
!     c-JS               mxlstu = max0(mxlstu,nnindx(i))                
                                                                        
 2800 continue 
!.....                                                                  
!.....if this is the last call to gneigh, then set newlst=0 to indicate 
!.....that the current list can be used and also increment ngtlst       
!.....                                                                  
      if ((i.eq.numnpp1).or.(numnpp1.lt.0.and.i.eq.numnp)) then 
         newlst = 0 
         ngtlst = ngtlst + 1 
         write(*,*) 'UPDATING NEIGHBORS' 
      end if 
!.....                                                                  
!.....compute the timing                                                
!.....                                                                  
!.....call timeused(icpu,io,isys,imem)                                  
!.....gnetim = (1.e-6)*float(icpu) - timin                              
!.....gnetmx = amax1(gnetmx,gnetim)                                     
!.....gnetmn = amin1(gnetmn,gnetim)                                     
                                                                        
!--   Vijay                                                             
      NumNeighbors(i) = totNeighbors 
!     print*, 'Atom: ', i, ' Neighbors: ', NumNeighbors(i)              
!--   Vijay                                                             
      return 
!.....                                                                  
!.....end of nmeth = 2                                                  
!.....                                                                  
      end subroutine gneigh 
!*********************************************************************  
!                                                                       
!     this subroutine computes the displacement of each particle since  
!     the last update of the neighbor list.                             
!     if the maximum displacement is more than 1/2 of dradn, then newlst
!     is set to flag the creation of a new neighbor list by gneigh      
!                                                                       
      subroutine chkdis(b,x) 
      use mod_global 
      implicit double precision (a-h,o-z) 
      dimension lstold(NEIMAX) 
      dimension b(ndf,*),x(nxdm,*) 
      logical firsttime 
      data firsttime /.true./ 
!     C--Jun Song: User defined distance-test whether atom moves too muc
      double precision DefDist(numnp) 
      integer nloop 
                                                                        
      DefDist=0.0d0 
!     C--JS: Reset UpdateNeigh(:)                                       
      UpdateNeigh=0 
!                                                                       
!     branch to appropriate neighbor method                             
!                                                                       
!     Ron's Changes                                                     
      if(firsttime) then 
         newlst=1 
         firsttime=.false. 
      else			 
         newlst = 0 
                                                                        
!     c--		return                                                       
      endif	 
                                                                        
      goto (1000,2000,3000,2000) nmeth 
 1000 return 
!                                                                       
!     nmeth = 2 by default                                              
!                                                                       
 2000 continue 
!                                                                       
!     treat the second call the same as the first in case defects have b
!     added                                                             
!                                                                       
!     also if newlst has been set to 1, do not override even if the part
!     have not moved                                                    
!                                                                       
      if (newlst.eq.1) goto 2500 
!                                                                       
!     compare the new positions with the old ones                       
!                                                                       
      do 2100 i = 1,numnp 
         if(IsRelaxed(i).ne.0) then 
            dis(1,i) = rold(1,i) - rv2(1,i,b,x) 
            dis(1,i) = dis(1,i) - perlen(1)*nint(dis(1,i)/perlen(1)) 
            dis(2,i) = rold(2,i) - rv2(2,i,b,x) 
            dis(2,i) = dis(2,i) - perlen(2)*nint(dis(2,i)/perlen(2)) 
            dis(3,i) = rold(3,i) - rv2(3,i,b,x) 
            dis(3,i) = dis(3,i) - perlen(3)*nint(dis(3,i)/perlen(3)) 
            rdyn(i) = dis(1,i)**2 + dis(2,i)**2 + dis(3,i)**2 
         endif 
 2100 continue 
!                                                                       
!     determine the maximum displacement and compare with  dradn        
!                                                                       
      drmax1 = 0.0 
      drmax2 = 0.0 
                                                                        
      do 2200 i = 1,numnp 
         if(IsRelaxed(i).ne.0) then 
            tmp = min(drmax1,rdyn(i)) 
            drmax1 = max(drmax1,rdyn(i)) 
            drmax2 = max(drmax2,tmp) 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
!     C--JS: set UpdateNeigh for i atom if exceed max allowable distance
!     c--this is to find the exact atom that leads to neighbr updating  
!     c--Note that here large neighbor list is used based on rctsq      
            DefDist(i)=2.0d0*dsqrt(rdyn(i)) 
            if(DefDist(i) .gt. dradn) then 
               UpdateNeigh(i)=1 
!     c--also set Update flag for all i's neighbors                     
               nloop=nnindx(i) 
               do j = 1,nnindx(i) 
                  UpdateNeigh(nnlst(i,j))=1 
                       enddo 
                    endif 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
                 endif 
 2200         continue 
              drmax = dsqrt(drmax1) + dsqrt(drmax2) 
              if (drmax.gt.dradn) goto 2500 
!                                                                       
!     if here the old neighbor list can be used                         
!                                                                       
              newlst = 0 
              return 
!-----------------                                                      
 2500         continue 
!                                                                       
!     if here, a new neighbor list is needed so store the current coordi
!                                                                       
                                                                        
              newlst = 1 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
!     C--Jun Song, only update rold for specific atoms                  
              if(SimStep .eq. 1) then 
                 do 2600 j = 1,3 
                    do 2600 i = 1,numnp 
 2600                  rold(j,i) = rv2(j,i,b,x) 
                 endif 
                                                                        
                 if(SimStep .gt. 1) then 
                    do i=1, numnp 
                       if(UpdateNeigh(i) .eq. 1) then 
                          do j=1, 3 
                             rold(j,i)=rv2(j,i,b,x) 
                          enddo 
                       endif 
                    enddo 
                 endif 
                                                                        
                 return 
!                                                                       
!     nmeth = 3                                                         
!                                                                       
 3000            stop 'method three not here.' 
                 end subroutine chkdis 
                                                                        
                                                                        
!     Get the neighborlist                                              
                                                                        
                                                                        
                                                                        
                                                                        
!.....&Modified: Wed Jan 18 10:10:47 1995 by foiles &                   
!***************************************************************        
      subroutine chkper 
      implicit double precision (a-h,o-z) 
      data ncalls/0/ 
      save nx0,ny0,nz0 
!                                                                       
!     check periodicity:                                                
!     allow only one periodic image in the x and y directions           
!     allow up to four in the z direction                               
!                                                                       
!     x and y:                                                          
!     if perlen is lt 2*dsqrt(rcutsq), then have troubles with periodici
!                                                                       
      rcutsq=CutoffR2(1) 
      rcut = dsqrt(rcutsq) 
      permin = 2.*rcut 
      istop = 0 
      do 100 i=1,2 
         if(perlen(i).ge.permin)go to 100 
         write(*,9230)i 
 9230    format('   periodicity is too short in the ',i2,'  direction') 
         write(*,9240)permin,perlen(i) 
 9240    format('   permin = ',e15.5,'  periodicity = ',e15.5) 
         istop = 1 
  100 continue 
      if(istop.eq.1)stop 
!                                                                       
!     z direction:                                                      
!     four levels                                                       
!                                                                       
!.....initial guess                                                     
      onlyoi=.true. 
      twoi=.false. 
      threei=.false. 
      fouri=.false. 
!                                                                       
      if(perlen(3).lt.2.*rcut)then 
         onlyoi=.false. 
         twoi=.true. 
      endif 
      if(onlyoi)go to 500 
      if(perlen(3).lt.rcut)threei=.true. 
      if(perlen(3).lt.2.*rcut/3.)fouri=.true. 
      if(perlen(3).lt.0.5*rcut)then 
         i=3 
         permin = 0.5*rcut 
         write(*,9230)i 
         write(*,9240)permin,perlen(i) 
         stop 
      endif 
  500 continue 
      nx = 1 
      ny = 1 
      nz = 1 
      if(fouri)nz = 4 
      if(threei.and..not.fouri)nz = 3 
      if(twoi.and..not.threei)nz = 2 
      if(onlyoi)nz = 1 
!                                                                       
!     print out the number of images used on the first call only        
!                                                                       
      if (ncalls.eq.0)then 
         write(*,9035)nx,ny,nz 
 9035    format(' # of periodic images   ',i10,5x,i10,5x,i10) 
         nx0 = nx 
         ny0 = ny 
         nz0 = nz 
      endif 
      ncalls = 1 
!                                                                       
!     print out notice if # of periodic images changes                  
      if(nx.ne.nx0.or.ny.ne.ny0.or.nz.ne.nz0)then 
         nx0 = nx 
         ny0 = ny 
         nz0 = nz 
         write(*,9040) 
 9040    format(//' ****** changing # of periodic images ') 
         write(*,9035)nx,ny,nz 
      endif 
      return 
      end subroutine chkper 
                                                                        
      END                                           
