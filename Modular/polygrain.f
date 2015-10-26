c**------------------------------------------------------------------
c** NearestBSite : get nearest Bravais site to the point in the model
c**
c**    Parmeters :
c**            u   (in) : Coordinates of the point
c**        nnear   (in) : Number of nearest sites wanted
c**   CheckModel   (in) : Flag, true if point must be in model
c**           xa  (out) : coordinates of the B-site
c**       igrain  (out) : Grain in which these sites are
c**
c**    Algorithm :
c**          for each grain
c**              find that B-site that is closest
c**              is it the smallest
c**          end do
c**
c--
      subroutine NearestBSite(u,nnear,CheckModel,xa,igrain)

      use mod_grain
      use mod_global
      implicit none

c--Variables Transferred
      double precision u(ndm), xa(3,1)
C90 ==> This variable should be removed and dynamic allocation used
      integer NNEAMAX
      parameter(NNEAMAX = 9)
      integer   igrain(1), nnear
      Logical CheckModel

c--Local Variables
      common/debugger/debug
      logical debug

      integer i, j, k,  k1, kt, ispt1
      double precision x(3,NNEAMAX), rmin(NNEAMAX), r(NNEAMAX),
     &                 rtmp, xtmp(3)

      logical change
c

      if(nnear .gt. NNEAMAX) then
         print*, 'NNEAMAX smaller than requested number neighs.'
         print*, 'Increase NNEAMAX in subroutine NearestBSite'
         print*, 'NNEAMAX =', NNEAMAX
         stop
      end if

      !Start with the first grain
      if(nnear.eq.1) then
         call NearBSiteInGrain1(u,1,CheckModel,x(1,1),r(1))
      else
         call NearBSiteInGrain(u,nnear,1,CheckModel,x(1,1),r(1))
      endif
      do k  = 1, nnear
         igrain(k) = 1
         rmin(k)  = r(k)
         do j = 1, 3
            xa(j,k) = x(j,k)
         end do
      end do

      !now for rest of the grains
      do i = 2, ngrains
         change = .false.
         if(nnear.eq.1) then
            call NearBSiteInGrain1(u,i,CheckModel,x(1,1),r(1))
         else
            call NearBSiteInGrain(u,nnear,i,CheckModel,x(1,1),r(1))
         endif
         do k = 1, nnear
            rtmp = rmin(k)
            do k1 = 1, nnear
               if (r(k1) .lt. rtmp) then
                  change = .true.
                  kt = k1
                  rtmp = r(k1)
               end if
            end do

            if(change) then
               do j = 1, 3
                  xtmp(j) = xa(j,k)   !flip atoms
                  xa(j,k) = x(j,kt)
                  x(j,kt) = xtmp(j)
               end do
               r(kt) = rmin(k)        !flip distances
               rmin(k) = rtmp
               igrain(k) = i          !set grain
            end if
         end do
      end do
      return
      end

c**------------------------------------------------------------------
c** NearBSiteInGrain : Get Bravais site nearest to given point in the grain
c**
c**      Parameters :
c**                u  (in) : coordinates of the point
c**            nnear  (in) : number of neartest B-sites requested
c**           igrain  (in) : the grain that is required
c**       CheckModel  (in) : Flag true if point must lie in the model
c**               xa (out) : coordinates of the B-sites
c**               rd (out) : distance square of near B-site
c**
c**      Algorithm :
c**           check if point is in the grain
c**           transform point to grain coordinates.
c**           ff (exterior point) then
c**              Find the nearest interior point
c**           end if
c**           find Cell and its eight neighbours
c**           for all B-sites in the cell,
c**              if (not interior) discard
c**              find the smallest
c**           end do
c**
c**      Notes :
c**           Does not change polygrain.
c**           This is a very messy sub.
c**
c--
      subroutine NearBSiteInGrain(u,nnear,igrain,CheckModel,xa,rd)

      use mod_grain
      use mod_global
      implicit none

      integer NNEAMAX
      parameter(NNEAMAX = 9)

c--Variables Transferred
      double precision u(ndm), xa(3,NNEAMAX), rd(NNEAMAX)
      integer   igrain, nnear
      logical CheckModel

c--Local Varibles

      logical PointInGrain, foundit, PointInPoly,
     &     pim
      double precision  xx, yy, une(ndm), tol, cx, cy, ul(ndm),
     &   ymin, zmin, uat(3), umin(3,NNEAMAX), xcorg, ycorg
      double precision,pointer :: cell(:,:)
      double precision  rmin(NNEAMAX), r, xc, yc,
     &                  rtmp, utmp(3), xtmp(3), xat(3)
      integer i, j, n, m, k, k1, ncell, nvt, icell,idummy
      parameter(tol=1.0d-6)

      if (ndm.eq.3) stop 'NearBSiteInGrain not ready for 3D'
      allocate(cell(3,grains(igrain)%ncell))
      do i = 1, ndm
         une(i)  = u(i)
      end do

      if( .not. PointInGrain(u,igrain)) then
         call FindNearestPoint(u, grains(igrain)%numvrts,
     &        grains(igrain)%grainarr,une,idummy)
      end if


      !load required data to local arrays
      nvt = grains(igrain)%numvrts
      ncell = grains(igrain)%ncell
      cx = grains(igrain)%dcell(1)
      cy = grains(igrain)%dcell(2)
      do i = 1, ncell
         do j = 1, 3
            cell(j,i) = grains(igrain)%cell(j,i)
         end do
      end do


      !find corrdinates of the point in local system
      ul(1) = grains(igrain)%rotmat(1)*( u(1) -
     &        grains(igrain)%refatom(1)) +
     &        grains(igrain)%rotmat(2)*( u(2) -
     &        grains(igrain)%refatom(2))
      ul(2) = -grains(igrain)%rotmat(2)*( u(1) -
     &        grains(igrain)%refatom(1)) +
     &        grains(igrain)%rotmat(1)*( u(2) -
     &        grains(igrain)%refatom(2))


      !shift origin to reference atom
      do i = 1, ndm
         une(i)  = une(i) - grains(igrain)%refatom(i)
      end do

      !get into grain coordinates
      xx =  une(1)*grains(igrain)%rotmat(1) +
     &       une(2)*grains(igrain)%rotmat(2)
      yy = -une(1)*grains(igrain)%rotmat(2) +
     &       une(2)*grains(igrain)%rotmat(1)

      !find which repeating cell
      m = int(xx*(1.0+tol)/cx)
      n = int(yy*(1.0+tol)/cy)
      if( xx.lt.tol .and. dabs(dble(m)*cx-xx).gt.tol) m = m - 1
      if( yy.lt.-tol .and. dabs(dble(n)*cy-yy).gt.tol)n = n - 1
      xcorg = dble(m)*cx
      ycorg = dble(n)*cy


      do i = 1, nnear
         rmin(i) = 1.d16
      end do
      do i = -2, 2, 1
         do j = -2 , 2, 1
            xc = xcorg + dble(i)*cx
            yc = ycorg + dble(j)*cy
            do icell = 1, ncell
               uat(1) = xc + cell(1,icell) + grains(igrain)%rfudgvec(1)
               uat(2) = yc + cell(2,icell) + grains(igrain)%rfudgvec(2)
               uat(3) = cell(3,icell)
               if(PointInPoly(uat,nvt,grains(igrain)%rotgrain(1,1)))
     $              then
                  uat(1) = uat(1) -  grains(igrain)%rfudgvec(1)
                  uat(2) = uat(2) -  grains(igrain)%rfudgvec(2)
                  xat(1) = grains(igrain)%rotmat(1)*uat(1)
     @                 -grains(igrain)%rotmat(2)*uat(2)
                  xat(2) = grains(igrain)%rotmat(2)*uat(1)
     @                 +grains(igrain)%rotmat(1)*uat(2)
                  xat(1) = xat(1) + grains(igrain)%refatom(1)
                  xat(2) = xat(2) + grains(igrain)%refatom(2)
                  xat(3)=uat(3)
                  r = (uat(1)-ul(1))**2 + (uat(2)-ul(2))**2
                  do k = 1, nnear
                     pim=.true.
                     if((r.le.rmin(k)).and.pim) then
                        rtmp  = rmin(k)
                        rmin(k) = r
                        r = rtmp
                        do k1 = 1, 3
                           utmp(k1) = umin(k1,k)
                           umin(k1,k) = uat(k1)
                           uat(k1) = utmp(k1)
                           xtmp(k1) = xa(k1,k)
                           xa(k1,k) = xat(k1)
                           xat(k1) = xtmp(k1)
                        end do
                     end if
                  end do
               end if
            end do
         end do
      end do

      foundit = .true.
      do i = 1, nnear
         if( rmin(k) .gt. 0.9e16) then
            foundit =  .false.
         end if
      end do

      if (foundit) then
          do i = 1, nnear
             rd(i) = rmin(i)
          end do
      else
          print*,'***Error : Could not locate all near atoms '
          stop
      end if
      deallocate(cell)
      return
      end

c**------------------------------------------------------------------
c** NearBSiteInGrain : Get Bravais site nearest to given point in the grain
c**
c**      Parameters :
c**                u  (in) : coordinates of the point
c**            nnear  (in) : number of neartest B-sites requested
c**           igrain  (in) : the grain that is required
c**       CheckModel  (in) : Flag true if point must lie in the model
c**               xa (out) : coordinates of the B-sites
c**               rd (out) : distance square of near B-site
c**
c**      Algorithm :
c**           check if point is in the grain
c**           transform point to grain coordinates.
c**           ff (exterior point) then
c**              Find the nearest interior point
c**           end if
c**           find Cell and its eight neighbours
c**           for all B-sites in the cell,
c**              if (not interior) discard
c**              find the smallest
c**           end do
c**
c**      Notes :
c**           Does not change polygrain.
c**           This is a very messy sub.
c**
c--
      subroutine NearBSiteInGrain1(u,igrain,CheckModel,xa,rd)

      use mod_grain
      use mod_material
      use mod_global
      implicit none

c--Variables Transferred
      double precision u(ndm), xa(3), rd
      integer   igrain, nnear
      logical CheckModel

c--Local Varibles

      logical PointInGrain
      double precision  une(ndm), tol
      integer idummy
      double precision a(3,3),b(3,3),bvec(3,3)
      parameter(tol=1.0d-6)
      common/debugger/debug
      logical debug

      if (ndm.eq.3) stop 'NearBSiteInGrain not ready for 3D'

      a=grains(igrain)%xlatvect
      b=material(grains(igrain)%matgrain)%bvec
      bvec=matmul(a,b)
      a=0.d0
      a(3,3)=1.d0
      a(1,1)=grains(igrain)%rotmat(1)
      a(2,2)=grains(igrain)%rotmat(1)
      a(1,2)=-grains(igrain)%rotmat(2)
      a(2,1)=grains(igrain)%rotmat(2)
      bvec=matmul(a,bvec)
      if(PointInGrain(u,igrain)) then
         call nearestbravais(bvec,grains(igrain)%refatom,u,xa,rd,
     &        CheckModel,igrain)
      else
         call FindNearestPoint(u, grains(igrain)%numvrts,
     &        grains(igrain)%grainarr,une,idummy)
         call nearestbravais(bvec,grains(igrain)%refatom,une,xa,rd,
     &        CheckModel,igrain)
         rd=sqrt((u(1)-xa(1))**2+(u(2)-xa(2))**2)
      end if

      end
c***********************************************************************
      subroutine nearestbravais(bvec,org,u,x,rd,CheckModel,igrain)
      implicit none
      common/debugger/debug
      logical debug
      double precision bvec(3,3),org(3),u(2),x(3),ZERO,ONE,BIG,rd,r2
     $     ,r2best,u3(3),binv(3,3),vec(3),p(3)
      integer i,j,l(3),i1,i2,i3,igrain
      parameter(ZERO=0.d0,ONE=1.0d0,BIG=1.e6)
      logical CheckModel,PointInGrain,pim
      u3(1)=u(1)-org(1)
      u3(2)=u(2)-org(2)
      u3(3)=ZERO-org(3)
      call inv33(bvec,binv)
      do i=1,3
         vec(i)=0
         do j=1,3
            vec(i)=vec(i)+binv(i,j)*u3(j)
         enddo
         if(vec(i).lt.ZERO) then
            l(i)=int(vec(i)-ONE)
         else
            l(i)=int(vec(i))
         endif
      enddo
      r2best=BIG
      do i1=l(1),l(1)+1
         do i2=l(2),l(2)+1
            do i3=l(3),l(3)+1
               p=org
               do i=1,3
                  p(i)=p(i)+bvec(i,1)*i1+bvec(i,2)*i2+bvec(i,3)*i3
               enddo
               pim=.true.
               if(pim) pim=PointInGrain(p,igrain)
               if(pim) then
                  r2=(p(1)-u(1))**2+(p(2)-u(2))**2
                  if(r2.lt.r2best)then
                     r2best=r2
                     x=p
                  endif
               endif
            enddo
         enddo
      enddo
      rd=sqrt(r2best)
      end subroutine nearestbravais


c**------------------------------------------------------------------
c** FindNearestPoint : Find the nearest point in the polygon to the
c**                    given point.
c**
c**      Parameters :
c**                uin   (in) : the given point
c**                nvt   (in) : number of vertices in the polygon
c**           nvt_phys   (in) : physical array dimension
c**               poly   (in) : the polygon
c**                  v  (out) : the nearest point
c**               iedge (out) : edge on which the nearest point is located
c**
c**      Algorithm :
c**             Scan over all edges in the polygon,
c**             find the nearest point on each edge,
c**             take the nearest of the nearest
c**
c**      Notes :
c**             Works only for 2D.
c**
c**      Changes :
c**            VBS, Aug 24, 1996 : Returns the edge on which the nearest
c**                                point sits
c**
c--
      subroutine FindNearestPoint(uin,nvt,poly,v,iedge)

      use mod_global
      implicit none

c--Variables Transferred
      integer nvt,iedge
      double precision uin(ndm), poly(ndm,*), v(ndm)

c--Local Variables
      integer ivt, ivt1, i
      double precision x, RIntersectRatio,vtmp(3),Dist2,dbest,dtemp
c
      dbest=1.e30
      do i=1,3
         vtmp(i)=0.d0
      enddo
c
      do ivt=1,nvt
         ivt1=ivt+1
         if (ivt1.gt.nvt) ivt1=1
c
c     intersection with edge(ivt)
c
         x = RIntersectRatio(poly(1,ivt), poly(1,ivt1), uin)
         if (x.lt.0.) then
            vtmp(1) = poly(1,ivt)
            vtmp(2) = poly(2,ivt)
         else if (x.gt.1.) then
            vtmp(1) = poly(1,ivt1)
            vtmp(2) = poly(2,ivt1)
         else
            vtmp(1) = poly(1,ivt) + x*(poly(1,ivt1)-poly(1,ivt))
            vtmp(2) = poly(2,ivt) + x*(poly(2,ivt1)-poly(2,ivt))
         endif
         dtemp=Dist2(vtmp,uin,2)
         if(dtemp.lt.dbest) then
            iedge = ivt
            dbest=dtemp
            v(1)=vtmp(1)
            v(2)=vtmp(2)
         endif
      enddo
      return
      end


c**------------------------------------------------------------------
c** RIntersectRatio : Given two points and a third point, return the
c**                  intersection ratio of the perpendicular
c**                  line form the third point to the line fromed
c**                  by the first two points
c**
c**     Parameters :
c**                x1 (in) : coordinates of the first point
c**                x2 (in) : coordinates of the second point
c**                 a (in) : coordinates of the third point
c**
c**
c**     Algorithm :
c**                  o x2
c**                  |                u  = x2 - x1
c**                  |                u1 = a  - x1
c**                  |
c**                  |                                  u1 . u
c**                 -+-------o a   RInterSectRatio =  ---------
c**                  |      /                           |u|^2
c**                  |     /
c**                  |    /
c**                  |   /
c**                  |  /
c**                  | /
c**                  |/
c**                  o
c**                  x1
c**
c**
c--
      double precision function RIntersectRatio(x1,x2,a)

      use mod_global
      implicit none

c--Variables Transferred
      double precision x1(ndm), x2(ndm), a(ndm)

c--Local Variables
      double precision u(ndm), u1(ndm), d, dot
      integer  i

      do i = 1, ndm
         u(i)  = x2(i) - x1(i)
         u1(i) =  a(i) - x1(i)
      end do

      d   = 0.0
      dot = 0.0

      do i = 1, ndm
         d   = d + u(i)*u(i)
         dot = dot + u(i)*u1(i)
      end do

      RIntersectRatio = dot/d

      return
      end

c**----------------------------------------------------------------
c** PointInGrain : Checks if given point is in given grain
c**
c**    Parameters :
c**               u (in)  : coordinates of the point
c**          igrain (in)  : grain number
c**
c**    Algorithm :
c**              Check using fudge vector
c**
c--
      logical function PointInGrain(u,igrain)
      use mod_grain
      use mod_global
      implicit none

c--Variables transferred
      double precision u(ndm)
      integer igrain

c--Local Variables
      logical PointInPoly
      double precision utmp(ndm)
      integer i

      do i = 1, ndm
         utmp(i) = u(i) + grains(igrain)%fudgevec(i)
      end do
      PointInGrain = PointInPoly(utmp(1:ndm), grains(igrain)%numvrts
     &     ,grains(igrain)%grainarr(1,1))

      return
      end

