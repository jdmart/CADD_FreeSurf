c** FEAPLIB:  These are a lot of handy utilities and other routines that
c**   don't fit neatly elsewhere in the code.  Most are self
c**   explanatory.  Here is a brief table of contents, but note that the
c**   actual order of the routines is different.
c**
c** data dumping routines:
c**    putout          : dump a double precision array to std. output
c**    intout          : dump an integer array to std. output
c** vector and matrix manipulations:
c**    distance        : 2-D distance between 2 points
c**    Dist2           : n-D distance-squared between 2 points
c**    dot             : n-D dot product
c**    cross_product   : 3-D cross product
c**    norm            : n-D vector normalization
c**    promul          : computes c=c+a*b, a is stiffness, c and b are vectors
c**    eq_inte         : equate two integer vectors
c**    addb            : do x=x+b
c**    subb            : do x=x-b
c**    maxsearch       : find the maximum absolute value in a vector
c**    store           : store a vector in another vector
c** handy geometric routines:
c**    FindArea        : find area of a quadrilateral
c**    PolygonArea     : find area of an n-sided polygon
c**    PointInPoly     : check if a point is in a polygon
c**    getaspect       : get the aspect ratio of a triangle
c**    area2x          : convert area coordinates to cartesian
c**    between0        : check if the point zero is between two points
c**    areaix          : get the area of an element
c** character string manipulations:
c**    freein          : parse free-form input lines
c**    downcase        : make a string all lower case
c**    pcomp           : compare to 4-character strings
c**    next            : find the next delimiter (comma or space) in a string
c** Feap things:
c**    blank           : blanks unconstrained components of a vector
c**    simplex         : shape function calculations
c**    dot11           : dot product of unconstrained components
c** QC specific, but called in many places:
c**    MapToRef        : map a point from deformed to reference configuration
c**    MapToCurrent    : map a point from reference to deformed configuration
c**    TwoMaxSum       : find sum of two maximum displacements
c**    MustUpdateStatus: check if status must be updated due to large displ's.
c**    RINFEIGEN       : find (rcut)*(max. eigenvalue of F)
c**    BALANC          : balance a matrix
c**    ELMHES          : get hessian form of matrix
c**    HQR             : get eigenvalues of a matrix
************************************************************************
      logical function pcomp(a,b)
      implicit double precision (a-h,o-z)
c
      character*4 a,b
c
      pcomp = .false.
      if (a.eq.b) pcomp = .true.
      return
      end
************************************************************************
      function next(last,input)
      implicit double precision (a-h,o-z)
c
      character*80 input
c
      do 10 i = last+1,80
      if ((input(i:i).eq.',').or.(input(i:i).eq.' ')) then
      next = i
      return
      end if
   10 continue
      next = 81
      return
      end
************************************************************************
      subroutine freein(input,ll,lr,int,real,isw)
      implicit double precision (a-h,o-z)
      character*80 input
      character*1 range1
      character*2 range2
      character*80 string
c
      string = input((ll+1):(lr-1))
      ldif = lr - ll - 1
      if (ldif.eq.0) then
      if (isw.eq.1) int = 0
      if (isw.eq.2) real = 0.
      else
      if (ldif.lt.10) then
      encode(1,'(i1)',range1) ldif
      if (isw.eq.1) decode(ldif,'(i'//range1//')',string) int
      if (isw.eq.2) decode(ldif,'(f'//range1//'.0)',string) real
      else
      encode(2,'(i2)',range2) ldif
      if (isw.eq.1) decode(ldif,'(i'//range2//')',string) int
      if (isw.eq.2) decode(ldif,'(f'//range2//'.0)',string) real
      end if
      end if
      return
      end
c**-----------------------------------------------------------------------
c** Dist2 : return distance square of between two points
c**
c**      Parameters :
c**              x, y (in) : two points
c**                 n (in) : dimensions
c--
      double precision function Dist2(x,y,n)
      implicit none

      integer n
      double precision x(n), y(n)

      Dist2=dot_product(x-y,x-y)
      end

************************************************************************
c** PointInPoly : Checks if given point is in given polygon
c**
c**    Parameters :
c**               u (in)  : coordinates of the point
c**            numv (in)  : number of vertices
c**         grainar (in)  : grain array (only one polygon)
c**
c**    Algorithm :
c**              The one by Joseph O'Rourke.
c**              Refer his book (look up in lib).
c**
c--

      logical function PointInPoly(u, numv, grainar)
      implicit none


c--Variables transferred
      double precision u(2), grainar(2,*)
      integer numv

c--Local variables
      common/debugger/debug
      logical debug
      integer NVMAX
      parameter(NVMAX = 100)
      integer icross, i, i1, j
      double  precision  x, grain(2,NVMAX)
      logical between0, OnVertex
      double precision tole
      parameter(tole = 1.0d-6)

c--Check data
      if (numv .gt. NVMAX) then
          print*, 'ERROR ** Too many vertices in PointInPoly'
          print*, 'Increase NVMAX',NVMAX,numv
          stop
      end if


c--Set up local variables
      do i = 1, numv
         do j = 1, 2
            grain(j,i) = grainar(j,i) - u(j)
         end do
      end do


c--For each edge, check if it crosses the ray
      icross = 0
      do i = 1, numv
         i1 = mod((i+numv-2), numv) + 1
         ! if edge straddles the xaxis ...
         if( ((grain(2,i) .gt. 0.0).and.(grain(2,i1).le.0.0)).or.
     &        ((grain(2,i1) .gt. 0.0).and.(grain(2,i).le.0.0)) ) then
            ! edge straddles ray, compute intersection
            x = (grain(1,i)*grain(2,i1)-grain(1,i1)*grain(2,i))/
     &           (grain(2,i1) - grain(2,i))
            if( x .gt. 0.0) icross = icross + 1
         else
c--A point on the boundary is considered in the polygon.
            if (between0(grain(1,i),grain(1,i1))) then
               PointInPoly = .TRUE.
               return
            end if
         endif
      end do


      !u is inside for odd number of crossings
      if( mod(icross,2) .eq. 1) then
         PointInPoly = .TRUE.
         return
      else
c--If all else fails make sure that the point is not one of the vertices
         OnVertex = .false.
         do i = 1, numv
           if( dabs(grain(1,i)).lt.tole .and.
     &           dabs(grain(2,i)).lt.tole)OnVertex = .true.
         enddo
         if(OnVertex) then
            PointInPoly = .TRUE.
         else
            PointInPoly = .FALSE.
         endif
         return
      endif

      end function PointInPoly

************************************************************************
      logical function between0(x,y)
c checks if the point zero fall on the closed segment between x and y
      implicit none
      double precision x(2),y(2),area2,zero,tol
      parameter(zero=0.d0, tol=1.e-9)
      between0=.false.
c--make sure the three points are colinear.
      if(abs(x(2)*y(1)-x(1)*y(2)).gt.tol) return
c--check betweenness.
      if (x(1).ne.y(1)) then
         between0=(((x(1).ge.zero).and.(y(1).le.zero)).or.
     @        ((x(1).le.zero).and.(y(1).ge.zero)))
      else
         between0=(((x(2).ge.zero).and.(y(2).le.zero)).or.
     @        ((x(2).le.zero).and.(y(2).ge.zero)))
      endif
      end function between0
************************************************************************
      logical function between(x1,x2,y)
      implicit none
      double precision x1(2),x2(2),y(2),RIntersectRatio,r1
      r1=RIntersectRatio(x1,x2,y)
      between=(r1.ge.0.d0).and.(r1.le.1.d0)
      end function between

************************************************************************
c** get the aspect ratio of a triangle
      subroutine getaspect(xtri,aspect)
      implicit double precision (a-h,o-z)
c
      dimension xtri(2,3)
c
      equi = 36.d0/dsqrt(3.d0)
      det = xtri(1,2)*xtri(2,3) - xtri(1,3)*xtri(2,2)
     1    + xtri(1,3)*xtri(2,1) - xtri(1,1)*xtri(2,3)
     2    + xtri(1,1)*xtri(2,2) - xtri(1,2)*xtri(2,1)
      area = det/2.d0
      perim = dsqrt((xtri(1,2) - xtri(1,1))**2
     1            + (xtri(2,2) - xtri(2,1))**2)
     2      + dsqrt((xtri(1,3) - xtri(1,2))**2
     3            + (xtri(2,3) - xtri(2,2))**2)
     4      + dsqrt((xtri(1,1) - xtri(1,3))**2
     5            + (xtri(2,1) - xtri(2,3))**2)
      aspect = equi*area/(perim*perim)
      return
      end
************************************************************************
      double precision function areaix(kount,ixnew,xnew,nxdm)
      implicit double precision (a-h,o-z)
      dimension ixnew(3,1),xnew(nxdm,1),xtri(2,3)
      do i = 1,3
         do j = 1,2
            ii = ixnew(i,kount)
            xtri(j,i) = xnew(j,ii)
         enddo
      enddo
      det = xtri(1,2)*xtri(2,3) - xtri(1,3)*xtri(2,2)
     1    + xtri(1,3)*xtri(2,1) - xtri(1,1)*xtri(2,3)
     2    + xtri(1,1)*xtri(2,2) - xtri(1,2)*xtri(2,1)
      areaix = det/2.d0
      return
      end
************************************************************************
      subroutine cross_product(x1,x2,x3)
c
c     finds cross prod of x1,x2 --> x3
c
      implicit none
      double precision x1(3),x2(3),x3(3)
      x3(1)=x1(2)*x2(3) - x1(3)*x2(2)
      x3(2)=x1(3)*x2(1) - x1(1)*x2(3)
      x3(3)=x1(1)*x2(2) - x1(2)*x2(1)
      end subroutine cross_product
************************************************************************
C     Logical function intri returns true when point p is inside
C     triangle with vertices x1,x2,x3 and transforms p coordinates
C     to area coordinates in s. ontri is set to true if point p
C     is on the triangle boundary.
      logical function intri(x1,x2,x3,p,s,ontri)
      implicit double precision (a-h,o-z)
      logical ontri
      parameter (tol=1.d-9)
      dimension x1(2),x2(2),x3(2),p(2),s(3)
C     Calculate area of triangle
      a=dabs(trarea(x1,x2,x3))
      if (a.eq.0.d0) then
         write(6,*)'***ERROR: Model contains coincident nodes.'
         write(6,*)x1(1),x1(2)
         write(6,*)x2(1),x2(2)
         write(6,*)x3(1),x3(2)
         stop
      endif
C     Calculate area of 3 new triangles formed by connecting p with
C     with vertices
      a1=dabs(trarea(p,x2,x3))
      a2=dabs(trarea(x1,p,x3))
      a3=dabs(trarea(x1,x2,p))
C     If p is inside triangle we should have a1+a2+a3 = a
      tola=tol*a
      intri=(dabs(a1+a2+a3-a).le.tola)
      ontri=.false.
      if (intri) then
         ontri=(a1.le.tola.or.a2.le.tola.or.a3.le.tola)
         s(1)=a1/a
         s(2)=a2/a
         s(3)=a3/a
      endif
      return
      end function intri
************************************************************************
C     Area of a triangle with vertices at x1,x2,x3
      function trarea(x1,x2,x3)
      implicit double precision (a-h,o-z)
      dimension x1(2),x2(2),x3(2)
      trarea=0.5*( x1(1)*(x2(2)-x3(2)) - x1(2)*(x2(1)-x3(1))
     &               + x2(1)*x3(2) - x3(1)*x2(2) )
      return
      end function trarea
************************************************************************
      subroutine inv33(a,b)
      implicit double precision (a-h,o-z)
      dimension a(3,3),b(3,3)
      det = det33(a)
      b(1,1) = (a(2,2)*a(3,3) - a(2,3)*a(3,2))/det
      b(1,2) = (a(1,3)*a(3,2) - a(1,2)*a(3,3))/det
      b(1,3) = (a(1,2)*a(2,3) - a(2,2)*a(1,3))/det
      b(2,1) = (a(2,3)*a(3,1) - a(2,1)*a(3,3))/det
      b(2,2) = (a(1,1)*a(3,3) - a(1,3)*a(3,1))/det
      b(2,3) = (a(2,1)*a(1,3) - a(1,1)*a(2,3))/det
      b(3,1) = (a(2,1)*a(3,2) - a(2,2)*a(3,1))/det
      b(3,2) = (a(3,1)*a(1,2) - a(1,1)*a(3,2))/det
      b(3,3) = (a(1,1)*a(2,2) - a(1,2)*a(2,1))/det
      return
      end subroutine inv33
************************************************************************
      double precision function det33(a)
      implicit double precision (a-h,o-z)
      dimension a(3,3)
      det33 = a(1,1)*a(2,2)*a(3,3) + a(2,1)*a(3,2)*a(1,3)
     1      + a(1,2)*a(2,3)*a(3,1) - a(3,1)*a(2,2)*a(1,3)
     2      - a(1,2)*a(2,1)*a(3,3) - a(3,2)*a(2,3)*a(1,1)
      return
      end function det33
************************************************************************
      subroutine linecoef(x1,x2,abc)
      implicit none
c
c     given two points on a line x1 and x2, find the equation of the
c     line:   abc(1)*x + abc(2)*y + abc(3) = 0
c
      double precision x1(2),x2(2),abc(3),bot
      if(x1(1).eq.0.d0.and.x2(1).eq.0.d0) then
         abc(1)=1.d0
         abc(2)=0.d0
         abc(3)=0.d0
      else if (x1(2).eq.0.d0.and.x2(2).eq.0.d0) then
         abc(1)=0.d0
         abc(2)=1.d0
         abc(3)=0.d0
      else
         bot=x1(1)*x2(2)-x2(1)*x1(2)
         abc(1)=(x1(2)-x2(2))/bot
         abc(2)=(-x1(1)+x2(1))/bot
         abc(3)=1.d0
      endif
      end subroutine linecoef

