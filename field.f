!     Qu modified on 02/01/2005---the displacement boundary conditions
!     are treated 
!     according to  anisotropic elasticity theory. The major
!     modification is
!     employed in subroutine kfiled_displ()
      subroutine ma01(id,x,ix,f,b,dr,db,input)
      use mod_global
      implicit  none
!     
      character input*80
      double precision b(ndf,*),x(nxdm,*),f(ndf,*),dr(*),db(*)
      integer id(ndf,*),ix(nen1,*)
!     
!     local variables
!     
      double precision btmp(ndf), lambda, mu
      integer i,j
      logical fixed

      call findAverageElasticConst(lambda, mu)

!     
!     applied displacement b.c.
!     
!JM      print *, '!!!!!!!!Calling k-field!!!!!!'
      print *, '!!!!!! NOT CALLING k-field !!!!!!'

 10   do i=1,numnp
!     fixed = .false.

!REMOVED THIS LOOP AND IF SINCE WE ARE NOT APPLYING A K-FIELD
!JM         do j = 1, ndf

!     fixed = (fixed.or.(id(j,i).eq.1))
!     enddo
!     if(fixed) then

!JM            if (id(j,i) .eq. 1) then 

!     print *, 'Calling k-field on  ', i, x(1,i), x(2,i)
!cc   --Only used when doing FEM, when b and f are obtained by scaling
!cc   --the current ones by time (that's why 1.0d0 is used below)

!JM               call kfield_displ(1.0d0, x(1, i), btmp, lambda, mu, 0.0)

!     print *, 'K_field =', i, btmp(1), btmp(2)
!cc   --End of Mod

!JM            endif

!JM         end do
         do j=1,ndf

!fix applied displacement to 0.0 for non-atom nodes
            if(id(j,i).eq.1) then
!JM               f(j,i)=btmp(j)
               f(j,i)=0.0
            endif

!JM         for atomic positions with prescribed displacements in y,
!           as chosen in mesh.f by id(i,j), apply a y displacement
!           in time via f(j,j)
            if(IsRelaxed(i).eq.1) then
               f(j,i) = 0.5
            endif

         enddo
      enddo
      end


      subroutine precrack(id,x,b,f, input)
      use mod_global
      implicit none
!     c
      integer id(ndf,*)
      double precision x(nxdm,*), b(ndf, *), f(ndf,*)
      double precision lambda, mu
!CC   --Hack parameter
      double precision yHShift
      integer i
      character*80 input
!     c
      integer lower,upper,next, idum
      double precision dflag, Kinit
      double precision btmp(ndf),temp,btmp1(ndf)
!     c
      lower = 4
      upper = next(lower,input)
      call freein(input,lower,upper,idum,Kinit,2)
      lower=upper
      upper = next(lower,input)
      call freein(input,lower,upper,idum,dflag,2)
!     c
!     c

      print *,'Time in precrack is', time, Kinit, dflag
      call findAverageElasticConst(lambda, mu)
      do i = 1, numnp
!CC   --JS Hack for H atoms
         yHShift=0.0d0
         if((atomSpecie(i).EQ.2).AND.(x(1,i).LT.0.1)) then
            if((x(2,i).GT.-0.2).AND.(x(2,i).LT.0.2)) then
               yHShift=0.1d0
            endif
            if((x(2,i).GT.1.0).AND.(x(2,i).LT.2.0)) then
               yHShift=-0.1d0
            endif
         endif
!CC   --END            
         if(dflag.eq.0.0) then
            time = Kinit
            call kfield_displ(Kinit, x(1, i), btmp, lambda, mu, yHShift)
            b(1,i)=btmp(1)
            b(2,i)=btmp(2) 
            if(id(1,i).eq.1) then
               f(1,i)=b(1,i)/time
            endif
            if(id(2,i).eq.1) then
               f(2,i)=b(2,i)/time
            endif
            if (id(1,i) .eq. 1 .or. id(2,i) .eq. 1) then 
!     print *, 'K_Field in prec', i, f(1,i), f(2,i), x(1,i), x(2,i)
            endif
         else
            call kfield_displ(Kinit, x(1, i), btmp1, lambda, mu,
     $           yHShift)
            temp=Kinit+dflag
            call kfield_displ(temp, x(1, i), btmp, lambda, mu, yHShift)
            b(1,i)=b(1,i)- btmp1(1)+ btmp(1)
            b(2,i)=b(2,i)- btmp1(2)+ btmp(2)
            if(id(1,i).eq.1) then
               f(1,i)=b(1,i)/(time)
            endif
            if(id(2,i).eq.1) then
               f(2,i)=b(2,i)/(time)
            endif
         endif

      enddo
      
!     c
      return
      end

!**********************************************************************
!**********************************************************************
!**********************************************************************

!     c**---------------------------------------------------------------
!     c**   pdelcalc : compute the pdelta curve
!     c**
      subroutine pdelcalc(f,dr,id,force,x)
      use mod_boundary
      use mod_global
      implicit none
!     c
      double precision f(ndf,*),dr(ndf,*),force,x(nxdm,*)
      integer id(ndf,*)
!     c
      integer i,j
!     c
      write(6,*) 'p-delta is not implemented'
      return
      end



      subroutine kfield_displ(K_u, xin, bout, lambda, mu, yshift)
      use mod_global
      use mod_crack
      implicit none
      double precision M_PI
      parameter (M_PI = 3.141592653589793d0)
!     C
!     C  This is valid for the hexagonal lattice only. Elastic constants
!     should be
!     C  exported from the code!
!     C
      double precision lambda
!     C      parameter (lambda = 0.5747d0) ! for 2D Hexagonal Al.
!     C	parameter (lambda = 0.3742) ! for 3D FCC Al.
!     C      parameter (lambda = 0.3784) ! for 3D FCC Al, Baskes and Daw
!     potential.
!     C      parameter (lambda = 0.7961) ! for 3D FCC Ni, Baskes and Daw
!     Potential
!     C
      double precision K_u, xin(nxdm), bout(ndf)
!     c
      double precision mu, nu, e, k, yshift
      double precision x, y, r, theta
      double precision K_I,K_II,ux1,ux2,uy1,uy2
      double precision s1x,s1y,s2x,s2y,p1x,p1y,p2x,p2y
      double precision q1x,q1y,q2x,q2y
      double precision xnu, xmu
      complex*16 s1,s2,p1,p2,q1,q2,b1,b2

      mu = (0.9581d0-lambda)/2.0d0 ! for 2D Hexagonal Al.
!      print *, 'mu = ', mu, ' lambda = ', lambda
!     C      mu = (0.7371 - lambda)/2.0d0 ! for 3D FCC Al.
!     C      mu = (0.7126 - lambda)/2.0d0 ! for 3D FCC Al, Baskes and
!     Daw potential.
!     C      mu = 0.3728 ! for 3D FCC Ni, Baskes and Daw Potential
      
      nu=lambda/2.0d0/(lambda+mu)
      e=mu*2.d0*(1.d0+nu)
      k=3.d0-4.d0*nu

!      print *, 'mu =  ', mu,' nu = ', nu, ' K_u = ',K_u
!     C      stop


      x = xin(1)-x0crack
      y = xin(2)-y0crack-yshift

      r = dsqrt(x*x+y*y)
      theta = datan2(y, x)
!     c Qu modification begins
      K_I=K_u
      K_II=0.d0
      open(unit=12,file='anisoEig.inp',status='old')
      read(12,*)s1x,s1y
      read(12,*)s2x,s2y
      read(12,*)p1x,p1y
      read(12,*)p2x,p2y
      read(12,*)q1x,q1y
      read(12,*)q2x,q2y
      close(12)
      s1=dcmplx(s1x,s1y)
      s2=dcmplx(s2x,s2y)
      p1=dcmplx(p1x,p1y)
      p2=dcmplx(p2x,p2y)
      q1=dcmplx(q1x,q1y)
      q2=dcmplx(q2x,q2y)	 
      b1=cdsqrt(dcos(theta)+s1*dsin(theta))
      b2=cdsqrt(dcos(theta)+s2*dsin(theta))
      call AnisoDispl(s1,s2,p1,p2,q1,q2,b1,b2,ux1,ux2,uy1,uy2)
      bout(1)=dsqrt(2.d0*r/M_PI)*(K_I*ux1+K_II*ux2)
      bout(2)=dsqrt(2.d0*r/M_PI)*(K_I*uy1+K_II*uy2)

!     bout(1)=0.0*1e6/(2.0*.1925*100e9)*xin(2)
!     bout(2)=0.01*xin(2) + 0.0*1e6/(2.0*.1925*100e9)*xin(1)

!     xmu=mu*100e9/1e30/1.602e-19
!     xnu=nu
!     bout(1)=K_u/xmu*dsqrt(r/(2.0*M_PI))*
!     &     dcos(theta/2.0)*(1.0-2.0*xnu+dsin(theta/2.0)*dsin(theta
!     /2.0))
!     bout(2)=K_u/xmu*dsqrt(r/(2.0*M_PI))*
!     &     dsin(theta/2.0)*(2.0-2.0*xnu-dcos(theta/2.0)*dcos(theta
!     /2.0))


!     c
!     bout(1)=K_u/2.d0/e*dsqrt(r/2.d0/M_PI)*(1.d0+nu)
!     &     *((2.d0*k-1.0d0)*dcos(theta/2.0d0)-dcos(1.5d0*theta))
!     bout(2)=K_u/2.d0/e*dsqrt(r/2.d0/M_PI)*(1.d0+nu)
!     &     *((2.d0*k+1.0d0)*dsin(theta/2.0d0)-dsin(1.5d0*theta))
!     c
!     c Qu modification ends
      return
      end


      subroutine AnisoDispl(s1,s2,p1,p2,q1,q2,b1,b2,ux1,ux2,uy1,uy2)
      implicit none
!     c normalized displacement fields around a crack tip in anisotropic
!     material
!     c    ux=ux/(sqrt(2*r/PI)/K_u)
!     c    ux1--x-component under mod_I load
!     c    ux2--x-component under mod_II load
!     c    uy1--y-component under mod_I load
!     c    uy2--y-component under mod_II load   
!     c
      complex*16 s1,s2,p1,p2,q1,q2,b1,b2
      double precision ux1,ux2,uy1,uy2

      ux1=dreal((s1*p2*b2-s2*p1*b1)/(s1-s2))
      ux2=dreal((p2*b2-p1*b1)/(s1-s2))
      uy1=dreal((s1*q2*b2-s2*q1*b1)/(s1-s2))
      uy2=dreal((q2*b2-q1*b1)/(s1-s2))

      end 

      subroutine findAverageElasticConst(lambda, mu)
      
      use mod_material
      implicit none
      double precision lambda, nu, mu, cEl(6,6)
!     C	type(bravaismat), dimension(:), pointer :: material
      
      cEl = material(1)%cc

      lambda = 1.0/15.0 * ( cEl(1,1) + cEl(3,3) + 5.0 * cEl(1,2) + 
     $     8.0 * cEl(1,3) - 4.0 * cEl(4,4) )
      
      mu = 1.0/30.0 * (7.0 * cEl(1,1) - 5.0 * cEl(1,2) + 
     $     2.0 * cEl(3,3) + 12.0 * cEl(4,4) - 4.0 * cEl(1,3) )  
      
      return
      end
      
      
