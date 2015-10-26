c     
c     $Id: disl_field.f,v 1.1.1.1 2003/03/12 20:09:00 shastry Exp $
c     
      subroutine disl_u(r0, burgers, th_e, th_s, r, u_out)
c     
c     calculate displacement field for dislocation using elastic
c     anisotropy
c     
c     1234567890123456789012345678901234567890123456789012345678901234
c     1         2         3         4         5         6        7
c     
      implicit none
      double precision pi, tol
      parameter (pi = 3.14159265358979323844d0)
      parameter (tol=1.0d-3)
      integer i,j
      logical :: aniso=.false.
      logical , save :: calledf=.false. 
      complex*16, save :: ad(3,3), p(3) 
      complex*16 eta
      double precision r0(3), burgers(3), th_e, th_s,
     &     r(3), u_out(3), u_temp(3), factor
c     
      double precision cc(6,6)
      double precision xe, xnu, xlambda, xmu
      integer i_elas
c     
      common /elastic/ xe, xnu, xlambda, xmu, cc, i_elas
c     
      character*80 error_message
      double precision b2, dx, dy, r2, b_ort_x, b_ort_y,
     &     x, y, theta, ux, uy, dxr, dyr
c     
      if(i_elas.ne.1) then
         error_message = 'dislocations: call fe_elastic first!'
         call error_handler(error_message)
      endif
c     
      b2 = burgers(1)*burgers(1)+burgers(2)*burgers(2)
      dx = r(1)-r0(1)
      dy = r(2)-r0(2)
      r2 = dx*dx+dy*dy
      if(r2.lt.tol) then
         u_out(1:3)=0.d0
         return
      endif
c     
c     
      if(aniso.eq..false.) then
c     
         b_ort_x = -burgers(2)
         b_ort_y = burgers(1)
         x =  dx*burgers(1) + dy*burgers(2)
         y =  dx*b_ort_x + dy*b_ort_y
c     
         theta = datan2(y,x)
         if( theta.le.th_e) then
            theta = theta + (pi-th_e)
         else
            theta = theta - (pi+th_e)
         endif
c     
         ux = (theta + x*y/b2/2.0d0/(1.0d0 - xnu)/r2)
         uy = -( (1.0d0-2.0d0*xnu)/4.0d0/(1.0d0-xnu)*dlog(r2)
     &        + (x*x-y*y)/b2/4.0d0/(1.0d0 - xnu)/r2 )
c     
         u_out(1) = (ux*burgers(1)-uy*burgers(2))/2.0d0/pi
         u_out(2) = (ux*burgers(2)+uy*burgers(1))/2.0d0/pi
c     
         dxr=-cos(th_s)*dx-sin(th_s)*dy
         dyr= sin(th_s)*dx-cos(th_s)*dy
         u_out(3)=burgers(3)/2.0d0/pi*datan2(dyr, dxr)
c     
c     
      elseif(aniso) then
c     
         if(calledf.eq..false.) then
            calledf=.true.
            call readEigfile(p,ad)
         endif
c     
c     
         if(burgers(1)*cos(th_s)+burgers(2)*sin(th_s)<0.0) then
            factor=-1.0
         else
            factor=1.0
         endif
c     
         x=cos(th_s)*dx+sin(th_s)*dy
         y=-sin(th_s)*dx+cos(th_s)*dy
c     
         u_out(1:3)=0.0
         do i = 1,3
            do j = 1,3
               eta = x + p(j)*y
               u_out(i) = u_out(i) - (0.5d0/pi)*aimag(factor*ad(i,j) 
     &              *(log(-eta)) )
            enddo
         enddo
c     
         u_temp(1) = cos(th_s)*u_out(1) - sin(th_s)*u_out(2)
         u_temp(2) = sin(th_s)*u_out(1) + cos(th_s)*u_out(2)
         u_out(1) = u_temp(1)
         u_out(2) = u_temp(2)
c     
c     
      endif
c     
      return
      end
c     
c     
c     
c     
c     
c     
      subroutine readEigfile(p,ad)
      implicit none
      complex*16 p(3),ad(3,3)
      integer i,j
      double precision aa1,aa2
      open(unit=12,file='anisoDis.inp',status='old')
      do i=1,3
         read(12,*) aa1, aa2
         p(i)=dcmplx(aa1,aa2)
c     write(*,*) p(i)
      enddo
      do i=1,3
         do j=1,3
            read(12,*) aa1, aa2
            ad(i,j)=dcmplx(aa1,aa2)
c     write(*,*) ad(i,j)
         enddo
      enddo
      close(12)
      return
      end
c     
c     
c     
c     
      subroutine readEigfile2(p,sfact)
      implicit none
      complex*16 p(3),ad(3,3),sfact(3,3,3)
      integer i,j,n
      double precision aa1,aa2
      open(unit=12,file='anisoDis.inp',status='old')
      do i=1,3
         read(12,*) aa1, aa2
         p(i)=dcmplx(aa1,aa2)
c     write(*,*) p(i)
      enddo
      do i=1,3
         do j=1,3
            read(12,*) aa1, aa2
            ad(i,j)=dcmplx(aa1,aa2)
c     write(*,*) ad(i,j)
         enddo
      enddo
      do i=1,3
         do j=1,3
            do n=1,3
               read(12,*) aa1, aa2
               sfact(i,j,n)=dcmplx(aa1,aa2)
            enddo
         enddo
      enddo
      close(12)
      return
      end
c     

c     
c     
c     
      subroutine disl_s(r0, burgers, r, s_out, th_s)
c     
c     2D Stress tensor
c     
      implicit none
      double precision pi
      parameter (pi = 3.14159265358979323844d0)
      double precision r0(3), burgers(3), r(3), s_out(3)
      double precision cc(6,6)
      double precision xe, xnu, xlambda, xmu
      integer i_elas
      common /elastic/ xe, xnu, xlambda, xmu, cc, i_elas
      double precision dx,dy,dxr, dyr, r2, factor, vp, sp
      character*80 error_message
      integer i,j,n,v
      logical :: aniso=.false.
      logical , save :: calledf=.false.
      complex*16, save :: sfact(3,3,3), p(3)
      complex*16 eta, dw
      double precision th_s, s_temp(3), rcore, b2, b1, r1
c     
c     1<=> s(1,1)
c     2<=> s(2,2)
c     3<=> s(1,2)
c     
      if(i_elas.ne.1) then
         error_message = 'dislocations: call fe_elastic first!'
         call error_handler(error_message)
      endif
      
      b2 = burgers(2)*burgers(2) + burgers(1)*burgers(1)
      b1 = dsqrt(b2)

      rcore = 2.d0*b1;

      dx = r(1)-r0(1)
      dy = r(2)-r0(2)
      r2 = dx*dx+dy*dy
c     Cutoff dislocation stresses at 2*burgers magnitude
c     *************************************************
c     This is not strictly correct, at the pad interface
c     but this is ok until the template is put in
      if (r2 < rcore**2) then;
         r1 = dsqrt(r2)
         dx = dx/r1*rcore
         dy = dy/r1*rcore
         r2 = rcore**2
      end if
c     *************************************************
c     
c     
      if(aniso.eq..false.) then
c     
c     
         factor = xmu/pi/2.0d0/(1.0d0-xnu)
c     
         vp = burgers(2)*dx-burgers(1)*dy
         sp = burgers(1)*dx+burgers(2)*dy
c     
         s_out(1) = (vp - 2.0d0*dx*dy*sp/r2)/r2 * factor
         s_out(2) = (vp + 2.0d0*dx*dy*sp/r2)/r2 * factor
         s_out(3) =  sp*(dx*dx-dy*dy)/r2/r2 * factor
c     
      elseif(aniso) then
c     
         if(calledf.eq..false.) then
            calledf=.true.
            call readEigfile2(p,sfact)
         endif
c     
         if(burgers(1)*cos(th_s)+burgers(2)*sin(th_s)<0.0) then
            factor=-1.0
         else
            factor=1.0
         endif
c     
         s_out(1:3) = 0.d0
         dxr=cos(th_s)*dx+sin(th_s)*dy
         dyr=-sin(th_s)*dx+cos(th_s)*dy
c     
         do v=1,3
            if(v.eq.1) then
               i=1
               j=1
            elseif(v.eq.2) then
               i=2
               j=2
            elseif(v.eq.3) then
               i=1
               j=2
            endif 
            do n=1,3
               eta = dxr + p(n)*dyr
               s_out(v) = s_out(v) - (0.5/pi) * aimag(
     &              factor*sfact(i,j,n) / eta )
            enddo
         enddo
c     
         s_temp(1) = cos(th_s)*cos(th_s)*s_out(1) + sin(th_s)*sin(th_s)
     $        *s_out(2) - 2*sin(th_s)*cos(th_s)*s_out(3)
         s_temp(2) = sin(th_s)*sin(th_s)*s_out(1) + cos(th_s)*cos(th_s)
     $        *s_out(2) + 2*sin(th_s)*cos(th_s)*s_out(3)
         s_temp(3) = cos(th_s)*sin(th_s)*(s_out(1) - s_out(2)) +
     $        (cos(th_s)*cos(th_s) - sin(th_s)*sin(th_s))*
     $        s_out(3)
c     
         s_out(1:3) = s_temp(1:3)
c     
      endif
      return
      end
c     
c     
c     
c     
c     
c     
      subroutine disl_displ(r, u_out)
      implicit none
      include 'disl_parameters.par'
      double precision r(3), u_out(3), u11(3)
c     
      character*80 error_message
      double precision u_temp(3)
      integer i, j
c     
      if(i_disl.ne.1) then
         error_message = 'disl_displ: call disl_setup first!'
         call error_handler(error_message)
      endif
c     
      do 1 i=1,3
 1       u_out(i)=0.0d0
c     
!     print *, 'Disl displ' , ndisl
      do i=1, ndisl
c$$$         if (elem_disl(i) .eq. 0) then 
c$$$            print *, "Calculating bucket dislocation displacment"
c$$$            print *, elem_disl(i), r_disl(1:2,i)
c$$$         end if
         call disl_u(r_disl(1,i), burgers(1,i), theta_e(i),
     &        theta_s(i), r, u_temp)
         do 2 j=1,3
 2          u_out(j)=u_out(j)+u_temp(j)
      enddo
      return;
      end
c     
c     
c     
      subroutine disl_stress(r, s_out)
      implicit none
      include 'disl_parameters.par'
      double precision r(3), s_out(3)
c     
      double precision s_temp(3)
      character*80 error_message
      integer i, j
c     
      if(i_disl.ne.1) then
         error_message = 'disl_stress: call disl_setup first!'
         call error_handler(error_message)
      endif
c     
      do 1 i=1,3
 1       s_out(i)=0.0d0
c     
      do i=1, ndisl
         if (elem_disl(i) > 0) then 
            call disl_s(r_disl(1,i), burgers(1,i), r, s_temp
     $           ,theta_s(i))
            do 2 j=1,3
 2             s_out(j)=s_out(j)+s_temp(j)
         endif
      enddo
      return
      end


