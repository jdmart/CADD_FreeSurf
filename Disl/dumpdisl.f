c$$$  subroutine dumpdisl(scale,logic,umag,b,ndf)
c$$$  use mod_dd_slip
c$$$  implicit none
c$$$  include 'disl_parameters.par'
c$$$  include 'fem_parameters.par'
c$$$  double precision xdef(2),u(3),umag,scale,bhat(3),s(3)
c$$$  integer logic,i,j,k1,k2,iel,node,ndf
c$$$  double precision xl(2,3),bl(3,3),b(ndf,*)
c$$$  logical intri,on,flag
c$$$  if(ndisl.lt.1) return
c$$$  write(logic,'('' VARIABLES = X Y UX UY UZ BX BY BZ'')')
c$$$  write(logic,*) 'ZONE'
c$$$  do i=1,ndisl
c$$$  iel=elem_disl(i)
c$$$  if(iel.gt.0) then
c$$$  do j=1,3
c$$$  node=iconn(j,iel)
c$$$  xl(1:2,j)=x0(1:2,node)
c$$$  call disl_displ(xl(1,j),u)
c$$$  bl(1:3,j)=b(1:3,imap(node))-u(1:3)
c$$$  enddo
c$$$  flag=intri(xl(1,1),xl(1,2),xl(1,3),r_disl(1,i),s,on)
c$$$  do k1=1,3
c$$$  bhat(k1)=0.
c$$$  do k2=1,3
c$$$  bhat(k1)=bhat(k1)+s(k2)*bl(k1,k2)
c$$$  enddo
c$$$  enddo
c$$$  else
c$$$  bhat=0
c$$$  endif
c$$$  call disl_displ(r_disl(1,i),u)
c$$$  u=u+bhat
c$$$  xdef(1:2)=r_disl(1:2,i)+umag*u(1:2)
c$$$  write(logic,1000) xdef,scale*u(1:3),burgers(1:3,i)
c$$$  end do
c$$$  1000   format(8e15.6)
c$$$  end

      module mod_disl_files
      type disl_files
      character(len=80):: fname
      end type disl_files
      
      end module mod_disl_files

      subroutine dumpdisl_vtk(scale, logic, umag, b, ndf)
      use mod_dd_slip
      implicit none

      include 'disl_parameters.par'
      include 'fem_parameters.par'
      double precision xdef(2),u(3),umag,scale,bhat(3),s(3)
      integer logic,i,j,k1,k2,iel,node,ndf, k
      integer islp, li, n, l
      double precision cphi, sphi, x, y, xstart, xend, ystart, yend
      double precision xl(2,3),bl(3,3),b(ndf,*)
      logical intri,on,flag
      double precision :: sd, xd(4,3)
      integer, dimension(:,:,:), allocatable :: lines
      double precision :: cosi, sini, ds, ds2, bsign
      ds = 50.0d0;
      ds2 = ds/2.d0

!! VTK files for dislocations 
      write(logic,fmt='(A)') '# vtk DataFile Version 2.0'
      write(logic,fmt='(A)') 'Dislocations in DD bending'
      write(logic,fmt='(A)') 'ASCII'
      write(logic,fmt='(A)') 'DATASET POLYDATA'
      write(logic,fmt='(A6,1x,I7,1x,A5)') 'POINTS',tot_disl*4, 'float'
      allocate(lines(2,3,tot_disl))

      lines = 0
      n = 0

      do islp = i, nslp
         do i = 1, ndis(islp)
            li = locphi(islp)
            cosi = cosphi(li)
            sini = sinphi(li)
            sd = sdis(i,islp)
            if (yendslp(islp) < 0.0d0) then 
               sd = -sd
            end if
            if (b_dd(i,islp) < 0.0d0) then 
               bsign = -1.0d0
            else
               bsign = 1.0d0
            end if
            xd = 0.0
            xd(1,1) = xslp(islp)+cosi*sd
            xd(1,2) = yslp(islp)+sini*sd
            xd(2,1) = xd(1,1) + bsign*ds*sini
            xd(2,2) = xd(1,2) - bsign*ds*cosi
            xd(3,1) = xd(1,1) - cosi*ds
            xd(3,2) = xd(1,2) - sini*ds
            xd(4,1) = xd(1,1) + cosi*ds
            xd(4,2) = xd(1,2) + sini*ds
            n = n + 1
            lines(1,1,n) = 2
            lines(1,2,n) = 4*(n-1) + 1 -1
            lines(1,3,n) = 4*(n-1) + 2 -1
            lines(2,1,n) = 2
            lines(2,2,n) = 4*(n-1) + 3 -1
            lines(2,3,n) = 4*(n-1) + 4 -1                 
            do k = 1, 4
               write(logic,fmt='(3(1X,E15.8))') (xd(k,l),l=1,3)
            end do
            write(logic,*)
         end do
      end do
      write(logic,*)
      write(logic,fmt='(A5,1X,I7,1X,I7)') 'LINES', 2*tot_disl, 6
     $     *tot_disl
      do i = 1, tot_disl
         do j = 1, 2
            do k = 1, 3
               write(logic, fmt='(1X,I7)',advance='no') lines(j,k,i)
            end do
            write(logic,*)
         end do
      end do
      deallocate(lines)


      end subroutine dumpdisl_vtk


      subroutine dumpdisl(scale,logic,umag,b,ndf)
      use mod_dd_slip
      implicit none
      include 'disl_parameters.par'
      include 'fem_parameters.par'
      double precision xdef(2),u(3),umag,scale,bhat(3),s(3)
      integer logic,i,j,k1,k2,iel,node,ndf
      integer islp, li
      double precision cphi, sphi, x, y, xstart, xend, ystart, yend
      double precision xl(2,3),bl(3,3),b(ndf,*)
      logical intri,on,flag
      if(ndisl.lt.1) return
      write(logic,'('' VARIABLES = X Y'')')

      write(logic,*) 'ZONE T= pos1'
      do islp = 1, nslp
         li = locphi(islp)
         cphi = cosphi(li)
         sphi = sinphi(li)
         xstart = xslp(islp)
         xend = xendslp(islp)
         ystart = yslp(islp)
         yend = yendslp(islp)
         if (li .eq. 1) then 
            do i = 1, ndis(islp)
               if (b_dd(i,islp) > 0) then 
                  if (yend > 0) then 
                     x = xstart + sdis(i,islp)*cphi
                     y = ystart + sdis(i,islp)*sphi
                  else
                     x = xstart - sdis(i,islp)*cphi
                     y = ystart - sdis(i,islp)*sphi
                  endif
                  write(logic,fmt='(2(1x,E15.9))') x,y
               endif
            end do
         endif
      end do

      write(logic,*) 'ZONE T= neg1'
      do islp = 1, nslp
         li = locphi(islp)
         cphi = cosphi(li)
         sphi = sinphi(li)
         xstart = xslp(islp)
         xend = xendslp(islp)
         ystart = yslp(islp)
         yend = yendslp(islp)
         if (li .eq. 1) then 
            do i = 1, ndis(islp)
               if (b_dd(i,islp) < 0) then 
                  if (yend > 0) then 
                     x = xstart + sdis(i,islp)*cphi
                     y = ystart + sdis(i,islp)*sphi
                  else
                     x = xstart - sdis(i,islp)*cphi
                     y = ystart - sdis(i,islp)*sphi
                  endif
                  write(logic,fmt='(2(1x,E15.9))') x,y
               endif
            end do
         endif
      end do


      write(logic,*) 'ZONE T= pos2'
      do islp = 1, nslp
         li = locphi(islp)
         cphi = cosphi(li)
         sphi = sinphi(li)
         xstart = xslp(islp)
         xend = xendslp(islp)
         ystart = yslp(islp)
         yend = yendslp(islp)
         if (li .eq. 3) then 
            do i = 1, ndis(islp)
               if (b_dd(i,islp) > 0) then 
                  if (yend > 0) then 
                     x = xstart + sdis(i,islp)*cphi
                     y = ystart + sdis(i,islp)*sphi
                  else
                     x = xstart - sdis(i,islp)*cphi
                     y = ystart - sdis(i,islp)*sphi
                  endif
                  write(logic,fmt='(2(1x,E15.9))') x,y
               endif
            end do
         endif
      end do


      write(logic,*) 'ZONE T= neg2'
      do islp = 1, nslp
         li = locphi(islp)
         cphi = cosphi(li)
         sphi = sinphi(li)
         xstart = xslp(islp)
         xend = xendslp(islp)
         ystart = yslp(islp)
         yend = yendslp(islp)
         if (li .eq. 2) then 
            do i = 1, ndis(islp)
               if (b_dd(i,islp) < 0) then 
                  if (yend > 0) then 
                     x = xstart + sdis(i,islp)*cphi
                     y = ystart + sdis(i,islp)*sphi
                  else
                     x = xstart - sdis(i,islp)*cphi
                     y = ystart - sdis(i,islp)*sphi
                  endif
                  write(logic,fmt='(2(1x,E15.9))') x,y
               endif
            end do
         endif
      end do


      write(logic,*) 'ZONE T= pos3'
      do islp = 1, nslp
         li = locphi(islp)
         cphi = cosphi(li)
         sphi = sinphi(li)
         xstart = xslp(islp)
         xend = xendslp(islp)
         ystart = yslp(islp)
         yend = yendslp(islp)
         if (li .eq. 3) then 
            do i = 1, ndis(islp)
               if (b_dd(i,islp) > 0) then 
                  x = xstart + sdis(i,islp)*cphi
                  y = ystart + sdis(i,islp)*sphi
               endif
               write(logic,fmt='(2(1x,E15.9))') x,y
            end do
         endif
      end do

      write(logic,*) 'ZONE T= neg3'
      do islp = 1, nslp
         li = locphi(islp)
         cphi = cosphi(li)
         sphi = sinphi(li)
         xstart = xslp(islp)
         xend = xendslp(islp)
         ystart = yslp(islp)
         yend = yendslp(islp)
         if (li .eq. 3) then 
            do i = 1, ndis(islp)
               if (b_dd(i,islp) < 0) then 
                  x = xstart + sdis(i,islp)*cphi
                  y = ystart + sdis(i,islp)*sphi
               endif
               write(logic,fmt='(2(1x,E15.9))') x,y
            end do
         endif
      end do
      end subroutine dumpdisl
