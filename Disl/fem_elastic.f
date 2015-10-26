c
c $Id: fem_elastic.f,v 1.1.1.1 2003/03/12 20:09:00 shastry Exp $
c
        subroutine fe_elastic(elas_const)
c
c123456789012345678901234567890123456789012345678901234567890123456789012
c         1         2         3         4         5         6         7
c
        implicit none
        double precision elas_const(6,6)
        double precision cc(6,6)
        double precision xe, xnu, xlambda, xmu
        integer i_elas
c
        common /elastic/ xe, xnu, xlambda, xmu, cc, i_elas
        data i_elas /0/
c
        double precision xinv1, xinv2
        integer i,j
        character*80 error_message
c
        if(i_elas.eq.0) then
          i_elas = 1
        else
          error_message = 'Second call to fe_elastic'
          call error_handler(error_message)
        endif
c
        do 1 j = 1,6
        do 1 i = 1,6
 1      cc(i,j)=elas_const(i,j)
c
c  Hirth and Lothe
c
ccHex        xinv1 = cc(1,1)+cc(2,2)+cc(3,3)+cc(4,4)+cc(4,4)+cc(5,5)+cc(5,5)
ccHex     &        +cc(6,6)+cc(6,6)
ccHex        xinv2 = cc(1,1)+cc(1,2)+cc(1,3)+cc(2,1)+cc(2,2)+cc(2,3)
ccHex     &        +cc(3,1)+cc(3,2)+cc(3,3)
c
ccHex        xlambda = (2.0d0*xinv2-xinv1)/15.0d0
ccHex        xmu = (3.0d0*xinv1-xinv2)/30.0d0
c
ccHex        xe = xmu*(2.0d0+xlambda/(xlambda+xmu))
ccHex        xnu = xlambda/2.0d0/(xlambda+xmu)

cc--JS: 2D plane strain elastic coefficient
     	xnu=(cc(1,1)+3.0d0*cc(1,2)-2.0d0*cc(6,6))/
     &	(4.0d0*cc(1,1)+4.0d0*cc(1,2))
        xe=(1.0d0+xnu)*(cc(1,1)-cc(1,2)+2.0d0*cc(6,6))/2.0d0
        xmu=xe/(2.0d0*(1.0d0+xnu))
        xlambda=xe*xnu/((1.0d0+xnu)*(1-2.0d0*xnu))
c$$$	e_dd = xe
c$$$        nu_dd = xnu
c$$$        mu_dd = xmu
	write(*,*)"xnu,xe,xlambda,xmu",xnu,xe,xlambda,xmu
c
        return
        end

c
c $Log: fem_elastic.f,v $
c Revision 1.1.1.1  2003/03/12 20:09:00  shastry
c vijay-   Initial import.
c
c Revision 1.2  2001/07/12 06:36:49  shilkrot
c Elastic constants are now passed to fem_setup and further to fe_elastic
c
c Revision 1.1  2001/06/25 02:38:09  shilkrot
c Sets up the stiffness matrix and the isotropic moduli
c
