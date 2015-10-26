c**---------------------------------------------------------------
c**  StrainEnergy : computes total strain energy in the mesh
c**
c**   Non-Obvious Parameters : NONE
c**
c**   Algorithm :-
c**        Loop over the repatoms and add appropriate nodal energy
c**        depending on local or nonlocal status.
c--
      double precision function StrainEnergy()
      use mod_global
      implicit none
      integer i
      double precision sed
      sed = 0.0
      do i = 1, numnp
         if(IsRelaxed(i).ne.0) then
c         if(IsRelaxed(i).ge.1) then
            sed = sed + energy(i)
         endif
      end do
      if(numnpp1.gt.numnp) sed=sed+energy(numnpp1)
      StrainEnergy = sed
      return
      end
************************************************************************
