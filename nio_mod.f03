      module nio_mod
!
!     This module supports the program nio.
!
!     -H. P. Hratchian, 2022.
!
!
!     USE Connections
!
      use mqc_general
      use mqc_molecule
      use mqc_gaussian
      use mqc_algebra2
      use mqc_algebra
      use iso_fortran_env
!
!     Variable Declarations
!
      implicit none
      integer,parameter::IOut=6
!
!
!     Module Procedures
!
      CONTAINS

!
!PROCEDURE determinantOverlap
      subroutine determinantOverlap(SMatrix,SMatrixMinusHalf,  &
        diffDensityEVals,diffDensityEVecs,CMOsKet,NOccKet,nBasis,  &
        detOverlap,nPlusOne,nMinusOne,iPlusOne,iMinusOne)
!
!     This routine evaluates the overlap between two determinants using the
!     density different natural orbitals and the MO coefiicients of the final
!     state (<CMOsKet>). The work below closely follows the derivation of Eq.
!     (24) in Harb and Hratchian, J. Chem. Phys. 154, 084104 (2021).
!
!
!     - H. P. Hratchian, 2022
!
!
!     Variable Declarations
!
      implicit none
      type(MQC_Variable),intent(in)::SMatrix,SMatrixMinusHalf,  &
        diffDensityEVals,diffDensityEVecs,CMOsKet
      integer(kind=int64),intent(in)::NOccKet,nBasis
      type(MQC_Variable),intent(out)::detOverlap
      integer(kind=int64),intent(out)::nPlusOne,nMinusOne,iPlusOne,  &
        iMinusOne
!
      integer(kind=int64)::i
      real(kind=real64),dimension(:),allocatable::vector
      type(MQC_Variable)::TMatrix
      type(MQC_Variable)::tmpMQCvar,tmpMQCvar1,tmpMQCvar2,tmpMQCvar3
      logical::DEBUG=.false.
!
!
!     Start by forming the V matrix, Eq. (7) from NIO polestrength JCP paper.
!

!hph+
      iPlusOne = 0
      iMinusOne = 0
      nPlusOne = 0
      nMinusOne = 0
      vector = diffDensityEVals
      do i = 1,nBasis
        if(vector(i).ge.0.999) then
          vector(i) = float(0)
          iPlusOne = i
          nPlusOne = nPlusOne + 1
        endIf
      endDo
      do i = 1,nBasis
        if(vector(i).le.-0.999) then
          vector(i) = float(0)
          iMinusOne = i
          nMinusOne = nMinusOne + 1
        endIf
      endDo
!hph-

      tmpMQCvar = MatMul(SMatrixMinusHalf,diffDensityEVecs)
      if(DEBUG) call tmpMQCvar%print(header='V')
      TMatrix = MatMul(Transpose(CMOsKet),MatMul(SMatrix,tmpMQCvar))
      if(DEBUG) then
        call TMatrix%print(header='TMatrix')
        call mqc_print(MatMul(Transpose(TMatrix),TMatrix),header='TMatrix(t).TMatrix')
      endIf
      tmpMQCvar = MQC_Variable_SubMatrix(TMatrix,newrange1=[1,nOccKet])

!hph+
      tmpMQCvar2 = vector
      tmpMQCvar1 = MatMul(MatMul(Transpose(tmpMQCvar),tmpMQCvar),tmpMQCvar2%diag())
!      tmpMQCvar1 = MatMul(MatMul(Transpose(tmpMQCvar),tmpMQCvar),diffDensityEVals%diag())
!hph-

      tmpMQCvar2 = MQC_Variable_UnitMatrix(nBasis)
      tmpMQCvar3 = tmpMQCvar2 - tmpMQCvar1
      if(DEBUG) then
        call tmpMQCvar2%print(header='unit matrix (nBasis)')
        call tmpMQCvar1%print(header='Tt.T.diagE')
        call tmpMQCvar3%print(header='I-T(occ)(t).T(occ).delta')
      endIf
      detOverlap = tmpMQCvar3%det()
      call detOverlap%print(header='detOverlap')
!
      return
      end subroutine determinantOverlap

!
!
!
      end module nio_mod
