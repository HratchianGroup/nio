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
      logical::DEBUG=.false.
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
!PROCEDURE categorizeDDNOs
      subroutine categorizeDDNOs(diffDensityEVals,infoDDNOs)
!
!     This routine analyzes a list of DDNOs to categorize DDNOs and, when
!     appropriate, identify relaxation pairs and attachment/detachement pairs.
!
!
!     - H. P. Hratchian, 2022
!
!
!     Variable Declarations
!
      implicit none
      type(MQC_Variable),intent(in)::diffDensityEVals
      type(MQC_Variable),intent(out)::infoDDNOs
      integer(kind=int64)::i,j,nDDNOs,nNonNullDDNOs
      integer(kind=int64),dimension(:,:),allocatable::myInfoDDNOs
      real(kind=real64),dimension(:),allocatable::vector
      logical::notNull
!
!
!     Most of the work below is done using fortran intrinsic data types.
!
      nDDNOs = Size(diffDensityEVals)
      Allocate(myInfoDDNOs(2,nDDNOs))
      myInfoDDNOs = 0
      nNonNullDDNOs = 0
      vector = diffDensityEVals
      do i = 1,nDDNOs
        notNull = .False.
        if(vector(i).ge.0.999) then
          notNull = .True.
          myInfoDDNOs(1,i) = 1
        elseIf(vector(i).le.-0.999) then
          notNull = .True.
          myInfoDDNOs(1,i) = -1
        elseIf(vector(i).ge.0.01) then
          notNull = .True.
          myInfoDDNOs(1,i) = 2
        elseIf(vector(i).le.-0.01) then
          notNull = .True.
          myInfoDDNOs(1,i) = -2
        endIf
        if(notNull) then
          nNonNullDDNOs = nNonNullDDNOs + 1
          if(myInfoDDNOs(2,i).eq.0) then
 jLoop1:    do j = 1,i-1
              if(ABS(vector(i)+vector(j)).le.(100*MQC_small)) then
                myInfoDDNOs(2,i) = j
                myInfoDDNOs(2,j) = i
                exit jLoop1
              endIf
            endDo jLoop1
          endIf
          if(myInfoDDNOs(2,i).eq.0) then
 jLoop2:    do j = i+1,nDDNOs
              if(ABS(vector(i)+vector(j)).le.(100*MQC_small)) then
                myInfoDDNOs(2,i) = j
                myInfoDDNOs(2,j) = i
                exit jLoop2
              endIf
            endDo jLoop2
          endIf
        endIf
      endDo
      write(*,*)
      write(*,*)' Hrant - nNonNullDDNOs = ',nNonNullDDNOs
      call mqc_print(6,myInfoDDNOs,header='Here is DDNO info')
      write(*,*)
      write(*,*)
!
      deallocate(vector)
!
      return
      end subroutine categorizeDDNOs

!
!PROCEDURE projectDDNOs
      subroutine projectDDNOs(infoDDNOs,DDNOs,SMatrix,CMOs,nOcc)
!
!     This routine finds all relaxation DDNO pairs and rotates them into a new
!     linear combination where one projects fully onto the initial state's
!     occupied sub-space and one projects fully onto the initial state's virtual
!     sub-space.
!
!
!     - H. P. Hratchian, 2022
!
!
!     Variable Declarations
!
      implicit none
      type(MQC_Variable),intent(in)::infoDDNOs,DDNOs,SMatrix,CMOs,nOcc
      type(MQC_Variable)::tmp
      real(kind=real64),dimension(:),allocatable::tmpVector1
      real(kind=real64),dimension(:,:),allocatable::tmpMatrix,occVirtPops
      integer(kind=int64)::no,i
!
!
!     Determine the projection of DDNOs onto MOs.
!
      tmpMatrix = MatMul(TRANSPOSE(CMOs),MatMul(SMatrix,DDNOs))
      call mqc_print(6,tmpMatrix,header='DDNOs projected into MO basis.')
      no = 8
      Allocate(tmpVector1(SIZE(DDNOs,1)))
      Allocate(occVirtPops(3,SIZE(DDNOs,1)))
      tmpVector1 = 0
      occVirtPops = 0
      do i = 1,SIZE(DDNOs,1)
        tmpVector1 = 0
        tmpVector1 = tmpMatrix(:no,i)
        occVirtPops(1,i) = dot_product(tmpVector1,tmpVector1)
        tmpVector1 = 0
        tmpVector1 = tmpMatrix((no+1):,i)
        occVirtPops(2,i) = dot_product(tmpVector1,tmpVector1)
        tmpVector1 = tmpMatrix(:,i)
        occVirtPops(3,i) = dot_product(tmpVector1,tmpVector1)
      endDo
      call mqc_print(6,occVirtPops,header='occ | virt | all pops...')
!
      return
      end subroutine projectDDNOs

!
!
!
      end module nio_mod
