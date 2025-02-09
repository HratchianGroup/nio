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
      iPlusOne = 0
      iMinusOne = 0
      nPlusOne = 0
      nMinusOne = 0
      vector = diffDensityEVals
      do i = 1,nBasis
        if(vector(i).ge.0.999) then
          vector(i) = dfloat(0)
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
      tmpMQCvar = MatMul(SMatrixMinusHalf,diffDensityEVecs)
      if(DEBUG) call tmpMQCvar%print(header='V')
      TMatrix = MatMul(Transpose(CMOsKet),MatMul(SMatrix,tmpMQCvar))
      if(DEBUG) then
        call TMatrix%print(header='TMatrix')
        call mqc_print(MatMul(Transpose(TMatrix),TMatrix),header='TMatrix(t).TMatrix')
      endIf
      tmpMQCvar = TMatrix%subMatrix(newrange1=[1,nOccKet])
      tmpMQCvar2 = vector
      tmpMQCvar1 = MatMul(MatMul(Transpose(tmpMQCvar),tmpMQCvar),tmpMQCvar2%diag())
      tmpMQCvar2 = MQC_Variable_UnitMatrix(nBasis)
      tmpMQCvar3 = tmpMQCvar2 - tmpMQCvar1
      if(DEBUG) then
        call tmpMQCvar2%print(header='unit matrix (nBasis)')
        call tmpMQCvar1%print(header='Tt.T.diagE')
        call tmpMQCvar3%print(header='I-T(occ)(t).T(occ).delta')
      endIf
      detOverlap = tmpMQCvar3%det()
!
      return
      end subroutine determinantOverlap

!
!PROCEDURE categorizeDDNOs
      subroutine categorizeDDNOs(diffDensityEVals,infoDDNOs,nRelaxationDDNOs)
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
      integer(kind=int64),intent(out)::nRelaxationDDNOs
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
      nRelaxationDDNOs = 0
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
          nRelaxationDDNOs = nRelaxationDDNOs + 1
        elseIf(vector(i).le.-0.01) then
          notNull = .True.
          myInfoDDNOs(1,i) = -2
          nRelaxationDDNOs = nRelaxationDDNOs + 1
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
      call mqc_print(myInfoDDNOs,6,header='Here is DDNO info')
      write(*,*)
      write(*,*)
!
      infoDDNOs = myInfoDDNOs
      deallocate(vector,myInfoDDNOs)

!hph+
      write(*,*)
      write(*,*)' Hrant - nRelaxationDDNOs: ',nRelaxationDDNOs
      write(*,*)
!hph-

!
      return
      end subroutine categorizeDDNOs

!
!PROCEDURE formOverlapMatrix
      function getOverlapMatrix(fileinfo) result(SMatrix)
!
!     This function returns the AO overlap matrix. The overlap is either
!     retrieved from the file provided by the calling program unit through
!     argument <fileinfo> or it is formed from other information available on
!     the file. The output of the function is the overlap matrix as an
!     MQC_Variable type (MQC_Algebra2 type).
!
!
!     - H. P. Hratchian, 2023
!
!
!     Variable Declarations
!
      implicit none
      integer(kind=int64)::i
      real(kind=real64)::tmpReal
      logical found,formedS
      type(mqc_gaussian_unformatted_matrix_file)::fileinfo
      type(MQC_Variable)::orthogonalBasis,SMatrix,sigma
!
!     Begin by setting <formedS> to FALSE and then asking for the overlap matrix
!     from the provided data file (<fileinfo>). If it is available on the file,
!     we're all done.
!
      formedS = .false.
      call fileinfo%getArray('overlap',mqcVarOut=SMatrix,foundOut=found)
      formedS = found
      if(formedS) return
!
!     The AO overlap matrix isn't on the file. Loop for the orthogonal basis
!     matrix. If that's available we can build the overlap matrix using it. If
!     it isn't on the file, then we die (for now, anyway). 
!
      call fileinfo%getArray('ORTHOGONAL BASIS',mqcVarOut=orthogonalBasis,foundOut=found)
      if(found) then
        sigma = MatMul(Transpose(orthogonalBasis),orthogonalBasis)
        do i = 1,Size(sigma,1)
          tmpReal = sigma%getVal([ i,i ])
          if(abs(tmpReal).gt.1.d-5) then
            tmpReal = float(1)/tmpReal**2
          else
            tmpReal = 0.d0
          endIf
          call sigma%put(tmpReal,[ i,i ])
        endDo
        SMatrix = MatMul(MatMul(orthogonalBasis,sigma),Transpose(orthogonalBasis))
        formedS = .true.
      else
        call mqc_error('formOverlapMatrix: Cannot complete operation with the data available.')
      endIf
      if(.not.formedS)  &
        call mqc_error('formOverlapMatrix: Failed to form the overlap matrix.')
!
      return
      end function getOverlapMatrix

!
!
!
      end module nio_mod
