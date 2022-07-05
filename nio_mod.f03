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
      tmpMQCvar = TMatrix%subMatrix(newrange1=[1,nOccKet])

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
      call mqc_print(6,myInfoDDNOs,header='Here is DDNO info')
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
!PROCEDURE projectDDNOs
      subroutine projectDDNOs(infoDDNOs,DDNOs,SMatrix,CMOs,nOcc,  &
        nRelaxationDDNOs,differenceDensity)
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
      type(MQC_Variable),intent(in)::infoDDNOs,DDNOs,SMatrix,CMOs
      integer(kind=int64),intent(in)::nOcc,nRelaxationDDNOs
      type(MQC_Variable)::differenceDensity
      type(MQC_Variable)::tmp1,tmp2
      real(kind=real64)::vectorNorm1,vectorNorm2
      real(kind=real64),dimension(:),allocatable::tmpVector1,tmpVector2,  &
        tmpVector3,tmpVector4
      real(kind=real64),dimension(:,:),allocatable::tmpMatrix,tmpMatrix2,  &
        tmpMatrix3,tmpMatrix4,occVirtPops,  &
        tmpVectors,intrinsicOverlap,intrinsicCMO
      integer(kind=int64)::i,j,iRelaxPartner,nRelaxFound
      integer(kind=int64)::no
!
      no = nOcc
      write(*,*)
      write(*,*)
      write(*,*)' Hrant - no = ',no
      write(*,*)
!
!
!     Determine the projection of DDNOs onto MOs.
!

!hph+
!      tmp1 = MatMul(TRANSPOSE(SMatrix),CMOs)
!      call tmp1%print(header='St.C')
!      tmp1 = MatMul(SMatrix,CMOs)
!      call tmp1%print(header='S.C')
!      tmp1 = MatMul(SMatrix,MatMul(CMOs,TRANSPOSE(CMOs)))
!      call tmp1%print(header='S.C.Ct = 1 ???')
!
!      call DDNOs%print(header='Input DDNOs in the AO Basis.')
!      tmp1 = MatMul(TRANSPOSE(CMOs),MatMul(SMatrix,DDNOs))
!      call tmp1%print(header='DDNOs in MO Basis')
!      tmp2 = MatMul(CMOs,tmp1)
!      call tmp2%print(header='DDNOs back-rotating to AO Basis')
!hph-


      tmpMatrix = MatMul(TRANSPOSE(CMOs),MatMul(SMatrix,DDNOs))
      call mqc_print(6,tmpMatrix,header='DDNOs projected into MO basis.')
      Allocate(tmpVector1(SIZE(DDNOs,1)),tmpVector2(SIZE(DDNOs,1)))
      Allocate(tmpVector3(SIZE(DDNOs,1)),tmpVector4(SIZE(DDNOs,1)))
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
!     Try to project DDNOs into occ or virt sub-space based on which the DDNO
!     already has majority overlap.
!
      Allocate(tmpMatrix2(SIZE(DDNOs,1),nRelaxationDDNOs))
      nRelaxFound = 0
      do i = 1,SIZE(DDNOs,2)
        if(Abs(Int(MQC_Variable_get_MQC(infoDDNOs,[1,i]))).eq.2) then
          nRelaxFound = nRelaxFound + 1
          tmpVector1 = tmpMatrix(:,i)
          if(occVirtPops(1,i).gt.occVirtPops(2,i)) then
            tmpVector1(no+1:) = 0.0
          else
            tmpVector1(1:no) = 0.0
          endIf
          vectorNorm1 = dot_product(tmpVector1,tmpVector1)
          if(abs(vectorNorm1).le.MQC_Small) then
            vectorNorm1 = 0.0
          else
            vectorNorm1 = 1.0/sqrt(vectorNorm1)
          endIf
          tmpMatrix2(:,nRelaxFound) = vectorNorm1*tmpVector1
          write(*,*)' Hrant - nRelaxFound,i = ',nRelaxFound,i
        endIf
      endDo
      call mqc_print(6,tmpMatrix2,header='Here are the relaxation DDNOs projected onto Occ/Virt sub-spaces.')
      call mqc_print(6,MatMul(TRANSPOSE(tmpMatrix2),tmpMatrix2),header='Relaxations DDNOs Occ/Virt Overlaps')


      return
      call mqc_error('hph-stop')




!
!     Build an overlap matrix of the NIOs with their occupied sub-spaces only.
!
      call infoDDNOs%print(header='infoDDNOs')
      tmpMatrix2 = tmpMatrix
      tmpMatrix2(no+1:,:) = float(0)
      call mqc_print(6,tmpMatrix2,header='Occ DDNO Matrix:')
      call mqc_dsyev_eigensystem_symmFull(tmpMatrix2,tmpVector1,tmpMatrix3)
      call mqc_print(6,tmpVector1,header='Occ DDNO E-Values')
      call mqc_print(6,tmpMatrix3,header='Occ DDNO E-Vectors')
      write(*,*)
      write(*,*)
      DeAllocate(tmpMatrix2)
      Allocate(tmpMatrix2(SIZE(DDNOs,1),nRelaxationDDNOs))
      nRelaxFound = 0
      do i = 1,SIZE(DDNOs,2)
        if(Abs(Int(MQC_Variable_get_MQC(infoDDNOs,[1,i]))).eq.2) then
          nRelaxFound = nRelaxFound + 1
          tmpVector1 = tmpMatrix(:,i)
          tmpVector1(no+1:) = 0.0
          vectorNorm1 = dot_product(tmpVector1,tmpVector1)
          if(abs(vectorNorm1).le.MQC_Small) then
            vectorNorm1 = 0.0
          else
            vectorNorm1 = 1.0/sqrt(vectorNorm1)
          endIf
          tmpMatrix2(:,nRelaxFound) = vectorNorm1*tmpVector1
          write(*,*)' Hrant - nRelaxFound,i = ',nRelaxFound,i
        endIf
      endDo
      call mqc_print(6,tmpMatrix2,header='Here are the relaxation DDNOs projected onto Occ sub-space.')
      call mqc_print(6,MatMul(TRANSPOSE(tmpMatrix2),tmpMatrix2),header='Relaxations DDNOs Occ Overlaps')

      do i = 2,nRelaxFound
        tmpVector1 = tmpMatrix2(:,i)
        do j = 1,i-1
          tmpVector2 = tmpMatrix2(:,j)
          if(dot_product(tmpVector1,tmpVector1).gt.mqc_small.and.  &
            dot_product(tmpVector2,tmpVector2).gt.mqc_small) then
            tmpVector1 = tmpVector1 - dot_product(tmpVector2,tmpVector1)*tmpVector2
            call mqc_normalizeVector(tmpVector1)
          endIf
        endDo
        tmpMatrix2(:,i) = tmpVector1
      endDo

      write(*,*)
      write(*,*)
      call mqc_print(6,tmpMatrix2,header='After GS -- Here are the relaxation DDNOs projected onto Occ sub-space.')
      call mqc_print(6,MatMul(TRANSPOSE(tmpMatrix2),tmpMatrix2),header='After GS -- Relaxations DDNOs Occ Overlaps')




      DeAllocate(tmpMatrix3)
      Allocate(tmpMatrix3(nRelaxationDDNOs,nRelaxationDDNOs))
      tmpMatrix3 = MatMul(TRANSPOSE(tmpMatrix2),tmpMatrix2)
      call mqc_dsyev_eigensystem_symmFull(tmpMatrix3,tmpVector1,tmpMatrix4)
      call mqc_print(6,tmpVector1,header='Overlap E-Values')
      call mqc_print(6,tmpMatrix4,header='Overlap E-Vectors')
      write(*,*)
      write(*,*)

      if(ALLOCATED(tmpMatrix)) DeAllocate(tmpMatrix)
      tmpMatrix = differenceDensity
      write(*,*)
      call mqc_print(6,MatMul(Transpose(tmpMatrix2),MatMul(tmpMatrix,tmpMatrix2)),  &
        header='vt.DeltaP.v')
      write(*,*)
      write(*,*)


      return





      tmpMatrix3 = float(0)
      do i = 1,SIZE(DDNOs,2)
        vectorNorm1 = dot_product(tmpMatrix2(:,i),tmpMatrix2(:,i))
        if(abs(vectorNorm1).lt.mqc_small) then
          vectorNorm1 = float(0)
        elseIf(vectorNorm1.ge.0.0) then
          vectorNorm1 = 1.0/sqrt(vectorNorm1)
        else
          call mqc_error('Hrant - vectorNorm1 < 0 ... WTF!!!')
        endIf
        tmpMatrix2(:,i) = vectorNorm1*tmpMatrix2(:,i)
      endDo
      call mqc_print(6,MatMul(Transpose(tmpMatrix2),tmpMatrix2),header='Occ-Occ Overlap Matrix:')

      call mqc_dsyev_eigensystem_symmFull(tmpMatrix2,tmpVector1,tmpMatrix3)
      call mqc_print(6,tmpVector1,header='Overlap E-Values')
      call mqc_print(6,tmpMatrix3,header='Overlap E-Vectors')
      write(*,*)
      write(*,*)

      return


!
!     Do the projections over relaxation pairs..
!
      call infoDDNOs%print(header='infoDDNOs')
      do i = 1,SIZE(DDNOs,2)
        if(Int(MQC_Variable_get_MQC(infoDDNOs,[1,i])).eq.2) then
          iRelaxPartner = Int(MQC_Variable_get_MQC(infoDDNOs,[2,i]))
          write(*,*)' Found a relaxation DDNO: ',i,'  Partner: ',iRelaxPartner
          tmpVector1 = 0
          tmpVector2 = 0
          tmpVector1(:no) = tmpMatrix(:no,i)
          tmpVector2(:no) = tmpMatrix(:no,iRelaxPartner)
          write(*,*)' i.iRelaxPartner = ',dot_product(tmpMatrix(:,i),tmpMatrix(:,iRelaxPartner))
          write(*,*)' vec1.vec2 (Occ)  = ',dot_product(tmpVector1,tmpVector2)
          vectorNorm1 = sqrt(dot_product(tmpVector1,tmpVector1))
          vectorNorm2 = sqrt(dot_product(tmpVector2,tmpVector2))
          write(*,*)' Angle(vec1.vec2) = ',ACos(dot_product(tmpVector1,tmpVector2)/(vectorNorm1*vectorNorm2))
          call mqc_print(6,tmpVector1,header='vec 1:')
          call mqc_print(6,tmpVector2,header='vec 2:')
          tmpVector3 = 1/sqrt(occVirtPops(1,i))*tmpVector1 +  &
            1/sqrt(occVirtPops(2,iRelaxPartner))*tmpVector2
          tmpVector1 = 0
          tmpVector2 = 0
          tmpVector1(no:) = tmpMatrix(no:,i)
          tmpVector2(no:) = tmpMatrix(no:,iRelaxPartner)
          write(*,*)' vec1.vec2 (Virt) = ',dot_product(tmpVector1,tmpVector2)
          tmpVector4 = tmpVector1 + tmpVector2
          call mqc_print(6,tmpVector3,header='vec 3:')
          call mqc_print(6,tmpVector4,header='vec 4:')
          write(*,*)
          write(*,*)' <d1` | d1 `> = ',dot_product(tmpVector3,tmpVector3)
          write(*,*)' <d2` | d2 `> = ',dot_product(tmpVector4,tmpVector4)
          write(*,*)' <d1` | d2 `> = ',dot_product(tmpVector3,tmpVector4)
          write(*,*)' <d2` | d1 `> = ',dot_product(tmpVector4,tmpVector3)
          write(*,*)
        endIf
      endDo
      DeAllocate(tmpVector2,tmpVector3,tmpVector4)
!
!     Try to do rotation of a relaxation DDNO pair.
!
      Allocate(tmpVector2(SIZE(DDNOs,1)))
      Allocate(tmpVectors(SIZE(DDNOs,1),2))
      call DDNOs%print(header='DDNOs')
      tmpMatrix = MQC_Variable_SubMatrix(DDNOs,newrange2=[2,2])
      tmpVector1 = RESHAPE(tmpMatrix,[SIZE(tmpMatrix)])
      call mqc_print(6,tmpMatrix,header='tmpMatrix')
      call mqc_print(6,tmpVector1,header='tmpVector1')
      tmpMatrix = MQC_Variable_SubMatrix(DDNOs,newrange2=[11,11])
      tmpVector2 = RESHAPE(tmpMatrix,[SIZE(tmpMatrix)])
      call mqc_print(6,tmpVector2,header='tmpVector2')
      Allocate(tmpVector3(SIZE(DDNOs,1)),tmpVector4(SIZE(DDNOs,1)))
      tmpVector3 = tmpVector1 - SQRT(occVirtPops(2,2))*tmpVector1
      tmpVector4 = tmpVector2 - occVirtPops(2,11)*tmpVector2
      call mqc_print(6,tmpVector3,header='tmpVector3')
      tmpVectors(:,1) = tmpVector3
      tmpVectors(:,2) = tmpVector4

      Allocate(intrinsicOverlap(SIZE(DDNOs,1),SIZE(DDNOs,1)),  &
        intrinsicCMO(SIZE(DDNOs,1),SIZE(DDNOs,1)))
      intrinsicOverlap = SMatrix
      intrinsicCMO = CMOs

      Call mqc_print(6,MatMul(TRANSPOSE(tmpVectors),MatMul(intrinsicOverlap,tmpVectors)),header='Vt.S.V')
      Call mqc_print(6,MatMul(TRANSPOSE(intrinsicCMO),MatMul(intrinsicOverlap,tmpVectors)),header='Ct.S.V')


      write(*,*)
      write(*,*)
!
      return
      end subroutine projectDDNOs

!
!
!
      end module nio_mod
