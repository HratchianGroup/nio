      program nio
!
!     This program carries out Natural Ionization Orbital analysis.
!
!     -H. P. Hratchian, 2022.
!
!
!     USE Connections
!
      use nio_mod
!
!     Variable Declarations
!
      implicit none
      integer(kind=int64)::nCommands,iPrint=1,i,nAtoms,nAtoms2,  &
        nBasis,nBasis2,nBasisUse,nBasisUse2,nEl1,nEl2,nElAlpha1,  &
        nElBeta1,nElAlpha2,NElBeta2,nPlusOneAlpha,nMinusOneAlpha,  &
        iPlusOneAlpha,iMinusOneAlpha,nPlusOneBeta,nMinusOneBeta,  &
        iPlusOneBeta,iMinusOneBeta,nPlusOne,nMinusOne,  &
        nRelaxationDDNOsAlpha,nRelaxationDDNOsBeta
      integer(kind=int64),allocatable,dimension(:)::tmpVectorInt
      real(kind=real64)::scfEnergy1,scfEnergy2,deltaSCFEnergy,temp
      real(kind=real64),dimension(3)::transitionDipole
      real(kind=real64),allocatable,dimension(:)::cart1,cart2,  &
        tmpVector
      real(kind=real64),allocatable,dimension(:,:)::tmpMatrix1,  &
        tmpMatrix2,tmpMatrix3
      character(len=512)::matrixFilename1,matrixFilename2,  &
        matrixFilenameOut
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile1,  &
        GMatrixFile2,GMatrixFileOut
      type(MQC_Variable)::DDNOsAlpha,DDNOsBeta,pDDNO,hDDNO,  &
        dipoleStrength,TDOverlapA,TDOverlapB,TDparticleHoleMag,  &
        oscillatorStrength
      type(MQC_Variable)::infoDDNOsAlpha,infoDDNOsBeta
      type(MQC_Variable)::SMatrixAO,SMatrixEVecs,SMatrixEVals,  &
        SMatrixAOHalf,SMatrixAOMinusHalf
      type(MQC_Variable)::PMatrixAlpha1,PMatrixBeta1,PMatrixTotal1,  &
        PMatrixAlpha2,PMatrixBeta2,PMatrixTotal2,diffDensityAlpha,  &
        diffDensityBeta,diffDensityAlphaEVecs,diffDensityAlphaEVals,  &
        diffDensityBetaEVecs,diffDensityBetaEVals,attachmentDensity,  &
        detachmentDensity,attachmentNumberAlpha,detachmentNumberAlpha,  &
        attachmentNumberBeta,detachmentNumberBeta,XnTot,RnTot,XnOO,  &
        XnOV,XnVO,XnVV,RnOO,RnOV,RnVO,RnVV
      type(MQC_Variable)::CAlpha1,CBeta1,CAlpha2,CBeta2,TAlpha,TBeta
      type(MQC_Variable)::dipoleAOx,dipoleAOy,dipoleAOz
      type(MQC_Variable)::tmpMQCvar,tmpMQCvar1,tmpMQCvar2,tmpMQCvar3,  &
        tmpMQCvar4,POcc,PVirt,PDeltaMO,PDeltaOO,PDeltaVV,PDeltaOV,  &
        PDeltaVO
      logical::doTestCode=.false.,doMatrixFileOut=.false.,isNIO,isDDNO
!
!     Format Statements
!
 1000 Format(1x,'Enter Program NIO.')
 1010 Format(1x,'Matrix File 1:      ',A,/,  &
             1x,'Matrix File 2:      ',A,/,  &
             1x,'Output Matrix File: ',A,/)
 1020 Format(1x,'Matrix File 1:      ',A,/,  &
             1x,'Matrix File 2:      ',A,/)
 1100 Format(1x,'nAtoms=',I4,3x,'nBasis   =',I4,3x,'nBasisUse=',I4,/,  &
             1x,'nEl1  =',I4,3x,'nElAlpha1=',I4,3x,'nElBeta  =',I4,/,  &
             1x,'nEl2  =',I4,3x,'nElAlpha2=',I4,3x,'nElBeta  =',I4,/)
 1200 Format(/,1x,'SCF Energies',/,  &
        3x,'SCF 1:     ',F15.8,' a.u.',/,  &
        3x,'SCF 2:     ',F15.8,' a.u.',/,  &
        3x,'Delta-SCF: ',F15.8,' a.u.',2x,'=',2x,F15.8,' eV',/,  &
        36x,'=',2x,F15.8,' cm^-1',/,  &
        36x,'=',2x,F15.8,' nm')
 1500 Format(/,1x,'NIO Polestrength = ',F9.6,/,  &
        3x,'Alpha Polestrength = ',F9.6,3x,'Beta Polestrength = ',F9.6)
 1600 Format(/,1x,'Overlap between Delta-SCF states = ',F9.6,/,  &
        3x,'Alpha Overlap = ',F9.6,3x,'Beta Overlap = ',F9.6)
 2000 Format(/,1x,'Attachment/Detachment Density Analysis',A,/,  &
        3x,'ALPHA: Attachement Number: ',F6.3,/,  &
        10x,'Detachment  Number: ',F6.3,/,  &
        3x,'BETA : Attachement Number: ',F6.3,/,  &
        10x,'Detachment  Number: ',F6.3)
 3000 Format(1x,'DDNO/NIO Excitation/Detachment Number',/,  &
        3x,'ALPHA: Attachment Number: ',I2,/,  &
        10x,'Detachment Number: ',I2,/,  &
        3x,'BETA : Attachment Number: ',I2,/,  &
        10x,'Detachment Number: ',I2)
 3100 Format(1x,'isNIO=',L1,3x,'isDDNO=',L1)
 4000 Format(/,1x,'Relaxation Number: ',F9.6)
 4010 Format(4x,'Alpha Contribution: ',F9.6,  &
        /,4x,'Beta Contribution: ',F9.6,/,/)
 7000 Format(/,1x,'Writing an output matrix file. Filename: ',A)
 8999 Format(/,1x,'END OF NIO PROGRAM')
 9000 Format(/,1x,'NIO has been compiled using an unsupported version of MQCPack.',/)
!
!
      write(IOut,1000)
      call mqc_version_print(iOut)
!
!     Do a check of the mqcPack version the program was built against to ensure
!     it's a supported version.
!
      if(.not.mqc_version_check(newerThanMajor=24,newerThanMinor=0,  &
        newerThanRevision=0)) then
        write(iOut,9000)
        goto 999
      endIf
!
!     Open the Gaussian matrix file and load the number of atomic centers.

      nCommands = command_argument_count()
      if(nCommands.lt.2)  &
        call mqc_error('Two input Gaussian matrix files must be provided in the command line.')
      call get_command_argument(1,matrixFilename1)
      call get_command_argument(2,matrixFilename2)
      call GMatrixFile1%load(matrixFilename1)
      call GMatrixFile2%load(matrixFilename2)
      if(nCommands.eq.3) then
        call get_command_argument(3,matrixFilenameOut)
        doMatrixFileOut = .true.
      endIf
      if(doMatrixFileOut) then
        write(IOut,1010) TRIM(matrixFilename1),TRIM(matrixFilename2),  &
          TRIM(matrixFilenameOut)
      else
        write(IOut,1020) TRIM(matrixFilename1),TRIM(matrixFilename2)
      endIf
!
!     Do some consistency checks and load the number of atoms, basis functions,
!     and linearly independent basis functions.
!
      nAtoms  = GMatrixFile1%getVal('nAtoms')
      nAtoms2 = GMatrixFile2%getVal('nAtoms')
      if(nAtoms.ne.nAtoms2) call mqc_error('nAtoms must be the same in the two matrix file!')
      nBasis  = GMatrixFile1%getVal('nBasis')
      nBasis2 = GMatrixFile2%getVal('nBasis')
      if(nBasis.ne.nBasis2) call mqc_error('nBasis must be the same in the two matrix file!')
      nBasisUse  = GMatrixFile1%getVal('nBasisUse')
      nBasisUse2 = GMatrixFile2%getVal('nBasisUse')
      if(nBasisUse.ne.nBasisUse2)  &
        call mqc_error('nBasisUse must be the same in the two matrix file!')
      nEl1      = GMatrixFile1%getVal('nElectrons')
      nElAlpha1 = GMatrixFile1%getVal('nAlpha')
      nElBeta1  = GMatrixFile1%getVal('nBeta')
      nEl2      = GMatrixFile2%getVal('nElectrons')
      nElAlpha2 = GMatrixFile2%getVal('nAlpha')
      nElBeta2  = GMatrixFile2%getVal('nBeta')
      write(IOut,1100) nAtoms,nBasis,nBasisUse,nEl1,nElAlpha1,nElBeta1,  &
        nEl2,nElAlpha2,nElBeta2
!
!     Another consistency check is to ensure the atomic coordinates are the same
!     for the two FAFs.
!
      cart1 = GMatrixFile1%getAtomCarts()
      cart2 = GMatrixFile2%getAtomCarts()
      Allocate(tmpVector(3*nAtoms))
      tmpVector = cart1-cart2
      if(MaxVal(abs(tmpVector)).gt.0.0001) then
        call mqc_print(tmpVector,iout=iOut,header='cart1 - cart2')
        call mqc_error('Coordinates of the two matrix files must be the same.')
      endIf
      DeAllocate(tmpVector)
!
!     Pull the energies from the two jobs and report the delta-E value in a few
!     units.
!
      scfEnergy1 = GMatrixFile1%getValReal('scfEnergy')
      scfEnergy2 = GMatrixFile2%getValReal('scfEnergy')
      deltaSCFEnergy = scfEnergy2-scfEnergy1
      write(iOut,1200) scfEnergy1,scfEnergy2,deltaSCFEnergy,  &
        deltaSCFEnergy*evPHartree,deltaSCFEnergy*cmM1PHartree,  &
        deltaSCFEnergy*evPHartree*nmPev
!
!     Load the atomic orbital overlap matrix and form S^(1/2) and S^(-1/2).
!
      SMatrixAO = getOverlapMatrix(GMatrixFile1)
      call SMatrixAO%eigen(SMatrixEVals,SMatrixEVecs)
      if(DEBUG) then
        call SMatrixAO%print(header='Overlap Matrix')
        call SMatrixEVals%print(header='S matrix eigen-values:')
        call mqc_print(MatMul(Transpose(SMatrixEVecs),SMatrixEVecs),header='SEVecs(t).SEVecs')
        call mqc_print(MatMul(MatMul(SMatrixEVecs,SMatrixEVals%diag()),TRANSPOSE(SMatrixEVecs)),6,'U.lambda.Ut')
      endIf
      tmpMQCvar = SMatrixEVals%rpower(0.5)
      SMAtrixAOHalf = MatMul(MatMul(SMatrixEVecs,tmpMQCvar%diag()),TRANSPOSE(SMatrixEVecs))
      tmpMQCvar = SMatrixEVals%rpower(-0.5)
      SMAtrixAOMinusHalf = MatMul(MatMul(SMatrixEVecs,tmpMQCvar%diag()),TRANSPOSE(SMatrixEVecs))
      if(DEBUG) then
        call SMAtrixAOHalf%print(header='S**(1/2)')
        call SMatrixAOMinusHalf%print(header='S**(-1/2)')
      endIf
!
!     Load the density matrices. The code below treats all systems as open
!     shell, so closed shell results are handled by copying the density matrix
!     from restricted calculations into alpha and beta density matrix arrays.
!
      call GMatrixFile1%getArray('ALPHA DENSITY MATRIX',mqcVarOut=PMatrixAlpha1)
      call GMatrixFile1%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=CAlpha1)
      if(GMatrixFile1%isUnrestricted()) then
        call GMatrixFile1%getArray('BETA DENSITY MATRIX',mqcVarOut=PMatrixBeta1)
        call GMatrixFile1%getArray('BETA MO COEFFICIENTS',mqcVarOut=CBeta1)
      else
        PMatrixBeta1  = PMatrixAlpha1
        CBeta1 = CAlpha1
      endIf
      PMatrixTotal1 = PMatrixAlpha1+PMatrixBeta1
      call GMatrixFile2%getArray('ALPHA DENSITY MATRIX',mqcVarOut=PMatrixAlpha2)
      call GMatrixFile2%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=CAlpha2)
      if(GMatrixFile2%isUnrestricted()) then
        call GMatrixFile2%getArray('BETA DENSITY MATRIX',mqcVarOut=PMatrixBeta2)
        call GMatrixFile2%getArray('BETA MO COEFFICIENTS',mqcVarOut=CBeta2)
      else
        PMatrixBeta2  = PMatrixAlpha2
        CBeta2 = CAlpha2
      endIf
      PMatrixTotal2 = PMatrixAlpha2+PMatrixBeta2
!
!     Form the difference density and construct the DDNOs, which are NIOs in
!     electron detachment cases.
!
      diffDensityAlpha = PMatrixAlpha2-PMatrixAlpha1
      diffDensityBeta  = PMatrixBeta2-PMatrixBeta1
      if(iPrint.ge.1.or.DEBUG) then
        call mqc_print(contraction(diffDensityAlpha,SMatrixAO),header='DP(alpha).S = ')
        call mqc_print(contraction(diffDensityBeta,SMatrixAO),header='DP(beta).S  = ')
      endIf
      tmpMQCvar = MatMul(SMatrixAOHalf,MatMul(diffDensityAlpha,SMatrixAOHalf))
      call tmpMQCvar%eigen(diffDensityAlphaEVals,diffDensityAlphaEVecs)
      tmpMQCvar = MatMul(SMatrixAOHalf,MatMul(diffDensityBeta,SMatrixAOHalf))
      call tmpMQCvar%eigen(diffDensityBetaEVals,diffDensityBetaEVecs)
      if(iPrint.ge.1.or.DEBUG) then
        call diffDensityAlphaEVals%print(header='Alpha Occupation Change Values')
        call diffDensityBetaEVals%print(header='Beta Occupation Change Value')
      endIf
!
!     Order the difference density eigenvalues and eigenvectors so that the
!     largest change occupation pair are at the HOMO-LUMO gap.
!
!     First, take care of the alpha orbitals.
      tmpVector = diffDensityAlphaEVals
      Allocate(tmpVectorInt(SIZE(tmpVector)))
      call sort(tmpVector(1:nElAlpha1),map=tmpVectorInt(1:nElAlpha1),  &
        sortListIn=.true.,reverse=.true.)
      call sort(tmpVector(nElAlpha1+1:),map=tmpVectorInt(nElAlpha1+1:),  &
        sortListIn=.true.,reverse=.true.)
      tmpVectorInt(nElAlpha1+1:) = tmpVectorInt(nElAlpha1+1:) + nElAlpha1
      diffDensityAlphaEVals = tmpVector
      tmpMatrix1 = diffDensityAlphaEVecs
      call mqc_matrixOrderedColumns_real(tmpMatrix1,tmpVectorInt)
      diffDensityAlphaEVecs = tmpMatrix1
!
!     Now, take care of the beta orbitals.
      tmpVector = diffDensityBetaEVals
      call sort(tmpVector(1:nElBeta1),map=tmpVectorInt(1:nElBeta1),  &
        sortListIn=.true.,reverse=.true.)
      call sort(tmpVector(nElBeta1+1:),map=tmpVectorInt(nElBeta1+1:),  &
        sortListIn=.true.,reverse=.true.)
      tmpVectorInt(nElBeta1+1:) = tmpVectorInt(nElBeta1+1:) + nElBeta1
      diffDensityBetaEVals = tmpVector
      tmpMatrix1 = diffDensityBetaEVecs
      call mqc_matrixOrderedColumns_real(tmpMatrix1,tmpVectorInt)
      diffDensityBetaEVecs = tmpMatrix1
!
!     Back transform the difference density natural orbitals to the AO basis.
!
      DDNOsAlpha = MatMul(SMatrixAOMinusHalf,diffDensityAlphaEVecs)
      DDNOsBeta  = MatMul(SMatrixAOMinusHalf,diffDensityBetaEVecs)

!hph+
!
!     This block of code carries out various numerical tests of theory I've been
!     working through.
!
      write(*,*)
      write(*,*)
      write(*,*)' Hrant - Tests...'
      write(*,*)
      call mqc_print(PMatrixAlpha1,iOut,header='P1_alpha')
      call mqc_print(MatMul(PMatrixAlpha1,MatMul(SMatrixAO,PMatrixAlpha1)),iOut,header='P.S.P')
      tmpMQCvar = MatMul(Transpose(CAlpha1),MatMul(PMatrixAlpha1,CAlpha1))
      call mqc_print(tmpMQCvar,iOut,header='Ct.P.C')
      tmpMQCvar = MatMul(Transpose(CAlpha1),MatMul(SMatrixAO,CAlpha1))
      call mqc_print(tmpMQCvar,iOut,header='Ct.S.C')
      tmpMQCvar = CAlpha1%subMatrix(newRange2=[1,nElAlpha1])
      call mqc_print(CAlpha1,iOut,header='C')
      call mqc_print(tmpMQCvar,iOut,header='C(occ)')
      call mqc_print(MatMul(tmpMQCvar,Transpose(tmpMQCvar)),iOut,header='C(occ).C(occ)T')
      call mqc_print(MatMul(Transpose(CAlpha1),CAlpha1),iOut,header='CT.C')




      tmpMQCvar  = MatMul(Transpose(CAlpha1),MatMul(SMatrixAO,PMatrixAlpha1))
      tmpMQCvar1 = MatMul(tmpMQCvar,MatMul(SMatrixAO,CAlpha1))
      call mqc_print(tmpMQCvar1,iOut,header='CT.S.P1ao.S.C')

      POcc = tmpMQCvar1
      call POcc%trimZero()
      PVirt = MQC_Variable_UnitMatrix(Size(POcc,1))-POcc
      call PVirt%trimZero()
      call mqc_print(POcc,iOut,header='POcc')
      call mqc_print(PVirt,iOut,header='PVirt')

      tmpMQCvar  = MatMul(Transpose(CAlpha1),MatMul(SMatrixAO,PMatrixAlpha2))
      tmpMQCvar2 = MatMul(tmpMQCvar,MatMul(SMatrixAO,CAlpha1))
      call mqc_print(tmpMQCvar2,iOut,header='CT.S.P2ao.S.C')
      tmpMQCvar  = MatMul(Transpose(CAlpha1),MatMul(SMatrixAO,(PMatrixAlpha2-PMatrixAlpha1)))
      PDeltaMO = MatMul(tmpMQCvar,MatMul(SMatrixAO,CAlpha1))
      call PDeltaMO%trimZero()
      call mqc_print(PDeltaMO,iOut,header='CT.S.Delta.S.C')


!
!     Build the excitation matrix and the relaxation matrix in the MO basis.
!
      call PDeltaMO%eigen(diffDensityAlphaEVals,diffDensityAlphaEVecs)
      call diffDensityAlphaEVals%print(header='Alpha Occupation Change Values MO Basis')
      if(Allocated(tmpMatrix1)) deAllocate(tmpMatrix1)
      if(Allocated(tmpMatrix2)) deAllocate(tmpMatrix2)
      if(Allocated(tmpMatrix3)) deAllocate(tmpMatrix3)
      Allocate(tmpMatrix1(nBasis,nBasis),tmpMatrix2(nBasis,nBasis),  &
        tmpMatrix3(nBasis,nBasis))
#ifdef PGI
      tmpMatrix1 = dfloat(0)
      tmpMatrix2 = dfloat(0)
#else
      tmpMatrix1 = float(0)
      tmpMatrix2 = float(0)
#endif
      do i = 1,nBasis
        if(abs(float(MQC_Variable_get_MQC(diffDensityAlphaEVals,[i]))).gt.  &
          (1.0-MQC_small)) then
          write(iOut,*)' i = ',i
          write(iOut,*)'    eVal = ',float(MQC_Variable_get_MQC(diffDensityAlphaEVals,[i]))
          write(iOut,*)
          tmpVector = diffDensityAlphaEVecs%column(i)
          tmpMatrix1 = tmpMatrix1 +  &
            mqc_outerProduct_real(tmpVector,tmpVector,  &
            float(MQC_Variable_get_MQC(diffDensityAlphaEVals,[i])))
        elseIf(abs(float(MQC_Variable_get_MQC(diffDensityAlphaEVals,[i]))).gt.  &
          MQC_small) then
          tmpVector = diffDensityAlphaEVecs%column(i)
          tmpMatrix2 = tmpMatrix2 +  &
            mqc_outerProduct_real(tmpVector,tmpVector,  &
            float(MQC_Variable_get_MQC(diffDensityAlphaEVals,[i])))
        endIf
      endDo
      XnTot = tmpMatrix1
      RnTot = tmpMatrix2
      XnOO = MatMul(POcc,MatMul(XnTot,POcc))
      XnOV = MatMul(POcc,MatMul(XnTot,PVirt))
      XnVO = MatMul(PVirt,MatMul(XnTot,POcc))
      XnVV = MatMul(PVirt,MatMul(XnTot,PVirt))
      if(DEBUG) then
        call XnOO%print(header='XnOO')
        call XnOV%print(header='XnOV')
        call XnVO%print(header='XnVO')
        call XnVV%print(header='XnVV')
      endIf
      write(iOut,*)
      write(iOut,*)' Tr(XnTot)= ',XnTot%trace()
      write(iOut,*)
      write(iOut,*)' Tr(XnOO) = ',XnOO%trace()
      write(iOut,*)' Tr(XnOV) = ',XnOV%trace()
      write(iOut,*)' Tr(XnVO) = ',XnVO%trace()
      write(iOut,*)' Tr(XnVV) = ',XnVV%trace()
      write(iOut,*)
      write(iOut,*)
      write(iOut,*)
      RnOO = MatMul(POcc,MatMul(RnTot,POcc))
      RnOV = MatMul(POcc,MatMul(RnTot,PVirt))
      RnVO = MatMul(PVirt,MatMul(RnTot,POcc))
      RnVV = MatMul(PVirt,MatMul(RnTot,PVirt))
      if(DEBUG) then
        call RnOO%print(header='RnOO')
        call RnOV%print(header='RnOV')
        call RnVO%print(header='RnVO')
        call RnVV%print(header='RnVV')
      endIf
      write(iOut,*)
      write(iOut,*)' Tr(RnTot)= ',RnTot%trace()
      write(iOut,*)
      write(iOut,*)' Tr(RnOO) = ',RnOO%trace()
      write(iOut,*)' Tr(RnOV) = ',RnOV%trace()
      write(iOut,*)' Tr(RnVO) = ',RnVO%trace()
      write(iOut,*)' Tr(RnVV) = ',RnVV%trace()
      write(iOut,*)
      write(iOut,*)
      write(iOut,*)


      write(iOut,*)
      write(iOut,*)
      write(iOut,*)' Hrant - TRY AGAIN...'
      write(iOut,*)
      write(iOut,*)' ALPHA'
      write(iOut,*)


      RnTot = partitionDensityMatrix(diffDensityAlphaEVecs,  &
        diffDensityAlphaEVals,'relaxation')
      RnOO = MatMul(POcc,MatMul(RnTot,POcc))
      RnOV = MatMul(POcc,MatMul(RnTot,PVirt))
      RnVO = MatMul(PVirt,MatMul(RnTot,POcc))
      RnVV = MatMul(PVirt,MatMul(RnTot,PVirt))
      write(iOut,*)' Tr(RnTot)= ',RnTot%trace()
      write(iOut,*)
      write(iOut,*)' Tr(RnOO) = ',RnOO%trace()
      write(iOut,*)' Tr(RnOV) = ',RnOV%trace()
      write(iOut,*)' Tr(RnVO) = ',RnVO%trace()
      write(iOut,*)' Tr(RnVV) = ',RnVV%trace()
      write(iOut,*)


      XnTot = partitionDensityMatrix(diffDensityAlphaEVecs,  &
        diffDensityAlphaEVals,'excitation')
      XnOO = MatMul(POcc,MatMul(XnTot,POcc))
      XnOV = MatMul(POcc,MatMul(XnTot,PVirt))
      XnVO = MatMul(PVirt,MatMul(XnTot,POcc))
      XnVV = MatMul(PVirt,MatMul(XnTot,PVirt))
      write(iOut,*)' Tr(XnTot)= ',XnTot%trace()
      write(iOut,*)
      write(iOut,*)' Tr(XnOO) = ',XnOO%trace()
      write(iOut,*)' Tr(XnOV) = ',XnOV%trace()
      write(iOut,*)' Tr(XnVO) = ',XnVO%trace()
      write(iOut,*)' Tr(XnVV) = ',XnVV%trace()
      write(iOut,*)

      write(iOut,*)
      write(iOut,*)' BETA'
      write(iOut,*)


      RnTot = partitionDensityMatrix(diffDensityAlphaEVecs,  &
        diffDensityAlphaEVals,'relaxation')
      RnOO = MatMul(POcc,MatMul(RnTot,POcc))
      RnOV = MatMul(POcc,MatMul(RnTot,PVirt))
      RnVO = MatMul(PVirt,MatMul(RnTot,POcc))
      RnVV = MatMul(PVirt,MatMul(RnTot,PVirt))
      write(iOut,*)' Tr(RnTot)= ',RnTot%trace()
      write(iOut,*)
      write(iOut,*)' Tr(RnOO) = ',RnOO%trace()
      write(iOut,*)' Tr(RnOV) = ',RnOV%trace()
      write(iOut,*)' Tr(RnVO) = ',RnVO%trace()
      write(iOut,*)' Tr(RnVV) = ',RnVV%trace()
      write(iOut,*)


      XnTot = partitionDensityMatrix(diffDensityAlphaEVecs,  &
        diffDensityAlphaEVals,'excitation')
      XnOO = MatMul(POcc,MatMul(XnTot,POcc))
      XnOV = MatMul(POcc,MatMul(XnTot,PVirt))
      XnVO = MatMul(PVirt,MatMul(XnTot,POcc))
      XnVV = MatMul(PVirt,MatMul(XnTot,PVirt))
      write(iOut,*)' Tr(XnTot)= ',XnTot%trace()
      write(iOut,*)
      write(iOut,*)' Tr(XnOO) = ',XnOO%trace()
      write(iOut,*)' Tr(XnOV) = ',XnOV%trace()
      write(iOut,*)' Tr(XnVO) = ',XnVO%trace()
      write(iOut,*)' Tr(XnVV) = ',XnVV%trace()
      write(iOut,*)


      goto 999




      PDeltaOO = MatMul(POcc,MatMul(PDeltaMO,POcc))
      PDeltaVV = MatMul(PVirt,MatMul(PDeltaMO,PVirt))
      PDeltaOV = MatMul(POcc,MatMul(PDeltaMO,PVirt))
      PDeltaVO = MatMul(PVirt,MatMul(PDeltaMO,POcc))

      call mqc_print(MatMul(POcc,MatMul(PDeltaMO,POcc)),iOut,header='Delta_oo')
      call mqc_print(MatMul(PVirt,MatMul(PDeltaMO,PVirt)),iOut,header='Delta_vv')

      call mqc_print(MatMul(POcc,MatMul(PDeltaMO,PVirt)),iOut,header='Delta_ov')
      call mqc_print(MatMul(PVirt,MatMul(PDeltaMO,POcc)),iOut,header='Delta_vo')

      write(iOut,*)
      write(iOut,*)' Tr(Delta_oo) = ',PDeltaOO%trace()
      write(iOut,*)' Tr(Delta_vv) = ',PDeltaVV%trace()
      write(iOut,*)' Tr(Delta_ov) = ',PDeltaOV%trace()
      write(iOut,*)' Tr(Delta_vo) = ',PDeltaVO%trace()
      write(iOut,*)
      write(iOut,*)
      write(iOut,*)' Back to the original NIO program.'
      write(iOut,*)
      write(iOut,*)

!hph      goto 999


!hph-

!
!     Form the polestrength (for detachment cases) or the N-1 overlap (for
!     excitation cases). At the end of this block, we decide if this is a
!     detachment (<isNIO>) or excitation job (<isDDNO>).
!
      call determinantOverlap(SMatrixAO,SMatrixAOMinusHalf,  &
        diffDensityAlphaEVals,diffDensityAlphaEVecs,CAlpha2,nElAlpha2,  &
        nBasis,TDOverlapA,nPlusOneAlpha,nMinusOneAlpha,iPlusOneAlpha,  &
        iMinusOneAlpha)
      call determinantOverlap(SMatrixAO,SMatrixAOMinusHalf,  &
        diffDensityBetaEVals,diffDensityBetaEVecs,CBeta2,nElBeta2,  &
        nBasis,TDOverlapB,nPlusOneBeta,nMinusOneBeta,iPlusOneBeta,  &
        iMinusOneBeta)
      isNIO  = (nPlusOneAlpha+nPlusOneBeta).lt.(nMinusOneAlpha+nMinusOneBeta)
      isDDNO = (nPlusOneAlpha+nPlusOneBeta).eq.  &
        (nMinusOneAlpha+nMinusOneBeta)
      write(iOut,3000) nPlusOneAlpha,nMinusOneAlpha,nPlusOneBeta,  &
        nMinusOneBeta
      if(iPrint.ge.1.or.DEBUG) write(iOut,3100) isNIO,isDDNO
      if(.not.(isNIO.xor.isDDNO))  &
        call mqc_error('Confused about NIO or DDNO job.')
      if(isNIO) write(iOut,1500) float(TDOverlapA*TDOverlapB),  &
        float(TDOverlapA),float(TDOverlapB)
!
!     If this is a DDNO job, calculate the overlap of the two determinants.
!
      if(isDDNO) then
        if(nElAlpha1.ne.nElAlpha2) then
          tmpMQCvar3 = 0.0
        else
          tmpMQCvar1 = MatMul(Transpose(CAlpha1%subMatrix(newrange2=[1,nElAlpha1])),  &
            MatMul(SMatrixAO,CAlpha2%subMatrix(newrange2=[1,nElAlpha2])))
          tmpMQCvar3 = tmpMQCvar1%det()
        endIf
        if(nElBeta1.ne.nElBeta2) then
          tmpMQCvar4 = 0.0
        else
          tmpMQCvar2 = MatMul(Transpose(CBeta1%subMatrix(newrange2=[1,nElBeta1])),  &
            MatMul(SMatrixAO,CBeta2%subMatrix(newrange2=[1,nElBeta2])))
          tmpMQCvar4 = tmpMQCvar2%det()
        endIf
        tmpMQCvar = tmpMQCvar3*tmpMQCvar4
        write(iOut,1600) float(tmpMQCvar),float(tmpMQCvar3),float(tmpMQCvar4)
      endIf
!
!     Carry out attachment/detachment density analysis.
!
      Allocate(tmpMatrix2(nBasis,nBasis2),tmpMatrix3(nBasis,nBasis2))
#ifdef PGI
      tmpMatrix2 = dfloat(0)
      tmpMatrix3 = dfloat(0)
#else
      tmpMatrix2 = float(0)
      tmpMatrix3 = float(0)
#endif
      do i = 1,nBasis
        tmpVector = DDNOsAlpha%column(i)
        tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
          float(MQC_Variable_get_MQC(diffDensityAlphaEVals,[i])))
        if(float(MQC_Variable_get_MQC(diffDensityAlphaEVals,[i])).ge.float(0)) then
          tmpMatrix2 = tmpMatrix2 + tmpMatrix1
        else
          tmpMatrix3 = tmpMatrix3 - tmpMatrix1
        endIf
      endDo
      attachmentDensity = tmpMatrix2
      detachmentDensity = tmpMatrix3
      attachmentNumberAlpha = contraction(attachmentDensity,SMatrixAO)
      detachmentNumberAlpha = contraction(detachmentDensity,SMatrixAO)
      tmpMatrix2 = float(0)
      tmpMatrix3 = float(0)
      do i = 1,nBasis
        tmpVector = DDNOsBeta%column(i)
        tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
          float(MQC_Variable_get_MQC(diffDensityBetaEVals,[i])))
        if(float(MQC_Variable_get_MQC(diffDensityBetaEVals,[i])).ge.float(0)) then
          tmpMatrix2 = tmpMatrix2 + tmpMatrix1
        else
          tmpMatrix3 = tmpMatrix3 - tmpMatrix1
        endIf
      endDo
      attachmentDensity = tmpMatrix2
      detachmentDensity = tmpMatrix3
      attachmentNumberBeta = contraction(attachmentDensity,SMatrixAO)
      detachmentNumberBeta = contraction(detachmentDensity,SMatrixAO)
      write(iOut,2000) ' (Original 1995 Model)',float(attachmentNumberAlpha),  &
        float(detachmentNumberAlpha),float(attachmentNumberBeta),  &
        float(detachmentNumberBeta)
!
!     Try promotion number a second time (Modified Model A). Alpha...
!
      tmpMQCvar1 = MatMul(Transpose(CAlpha1%subMatrix(newrange2=[1,nElAlpha1])),  &
        MatMul(SMatrixAO,DDNOsAlpha))
      tmpMQCvar2 = MatMul(CAlpha1%subMatrix(newrange2=[1,nElAlpha1]),tmpMQCvar1)
      tmpMQCvar1 = MatMul(Transpose(CAlpha1%subMatrix(newrange2=[nElAlpha1+1,nBasisUse])),  &
        MatMul(SMatrixAO,DDNOsAlpha))
      tmpMQCvar3 = MatMul(CAlpha1%subMatrix(newrange2=[nElAlpha1+1,nBasisUse]),tmpMQCvar1)
      tmpMatrix2 = float(0)
      tmpMatrix3 = float(0)
      do i = 1,nBasis
        if(float(MQC_Variable_get_MQC(diffDensityAlphaEVals,[i])).ge.float(0)) then
          tmpVector = tmpMQCvar3%column(i)
          tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
            float(MQC_Variable_get_MQC(diffDensityAlphaEVals,[i])))
          tmpMatrix2 = tmpMatrix2 + tmpMatrix1
        else
          tmpVector = tmpMQCvar2%column(i)
          tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
            float(MQC_Variable_get_MQC(diffDensityAlphaEVals,[i])))
          tmpMatrix3 = tmpMatrix3 - tmpMatrix1
        endIf
      endDo
      attachmentDensity = tmpMatrix2
      detachmentDensity = tmpMatrix3
      attachmentNumberAlpha = contraction(attachmentDensity,SMatrixAO)
      detachmentNumberAlpha = contraction(detachmentDensity,SMatrixAO)
!
!     Promotion Number Model A...beta...
!
      tmpMQCvar1 = MatMul(Transpose(CBeta1%subMatrix(newrange2=[1,nElBeta1])),  &
        MatMul(SMatrixAO,DDNOsBeta))
      tmpMQCvar2 = MatMul(CBeta1%subMatrix(newrange2=[1,nElBeta1]),tmpMQCvar1)
      tmpMQCvar1 = MatMul(Transpose(CBeta1%subMatrix(newrange2=[nElBeta1+1,nBasisUse])),  &
        MatMul(SMatrixAO,DDNOsBeta))
      tmpMQCvar3 = MatMul(CBeta1%subMatrix(newrange2=[nElBeta1+1,nBasisUse]),tmpMQCvar1)
      tmpMatrix2 = float(0)
      tmpMatrix3 = float(0)
      do i = 1,nBasis
        if(float(MQC_Variable_get_MQC(diffDensityBetaEVals,[i])).ge.float(0)) then
          tmpVector = tmpMQCvar3%column(i)
          tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
            float(MQC_Variable_get_MQC(diffDensityBetaEVals,[i])))
          tmpMatrix2 = tmpMatrix2 + tmpMatrix1
        else
          tmpVector = tmpMQCvar2%column(i)
          tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
            float(MQC_Variable_get_MQC(diffDensityBetaEVals,[i])))
          tmpMatrix3 = tmpMatrix3 - tmpMatrix1
        endIf
      endDo
      attachmentDensity = tmpMatrix2
      detachmentDensity = tmpMatrix3
      attachmentNumberBeta = contraction(attachmentDensity,SMatrixAO)
      detachmentNumberBeta = contraction(detachmentDensity,SMatrixAO)
      write(iOut,2000) ' (Modified Model A)',float(attachmentNumberAlpha),  &
        float(detachmentNumberAlpha),float(attachmentNumberBeta),  &
        float(detachmentNumberBeta)
!
!     Try promotion number a third time...alpha.
!

      tmpMQCvar1 = MatMul(Transpose(CAlpha1%subMatrix(newrange2=[1,nElAlpha1])),  &
        MatMul(SMatrixAO,DDNOsAlpha))
      tmpMQCvar3 = MatMul(CAlpha1%subMatrix(newrange2=[1,nElAlpha1]),tmpMQCvar1)
      tmpMQCvar1 = MatMul(Transpose(CAlpha1%subMatrix(newrange2=[nElAlpha1+1,nBasisUse])),  &
        MatMul(SMatrixAO,DDNOsAlpha))
      tmpMQCvar4 = MatMul(CAlpha1%subMatrix(newrange2=[nElAlpha1+1,nBasisUse]),tmpMQCvar1)
      tmpMatrix2 = float(0)
      tmpMatrix3 = float(0)
      do i = 1,nBasis
        tmpVector = tmpMQCvar4%column(i)
        tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
          float(MQC_Variable_get_MQC(diffDensityAlphaEVals,[i])))
        tmpMatrix2 = tmpMatrix2 + tmpMatrix1
!
        tmpVector = tmpMQCvar3%column(i)
        tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
          float(MQC_Variable_get_MQC(diffDensityAlphaEVals,[i])))
        tmpMatrix3 = tmpMatrix3 - tmpMatrix1
      endDo
      attachmentDensity = tmpMatrix2
      detachmentDensity = tmpMatrix3
      attachmentNumberAlpha = contraction(attachmentDensity,SMatrixAO)
      detachmentNumberAlpha = contraction(detachmentDensity,SMatrixAO)
!
!     Try promotion number a third time...beta.
!
      tmpMQCvar1 = MatMul(Transpose(CBeta1%subMatrix(newrange2=[1,nElBeta1])),  &
        MatMul(SMatrixAO,DDNOsBeta))
      tmpMQCvar3 = MatMul(CBeta1%subMatrix(newrange2=[1,nElBeta1]),tmpMQCvar1)
      tmpMQCvar1 = MatMul(Transpose(CBeta1%subMatrix(newrange2=[nElBeta1+1,nBasisUse])),  &
        MatMul(SMatrixAO,DDNOsBeta))
      tmpMQCvar4 = MatMul(CBeta1%subMatrix(newrange2=[nElBeta1+1,nBasisUse]),tmpMQCvar1)
      tmpMatrix2 = float(0)
      tmpMatrix3 = float(0)
      do i = 1,nBasis
        tmpVector = tmpMQCvar4%column(i)
        tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
          float(MQC_Variable_get_MQC(diffDensityBetaEVals,[i])))
        tmpMatrix2 = tmpMatrix2 + tmpMatrix1
!
        tmpVector = tmpMQCvar3%column(i)
        tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
          float(MQC_Variable_get_MQC(diffDensityBetaEVals,[i])))
        tmpMatrix3 = tmpMatrix3 - tmpMatrix1
      endDo
      attachmentDensity = tmpMatrix2
      detachmentDensity = tmpMatrix3
      attachmentNumberBeta = contraction(attachmentDensity,SMatrixAO)
      detachmentNumberBeta = contraction(detachmentDensity,SMatrixAO)
      write(iOut,2000) ' (Modified Model B)',float(attachmentNumberAlpha),  &
        float(detachmentNumberAlpha),float(attachmentNumberBeta),  &
        float(detachmentNumberBeta)
!
!     Compute the transition dipole and dipole strength for DDNO jobs.
!
      if(isDDNO) then
        call GMatrixFile1%getArray('Dipole Integrals',  &
          mqcVarOut=dipoleAOx,arraynum=1)
        call GMatrixFile1%getArray('Dipole Integrals',  &
          mqcVarOut=dipoleAOy,arraynum=2)
        call GMatrixFile1%getArray('Dipole Integrals',  &
          mqcVarOut=dipoleAOz,arraynum=3)
        if(DEBUG) then
          call mqc_print(contraction(PMatrixTotal1,dipoleAOx),header='P1(total).dipoleX')
          call mqc_print(contraction(PMatrixTotal1,dipoleAOy),header='P1(total).dipoleY')
          call mqc_print(contraction(PMatrixTotal1,dipoleAOz),header='P1(total).dipoleZ')
          call mqc_print(contraction(PMatrixTotal2,dipoleAOx),header='P2(total).dipoleX')
          call mqc_print(contraction(PMatrixTotal2,dipoleAOy),header='P2(total).dipoleY')
          call mqc_print(contraction(PMatrixTotal2,dipoleAOz),header='P2(total).dipoleZ')
        endIf
        if(iPlusOneAlpha.gt.0) then
          pDDNO = DDNOsAlpha%column(iPlusOneAlpha)
        elseIf(iPlusOneBeta.gt.0) then
          pDDNO = DDNOsBeta%column(iPlusOneBeta)
        else
          call mqc_error('No particle DDNO located.')
        endIf
        if(iMinusOneAlpha.gt.0) then
          hDDNO = DDNOsAlpha%column(iMinusOneAlpha)
        elseIf(iMinusOneBeta.gt.0) then
          hDDNO = DDNOsBeta%column(iMinusOneBeta)
        else
          call mqc_error('No hole DDNO located.')
        endIf
        if(iPrint.ge.1.or.DEBUG) then
          call pDDNO%print(header='particle DDNO')
          call hDDNO%print(header='hole DDNO')
        endIf
        if(DEBUG) then
          call mqc_print_scalar_real(float(dot_product(pDDNO,MQC_Variable_MatrixVector(SMatrixAO,pDDNO))),6,header='<p|p>')
          call mqc_print_scalar_real(float(dot_product(hDDNO,MQC_Variable_MatrixVector(SMatrixAO,hDDNO))),6,header='<h|h>')
          call mqc_print_scalar_real(float(dot_product(hDDNO,MQC_Variable_MatrixVector(SMatrixAO,pDDNO))),6,header='<h|p>')
          call mqc_print_scalar_real(float(dot_product(pDDNO,MQC_Variable_MatrixVector(SMatrixAO,hDDNO))),6,header='<p|h>')
        endIf
        if(DEBUG) then
          tmpMQCvar = MQC_Variable_MatrixVector(dipoleAOy,hDDNO)
          call tmpMQCvar%print(header='dipoleAOy.hDDNO')
          tmpMQCvar = MQC_Variable_MatrixVector(dipoleAOy,pDDNO)
          call tmpMQCvar%print(header='dipoleAOy.pDDNO')
          tmpMQCvar = dot_product(pDDNO,MQC_Variable_MatrixVector(dipoleAOy,hDDNO))
          call tmpMQCvar%print(header='same mu^y')
          tmpMQCvar = dot_product(hDDNO,MQC_Variable_MatrixVector(dipoleAOy,pDDNO))
          call tmpMQCvar%print(header='flipped mu^y')
        endIf
        transitionDipole(1) =  dot_product(pDDNO,MQC_Variable_MatrixVector(dipoleAOx,hDDNO))
        transitionDipole(2) =  dot_product(pDDNO,MQC_Variable_MatrixVector(dipoleAOy,hDDNO))
        transitionDipole(3) =  dot_product(pDDNO,MQC_Variable_MatrixVector(dipoleAOz,hDDNO))
        call mqc_print(transitionDipole,6,header='Transition Dipole Moment',blank_at_top=.true.)
        TDparticleHoleMag = dot_product(transitionDipole,transitionDipole)
        call TDparticleHoleMag%print(header='TDparticleHoleMag = ')
        if(DEBUG) then
          call TDparticleHoleMag%print(header='Transition Dipole contribution to the Dipole Strength =')
          dipoleStrength = TDOverlapA*TDOverlapB*TDparticleHoleMag
          call dipoleStrength%print(header='OLD Dipole Strength (au) =')
          dipoleStrength = TDOverlapA*TDOverlapA*TDOverlapB*TDOverlapB*TDparticleHoleMag
          call dipoleStrength%print(header='NEW Dipole Strength (au) =')
        else
          dipoleStrength = TDOverlapA*TDOverlapA*TDOverlapB*TDOverlapB*TDparticleHoleMag
          call dipoleStrength%print(header='Dipole Strength (au) =')
        endIf
        temp = TDOverlapA*TDOverlapA*TDOverlapB*TDOverlapB
        temp = float(1) - temp
        write(iOut,4000) temp
        write(iOut,4010) float(TDOverlapA*TDOverlapA),float(TDOverlapB*TDOverlapB)
        oscillatorStrength = float(2)*deltaSCFEnergy/float(3)
        oscillatorStrength = oscillatorStrength*dipoleStrength
        call dipoleStrength%print(header='dipoleStrength')
        call oscillatorStrength%print(header='Oscillator Strength  =')
      endIf

      if(isNIO.or..not.doTestCode) goto 998

!hph+
      goto 999
!hph-


!hph+
      write(*,*)
      write(*,*)
      write(*,*)' Hrant - Calling categorization routine...'
      call categorizeDDNOs(diffDensityAlphaEVals,infoDDNOsAlpha,nRelaxationDDNOsAlpha)
      call categorizeDDNOs(diffDensityBetaEVals,infoDDNOsBeta,nRelaxationDDNOsBeta)
      write(*,*)' Hrant - Back from categorization routine!'
      write(*,*)
      write(*,*)
      write(*,*)' Hrant - Calling projectDDNOs for ALPHA spin...'
      call projectDDNOs(infoDDNOsAlpha,DDNOsAlpha,SMatrixAO,CAlpha1,  &
        nElAlpha1,nRelaxationDDNOsAlpha,diffDensityAlpha)
      write(*,*)
      write(*,*)
      write(*,*)' Hrant - Calling projectDDNOs for BETA  spin...'
      call projectDDNOs(infoDDNOsBeta,DDNOsBeta,SMatrixAO,CBeta1,  &
        nElBeta1,nRelaxationDDNOsBeta,diffDensityBeta)
      write(*,*)

      goto 999

      write(*,*)
      write(*,*)' Hrant - Calling projectDDNOs for BETA spin...'
      call projectDDNOs(infoDDNOsBeta,DDNOsBeta,SMatrixAO,CBeta1,  &
        nElBeta1,nRelaxationDDNOsBeta,diffDensityBeta)
      write(*,*)
      write(*,*)' Hrant - Back from calling projectDDNOs!'
      write(*,*)
      write(*,*)
      
!hph-

  998 Continue

!
!     If requested, write results to an output Gaussian matrix file.
!
      if(doMatrixFileOut) then
        write(iOut,7000) TRIM(matrixFilenameOut)
        GMatrixFileOut = GMatrixFile1
        call GMatrixFileOut%create(TRIM(matrixFilenameOut))
!
!       Basis set info...
        call GMatrixFile1%getArray('SHELL TO ATOM MAP',mqcVarOut=tmpMQCvar)
        call GMatrixFileOut%writeArray2('SHELL TO ATOM MAP',tmpMQCvar)
        call GMatrixFile1%getArray('SHELL TYPES',mqcVarOut=tmpMQCvar)
        call GMatrixFileOut%writeArray2('SHELL TYPES',tmpMQCvar)
        call GMatrixFile1%getArray('NUMBER OF PRIMITIVES PER SHELL',mqcVarOut=tmpMQCvar)
        call GMatrixFileOut%writeArray2('NUMBER OF PRIMITIVES PER SHELL',tmpMQCvar)
        call GMatrixFile1%getArray('PRIMITIVE EXPONENTS',mqcVarOut=tmpMQCvar)
        call GMatrixFileOut%writeArray2('PRIMITIVE EXPONENTS',tmpMQCvar)
        call GMatrixFile1%getArray('CONTRACTION COEFFICIENTS',mqcVarOut=tmpMQCvar)
        call GMatrixFileOut%writeArray2('CONTRACTION COEFFICIENTS',tmpMQCvar)
        call GMatrixFile1%getArray('P(S=P) CONTRACTION COEFFICIENTS',mqcVarOut=tmpMQCvar)
        call GMatrixFileOut%writeArray2('P(S=P) CONTRACTION COEFFICIENTS',tmpMQCvar)
        call GMatrixFile1%getArray('COORDINATES OF EACH SHELL',mqcVarOut=tmpMQCvar)
        call GMatrixFileOut%writeArray2('COORDINATES OF EACH SHELL',tmpMQCvar)
!
!       DDNO eigenvectors and eigenvalues...
        call GMatrixFileOut%writeArray2('ALPHA ORBITAL ENERGIES',diffDensityAlphaEVals)
        call GMatrixFileOut%writeArray2('BETA ORBITAL ENERGIES',diffDensityBetaEVals)
        call GMatrixFileOut%writeArray2('ALPHA MO COEFFICIENTS',DDNOsAlpha)
        call GMatrixFileOut%writeArray2('BETA MO COEFFICIENTS',DDNOsBeta)
!
!       Close out the matrix file.
        call GMatrixFileOut%closeFile()
      endIf
!
  999 Continue
      write(iOut,8999)
      end program nio
