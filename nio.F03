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
      integer(kind=int64)::nCommands,iPrint=0,i,nAtoms,nAtoms2,  &
        nBasis,nBasis2,nBasisUse,nBasisUse2,nEl1,nEl2,nElAlpha1,  &
        nElBeta1,nElAlpha2,NElBeta2,nPlusOneAlpha,nMinusOneAlpha,  &
        iPlusOneAlpha,iMinusOneAlpha,nPlusOneBeta,nMinusOneBeta,  &
        iPlusOneBeta,iMinusOneBeta,nPlusOne,nMinusOne,  &
        nRelaxationDDNOsAlpha,nRelaxationDDNOsBeta
      integer(kind=int64),allocatable,dimension(:)::tmpVectorInt
      real(kind=real64)::scfEnergy1,scfEnergy2,deltaSCFEnergy,  &
        totalEnergy1,totalEnergy2,deltaTotalEnergy,temp
      real(kind=real64),dimension(3)::transitionDipole,particleXYZ(3),  &
        holeXYZ(3)
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
        PMatrixAlpha2,PMatrixBeta2,PMatrixTotal2,PMatrixCorrAlpha1,  &
        PMatrixCorrBeta1,PMatrixCorrTotal1,PMatrixCorrAlpha2,  &
        PMatrixCorrBeta2,PMatrixCorrTotal2,PMatrixAlpha1MO,  &
        PMatrixBeta1MO,PMatrixAlpha2MO,PMatrixBeta2MO,diffDensityAlpha,  &
        diffDensityBeta,diffDensityAlphaEVecs,diffDensityAlphaEVals,  &
        diffDensityBetaEVecs,diffDensityBetaEVals,diffDensityAlphaMO,  &
        diffDensityBetaMO,diffDensityAlphaEVecsMO,diffDensityAlphaEValsMO,  &
        diffDensityBetaEVecsMO,diffDensityBetaEValsMO,attachmentDensity,  &
        detachmentDensity,attachmentNumberAlpha,detachmentNumberAlpha,  &
        attachmentNumberBeta,detachmentNumberBeta,excitationMatrixAlpha,  &
        relaxationMatrixAlpha,excitationMatrixBeta,relaxationMatrixBeta,  &
        excitationNumberAlphaO2V,excitationNumberAlphaV2O,  &
        excitationNumberBetaO2V,excitationNumberBetaV2O,  &
        relaxationNumberAlphaO2V,relaxationNumberAlphaV2O,  &
        relaxationNumberBetaO2V,relaxationNumberBetaV2O
      type(MQC_Variable)::CAlpha1,CBeta1,CAlpha2,CBeta2,TAlpha,TBeta
      type(MQC_Variable)::dipoleAOx,dipoleAOy,dipoleAOz
      type(MQC_Variable)::tmpMQCvar,tmpMQCvar1,tmpMQCvar2,tmpMQCvar3,  &
        tmpMQCvar4
      logical::doTestCode=.false.,doMatrixFileOut=.false.,isNIO,isDDNO,  &
        found,postSCF
!hph
      type(MQC_Variable)::POccAlpha1,POccBeta1,PDeltaMO,PDeltaOO,PDeltaVV,  &
        PDeltaOV,PDeltaVO
!
!     Format Statements
!
 1000 Format(1x,'Enter Program NIO'/,  &
             1x,'Version 2025.09.10.00')
 1010 Format(1x,'Matrix File 1:      ',A,/,  &
             1x,'Matrix File 2:      ',A,/,  &
             1x,'Output Matrix File: ',A,/)
 1020 Format(1x,'Matrix File 1:      ',A,/,  &
             1x,'Matrix File 2:      ',A,/)
 1100 Format(1x,'nAtoms=',I4,3x,'nBasis   =',I4,3x,'nBasisUse=',I4,/,  &
             1x,'nEl1  =',I4,3x,'nElAlpha1=',I4,3x,'nElBeta  =',I4,/,  &
             1x,'nEl2  =',I4,3x,'nElAlpha2=',I4,3x,'nElBeta  =',I4,/)
 1200 Format(/,1x,A,' Energies',/,  &
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
        1x,'Citation: ',A,/,  &
        3x,'ALPHA: Attachement Number: ',F6.3,/,  &
        10x,'Detachment  Number: ',F6.3,/,  &
        3x,'BETA : Attachement Number: ',F6.3,/,  &
        10x,'Detachment  Number: ',F6.3)
 2100 Format(/,1x,'Excitation/Relaxation Density Analysis',/,  &
        3x,'ALPHA: Excitation Population o->v: ',F6.3,/,  &
        3x,'       Excitation Population v->o: ',F6.3,/,  &
        3x,'BETA : Excitation Population o->v: ',F6.3,/,  &
        3x,'       Excitation Population v->o: ',F6.3,/,  &
        3x,'TOTAL: Excitation Population:      ',F6.3,/,/,  &
        3x,'ALPHA: Relaxation Population o->v: ',F6.3,/,  &
        3x,'       Relaxation Population v->o: ',F6.3,/,  &
        3x,'BETA : Relaxation Population o->v: ',F6.3,/,  &
        3x,'       Relaxation Population v->o: ',F6.3,/,  &
        3x,'TOTAL: Relaxation Population:      ',F6.3,/)
 3000 Format(/,1x,'DDNO/NIO Excitation/Detachment Number',/,  &
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
      if(.not.mqc_version_check(newerThanMajor=25,newerThanMinor=9,  &
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
!     Load the atomic orbital overlap matrix and form S^(1/2) and S^(-1/2).
!
      SMatrixAO = getOverlapMatrix(GMatrixFile1)
      call SMatrixAO%eigen(SMatrixEVals,SMatrixEVecs)
      if(DEBUG) then
        call SMatrixAO%print(header='Overlap Matrix')
        call SMatrixEVals%print(header='S matrix eigen-values:')
        call mqc_print(MatMul(Transpose(SMatrixEVecs),SMatrixEVecs),  &
          header='SEVecs(t).SEVecs')
        call mqc_print(MatMul(MatMul(SMatrixEVecs,SMatrixEVals%diag()),  &
          TRANSPOSE(SMatrixEVecs)),6,'U.lambda.Ut')
      endIf
      tmpMQCvar = SMatrixEVals%rpower(0.5)
      SMAtrixAOHalf = MatMul(MatMul(SMatrixEVecs,tmpMQCvar%diag()),  &
        TRANSPOSE(SMatrixEVecs))
      tmpMQCvar = SMatrixEVals%rpower(-0.5)
      SMAtrixAOMinusHalf = MatMul(MatMul(SMatrixEVecs,tmpMQCvar%diag()),  &
        TRANSPOSE(SMatrixEVecs))
      if(DEBUG) then
        call SMAtrixAOHalf%print(header='S**(1/2)')
        call SMatrixAOMinusHalf%print(header='S**(-1/2)')
      endIf
!
!     Load the SCF densities and MO coefficient matrices (in the AO basis). The
!     code below treats all systems as open shell, so closed shell results are
!     handled by copying the density matrix from restricted calculations into
!     alpha and beta density matrix arrays.
!

!hph+
!      call GMatrixFile1%getArray('ALPHA MP2 DENSITY MATRIX',mqcVarOut=PMatrixAlpha1)
!      call GMatrixFile1%getArray('ALPHA SCF DENSITY MATRIX',mqcVarOut=tmpMQCvar)
!      PMatrixAlpha1 = PMatrixAlpha1-tmpMQCvar
!      call GMatrixFile1%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=CAlpha1)
!      if(GMatrixFile1%isUnrestricted()) then
!        call GMatrixFile1%getArray('BETA MP2 DENSITY MATRIX',mqcVarOut=PMatrixBeta1)
!        call GMatrixFile1%getArray('BETA SCF DENSITY MATRIX',mqcVarOut=tmpMQCvar)
!        PMatrixBeta1 = PMatrixBeta1-tmpMQCvar
!        call GMatrixFile1%getArray('BETA MO COEFFICIENTS',mqcVarOut=CBeta1)
!      else
!        PMatrixBeta1  = PMatrixAlpha1
!        CBeta1 = CAlpha1
!      endIf
!      PMatrixTotal1 = PMatrixAlpha1+PMatrixBeta1
!      call GMatrixFile2%getArray('ALPHA MP2 DENSITY MATRIX',mqcVarOut=PMatrixAlpha2)
!      call GMatrixFile2%getArray('ALPHA SCF DENSITY MATRIX',mqcVarOut=tmpMQCvar)
!      PMatrixAlpha2 = PMatrixAlpha2-tmpMQCvar
!      call GMatrixFile2%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=CAlpha2)
!      if(GMatrixFile2%isUnrestricted()) then
!        call GMatrixFile2%getArray('BETA MP2 DENSITY MATRIX',mqcVarOut=PMatrixBeta2)
!        call GMatrixFile2%getArray('BETA SCF DENSITY MATRIX',mqcVarOut=tmpMQCvar)
!        PMatrixBeta2 = PMatrixBeta2-tmpMQCvar
!        call GMatrixFile2%getArray('BETA MO COEFFICIENTS',mqcVarOut=CBeta2)
!      else
!        PMatrixBeta2  = PMatrixAlpha2
!        CBeta2 = CAlpha2
!      endIf
!      PMatrixTotal2 = PMatrixAlpha2+PMatrixBeta2
!hph-

      call GMatrixFile1%getArray('ALPHA SCF DENSITY MATRIX',mqcVarOut=PMatrixAlpha1)
      call GMatrixFile1%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=CAlpha1)
      if(GMatrixFile1%isUnrestricted()) then
        call GMatrixFile1%getArray('BETA SCF DENSITY MATRIX',mqcVarOut=PMatrixBeta1)
        call GMatrixFile1%getArray('BETA MO COEFFICIENTS',mqcVarOut=CBeta1)
      else
        PMatrixBeta1  = PMatrixAlpha1
        CBeta1 = CAlpha1
      endIf
      PMatrixTotal1 = PMatrixAlpha1+PMatrixBeta1
      call GMatrixFile2%getArray('ALPHA SCF DENSITY MATRIX',mqcVarOut=PMatrixAlpha2)
      call GMatrixFile2%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=CAlpha2)
      if(GMatrixFile2%isUnrestricted()) then
        call GMatrixFile2%getArray('BETA SCF DENSITY MATRIX',mqcVarOut=PMatrixBeta2)
        call GMatrixFile2%getArray('BETA MO COEFFICIENTS',mqcVarOut=CBeta2)
      else
        PMatrixBeta2  = PMatrixAlpha2
        CBeta2 = CAlpha2
      endIf
      PMatrixTotal2 = PMatrixAlpha2+PMatrixBeta2
!
!     Look for the total density and see if we think we're looking at a Post-SCF
!     job. If we are, then we fill in the correlated density matrices.
!
      call GMatrixFile1%getArray('ALPHA DENSITY MATRIX',  &
        mqcVarOut=PMatrixCorrAlpha1,foundOut=found)
      postSCF = found

      write(iOut,*)' Hrant - postSCF = ',postSCF
      if(postSCF) then
        call GMatrixFile2%getArray('ALPHA DENSITY MATRIX',  &
          mqcVarOut=PMatrixCorrAlpha2,foundOut=found)
        postSCF = found
      endIf
      if(postSCF) then
        if(GMatrixFile1%isUnrestricted()) then
          call GMatrixFile1%getArray('BETA DENSITY MATRIX',mqcVarOut=PMatrixCorrBeta1)
        else
          PMatrixCorrBeta1  = PMatrixCorrAlpha1
        endIf
        if(GMatrixFile2%isUnrestricted()) then
          call GMatrixFile2%getArray('BETA DENSITY MATRIX',mqcVarOut=PMatrixCorrBeta2)
        else
          PMatrixCorrBeta2  = PMatrixCorrAlpha2
        endIf
        PMatrixCorrTotal1 = PMatrixCorrAlpha1+PMatrixCorrBeta1
        PMatrixCorrTotal2 = PMatrixCorrAlpha2+PMatrixCorrBeta2
      endIf



!hph Hrant - I left off here on building up post-SCF logic and data collection.





!
!     Pull the energies from the two jobs and report the delta-E value in a few
!     units.
!
      scfEnergy1 = GMatrixFile1%getValReal('scfEnergy')
      scfEnergy2 = GMatrixFile2%getValReal('scfEnergy')
      deltaSCFEnergy = scfEnergy2-scfEnergy1
      write(iOut,1200) 'SCF',scfEnergy1,scfEnergy2,deltaSCFEnergy,  &
        deltaSCFEnergy*evPHartree,deltaSCFEnergy*cmM1PHartree,  &
        mqc_float(10000000)*(planck*sLight)/(deltaSCFEnergy*jPHartree)
      totalEnergy1 = GMatrixFile1%getValReal('energy')
      totalEnergy2 = GMatrixFile2%getValReal('energy')
      deltaTotalEnergy = totalEnergy2-totalEnergy1
      write(iOut,1200) 'Total',totalEnergy1,totalEnergy2,  &
        deltaTotalEnergy,deltaTotalEnergy*evPHartree,  &
        deltaTotalEnergy*cmM1PHartree,  &
        mqc_float(10000000)*(planck*sLight)/(deltaTotalEnergy*jPHartree)
!
!     Form the difference density and construct the DDNOs, which are NIOs in
!     electron detachment cases.
!
      diffDensityAlpha = PMatrixAlpha2-PMatrixAlpha1
      diffDensityBeta  = PMatrixBeta2-PMatrixBeta1
      if(iPrint.ge.1.or.DEBUG) then
        call mqc_print(contraction(diffDensityAlpha,SMatrixAO),  &
          header='DP(alpha).S = ')
        call mqc_print(contraction(diffDensityBeta,SMatrixAO),  &
          header='DP(beta).S  = ')
      endIf
      tmpMQCvar = MatMul(SMatrixAOHalf,MatMul(diffDensityAlpha,  &
        SMatrixAOHalf))
      if(DEBUG) call tmpMQCvar%print(iOut,header='alpha S^{1/2}PS^{1/2}')
      call tmpMQCvar%eigen(diffDensityAlphaEVals,diffDensityAlphaEVecs)
      tmpMQCvar = MatMul(SMatrixAOHalf,MatMul(diffDensityBeta,  &
        SMatrixAOHalf))
      if(DEBUG) call tmpMQCvar%print(iOut,header='beta  S^{1/2}PS^{1/2}')
      call tmpMQCvar%eigen(diffDensityBetaEVals,diffDensityBetaEVecs)
!
!     Order the AO difference density eigenvectors and eigenvalues so that the
!     largest change occupation pair is at the HOMO-LUMO gap.
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
      if(iPrint.ge.1.or.DEBUG) then
        call diffDensityAlphaEVals%print(header='Alpha Occupation Change Values')
        call diffDensityBetaEVals%print(header='Beta Occupation Change Value')
      endIf
!
!     Back transform the difference density natural orbitals to the AO basis.
!
      DDNOsAlpha = MatMul(SMatrixAOMinusHalf,diffDensityAlphaEVecs)
      DDNOsBeta  = MatMul(SMatrixAOMinusHalf,diffDensityBetaEVecs)
!
!     Form the densities, difference densities, and DDNOs in the initial state
!     MO basis.
!
      tmpMQCvar  = MatMul(Transpose(CAlpha1),MatMul(SMatrixAO,  &
        PMatrixAlpha1))
      PMatrixAlpha1MO = MatMul(tmpMQCvar,MatMul(SMatrixAO,CAlpha1))
      tmpMQCvar  = MatMul(Transpose(CBeta1),MatMul(SMatrixAO,  &
        PMatrixBeta1))
      PMatrixBeta1MO = MatMul(tmpMQCvar,MatMul(SMatrixAO,CBeta1))
      tmpMQCvar  = MatMul(Transpose(CAlpha1),MatMul(SMatrixAO,diffDensityAlpha))
      diffDensityAlphaMO = MatMul(tmpMQCvar,MatMul(SMatrixAO,CAlpha1))
      call diffDensityAlphaMO%trimZero()
      tmpMQCvar  = MatMul(Transpose(CBeta1),MatMul(SMatrixAO,diffDensityBeta))
      diffDensityBetaMO = MatMul(tmpMQCvar,MatMul(SMatrixAO,CBeta1))
      call diffDensityBetaMO%trimZero()
      call diffDensityAlphaMO%eigen(diffDensityAlphaEValsMO,  &
        diffDensityAlphaEVecsMO)
      call diffDensityBetaMO%eigen(diffDensityBetaEValsMO,  &
        diffDensityBetaEVecsMO)
      if(iPrint.ge.1.or.DEBUG) then
        call diffDensityAlphaEValsMO%print(  &
          header='Alpha Occupation Change Values, MO')
        call diffDensityBetaEValsMO%print(  &
          header='Beta Occupation Change Value, MO')
      endIf
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
      if(isNIO.eqv.isDDNO)  &
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
        call tmpMQCvar%trimZero()
        call tmpMQCvar3%trimZero()
        call tmpMQCvar4%trimZero()
        write(iOut,1600) float(tmpMQCvar),float(tmpMQCvar3),float(tmpMQCvar4)
      endIf
!
!     Carry out attachment/detachment density analysis.
!
      tmpMQCvar = partitionDensityMatrix(diffDensityAlphaEVecsMO,  &
        diffDensityAlphaEValsMO,'attachment')
      attachmentNumberAlpha = tmpMQCvar%trace()
      tmpMQCvar = partitionDensityMatrix(diffDensityAlphaEVecsMO,  &
        diffDensityAlphaEValsMO,'detachment')
      detachmentNumberAlpha = tmpMQCvar%trace()
      tmpMQCvar = partitionDensityMatrix(diffDensityBetaEVecsMO,  &
        diffDensityBetaEValsMO,'attachment')
      attachmentNumberBeta = tmpMQCvar%trace()
      tmpMQCvar = partitionDensityMatrix(diffDensityBetaEVecsMO,  &
        diffDensityBetaEValsMO,'detachment')
      detachmentNumberBeta = tmpMQCvar%trace()
      write(iOut,2000) ' (Original 1995 Model)',  &
        'Head-Gordon, Grana, Maurice, and White, J. Phys. Chem. 99, 14261-14270 (1995).',  &
        float(attachmentNumberAlpha),  &
        float(detachmentNumberAlpha),float(attachmentNumberBeta),  &
        float(detachmentNumberBeta)
!
!     Evaluate the excitation number.
!
      call MO_partition(diffDensityAlphaMO,PMatrixAlpha1MO,  &
        traceXvv=attachmentNumberAlpha)
      call MO_partition(diffDensityAlphaMO,PMatrixAlpha1MO,  &
        traceXoo=detachmentNumberAlpha)
      call MO_partition(diffDensityBetaMO,PMatrixBeta1MO,  &
        traceXvv=attachmentNumberBeta)
      call MO_partition(diffDensityBetaMO,PMatrixBeta1MO,  &
        traceXoo=detachmentNumberBeta)
      write(iOut,2000) ' (Excitation Number Model)',  &
        'Barca, Gilbert, and Gill, J. Chem. Theory Comput. 14, 9-13 (2018).',  &
        float(attachmentNumberAlpha),  &
        -float(detachmentNumberAlpha),float(attachmentNumberBeta),  &
        -float(detachmentNumberBeta)
!
!     Evaluate the Excitation/Relaxation model.
!
      tmpMQCvar = partitionDensityMatrix(diffDensityAlphaEVecsMO,  &
        diffDensityAlphaEValsMO,'excitation')
      excitationMatrixAlpha = tmpMQCvar
      tmpMQCvar = partitionDensityMatrix(diffDensityAlphaEVecsMO,  &
        diffDensityAlphaEValsMO,'relaxation')
      relaxationMatrixAlpha = tmpMQCvar
      tmpMQCvar = partitionDensityMatrix(diffDensityBetaEVecsMO,  &
        diffDensityBetaEValsMO,'excitation')
      excitationMatrixBeta = tmpMQCvar
      tmpMQCvar = partitionDensityMatrix(diffDensityBetaEVecsMO,  &
        diffDensityBetaEValsMO,'relaxation')
      relaxationMatrixBeta = tmpMQCvar
      call MO_partition(excitationMatrixAlpha,PMatrixAlpha1MO,  &
        traceXoo=excitationNumberAlphaO2V,  &
        traceXvv=excitationNumberAlphaV2O)
      call MO_partition(excitationMatrixBeta,PMatrixBeta1MO,  &
        traceXoo=excitationNumberBetaO2V,  &
        traceXvv=excitationNumberBetaV2O)
      call MO_partition(relaxationMatrixAlpha,PMatrixAlpha1MO,  &
        traceXoo=relaxationNumberAlphaO2V,  &
        traceXvv=relaxationNumberAlphaV2O)
      call MO_partition(relaxationMatrixBeta,PMatrixBeta1MO,  &
        traceXoo=relaxationNumberBetaO2V,  &
        traceXvv=relaxationNumberBetaV2O)
      write(iOut,2100) float(excitationNumberAlphaO2V),  &
        float(excitationNumberAlphaV2O),  &
        float(excitationNumberBetaO2V),  &
        float(excitationNumberBetaV2O),  &
        float(excitationNumberAlphaV2O)+float(excitationNumberBetaV2O),  &
        float(relaxationNumberAlphaO2V),  &
        float(relaxationNumberAlphaV2O),  &
        float(relaxationNumberBetaO2V),  &
        float(relaxationNumberBetaV2O),  &
        float(relaxationNumberAlphaV2O)+float(relaxationNumberBetaV2O)

!hph+
!!
!!     Try promotion number a second time (Modified Model A). Alpha...
!!
!      tmpMQCvar1 = MatMul(Transpose(CAlpha1%subMatrix(newrange2=[1,nElAlpha1])),  &
!        MatMul(SMatrixAO,DDNOsAlpha))
!      tmpMQCvar2 = MatMul(CAlpha1%subMatrix(newrange2=[1,nElAlpha1]),tmpMQCvar1)
!      tmpMQCvar1 = MatMul(Transpose(CAlpha1%subMatrix(newrange2=[nElAlpha1+1,nBasisUse])),  &
!        MatMul(SMatrixAO,DDNOsAlpha))
!      tmpMQCvar3 = MatMul(CAlpha1%subMatrix(newrange2=[nElAlpha1+1,nBasisUse]),tmpMQCvar1)
!      tmpMatrix2 = float(0)
!      tmpMatrix3 = float(0)
!      do i = 1,nBasis
!        if(float(MQC_Variable_get_MQC(diffDensityAlphaEVals,[i])).ge.float(0)) then
!          tmpVector = tmpMQCvar3%column(i)
!          tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
!            float(MQC_Variable_get_MQC(diffDensityAlphaEVals,[i])))
!          tmpMatrix2 = tmpMatrix2 + tmpMatrix1
!        else
!          tmpVector = tmpMQCvar2%column(i)
!          tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
!            float(MQC_Variable_get_MQC(diffDensityAlphaEVals,[i])))
!          tmpMatrix3 = tmpMatrix3 - tmpMatrix1
!        endIf
!      endDo
!      attachmentDensity = tmpMatrix2
!      detachmentDensity = tmpMatrix3
!      attachmentNumberAlpha = contraction(attachmentDensity,SMatrixAO)
!      detachmentNumberAlpha = contraction(detachmentDensity,SMatrixAO)
!!
!!     Promotion Number Model A...beta...
!!
!      tmpMQCvar1 = MatMul(Transpose(CBeta1%subMatrix(newrange2=[1,nElBeta1])),  &
!        MatMul(SMatrixAO,DDNOsBeta))
!      tmpMQCvar2 = MatMul(CBeta1%subMatrix(newrange2=[1,nElBeta1]),tmpMQCvar1)
!      tmpMQCvar1 = MatMul(Transpose(CBeta1%subMatrix(newrange2=[nElBeta1+1,nBasisUse])),  &
!        MatMul(SMatrixAO,DDNOsBeta))
!      tmpMQCvar3 = MatMul(CBeta1%subMatrix(newrange2=[nElBeta1+1,nBasisUse]),tmpMQCvar1)
!      tmpMatrix2 = float(0)
!      tmpMatrix3 = float(0)
!      do i = 1,nBasis
!        if(float(MQC_Variable_get_MQC(diffDensityBetaEVals,[i])).ge.float(0)) then
!          tmpVector = tmpMQCvar3%column(i)
!          tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
!            float(MQC_Variable_get_MQC(diffDensityBetaEVals,[i])))
!          tmpMatrix2 = tmpMatrix2 + tmpMatrix1
!        else
!          tmpVector = tmpMQCvar2%column(i)
!          tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
!            float(MQC_Variable_get_MQC(diffDensityBetaEVals,[i])))
!          tmpMatrix3 = tmpMatrix3 - tmpMatrix1
!        endIf
!      endDo
!      attachmentDensity = tmpMatrix2
!      detachmentDensity = tmpMatrix3
!      attachmentNumberBeta = contraction(attachmentDensity,SMatrixAO)
!      detachmentNumberBeta = contraction(detachmentDensity,SMatrixAO)
!      write(iOut,2000) ' (Modified Model A)',float(attachmentNumberAlpha),  &
!        float(detachmentNumberAlpha),float(attachmentNumberBeta),  &
!        float(detachmentNumberBeta)
!!
!!     Try promotion number a third time...alpha.
!!
!
!      tmpMQCvar1 = MatMul(Transpose(CAlpha1%subMatrix(newrange2=[1,nElAlpha1])),  &
!        MatMul(SMatrixAO,DDNOsAlpha))
!      tmpMQCvar3 = MatMul(CAlpha1%subMatrix(newrange2=[1,nElAlpha1]),tmpMQCvar1)
!      tmpMQCvar1 = MatMul(Transpose(CAlpha1%subMatrix(newrange2=[nElAlpha1+1,nBasisUse])),  &
!        MatMul(SMatrixAO,DDNOsAlpha))
!      tmpMQCvar4 = MatMul(CAlpha1%subMatrix(newrange2=[nElAlpha1+1,nBasisUse]),tmpMQCvar1)
!      tmpMatrix2 = float(0)
!      tmpMatrix3 = float(0)
!      do i = 1,nBasis
!        tmpVector = tmpMQCvar4%column(i)
!        tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
!          float(MQC_Variable_get_MQC(diffDensityAlphaEVals,[i])))
!        tmpMatrix2 = tmpMatrix2 + tmpMatrix1
!!
!        tmpVector = tmpMQCvar3%column(i)
!        tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
!          float(MQC_Variable_get_MQC(diffDensityAlphaEVals,[i])))
!        tmpMatrix3 = tmpMatrix3 - tmpMatrix1
!      endDo
!      attachmentDensity = tmpMatrix2
!      detachmentDensity = tmpMatrix3
!      attachmentNumberAlpha = contraction(attachmentDensity,SMatrixAO)
!      detachmentNumberAlpha = contraction(detachmentDensity,SMatrixAO)
!!
!!     Try promotion number a third time...beta.
!!
!      tmpMQCvar1 = MatMul(Transpose(CBeta1%subMatrix(newrange2=[1,nElBeta1])),  &
!        MatMul(SMatrixAO,DDNOsBeta))
!      tmpMQCvar3 = MatMul(CBeta1%subMatrix(newrange2=[1,nElBeta1]),tmpMQCvar1)
!      tmpMQCvar1 = MatMul(Transpose(CBeta1%subMatrix(newrange2=[nElBeta1+1,nBasisUse])),  &
!        MatMul(SMatrixAO,DDNOsBeta))
!      tmpMQCvar4 = MatMul(CBeta1%subMatrix(newrange2=[nElBeta1+1,nBasisUse]),tmpMQCvar1)
!      tmpMatrix2 = float(0)
!      tmpMatrix3 = float(0)
!      do i = 1,nBasis
!        tmpVector = tmpMQCvar4%column(i)
!        tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
!          float(MQC_Variable_get_MQC(diffDensityBetaEVals,[i])))
!        tmpMatrix2 = tmpMatrix2 + tmpMatrix1
!!
!        tmpVector = tmpMQCvar3%column(i)
!        tmpMatrix1 = mqc_outerProduct_real(tmpVector,tmpVector,  &
!          float(MQC_Variable_get_MQC(diffDensityBetaEVals,[i])))
!        tmpMatrix3 = tmpMatrix3 - tmpMatrix1
!      endDo
!      attachmentDensity = tmpMatrix2
!      detachmentDensity = tmpMatrix3
!      attachmentNumberBeta = contraction(attachmentDensity,SMatrixAO)
!      detachmentNumberBeta = contraction(detachmentDensity,SMatrixAO)
!      write(iOut,2000) ' (Modified Model B)',float(attachmentNumberAlpha),  &
!        float(detachmentNumberAlpha),float(attachmentNumberBeta),  &
!        float(detachmentNumberBeta)
!
!hph-


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
          call mqc_print_scalar_real(float(dot_product(pDDNO,MQC_Variable_MatrixVector(SMatrixAO,pDDNO))),6,header='(p|p)')
          call mqc_print_scalar_real(float(dot_product(hDDNO,MQC_Variable_MatrixVector(SMatrixAO,hDDNO))),6,header='(h|h)')
          call mqc_print_scalar_real(float(dot_product(hDDNO,MQC_Variable_MatrixVector(SMatrixAO,pDDNO))),6,header='(h|p)')
          call mqc_print_scalar_real(float(dot_product(pDDNO,MQC_Variable_MatrixVector(SMatrixAO,hDDNO))),6,header='(p|h)')
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
        particleXYZ(1) =  dot_product(pDDNO,MQC_Variable_MatrixVector(dipoleAOx,pDDNO))
        particleXYZ(2) =  dot_product(pDDNO,MQC_Variable_MatrixVector(dipoleAOy,pDDNO))
        particleXYZ(3) =  dot_product(pDDNO,MQC_Variable_MatrixVector(dipoleAOz,pDDNO))
        holeXYZ(1) =  dot_product(hDDNO,MQC_Variable_MatrixVector(dipoleAOx,hDDNO))
        holeXYZ(2) =  dot_product(hDDNO,MQC_Variable_MatrixVector(dipoleAOy,hDDNO))
        holeXYZ(3) =  dot_product(hDDNO,MQC_Variable_MatrixVector(dipoleAOz,hDDNO))
        call mqc_print(particleXYZ,6,header='<x> particle',blank_at_top=.true.)
        call mqc_print(holeXYZ,6,header='<x> hole',blank_at_top=.true.)
      endIf

      if(isNIO.or..not.doTestCode) goto 998

!hph+
      goto 998
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

      goto 998

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
!       HPH -- Do a test to see if NBasis.ne.NBasisUse and then print out
!       dimensions of the arrays being written in the lines below.
!
        if(nBasis.ne.nBasisUse) then
          write(iOut,*)
          write(iOut,*)' nBasis != NBasisUse'
          write(iOut,*)'            Dim(DDNOsAlpha)'
          write(iOut,*)'                nRows = ',Size(DDNOsAlpha,1)
          write(iOut,*)'                nCols = ',Size(DDNOsAlpha,2)
          write(iOut,*)'            Dim(DDNOsBeta)'
          write(iOut,*)'                nRows = ',Size(DDNOsBeta,1)
          write(iOut,*)'                nCols = ',Size(DDNOsBeta,2)
          write(iOut,*)'            Dim(diffDensityAlphaEVals)'
          write(iOut,*)'                nRows = ',Size(diffDensityAlphaEVals)
          write(iOut,*)'            Dim(diffDensityBetaEVals)'
          write(iOut,*)'                nRows = ',Size(diffDensityBetaEVals)
          tmpVector = diffDensityAlphaEVals
          write(iOut,*)'            Max pruned alpha eval = ',MaxVal(ABS(tmpVector(nBasisUse+1:)))
          diffDensityAlphaEVals = tmpVector(1:nBasisUse)
          tmpVector = diffDensityBetaEVals
          write(iOut,*)'            Max pruned Beta eval  = ',MaxVal(ABS(tmpVector(nBasisUse+1:)))
          diffDensityBetaEVals = tmpVector(1:nBasisUse)
          tmpMatrix1 = DDNOsAlpha
          DDNOsAlpha = tmpMatrix1(:,1:nBasisUse)
          tmpMatrix1 = DDNOsBeta
          DDNOsBeta = tmpMatrix1(:,1:nBasisUse)
          write(iOut,*)
          write(iOut,*)' AFTER pruning arrays'
          write(iOut,*)'            Dim(DDNOsAlpha)'
          write(iOut,*)'                nRows = ',Size(DDNOsAlpha,1)
          write(iOut,*)'                nCols = ',Size(DDNOsAlpha,2)
          write(iOut,*)'            Dim(DDNOsBeta)'
          write(iOut,*)'                nRows = ',Size(DDNOsBeta,1)
          write(iOut,*)'                nCols = ',Size(DDNOsBeta,2)
          write(iOut,*)'            Dim(diffDensityAlphaEVals)'
          write(iOut,*)'                nRows = ',Size(diffDensityAlphaEVals)
          write(iOut,*)'            Dim(diffDensityBetaEVals)'
          write(iOut,*)'                nRows = ',Size(diffDensityBetaEVals)
          write(iOut,*)
        endIf
!
!       DDNO eigenvectors and eigenvalues...
        call GMatrixFileOut%writeArray2('ALPHA ORBITAL ENERGIES',diffDensityAlphaEVals)
        call GMatrixFileOut%writeArray2('BETA ORBITAL ENERGIES',diffDensityBetaEVals)
        call GMatrixFileOut%writeArray2('ALPHA MO COEFFICIENTS',DDNOsAlpha)
        call GMatrixFileOut%writeArray2('BETA MO COEFFICIENTS',DDNOsBeta)
        call GMatrixFileOut%writeArray2('ALPHA SCF DENSITY MATRIX',PMatrixAlpha1,storage='l')
        call GMatrixFileOut%writeArray2('BETA SCF DENSITY MATRIX',PMatrixBeta1,storage='l')
!hph        call GMatrixFileOut%writeArray2('BETA SCF DENSITY MATRIX',PMatrixBeta1)
!
!       Close out the matrix file.
        call GMatrixFileOut%closeFile()
      endIf
!
  999 Continue
      write(iOut,8999)
      end program nio
