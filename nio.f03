INCLUDE 'nio_mod.f03'
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
      integer(kind=int64)::nCommands,i,nAtoms,nAtoms2,  &
        nBasis,nBasis2,nBasisUse,nBasisUse2,nEl1,nEl2,nElAlpha1,  &
        nElBeta1,nElAlpha2,NElBeta2,nPlusOneAlpha,nMinusOneAlpha,  &
        iPlusOneAlpha,iMinusOneAlpha,nPlusOneBeta,nMinusOneBeta,  &
        iPlusOneBeta,iMinusOneBeta,nPlusOne,nMinusOne,  &
        nRelaxationDDNOsAlpha,nRelaxationDDNOsBeta
      real(kind=real64),dimension(3)::transitionDipole
      character(len=512)::matrixFilename1,matrixFilename2
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile1,  &
        GMatrixFile2,GMatrixFileOut
      type(MQC_Variable)::DDNOsAlpha,DDNOsBeta,pDDNO,hDDNO,  &
        dipoleStrength,TDOverlapA,TDOverlapB,TDparticleHoleMag
      type(MQC_Variable)::infoDDNOsAlpha,infoDDNOsBeta
      type(MQC_Variable)::SMatrixAO,SMatrixEVecs,SMatrixEVals,  &
        SMatrixAOHalf,SMatrixAOMinusHalf
      type(MQC_Variable)::PMatrixAlpha1,PMatrixBeta1,PMatrixTotal1,  &
        PMatrixAlpha2,PMatrixBeta2,PMatrixTotal2,diffDensityAlpha,  &
        diffDensityBeta,diffDensityAlphaEVecs,diffDensityAlphaEVals,  &
        diffDensityBetaEVecs,diffDensityBetaEVals
      type(MQC_Variable)::CAlpha1,CBeta1,CAlpha2,CBeta2,TAlpha,TBeta
      type(MQC_Variable)::dipoleAOx,dipoleAOy,dipoleAOz
      type(MQC_Variable)::tmpMQCvar,tmpMQCvar1,tmpMQCvar2,tmpMQCvar3,  &
        tmpMQCvar4
      logical::doTestCode=.False.,isNIO,isDDNO
!
!     Format Statements
!
 1000 Format(1x,'Enter Program NIO.')
 1010 Format(3x,'Matrix File 1: ',A,/,  &
             3x,'Matrix File 2: ',A,/)
 1100 Format(1x,'nAtoms=',I4,3x,'nBasis=',I4,3x,'nBasisUse=',I4,/,  &
        1x,'nEl1=',I4,3x,'nElAlpha1=',I4,3x,'nElBeta=',I4,/,  &
        1x,'nEl2=',I4,3x,'nElAlpha2=',I4,3x,'nElBeta=',I4,/)
 1200 Format(1x,'Atomic Coordinates (Angstrom)')
 1210 Format(3x,I3,2x,A2,5x,F7.4,3x,F7.4,3x,F7.4)
 1300 Format(1x,'Nuclear Repulsion Energy = ',F20.6)
 1500 Format(/,1x,'Overlap between Delta-SCF states = ',F9.6,/,  &
        3x,'Alpha Overlap = ',F9.6,3x,'Beta Overlap = ',F9.6,/)
 2000 Format(/,1x,'nPlusOneAlpha=',I2,3x,'nMinusAlpha=',I2,/,  &
        1x,'nPlusOneBeta =',I2,3x,'nMinusBeta =',I2,/,  &
        1x,'isNIO=',L1,3x,'isDDNO=',L1)
 8999 Format(/,1x,'END OF NIO PROGRAM')
!
!
      write(IOut,1000)
!
!     Open the Gaussian matrix file and load the number of atomic centers.

      nCommands = command_argument_count()
      if(nCommands.ne.2)  &
        call mqc_error('Two input Gaussian matrix files must be provided in the command line.')
      call get_command_argument(1,matrixFilename1)
      call get_command_argument(2,matrixFilename2)
      call GMatrixFile1%load(matrixFilename1)
      call GMatrixFile2%load(matrixFilename2)
      write(IOut,1010) TRIM(matrixFilename1),TRIM(matrixFilename2)
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


!hph+
!      tmpMQCvar = MQC_Gaussian_Unformatted_Matrix_Get_Value_Real(GMatrixFile1,'gaussian scalars')
!hph-



!
!     Load the atomic orbital overlap matrix and form S^(1/2) and S^(-1/2).
!
      call GMatrixFile1%getArray('OVERLAP',mqcVarOut=SMatrixAO)
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
      call mqc_print(contraction(diffDensityAlpha,SMatrixAO),header='P(alpha).S')
      call mqc_print(contraction(diffDensityBeta,SMatrixAO),header='P(beta).S')
      tmpMQCvar = MatMul(SMatrixAOHalf,MatMul(diffDensityAlpha,SMatrixAOHalf))
      call tmpMQCvar%eigen(diffDensityAlphaEVals,diffDensityAlphaEVecs)
      DDNOsAlpha = MatMul(SMatrixAOMinusHalf,diffDensityAlphaEVecs)
      tmpMQCvar = MatMul(SMatrixAOHalf,MatMul(diffDensityBeta,SMatrixAOHalf))
      call tmpMQCvar%eigen(diffDensityBetaEVals,diffDensityBetaEVecs)
      DDNOsBeta  = MatMul(SMatrixAOMinusHalf,diffDensityBetaEVecs)
      call diffDensityAlphaEVals%print(header='Alpha Occ Change EVals')
      call diffDensityBetaEVals%print(header='Beta Occ Change EVals')
!
!     Form the polestrength (for detachment cases) or the N-1 overlap (for
!     excitation cases). At the end of this block, we decide if this is a
!     detachment (<isNIO>) or excitation job (<isDDNO>).
!
      call determinantOverlap(SMatrixAO,SMatrixAOMinusHalf,  &
        diffDensityAlphaEVals,diffDensityAlphaEVecs,CAlpha2,nElAlpha2,  &
        nBasis,TDOverlapA,nPlusOneAlpha,nMinusOneAlpha,iPlusOneAlpha,  &
        iMinusOneAlpha)
      call TDOverlapA%print(header='determinant overlap for ALPHA')
      call determinantOverlap(SMatrixAO,SMatrixAOMinusHalf,  &
        diffDensityBetaEVals,diffDensityBetaEVecs,CBeta2,nElBeta2,  &
        nBasis,TDOverlapB,nPlusOneBeta,nMinusOneBeta,iPlusOneBeta,  &
        iMinusOneBeta)
      call TDOverlapB%print(header='determinant overlap for BETA ')
      isNIO  = ((nPlusOneAlpha+nPlusOneBeta).eq.0.and.  &
        (nMinusOneAlpha+nMinusOneBeta).eq.1)
      isDDNO = ((nPlusOneAlpha+nPlusOneBeta).eq.1.and.  &
        (nMinusOneAlpha+nMinusOneBeta).eq.1)
      write(iOut,2000) nPlusOneAlpha,nMinusOneAlpha,nPlusOneBeta,  &
        nMinusOneBeta,isNIO,isDDNO
!
!     If this is a DDNO job, calculate the overlap of the two determinants.
!
      if(isDDNO) then
        tmpMQCvar1 = MatMul(Transpose(CAlpha1%subMatrix(newrange2=[1,nElAlpha1])),  &
          MatMul(SMatrixAO,CAlpha2%subMatrix(newrange2=[1,nElAlpha2])))
        tmpMQCvar2 = MatMul(Transpose(CBeta1),MatMul(SMatrixAO,CBeta2))
        tmpMQCvar3 = tmpMQCvar1%det()
        tmpMQCvar4 = tmpMQCvar2%det()
        tmpMQCvar = tmpMQCvar3*tmpMQCvar4
        write(iOut,1500) float(tmpMQCvar),float(tmpMQCvar3),float(tmpMQCvar4)
      endIf
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
        write(*,*)' iPlusOneAlpha  = ',iPlusOneAlpha
        write(*,*)' iMinusOneAlpha = ',iMinusOneAlpha
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
        call pDDNO%print(header='particle DDNO')
        call hDDNO%print(header='hole DDNO')
        if(DEBUG) then
          call mqc_print_scalar_real(6,float(dot_product(pDDNO,MQC_Variable_MatrixVector(SMatrixAO,pDDNO))),header='<p|p>')
          call mqc_print_scalar_real(6,float(dot_product(hDDNO,MQC_Variable_MatrixVector(SMatrixAO,hDDNO))),header='<h|h>')
          call mqc_print_scalar_real(6,float(dot_product(hDDNO,MQC_Variable_MatrixVector(SMatrixAO,pDDNO))),header='<h|p>')
          call mqc_print_scalar_real(6,float(dot_product(pDDNO,MQC_Variable_MatrixVector(SMatrixAO,hDDNO))),header='<p|h>')
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
        call mqc_print(6,transitionDipole,header='Transition Dipole Moment')
        TDparticleHoleMag = dot_product(transitionDipole,transitionDipole)
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
      endIf


      if(.not.doTestCode) goto 999


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






!
!     Write results to the output Gaussian matrix file.
!
      GMatrixFileOut = GMatrixFile1
      call GMatrixFileOut%create('new.mat')
      write(*,*)
      write(*,*)' Hrant - filename = ',TRIM(GMatrixFileOut%filename)
!
!     Basis set info...
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
!     DDNO eigenvectors and eigenvalues...
      call GMatrixFileOut%writeArray2('ALPHA ORBITAL ENERGIES',diffDensityAlphaEVals)
      call GMatrixFileOut%writeArray2('BETA ORBITAL ENERGIES',diffDensityBetaEVals)
      call GMatrixFileOut%writeArray2('ALPHA MO COEFFICIENTS',DDNOsAlpha)
      call GMatrixFileOut%writeArray2('BETA MO COEFFICIENTS',DDNOsBeta)
!
!     Close out the matrix file.
      call GMatrixFileOut%closeFile()
!
  999 Continue
      write(iOut,8999)
      end program nio
