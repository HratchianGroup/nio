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
        nBasis,nBasis2,nBasisUse,nBasisUse2
      character(len=512)::matrixFilename1,matrixFilename2
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile1,GMatrixFile2
      type(MQC_Variable)::SMatrixAO,SMatrixEVecs,SMatrixEVals,  &
        SMatrixAOHalf,SMatrixAOMinusHalf
      type(MQC_Variable)::tmpMQCvar
      type(MQC_Variable)::PMatrixAlpha1,PMatrixBeta1,PMatrixTotal1,  &
        PMatrixAlpha2,PMatrixBeta2,PMatrixTotal2,diffDensityAlpha,  &
        diffDensityBeta,diffDensityAlphaEVecs,diffDensityAlphaEVals,  &
        diffDensityBetaEVecs,diffDensityBetaEVals
!
!     Format Statements
!
 1000 Format(1x,'Enter Program NIO.')
 1010 Format(3x,'Matrix File 1: ',A,/,  &
             3x,'Matrix File 2: ',A,/)
 1100 Format(1x,'nAtoms=',I4,3x,'nBasis=',I4,3x,'nBasisUse=',I4)
 1200 Format(1x,'Atomic Coordinates (Angstrom)')
 1210 Format(3x,I3,2x,A2,5x,F7.4,3x,F7.4,3x,F7.4)
 1300 Format(1x,'Nuclear Repulsion Energy = ',F20.6)
 8999 Format(/,1x,'END OF TEST PROGRAM scfEnergyTerms.')
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
      write(IOut,1100) nAtoms,nBasis,nBasisUse
      write(iOut,*)
      write(iOut,*)' Hrant - the unit number for the first  matrix file is: ',GMatrixFile1%UnitNumber
      write(iOut,*)' Hrant - the unit number for the second matrix file is: ',GMatrixFile2%UnitNumber
      write(iOut,*)
!
!     Load the atomic orbital overlap matrix and form S^(1/2) and S^(-1/2).
!
      call GMatrixFile1%getArray('OVERLAP',mqcVarOut=SMatrixAO)
      call SMatrixAO%print(header='Overlap Matrix')
      call SMatrixAO%eigen(SMatrixEVals,SMatrixEVecs)
      call SMatrixEVals%print(header='S matrix eigen-values:')

      call mqc_print(MatMul(Transpose(SMatrixEVecs),SMatrixEVecs),header='SEVecs(t).SEVecs')
      call mqc_print(MatMul(MatMul(SMatrixEVecs,SMatrixEVals%diag()),TRANSPOSE(SMatrixEVecs)),6,'U.lambda.Ut')

      tmpMQCvar = SMatrixEVals%rpower(0.5)
      SMAtrixAOHalf = MatMul(MatMul(SMatrixEVecs,tmpMQCvar%diag()),TRANSPOSE(SMatrixEVecs))
      tmpMQCvar = SMatrixEVals%rpower(-0.5)
      SMAtrixAOMinusHalf = MatMul(MatMul(SMatrixEVecs,tmpMQCvar%diag()),TRANSPOSE(SMatrixEVecs))

      call SMAtrixAOHalf%print(header='S**(1/2)')
      call SMatrixAOMinusHalf%print(header='S**(-1/2)')

!
!     Load the density matrices. The code below treats all systems as open
!     shell, so closed shell results are handled by copying the density matrix
!     from restricted calculations into alpha and beta density matrix arrays.
!
      call GMatrixFile1%getArray('ALPHA DENSITY MATRIX',mqcVarOut=PMatrixAlpha1)
      if(GMatrixFile1%isUnrestricted()) then
        call GMatrixFile1%getArray('BETA DENSITY MATRIX',mqcVarOut=PMatrixBeta1)
      else
        PMatrixBeta1  = PMatrixAlpha1
      endIf
      PMatrixTotal1 = PMatrixAlpha1+PMatrixBeta1
      call GMatrixFile2%getArray('ALPHA DENSITY MATRIX',mqcVarOut=PMatrixAlpha2)
      if(GMatrixFile2%isUnrestricted()) then
        call GMatrixFile2%getArray('BETA DENSITY MATRIX',mqcVarOut=PMatrixBeta2)
      else
        PMatrixBeta2  = PMatrixAlpha2
      endIf
      PMatrixTotal2 = PMatrixAlpha2+PMatrixBeta2
!
!     Form the difference density and construct the NIOs.
!
      diffDensityAlpha = PMatrixAlpha2-PMatrixAlpha1
      diffDensityBeta  = PMatrixBeta2-PMatrixBeta1

      call mqc_print(contraction(diffDensityAlpha,SMatrixAO),header='P(alpha).S')
      call mqc_print(contraction(diffDensityBeta,SMatrixAO),header='P(beta).S')

      tmpMQCvar = MatMul(SMatrixAOHalf,MatMul(diffDensityAlpha,SMatrixAOHalf))
      call tmpMQCvar%eigen(diffDensityAlphaEVals,diffDensityAlphaEVecs)
      tmpMQCvar = MatMul(SMatrixAOHalf,MatMul(diffDensityBeta,SMatrixAOHalf))
      call tmpMQCvar%eigen(diffDensityBetaEVals,diffDensityBetaEVecs)

      call diffDensityAlphaEVals%print(header='Alpha Occ Change EVals')
      call diffDensityBetaEVals%print(header='Beta Occ Change EVals')

      goto 999


!
  999 Continue
      write(iOut,8999)
      end program nio
