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
      integer(kind=int64)::nCommands,i,j,k1,k2,nAtoms,nAtoms2,nAt3,  &
        nBasis,nBasis2,nBasisUse,nBasisUse2
      integer(kind=int64),dimension(:),allocatable::atomicNumbers
      real(kind=real64)::Vnn,Escf
      real(kind=real64),dimension(3)::tmp3Vec
      real(kind=real64),dimension(:),allocatable::cartesians,cartesians2
      real(kind=real64),dimension(:,:),allocatable::distanceMatrix
      character(len=512)::matrixFilename1,matrixFilename2,tmpString
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFile1,GMatrixFile2
      type(MQC_Variable)::tmpMQCvar
      type(MQC_Variable)::nEalpha,nEbeta,nEtot,KEnergy,VEnergy,OneElEnergy,  &
        TwoElEnergy,scfEnergy
      type(MQC_Variable)::SMatrixAO,TMatrixAO,VMatrixAO,HCoreMatrixAO,  &
        FMatrixAlpha,FMatrixBeta,PMatrixAlpha1,PMatrixBeta1,PMatrixTotal1,  &
        PMatrixAlpha2,PMatrixBeta2,PMatrixTotal2,  &
        ERIs,JMatrixAlpha,KMatrixAlpha
      type(MQC_R4Tensor)::tmpR4
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
!     Figure out nAt3, then allocate memory for key arrays.
!
      nAt3 = 3*nAtoms
      Allocate(cartesians(NAt3),cartesians2(NAt3),atomicNumbers(NAtoms))
!
!     Load the atomic numbers and Cartesian coordinates.
!
      atomicNumbers = GMatrixFile1%getAtomicNumbers()
      cartesians = GMatrixFile1%getAtomCarts()
      cartesians = cartesians*angPBohr
      cartesians2 = GMatrixFile2%getAtomCarts()
      cartesians2 = cartesians2*angPBohr
!
!     Load the atomic orbital overlap matrix and the density matrices.
!
      call GMatrixFile1%getArray('OVERLAP',mqcVarOut=SMatrixAO)
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
      







!hph+      
!!
!!     Load up a few matrices from the matrix file.
!!
!      call GMatrixFile%getArray('KINETIC ENERGY',mqcVarOut=TMatrixAO)
!      call GMatrixFile%getArray('CORE HAMILTONIAN ALPHA',mqcVarOut=HCoreMatrixAO)
!      call GMatrixFile%getArray('ALPHA FOCK MATRIX',mqcVarOut=FMatrixAlpha)
!      if(GMatrixFile%isUnrestricted()) then
!        call GMatrixFile%getArray('BETA FOCK MATRIX',mqcVarOut=FMatrixBeta)
!      else
!        FMatrixBeta  = FMatrixAlpha
!      endIf
!      VMatrixAO = HCoreMatrixAO-TMatrixAO
!
!!hph+
!      tmpMQCvar = FMatrixAlpha-HCoreMatrixAO
!      call tmpMQCvar%print(header='F-H')
!      MQC_Gaussian_DEBUGHPH = .True.
!      call GMatrixFile%getArray('REGULAR 2E INTEGRALS',mqcVarOut=ERIs)
!      call ERIs%print(IOut,' ERIs=')
!      write(*,*)' 1,2,2,1 = ',float(ERIs%getVal([1,2,2,1]))
!      write(*,*)' 1,2,2,2 = ',float(ERIs%getVal([1,2,2,2]))
!      tmpMQCvar = HCoreMatrixAO
!      write(iOut,*)
!      write(iOut,*)' Forming Coulomb matrix...'
!      call formCoulomb(Int(GMatrixFile%getVal('nbasis')),PMatrixAlpha,ERIs,tmpMQCvar,initialize=.true.)
!      call formCoulomb(Int(GMatrixFile%getVal('nbasis')),PMatrixBeta,ERIs,tmpMQCvar,initialize=.false.)
!      JMatrixAlpha = tmpMQCvar
!      write(iOut,*)
!      write(iOut,*)' Forming Exchange matrix...'
!      call formExchange(Int(GMatrixFile%getVal('nbasis')),PMatrixAlpha,ERIs,tmpMQCvar,initialize=.false.)
!      KMatrixAlpha = tmpMQCvar
!      call PMatrixAlpha%print(header='PAlpha')
!      call JMatrixAlpha%print(header='JAlpha')
!      call KMatrixAlpha%print(header='KAlpha')
!      call tmpMQCvar%print(header='temp fock matrix without H',blankAtTop=.true.)
!      tmpMQCvar = tmpMQCvar + HCoreMatrixAO
!      call tmpMQCvar%print(header='temp fock matrix with H',blankAtTop=.true.)
!      call FMatrixAlpha%print(header='fock matrix from Gaussian')
!
!!      call mqc_error('Hrant - STOP')
!!hph-
!
!
!!
!!     Calculate the number of electrons using <PS>.
!!
!      nEalpha = Contraction(PMatrixAlpha,SMatrixAO)
!      nEbeta  = Contraction(PMatrixBeta,SMatrixAO)
!      nEtot   = Contraction(PMatrixTotal,SMatrixAO)
!      call nEalpha%print(IOut,' <P(Alpha)S>=')
!      call nEbeta%print(IOut,' <P(Beta )S>=')
!      call nEtot%print(IOut,' <P(Total)S>=')
!!
!!     Calculate the 1-electron energy and component pieces of the 1-electron
!!     energy. Also, calculate the 2-electron energy.
!!
!      KEnergy     = Contraction(PMatrixTotal,TMatrixAO)
!      VEnergy     = Contraction(PMatrixTotal,VMatrixAO)
!      OneElEnergy = Contraction(PMatrixTotal,HCoreMatrixAO)
!      TwoElEnergy = Contraction(PMatrixTotal,FMatrixAlpha)
!      if(GMatrixFile%isUnrestricted()) then
!        tmpMQCvar = Contraction(PMatrixTotal,FMatrixBeta)
!        TwoElEnergy = TwoElEnergy + tmpMQCvar
!      endIf
!      TwoElEnergy = TwoElEnergy - OneElEnergy
!      TwoElEnergy = MQC(0.5)*TwoElEnergy
!      call KEnergy%print(IOut,' <P.K> = ')
!      call VEnergy%print(IOut,' <P.V> =')
!      call OneElEnergy%print(IOut,' <P.H> = ')
!      call TwoElEnergy%print(IOut,' <P.F>-<P.H> = ')
!!
!!     Load the atommic numbers and Cartesian coordinates into our intrinsic
!!     arrays.
!!
!      atomicNumbers = GMatrixFile%getAtomicNumbers()
!      cartesians = GMatrixFile%getAtomCarts()
!      cartesians = cartesians*angPBohr
!!
!!     Print out the atomic numbers and Cartesian coordiantes for each atomic
!!     center.
!!
!      write(IOut,1200)
!      do i = 1,NAtoms
!        j = 3*(i-1)
!        write(IOut,1210) i,mqc_element_symbol(atomicNumbers(i)),  &
!          cartesians(j+1),cartesians(j+2),cartesians(j+3)
!      endDo
!!
!!     Form the distance matrix between atomic centers.
!!
!      Allocate(distanceMatrix(nAtoms,nAtoms))
!      do i = 1,nAtoms-1
!        distanceMatrix(i,i) = float(0)
!        k1 = 3*(i-1)+1
!        do j = i+1,NAtoms
!          k2 = 3*(j-1)+1
!          tmp3Vec = cartesians(k1:k1+2)-cartesians(k2:k2+2)
!          distanceMatrix(i,j) = sqrt(dot_product(tmp3Vec,tmp3Vec))
!          distanceMatrix(j,i) = distanceMatrix(i,j)
!        endDo
!      endDo
!!
!!     Calculate the nuclear-nuclear repulsion energy.
!!
!      distanceMatrix = distanceMatrix/angPBohr
!      Vnn = float(0)
!      do i = 1,NAtoms-1
!        do j = i+1,NAtoms
!          Vnn = Vnn + float(atomicNumbers(i)*atomicNumbers(j))/distanceMatrix(i,j)
!        endDo
!      endDo
!      write(iOut,1300) Vnn
!!
!!     Put things together and report the SCF energy.
!!
!      scfEnergy = oneElEnergy + twoElEnergy
!      scfEnergy = scfEnergy + MQC(Vnn)
!      call scfEnergy%print(IOut,' SCF Energy = ')
!hph-


!
  999 Continue
      write(iOut,8999)
      end program nio
