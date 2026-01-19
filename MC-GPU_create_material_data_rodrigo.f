C  ****  Files included to simplify compilation.
      INCLUDE 'penelope.f'
      INCLUDE 'rita.f'
	  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   This program was adapted from Andreu Badal program to generate     C
C   the material tables to MC-GPU MC code from PENELOPE.               C
C   We updated the functions to work with PENELOPE 2018 version.       C 
C                              Rodrigo Massera, 2019-09-24             C
C                                                                      C
C   This program reads a PENELOPE 2018 material file and outputs a     C
C   table with photon interaction mean free paths (MFP), and data for  C
C   Rayleigh and Compton interaction sampling.                         C
C   In PENELOPE 2006 for Rayleigh scattering the X unity was used.     C
C   However, in newer versions, the Q/(mec2) is the new standard.      C
C   The relation bewteen them is given by eq 2.14 in PENELOPE 2014     C
C   manual:                                                            C
C    X =  20.6074 x (q/mec)                                            C
C	Where me is the electron rest mass and c the speed of light.       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   ORIGINAL VERSION:                                                  C
C   This program reads a PENELOPE 2006 material file and outputs a     C
C   table with photon interaction mean free paths (MFP), and data for  C
C   Rayleigh and Compton interaction sampling.                         C
C                                                                      C
C   While the PENELOPE database is linearly interpolated in LOG-LOG,   C
C   the energy grid in the output table is a linear, ie, has equally   C
C   spaced energy bins. A small bin width is required to allow direct  C
C   linear interpolation of the MFP, avoiding the LOG computation.     C
C                                                                      C
C   This source code is based on PENELOPE's "tables.f".                C
C                                                                      C 
C                              Andreu Badal, 2009-03-31                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C --Copyright notice from tables.f and penelope.f:                     C
C                                                                      C
C  PENELOPE/PENGEOM (version 2018)                                     C
C  Copyright (c) 2001-2018                                             C
C  Universitat de Barcelona                                            C
C                                                                      C
C  Permission to use, copy, modify, distribute and sell this software  C
C  and its documentation for any purpose is hereby granted without     C
C  fee, provided that the above copyright notice appears in all        C
C  copies and that both that copyright notice and this permission      C
C  notice appear in all supporting documentation. The Universitat de   C
C  Barcelona makes no representations about the suitability of this    C
C  software for any purpose. It is provided 'as is' without express    C
C  or implied warranty.                                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C  *********************************************************************
C                       MAIN PROGRAM
C  *********************************************************************

      USE PENELOPE_mod
      USE PENERROR_mod
	  
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*20 MFNAME, OUTNAME 
      CHARACTER*62 NAME
	  CHARACTER PMFILE*20
	  DIMENSION PMFILE(MAXMAT)
      
C  ****  Auxiliary arrays.
      DIMENSION E_MFP(6)
      PARAMETER (MAX_ENERGY_BINS=60005)
      DIMENSION PMAX_linear_energy(MAX_ENERGY_BINS)
            
C  ****  Simulation parameters.
C  ****  PARAMETER (MAXMAT=10)
C  ****  COMMON/CSIMPA/EABS(3,MAXMAT),C1(MAXMAT),C2(MAXMAT),WCC(MAXMAT),
C  ****  1  WCR(MAXMAT)
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)

C  ****  Penelope energy grid and Rayleigh sampling data:
C  ****  PARAMETER (NEGP=200)
C  ****      PARAMETER (NP=128,NPM1=NP-1)
C  ****      COMMON/CGRA/XCO(NP,MAXMAT),PCO(NP,MAXMAT),ACO(NP,MAXMAT),
C  ****     1  BCO(NP,MAXMAT),PMAX(NEGP,MAXMAT),ITLCO(NP,MAXMAT),
C  ****     2  ITUCO(NP,MAXMAT)

C  ****  Photon simulation tables.
      COMMON/CGIMFP/SGRA(MAXMAT,NEGP),SGCO(MAXMAT,NEGP),
     1  SGPH(MAXMAT,NEGP),SGPP(MAXMAT,NEGP),SGAUX(MAXMAT,NEGP)

C  ****  Rayleigh scattering.
      PARAMETER (NQ=250,NEX=1024)
C      COMMON/CGRA00/FACTE,Q2MAX,MM,MOM
      COMMON/CGRA01/FF(MAXMAT,NQ),ERA(NEX),XSRA(MAXMAT,NEX),
     1    IED(NEGP),IEU(NEGP),NE
      COMMON/CGRA02/QQ(NQ),AR(MAXMAT,NQ),BR(MAXMAT,NQ),CR(MAXMAT,NQ),
     1              DR(MAXMAT,NQ),FF0(MAXMAT),QQM
      PARAMETER (NP=150,NPM1=NP-1)
      COMMON/CGRA03/QRA(NP,MAXMAT),PRA(NP,MAXMAT),DPRA(NP,MAXMAT),
     1  ARA(NP,MAXMAT),BRA(NP,MAXMAT),PMAX(NEGP,MAXMAT),
     2  ITLRA(NP,MAXMAT),ITURA(NP,MAXMAT)

C  ****  Energy grid and interpolation constants for the current energy.
      COMMON/CEGRID/EMIN,EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE

C  ****  Compton scattering.
      PARAMETER (NOCO=512)
      COMMON/CGCO/FCO(MAXMAT,NOCO),UICO(MAXMAT,NOCO),FJ0(MAXMAT,NOCO),
     1  PTRSH(MAXMAT,NOCO),KZCO(MAXMAT,NOCO),KSCO(MAXMAT,NOCO),
     2  NOSCCO(MAXMAT)
	 
	 !(Q/mec)2 to X2 
	 PARAMETER (Q2X=424.66493475999994D0)



      WRITE(6,*)" "
      WRITE(6,*)" "
      WRITE(6,*)"         ***********************************"
      WRITE(6,*)"         *** MC-GPU_create_material_data ***"
      WRITE(6,*)"         ***********************************"
      WRITE(6,*)" "
      WRITE(6,*)" "
      WRITE(6,*)"    Creating a material input file for MC-GPU."
      WRITE(6,*)
     &"    This program reads a PENELOPE 2018 material file and outputs"
      WRITE(6,*)
     & "    a table with photon interaction mean free paths (MFP) and"
      WRITE(6,*)
     & "    data for Rayleigh and Compton interaction sampling."
C
C  ****  Parameters (to tabulate the complete energy range and to switch
C        soft interactions off).
C

C
C  ****  Material data file.
C
      WRITE(6,'(a)') '  '
      WRITE(6,'(a)') '  -- Enter the energy range to tabulate: '//
     &           ' Emin, Emax (eg, 5000  125000):'
      READ(5,*) EMIN, EMAX
      WRITE(6,'(a)') '  -- Enter the number of energy bins (eg, 8192):'    ! 8192 = 2^13
      READ(5,*) NBINS
      DE=(EMAX-EMIN)/DBLE(NBINS-1) !added -1 to include the max value
      WRITE(6,'(a,1pe17.10)')
     &      '      - Energy bin width set to (EMAX-EMIN)/NBINS = ',DE
      WRITE(6,'(a)') '  -- Enter the name of the PENELOPE 2006'//
     &               ' material data file (eg, water.mat):'
      READ(5,'(A20)') MFNAME
      WRITE(6,'(a)') '  -- Enter the name of the output data file'//
     &               'for MC-GPU (eg, water.mcgpu)...'
      READ(5,'(A20)') OUTNAME
      WRITE(6,'('' Material data file:  '', A40)') MFNAME
      WRITE(6,'(a)') '  '
      WRITE(6,'(a)') 'Processing material data. Please, wait...'


      ! -- Initializing PENELOPE with the material information:
      !    Tabulate the material tables between the input maximum and minimum energies.
      DO M=1,MAXMAT
        EABS(1,M) = EMIN
        EABS(2,M) = EMIN
        EABS(3,M) = EMIN
        C1(M)     =  0.0D0
        C2(M)     =  0.0D0
        WCC(M)    =  0.0D0
        WCR(M)    =-10.0D0
      ENDDO
	  
	  !starts PENELOPE to load the coefficients
	  PMFILE(1)=MFNAME
      OPEN(16,FILE='material.dat')
	  !INFO can be changed from 0 (none) to 5 (max verbose)
	  INFO =5 
	  CALL PEINIT(EMAX,1,16,INFO,PMFILE)
	  CLOSE(16)
	  
	  write(*,*) 'done'

      ! -- Re-open the material file and read the material name (2nd line):
      OPEN(11,FILE=MFNAME)
      READ(11,'(A55)') NAME
      READ(11,'(11X,a62)') NAME
      CLOSE(11)
 

C  ****  Calculate photon mean free paths:
C    ** Function PHMFP returns the mean free path,MFP, [cm] for the input energy, kind of particle,
C    **  material number (from input file), and kind of interaction.                         
C    ** The cross section is found dividing the inverse MFP by the molar volume [atoms/cm^3].
C
      WRITE(6,*)'  '
c      WRITE(6,*)'====================================================='
c      write(6,*)' PENELOPEs function PHMFP returns the mean free'//
c     &          ' path [cm] for the input energy, kind of particle, '//
c     &          ' material number, and kind of interaction: '
c      write(6,*)'    MFP=PHMFP(E,KPAR,M,ICOL)'
c      write(6,*)' The cross section is found dividing the inverse'//
c     &          ' mean free path (=attenuation coefficient) by the'//
c     &          ' molar volume [atoms/cm^3].'
c      write(6,*)'    XS=(1.0D0/MFP)/VMOL(M)'
c     WRITE(6,*)'====================================================='
c      WRITE(6,*) '  '      

      ! Set mat number and particle:
      M=1                ! Use first material defined in the input material file
      KPAR = 2           ! Select photons (1=electron, 2=photon, 3=positron)
      
      ! -- Open output file:
      OPEN(1, FILE=OUTNAME)

      ! -- Write file header:
      WRITE(1,'(a)')'#[MATERIAL DEFINITION FOR MC-GPU: interaction'//
     &      ' mean free path and sampling data from PENELOPE 2018]'
      WRITE(1,'(a)')'#[MATERIAL NAME]'
      WRITE(1,1001) NAME
 1001 format('# ',a)
      WRITE(1,'(a)')'#[NOMINAL DENSITY (g/cm^3)]'
      WRITE(1,1002) RHO(M)
 1002 format('# ',f12.8)
      WRITE(1,'(a)')'#[NUMBER OF DATA VALUES]'
      WRITE(1,1003) NBINS
 1003 format('# ',I6)
      WRITE(1,'(a)')'#[MEAN FREE PATHS (cm)'//
     &        ' (ie, average distance between interactions)]'
      WRITE(1,'(a)') '#[Energy (eV)     | Rayleigh        |'//
     &                ' Compton         | Photoelectric   |'//           
     &                ' TOTAL (+pair prod) (cm) |'//            !  &  ' Pair-production | TOTAL (cm) |'//
     &                ' Rayleigh: max cumul prob F^2]'


ccccc *** MEAN FREE PATH DATA (and Rayleigh cumulative prob) **********

      ! -- Re-calculate the maximum Rayleigh cumulative probability for each linear energy bin instead of the PENELOPE grid:
      call GRAaI_linear_energy(M, NBINS, EMIN, DE, PMAX_linear_energy)


      do i = 1, NBINS

        E = EMIN + (i-1)*DE             ! Set bin energy

        IF(E.LT.EABS(KPAR,M).OR.E.GT.EMAX) THEN
          WRITE(6,*) '!!ERROR!! Energy outside the table interval!',
     &               ' #bin, E = ', i, E
          STOP 'ERROR!'
        ENDIF

        E_MFP(1) = E                    ! Store the bin energy
        E_MFP(2) = PHMFP(E,KPAR,M,1)    ! Store the bin MFPs: (1) Rayleigh
        E_MFP(3) = PHMFP(E,KPAR,M,2)    ! Store the bin MFPs: (2) Compton
        E_MFP(4) = PHMFP(E,KPAR,M,3)    ! Store the bin MFPs: (3) photoelectric
        E_MFP(5) = PHMFP(E,KPAR,M,4)    ! Store the bin MFPs: (4) pair production

        E_MFP(6) = 1.0/E_MFP(2)+1.0/E_MFP(3)+1.0/E_MFP(4)+1.0/E_MFP(5)
        E_MFP(6) = 1.0/E_MFP(6)         ! Store the bin total MFP

        write(1,'(6(1x,1pe17.10))') E_MFP(1), E_MFP(2), E_MFP(3),    ! Write MFP table to external file
     &                 E_MFP(4), E_MFP(6), PMAX_linear_energy(i)     ! Write the Rayleigh cumulative probability for the energy bin

          ! E_MFP(5) --> Pair production MFP is not written bc it is not used in the simulation, but it is included in the TOTAL MFP
     

      enddo


ccccc *** RAYLEIGH DATA ***********************************************

      ! -- Rayleigh sampling data header:
      WRITE(1,'(a)')'#[RAYLEIGH INTERACTIONS (RITA sampling '//
     &              ' of atomic form factor from EPDL database)]'
      WRITE(1,'(a)')
     & '#[DATA VALUES TO SAMPLE SQUARED MOLECULAR FORM FACTOR (F^2)]'
      WRITE(1,1003) NP
      WRITE(1,'(a)')
	  !here a major modification: PENELOPE 2006 USED X INSTEAD OF Q FROM 2018
	  !we need to convert these values
	  !CGRA->CGRA03
	  !X-> QRA
	  !P-> PRA
	  !A -> ARA
	  !B -> BRA
	  !ITL -> ILTRA
	  !ITU -> ITURA
     & '#[SAMPLING DATA FROM COMMON/CGRA/: X, P, A, B, ITL, ITU]'   ! X == momentum transfer data value (adaptive grid), tabulated from the minimum to the maximum possible momentum transfers
                                                                     ! P == squared Molecular Form Factor cumulative prob at this X (adaptive grid)
                                                                     ! A & B == RITA sampling parameters
                                                                     ! ITL & ITU == lower and upper limits to speed binary search
	  																 
																	 
      do i = 1, NP
        write(1,5555) QRA(i,M)*Q2X, PRA(i,M), 
     1                ARA(i,M),BRA(i,M), ITLRA(i,M), ITURA(i,M)
      enddo
5555  format(4(1x,1pe17.10),1x,i4,1x,i4)


ccccc *** COMPTON DATA ************************************************
 
      ! -- Compton sampling data header:
      WRITE(1,'(a)')
     &  '#[COMPTON INTERACTIONS (relativistic impulse model with'//
     &  ' approximated one-electron analytical profiles)]'
      WRITE(1,'(a)')'#[NUMBER OF SHELLS]'
      WRITE(1,1003) NOSCCO(M)
      WRITE(1,'(a)')'#[SHELL INFORMATION FROM COMMON/CGCO/:'//      ! FCO == equivalent number of electrons in the shell?? (eq. 2.36 penelope 2008)
     &              ' FCO, UICO, FJ0, KZCO, KSCO]'                   ! UICO == shell ionization energy
                                                                     ! FJ0 == one-electron shell profile at p_z=0 (eq. 2.54, page 72, penelope 2008)
                                                                     ! KZCO == element that "owns" the shell??
                                                                     ! KSCO == atomic shell number, ie, atomic transition line
                                                                     ! NOSCCO == number of shells, after grouping 
      do i = 1, NOSCCO(M)
        write(1,5107) FCO(M,i), UICO(M,i), FJ0(M,i),
     &                KZCO(M,i), KSCO(M,i)
      enddo
 5107 format(3(1X,E16.8),2(1X,I4))


      WRITE(1,'(a)')' '
      CLOSE(1)

      WRITE(6,'(a)')
     & '*** Material file correctly generated. Have a nice simulation!'
      WRITE(6,*)' '
	  
	  !!!!!DEBUG!!!
	  !call GRAaD_X
	  !!!DEBUG!!!!!
     
      END
	  
	  
C  *********************************************************************
C        Code based on PENELOPE's subroutine:  SUBROUTINE GRAaR 
C        GraaI was depreciated in PENELOPE newer versions
C  *********************************************************************
      SUBROUTINE GRAaI_linear_energy(M, nbins, emin2, de, PMAX_linear_e)
	  
      USE PENELOPE_mod
      USE PENERROR_mod
	  
	  
	  IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER*2 LASYMB
      PARAMETER (REV=5.10998928D5)  ! Electron rest energy (eV)
      PARAMETER (RREV=1.0D0/REV)
	  
	  
C  ****  Energy grid and interpolation constants for the current energy.
      COMMON/CEGRID/EMIN,EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
	  
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
	 
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RSCR(99),ETA(99),EB(99,30),
     1  ALW(99,30),CP0(99,30),IFI(99,30),IKS(99,30),NSHT(99),LASYMB(99)
	 

	  PARAMETER (NM=512)
      COMMON/CRITA/QTI(NM),PACI(NM),DPACI(NM),AI(NM),BI(NM),
     1             ITLI(NM),ITUI(NM),NPM1I
	 
	  PARAMETER (NP=150)
      COMMON/CGRA03/QRA(NP,MAXMAT),PRA(NP,MAXMAT),DPRA(NP,MAXMAT),
     1  ARA(NP,MAXMAT),BRA(NP,MAXMAT),PMAX(NEGP,MAXMAT),
     2  ITLRA(NP,MAXMAT),ITURA(NP,MAXMAT)
	 
      PARAMETER (NIP=51)
	  PARAMETER (NQ=250,NEX=1024)
      DIMENSION QI(NIP),FUN(NIP),SUM(NIP)
	  COMMON/CGRA00/FACTE,Q2MAX,MM,MOM
	  COMMON/CGRA02/QQ(NQ)
	       
	  DIMENSION Q(NQ),F(NQ),FFI(NQ),ER(NEX),A(NQ),B(NQ),C(NQ),D(NQ)
	  
      EXTERNAL GRAaF2
	  
	  !! Dimension output array:
      PARAMETER (MAX_ENERGY_BINS=50000)
      DIMENSION PMAX_linear_e(MAX_ENERGY_BINS)
	  
	  MM=M
      Q2MIN=0.0D0
      Q2MAX=0.0D0
      NPT=NP
      NU=NPT/4
	  
	  
	  !q2 max debug
C	  Q2MAX = 9000000.0D0
      !routine to calculate Q2MAX
	  DO I=2,NQ
        IF(GRAaF2(EXP(QQ(I))).GT.1.0D-35) Q2MAX=EXP(QQ(I-1))
      ENDDO
	  
	  !debug
C	  WRITE(*,*) Q2MAX
	  
	  !Function GraaD1 no longer exists, changed to GraaF2
	  
C     DO I=2,NQ
C	    Q(I) = SQRT(EXP(QQ(I)))
C        IF(GRAaF2(Q(I)**2).GT.1.0D-35) Q2MAX=Q(I-1)**2
C      ENDDO
	  !debug 
	  
	  
	  CALL RITAI0(GRAaF2,Q2MIN,Q2MAX,NPT,NU,ERRM,0,IER)
	  !very important to work
	  NPI=NPM1I+1
	  
	  
C  ****  Upper limit of the X2 interval for the PENELOPE grid energies.
C
        DO IE=1,nbins
        QM=2.0D0*(emin2+(IE-1)*de)*RREV
        Q2M=QM*QM
        IF(Q2M.GT.QTI(1)) THEN
          IF(Q2M.LT.QTI(NP)) THEN
            I=1
            J=NPI
    1       IT=(I+J)/2
            IF(Q2M.GT.QTI(IT)) THEN
              I=IT
            ELSE
              J=IT
            ENDIF
            IF(J-I.GT.1) GO TO 1
			
C
            Q1=QTI(I)
            Q2=Q2M
            DQ=(Q2-Q1)/DBLE(NIP-1)
            DO K=1,NIP
              QI(K)=Q1+DBLE(K-1)*DQ
              TAU=(QI(K)-QTI(I))/(QTI(I+1)-QTI(I))
              CON1=2.0D0*BI(I)*TAU
              CI=1.0D0+AI(I)+BI(I)
              CON2=CI-AI(I)*TAU

              IF(ABS(CON1).GT.1.0D-16*ABS(CON2)) THEN
                ETAP=CON2*(1.0D0-SQRT(1.0D0-2.0D0*TAU*CON1/CON2**2))
     1              /CON1
              ELSE
                ETAP=TAU/CON2
              ENDIF
              FUN(K)=DPACI(I)*(1.0D0+(AI(I)+BI(I)*ETAP)*ETAP)**2
     1              /((1.0D0-BI(I)*ETAP*ETAP)*CI*(QTI(I+1)-QTI(I)))
            ENDDO
			!SIMPSU routine was changed to SLAG6 for integration
            CALL SLAG6(DQ,FUN,SUM,NIP)
            IF(IRETRN.NE.0) RETURN
            PMAX_linear_e(IE)=PACI(I)+SUM(NIP)
          ELSE
           PMAX_linear_e(IE)=1.0D0
          ENDIF
        ELSE
          PMAX_linear_e(IE)=PACI(1)
        ENDIF
      ENDDO
	  
	  
      RETURN
      END
	  
C  ********************************************************************* 
C  DEBUG FUNCTION
C                       FUNCTION GRAaF2X
C  *********************************************************************
      FUNCTION GRAaF2X(X2)
C
C  Squared molecular form factor, as a 
C  function of 20.6070**2*(Q*SL/REV)**2.
C
      USE PENELOPE_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NQ=250)
	  PARAMETER (X2Q=0.002354797672581918D0)
      COMMON/CGRA00/FACTE,Q2MAX,MM,MOM
      COMMON/CGRA02/QQ(NQ),AR(MAXMAT,NQ),BR(MAXMAT,NQ),CR(MAXMAT,NQ),
     1              DR(MAXMAT,NQ),FF0(MAXMAT),QQM
C
	  Q2 = X2*X2Q
      IF(Q2.LT.1.0D-9) THEN
        GRAaF2X=FF0(MM)
      ELSE IF(Q2.GT.QQM) THEN
        GRAaF2X=0.0D0
      ELSE
        QL=LOG(Q2)
        CALL FINDI(QQ,QL,NQ,I)
        F2=AR(MM,I)+QL*(BR(MM,I)+QL*(CR(MM,I)+QL*DR(MM,I)))
        GRAaF2X=EXP(F2)
      ENDIF
      RETURN
      END

C  *********************************************************************
C        DEBUG SUBROUTINE
C        Code based on PENELOPE's subroutine:  SUBROUTINE GRAaR 
C        GraaI was depreciated in PENELOPE newer versions
C  *********************************************************************
      SUBROUTINE GRAaD_X!(QRAX, PRAX, ARAX, BRAX, ILTRAX, ITURAX)
	  
	  
      USE PENELOPE_mod
      USE PENERROR_mod
	  IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
	  
C  ****  Random sampling (RITA).
      PARAMETER (NM=512)
	  PARAMETER (NQ=250)
      COMMON/CRITA/QRAX(NM),PRAX(NM),DPRAX(NM),ARAX(NM),BRAX(NM),
     1             ITLRAX(NM),ITURAX(NM),NPM1I
C
      COMMON/CGRA00/FACTE,Q2MAX,MM,MOM
	  COMMON/CGRA02/QQ(NQ)
      EXTERNAL  GRAaF2X
	  PARAMETER (Q2X=424.66493475999994D0)
	  PARAMETER (NP=150)
	  
	  
	  
	  MM=1
      X2MIN=0.0D0
      X2MAX=0.0D0
      NPT=NP
      NU=NPT/4

      X2MAX = 9000000.0D0*Q2X
	  
	  WRITE(*,*) GRAaF2X(10D0)
	  
	  
	  
	  CALL RITAI0(GRAaF2X,X2MIN,X2MAX,NPT,NU,ERRM,0,IER)
	  
	  do i = 1, NP
        write(*,*) QRAX(i), PRAX(i), 
     1                ARAX(i),BRAX(i), ITLRAX(i), ITURAX(i)
      enddo
	  
	  RETURN
      END
	  
	  
      
    

    