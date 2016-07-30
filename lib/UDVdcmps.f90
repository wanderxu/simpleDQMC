!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Two subroutines used to calculate the UDV decomposition for arbitrary N*M matrix with N >= M.
!              Including the subroutines for real and complex cases. The UDV is actually the QR decomposition and we just put 
!               the norm of all the new orthogonal vectors in U, to the diagonal elements of D matrix. 
!
!          Need to use the NAG mathematical Library.
!
! COMMENT: Common file.  
! AUTHOR:  Yuan-Yao He
! DATE:    2014-12-01
! PURPOSE: Different subroutines are introduced as following:
!
!    UDVDcmps_R --> Subroutine to calculate the UDV decomposition for real N*M matrix;
!    UDVDcmps_C --> Subroutine to calculate the UDV decomposition for complex N*M matrix.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine ____________________________________________________________________________________
!__________________________________________________________________________________________________________________________________________
	subroutine UDVdcmps_R(ND1, ND2, A, U, D, V)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  UDVdcmps_R(ND1, ND2, A, U, D, V)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to create the UDV decomposition of N*M real matrix of A.
! KEYWORDS: Calculate the UDV decomposition, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2014-12-01
! DESCRIPTION:
!
!     Calculate the UDV decomposition, real version.
!
!     Input: ND1 --> Dimension of A matrix;
!            ND2 --> Dimension of A matrix;
!            A   --> Input real N*M matrix (N>M);
!
!     Outpt: U --> Result output orthogonal vectors;
!            D --> A vector with diaognal elements;
!            V --> Upper triangular matrix.
!
! MODULES USED: RealPrecsn, StdInOutSt, MatrixMult, IniDiagMat
! SUBROUTINES CALLED: (1) NAG Library: F01QCF, F01QEF
!                     (2) MatrCmpr_R, MatrMult_R
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		integer ND1             ! Dimension of A matrix A(ND1, ND2)
		integer ND2             ! Dimension of A matrix A(ND1, ND2)
		integer NCon            ! Integer indicating whether to perform double-check for the UDV decomposition
		real(8) A(ND1, ND2)    ! Input A(ND1, ND2) matrix, Result: A=UDV
		real(8) U(ND1, ND2)    ! Result Output U(ND1, ND2) matrix
		real(8) D(ND2)         ! Result Output D(ND2, ND2) matrix --> Only output the diagonal matrix elements
		real(8) V(ND2, ND2)    ! Result Output V(ND1, ND2) matrix --> An upper triangular matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer I1     ! Loop integer 
		integer I2     ! Loop integer
		integer I3     ! Loop integer
		integer K

		integer IMax
		integer Ifail
				
		real(8) XMax
		real(8) XMean
		real(8) Z
				
		integer, allocatable :: IVPT(:)
		integer, allocatable :: IVPTM1(:)
		real(8), allocatable :: XNorm(:)
		real(8), allocatable :: VHelp(:)
		real(8), allocatable :: Theta(:)
		real(8), allocatable :: Work (:)
				
		real(8), allocatable :: Tmp(:, :)
		real(8), allocatable :: V1(:, :)
				
		real(8), allocatable :: Test(:, :)
		real(8), allocatable :: Test1(:, :)
		real(8), allocatable :: Test2(:, :)
!______________________________________________________________________________________________________________	  
!__________________________________ Allocate Array and Initializations ________________________________________
!______________________________________________________________________________________________________________
		allocate(XNorm (ND2))
		allocate(VHelp (ND2))
		allocate(IVPT  (ND2))
		allocate(IVPTM1(ND2))
		allocate(Work  (ND2))
		allocate(Theta (ND2))
		allocate(Tmp(ND1, ND2))
		allocate(V1 (ND2, ND2))
		V1 = 0.0_8
				
		U = 0.0_8
		D = 0.0_8
		V = 0.0_8
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Creating matrix _________________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!________________ 1. Calculate the summation of absolute values of every column in A  ___________
!________________________________________________________________________________________________			
		do I1 = 1, ND2
			XNorm(I1) = 0.0_8
			do I2 = 1, ND1
				XNorm(I1) = XNorm(I1) + abs(A(I2, I1))
			enddo
		enddo
		do I1 = 1, ND2
			VHelp(I1) = XNorm(I1)
		enddo
!________________________________________________________________________________________________	  
!________________ 2. Pivoting according to absolute value summation in VHelp ____________________
!________________________________________________________________________________________________					
		do I1 = 1, ND2
			XMax = 0.0_8
			do I2 = 1, ND2
				if(VHelp(I2) .gt. XMax) then
					IMax = I2
					XMax = VHelp(I2)
				end if 
			enddo
			VHelp(IMax) = -1.0_8
			IVPTM1(IMax) = I1
			IVPT(I1) = IMax
		enddo
!________________________________________________________________________________________________	  
!________________ 3. Perform a scailing for the A matrix by the XNorm ___________________________
!________________________________________________________________________________________________				
		do I1 = 1, ND2
			Z = 1.0_8 / XNorm(IVPT(I1))
			K = IVPT(I1)
			do I2 = 1, ND1
				Tmp(I2, I1) = A(I2, K) * Z
			enddo
		enddo
!________________________________________________________________________________________________	  
!________________ 4. Call NAG subroutine for QR decomposition calculations ______________________
!________________________________________________________________________________________________				
		IFail = 0
		call F01QCF(ND1, ND2, Tmp, ND1, Theta, IFail)
!________________________________________________________________________________________________	  
!________________ 5. Get the D vector as diagonal elements of QR output Tmp matrix ______________
!________________________________________________________________________________________________				
		do I1 = 1, ND2
			D(I1) = abs(Tmp(I1, I1))
		enddo
!________________________________________________________________________________________________	  
!________________ 6. Keep the upper triangular part of Tmp matrix for restoring R matrix ________
!________________________________________________________________________________________________				
		do I1 = 1, ND2
			Z = 1.0_8 / D(I1)
			do I2 = I1, ND2
				V1(I1, I2) = Tmp(I1, I2) * Z
			enddo
		enddo
!________________________________________________________________________________________________	  
!________________ 7. Restore the Q matrix from the Tmp matrix in QR decomposition _______________
!________________________________________________________________________________________________				
		IFail = 0
		call F01QEF("Separate", ND1, ND2, ND2, Tmp, ND1, Theta, Work, IFail)
		do I1 = 1, ND1
			do I2 = 1, ND2
				U(I1, I2) = Tmp(I1, I2)
			enddo
		enddo
!________________________________________________________________________________________________	  
!________________ 8. Rescailing for the D vector and V1 matrix by XNorm factors _________________
!________________________________________________________________________________________________					
		do I1 = 1, ND2
			D(I1) = D(I1) * XNorm(IVPT(I1))
		enddo
		do I1 = 1, ND2-1
			Z = 1.0_8 / XNorm(IVPT(I1))
			do I2 = I1+1, ND2
				V1(I1, I2) = V1(I1, I2) * XNorm(IVPT(I2)) * Z
			enddo
		enddo
!________________________________________________________________________________________________	  
!________________ 9. Get the Upper triangulare V matrix in the UDV decomposition ________________
!________________________________________________________________________________________________				
		do I2 = 1, ND2
			do I1 = 1, ND2
				V(I1, I2) = V1(I1, IVPTM1(I2))
			enddo
		enddo
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		deallocate(XNorm )
		deallocate(VHelp )
		deallocate(IVPT  )
		deallocate(IVPTM1)
		deallocate(Work  )
		deallocate(Theta )
		deallocate(Tmp   )
		deallocate(V1    )
				
	end subroutine
!__________________________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ______________________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine ____________________________________________________________________________________
!__________________________________________________________________________________________________________________________________________
	subroutine UDVdcmps_C(ND1, ND2, A, U, D, V)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
! PROGRAM:  UDVdcmps_C(ND1, ND2, A, U, D, V)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to create the UDV decomposition of N*M real matrix of A.
! KEYWORDS: Calculate the UDV decomposition, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2014-12-01
! DESCRIPTION:
!
!     Calculate the UDV decomposition, complex version.
!
!     Input: ND1 --> Dimension of A matrix;
!            ND2 --> Dimension of A matrix;
!            A   --> Input real N*M matrix (N>M);
!
!     Outpt: U --> Result output orthogonal vectors;
!            D --> A vector with diaognal elements;
!            V --> Upper triangular matrix.
!            NCon --> Show if the calculation is good or not.
!
! MODULES USED: RealPrecsn, StdInOutSt, MatrixMult, IniDiagMat
! SUBROUTINES CALLED: (1) NAG Library: F01RCF, F01REF
!                     (2) MatrCmpr_C, MatrMult_C
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		integer ND1                ! Dimension of A matrix A(ND1, ND2)
		integer ND2                ! Dimension of A matrix A(ND1, ND2)
		integer NCon               ! Integer indicating whether to perform double-check for the UDV decomposition
		complex(8) A(ND1, ND2)    ! Input A(ND1, ND2) matrix, Result: A=UDV
		complex(8) U(ND1, ND2)    ! Result Output U(ND1, ND2) matrix
		real(8) D(ND2)            ! Result Output D(ND2, ND2) matrix --> Only output the diagonal matrix elements
		complex(8) V(ND2, ND2)    ! Result Output V(ND1, ND2) matrix --> An upper triangular matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer I1     ! Loop integer 
		integer I2     ! Loop integer
		integer I3     ! Loop integer
				
		integer Ifail
				
		real(8) DetV
		real(8) X
				
		real(8) XMax
		real(8) XMean
				
		complex(8) Z

		complex(8), allocatable :: Theta(:)
		complex(8), allocatable :: Work (:)
				
		complex(8), allocatable :: Tmp(:, :)
				
		complex(8), allocatable :: Test(:, :)
		complex(8), allocatable :: Test1(:, :)
		complex(8), allocatable :: Test2(:, :)
!______________________________________________________________________________________________________________	  
!__________________________________ Allocate Array and Initializations ________________________________________
!______________________________________________________________________________________________________________
		allocate(  Tmp(ND1, ND2))
		allocate( Work(ND2))
		allocate(Theta(ND2))
		
		U = cmplx(0.0_8, 0.0_8, 8)
		D = cmplx(0.0_8, 0.0_8, 8)
		V = cmplx(0.0_8, 0.0_8, 8)
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Creating matrix _________________________________________
!______________________________________________________________________________________________________________
!________________________________________________________________________________________________	  
!________________ 1. Call NAG subroutine for QR decomposition calculations ______________________
!________________________________________________________________________________________________
		Tmp   = A
		IFail = 0
		call F01RCF(ND1, ND2, Tmp, ND1, Theta, IFail)
!________________________________________________________________________________________________	  
!________________ 2. Keep the upper triangular part of Tmp matrix for restoring R matrix ________
!____________________________ and the determinant of R matrix ___________________________________
!________________________________________________________________________________________________
		do I1 = 1, ND2
			do I2 = I1, ND2
				V(I1, I2) = Tmp(I1, I2)
			enddo
		enddo
		DetV = 1.0_8
		do I2 = 1, ND2
			DetV = DetV * dble(Tmp(I2, I2))
		enddo
!________________________________________________________________________________________________	  
!________________ 4. Restore the Q matrix from the Tmp matrix in QR decomposition _______________
!________________________________________________________________________________________________				
		call F01REF('Separate', ND1, ND2, ND2, Tmp, ND1, Theta, Work, IFail)
		do I1 = 1, ND1
			do I2 = 1, ND2
				U(I1, I2) = Tmp(I1, I2)
			enddo
		enddo
!________________________________________________________________________________________________	  
!________________ 5. Rescailing for the D vector and V1 matrix by XNorm factors _________________
!________________________________________________________________________________________________					
		if( dble(DetV) .lt. 0.0_8 ) then
			do I1 = 1, ND1
				U(I1, 1) = -U(I1, 1)
			enddo
			do I2 = 1, ND2
				V(1, I2) = -V(1, I2)
			enddo
		end if
!________________________________________________________________________________________________	  
!________________ 6. Get the Upper triangulare V matrix in the UDV decomposition ________________
!________________________________________________________________________________________________				
		do I2 = 1, ND2
			D(I2) = abs(dble(V(I2, I2)))
		enddo
		do I1 = 1, ND2
			Z = cmplx(1.0_8, 0.0_8, 8) / cmplx( D(I1), 0.0_8, 8 )
			do I2 = I1, ND2
				V(I1, I2) = V(I1, I2) * Z
			enddo
		enddo
!______________________________________________________________________________________________________________	  
!___________________________________________ Deallocate the arrays ____________________________________________
!______________________________________________________________________________________________________________
		deallocate(Work)
		deallocate(Theta)
		deallocate(Tmp)
				
	end subroutine
!__________________________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ______________________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$    
