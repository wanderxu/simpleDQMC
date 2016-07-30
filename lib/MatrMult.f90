!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Two subroutines used to construct the interface for matrix multiplication, including the real matrices
!              multiplication and complex matrices multiplication. The kernel subroutine used to perform the matrix multiplication is
!              subroutines from the BLAS mathematical library. 
!
!          Need BLAS mathematical Library.
!
! COMMENT: Common file.  
! AUTHOR:  Yuan-Yao He
! DATE:    2014-12-01
! PURPOSE: Different subroutines are introduced as following:
!             
!    MatrMult_R --> Subroutine to perform the real matrices multiplication as C=A*B;
!    MatrMult_C --> Subroutine to perform the complex matrices multiplication as C=A*B.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine ____________________________________________________________________________________
!__________________________________________________________________________________________________________________________________________
	subroutine MatrMult_R(N, M, K, A, B, C)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  MatrMult_R(N, M, K, A, B, C)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the real matrices multiplization  as: C = A * B
! KEYWORDS: Matrice multiplication, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2014-11-17
! DESCRIPTION:
!
!     Matrice multiplication, real version. 
!
!     Input: N --> Dimension of A matrix as A(N, M);
!            N --> Dimension of A matrix as A(N, M), and B matrix as B(M, K);
!            K --> Dimension of B matrix as B(M, K);
!            A --> Input matrix A(N, M);
!            B --> Input matrix B(M, K).
!             
!     Outpt: C --> Result output matrix as C(N, K).
!
! MODULES USED: RealPrecsn
! SUBROUTINES CALLED: (1) BLAS Library: DGEMM
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		integer N            ! Dimension of A matrix A(N, M)
		integer M            ! Dimension of A matrix A(N, M)
		integer K            ! Dimension of B matrix B(M, K) 
		real(8) A(N, M)     ! Input A matrix
		real(8) B(M, K)     ! Input B matrix
		real(8) C(N, K)     ! Ouput C matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer LDA   ! leading dimension of A matrix
		integer LDB   ! leading dimension of B matrix
		integer LDC   ! leading dimension of C matrix
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of the matrix multiplication _______________________________
!______________________________________________________________________________________________________________
		LDA = N
		LDB = M
		LDC = N
				
		call DGEMM("N", "N", N, K, M, 1.0_8, A, LDA, B, LDB, 0.0_8, C, LDC) ! Call BLAS subroutine
		
	end subroutine
!__________________________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ______________________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine ____________________________________________________________________________________
!__________________________________________________________________________________________________________________________________________
	subroutine MatrMult_C(N, M, K, A, B, C)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  MatrMult_C(N, M, K, A, B, C)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the complex matrices multiplization  as: C = A * B
! KEYWORDS: Matrice multiplication, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2014-11-17
! DESCRIPTION:
!
!     Matrice multiplication, complex version. 
!
!     Input: N --> Dimension of A matrix as A(N, M);
!            N --> Dimension of A matrix as A(N, M), and B matrix as B(M, K);
!            K --> Dimension of B matrix as B(M, K);
!            A --> Input matrix A(N, M);
!            B --> Input matrix B(M, K).
!             
!     Outpt: C --> Result output matrix as C(N, K).
!
! MODULES USED: RealPrecsn
! SUBROUTINES CALLED: (1) BLAS Library: ZGEMM
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		integer N            ! Dimension of A matrix A(N, M)
		integer M            ! Dimension of A matrix A(N, M)
		integer K            ! Dimension of B matrix B(M, K) 
		complex(8) A(N, M)  ! Input A matrix
		complex(8) B(M, K)  ! Input B matrix
		complex(8) C(N, K)  ! Ouput C matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer LDA   ! leading dimension of A matrix
		integer LDB   ! leading dimension of B matrix
		integer LDC   ! leading dimension of C matrix
		
		complex(8) Ztp1  ! Complex temporary number used in calculations
		complex(8) Ztp2  ! Complex temporary number used in calculations
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of the matrix multiplication _______________________________
!______________________________________________________________________________________________________________
		LDA = N
		LDB = M
		LDC = N
		Ztp1 = cmplx(1.0_8, 0.0_8, 8)
		Ztp2 = cmplx(0.0_8, 0.0_8, 8)
				
		call ZGEMM("N", "N", N, K, M, Ztp1, A, LDA, B, LDB, Ztp2, C, LDC) ! Call BLAS subroutine
		
	end subroutine
!__________________________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ______________________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$