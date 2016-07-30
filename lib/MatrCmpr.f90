!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Two subroutines used to construct compare two different matrices and calculate the difference of these two
!              matrices. Calculate the largest difference element and calculate the average difference for all elements.
! COMMENT: Common file.  
! AUTHOR:  Yuan-Yao He
! DATE:    2014-12-01
! PURPOSE: Different subroutines are introduced as following:
!
!    DiagMatr_R --> Subroutine to create real diagonal matrix;
!    DiagMatr_C --> Subroutine to create complex diagonal matrix.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine ____________________________________________________________________________________
!__________________________________________________________________________________________________________________________________________
	subroutine MatrCmpr_R(N, M, A, B, XMax, XMean)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  MatrCmpr_R(N, M, A, B, XMax, XMean)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to compare two real different matrices and calculate their differences in elements.
! KEYWORDS: Compare two matrices, real version.
! AUTHOR:   Yuan-Yao He
! TIME:     2014-12-01
! DESCRIPTION:
!
!     Compare two matrices, real version. 
!
!     Input: N --> Dimension of A, B matrix as A(N, M) and B(N, M);
!            M --> Dimension of A, B matrix as A(N, M) and B(N, M);
!            A --> Input matrix A(N, M);
!            B --> Input matrix A(N, M);
!     
!     Outpt: XMax  --> The largest absolute difference;
!            XMean --> The average difference in all elements.
!
! MODULES USED: RealPrecsn
! SUBROUTINES CALLED:  (none)
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		integer N           ! Dimension of A, B matrices
		integer M           ! Dimension of A, B matrices
		real(8) XMax       ! Maximum difference for matrix elements in A, B matrices
		real(8) XMean      ! Average difference for matrix elements in A, B matrices
		real(8) A(N, M)    ! Input A matrix
		real(8) B(N, M)    ! Input B matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer I1      ! Loop integer 
		integer I2      ! Loop integer
		real(8) Diff   ! Difference between the elements in A and B
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Creating matrix _________________________________________
!______________________________________________________________________________________________________________
		XMax  = 0.0_8
		XMean = 0.0_8
				
		do I1 = 1, N
			do I2 = 1, M
				Diff = abs(B(I1, I2) - A(I1, I2))
				if(Diff .gt. XMax) then
					XMax = Diff
				end if
				XMean = XMean + Diff
			enddo
		enddo
				
		XMean = XMean / (M * N * 1.0_8)
		
	end subroutine
!__________________________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ______________________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine ____________________________________________________________________________________
!__________________________________________________________________________________________________________________________________________
	subroutine MatrCmpr_C(N, M, A, B, XMax, XMean)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
! PROGRAM:  MatrCmpr_C(N, M, A, B, XMax, XMean)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to compare two different complex matrices and calculate their differences in elements.
! KEYWORDS: Compare two matrices, complex version.
! AUTHOR:   Yuan-Yao He
! TIME:     2014-12-01
! DESCRIPTION:
!
!     Compare two matrices, complex version. 
!
!     Input: N --> Dimension of A, B matrix as A(N, M) and B(N, M);
!            M --> Dimension of A, B matrix as A(N, M) and B(N, M);
!            A --> Input matrix A(N, M);
!            B --> Input matrix A(N, M);
!     
!     Outpt: XMax  --> The largest absolute difference;
!            XMean --> The average difference in all elements.
!
! MODULES USED: RealPrecsn
! SUBROUTINES CALLED:  (none)
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 			
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		integer N              ! Dimension of A, B matrices
		integer M              ! Dimension of A, B matrices
		real(8) XMax          ! Maximum difference for matrix elements in A, B matrices
		real(8) XMean         ! Average difference for matrix elements in A, B matrices
		complex(8) A(N, M)    ! Input A matrix
		complex(8) B(N, M)    ! Input B matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer I1      ! Loop integer 
		integer I2      ! Loop integer
		real(8) Diff   ! Difference between the elements in A and B
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of Creating matrix _________________________________________
!______________________________________________________________________________________________________________			
		XMax  = 0.0_8
		XMean = 0.0_8
				
		do I1 = 1, N
			do I2 = 1, M
				Diff = sqrt( real( (B(I1, I2) - A(I1, I2)) * conjg(B(I1, I2) - A(I1, I2)) ) )
				if(Diff .gt. XMax) then
					XMax = Diff
				end if
				XMean = XMean + Diff
			enddo
		enddo
				
		XMean = XMean / (M * N * 1.0_8)
		
	end subroutine
!__________________________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ______________________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 