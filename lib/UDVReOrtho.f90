!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PROGRAM: Subroutines used to carry out the reorthogonalization process for both complex and real matrices, by applying the UDV
!                decomposition.
! COMMENT: Performing the reorthogonalization process. 
! AUTHOR:  Yuan-Yao He
! DATE:    2016-05-03
! PURPOSE: Different subroutines are introduced as following:
!             
!   UDVOrthg_C --> Subroutine to perform the UDV decomposition for complex matrix;
!   UDVOrthg_R --> Subroutine to perform the UDV decomposition for real    matrix.
!             
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine ____________________________________________________________________________________
!__________________________________________________________________________________________________________________________________________
	subroutine UDVOrthg_C(ND1, ND2, A) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  UDVOrthg_C(ND1, ND2, A)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the UDV decomposition for input A matrix as: A=UDV and give U to A: A=U.
!                 The UDV decomposition is actually a QR decomposition as reorthogonal process.
! KEYWORDS: UDV decomposition for A matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2016-05-03
! DESCRIPTION:
!
!     UDV decomposition for A matrix as A=UDV and finally A=U.
!         This decomposition is done by calling subroutines defined in UDVDecomps module.
!
!     Input:  (none)   Output: (none)
!
! MODULES USED: RealPrecsn
! SUBROUTINES CALLED: UDVdcmps_C
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		integer ND1               ! dimension of A(ND1, ND2)
		integer ND2               ! dimension of A(ND1, ND2)
		complex(8) A(ND1, ND2)   ! Input complex matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer I1     ! Loop integer
		integer I2     ! Loop integer
		
		complex(8), allocatable ::  D(:)      ! The D vector in UDV decomposition
		complex(8), allocatable ::  U(:, :)   ! The U matrix in UDV decomposition
		complex(8), allocatable ::  V(:, :)   ! The V matrix in UDV decomposition
		complex(8), allocatable :: AT(:, :)   ! The transpose matrix of A matrix
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of UDV decomposition _______________________________________
!______________________________________________________________________________________________________________
!______________________________________________________________________________________________	  
!_____________________ 0. For the ND1 > ND2 case ______________________________________________
!______________________________________________________________________________________________
		if(ND1 .gt. ND2) then		
			allocate(U(ND1, ND2))
			allocate(D(ND2))
			allocate(V(ND2, ND2))
			call UDVdcmps_C(ND1, ND2, A, U, D, V)
			A = U
			deallocate(D)
			deallocate(U)
			deallocate(V)
!______________________________________________________________________________________________	  
!_____________________ 1. For the ND1 < ND2 case ______________________________________________
!______________________________________________________________________________________________
		else
			allocate(U(ND2, ND1))
			allocate(D(ND1))
			allocate(V(ND1, ND1))
			allocate(AT(ND2, ND1))
			do I1 = 1, ND1
				do I2 = 1, ND2
					AT(I2, I1) = A(I1, I2)
				enddo
			enddo
			call UDVdcmps_C(ND2, ND1, AT, U, D, V)
			do I1 = 1, ND1
				do I2 = 1, ND2
					A(I1, I2) = U(I2, I1)
				enddo
			enddo
			deallocate(D)
			deallocate(U)
			deallocate(V)
			deallocate(AT)
		end if
		
   end subroutine
!__________________________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ______________________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   
   
   
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine ____________________________________________________________________________________
!__________________________________________________________________________________________________________________________________________
	subroutine UDVOrthg_R(ND1, ND2, A) 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
! PROGRAM:  UDVOrthg_R(ND1, ND2, A)
! TYPE:     subroutine
! PURPOSE:  This Subroutine is used to perform the UDV decomposition for input A matrix as: A=UDV and give U to A: A=U.
!                 The UDV decomposition is actually a QR decomposition as reorthogonal process.
! KEYWORDS: UDV decomposition for A matrix.
! AUTHOR:   Yuan-Yao He
! TIME:     2016-05-03
! DESCRIPTION:
!
!     UDV decomposition for A matrix as A=UDV and finally A=U.
!         This decomposition is done by calling subroutines defined in UDVDecomps module.
!
!     Input:  (none)   Output: (none)
!
! MODULES USED: RealPrecsn
! SUBROUTINES CALLED: UDVdcmps_R
! END PROGRAM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!______________________________________________________________________________________________________________	  
!_________________________________________ Modules used in this subroutine ____________________________________
!______________________________________________________________________________________________________________ 
		implicit none
!______________________________________________________________________________________________________________	  
!_________________________________________ All Input and Output Quantities ____________________________________
!______________________________________________________________________________________________________________
		integer ND1               ! dimension of A(ND1, ND2)
		integer ND2               ! dimension of A(ND1, ND2)
		real(8) A(ND1, ND2)      ! Input real matrix
!______________________________________________________________________________________________________________	  
!______________________________ Temperory Quantities used in the calculations _________________________________
!______________________________________________________________________________________________________________
		integer I1     ! Loop integer
		integer I2     ! Loop integer
		
		real(8), allocatable ::  D(:)      ! The D vector in UDV decomposition
		real(8), allocatable ::  U(:, :)   ! The U matrix in UDV decomposition
		real(8), allocatable ::  V(:, :)   ! The V matrix in UDV decomposition
		real(8), allocatable :: AT(:, :)   ! The transpose matrix of A matrix
!______________________________________________________________________________________________________________	  
!_______________________________ Main calculations of UDV decomposition _______________________________________
!______________________________________________________________________________________________________________
!______________________________________________________________________________________________	  
!_____________________ 0. For the ND1 > ND2 case ______________________________________________
!______________________________________________________________________________________________
		if(ND1 .gt. ND2) then		
			allocate(U(ND1, ND2))
			allocate(D(ND2))
			allocate(V(ND2, ND2))
			call UDVdcmps_R(ND1, ND2, A, U, D, V)
			A = U
			deallocate(D)
			deallocate(U)
			deallocate(V)
!______________________________________________________________________________________________	  
!_____________________ 1. For the ND1 < ND2 case ______________________________________________
!______________________________________________________________________________________________
		else
			allocate(U(ND2, ND1))
			allocate(D(ND1))
			allocate(V(ND1, ND1))
			allocate(AT(ND2, ND1))
			do I1 = 1, ND1
				do I2 = 1, ND2
					AT(I2, I1) = A(I1, I2)
				enddo
			enddo
			call UDVdcmps_R(ND2, ND1, AT, U, D, V)
			do I1 = 1, ND1
				do I2 = 1, ND2
					A(I1, I2) = U(I2, I1)
				enddo
			enddo
			deallocate(D)
			deallocate(U)
			deallocate(V)
			deallocate(AT)
		end if

	end subroutine
!__________________________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ______________________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$