!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! PROGRAM: A few subroutines from slatec which were modified to solve eigensystems in general complex Hermitian arithmetics.
! COMMENT: to learn about netlib, send an otherwise empty message to netlib@research.att.com
!          containing 'send index' in the subject header)
!          The WWW address of netlib is http://netlib.att.com/netlib/search.html
! AUTHOR:  Yuan-Yao He  -->  661
! DATE:    2014-03-22
! PURPOSE: Different subroutines are introduced as following:
!             
!   HermDiag1(NM, N, AR, AI, W, IERR)         --> Diagonalizing complex Hermitian Matrix AR(N, N) + i * AI(N, N), Only EigVal;
!   HermDiag2(NM, N, AR, AI, W, ZR, ZI, IERR) --> Diagonalizing complex Hermitian Matrix AR(N, N) + i * AI(N, N), Both EigVal and EigVec;
!   HermDiag3(NM, N, Amat, W, Zmat )          --> Diagonalizing complex Hermitian Matrix Amat(N,N), Both EigVal and EigVec;
!   HTRIDI(NM,N,AR,AI,D,E,TAU)                --> Reduce a complex Hermitian matrix to a real symmetric tridiagonal matrix;
!   RTQL2(NM,N,D,E,Z,IERR)                    --> Compute the eigenvalues and eigenvectors of a real symmetric tridiagonal matrix;
!   HTRIBK(NM,N,AR,AI,TAU,M,ZR,ZI)            --> Form the eigenvectors of a complex Hermitian matrix from the eigenvectors of 
!                                                              the real symmetric tridiagonal matrix output from HTRIDI.  
!   PyThagC                                   --> sqrt(a** 2+b**2) without overflow or destructive underflow
!
! SUBROUTINES CALLED:  Function  : PyThagC
! PURPOSE: Get all the Eigenvalues and Eigenvectors for general Complex Hermitian Matrix
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine ____________________________________________________________________________________
!__________________________________________________________________________________________________________________________________________   
	subroutine HermDiag1(NM, N, AR, AI, W, IERR)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
! PROGRAM:  HermDiag1(NM, N, AR, AI, W, IERR)
! TYPE:     subroutine
! PURPOSE:  Compute the eigenvalues and, optionally, the eigenvectors of a complex Hermitian matrix.
! LIBRARY:  SLATEC (EISPACK)
! CATEGORY: D4A3
! KEYWORDS: EIGENVALUES, EIGENVECTORS, EISPACK
! AUTHOR:   Smith, B. T., et al.
! DESCRIPTION
!
!     This subroutine calls the recommended sequence of subroutines from the eigensystem subroutine package (EISPACK)
!     to find the eigenvalues and eigenvectors (if desired) of a COMPLEX HERMITIAN matrix.
!
!     On INPUT
!
!        NM: must be set to the row dimension of the two-dimensional array parameters, AR, AI, ZR and ZI, as declared in the
!            calling program dimension statement.  NM is an INTEGER variable.
!        N:  is the order of the matrix A=(AR,AI).  N is an INTEGER variable.  N must be less than or equal to NM.
!        AR and AI contain the real and imaginary parts, respectively,of the complex Hermitian matrix.  AR and AI are two-dimensional
!            REAL arrays, dimensioned AR(NM,N) and AI(NM,N).
!
!     On OUTPUT
!
!        W contains the eigenvalues in ascending order. W is a one-dimensional REAL array, dimensioned W(N).
!        ZR and ZI contain the real and imaginary parts, respectively, of the eigenvectors.  ZR and ZI are two-dimensional REAL 
!              arrays, dimensioned ZR(NM,N) and ZI(NM,N).
!        IERR is an INTEGER flag set to
!               Zero       for normal return,
!               10*N       if N is greater than NM,
!               J          if the J-th eigenvalue has not been determined after a total of 30 iterations. The eigenvalues should be 
!                                correct for indices 1, 2, ..., IERR-1, but no eigenvectors are computed.
!        FV1 and FV2 are one-dimensional REAL arrays used for temporary storage, dimensioned FV1(N) and FV2(N).
!        FM1 is a two-dimensional REAL array used for temporary storage, dimensioned FM1(2,N).

!     Questions and comments should be directed to B. S. Garbow, APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
! REFERENCES:  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,1976.
! ROUTINES CALLED  HTRIBK, HTRIDI, TQL2, TQLRAT
! REVISION HISTORY  (YYMMDD)
!      760101  DATE WRITTEN
!      890831  Modified array declarations.  (WRB)
!      890831  REVISION DATE from Version 3.2
!      891214  Prologue converted to Version 4.0 format.  (BAB)
!      920501  Reformatted the REFERENCES section.  (WRB)
! END PROLOGUE  CH
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
		implicit none

		integer NM,N,IERR,i
		real(8) W(N)
		real(8) AR(NM,N),AI(NM,N)
		real(8), allocatable :: FV1(:)
		real(8), allocatable :: FM1(:,:)
		real(8), allocatable :: ZR(:, :)
		real(8), allocatable :: ZI(:, :)

		allocate(FV1(N))
		allocate(FM1(2, N))
		allocate(ZR(NM, N))
		allocate(ZI(NM, N))
	  
		call HTRIDI(NM,N,AR,AI,W,FV1,FM1)

!%%%%%%%%%%%%%%%%%%%% FIND BOTH EIGENVALUES AND EIGENVECTORS%%%%%%%%%%%%%%%%%%
		ZR(1:N,1:N)=0.0_8
		do i=1,N
			ZR(i,i)=1.0_8
		enddo
		CALL  TQL2(NM,N,W,FV1,ZR,IERR)
		CALL  HTRIBK(NM,N,AR,AI,FM1,N,ZR,ZI)
	  
		deallocate(FV1)
		deallocate(FM1)
		deallocate(ZR)
		deallocate(ZI)

	end subroutine
!__________________________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ______________________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine ____________________________________________________________________________________
!__________________________________________________________________________________________________________________________________________   
	subroutine HermDiag2(NM, N, AR, AI, W, ZR, ZI, IERR)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
! PROGRAM:  HermDiag2(NM, N, AR, AI, W, ZR, ZI, IERR)
! TYPE:     subroutine
! PURPOSE:  Compute the eigenvalues and, optionally, the eigenvectors of a complex Hermitian matrix.
! LIBRARY:  SLATEC (EISPACK)
! CATEGORY: D4A3
! KEYWORDS: EIGENVALUES, EIGENVECTORS, EISPACK
! AUTHOR:   Smith, B. T., et al.
! DESCRIPTION
!
!     This subroutine calls the recommended sequence of subroutines from the eigensystem subroutine package (EISPACK)
!     to find the eigenvalues and eigenvectors (if desired) of a COMPLEX HERMITIAN matrix.
!
!     On INPUT
!
!        NM: must be set to the row dimension of the two-dimensional array parameters, AR, AI, ZR and ZI, as declared in the
!            calling program dimension statement.  NM is an INTEGER variable.
!        N:  is the order of the matrix A=(AR,AI).  N is an INTEGER variable.  N must be less than or equal to NM.
!        AR and AI contain the real and imaginary parts, respectively,of the complex Hermitian matrix.  AR and AI are two-dimensional
!            REAL arrays, dimensioned AR(NM,N) and AI(NM,N).
!
!     On OUTPUT
!
!        W contains the eigenvalues in ascending order. W is a one-dimensional REAL array, dimensioned W(N).
!        ZR and ZI contain the real and imaginary parts, respectively, of the eigenvectors.  ZR and ZI are two-dimensional REAL 
!              arrays, dimensioned ZR(NM,N) and ZI(NM,N).
!        IERR is an INTEGER flag set to
!               Zero       for normal return,
!               10*N       if N is greater than NM,
!               J          if the J-th eigenvalue has not been determined after a total of 30 iterations. The eigenvalues should be 
!                                correct for indices 1, 2, ..., IERR-1, but no eigenvectors are computed.
!        FV1 and FV2 are one-dimensional REAL arrays used for temporary storage, dimensioned FV1(N) and FV2(N).
!        FM1 is a two-dimensional REAL array used for temporary storage, dimensioned FM1(2,N).

!     Questions and comments should be directed to B. S. Garbow, APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
! REFERENCES:  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,1976.
! ROUTINES CALLED  HTRIBK, HTRIDI, TQL2, TQLRAT
! REVISION HISTORY  (YYMMDD)
!      760101  DATE WRITTEN
!      890831  Modified array declarations.  (WRB)
!      890831  REVISION DATE from Version 3.2
!      891214  Prologue converted to Version 4.0 format.  (BAB)
!      920501  Reformatted the REFERENCES section.  (WRB)
! END PROLOGUE  CH
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
		implicit none

		integer NM,N,IERR,i
		real(8) W(N)
		real(8) ZR(NM,N),ZI(NM,N)
		real(8) AR(NM,N),AI(NM,N)
		real(8),allocatable :: FV1(:)
		real(8),allocatable :: FM1(:,:)

		allocate(FV1(N))
		allocate(FM1(2,N))
	  
		call HTRIDI(NM,N,AR,AI,W,FV1,FM1)

!%%%%%%%%%%%%%%%%%%%% FIND BOTH EIGENVALUES AND EIGENVECTORS%%%%%%%%%%%%%%%%%%
		ZR(1:N,1:N)=0.0_8
		do i=1,N
			ZR(i,i)=1.0_8
		enddo
		CALL  TQL2(NM,N,W,FV1,ZR,IERR)
		CALL  HTRIBK(NM,N,AR,AI,FM1,N,ZR,ZI)
	  
		deallocate(FV1)
		deallocate(FM1)

	end subroutine
!__________________________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ______________________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine ____________________________________________________________________________________
!__________________________________________________________________________________________________________________________________________   
	subroutine HermDiag3(NM, N, Amat, W, Zmat )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
! PROGRAM:  HermDiag2(NM, N, Amat, W, Zmat )
! TYPE:     subroutine
! PURPOSE:  Compute the eigenvalues and, optionally, the eigenvectors of a complex Hermitian matrix.
! LIBRARY:  SLATEC (EISPACK)
! CATEGORY: D4A3
! KEYWORDS: EIGENVALUES, EIGENVECTORS, EISPACK
! AUTHOR:   Smith, B. T., et al.
! DESCRIPTION
!
!     This subroutine calls the recommended sequence of subroutines from the eigensystem subroutine package (EISPACK)
!     to find the eigenvalues and eigenvectors (if desired) of a COMPLEX HERMITIAN matrix.
!
!     On INPUT
!
!        NM: must be set to the row dimension of the two-dimensional array parameters, Amat, Zmat, AR, AI, ZR, ZI as declared in the
!            calling program dimension statement.  NM is an INTEGER variable.
!        N:  is the order of the matrix A=(AR,AI).  N is an INTEGER variable.  N must be less than or equal to NM.
!        Amat: the complex Hermitian matrix.  Amat are two-dimensional COMPLEX arrays, dimensioned Amat(NM,N)
!
!     On OUTPUT
!
!        W contains the eigenvalues in ascending order. W is a one-dimensional REAL array, dimensioned W(N).
!        Zmat: the eigenvectors.  ZR and ZI are two-dimensional COMPLEX arrays, dimensioned Zmat(NM,N).
!
!        AR and AI are real and imaginary part of Amat
!        ZR and ZI are real and imaginary part of Zmat

!        FV1 and FV2 are one-dimensional REAL arrays used for temporary storage, dimensioned FV1(N) and FV2(N).
!        FM1 is a two-dimensional REAL array used for temporary storage, dimensioned FM1(2,N).

!     Questions and comments should be directed to B. S. Garbow, APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
! REFERENCES:  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,1976.
! ROUTINES CALLED  HTRIBK, HTRIDI, TQL2, TQLRAT
! REVISION HISTORY  (YYMMDD)
!      760101  DATE WRITTEN
!      890831  Modified array declarations.  (WRB)
!      890831  REVISION DATE from Version 3.2
!      891214  Prologue converted to Version 4.0 format.  (BAB)
!      920501  Reformatted the REFERENCES section.  (WRB)
! END PROLOGUE  CH
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
		implicit none

		integer NM,N,IERR,i, j
		real(8) W(N)
        complex(8) Amat(NM,N), Zmat(NM,N)
		real(8),allocatable :: FV1(:)
		real(8),allocatable :: FM1(:,:)
        real(8),allocatable :: AR(:,:), AI(:,:)
        real(8),allocatable :: ZR(:,:), ZI(:,:)

		allocate(FV1(N))
		allocate(FM1(2,N))
        allocate( AR(NM,N),AI(NM,N) )
        allocate( ZR(NM,N),ZI(NM,N) )

       do j = 1, N
           do i = 1, NM
               AR(i,j) =  real( Amat(i,j) )
               AI(i,j) = aimag( Amat(i,j) )
           end do
       end do
	  
		call HTRIDI(NM,N,AR,AI,W,FV1,FM1)

!%%%%%%%%%%%%%%%%%%%% FIND BOTH EIGENVALUES AND EIGENVECTORS%%%%%%%%%%%%%%%%%%
		ZR(1:N,1:N)=0.0_8
		do i=1,N
			ZR(i,i)=1.0_8
		enddo
		CALL  TQL2(NM,N,W,FV1,ZR,IERR)
		CALL  HTRIBK(NM,N,AR,AI,FM1,N,ZR,ZI)

       do j = 1, N
           do i = 1, NM
               Zmat(i,j) = cmplx( ZR(i,j), ZI(i,j), 8 )
           end do
       end do

		deallocate(ZI,ZR)
		deallocate(AI,AR)
		deallocate(FM1)
		deallocate(FV1)

	end subroutine
!__________________________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ______________________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine ____________________________________________________________________________________
!__________________________________________________________________________________________________________________________________________   
	subroutine HTRIDI(NM,N,AR,AI,D,E,TAU)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
! PROGRAM:  HTRIDI(NM,N,AR,AI,D,E,TAU)
! TYPE:     subroutine
! PURPOSE:  Reduce a complex Hermitian matrix to a real symmetric tridiagonal matrix using unitary similarity transformations.
! LIBRARY:  slatec (eispack)
! CATEGORY: D4C1B1
! KEYWORDS: eigenvalues, eigenvectors, eispack
! AUTHOR:   Smith, B. T., et al.
! DESCRIPTION
!
!     This subroutine is a translation of a complex analogue of the ALGOL procedure TRED1, NUM. MATH. 11, 181-195(1968)
!        by Martin, Reinsch, and Wilkinson.HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     This subroutine reduces a COMPLEX HERMITIAN matrix to a real symmetric tridiagonal matrix using unitary similarity transformations.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional array parameters, AR and AI, as declared in the calling
!             program dimension statement.  NM is an INTEGER variable.
!        N is the order of the matrix A=(AR,AI).  N is an INTEGER variable. N must be less than or equal to NM.
!        AR and AI contain the real and imaginary parts, respectively,of the complex Hermitian input matrix.  Only the lower
!            triangle of the matrix need be supplied.  AR and AI are two-dimensional REAL arrays, dimensioned AR(NM,N) and AI(NM,N).
!
!     On OUTPUT
!
!        AR and AI contain some information about the unitary transformations used in the reduction in the strict lower triangle
!             of AR and the full lower triangle of AI.  The rest of the matrices are unaltered.
!        D contains the diagonal elements of the real symmetric tridiagonal matrix. D is a one-dimensional REAL array, dimensioned D(N).
!        E contains the subdiagonal elements of the real tridiagonal matrix in its last N-1 positions.  E(1) is set to zero. E is a 
!            one-dimensional REAL array, dimensioned E(N).
!        TAU contains further information about the transformations. TAU is a one-dimensional REAL array, dimensioned TAU(2,N).
!
!        Calls PyThagC(A,B) for sqrt(A**2+B**2).
!
!     Questions and comments should be directed to B. S. Garbow, APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
! REFERENCES:  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                   system Routines - EISPACK Guide, Springer-Verlag,1976.
! SUBROUTINES CALLED:  PyThagC
! REVISION HISTORY  (YYMMDD)
!      760101  DATE WRITTEN
!      890831  Modified array declarations.  (WRB)
!      890831  REVISION DATE from Version 3.2
!      891214  Prologue converted to Version 4.0 format.  (BAB)
!      920501  Reformatted the REFERENCES section.  (WRB)
! END PROLOGUE  HTRIDI
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
		implicit none

		integer NM,N,I,J,K,L,II,JP1
		real(8) F,G,H,FI,GI,HH,SI,SCALE
		real(8) PyThagC
		real(8) AR(NM,N),AI(NM,N)
		real(8) TAU(2,N)
		real(8) D(N),E(N)

! ***FIRST EXECUTABLE STATEMENT  HTRIDI
		TAU(1,N)=1.0_8
		TAU(2,N)=0.0_8
!
		DO 100 I=1,N
100   D(I)=AR(I,I)
!     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
		DO 300 II=1,N
			I=N+1-II
			L=I-1
			H=0.0_8
			SCALE=0.0_8
			IF (L.LT.1) GO TO 130
!     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
			DO 120 K=1,L
120      SCALE=SCALE+ABS(AR(I,K))+ABS(AI(I,K))
!
			IF (SCALE.NE.0.0_8) GO TO 140
			TAU(1,L)=1.0_8
			TAU(2,L)=0.0_8
130      E(I)=0.0_8
			GO TO 290
!
140      DO 150 K =1,L
				AR(I,K)=AR(I,K)/SCALE
				AI(I,K)=AI(I,K)/SCALE
				H=H+AR(I,K)*AR(I,K)+AI(I,K)*AI(I,K)
150      CONTINUE
!
			G=SQRT(H)
			E(I)=SCALE * G
			F=PyThagC(AR(I,L),AI(I,L))
!     .......... FORM NEXT DIAGONAL ELEMENT OF MATRIX T ..........
			IF (F.EQ.0.0_8) GO TO 160
			TAU(1,L)=(AI(I,L)*TAU(2,I)-AR(I,L)*TAU(1,I))/F
			SI=(AR(I,L)*TAU(2,I)+AI(I,L)*TAU(1,I))/F
			H=H+F*G
			G=1.0_8+G/F
			AR(I,L)=G*AR(I,L)
			AI(I,L)=G*AI(I,L)
			IF (L.EQ.1) GO TO 270
			GO TO 170
160      TAU(1,L)=-TAU(1,I)
			SI=TAU(2,I)
			AR(I,L)=G
170      F=0.0_8
!
			DO 240 J =1,L
				G=0.0_8
				GI=0.0_8
!     .......... FORM ELEMENT OF A*U ..........
				DO 180 K=1,J
					G=G+AR(J,K)*AR(I,K)+AI(J,K)*AI(I,K)
					GI=GI-AR(J,K)*AI(I,K)+AI(J,K)*AR(I,K)
180         CONTINUE
!
				JP1=J+1
				IF (L.LT.JP1) GO TO 220
!
				DO 200 K=JP1, L
					G=G+AR(K,J)*AR(I,K)-AI(K,J)*AI(I,K)
					GI=GI-AR(K,J)*AI(I,K)-AI(K,J)*AR(I,K)
200         CONTINUE
!     .......... FORM ELEMENT OF P ..........
220         E(J)=G/H
				TAU(2,J)=GI/H
				F=F+E(J)*AR(I,J)-TAU(2,J)*AI(I,J)
240      CONTINUE
!
			HH=F/(H + H)
!     .......... FORM REDUCED A ..........
			DO 260 J=1,L
				F=AR(I,J)
				G=E(J)- HH * F
				E(J)=G
				FI=-AI(I,J)
				GI=TAU(2,J)-HH*FI
				TAU(2,J)=-GI
!
				DO 260 K=1,J
					AR(J,K)=AR(J,K)-F*E(K)-G * AR(I,K)+FI*TAU(2,K)+GI*AI(I,K)
					AI(J,K)=AI(J,K)-F*TAU(2,K)-G*AI(I,K)-FI*E(K)-GI*AR(I,K)
260      CONTINUE
!
270      DO 280 K = 1, L
				AR(I,K) = SCALE * AR(I,K)
				AI(I,K) = SCALE * AI(I,K)
280      CONTINUE
!
			TAU(2,L)=-SI
290      HH=D(I)
			D(I)=AR(I,I)
			AR(I,I)=HH
			AI(I,I)=SCALE*sqrt(H)
300   CONTINUE
	  
	end subroutine
!__________________________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ______________________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine ____________________________________________________________________________________
!__________________________________________________________________________________________________________________________________________
	subroutine TQL2(NM,N,D,E,Z,IERR)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
! PROGRAM:  TQL2(NM,N,D,E,Z,IERR)
! TYPE:     subroutine
! PURPOSE:  Compute the eigenvalues and eigenvectors of symmetric tridiagonal matrix.
! LIBRARY:  SLATEC (EISPACK)
! CATEGORY: D4A5, D4C2A
! KEYWORDS: EIGENVALUES, EIGENVECTORS, EISPACK
! AUTHOR:   Smith, B. T., et al.
! DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure TQL2, NUM. MATH. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
!          Wilkinson. HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
!
!     This subroutine finds the eigenvalues and eigenvectors of a SYMMETRIC TRIDIAGONAL matrix by the QL method. The eigenvectors 
!          of a FULL SYMMETRIC matrix can also be found if  TRED2  has been used to reduce this full matrix to tridiagonal form.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional array parameter, Z, as declared in the calling program
!            dimension statement.  NM is an INTEGER variable.
!        N is the order of the matrix.  N is an INTEGER variable. N must be less than or equal to NM.
!        D contains the diagonal elements of the symmetric tridiagonal matrix.  D is a one-dimensional REAL array, dimensioned D(N).
!        E contains the subdiagonal elements of the symmetric tridiagonal matrix in its last N-1 positions.  E(1) is
!            arbitrary.  E is a one-dimensional REAL array, dimensioned E(N).
!        Z contains the transformation matrix produced in the reduction by  TRED2, if performed.  If the eigenvectors of the 
!            tridiagonal matrix are desired, Z must contain the identity matrix.  Z is a two-dimensional REAL array, dimensioned Z(NM,N).
!
!      On Output
!
!        D contains the eigenvalues in ascending order.  If an error exit is made, the eigenvalues are correct but unordered
!             for indices 1, 2, ..., IERR-1.
!        E has been destroyed.
!        Z contains orthonormal eigenvectors of the symmetric tridiagonal (or full) matrix.  If an error exit is made, Z contains 
!             the eigenvectors associated with the stored eigenvalues.
!        IERR is an INTEGER flag set to
!               Zero       for normal return,
!               J          if the J-th eigenvalue has not been determined after 30 iterations.
!
!        Calls PyThagC(A,B) for sqrt(A**2+B**2).
!
!     Questions and comments should be directed to B. S. Garbow,APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
! REFERENCES:  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-system 
!                      Routines - EISPACK Guide, Springer-Verlag,1976.
! SUBROUTINES CALLED:  PyThagC
! REVISION HISTORY  (YYMMDD)
!     760101  DATE WRITTEN
!     890831  Modified array declarations.  (WRB)
!     890831  REVISION DATE from Version 3.2
!     891214  Prologue converted to Version 4.0 format.  (BAB)
!     920501  Reformatted the REFERENCES section.  (WRB)
! END PROLOGUE  TQL2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
		implicit none

		integer NM,N,IERR
		integer I,J,K,L,M,II,L1,L2,MML
		real(8) B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
		real(8) PyThagC
		real(8) Z(NM,N)
		real(8) D(N),E(N)
		
!***FIRST EXECUTABLE STATEMENT  TQL2
		IERR=0
		IF (N .EQ. 1) GO TO 1001
!
		DO 100 I = 2, N
100   E(I-1) = E(I)
!
		F = 0.0_8
		B = 0.0_8
		E(N) = 0.0_8
!
		DO 240 L = 1, N
			J = 0
			H = ABS(D(L)) + ABS(E(L))
			IF (B .LT. H) B = H
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
			DO 110 M = L, N
				IF (B + ABS(E(M)) .EQ. B) GO TO 120
!     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
110      CONTINUE
!
120      IF (M .EQ. L) GO TO 220
130      IF (J .EQ. 30) GO TO 1000
			J = J + 1
!     .......... FORM SHIFT ..........
			L1 = L + 1
			L2 = L1 + 1
			G = D(L)
			P = (D(L1) - G) / (2.0_8 * E(L))
			R = PyThagC(P,1.0_8)
			D(L) = E(L) / (P + SIGN(R,P))
			D(L1) = E(L) * (P + SIGN(R,P))
			DL1 = D(L1)
			H = G - D(L)
			IF (L2 .GT. N) GO TO 145
!
			DO 140 I = L2, N
140      D(I) = D(I) - H
!
145      F = F + H
!     .......... QL TRANSFORMATION ..........
			P = D(M)
			C = 1.0_8
			C2 = C
			EL1 = E(L1)
			S = 0.0_8
			MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
			DO 200 II = 1, MML
				C3 = C2
				C2 = C
				S2 = S
				I = M - II
				G = C * E(I)
				H = C * P
				IF (ABS(P) .LT. ABS(E(I))) GO TO 150
				C = E(I) / P
				R = SQRT(C*C+1.0_8)
				E(I+1) = S * P * R
				S = C / R
				C = 1.0_8 / R
				GO TO 160
150         C = P / E(I)
				R = SQRT(C*C+1.0_8)
				E(I+1) = S * E(I) * R
				S = 1.0_8 / R
				C = C * S
160         P = C * D(I) - S * G
				D(I+1) = H + S * (C * G + S * D(I))
!     .......... FORM VECTOR ..........
				DO 180 K = 1, N
					H = Z(K,I+1)
					Z(K,I+1) = S * Z(K,I) + C * H
					Z(K,I) = C * Z(K,I) - S * H
180         CONTINUE
!
200      CONTINUE
!
			P = -S * S2 * C3 * EL1 * E(L) / DL1
			E(L) = S * P
			D(L) = C * P
			IF (B + ABS(E(L)) .GT. B) GO TO 130
220      D(L) = D(L) + F
240   CONTINUE
!     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
		DO 300 II = 2, N
			I = II - 1
			K = I
			P = D(I)
!
			DO 260 J = II, N
				IF (D(J) .GE. P) GO TO 260
				K = J
				P = D(J)
260      CONTINUE
!		
			IF (K .EQ. I) GO TO 300
			D(K) = D(I)
			D(I) = P
!
			DO 280 J = 1, N
				P = Z(J,I)
				Z(J,I) = Z(J,K)
				Z(J,K) = P
280      CONTINUE
!
300   CONTINUE
!
		GO TO 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
1000  IERR=L
1001  RETURN

	end subroutine
!__________________________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ______________________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin subroutine ____________________________________________________________________________________
!__________________________________________________________________________________________________________________________________________
	subroutine HTRIBK(NM,N,AR,AI,TAU,M,ZR,ZI)      
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
! PROGRAM:  HTRIBK(NM,N,AR,AI,TAU,M,ZR,ZI)
! TYPE:     subroutine
! PURPOSE:  Form the eigenvectors of a complex Hermitian matrix from the eigenvectors of a real symmetric tridiagonal matrix
!              output from HTRIDI.
! LIBRARY:  SLATEC (EISPACK)
! CATEGORY: D4C4
! KEYWORDS: EIGENVALUES, EIGENVECTORS, EISPACK
! AUTHOR:   Smith, B. T., et al.
! DESCRIPTION
!
!     This subroutine is a translation of a complex analogue of the ALGOL procedure TRBAK1, NUM. MATH. 11, 181-195(1968)
!          by Martin, Reinsch, and Wilkinson. HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     This subroutine forms the eigenvectors of a COMPLEX HERMITIAN matrix by back transforming those of the corresponding
!          real symmetric tridiagonal matrix determined by  HTRIDI.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional array parameters, AR, AI, ZR, and ZI, as declared in the
!            calling program dimension statement.  NM is an INTEGER variable.
!        N is the order of the matrix.  N is an INTEGER variable. N must be less than or equal to NM.
!        AR and AI contain some information about the unitary transformations used in the reduction by  HTRIDI  in the strict 
!            lower triangle of AR and the full lower triangle of AI.  The remaining upper parts of the matrices are arbitrary.
!            AR and AI are two-dimensional REAL arrays, dimensioned AR(NM,N) and AI(NM,N).
!        TAU contains further information about the transformations. TAU is a one-dimensional REAL array, dimensioned TAU(2,N).
!        M is the number of eigenvectors to be back transformed. M is an INTEGER variable.
!        ZR contains the eigenvectors to be back transformed in its first M columns.  The contents of ZI are immaterial.  
!            ZR and ZI are two-dimensional REAL arrays, dimensioned ZR(NM,M) and ZI(NM,M).
!
!     On OUTPUT
!
!        ZR and ZI contain the real and imaginary parts, respectively, of the transformed eigenvectors in their first M columns.
!
!     Note that the last component of each returned vector is real and that vector Euclidean norms are preserved.
!
!     Questions and comments should be directed to B. S. Garbow, APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
! REFERENCES:  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-system
!                      Routines - EISPACK Guide, Springer-Verlag,1976.
! SUBROUTINES CALLED:  (NONE)
! REVISION HISTORY  (YYMMDD)
!      760101  DATE WRITTEN
!      890831  Modified array declarations.  (WRB)
!      890831  REVISION DATE from Version 3.2
!      891214  Prologue converted to Version 4.0 format.  (BAB)
!      920501  Reformatted the REFERENCES section.  (WRB)
! END PROLOGUE  HTRIBK
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
		implicit none

		integer NM,N,I,J,K,L,M
	   real(8) H,S,SI
	   real(8) AR(NM,N),AI(NM,N),ZR(NM,N),ZI(NM,N)
	   real(8) TAU(2,N)

!***FIRST EXECUTABLE STATEMENT  HTRIBK
		IF (M .EQ. 0) GO TO 200
!     .......... TRANSFORM THE EIGENVECTORS OF THE REAL SYMMETRIC TRIDIAGONAL MATRIX TO THOSE OF THE HERMITIAN TRIDIAGONAL MATRIX. ..........
		DO 50 K = 1, N
!
			DO 50 J = 1, M
				ZI(K,J) = -ZR(K,J) * TAU(2,K)
				ZR(K,J) = ZR(K,J) * TAU(1,K)
50    CONTINUE
!
		IF (N .EQ. 1) GO TO 200
!     .......... RECOVER AND APPLY THE HOUSEHOLDER MATRICES ..........
		DO 140 I = 2, N
			L = I - 1
			H = AI(I,I)
			IF (H .EQ. 0.0_8) GO TO 140
!
			DO 130 J = 1, M
				S = 0.0_8
				SI = 0.0_8
!
				DO 110 K = 1, L
					S = S + AR(I,K) * ZR(K,J) - AI(I,K) * ZI(K,J)
					SI = SI + AR(I,K) * ZI(K,J) + AI(I,K) * ZR(K,J)
110         CONTINUE
!     .......... DOUBLE DIVISIONS AVOID POSSIBLE UNDERFLOW ..........
				S = (S / H) / H
				SI = (SI / H) / H
!
				DO 120 K = 1, L
					ZR(K,J) = ZR(K,J) - S * AR(I,K) - SI * AI(I,K)
					ZI(K,J) = ZI(K,J) - SI * AR(I,K) + S * AI(I,K)
120         CONTINUE
130      CONTINUE
140   CONTINUE
200   RETURN

	end subroutine
!__________________________________________________________________________________________________________________________________________  
!____________________________________ End subroutine ______________________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!____________________________________ Begin Function ______________________________________________________________________________________
!__________________________________________________________________________________________________________________________________________
	function PyThagC(a, b)
!	find sqrt(a** 2+b**2) without overflow or destructive underflow

		!use RealPrecsn
		implicit none

		real(8) PyThagC, a, b, p, r, s, t

		p = max(abs(a), abs(b))
	
		if(p == 0.0_8) then
			PyThagC = p
			return
		end if

		r = (min(abs(a), abs(b)) / p)**2
		do while (4 + r /= 4)
			s = r / (4.0_8 + r)
			t = 1 + 2 * s
			p = t * p
			r = (s / t)**2 * r
		enddo

		PyThagC = p

	end function
!__________________________________________________________________________________________________________________________________________  
!____________________________________ End Funtion _________________________________________________________________________________________
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
