program in_out
	use Linear_system
	implicit none
	
!	----------------------------------------------------
!	|	Program solves linear equation system		   |
!	|												   |
!	|	Method types and keys:						   |
!	|												   |
!	|		Simple Gauss method -- Gauss			   |
!	|		Jordan method -- Jordan 				   |
!	|		Modified Gauss method -- ModGauss	   |
!	|												   |
!	----------------------------------------------------

	real(8), dimension(:,:), allocatable :: A							!system amtrix
	real(8), dimension(:), allocatable :: B, X							!column vector and solution vector
	integer Sys_size 													!square matrix size
	integer i,j,k														!string
	character*14 method_type											!Method key
	
	
	call getarg(1, method_type)											!Console method key input
	
	do while(method_type/='Gauss' .and. method_type/='Jordan' .and. method_type/='ModGauss')	!Corectness check
		write(*,*) 'Wrong input! Try Again!', method_type
		read(*,*) method_type
	enddo
	
	
	open ( unit=1, file='data.dat', status='old', action='read' )
	open ( unit=3, file='result.dat',status='REPLACE')

	read(1,'(2x,I6)') Sys_size											!read sharp symbol and matrix size
	
	allocate(A(Sys_size,Sys_size))
	allocate(B(Sys_size),X(Sys_size))
		
	read(1,*) A															!read system matrix
	A=TRANSPOSE(A)
	
	do i=1, Sys_size
		read(1,*) B(i)													!read column vector
	enddo
	
	call Lin_sys_solver(A,B,X,Sys_size,method_type)								!call subroutine to choose the method 
	
	write(3,55) Sys_size												!print system size
	55 format('#',I5)													
	
	do j=1,Sys_size
				write(3,*) X(j)											!write solution vector
	enddo
	
	write(*,*) 'Norm of Residual vector is:', norm2( matmul(A,X)-B )
	
	close(1)
	close(3)
	
	deallocate(A,B,X)
	
end program in_out
