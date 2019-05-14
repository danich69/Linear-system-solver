module Linear_system
	implicit none
	contains
	
	subroutine Lin_sys_solver(a,b,x,s_size,method_type)
	
		integer :: s_size, i, k, num_bufer,j
		integer, dimension(s_size) ::  order							! order - vector of actual x indexses, is used for column permutations
		integer, dimension(2) :: max_coords								! coords of maxloc in current matrix
		real(8), dimension(s_size+1) :: row_bufer					
		real(8), dimension(s_size, s_size) :: a
		real(8), dimension(s_size, s_size+1) :: new
		real(8), dimension(s_size) :: column_bufer, b, x
		character*14 method_type
				
		new=a
		new(:,s_size+1)=b
		
		do i=1, s_size
			order(i)=i
		enddo
			
		if(method_type=='ModGauss') then								!changes for 
				
			do k=1, s_size  											!determine the maximum element and using row and column permutations put it on the diagonal
				max_coords(:)=MAXLOC(abs(new(k:s_size,k:s_size)))+(k-1)
				if(k /= max_coords(1)) then								!row permutation 
					row_bufer=new(k,:)
					new(k,:)=new(max_coords(1),:)
					new(max_coords(1),:)=row_bufer
				endif
			
				if(k /= max_coords(2)) then								!column permutation
					column_bufer=new(:,k)
					new(:,k)=new(:,max_coords(2))
					new(:,max_coords(2))=column_bufer
				
					num_bufer=order(k)									!remember which columns were permutated, as it's a permutation of solution vector coordinates
					order(k)=order(max_coords(2))							
					order(max_coords(2))=num_bufer
				endif
			enddo			
		
		endif
		
            write(*,*) new
            write(*,*) 
		
		do k=1, s_size
!            if(method_type=='ModGauss') then
!                max_coords(:)=MAXLOC(abs(new(k:s_size,k:s_size)))+(k-1)
!				if(k /= max_coords(1)) then								!row permutation 
!					row_bufer=new(k,:)
!					new(k,:)=new(max_coords(1),:)
!					new(max_coords(1),:)=row_bufer
!				endif
			
!				if(k /= max_coords(2)) then								!column permutation
!					column_bufer=new(:,k)
!					new(:,k)=new(:,max_coords(2))
!					new(:,max_coords(2))=column_bufer
				
!					num_bufer=order(k)									!remember which columns were permutated, as it's a permutation of solution vector coordinates
!					order(k)=order(max_coords(2))							
!					order(max_coords(2))=num_bufer
!				endif
!            write(*,*) new
!            write(*,*) 
!            endif
            
			do i=k+1,s_size
			
				if( abs(new(k, k)) < 0.0001 ) then
					write(*,*) 'Division by a number close to 0. Row number ',k  !check weather program devide by "0" or not
				endif
            write(*,*) new
            write(*,*) 
				
				new(i,:)=new(i,:)-new(k,:)*new(i,k)/new(k,k)			! 1-st step: subtract the first row multiplied by the coefficient 
				
				if(abs(new(i, i)) < 0.0001) then
					write(*,*) 'Division by a number close to 0. Row number ',i  !check weather program devide by "0" or not
					
				endif
				new(i,:)=new(i,:)/new(i,i)								! normalize new row
			enddo
		enddo
		
		
		if(method_type=='Jordan') then
			do i=1, s_size
				do k=1, i-1
				new(i,k)=0												!zeroing of elements which should be 0, but in fact can be very small (exmpl - 1.75834E-8)
				enddo
			enddo
		
			do k=s_size-1, 1, -1
				do i=k+1, s_size
					new(k,:)=new(k,:)-new(i,:)*new(k,i)/new(i,i)		!bring the triangular matrix to a diagonal view
				enddo
			enddo
				
			if(new(1,1)/=1) then
				new(1,:)=new(1,:)/new(1,1)								! normalize first row		
			endif
			
			x=new(:,s_size+1)
		
		else
		
			if(new(1,1)/=1) then
				new(1,:)=new(1,:)/new(1,1)								! normalize first row		
			endif		
		
			x(order(s_size))=new(s_size,s_size+1)						! calculation of x(n)
		
			do i=s_size-1, 1, -1
				x(order(i))=new(i,s_size+1)
				
				do j=s_size,i+1, -1
					x(order(i))=x(order(i))-new(i,j)*x(order(j))		! 2-nd step: calculation of resulting vector remembering column permutations
				enddo
			
			enddo
		endif
		
	end subroutine Lin_sys_solver
	
end module Linear_system
