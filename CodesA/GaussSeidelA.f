      program GS !programme principal
        implicit none
        interface !debut interface
          subroutine matriceCarre(A,carre)
            implicit none
            real,dimension(:,:), intent(in) :: A
            logical,intent(out) :: carre
          endsubroutine matriceCarre

          subroutine diagonaleDom(A,diag)
            implicit none
            real,dimension(:,:), intent(in) :: A
            logical,intent(out) :: diag
          endsubroutine diagonaleDom

          subroutine GaussSeidelFunction(A,B,Eps,x,m,n)
            implicit none
            integer,intent(in) :: n
            real, intent(in) :: Eps
            real, dimension(:,:), intent(in) :: A
            real, dimension(:),intent(in) :: B
            real,dimension(:), intent(inout) :: x
            integer, intent(out) :: m
          endsubroutine GaussSeidelFunction
          
          subroutine MatHilbert(A,w,n)
            implicit none
            real,dimension(:,:),intent(out) :: A
            integer,intent(in) :: n
            real,intent(in) :: w
            integer :: i,j
          endsubroutine MatHilbert
          
          subroutine CalculB(A,B,n)
            real,dimension(:,:),intent(in) :: A
            real,dimension(:), intent(out) :: B
            integer, intent(in) :: n
            integer :: i
          endsubroutine CalculB
          
        endinterface!fin interface
        real, dimension(:,:),allocatable :: A
        real,dimension(:),allocatable :: B
        real, dimension(:),allocatable :: x
        real :: Eps,w
        integer :: m,n,dim
        logical :: carre,diag
	 !lit les paramètres
        print*, "Ecrivez la taille de la matrice que vous souhaitez"
        read(*,*) dim
        

        !allocation
        allocate(A(dim,dim))
        allocate(B(dim))
        allocate(x(dim))
        
        print*,"Ecrivez le w que vous souhaitez"
        read(*,*) w
        

        call MatHilbert(A,w,dim)
        call CalculB(A,B,dim)
        
        print*,"Ecrivez votre epsilon"
        read(*,*) Eps
        
        print*,"Ecrivez le nombre d'itérations que vous souhaitez"
        read(*,*) n
        


        call matriceCarre(A,carre)
        call diagonaleDom(A,diag)
     	 
	 !si à diag dominante et matrice carre
        if (diag .and. carre) then
          call GaussSeidelFunction(A,B,Eps,x,m,n)!appel de la subroutine Gauss Seidel
          print*,m,x
        else
          print*,carre,diag
        endif

        deallocate(x)
        deallocate(A)
        deallocate(B)



      endprogram GS!fin programme principale

      subroutine GaussSeidelFunction(A,B,Eps,x0,m,n)!début subroutine GaussSeidel
        implicit none
        interface
          !appel des subroutines dans l'interface
          subroutine inverse_matrice(A, A_inv, n)
            implicit none
            real,dimension(:,:), intent(in) :: A
            real, dimension(:,:),intent(out) :: A_inv
            real :: augmented(2*n, n)
            integer :: i, j, k
            real :: factor
            integer, intent(in) ::n
          endsubroutine inverse_matrice
          
          subroutine diagonale(A,D)
      	     implicit none
      	     real,dimension(:,:),intent(in) :: A
      	     real,dimension(:,:),intent(out) :: D
      	     integer,dimension(:),allocatable :: dims1
      	     integer :: i,dims
      	   endsubroutine diagonale
      	   
      	   subroutine triangulaireInf(A,L)
     	     implicit none
     	     real,dimension(:,:),intent(in) :: A
     	     real,dimension(:,:),intent(out) :: L
     	     integer :: i,j,dims
     	   endsubroutine triangulaireInf
     	   
     	   subroutine triangulaireSupp(A,U)
            implicit none
            real,dimension(:,:),intent(in) :: A
            real,dimension(:,:),intent(out) :: U
            integer :: i,j,dims
          endsubroutine triangulaireSupp
          
          subroutine matrice_multiplication(A, B, C, n)
            real,dimension(:,:), intent(in) :: A, B
            real,dimension(:,:), intent(out) :: C
            integer :: i, j, k
            integer, intent(in) ::n
          endsubroutine matrice_multiplication
          
          subroutine matrice_vect_multiplication(A, b, c, n)
            real,dimension(:,:), intent(in) :: A
            real,dimension(:), intent(in) :: b
            real,dimension(:), intent(out) :: c
            integer :: i, j
            integer, intent(in) ::n
          endsubroutine matrice_vect_multiplication
        endinterface
     	   
      	   
        integer, intent(in) :: n
        real,intent(in) :: Eps
        real,dimension(:,:), intent(in) :: A
        real,dimension(:),intent(in) :: B
        real,dimension(:),intent(inout) :: x0
        integer, intent(out) :: m
        real :: erreur
        integer ::i,j,iter,matdim
        real,dimension(:,:),allocatable :: D, L, U, F, DL, DL_inv
        real,dimension(:),allocatable :: G,x1

        
        
        !initialisation
        matdim=size(A,1)
        x0=0.0
        m=0
        erreur=Eps+1.0
        
        
                
        allocate(x1(matdim))
        allocate(F(matdim,matdim))
        allocate(G(matdim))
        allocate(D(matdim,matdim))
        allocate(U(matdim,matdim))
        allocate(L(matdim,matdim))
        allocate(DL(matdim,matdim))
        allocate(DL_inv(matdim,matdim))
        
        
        call diagonale(A,D) !calcul diagonale
        call triangulaireInf(A,L) ! Calcul triangulaire Inf
        call triangulaireSupp(A,U)! Calcul triangulaire supp
        
        
        
        DL = D + L !D+L
    
        call inverse_matrice(DL, DL_inv, matdim)!calcul de DL_inv
        call matrice_multiplication(DL_inv, U, F, matdim)!calcul de F=DL_inv*U
        F = -F !F=-F
        
        call matrice_vect_multiplication(DL_inv, B, G, matdim)!Calcul G=DL_inv*B
	 !itération X_k+1=FX_k+G
	 x1=x0
        do iter=1,n
          if (erreur>Eps) then      
            m=m+1
            call matrice_vect_multiplication(F, x0, x1, matdim)
            x1=x1+G
            erreur=maxval(abs(x1-x0))
            print*,"erreur :"
            print*,erreur
          endif
          x0=x1
        enddo
        
        
        
        
     	
        
        if (erreur>Eps) then
          print *, "Pas de convergence"
        endif
	 
	 
	 
        deallocate(x1)
        deallocate(F)
        deallocate(G)
 	 deallocate(D)
 	 deallocate(U)
 	 deallocate(L)
 	 deallocate(DL)
 	 deallocate(DL_inv)      
      endsubroutine GaussSeidelFunction


      subroutine matriceCarre(A,carre)
        implicit none
        real,dimension(:,:),intent(in) :: A
        integer, dimension(:), allocatable : : dims
        logical,intent(out) :: carre

        dims=shape(A)
        carre=.true.!initialisation matrice carre vraie
        
        if (dims(1) /= dims(2)) then
                carre=.false.!si pas vérifiée alors faux
        endif
      endsubroutine matriceCarre


      subroutine diagonaleDom(A,diag)
        implicit none
        real,dimension(:,:),intent(in) :: A
        logical,intent(out) :: diag
        integer, dimension(:), allocatable :: dims1
        integer :: i,j,dims 
        real :: somme
        
        dims1=shape(A)
        dims= dims1(1)
        diag=.true.!initialisation diagonale Dominante vraie

        do i=1,dims
          somme=0.0
          do j=1,dims
            if (i/=j) then
              somme=somme+abs(A(i,j))
            endif
          enddo
          if (abs(A(i,i)) < somme) then!si pas vérifiée alors faux
            diag=.false.
          endif
        enddo
      endsubroutine diagonaleDom
      
      
      subroutine diagonale(A,D)
      	implicit none
      	real,dimension(:,:),intent(in) :: A
      	real,dimension(:,:),intent(out) :: D
      	integer,dimension(:),allocatable :: dims1
      	integer :: i,dims
      	
      	dims1=shape(A)
      	dims=dims1(1)
      	
      	do i=1,dims
      	  D(i,i)=A(i,i)!calcul diagonale
      	enddo
      	
      endsubroutine diagonale
     
     
      subroutine triangulaireInf(A,L)
     	 implicit none
     	 real,dimension(:,:),intent(in) :: A
     	 real,dimension(:,:),intent(out):: L
     	 integer :: i,j,dims
     	 
     	 dims = size(A, 1)!Dimension de A
     	 L=0.0
     	 
     	 do i=1,dims
     	   do j=1,dims
     	     if (i>j) then
     	       L(i,j)=A(i,j)!Calcul triangulaire Inférieure
     	       
     	     endif
     	   enddo
        enddo
     	
       
     	
      endsubroutine triangulaireInf
      
      subroutine triangulaireSupp(A,U)
        implicit none
        real,dimension(:,:),intent(in) :: A
        real,dimension(:,:),intent(out):: U
        integer :: i,j,dims
	 
	 
	 dims = size(A, 1)!Dimension de A
	
        U=0.0
        
        do i=1,dims
     	   do j=1,dims
     	     if (i<j) then
     	       U(i,j)=A(i,j) !Calcul triangulaire suppérieure
     	     endif
     	   enddo
        enddo
       
      endsubroutine triangulaireSupp
      
      subroutine inverse_matrice(A, A_inv, n)
      	 implicit none
        real,dimension(:,:), intent(in) :: A
        integer, intent(in) ::n
        real, dimension(:,:),intent(out) :: A_inv
        real :: augmented(2*n, n)
        integer :: i, j, k
        real :: factor

          ! Construire la matrice augmentée [A | I]
        augmented(1:n, 1:n) = A
        augmented(1:n, n+1:2*n) = 0.0
        do i = 1, n
          augmented(i, n+i) = 1.0
        end do

        ! Gauss-Jordan : mise en forme échelonnée
        do i = 1, n
          ! Normalisation de la ligne pivot
          factor = augmented(i, i)
          do j = 1, 2*n
            augmented(i,j) = augmented(i,j) / factor
          end do

            ! Élimination des autres lignes
          do k = 1, n
            if (k /= i) then
              factor = augmented(k, i)
              do j = 1, 2*n
                augmented(k,j)=augmented(k,j) -factor*augmented(i, j)
              end do
            end if
          end do
        end do

          ! Extraction de l'inverse
        A_inv = augmented(1:n, n+1:2*n)
      end subroutine inverse_matrice

      ! Sous-programme pour le produit de deux matrices
      subroutine matrice_multiplication(A, B, C, n)
        real,dimension(:,:), intent(in) :: A, B
        integer, intent(in) :: n
        real,dimension(:,:), intent(out) :: C
        integer :: i, j, k

        C = 0.0
        do i = 1, n
          do j = 1, n
            do k = 1, n
              C(i, j) = C(i, j) + A(i, k) * B(k, j)
            end do
          end do
        end do
      end subroutine matrice_multiplication

      ! Sous-programme pour le produit matrice-vecteur
      subroutine matrice_vect_multiplication(A, b, c, n)
        real,dimension(:,:), intent(in) :: A
        real,dimension(:), intent(in) :: b
        real,dimension(:), intent(out) :: c
        integer, intent(in) ::n
        integer :: i, j

        c = 0.0
        do i = 1, n
          do j = 1, n
            c(i) = c(i) + A(i, j) * b(j)
          end do
        end do
      end subroutine matrice_vect_multiplication
      !calcul de la matrice de Hilbert
      subroutine MatHilbert(A,w,n)
        implicit none
        real,dimension(:,:),intent(out) :: A
        integer,intent(in) :: n
        real,intent(in) :: w
        integer :: i,j
        
        do i=1,n
          do j=1,n
            A(i,j)=1.0/real(i+j-1)
          enddo
        enddo
    	
    	do i=1,n
    	  A(i,i)=A(i,i)+w
    	enddo
      endsubroutine MatHilbert
      !Caclul du vecteur B
      subroutine CalculB(A,B,n)
        real,dimension(:,:),intent(in) :: A
        real,dimension(:), intent(out) :: B
        integer, intent(in) :: n
        integer :: i
        
        do i=1,n
          B(i)=sum(A(i,:))
        enddo
      endsubroutine CalculB
      
      
      
      
           	











                
