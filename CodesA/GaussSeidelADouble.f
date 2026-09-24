      program GS!début programme principal
        implicit none
        integer, parameter :: dp = kind(0.d0)!Paramètre pour mettre en double précision
        interface!Début interface pour appeler subroutine
          subroutine matriceCarre(A,carre)
            implicit none
            integer, parameter :: dp = kind(0.d0)
            real(dp),dimension(:,:), intent(in) :: A(:,:)
            logical,intent(out) :: carre
          endsubroutine matriceCarre

          subroutine diagonaleDom(A,diag)
            implicit none
            integer, parameter :: dp = kind(0.d0)
            real(dp),dimension(:,:), intent(in) :: A
            logical,intent(out) :: diag
          endsubroutine diagonaleDom

          subroutine GaussSeidelFunction(A,B,Eps,x,m,n)
            implicit none
            integer, parameter :: dp = kind(0.d0)
            integer,intent(in) :: n
            real(dp), intent(in) :: Eps
            real(dp), dimension(:,:), intent(in) :: A
            real(dp), dimension(:),intent(in) :: B
            real(dp),dimension(:), intent(inout) :: x
            integer, intent(out) :: m
          endsubroutine GaussSeidelFunction
          
          subroutine MatHilbert(A,w,n)
            implicit none
            integer, parameter :: dp = kind(0.d0)
            real(dp),dimension(:,:),intent(out) :: A
            integer,intent(in) :: n
            real(dp),intent(in) :: w
            integer :: i,j
          endsubroutine MatHilbert
          
          subroutine CalculB(A,B,n)
            implicit none
            integer, parameter :: dp = kind(0.d0)
            real(dp),dimension(:,:),intent(in) :: A
            real(dp),dimension(:), intent(out) :: B
            integer, intent(in) :: n
            integer :: i
          endsubroutine CalculB
          
        endinterface!fin interface
        real(dp), dimension(:,:),allocatable :: A
        real(dp),dimension(:),allocatable :: B
        real(dp), dimension(:),allocatable :: x
        real(dp) :: Eps,w
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
     	 
	 !si à diagonale dominate et matrice carré
        if (diag .and. carre) then
          call GaussSeidelFunction(A,B,Eps,x,m,n)!appel de la subroutine Gauss Siedel
          print*,m,x
        else
          print*,carre,diag
        endif

        deallocate(x)
        deallocate(A)
        deallocate(B)



      endprogram GS!fin programme principal

      subroutine GaussSeidelFunction(A,B,Eps,x0,m,n)
        implicit none
        integer, parameter :: dp = kind(0.d0)
        interface
          !appel des subroutines dans l'interface
          subroutine inverse_matrice(A, A_inv, n)
            implicit none
            integer, parameter :: dp = kind(0.d0)
            real(dp),dimension(:,:), intent(in) :: A
            real(dp), dimension(:,:),intent(out) :: A_inv
            real(dp) :: augmented(2*n, n)
            integer :: i, j, k
            real(dp) :: factor
            integer, intent(in) ::n
          endsubroutine inverse_matrice
          
          subroutine diagonale(A,D)
      	     implicit none
      	     integer, parameter :: dp = kind(0.d0)
      	     real(dp),dimension(:,:),intent(in) :: A
      	     real(dp),dimension(:,:),intent(out) :: D
      	     integer,dimension(:),allocatable :: dims1
      	     integer :: i,dims
      	   endsubroutine diagonale
      	   
      	   subroutine triangulaireInf(A,L)
     	     implicit none
     	     integer, parameter :: dp = kind(0.d0)
     	     real(dp),dimension(:,:),intent(in) :: A
     	     real(dp),dimension(:,:),intent(out) :: L
     	     integer :: i,j,dims
     	   endsubroutine triangulaireInf
     	   
     	   subroutine triangulaireSupp(A,U)
            implicit none
            integer, parameter :: dp = kind(0.d0)
            real(dp),dimension(:,:),intent(in) :: A
            real(dp),dimension(:,:),intent(out) :: U
            integer :: i,j,dims
          endsubroutine triangulaireSupp
          
          subroutine matrice_multiplication(A, B, C, n)
            implicit none
            integer, parameter :: dp = kind(0.d0)
            real(dp),dimension(:,:), intent(in) :: A, B
            real(dp),dimension(:,:), intent(out) :: C
            integer :: i, j, k
            integer, intent(in) ::n
          endsubroutine matrice_multiplication
          
          subroutine matrice_vect_multiplication(A, b, c, n)
            implicit none
            integer, parameter :: dp = kind(0.d0)
            real(dp),dimension(:,:), intent(in) :: A
            real(dp),dimension(:), intent(in) :: b
            real(dp),dimension(:), intent(out) :: c
            integer :: i, j
            integer, intent(in) ::n
          endsubroutine matrice_vect_multiplication
        endinterface
     	   
      	   
        integer, intent(in) :: n
        real(dp),intent(in) :: Eps
        real(dp),dimension(:,:), intent(in) :: A
        real(dp),dimension(:),intent(in) :: B
        real(dp),dimension(:),intent(inout) :: x0
        integer, intent(out) :: m
        real(dp) :: erreur
        integer ::i,j,iter,matdim
        real(dp),dimension(:,:),allocatable :: D, L, U, F, DL, DL_inv
        real(dp),dimension(:),allocatable :: G,x1

        
        
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
        
        
        call diagonale(A,D)!calcul diagonale
        call triangulaireInf(A,L)!calul de matrice Inférieure
        call triangulaireSupp(A,U)!calcul de matrice suppérieure
        
        
        
        DL = D + L!D+L
    
        call inverse_matrice(DL, DL_inv, matdim)!calcul de DL_inv
        call matrice_multiplication(DL_inv, U, F, matdim)!Calcul de F=DL_inv*U
        F = -F!F=-F
        
        call matrice_vect_multiplication(DL_inv, B, G, matdim)!calcul de G=DL_inv*B
	 !iteration X_k+1=FX_k+G
	 x1=x0
        do iter=1,n
          if (erreur>Eps) then      
            m=m+1
            call matrice_vect_multiplication(F, x0, x1, matdim)
            x1=x1+G
            erreur=maxval(abs(x1-x0))
            print*,"Erreur:"
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
        integer, parameter :: dp = kind(0.d0)
        real(dp),dimension(:,:),intent(in) :: A
        integer, dimension(:), allocatable : : dims
        logical,intent(out) :: carre

        dims=shape(A)
        carre=.true.!initialisation matrice carre vraie
        
        if (dims(1) /= dims(2)) then!si pas carre alors faux
                carre=.false.
        endif
      endsubroutine matriceCarre


      subroutine diagonaleDom(A,diag)
        implicit none
        integer, parameter :: dp = kind(0.d0)
        real(dp),dimension(:,:),intent(in) :: A
        logical,intent(out) :: diag
        integer, dimension(:), allocatable :: dims1
        integer :: i,j,dims 
        real(dp) :: somme
        
        dims1=shape(A)
        dims= dims1(1)
        diag=.true.!initialisation diagonale dominante vraie

        do i=1,dims
          somme=0.0
          do j=1,dims
            if (i/=j) then!si pas respéctée alors faux
              somme=somme+abs(A(i,j))
            endif
          enddo
          if (abs(A(i,i)) < somme) then
            diag=.false.
          endif
        enddo
      endsubroutine diagonaleDom
      
      
      subroutine diagonale(A,D)
      	implicit none
      	integer, parameter :: dp = kind(0.d0)
      	real(dp),dimension(:,:),intent(in) :: A
      	real(dp),dimension(:,:),intent(out) :: D
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
     	 integer, parameter :: dp = kind(0.d0)
     	 real(dp),dimension(:,:),intent(in) :: A
     	 real(dp),dimension(:,:),intent(out):: L
     	 integer :: i,j,dims
     	 
     	 dims = size(A, 1)
     	 L=0.0
     	 
     	 do i=1,dims
     	   do j=1,dims
     	     if (i>j) then!calcul triangulaire inférieure
     	       L(i,j)=A(i,j)
     	       
     	     endif
     	   enddo
        enddo
     	
       
     	
      endsubroutine triangulaireInf
    
      subroutine triangulaireSupp(A,U)
        implicit none
        integer, parameter :: dp = kind(0.d0)
        real(dp),dimension(:,:),intent(in) :: A
        real(dp),dimension(:,:),intent(out):: U
        integer :: i,j,dims
	 
	 
	 dims = size(A, 1)
	
        U=0.0
        
        do i=1,dims
     	   do j=1,dims
     	     if (i<j) then!calcul triangulaire supérieure
     	       U(i,j)=A(i,j)
     	     endif
     	   enddo
        enddo
       
      endsubroutine triangulaireSupp
      
      subroutine inverse_matrice(A, A_inv, n)
      	 implicit none
      	 integer, parameter :: dp = kind(0.d0)
        real(dp),dimension(:,:), intent(in) :: A
        integer, intent(in) ::n
        real(dp), dimension(:,:),intent(out) :: A_inv
        real(dp) :: augmented(2*n, n)
        integer :: i, j, k
        real(dp) :: factor

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
      	 implicit none
      	 integer, parameter :: dp = kind(0.d0)
        real(dp),dimension(:,:), intent(in) :: A, B
        integer, intent(in) :: n
        real(dp),dimension(:,:), intent(out) :: C
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
      	 implicit none
      	 integer, parameter :: dp = kind(0.d0)
        real(dp),dimension(:,:), intent(in) :: A
        real(dp),dimension(:), intent(in) :: b
        real(dp),dimension(:), intent(out) :: c
        integer, intent(in) ::n
        integer :: i, j

        c = 0.0
        do i = 1, n
          do j = 1, n
            c(i) = c(i) + A(i, j) * b(j)
          end do
        end do
      end subroutine matrice_vect_multiplication
      
      subroutine MatHilbert(A,w,n)!calcul matrice de Hilbert
        implicit none
        integer, parameter :: dp = kind(0.d0)
        real(dp),dimension(:,:),intent(out) :: A
        integer,intent(in) :: n
        real(dp),intent(in) :: w
        integer :: i,j
        
        do i=1,n
          do j=1,n
            A(i,j)=1.0/dble(i+j-1)
          enddo
        enddo
    	
    	do i=1,n
    	  A(i,i)=A(i,i)+w
    	enddo
      endsubroutine MatHilbert
      
      subroutine CalculB(A,B,n)!calcul du vecteur B
      	 implicit none
      	 integer, parameter :: dp = kind(0.d0)
        real(dp),dimension(:,:),intent(in) :: A
        real(dp),dimension(:), intent(out) :: B
        integer, intent(in) :: n
        integer :: i
        
        do i=1,n
          B(i)=sum(A(i,:))
        enddo
      endsubroutine CalculB
      
      
      
      
           	











                
