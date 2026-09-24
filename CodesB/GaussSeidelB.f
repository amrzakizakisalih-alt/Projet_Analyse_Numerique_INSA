      program GS!début programme principal
        implicit none
        interface!debut interface pour appeler subroutine
          subroutine matriceCarre(A,carre)
            implicit none
            real,dimension(:,:), intent(in) :: A
            logical,intent(out) :: carre
          endsubroutine matriceCarre

          subroutine diagonale(A,diag)
            implicit none
            real,dimension(:,:), intent(in) :: A
            logical,intent(out) :: diag
          endsubroutine diagonale

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
	 ! lit les paramètres
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
        call diagonale(A,diag)
     
	 !si matrice à diagonale dominante et carre alors on appelle la subroutine GaussSeidel
        if (diag .and. carre) then
          call GaussSeidelFunction(A,B,Eps,x,m,n)
          print*,m,x
        else
          print*,carre,diag
        endif
        
        deallocate(A)
        deallocate(B)
        deallocate(x)

        



      endprogram GS

      subroutine GaussSeidelFunction(A,B,Eps,x,m,n)!début de la subroutine GaussSeidel
        implicit none
        integer, intent(in) :: n
        real,intent(in) :: Eps
        real,dimension(:,:), intent(in) :: A
        real,dimension(:),intent(in) :: B
        real,dimension(:),intent(inout) :: x
        integer, intent(out) :: m
        real, dimension(:), allocatable :: y
        real :: erreur
        integer :: tailleB,i,j,iter

        tailleB=size(B)
        x=0.0
        m=0
        erreur=Eps+1.0
                
        allocate(y(tailleB))
	 !itération sans décomposition de la matrice
        do iter=1,n
          
          if (erreur>Eps) then      
            m=m+1

            do i=1,tailleB   
              y(i)=B(i)                         
              do j=1,i-1
                y(i)=y(i)-A(i,j)*x(j)
              enddo                 
              do j=i+1,tailleB
                y(i)=y(i)-A(i,j)*x(j)
              enddo
                                        
              y(i)=y(i)/A(i,i)
            enddo
            erreur=maxval(abs(y-x))
            x=y
          endif
        enddo
        
        if (erreur>Eps) then
          print *, "Pas de convergence"
        endif

        deallocate(y)
        
      endsubroutine GaussSeidelFunction
      !subroutine matriceCarre
      subroutine matriceCarre(A,carre)
        implicit none
        real,dimension(:,:),intent(in) :: A
        integer, dimension(:), allocatable : : dims
        logical,intent(out) :: carre

        dims=shape(A)
        
        carre=.true.!carre initialisé à vrai
        
        

        if (dims(1) /= dims(2)) then!si matrice non carré alors faux
                carre=.false.
        endif
      endsubroutine matriceCarre
     
      subroutine diagonale(A,diag)
        implicit none
        real,dimension(:,:),intent(in) :: A
        logical,intent(out) :: diag
        integer, dimension(:), allocatable :: dims1
        integer :: i,j,dims 
        real :: somme
        
        dims1=shape(A)
        dims= dims1(1)
        
        

        diag=.true. !diagonale initialisé à vrai

        do i=1,dims
          somme=0.0
          do j=1,dims
            if (i/=j) then
              somme=somme+abs(A(i,j)) !si pas à diagonale dominante alors faux
            endif
          enddo
          if (abs(A(i,i)) < somme) then
            diag=.false.
          endif
        enddo
      endsubroutine diagonale
      !Calcul de la matrice de Hilbert
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
      !Calcul du vecteur B
      subroutine CalculB(A,B,n)
        real,dimension(:,:),intent(in) :: A
        real,dimension(:), intent(out) :: B
        integer, intent(in) :: n
        integer :: i
        
        do i=1,n
          B(i)=sum(A(i,:))
        enddo
      endsubroutine CalculB  











                
