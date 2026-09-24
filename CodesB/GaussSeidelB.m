	clc
	clear
	format long
	n=input('Donnez la taille de la matrice n');
	w=input('Donnez le omega');
	Eps=input('Donnez le epsilon');
	iter=input('Donner le nombre d''iterations');
	A = hilb(n)+ w*eye(n);
	B=sum(A,2);
	
	%vérification que la matrice soit bien carrée
	
	[M,N]=size(A);
	if M~=N
		disp('Pas de solution, mettez une matrice carrée')
		return
	end
	
	%Calcul du déterminant, il faut qu'il soit différent de 0
	
	D=det(A)
	if (D==0)
		disp('Pas de solution, le déterminant est égal à 0')
		return
	end
	
	%vérification que la matrice soit bien à diagonale dominante
	
	for i=1:N
		SommeI= sum(abs(A(i, 1:i-1)))+sum(abs(A(i, i+1:N)));
		if abs(A(i,i))<SommeI
			disp('la matrice n''est pas à diagonale dominante')
			return
		end
	end
	
	%début Gauss Seidel
	
	x=zeros(N,1);
	m=0;
	erreur= Eps+1;
	
	for iter=1:n
		if erreur> Eps
			m=m+1
			for i=1:N
				y(i)=B(i);
				for j =1:i-1
					y(i)= y(i)-A(i,j)*x(j);
				end
				for j=i+1:N
					y(i)=y(i)-A(i,j)*x(j);
				end
				y(i)=y(i)/A(i,i);
			end
			erreur=max(abs(y-x));
			x=y;
		end
	end
	
	disp(x)
	
	% Avertissement convergence
	if erreur>Eps 
		disp('Il n''y a pas convergence, le nombre d''itération n''est peut-être pas suffisant')
	end
	

