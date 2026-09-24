	clc
	clear
	format long
	n=input('Donnez la taille de la matrice n');
	w=input('Donnez le omega');
	eps=input('Donnez le epsilon');
	iter=input('Donner le nombre d''itérations');
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
	
	x0=zeros(N,1);
	m=0;
	erreur= eps+1;


	% Décomposer A en D, L, U
	D = diag(diag(A));           % Partie diagonale
	L = tril(A, -1);             % Partie strictement inférieure
	U = triu(A, 1);              % Partie strictement supérieure

	% Calcul de F et G pour Gauss-Seidel
	F = -(D + L) \ U;            % F = -(D+L)^(-1) * U
	G = (D + L) \ B;             % G = (D+L)^(-1) * B

	% Itération
	x1 = x0;
	for i = 1:n
	    m=m+1
	    x1 = F * x0 + G;         % Nouvelle approximation
	    erreur = norm(x1 - x0);     % Résidu relatif
	    if erreur <= eps
		 disp('Convergence atteinte en %d itérations.\n', m);
		 break;
	    end
	    x0 = x1;                 % Mise à jour de l'itération précédente
	end

	% Résultat final
	disp('Solution approximative :');
	disp(x1);
	
	if erreur> eps
		disp('Il n''y a pas convergence, le nombre d''itération n''est peut-être pas suffisant')
	end

