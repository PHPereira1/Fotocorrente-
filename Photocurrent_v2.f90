module var
!---------------------------------------------------
!                 input
!-----------------------------------------------------
integer, parameter   :: tn=2**16 ! numero de pontos na rede do tempo
real*8, parameter    :: pi = 3.14159265359d0
real*8, parameter    :: me = 9.10938356E-31 ![kg]
real*8, parameter    :: hbar = 1.054571817E-34 ![J][s]
real*8, parameter    :: e_charge = 1.602176634E-19 ![C]

real*8, parameter    :: eo_InGaAs = 14.3d0 !constante dieletrica InGaAs
real*8, parameter    :: eo_InAlAs = 12.7d0 !constante dieletrica InAlAs
real*8, parameter    :: m_InGaAs = 0.0436d0*me 
real*8, parameter    :: m_InAlAs = 0.0836d0*me 
real*8, parameter    :: nm = 1E-9 ![m]
real*8, parameter    :: Vo=0.503d0*e_charge ![J]
						
					 !fonte:DOI: 10.1103/PhysRevB.67.085318	
real*8, parameter    :: Nonparabw = 1.3E-18
real*8, parameter    :: Ew = hbar*hbar/(2*m_InGaAs*Nonparabw)    ! gap energy of InGaAs in J!
real*8, parameter    :: Eb = Ew*m_InAlAs/m_InGaAs                ! gap energy of InAlAs in J !


real*8, parameter    :: a0_InGaAs = 0.53d0*eo_InGaAs/m_InGaAs 
real*8, parameter    :: a0_InAlAs = 0.53d0*eo_InAlAs/m_InAlAs 


real*8, parameter    :: tol=1.0d-4 ! tolerancia

complex*8, parameter :: ic = (0.d0,1.d0)
real*8, allocatable  :: v(:),  x(:), m(:)
END MODULE var
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
program split
use var
implicit none
integer        :: j, jn, i, En, ne, x_cut, num_QB_left, num_QB_right, n, nn
real*8         :: ener,dE, E_initial, E_final, E0, E0a, fs, start, finish, dx
real*8 	       :: m_qw, m_qb, efs_InGaAs, efs_InAlAs, a0, ry, Fd, soma
real*8         :: L0_left, L_qb_left, L_qw_left, L_bm_left, L_left, L_cqw, L_qb_right, L_qw_right, L_bm_right, L0_right, L_right
double complex :: aux1, beta, Je, Jd, trap
double complex :: P_inicial(2,2), P_final(2,2), Dn(2,2), T(2,2), AUX(2,2), AUX2(2,2), T0(2,2), Ti(2,2)
double complex :: W(2,2), W0(2,2), Wi(2,2), Wf(2,2), F(2), beta_minus, beta_plus, cte, ctecur_left, ctecur_right
character*1024   :: char_QBl, char_Lqw_l, char_Lqb_l, char_cqw, char_QBr, char_Lqw_r, char_Lqb_r, SL
character (len = 50) :: SL_info, SL_folder

real*8, allocatable  		 :: energy(:)
double complex, allocatable  :: psi(:), prod(:), ki(:), psi0(:,:) 

m_qw = m_InGaAs
m_qb = m_InAlAs

! ------------------------------------------------------------------
!                 dimension of the heterostructure [m]
! ------------------------------------------------------------------
! quantum bragg mirror to left of central quantum well
num_QB_left  = 5		!Number of well/barrier of the Bragg mirror
L_qw_left    = 2.0d0 	!Width of lateral quantum well [nm]
L_qb_left    = 7.0d0 	!Width of quantum barrier [nm]

! central quantum well
L_cqw        = 2.5d0 	!Width of central quantum well [nm]

! quantum bragg mirror to right of central quantum well
num_QB_right = 1 		!Number of well/barrier of the Bragg mirror
L_qw_right   = 2.0d0 	!Width of lateral quantum well [nm]
L_qb_right   = 7.0d0 	!Width of quantum barrier [nm]



dx           = 0.1*nm   ! step size dx
nn           = 50       ! The number of states you want to calculate
!------------------------------------------------------------------

! quantum bragg mirror to left of central quantum well
write(char_QBl,'(I1)') num_QB_left
write(char_Lqw_l,'(f4.2)')L_qw_left
L_qw_left    = L_qw_left*nm 	!Width of lateral quantum well
write(char_Lqb_l,'(f4.2)')L_qb_left
L_qb_left    = L_qb_left*nm 	!Width of quantum barrier
L_bm_left    = -num_QB_left*(L_qb_left+L_qw_left)
L0_left      = -50*nm		!Thickness of the first barrier 
L_left       = (L_BM_left + L0_left)


! central quantum well
write(char_cqw,'(f4.2)')L_cqw
L_cqw        = L_cqw*nm 	!Width of central quantum well

! quantum bragg mirror to right of central quantum well
write(char_QBr,'(I1)') num_QB_right
write(char_Lqb_r,'(f4.2)')L_qb_right
L_qb_right   = L_qb_right*nm 	!Width of quantum barrier
write(char_Lqw_r,'(f4.2)')L_qw_right
L_qw_right   = L_qw_right*nm 	!Width of lateral quantum well
L_bm_right   = num_QB_right*(L_qb_right+L_qw_right)
L0_right     = -L_left - L_cqw - L_BM_right !Thickness of the last barrier
L_right      = L_cqw + L_BM_right + L0_right + L_cqw



n			 = int((L_right-L_left)/dx) + 2

allocate (v(n),  x(n), m(n), psi(n), prod(n), ki(n), psi0(n,nn), energy(nn))
!------------------------------------------------------------------

!------------------------------------------------------------------
! Creating the filename with info of the suprlattice
SL_info = trim(char_QBl)//"x("//trim(char_Lqw_l)//','//trim(char_Lqb_l)//')-'//trim(char_cqw)//'-'&
&//trim(char_QBr)//'x('//trim(char_Lqw_r)//','//trim(char_Lqb_r)//')' 

! Creating the filename of the folder of the results 
SL_folder = trim(char_QBl)//"x"//trim(char_Lqw_l)//"+"//trim(char_Lqb_l)//"-"//trim(char_cqw)//'-'//&
&trim(char_QBr)//"x"//trim(char_Lqw_r)//"+"//trim(char_Lqb_r)

! Making the folder where the results will be saved 
call system('mkdir -p Structures/'//adjustl(trim(SL_folder)))

open(26, file='Structures/'//adjustl(trim(SL_folder))//'/wavefunction_SL_'//trim(SL_info)//'.txt')
open(7,  file='Structures/'//adjustl(trim(SL_folder))//'/Transmission_SL_'//trim(SL_info)//'.txt') 
open(8,  file='Structures/'//adjustl(trim(SL_folder))//'/Photocurrent_SL_'//trim(SL_info)//'.txt')
open(25, file='Structures/'//adjustl(trim(SL_folder))//'/Potencial_SL_'//trim(SL_info)//'.txt') 
open(28, file='Structures/'//adjustl(trim(SL_folder))//'/OscStr_SL_'//trim(SL_info)//'.txt')




do j=1,n      
	  x(j)=((j-1)*dx + L_left) ! vetor x da super-rede
end do

call mass(m_qw, m_qb, num_qb_left, L_qw_left, L_qb_left, L_cqw, num_qb_right, L_qw_right, L_qb_right, n)
fs= 0*efs_InGaAs
call potencial(fs, num_qb_left, L_qw_left, L_qb_left, L_cqw, num_qb_right, L_qw_right, L_qb_right, x_cut, n)
  
do j=1,n
	write(25,*) x(j)/nm, 1E3*v(j)/e_charge!, m(j)/me!
end do	
!---------------------------------------------------------------------
!   Computing the wavefunction of the ground state using the Numerov method
!---------------------------------------------------------------------
call cpu_time(start)

E_initial = 0.d0   ![meV] 
E_final   = 800.d0 ![meV] 
dE        = 0.1    ![meV]
En        = (E_final - E_initial)/dE
Ener      = E_initial


call wavefunction_numerov(num_qb_left, L_qw_left, L_qb_left, L_cqw, num_qb_right, L_qw_right, L_qb_right, x_cut, dx, & 
                          E_initial, E_final, dE, En, ne, Ener, psi0, energy, n, nn)

psi(:) = psi0(:,1) ! wavefucntion of the ground state
E0 = Energy(1)![J] energy of the ground state

!	---------------------------------------------------------------
!   Escrevendo as funcoes de onda do eletron
do i=1,n
	write(26,*) real(x(i)/nm), ((real(psi0(i,j)**2) + imag(psi0(i,j)**2))*400*nm + 1E3*Energy(j)/e_charge, j=1,ne)
end do
!	------------------------------------------------------

!------------------------------------------------------------------------
!     Calculanting the oscillator strength
call oscillator_strength(psi0, energy, ne, dx, n)
!------------------------------------------------------------------------

!	---------------------------------------------------------------------

Fd = 1*1E5*e_charge ! Amplitude of the external electric field [J]/[M]

!---------------------------------------------------------------------
!   Computing the photovoltaic photocurrent spectrum using the coherent carrier propagation in the continuum 
!   doi: https://doi.org/10.1103/PhysRevB.60.R13993
!---------------------------------------------------------------------

E_initial = 200.d0 ![meV] 
E_final   = 500.d0 ![meV] 
dE        = 0.1    ![meV]
En        = (E_final - E_initial)/dE
do  jn=1,En     
	Ener= (E_initial + jn*dE) ![meV]
	Ener= Ener*e_charge*1E-3  ![J]
	
	! non-parabolicity 	  
    m_qw = m_InGaAs*(1.d0 + ener/Ew)
    m_qb = m_InAlAs*(1.d0  + (E0+ ener-vo)/Eb)	  
    call mass(m_qw, m_qb, num_qb_left, L_qw_left, L_qb_left, L_cqw, num_qb_right, L_qw_right, L_qb_right, n) 	
	
	T0 (1,1) = 1.d0 
	T0 (1,2) = 0.d0
	T0 (2,1) = 0.d0
	T0 (2,2) = 1.d0
	
	Wf (1,1) = 0.d0 
	Wf (1,2) = 0.d0
	Wf (2,1) = 0.d0
	Wf (2,2) = 0.d0
	
	W (1,1) = 1.d0 
	W (1,2) = 0.d0
	W (2,1) = 0.d0
	W (2,2) = 1.d0
	
	W0 (1,1) = 1.d0 
	W0 (1,2) = 0.d0
	W0 (2,1) = 0.d0
	W0 (2,2) = 1.d0
	
	F(1) = 0.d0
	F(2) = 0.d0
	
	do i=n-1,1, -1
		aux1 = (2*m(i)*(Ener+E0-v(i)))/hbar/hbar
		ki(i) = cdsqrt(aux1)

		aux1 = (2*m(i+1)*(Ener+E0-v(i+1)))/hbar/hbar
		ki(i+1) = cdsqrt(aux1)	
		

		!Matriz propagacao para i=1
		P_inicial(1,1) = zexp(-ic*ki(i+1)*(x(i))) 
		P_inicial(1,2) = 0.d0
		P_inicial(2,1) = 0.d0
		P_inicial(2,2) = zexp(+ic*ki(i+1)*(x(i)))
		
		P_final(1,1) = zexp(+ic*ki(i)*(x(i))) 
		P_final(1,2) = 0.d0
		P_final(2,1) = 0.d0
		P_final(2,2) = zexp(-ic*ki(i)*(x(i)))		
	
		!Matriz descontinuidade para i=1 
		beta = (ki(i)*m(i+1))/(ki(i+1)*m(i))
		Dn(1,1) = 0.5d0*(1 + beta)
		Dn(1,2) = 0.5d0*(1 - beta)
		Dn(2,1) = 0.5d0*(1 - beta)
		Dn(2,2) = 0.5d0*(1 + beta)
		
		
		Call  matrix_mult(P_inicial,Dn,AUX)
		Call  matrix_mult(AUX,P_final,Ti)
		Call matrix_mult(T0,Ti,T)		
		T0 = T	
		
!		----------------------------------------------------
!		Calculo da interacao eletron-foton na interface i
!       ----------------------------------------------------
		cte = (m(i)*Fd)/(2.d0*ic*hbar*hbar*ki(i))
		beta_plus  = -cte*(zexp(+ic*ki(i)*x(i))*x(i)*psi(i)-zexp(+ic*ki(i)*x(i+1))*x(i+1)*psi(i+1))*dx
		beta_minus = +cte*(zexp(-ic*ki(i)*x(i))*x(i)*psi(i)-zexp(-ic*ki(i)*x(i+1))*x(i+1)*psi(i+1))*dx
		
		F(1) = F(1) + (T(1,1)*beta_minus + T(1,2)*beta_plus)
		F(2) = F(2) + (T(2,1)*beta_minus + T(2,2)*beta_plus)  	
				
	end do
	
!	------------------------------------------------------------
!	Calculo da fotocorrente em funcao da energia do foton (E-E0)
!   ------------------------------------------------------------

ctecur_left  = (hbar*ki(1))/m(1)
ctecur_right = (hbar*ki(n))/m(n)

Je = ctecur_left*Fd*(-F(2)/T(2,2)*dconjg(-F(2)/T(2,2)))
Jd = ctecur_right*Fd*(F(1)-F(2)*T(1,2)/T(2,2))*dconjg((F(1)-F(2)*T(1,2)/T(2,2)))	

write(8,*)1E3*Ener/e_charge,real(Jd-Je)!/6.966529605e-23	
write(7,*) 1E3*Ener/e_charge, real(1/cdabs(T(1,1))) !real(T(1,1))!

enddo
call cpu_time(finish)
write(*,*) finish-start
open(30, file='fimprograma.txt')
stop
end

!----------------------------------------------------------------
!		Multiple matrix 2x2
!----------------------------------------------------------------

subroutine matrix_mult(A,B,C)
	implicit none
	double complex :: A(2,2), B(2,2), C(2,2)
	
	C(1,1) = A(1,1)*B(1,1) + A(1,2)*B(2,1)
	C(1,2) = A(1,1)*B(1,2) + A(1,2)*B(2,2)
	C(2,1) = A(2,1)*B(1,1) + A(2,2)*B(2,1)
	C(2,2) = A(2,1)*B(1,2) + A(2,2)*B(2,2)	
	return
end subroutine matrix_mult
	  
	  
!---------------------------------------------------------------
!          Rotina para o calculo das funcoes de onda utilizando o metodo numerov 
!---------------------------------------------------------------  
  
 subroutine wavefunction_numerov(num_qb_left, L_qw_left, L_qb_left, L_cqw, num_qb_right, L_qw_right, L_qb_right, x_cut, dx, &
	 							Ei, Ef, dE, jn, ne, Ener, psi0, energy, n, nn)
														
	 use var
	 Implicit none
	 Integer :: j, x_cut, jn, k, ne, num_qb_left, num_qb_right, En, n, nn			
	 Real*8  :: Ei, Ef, dE, Ener, Ener_Old, Ener_New, Energy(nn), mult1, m_qw, m_qb, soma 
	 real*8  :: dx, L_qw_left, L_qb_left, L_cqw, L_qw_right, L_qb_right
	 Double complex :: psi(n), psi_Old, psi_New, prod(n), trap, psi0(n,nn) 
	 
	 
	 ne = 0
     m_qw = m_InGaAs*(1.d0 + ener/Ew)
     m_qb = m_InAlAs*(1.d0  + (ener-vo)/Eb)
	 call mass(m_qw, m_qb, num_qb_left, L_qw_left, L_qb_left, L_cqw, num_qb_right, L_qw_right, L_qb_right, n) 	
	 call numerov(Ener, psi, dx, n)
	 psi_Old = psi(n)
	 Ener_old = Ener

	 do  j=1,jn   ! numero de autovalores do problema 
	    psi_Old = psi(n)
!   calculando a funcao de onda do eletron para energia Ener  
	    Ener = Ei + dE*j!  
	    Ener = Ener*e_charge*1E-3  ![J]
        m_qw = m_InGaAs*(1.d0 + ener/Ew)
        m_qb = m_InAlAs*(1.d0  + (ener-vo)/Eb)
        call mass(m_qw, m_qb, num_qb_left, L_qw_left, L_qb_left, L_cqw, num_qb_right, L_qw_right, L_qb_right, n)	
	    call numerov(Ener, psi, dx, n)
	    psi_New = psi(n)
		
!   condicao necessaria para o "do while" operar (procure a regiao em que a multiplicacao passou por zero)	
	    mult1 = real(psi_old*psi_New)
	    if(mult1.le.0.d0)then
	    	Ener_new = Ener
	    	k = 0
	    15  continue 
			
!	------------------------------------------------------
!	calculando o funcao de onda do eletron com energia Ener_Old
	    Ener = Ener_old	  
        m_qw = m_InGaAs*(1.d0 + ener/Ew)
        m_qb = m_InAlAs*(1.d0  + (ener-vo)/Eb)
   	    call mass(m_qw, m_qb, num_qb_left, L_qw_left, L_qb_left, L_cqw, num_qb_right, L_qw_right, L_qb_right, n) 		
	    call numerov(Ener, psi, dx, n)
	    psi_old = psi(n)
		
!	------------------------------------------------------
!	calculando o funcao de onda do eletron com energia Ener_New
	    Ener = Ener_new	  
        m_qw = m_InGaAs*(1.d0 + ener/Ew)
        m_qb = m_InAlAs*(1.d0 + (ener-vo)/Eb)
   	    call mass(m_qw, m_qb, num_qb_left, L_qw_left, L_qb_left, L_cqw, num_qb_right, L_qw_right, L_qb_right, n) 	
	    call numerov(Ener, psi, dx, n)
	    psi_new = psi(n)
!	------------------------------------------------------
!	condicao para saber se as psiN (ener_old e ener_new) passar por zero ou nao
	    mult1 = real(psi_old*psi_new)
!	------------------------------------------------------	

!	------------------------------------------------------	
!	Inicio do metodo de bissecao para encontrar os autovalores da funcao de onda		
	    if (mult1.lt.0) then
!caso a condicao nao tenha passado por zero, o valor ener sera a media entre ener_old e ener_new	
	    	Ener = (Ener_old + Ener_new)/2.0
	    end if
!	------------------------------------------------------
!	calculando o funcao de onda para nova energia Ener	  
		m_qw = m_InGaAs*(1.d0 + ener/Ew)
		m_qb = m_InAlAs*(1.d0  + (ener-vo)/Eb)
		call mass(m_qw, m_qb, num_qb_left, L_qw_left, L_qb_left, L_cqw, num_qb_right, L_qw_right, L_qb_right, n)  	
	    call numerov(Ener, psi, dx, n)
	    psi_New = psi(n)

!	------------------------------------------------------
!	condicao para saber se as psiN (ener_old e ener_new) passou por zero ou nao	
	    mult1 = real(psi_old*psi_New)
!	------------------------------------------------------	
	    if (mult1 .lt. 0) then
	    	Ener_new = Ener 
	    	else
	    	Ener_old = Ener
	    end if
	    k = k + 1 !contador de vezes que o metodo da bissecao foi utilizado
	    if (k.eq.100) go to 16 !evitando que o loop se torne infinito 
!		Caso a amplitude a funcao de onda em psi(n) for menor que a tolerancia, o problema foi resolvido 	
	    if (abs(psi_new) .gt. tol) goto 15 
	
 16  	continue
 ! Quando ocorre o loop infinito significa que a funcao de onda está explodindo no final da estrutura
!Para contornar esse problema, o calculo da funcao de onda precisa ser dividido em duas partes: um loop que vai do 
!inicio da estrutura i = 0 até no meio da estrutura n/2 e outro loop saindo do final da estrutura 
! no ponto n ate no meio da estrutura n/2
  		ne = ne + 1
   		if(Ener.lt.vo)then
	        m_qw = m_InGaAs*(1.d0 + ener/Ew)
	        m_qb = m_InAlAs*(1.d0  + (ener-vo)/Eb)
	   		call mass(m_qw, m_qb, num_qb_left, L_qw_left, L_qb_left, L_cqw, num_qb_right, L_qw_right, L_qb_right, n)
	    	call funcWave_Numerov2(Ener, dx, x_cut, psi, ne, n) 
	    end if
	
	    !write(*,*)"The energy is", 1E3*Ener/e_charge
	    psi0(:,ne) =  psi(:)
	    Energy(ne) =  Ener
		
	    if (ne.eq.nn) go to 17 ! se o numero de autoestados encontrados (ne) for igual ao numero de autoestado desejados(nn), finalize o loop
! Aqui estou calculando novamente a funcao de onda completa pois quando uso o numerov2 eu coloco psi nas bordas é zero e com isso 
! a condicao mult1 de calculo da funcao de onda nao é sastifeita
		Ener = Ei + dE*j
		m_qw = m_InGaAs*(1.d0 + ener/Ew)
		m_qb = m_InAlAs*(1.d0 + (ener-vo)/Eb)
		call mass(m_qw, m_qb, num_qb_left, L_qw_left, L_qb_left, L_cqw, num_qb_right, L_qw_right, L_qb_right, n)  		
	    call numerov(Ener, psi, dx, n)
	  	end if 	
	end do
17 continue  
!---------------------------------------------------------------
!	Normalizando as funcoes de onda	
	do j=1,ne
	    psi(:) = psi0(:,j)
	    call norm_psi(psi, 1, n, dx, n)
	    psi0(:,j) = psi(:)
	end do	
	
	return 
end subroutine wavefunction_numerov	


!----------------------------------------------------------------
!		Integracao da equacao de Schrodinger pelo metodo Numerov
!----------------------------------------------------------------

	subroutine numerov(Ener,psi, dx, n)
	use var
	implicit none
	double complex :: psi(n), aux1, ki(n), a, b, c
	real*8         :: ener, beta, am, dx
	integer 	   :: i, n
	
	
	psi(1) = 0.d0
	psi(2) = 0.1d0
	beta = dx*dx/12.d0
	
	do i=2,n-1
		
		am = (m(i-1) + m(i))/2.d0 ! média das massas do eletron
		aux1 = (2*am*(Ener-v(i-1)))/hbar/hbar
		ki(i-1) = cdsqrt(aux1)
		
		am = (m(i) + m(i+1))/2.d0 ! média das massas do eletron
		aux1 = (2*am*(Ener-v(i)))/hbar/hbar
		ki(i) = cdsqrt(aux1)
		
		am = (m(i+1) + m(i))/2.d0 ! média das massas do eletron
		aux1 = (2*am*(Ener-v(i+1)))/hbar/hbar
		ki(i+1) = cdsqrt(aux1)
		
		a = 2.d0 - 10.d0*beta*ki(i)**2
		b = 1.d0 + beta*ki(i-1)**2
		c = 1.d0 + beta*ki(i+1)**2
		
		psi(i+1) = (a*psi(i) - b*psi(i-1))/c			
	end do
	
	return
	end subroutine	numerov
	
	
	!----------------------------------------------------------------
	!Integracao da equacao de Schrodinger pelo metodo Numerov para 
	! quando há explosao da funcao de onda no extremo 
	!----------------------------------------------------------------	
	subroutine funcWave_Numerov2(Ener, dx, x_cut, psi, ne, n)
	use var
	implicit none
	double complex :: psi(n), psi_left(n), psi_right(n), aux1, ki(n), a, b, c
	real*8         :: ener, beta, am, k, dx
	integer 	   :: i, ne, x_cut, n
	
	psi_left(1) = 0.d0
	psi_left(2) = 0.1d0
	beta = dx*dx/12.d0
	

	! calculando a funcao de onda do inicio da estrutura até o meio da estrutura
	do i=2, x_cut + 1
	
		am = (m(i-1) + m(i))/2.d0 ! média das massas do eletron no ponto i-1
		aux1 = (2*am*(Ener-v(i-1)))/hbar/hbar
		ki(i-1) = cdsqrt(aux1)
		
		am = (m(i) + m(i+1))/2.d0 ! média das massas do eletron no ponto i
		aux1 = (2*am*(Ener-v(i)))/hbar/hbar
		ki(i) = cdsqrt(aux1)
	
		am = (m(i+1) + m(i+2))/2.d0 ! média das massas do eletron no ponto i+1
		aux1 = (2*am*(Ener-v(i+1)))/hbar/hbar
		ki(i+1) = cdsqrt(aux1)
	
		a = 2.d0 - 10.d0*beta*ki(i)**2
		b = 1.d0 + beta*ki(i-1)**2
		c = 1.d0 + beta*ki(i+1)**2
	
		psi_left(i+1) = (a*psi_left(i) - b*psi_left(i-1))/c
		
	end do

	psi_right(n) = 0.d0
	psi_right(n-1) = 0.1d0

	! calculando a funcao de onda do final da estrutura até o meio da estrutura
	do i=n-1,x_cut + 1, -1 
	
		am = (m(i-2) + m(i-1))/2.d0 ! média das massas do eletron no ponto i-1
		aux1 = (2*am*(Ener-v(i-1)))/hbar/hbar
		ki(i-1) = cdsqrt(aux1)
	
		am = (m(i) + m(i-1))/2.d0 ! média das massas do eletron no ponto i
		aux1 = (2*am*(Ener-v(i)))/hbar/hbar
		ki(i) = cdsqrt(aux1)
	
		am = (m(i) + m(i+1))/2.d0 ! média das massas do eletron no ponto i+1
		aux1 = (2*am*(Ener-v(i+1)))/hbar/hbar
		ki(i+1) = cdsqrt(aux1)
	
		a = 2.d0 - 10.d0*beta*ki(i)**2
		b = 1.d0 + beta*ki(i-1)**2
		c = 1.d0 + beta*ki(i+1)**2
	
		psi_right(i-1) = (a*psi_right(i) - c*psi_right(i+1))/b	
				
	end do
	k = ne + 1
	
	!A funcao de onda a direita apresenta uma paridade inversa que a funcao de onda da esquerda
	!para resolver o problema estou multiplicando por -1 as funcoes de onda impares
	do i=n-1,x_cut + 1, -1 
		psi_right(i-1) = (-1)**k*psi_right(i-1)
	end do
	
	!normalizando a funcao de onda da esquerda com a amplitude da funcao de onda da direita
	do i = 1,x_cut
		psi_left(i) = psi_left(i)*psi_right(x_cut)/psi_left(x_cut)
	end do
	
	!juntando as funcoes de onda esquerda e direita em apenas uma funcao de onda
	do i = 1,x_cut
		psi(i) = psi_left(i)
	end do
	
	do i = x_cut,n 
		psi(i) = psi_right(i)
	end do
	
	return 
	end subroutine funcWave_Numerov2
	
	
	!----------------------------------------------------------------
	!		Normalizando as funcoes de onda
	!----------------------------------------------------------------
	
		subroutine norm_psi(psi, ni, nf, dx, n)
			use var
			implicit none
			double complex :: psi(n), prod(n), trap
			real*8         :: norm, dx 
			integer 	   :: ni, nf, j, n
	 	
			do j = ni, nf
				prod(j) = conjg(psi(j))*psi(j)
			end do
		
			norm = trap(prod, ni, nf, dx, n)
		
			do j = ni, nf
				psi(j) = psi(j)/sqrt(norm)
			end do
		return 
		end subroutine norm_psi
	
	
	!---------------------------------------------------------------------
	!       Integracao - usando o metodo do trapezio 
	!---------------------------------------------------------------------
		Function trap(d, ni, nf, dx, n)
			use var
			Implicit none
			Integer        :: ni, nf, j, n
			double complex :: d(n), trap
			real*8		   :: dx			
		
			trap = 0

			do j=ni+2,(nf-1)
				trap = trap +d(j)*dx
			end do
		
			trap = trap+0.5*(d(ni)+d(nf))*dx
			return
		end function trap
		
		
!----------------------------------------------------------------
!		Calculando a forca de oscilador
!----------------------------------------------------------------
	subroutine oscillator_strength(psi0, energy, ne, dx, n)
		use var
		implicit none
		integer 	   :: i, ne, n
		double complex :: psi0(n,ne), prod(n), trap
		real*8         :: norm, energy(ne), f(ne),dE, soma, dx
				
		do i = 1, ne  
			prod(:) = conjg(psi0(:,i))*x(:)*psi0(:,1)
			dE = (Energy(i) - Energy(1))
			soma = trap(prod, 1, n, dx, n)
			f(i) = 2*me*dE*abs(soma)**2/hbar/hbar
			write(28,*) 1E3*dE/e_charge, f(i)
		end do
		
	return 
	end subroutine oscillator_strength		
	

!---------------------------------------------------------------
!                    mass
!---------------------------------------------------------------
	subroutine mass(m_qw, m_qb, num_qb_left, L_qw_left, L_qb_left, L_cqw, num_qb_right, L_qw_right, L_qb_right, n)  
		use var
	    Implicit none
	  	Integer :: i, j, k, num_qb_left, num_qb_right, n_qb_left, n_qb_right, n
	    Real*8  :: m_qw, m_qb, L_qw_left, L_qb_left, L_cqw, L_qw_right, L_qb_right
	    Real*8  :: L0, L0_left, L_bm_left
		
		L_bm_left = num_QB_left*(L_qb_left+L_qw_left)    ! Comprimento do espelho de Bragg a esquerda do poco quantico de defeito
		L0        =  L_bm_left							 ! L0 é a comprimento do ultimo poco quantico da super-rede	
		L0_left   = L0                                   ! Guardando esse valor para ser usado no looping da construcao da super-rede


!----------------------------------------------------------------------
!		Escrevendo a massa do eletron da primeira barreira quantica a esquerda da amostra
!----------------------------------------------------------------------		
		
		do j = 1, n	
			if(x(j).le.-L0)then
				m(j) = m_qb
			else 
				go to 10
			end if
		end do		
10      continue	
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!	    Escrevendo a massa do eletron da super-rede a esquerda do poco quantico central
!		Aqui a repeticao vai ser poco/barreira
!----------------------------------------------------------------------
		k = 1         ! inteiro que vai ser usada para determinar o fim da repeticao poco/barreira
		n_qb_left = 0 !numero de repeticao do conjunto poco/barreira indo de 0 até 5
		
		do i = j, n
			L0 = L0_left - n_qb_left*(L_qb_left+L_qw_left)
			if(x(i).ge.-L0.and.x(i).le.-L0+L_qw_left)then
				m(i) = m_qw
			else if(x(i).ge.-L0+L_qw_left.and.x(i).le.-L0+L_qw_left+L_qb_left)then
				m(i) = m_qb
			else
			m(i) = m_qw!m_qb
			n_qb_left = n_qb_left + 1 ! adicione uma repeticao poco/barreira no potencial
			k = k+1
			end if	
			if(k.gt.num_qb_left) go to 20 ! se o numero de repeticao (k) for igual a numero fornecido no input (um_qb_left), finalize o looping
		end do
20      continue

!----------------------------------------------------------------------
		
		
!----------------------------------------------------------------------
!		Escrevendo a massa do eletron do poco quantico central
!----------------------------------------------------------------------		

		do j = i, n
			if(x(j).ge.0.d0.and.x(j).le.L_cqw)then
				m(j) = m_qw
			else
			go to 30
		end if
	    end do		
30		continue 
!----------------------------------------------------------------------



!----------------------------------------------------------------------
!	    Escrevendo a massa do eletron da super-rede a direita do poco quantico central
!		Aqui a repeticao vai ser barreira/poco 
!----------------------------------------------------------------------		
		k = 1
		n_qb_right = 0
		do i=j,n
		L0 = L_cqw + n_qb_right*(L_qw_right+L_qb_right)
		if(x(i).ge.L0.and.x(i).le.L0+L_qb_right)then
			m(i) = m_qb
		else if(x(i).ge.L0+L_qb_right.and.x(i).le.L0+L_qb_right+L_qw_right)then
			m(i) = m_qw
		else
			m(i) = m_qb!m_qw
			n_qb_right = n_qb_right + 1
			k = k+1	
		end if
		if(k.gt.num_qb_right) go to 40	
		end do
40		continue

!----------------------------------------------------------------------



!----------------------------------------------------------------------
!		Escrevendo a massa do eletron da primeira barreira quantica a esquerda da amostra
!----------------------------------------------------------------------	
		do j=i,n
			m(j) = m_qb
	    end do
!----------------------------------------------------------------------		
	  return
	end subroutine mass
	
!---------------------------------------------------------------
!                     potencial 
!---------------------------------------------------------------
	subroutine potencial(fs, num_qb_left, L_qw_left, L_qb_left, L_cqw, num_qb_right, L_qw_right, L_qb_right, x_cut, n)
		use var
	    Implicit none
	    Integer :: j, i, x_cut, k, num_qb_left, num_qb_right, n_qb_left, n_qb_right, n	
	    Real*8  :: fs, L_cqw, L_qb,L_lqw, L_qw_left, L_qb_left, L_qw_right, L_qb_right, L_left, L0
	    Real*8  :: L0_left, L_bm_left


		L_bm_left = num_QB_left*(L_qb_left+L_qw_left)    ! Comprimento do espelho de Bragg a esquerda do poco quantico de defeito
		L0        =  L_bm_left							 ! L0 é a comprimento do ultimo poco quantico da super-rede	
		L0_left   = L0                                   ! Guardando esse valor para ser usado no looping da construcao da super-rede


!----------------------------------------------------------------------
!		Escrevendo o potencial da primeira barreira quantica a esquerda da amostra
!----------------------------------------------------------------------		
		
		do j = 1, n	
			if(x(j).le.-L0)then
				v(j) = Vo
			else 
				go to 10
			end if
		end do		
10      continue	
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!	    Escrevendo o potencial da super-rede a esquerda do poco quantico central
!		Aqui a repeticao vai ser poco/barreira
!----------------------------------------------------------------------
		k = 1         ! inteiro que vai ser usada para determinar o fim da repeticao poco/barreira
		n_qb_left = 0 !numero de repeticao do conjunto poco/barreira indo de 0 até 5
		
		do i = j, n
			L0 = L0_left - n_qb_left*(L_qb_left+L_qw_left)
			if(x(i).ge.-L0.and.x(i).le.-L0+L_qw_left)then
				v(i) = 0
			else if(x(i).ge.-L0+L_qw_left.and.x(i).le.-L0+L_qw_left+L_qb_left)then
				v(i) = vo
			else
				v(i) = 0!vo
				n_qb_left = n_qb_left + 1 ! adicione uma repeticao poco/barreira no potencial
				k = k+1
			end if	
			if(k.gt.num_qb_left) go to 20 ! se o numero de repeticao (k) for igual a numero fornecido no input (um_qb_left), finalize o looping
		end do
20      continue

!----------------------------------------------------------------------
				
!----------------------------------------------------------------------
!		Escrevendo o potencial do poco quantico central
!----------------------------------------------------------------------		

		do j = i, n
			if(x(j).ge.0.d0.and.x(j).le.L_cqw)then
				v(j) = 0
			else
			go to 30
		end if
	    end do		
30		continue 

!----------------------------------------------------------------------

	    x_cut = j ! inteiro x_cut para ser usado para calcular a funcao de onda pelo Numerov

!----------------------------------------------------------------------
!	    Escrevendo o potencial da super-rede a direita do poco quantico central
!		Aqui a repeticao vai ser barreira/poco 
!----------------------------------------------------------------------		
		k = 1
		n_qb_right = 0
		do i=j,n
		L0 = L_cqw + n_qb_right*(L_qw_right+L_qb_right)
		if(x(i).ge.L0.and.x(i).le.L0+L_qb_right)then
			v(i) = vo
		else if(x(i).ge.L0+L_qb_right.and.x(i).le.L0+L_qb_right+L_qw_right)then
			v(i) = 0.d0
		else
			v(i) = vo!0.d0
			n_qb_right = n_qb_right + 1
			k = k+1	
		end if
		if(k.gt.num_qb_right) go to 40	
		end do
40		continue

!----------------------------------------------------------------------

!----------------------------------------------------------------------
!		Escrevendo o potencial da primeira barreira quantica a esquerda da amostra
!----------------------------------------------------------------------	
		do j=i,n
			v(j) = vo
	    end do
!----------------------------------------------------------------------			 

	      return
	    end subroutine potencial
	
