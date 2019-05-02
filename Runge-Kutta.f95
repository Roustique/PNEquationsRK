module RungeKuttaPNEquations
use my_prec
implicit none
real(mp), parameter :: c=1.0 !Speed of light 
real(mp), parameter :: G=1.47662e-1 !Gravitational Parameter
contains

subroutine EqMotPostNewt(x,v,Mass,k)
real(mp), intent(in), dimension(2) :: x, v
real(mp), intent(in) :: Mass
real(mp), intent(out), dimension(4) :: k
real(mp) :: C1, C2, r

    r=sqrt(x(1)**2+x(2)**2)
    
    C1=-(c*c+v(1)*v(1)+v(2)*v(2)-4*G*Mass/r)
    C2=4*(v(1)*x(1)+v(2)*x(2))

    k(1)=v(1)
    k(2)=v(2)
    k(3)=(C1*x(1)+C2*v(1))*G*Mass/(c*c*r**3)
    k(4)=(C1*x(2)+C2*v(2))*G*Mass/(c*c*r**3)

end subroutine EqMotPostNewt

subroutine EqMotNewt(x,v,Mass,k)
real(mp), intent(in), dimension(2) :: x, v
real(mp), intent(in) :: Mass
real(mp), intent(out), dimension(4) :: k
real(mp) :: gradphi, r

	r=sqrt(x(1)**2+x(2)**2)

	gradphi=-G*Mass/r**3

    k(1)=v(1)
	k(2)=v(2)
	k(3)=gradphi*x(1)
	k(4)=gradphi*x(2)

end subroutine EqMotNewt

subroutine RKIntegrate(step, amount, Coord0, Velocity0, Mass, Coord, Velocity, Energy, Angmom, Key)
real(mp), intent(in) :: step, Mass
integer(4), intent(in) :: amount
real(mp), intent(in), dimension(2) :: Coord0, Velocity0   !Initial Values 
real(mp), intent(out), dimension(0:amount,2) :: Coord, Velocity
real(mp), intent(out), dimension(0:amount) :: Energy, Angmom
real(mp) :: r  !r is module of Coord vector
real(mp) :: v2 !v2 is scalar square of velocity
real(mp), dimension(4,4) :: k  !RK4 coefficients
integer(4) :: i
character(2), intent(in) :: Key

    Coord(0,:)=Coord0(:)
    Velocity(0,:)=Velocity0(:)
    r=sqrt(sum(Coord0**2))
	v2=sum(Velocity(0,:)**2)

	if (Key == "-P") then

        Energy(0)=v2/2-G*Mass/r+(1/c**2)*(3/8*v2**2+G*Mass/2/r*3*v2+(G*Mass/r)**2/2)
		Angmom(0)=(Coord0(1)*Velocity0(2)-Coord0(2)*Velocity0(1))*(1+v2/2/c**2+3*G*Mass/r/c**2)

	elseif (Key == "-N") then

		Energy(0)=v2/2-G*Mass/r
		Angmom(0)=Coord0(1)*Velocity0(2)-Coord0(2)*Velocity0(1)

	endif
    
    do i=1,amount
       
		if (Key=="-P") then

        	call EqMotPostNewt(Coord(i-1,:),Velocity(i-1,:),Mass,k(1,:))
        	call EqMotPostNewt(Coord(i-1,:)+step*k(1,:2)/2,Velocity(i-1,:)+step*k(1,3:)/2,Mass,k(2,:))
        	call EqMotPostNewt(Coord(i-1,:)+step*k(2,:2)/2,Velocity(i-1,:)+step*k(2,3:)/2,Mass,k(3,:))
        	call EqMotPostNewt(Coord(i-1,:)+step*k(3,:2),Velocity(i-1,:)+step*k(3,3:),Mass,k(4,:))

		elseif (Key=="-N") then

			call EqMotNewt(Coord(i-1,:),Velocity(i-1,:),Mass,k(1,:))
        	call EqMotNewt(Coord(i-1,:)+step*k(1,:2)/2,Velocity(i-1,:)+step*k(1,3:)/2,Mass,k(2,:))
        	call EqMotNewt(Coord(i-1,:)+step*k(2,:2)/2,Velocity(i-1,:)+step*k(2,3:)/2,Mass,k(3,:))
        	call EqMotNewt(Coord(i-1,:)+step*k(3,:2),Velocity(i-1,:)+step*k(3,3:),Mass,k(4,:))

		else

			stop("Wrong Argument! (Use '-N' or '-P')")

		endif

		Coord(i,:)=Coord(i-1,:)+step/6*(k(1,:2)+2*k(2,:2)+2*k(3,:2)+k(4,:2))
        Velocity(i,:)=Velocity(i-1,:)+step/6*(k(1,3:)+2*k(2,3:)+2*k(3,3:)+k(4,3:))
 
        r=sqrt(sum(Coord(i,:)**2))
		v2=sum(Velocity(i,:)**2)

		if (Key == "-P") then

        	Energy(i)=v2/2-G*Mass/r+(1/c**2)*(3/8*v2**2+G*Mass/2/r*3*v2+(G*Mass/r)**2/2)
			Angmom(i)=(Coord(i,1)*Velocity(i,2)-Coord(i,2)*Velocity(i,1))*(1+v2/2/c**2+3*G*Mass/r/c**2)

		elseif (Key == "-N") then

			Energy(i)=v2/2-G*Mass/r
			Angmom(i)=Coord(i,1)*Velocity(i,2)-Coord(i,2)*Velocity(i,1)

		endif

    enddo

end subroutine RKIntegrate


end module RungeKuttaPNEquations
