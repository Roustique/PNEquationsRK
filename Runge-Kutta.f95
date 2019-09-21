module RungeKuttaPNEquations
use my_prec
implicit none
real(mp), parameter :: c=1.0 !Speed of light
real(mp), parameter :: G=1.0 !Gravitational Parameter
contains

subroutine EqMot(x,v,Mass,k,Key)
real(mp), intent(in), dimension(2) :: x, v
real(mp), intent(in) :: Mass
real(mp), intent(out), dimension(4) :: k
real(mp) :: phi, gradphi, C1, C2, C3, C4, A, B, r, v2, rv
character(3), intent(in) :: Key

    r=sqrt(x(1)**2 + x(2)**2)
    v2 = v(1)**2 + v(2)**2
    phi = - G * Mass / r
    rv = x(1) * v(1) + x(2) * v(2)

    k(1) = v(1)
    k(2) = v(2)

    if (Key == "-NT") then
        C1 = 0.0_mp
        C2 = 0.0_mp
        C3 = 0.0_mp
        C4 = 0.0_mp
    elseif (Key == "-GS") then
        C1 = 2.0_mp
        C2 = 2.0_mp
        C3 = -3.0_mp
        C4 = 2.0_mp
    elseif ((Key == "-GH").or.(Key == "-FF")) then
        C1 = 4.0_mp
        C2 = 1.0_mp
        C3 = 0.0_mp
        C4 = 4.0_mp
    elseif (key == "-F0") then
        C1 = 3.0_mp
        C2 = 1.0_mp
        C3 = 0.0_mp
        C4 = 4.0_mp
    end if

    A = phi / r**2 * (1 + C1 * phi / c**2 + C2 * v2 / c**2 + C3 * rv ** 2 / c**2 / r**2)
    B = - C4 * phi / c**2 / r**2 * rv

    k(3) = A * x(1) + B * v(1)
    k(4) = A * x(2) + B * v(2)

end subroutine EqMot

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
character(3), intent(in) :: Key

    Coord(0,:)=Coord0(:)
    Velocity(0,:)=Velocity0(:)
    r=sqrt(sum(Coord0**2))
	v2=sum(Velocity(0,:)**2)

	if ((Key == "-GS").or.(Key == "-GH").or.(Key == "-F0").or.(Key == "-FF")) then

        Energy(0)=v2/2-G*Mass/r+(1/c**2)*(3/8*v2**2+G*Mass/2/r*3*v2+(G*Mass/r)**2/2)
		Angmom(0)=(Coord0(1)*Velocity0(2)-Coord0(2)*Velocity0(1))*(1+v2/2/c**2+3*G*Mass/r/c**2)

	elseif (Key == "-NT") then

		Energy(0)=v2/2-G*Mass/r
		Angmom(0)=Coord0(1)*Velocity0(2)-Coord0(2)*Velocity0(1)

	endif

    do i=1,amount

        call EqMot(Coord(i-1,:), Velocity(i-1,:), Mass, k(1,:), Key)
        call EqMot(Coord(i-1,:) + step * k(1,:2) / 2, Velocity(i-1,:) + step * k(1,3:) / 2, Mass, k(2,:), Key)
        call EqMot(Coord(i-1,:) + step * k(2,:2) / 2, Velocity(i-1,:) + step * k(2,3:) / 2, Mass, k(3,:), Key)
        call EqMot(Coord(i-1,:) + step * k(3,:2), Velocity(i-1,:) + step * k(3,3:), Mass, k(4,:), Key)

		Coord(i,:)=Coord(i-1,:)+step/6*(k(1,:2)+2*k(2,:2)+2*k(3,:2)+k(4,:2))
        Velocity(i,:)=Velocity(i-1,:)+step/6*(k(1,3:)+2*k(2,3:)+2*k(3,3:)+k(4,3:))

        r=sqrt(sum(Coord(i,:)**2))
		v2=sum(Velocity(i,:)**2)

		if ((Key == "-GS").or.(Key == "-GH").or.(Key == "-F0").or.(Key == "-FF")) then

        	Energy(i)=v2/2-G*Mass/r+(1/c**2)*(3/8*v2**2+G*Mass/2/r*3*v2+(G*Mass/r)**2/2)
			Angmom(i)=(Coord(i,1)*Velocity(i,2)-Coord(i,2)*Velocity(i,1))*(1+v2/2/c**2+3*G*Mass/r/c**2)

		elseif (Key == "-NT") then

			Energy(i)=v2/2-G*Mass/r
			Angmom(i)=Coord(i,1)*Velocity(i,2)-Coord(i,2)*Velocity(i,1)

		endif

    enddo

end subroutine RKIntegrate


end module RungeKuttaPNEquations
