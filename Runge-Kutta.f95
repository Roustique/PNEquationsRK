module RungeKuttaPNEquations
use my_prec
implicit none
contains

subroutine EqMot(x,v,GM,c,k)
real(mp), intent(in), dimension(2) :: x, v
real(mp), intent(in) :: GM, c
real(mp), intent(out), dimension(4) :: k
real(mp) :: C1, C2, r

    r=sqrt(x(1)**2+x(2)**2)
    
    C1=-(c*c+v(1)*v(1)+v(2)*v(2)-4.0_mp*GM/r)
    C2=4*(v(1)*x(1)+v(2)*x(2))

    k(1)=v(1)
    k(2)=v(2)
    k(3)=(C1*x(1)+C2*v(1))*GM/(c**2.0_mp*r**3.0_mp)
    k(4)=(C1*x(2)+C2*v(2))*GM/(c**2.0_mp*r**3.0_mp)

end subroutine EqMot

subroutine RKIntegrate(step, amount, Coord0, Velocity0, GM, Coord, Velocity, Energy)
real(mp), intent(in) :: step, GM
integer(4), intent(in) :: amount
real(mp), intent(in), dimension(2) :: Coord0, Velocity0   !Initial Values 
real(mp), intent(out), dimension(0:amount,2) :: Coord, Velocity
real(mp), intent(out), dimension(0:amount) :: Energy
real(mp), parameter :: c=3.20727117779e+6 !Speed of light 
real(mp) :: r  !r is module of Coord vector
real(mp), dimension(4,4) :: k  !RK4 coefficients
integer(4) :: i

    Coord(0,:)=Coord0(:)
    Velocity(0,:)=Velocity0(:)
    r=sqrt(Coord0(1)**2+Coord0(2)**2)
    Energy(0)=sum(Velocity(0,:)**2)/2.0_mp-GM/r!+(GM/r)**2/(2*c**2)
    
    do i=1,amount
       
        call EqMot(Coord(i-1,:),Velocity(i-1,:),GM,c,k(1,:))
        call EqMot(Coord(i-1,:)+step*k(1,:2)/2.0_mp,Velocity(i-1,:)+step*k(1,3:)/2.0_mp,GM,c,k(2,:))
        call EqMot(Coord(i-1,:)+step*k(2,:2)/2.0_mp,Velocity(i-1,:)+step*k(2,3:)/2.0_mp,GM,c,k(3,:))
        call EqMot(Coord(i-1,:)+step*k(3,:2),Velocity(i-1,:)+step*k(3,3:),GM,c,k(4,:))

		Coord(i,:)=Coord(i-1,:)+step/6*(k(1,:2)+2*k(2,:2)+2*k(3,:2)+k(4,:2))
        Velocity(i,:)=Velocity(i-1,:)+step/6*(k(1,3:)+2*k(2,3:)+2*k(3,3:)+k(4,3:))
 
        r=sqrt(Coord(i,1)**2+Coord(i,2)**2)

        Energy(i)=sum(Velocity(i,:)**2)/2.0_mp-GM/r!+(GM/r)**2/(2*c**2)

    enddo

end subroutine RKIntegrate

!subroutine RKIntegrate2()

!end subroutine RKIntegrate2()

end module RungeKuttaPNEquations
