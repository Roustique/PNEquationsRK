Program main
use my_prec
use RungeKuttaPNEquations
implicit none
real(mp), dimension(:,:), allocatable :: Coord, Velocity
real(mp) :: step, GM, E0
integer(4) :: amount, i
real(mp), dimension(2) :: Coord0, Velocity0

read(*,*)step, amount
read(*,*)Coord0(1), Coord0(2), Velocity0(1), Velocity0(2)
read(*,*)GM

allocate(Coord(0:amount,2),Velocity(0:amount,2))

call RKIntegrate(step, amount, Coord0, Velocity0, GM, Coord, Velocity)

E0=(Velocity(0,1)**2+Velocity(0,2)**2)/2-GM/sqrt(Coord(0,1)**2+Coord(0,2)**2)

write(*,*)"x, y, vx, vy, E"
do i=0,amount
   write(*,*)Coord(i,1)*2.95e+11,", ", Coord(i,2)*2.95e+11,", ", Velocity(i,1)*2.95e+11/3.156e+7,", ",&
   Velocity(i,2)*2.95e+11/3.156e+7, ", ",((Velocity(i,1)**2+Velocity(i,2)**2)/2-GM/sqrt(Coord(i,1)**2+Coord(i,2)**2))/E0   
enddo

end

