Program main
use my_prec
use RungeKuttaPNEquations
implicit none
real(mp), dimension(:,:), allocatable :: Coord, Velocity
real(mp), dimension(:), allocatable :: Energy
real(mp) :: step, GM, E0
integer(4) :: amount, i
real(mp), dimension(2) :: Coord0, Velocity0
character(2) :: Key

call getarg(1,Key)

open(1,file="data.dat")

 read(1,*)step, amount
 read(1,*)Coord0(1), Coord0(2), Velocity0(1), Velocity0(2)
 read(1,*)GM

close(1)

Coord0=Coord0
Velocity0=Velocity0

allocate(Coord(0:amount,2),Velocity(0:amount,2),Energy(0:amount))

call RKIntegrate(step, amount, Coord0, Velocity0, GM, Coord, Velocity, Energy, Key)

open(2,file="result.dat")

write(2,*)"x, y, vx, vy, E"
do i=0,amount
   write(2,*)Coord(i,1)*2.95e+5,", ", Coord(i,2)*2.95e+5,", ", Velocity(i,1)*2.95e+5,", ",&
   Velocity(i,2)*2.95e+5, ", ",Energy(i)
enddo

close(2)

end

