Program main
use my_prec
use RungeKuttaPNEquations
implicit none
real(mp), dimension(:,:), allocatable :: Coord, Velocity
real(mp), dimension(:), allocatable :: Energy, Angmom
real(mp) :: step, Mass, E0
integer(4) :: amount, i
real(mp), dimension(2) :: Coord0, Velocity0
character(3) :: Key

call getarg(1,Key)

if ((Key /= "-NT").and.(Key /= "-GS").and.(Key /= "-GH").and.(Key /= "-F0").and.(Key /= "-FF")) then
    write(*,*)"The key is missing or incorrect"
    write(*,*)"Write down one of the following keys:"
    write(*,*)"-NT for Newtonian EoM"
    write(*,*)"-GS for GR EoM in Schwarzschild coordinates"
    write(*,*)"-GH for GR EoM in harmonic coordinates"
    write(*,*)"-F0 for FGT EoM with K=0"
    write(*,*)"-FF for FGT EoM with K=1/2"
    write(*,*)"or anything to quit"
    read(*,*)Key
    if ((Key /= "-NT").and.(Key /= "-GS").and.(Key /= "-GH").and.(Key /= "-F0").and.(Key /= "-FF")) then
        call exit(1)
    end if
end if

open(1,file="./input/data.dat")

read(1,*)step, amount
read(1,*)Coord0(1), Coord0(2), Velocity0(1), Velocity0(2)
read(1,*)Mass

close(1)

Coord0=Coord0
Velocity0=Velocity0

allocate(Coord(0:amount,2),Velocity(0:amount,2),Energy(0:amount),Angmom(0:amount))

call RKIntegrate(step, amount, Coord0, Velocity0, Mass, Coord, Velocity, Energy, Angmom, Key)

open(2,file="./output/result" // Key // ".dat")

write(2,*)"x, y, vx, vy, E, L"
do i=0,amount
    write(2,*)Coord(i,1),", ", Coord(i,2),", ", Velocity(i,1),", ",&
    Velocity(i,2), ", ",Energy(i), ", ", Angmom(i)
enddo

close(2)

end

