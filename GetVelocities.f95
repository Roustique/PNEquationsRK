function F(v, phi, c, e)
use my_prec
implicit none

	real(mp), dimension(2) :: F, v
	real(mp) :: phi, c, e

	F(1) = v(1)**4 - v(2)**4 + 4 * (c**2 / 3 + phi * (1 + e)) * v(1)**2 + &
		4 * (c**2 / 3 - phi * (1 - e)) * v(2)**2 + 16 / 3 * e * phi * (phi - c**2)
	F(2) = 1 / (1 + e) * v(1)**3 - 1 / (1 - e) * v(2)**3 + &
		(3 * phi + c**2 / (1 + e)) * v(1) - (3 * phi + c**2 / (1 - e)) * v(2)

end function F


function J(v, phi, c, e)
use my_prec
implicit none

	real(mp), dimension(2,2) :: J
	real(mp), dimension(2) :: v
	real(mp) :: phi, c, e

	J(1,1) =   4 * v(1)**3 + 8 * (c**2 / 3 + phi * (1 + e)) * v(1)
	J(1,2) = - 4 * v(2)**3 + 8 * (c**2 / 3 - phi * (1 - e)) * v(2)
	J(2,1) =   3 / (1 + e) * v(1)**2 + (3 * phi + c**2 / (1 + e))
	J(2,2) = - 3 / (1 - e) * v(2)**2 - (3 * phi + c**2 / (1 - e))

end function J


Program getvel
use my_prec
implicit none

	interface

		function F(v, phi, c, e)
		use my_prec

			real(mp), dimension(2) :: F, v
			real(mp) :: phi, c, e

		end function F


		function J(v, phi, c, e)
		use my_prec

			real(mp), dimension(2,2) :: J
			real(mp), dimension(2) :: v
			real(mp) :: phi, c, e

		end function J

	end interface

	real(mp), dimension(2,2) :: A
	real(mp), dimension(2) :: B, v, vold
	real(mp) :: phi, e, Det, Det1, Det2, M, majsemiax
	real(mp), parameter :: c = 1.0
	real(mp), parameter :: eps = 1e-15

	open(1, file = './input/initvel.dat')
	read(1,*) majsemiax, e, M
	close(1)

	phi = M * 1.47662e-1 / (majsemiax * (1 - e**2))

	v(1) = sqrt(phi) * (1 + e)
	v(2) = sqrt(phi) * (1 - e)

	write(*,*)"Newtonian velocities:"
	write(*,*) v
	write(*,*)"Iterations finding PN velocities:"

	do while (norm2(v - vold) > eps)

		A = J(v, phi, c, e)
		B = matmul(A, v) - F(v, phi, c, e)

		vold = v

		!Det  = A(1,1) * A(2,2) - A(2,1) * A(1,2)
		!Det1 =   B(1) * A(2,2) -   B(2) * A(1,2)
		!Det2 = A(1,1) *   B(2) - A(2,1) *   B(1)

		!v(1) = Det1 / Det
		!v(2) = Det2 / Det

		A(1,2) = A(1,2) / A(1,1)
		B(1) = B(1) / A(1,1)
		A(1,1) = 1

		A(2,1) = 0
		A(2,2) = A(2,2) - A(1,2) * A(2,1)
		B(2) = B(2) - B(1) * A(2,1)

		B(2) = B(2) / A(2,2)
		A(2,2) = 1

		v(2) = B(2)
		v(1) = B(1) + A(1,2) * B(2)

		write(*,*) v

	enddo

	open(2, file = './output/perivel.dat')
	write(2,*) v(1)
	close(2)

end program getvel
