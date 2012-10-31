
! Bulirsch-Stoer Integrator

MODULE BSINT
	USE PHYSICS

CONTAINS


	SUBROUTINE RK4(posA, velA, massA, posB, velB, massB, posC, velC, massC, dt)
		IMPLICIT NONE
		DOUBLE PRECISION, DIMENSION(3) :: posA, velA, posB, velB, posC, velC
		DOUBLE PRECISION :: massA, massB, massC, dt
		DOUBLE PRECISION, DIMENSION(3) :: k0, j0, k1, j1, k2, j2, k3, j3, accel

		k0(:) = velC(:)*dt
		CALL Gravity(posA, massA, posB, massB, posC(:), massC, accel)
		j0(:) = accel(:)*dt
		
		k1(:) = (velC(:)+j0(:)/2D0)*dt
		CALL Gravity(posA, massA, posB, massB, posC(:)+k0(:)/2D0, massC, accel)
		j1(:) = accel(:)*dt
		
		k2(:) = (velC(:)+j1(:)/2D0)*dt
		CALL Gravity(posA, massA, posB, massB, posC(:)+k1(:)/2D0, massC, accel)
		j2(:) = accel(:)*dt
		
		k3(:) = (velC(:)+j2(:))*dt
		CALL Gravity(posA, massA, posB, massB, posC(:)+k2(:), massC, accel)
		j3(:) = accel(:)*dt
		
		posC(:) = posC(:) + (k0(:) + 2D0*k1(:) + 2D0*k2(:) + k3(:))/6D0
		velC(:) = velC(:) + (j0(:) + 2D0*j1(:) + 2D0*j2(:) + j3(:))/6D0

	END SUBROUTINE RK4

	SUBROUTINE Euler2(posA, velA, massA, posB, velB, massB, posC, velC, massC, dt)
		IMPLICIT NONE
		INTEGER :: ii
		DOUBLE PRECISION, DIMENSION(3) :: posA, velA, posB, velB, posC, velC
		DOUBLE PRECISION, DIMENSION(3) :: accel, dist1, dist2
		DOUBLE PRECISION :: grav1, grav2, massA, massB, massC, dt
		
		CALL Gravity(posA, massA, posB, massB, posC, massC, accel)

		! integrate velocity
		velC(1) = velC(1) + accel(1)*dt
		velC(2) = velC(2) + accel(2)*dt
		velC(3) = velC(3) + accel(3)*dt

		! integrate position
		posC(1) = posC(1) + velC(1)*dt
		posC(2) = posC(2) + velC(2)*dt
		posC(3) = posC(3) + velC(3)*dt

	END SUBROUTINE Euler2

	SUBROUTINE BS(posA, velA, massA, posB, velB, massB, posC, velC, massC, dt, steps)
		INTEGER :: n, i, j, iter, steps
		PARAMETER (iter = 4)
		DOUBLE PRECISION, DIMENSION(3) :: posA, velA, posB, velB, posC, velC
		DOUBLE PRECISION :: massA, massB, massC, dt
		DOUBLE PRECISION, DIMENSION(3) :: vel_n, vel_nn, vel_old, accel, factor, ret
		DOUBLE PRECISION, DIMENSION(3) :: pos_n, pos_nn, pos_old
		DOUBLE PRECISION :: h
		DOUBLE PRECISION, DIMENSION(3,iter) :: res, res2, vel_all, pos_all, ans
		! START  integration at different stepsizes
		DO n=1,iter
			pos_n(:) = posC(:)
			pos_nn(:) = posC(:)
			pos_old(:) = posC(:)
			vel_n(:) = velC(:)
			vel_nn(:) = velC(:)
			vel_old(:) = velC(:)
			h = 0.5D0*dt/(n*2D0)
			! TODO: co-evolve positions for stability
			CALL Gravity(posA, massA, posB, massB, pos_old, massC, accel)
			vel_nn(:) = vel_old(:) + h*accel(:)
			pos_nn(:) = pos_old(:) + h*vel_old(:)
			CALL Gravity(posA, massA, posB, massB, pos_nn, massC, accel)
			DO i=1,2*n-1
				vel_n(:) = vel_n(:) + 2D0*h*accel(:)
				pos_n(:) = pos_n(:) + 2D8*h*vel_nn(:)
				CALL Gravity(posA, massA, posB, massB, pos_n, massC, accel)
				vel_nn(:) = vel_nn(:) + 2D0*h*accel(:)
				pos_nn(:) = pos_nn(:) + 2D0*h*vel_n(:)
				CALL Gravity(posA, massA, posB, massB, pos_nn, massC, accel)
			ENDDO
			vel_n(:) = vel_n(:) + 2D0*h*accel(:)
			pos_n(:) = pos_n(:) + 2D0*h*vel_nn(:)
			CALL Gravity(posA, massA, posB, massB, pos_n, massC, accel)
			vel_nn(:) = vel_nn(:) + h*accel(:)
			pos_nn(:) = pos_nn(:) + h*vel_n(:)
			res2(:,n) = (vel_n(:)+vel_nn(:))*0.5D0
			res(:,n) = (pos_n(:)+pos_nn(:))*0.5D0
		ENDDO

		! START  rational extrapolation
		DO i=1,iter
			pos_all(:,i) = posC(:)
			vel_all(:,i) = velC(:)
			ans(:,i) = 0D0
		ENDDO

		DO n=1,iter
			DO i=1,iter-n
				factor(:) = 1D0 - (res(:,i+1)-res(:,i))/(res(:,i+1)-pos_all(:,i+1))
				! check
				DO j=1,3
					IF (ABS(res(j,i+1)-pos_all(j,i+1)) .LT. 1D-9) THEN
						factor(j) = 1D0
					ENDIF
				ENDDO
				factor(:) = factor(:)*((i+n+1D0)/(i+1D0))**2D0 - 1D0
				ans(:,i) = res(:,i+1) + (res(:,i+1)-res(:,i))/factor(:)
			ENDDO
			DO i=1,iter
				pos_all(:,i) = res(:,i)
				!res(:,i) = ans(:,i)
			ENDDO
		ENDDO
		velC(:) = res2(:,iter)
		posC(:) = res(:,iter)
	END SUBROUTINE BS

END MODULE BSINT
