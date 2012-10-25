
! Bulirsch-Stoer Integrator

MODULE BSINT
	USE PHYSICS

CONTAINS

	SUBROUTINE Gravity(posA, massA, posB, massB, posC, massC, accel)
		IMPLICIT NONE
		DOUBLE PRECISION, DIMENSION(3) :: posA, velA, posB, velB, posC, velC
		DOUBLE PRECISION, DIMENSION(3) :: accel, dist1, dist2
		DOUBLE PRECISION :: massA, massB, massC, gravA, gravB

		dist1(1) = posC(1) - posA(1)
		dist1(2) = posC(2) - posA(2)
		dist1(3) = posC(3) - posA(3)
		dist2(1) = posC(1) - posB(1)
		dist2(2) = posC(2) - posB(2)
		dist2(3) = posC(3) - posB(3)

		! solve for gravity
		gravA = -G*massA/(VecMag(dist1))**3D0
		gravB = -G*massB/(VecMag(dist2))**3D0

		accel(1) = gravA*dist1(1) + gravB*dist2(1)
		accel(2) = gravA*dist1(2) + gravB*dist2(2)
		accel(3) = gravA*dist1(3) + gravB*dist2(3)
	END SUBROUTINE Gravity

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


	SUBROUTINE BS(posA, velA, massA, posB, velB, massB, posC, velC, massC, dt, steps, ret)
		INTEGER :: n, i, iter, steps
		PARAMETER (iter = 3)
		DOUBLE PRECISION, DIMENSION(3) :: posA, velA, posB, velB, posC, velC
		DOUBLE PRECISION :: massA, massB, massC, dt
		DOUBLE PRECISION, DIMENSION(3) :: pos_n, pos_nn, pos_old, accel, factor, ret
		DOUBLE PRECISION :: h
		DOUBLE PRECISION, DIMENSION(3,iter) :: res, pos_all, ans
		! START  integration at different stepsizes
		DO n=1,iter
			pos_n(:) = 0D0 ! CHANGE
			pos_nn(:) = 0D0
			pos_old(:) = 0D0
			h = 0.5D0/(n*steps*1D0)
			! TODO: co-evolve positions for stability
			CALL Gravity(posA, massA, posB, massB, pos_old, massC, accel)
			pos_nn(:) = pos_old(:) + h*accel(:)
			
			DO i=1,2*n*steps-1,2
				CALL Gravity(posA, massA, posB, massB, pos_nn, massC, accel)
				pos_n(:) = pos_n(:) + 2D0*h*accel(:)
				CALL Gravity(posA, massA, posB, massB, pos_n, massC, accel)
				pos_nn(:) = pos_nn(:) + 2D0*h*accel(:)
			ENDDO
			CALL Gravity(posA, massA, posB, massB, pos_n, massC, accel)
			pos_nn(:) = pos_nn(:) - h*accel(:)
			res(:,n) = (pos_n(:)+pos_nn(:))*0.5D0
		ENDDO
		

		! START  rational extrapolation
		DO i=1,iter
			pos_all(:,i) = 0D0
			ans(:,i) = 0D0
		ENDDO

		DO n=1,iter
			DO i=1,iter-n
				factor(:) = 1D0 - (res(:,i+1)-res(:,i))/(res(:,i+1)-pos_all(:,i+1))
				factor(:) = factor(:)*((i+n+1D0)/(i+1D0))**2D0 - 1D0
				ans(:,i) = res(:,i+1) + (res(:,i+1)-res(:,i))/factor(:)
			ENDDO
			DO i=1,iter
				pos_all(:,i) = res(:,i)
				res(:,i) = ans(:,i)
			ENDDO
		ENDDO
		ret(:) = ans(:,1)
	END SUBROUTINE BS

END MODULE BSINT
