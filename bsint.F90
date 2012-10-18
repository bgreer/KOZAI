
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

	! testing integrator
	SUBROUTINE Euler(posA, velA, massA, posB, velB, massB, dt)
		IMPLICIT NONE
		INTEGER :: ii
		DOUBLE PRECISION, DIMENSION(3) :: posA, velA, posB, velB
		DOUBLE PRECISION, DIMENSION(3) :: accel, dist
		DOUBLE PRECISION :: grav, massA, massB, dt

		dist(1) = posB(1) - posA(1)
		dist(2) = posB(2) - posA(2)
		dist(3) = posB(3) - posA(3)

		! solve for gravity
		grav = -G*massA/(VecMag(dist))**3D0
		accel(1) = grav*dist(1)
		accel(2) = grav*dist(2)
		accel(3) = grav*dist(3)

		! integrate velocity
		velB(1) = velB(1) + accel(1)*dt
		velB(2) = velB(2) + accel(2)*dt
		velB(3) = velB(3) + accel(3)*dt

		! integrate position
		posB(1) = posB(1) + velB(1)*dt
		posB(2) = posB(2) + velB(2)*dt
		posB(3) = posB(3) + velB(3)*dt

	END SUBROUTINE Euler

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


	SUBROUTINE BS(steps)
		INTEGER :: n, i, iter, steps
		PARAMETER (iter = 3)
		DOUBLE PRECISION :: yn, xn, ynn, xnn, x_old, y_old, h
		DOUBLE PRECISION :: factor, ret
		DOUBLE PRECISION, DIMENSION(iter) :: res, yold, ans
		! START  integration at different stepsizes
		DO n=1,iter
			yn = 0D0
			xn = 0D0
			ynn = 0D0
			xnn = 0D0
			y_old = 0D0
			x_old = 0D0
			h = 0.5D0/(n*steps*1D0)
			yn = y_old
			xn = x_old
			ynn = y_old + h!*f(x_old,y_old)
			xnn = x_old + h
			DO i=1,2*n*steps-1,2
				yn = yn + 2.0*h!*f(xnn,ynn)
				xn = xn + 2.0*h
				ynn = ynn + 2.0*h!*f(xn,yn)
				xnn = xnn + 2.0*h
			ENDDO
			ynn = ynn - h!*f(xn,yn)
			res(n) = (yn+ynn)*0.5
		ENDDO
		

		! START  rational extrapolation
!	for (i=0; i<iter; i++) 
!	{
!		yold[i] = 0.0;
!		ans[i] = 0.0;
!	}
!
!	for (n=1; n<iter; n++)
!	{
!		for (i=0; i<iter-n; i++)
!		{
!			factor = 1.0 - (res[i+1]-res[i])/(res[i+1]-yold[i+1]);
!			factor = factor*pow((i+n+1.)/(i+1.), 2.0) - 1.0;
!			ans[i] = res[i+1] + (res[i+1]-res[i])/factor;
!		}
!		for (i=0; i<iter; i++)
!		{
!			yold[i] = res[i];
!			res[i] = ans[i];
!		}
!	}
!	ret = ans[0];
!	free(res);
!	free(yold);
!	free(ans);
!	return ret;

	END SUBROUTINE BS

END MODULE BSINT
