

MODULE PHYSICS

	! G in AU, years, solar masses
	DOUBLE PRECISION :: G = 39.4812494562
	DOUBLE PRECISION :: PI = 3.14159265359

CONTAINS

	! measure current eccentricity
	SUBROUTINE Eccentricity(posA, velA, massA, posB, velB, massB, ecc)
		IMPLICIT NONE
		DOUBLE PRECISION, DIMENSION(3) :: posA, velA, posB, velB, angmom, dist, relvel
		DOUBLE PRECISION :: ecc, spm, en, h, mu, massA, massB

		dist(1) = posB(1) - posA(1)
		dist(2) = posB(2) - posA(2)
		dist(3) = posB(3) - posA(3)

		relvel(1) = velB(1) - velA(1)
		relvel(2) = velB(2) - velA(2)
		relvel(3) = velB(3) - velA(3)

		! reduced mass
		spm = massA*massB/(massA+massB)
		! specific orbital energy
		en = (0.5D0 * massB * (VecMag(relvel))**2D0 - G * massB * massA / VecMag(dist))/spm
		! gravitational parameter
		mu = G*(massA+massB)
		! specific angular momentum
		angmom(1) = dist(2)*relvel(3) - dist(3)*relvel(2)
		angmom(2) = dist(3)*relvel(1) - dist(1)*relvel(3)
		angmom(3) = dist(1)*relvel(2) - dist(2)*relvel(1)
		h = VecMag(angmom)
		h = angmom(3)

		ecc = SQRT(1D0 + 2D0*en*h*h/(mu*mu))
		
	END SUBROUTINE Eccentricity

	! analytic solution for driver position
	SUBROUTINE DriverPos(pos, time, period)
		IMPLICIT NONE
		DOUBLE PRECISION, DIMENSION(3) :: pos
		DOUBLE PRECISION :: time, period, x

		x = 2D0*PI*time/period

		pos(1) = 5.2D0*COS(x)*COS(0.1)
		pos(2) = -5.2D0*SIN(x)
		pos(3) = 5.2D0*COS(x)*SIN(0.1)

	END SUBROUTINE DriverPos

	! determine if a bodyA is escaping another bodyB
	LOGICAL FUNCTION Escaping(posA, velA, massA, posB, velB, massB)
		DOUBLE PRECISION, DIMENSION(3) :: posA, velA, posB, velB, dist, relvel
		DOUBLE PRECISION :: massA, massB, KE, PE
		
		dist(1) = posB(1) - posA(1)
		dist(2) = posB(2) - posA(2)
		dist(3) = posB(3) - posA(3)

		relvel(1) = velB(1) - velA(1)
		relvel(2) = velB(2) - velA(2)
		relvel(3) = velB(3) - velA(3)

		KE = 0.5D0 * massA * (VecMag(relvel))**2D0
		PE = G * massB * massA / VecMag(dist)

		IF (KE .GT. PE) THEN
			Escaping = .TRUE.
		ELSE
			Escaping = .FALSE.
		ENDIF
	END FUNCTION Escaping

	! have two bodies come close enough to collide?
	LOGICAL FUNCTION HasCollided(posA, posB)
		DOUBLE PRECISION, DIMENSION(3) :: posA, posB, dist
		dist(1) = posB(1) - posA(1)
		dist(2) = posB(2) - posA(2)
		dist(3) = posB(3) - posA(3)
		IF (VecMag(dist) .LT. 1d-5) THEN
			HasCollided = .TRUE.
		ELSE
			HasCollided = .FALSE.
		ENDIF
	END FUNCTION HasCollided

	! magnitude of a 3-vector
	DOUBLE PRECISION FUNCTION VecMag(vector)
		DOUBLE PRECISION, DIMENSION(3) :: vector
		VecMag = SQRT(vector(1)*vector(1) + vector(2)*vector(2) + vector(3)*vector(3))
	END FUNCTION VecMag

END MODULE PHYSICS
