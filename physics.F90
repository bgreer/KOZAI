

MODULE PHYSICS

	! G in AU, years, solar masses
	DOUBLE PRECISION :: G = 39.4812494562
	DOUBLE PRECISION :: PI = 3.14159265359

CONTAINS

	SUBROUTINE Gravity(posA, massA, posB, massB, posC, massC, accel)
		IMPLICIT NONE
		DOUBLE PRECISION, DIMENSION(3) :: posA, velA, posB, velB, posC, velC
		DOUBLE PRECISION, DIMENSION(3) :: accel, dist1, dist2, dist3
		DOUBLE PRECISION :: massA, massB, massC, gravA, gravB, gravC

		dist1(:) = posC(:) - posA(:)
		dist2(:) = posC(:) - posB(:)
		dist3(:) = posB(:) - posA(:)

		! solve for gravity
		gravA = -G*massA/(VecMag(dist1))**3D0
		gravB = -G*massB/(VecMag(dist2))**3D0
		gravC = -G*massB/(VecMag(dist3))**3D0

		accel(:) = gravA*dist1(:) + gravB*dist2(:) + gravC*dist3(:)
	END SUBROUTINE Gravity

	! measure current inclination relative to.. uh.. z-axis?
	SUBROUTINE Inclination(posB, velB, i)
		IMPLICIT NONE
		DOUBLE PRECISION, DIMENSION(3) :: posB, velB, angmom
		DOUBLE PRECISION :: i

		angmom(1) = posB(2)*velB(3) - posB(3)*velB(2)
		angmom(2) = posB(3)*velB(1) - posB(1)*velB(3)
		angmom(3) = posB(1)*velB(2) - posB(2)*velB(1)
		! A dot B = ABcos(theta)
		i = ACOS(-angmom(3)/VecMag(angmom))
	END SUBROUTINE Inclination

	! new because im dumb
	SUBROUTINE Eccentricity(posA, velA, massA, posB, velB, massB, ecc)
		IMPLICIT NONE
		DOUBLE PRECISION, DIMENSION(3) :: posA, velA, posB, velB, angmom, dist, relvel, eccvec
		DOUBLE PRECISION :: ecc, spm, en, h, mu, massA, massB, r

		dist(1) = posB(1) - posA(1)
		dist(2) = posB(2) - posA(2)
		dist(3) = posB(3) - posA(3)
		relvel(1) = velB(1) - velA(1)
		relvel(2) = velB(2) - velA(2)
		relvel(3) = velB(3) - velA(3)
		r = VecMag(dist)

		! reduced mass
		spm = massA*massB/(massA+massB)
		! gravitational parameter
		mu = G*(massA+massB)
		! specific angular momentum
		angmom(1) = dist(2)*relvel(3) - dist(3)*relvel(2)
		angmom(2) = dist(3)*relvel(1) - dist(1)*relvel(3)
		angmom(3) = dist(1)*relvel(2) - dist(2)*relvel(1)
		angmom(:) = angmom(:)

		eccvec(1) = relvel(2)*angmom(3) - relvel(3)*angmom(2)
		eccvec(2) = relvel(3)*angmom(1) - relvel(1)*angmom(3)
		eccvec(3) = relvel(1)*angmom(2) - relvel(2)*angmom(1)

		eccvec(:) = eccvec(:) / mu

		eccvec(1) = eccvec(1) - dist(1)/r
		eccvec(2) = eccvec(2) - dist(2)/r
		eccvec(3) = eccvec(3) - dist(3)/r

		ecc = VecMag(eccvec)
	END SUBROUTINE Eccentricity

	! analytic solution for driver position
	SUBROUTINE DriverPos(pos, time, ecc, sma, inc, massA, massB)
		IMPLICIT NONE
		DOUBLE PRECISION, DIMENSION(3) :: pos
		DOUBLE PRECISION :: time, x, ecc, sma, inc, period, massA, massB

		period = SQRT(4D0*PI*PI*sma*sma*sma/(G*(massA+massB)))
		x = 2D0*PI*time/period

		pos(1) = sma*COS(x)*COS(inc)
		pos(2) = -sma*SIN(x)
		pos(3) = sma*COS(x)*SIN(inc)

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
		IF (VecMag(dist) .LT. 1d-2) THEN
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

	! Initialize an orbit based on orbital parameters
	! Stolen from Sam, but I'll go ahead and take credit.
	SUBROUTINE SetOrbit(massA,massB,a,e,i,pos,vel)
		IMPLICIT NONE
		DOUBLE PRECISION :: massA, massB, a, p, e, i, o, w, v, GMU
		DOUBLE PRECISION, DIMENSION(3) :: pos, vel
		DOUBLE PRECISION, DIMENSION(3,1) :: rpqw
		DOUBLE PRECISION, DIMENSION(3,1) ::  vpqw 
		DOUBLE PRECISION, DIMENSION(3,3) :: pqw2xyz ! Direction Cosine Matrix
		DOUBLE PRECISION, DIMENSION(3,1) :: rxyz !radius in inertial frame
		DOUBLE PRECISION, DIMENSION(3,1) :: vxyz !velocity in inertial frame

		! maybe these will be needed later?
		o = 0D0 ! right ascension of something
		w = 0D0 ! argument of periapsis
		v = 0D0 ! true anomaly

		GMU = G*(massA+massB)
		p = a*(1D0-e*e)

		!Finds Rp, Rq, Rw of orbit frame
		rpqw(1,1)=(p*cos(v))/(1+e*cos(v)) 
		rpqw(2,1)=(p*sin(v))/(1+e*cos(v))
		rpqw(3,1)=0D0

		!Finds Vp, Vq, Vw of orbit frame
		vpqw(1,1)=-(sqrt(GMU/p))*sin(v)
		vpqw(2,1)=(sqrt(GMU/p))*(e+cos(v))
		vpqw(3,1)=0D0

		!Finds direction cosine matrix knowing 3D orbit values o,w,i
		pqw2xyz(1,1)=cos(o)*cos(w)-sin(o)*sin(w)*cos(i)
		pqw2xyz(1,2)=-cos(o)*sin(w)-sin(o)*cos(w)*cos(i)
		pqw2xyz(1,3)=sin(o)*sin(i)
		pqw2xyz(2,1)=sin(o)*cos(w)+cos(o)*sin(w)*cos(i)
		pqw2xyz(2,2)=-sin(o)*sin(w)+cos(o)*cos(w)*cos(i)
		pqw2xyz(2,3)=-cos(o)*sin(i)
		pqw2xyz(3,1)=sin(w)*sin(i)
		pqw2xyz(3,2)=cos(w)*sin(i)
		pqw2xyz(3,3)=cos(i)

		!Goes from Orbit frame to inertial frame
		rxyz=matmul(pqw2xyz,rpqw)
		vxyz=matmul(pqw2xyz,vpqw)

		! Output
		pos(:)=rxyz(:,1)
		vel(:)=vxyz(:,1)
		
	END SUBROUTINE SetOrbit

END MODULE PHYSICS
