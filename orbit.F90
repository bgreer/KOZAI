PROGRAM ORBIT
	USE omp_lib
	USE BSINT
	USE PHYSICS
	IMPLICIT NONE

	! misc variables
	INTEGER :: ii, ij, ik, id
	DOUBLE PRECISION :: ecc, inc

	! number of systems to evolve simultaneously
	INTEGER :: nsystems
	PARAMETER(nsystems=1)
	
	! body parameters
	! mass in solar mass
	! position in AU
	! velocity in AU/yr
	DOUBLE PRECISION, DIMENSION(nsystems) :: mass_a, mass_b, mass_c
	DOUBLE PRECISION, DIMENSION(nsystems) :: ecc_b, ecc_c, sma_b, sma_c, inc_b, inc_c
	DOUBLE PRECISION, DIMENSION(3,nsystems) :: pos_a, pos_b, pos_c
	DOUBLE PRECISION, DIMENSION(3,nsystems) :: vel_a, vel_b, vel_c

	! integration parameters
	! time in years
	DOUBLE PRECISION, DIMENSION(nsystems) :: currtime, dt, endtime
	! 0 = ok, 1 = escaped, 2 = collided
	INTEGER, DIMENSION(nsystems) :: stat

	WRITE(*,'(A)') "Run Parameters:"
	!$OMP PARALLEL
	id = OMP_GET_THREAD_NUM()
	IF (id .EQ. 0) THEN
		WRITE(*,'(A,I4)') "  Number of threads = ",OMP_GET_NUM_THREADS()
	ENDIF
	!$OMP END PARALLEL
	WRITE(*,'(A,I7,A)') "  ",nsystems," Systems"

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! initialize a bunch of systems based on orbit parameters
	DO ii=1,nsystems
		mass_a(ii) = 1D0 ! sun
		
		mass_b(ii) = 1D-3 ! driver
		ecc_b(ii) = 0.2D0 ! eccentricity
		sma_b(ii) = 10D0 ! semimajor axis
		inc_b(ii) = 2.0D0 ! inclination in radians?

		mass_c(ii) = 1D-6 ! test particle
		ecc_c(ii) = 0.5D0
		sma_c(ii) = 1D0
		inc_c(ii) = 0D0
	ENDDO
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! reparameterize in terms of cartesian pos/vel
	DO ii=1,nsystems
		! sun does not move
		pos_a(1,ii) = 0D0
		pos_a(2,ii) = 0D0
		pos_a(3,ii) = 0D0
		vel_a(1,ii) = 0D0
		vel_a(2,ii) = 0D0
		vel_a(3,ii) = 0D0
		
		! set orbit for test particle c
		CALL SetOrbit(mass_a(ii),mass_c(ii),sma_c(ii),ecc_c(ii),inc_c(ii),pos_c,vel_c)

		currtime(ii) = 0D0
		dt(ii) = 0.0003D0
		endtime(ii) = 1D6
		stat(ii) = 0
	ENDDO

	WRITE(*,'(A)') "Systems Initialized, Beginning Integration."

	! evolve the systems in parallel
	!$OMP PARALLEL DO PRIVATE(ii,ij) SCHEDULE(DYNAMIC,1)
	DO ii=1,nsystems
		ij = 0
		DO WHILE (currtime(ii) .LT. endtime(ii) .AND. stat(ii) .EQ. 0)

			! decide timestep?

			CALL DriverPos(pos_b(:,ii),currtime(ii),ecc_b(ii),sma_b(ii),inc_b(ii),mass_a(ii),mass_b(ii))
			! integrate
			CALL RK4(pos_a(:,ii), vel_a(:,ii), mass_a(ii), &
					pos_b(:,ii), vel_b(:,ii), mass_b(ii), &
					pos_c(:,ii), vel_c(:,ii), mass_c(ii), dt(ii))

			! check for escape
			IF (Escaping(pos_c(:,ii), vel_c(:,ii), mass_c(ii), &
					pos_a(:,ii), vel_a(:,ii), mass_a(ii))) THEN
!				WRITE(*,'(A,I6,A)') "--System ",ii," Escaped"
!				stat(ii) = 1
			ENDIF

			! check for collision
			IF (HasCollided(pos_c(:,ii), pos_a(:,ii))) THEN
				stat(ii) = 2
			ENDIF

			IF (ii .EQ. 1 .AND. ij .EQ. 10000) THEN
				CALL Eccentricity(pos_a(:,ii), vel_a(:,ii), mass_a(ii), &
					pos_c(:,ii), vel_c(:,ii), mass_c(ii), ecc)
				CALL Inclination(pos_c(:,ii), vel_c(:,ii), inc)
				!WRITE(*,'(E15.6,7E15.6)') currtime(ii), pos_b(:,ii), pos_c(:,ii), ecc
				WRITE(*,'(6E15.6)') currtime(ii), ecc, 180D0-inc*180D0/3.14159D0, pos_c(:,ii)
				ij = 0
			ENDIF
			currtime(ii) = currtime(ii) + dt(ii)
			ij = ij + 1
		ENDDO
	ENDDO
	!$OMP END PARALLEL DO

	WRITE(*,'(A)') "Integration Complete"

	! output stats
	

END PROGRAM

