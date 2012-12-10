PROGRAM ORBIT
	USE omp_lib
	USE ParseInput
	USE BSINT
	USE PHYSICS
	IMPLICIT NONE

	! misc variables
	INTEGER :: ii, ij, ik, id, argcount
	DOUBLE PRECISION :: ecc, inc
	REAL :: rand
	INTEGER :: values(1:8), k
	INTEGER, DIMENSION(:), ALLOCATABLE :: seed
	
	! body parameters
	! mass in solar mass
	! position in AU
	! velocity in AU/yr
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mass_a, mass_b, mass_c
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ecc_b, ecc_c, sma_b, sma_c, inc_b, inc_c
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: pos_a, pos_b, pos_c
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: vel_a, vel_b, vel_c

	! used for determining kozai resonance
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ecc_max, ecc_min, inc_max, inc_min

	! integration parameters
	! time in years
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: currtime, dt
	! 0 = ok, 1 = escaped, 2 = collided
	INTEGER, DIMENSION(:), ALLOCATABLE :: stat

	! Parse command-line args
	CALL SetDefaults()
	CALL ReadCommandLine()
	CALL PrintDetails()

	! Init RNG on current time
	CALL DATE_AND_TIME(values = values)
	CALL RANDOM_SEED(size=k)
	ALLOCATE(seed(1:k))
	seed(:) = values(8)
	CALL RANDOM_SEED(put=seed)
	DEALLOCATE(seed)

	! Allocate space for systems
	ALLOCATE(mass_a(nsystems))
	ALLOCATE(mass_b(nsystems))
	ALLOCATE(mass_c(nsystems))
	ALLOCATE(ecc_b(nsystems))
	ALLOCATE(ecc_c(nsystems))
	ALLOCATE(sma_b(nsystems))
	ALLOCATE(sma_c(nsystems))
	ALLOCATE(inc_b(nsystems))
	ALLOCATE(inc_c(nsystems))
	ALLOCATE(pos_a(3,nsystems))
	ALLOCATE(pos_b(3,nsystems))
	ALLOCATE(pos_c(3,nsystems))
	ALLOCATE(vel_a(3,nsystems))
	ALLOCATE(vel_b(3,nsystems))
	ALLOCATE(vel_c(3,nsystems))
	ALLOCATE(currtime(nsystems))
	ALLOCATE(dt(nsystems))
	ALLOCATE(stat(nsystems))
	ALLOCATE(ecc_max(2,nsystems))
	ALLOCATE(ecc_min(2,nsystems))
	ALLOCATE(inc_max(2,nsystems))
	ALLOCATE(inc_min(2,nsystems))
	ecc_max(1,:) = 0D0
	ecc_max(2,:) = -1D0
	ecc_min(1,:) = 1D0
	ecc_min(2,:) = 0D0
	inc_max = 0D0
	inc_min = 0D0

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! initialize a bunch of systems based on orbit parameters
	DO ii=1,nsystems
		mass_a(ii) = 1D0 ! sun

		CALL RANDOM_NUMBER(rand)
		mass_b(ii) = rand*(driver_mass_max-driver_mass_min)+driver_mass_min ! driver
		CALL RANDOM_NUMBER(rand)
		ecc_b(ii) = rand*(driver_ecc_max-driver_ecc_min)+driver_ecc_min
		CALL RANDOM_NUMBER(rand)
		sma_b(ii) = rand*(driver_sma_max-driver_sma_min)+driver_sma_min
		CALL RANDOM_NUMBER(rand)
		inc_b(ii) = rand*(driver_inc_max-driver_inc_min)+driver_inc_min

		CALL RANDOM_NUMBER(rand)
		mass_c(ii) = rand*(test_mass_max-test_mass_min)+test_mass_min ! test particle
		CALL RANDOM_NUMBER(rand)
		ecc_c(ii) = rand*(test_ecc_max-test_ecc_min)+test_ecc_min
		CALL RANDOM_NUMBER(rand)
		sma_c(ii) = rand*(test_sma_max-test_sma_min)+test_sma_min
		CALL RANDOM_NUMBER(rand)
		inc_c(ii) = rand*(test_inc_max-test_inc_min)+test_inc_min
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
		CALL SetOrbit(mass_a(ii),mass_c(ii),sma_c(ii),ecc_c(ii),inc_c(ii),pos_c(:,ii),vel_c(:,ii))

		currtime(ii) = 0D0
		dt(ii) = 0.0003D0
		stat(ii) = 0
	ENDDO
	WRITE(*,'(A)') "Systems Initialized, Beginning Integration."

	! evolve the systems in parallel
	!$OMP PARALLEL DO PRIVATE(ii,ij,ecc,inc) SCHEDULE(DYNAMIC,1)
	DO ii=1,nsystems
		ij = 0
		DO WHILE (currtime(ii) .LT. endtime .AND. stat(ii) .EQ. 0)

			! decide timestep?

			CALL DriverPos(pos_b(:,ii),currtime(ii),ecc_b(ii),sma_b(ii),inc_b(ii),mass_a(ii),mass_b(ii))
			! integrate
			CALL RK4(pos_a(:,ii), vel_a(:,ii), mass_a(ii), &
					pos_b(:,ii), vel_b(:,ii), mass_b(ii), &
					pos_c(:,ii), vel_c(:,ii), mass_c(ii), dt(ii))
			! check for escape
			IF (Escaping(pos_c(:,ii), vel_c(:,ii), mass_c(ii), &
					pos_a(:,ii), vel_a(:,ii), mass_a(ii))) THEN
				WRITE(*,'(A,I6,A)') "--System ",ii," Escaped"
				stat(ii) = 1
			ENDIF

			! check for collision
			IF (HasCollided(pos_c(:,ii), pos_a(:,ii))) THEN
				stat(ii) = 2
			ENDIF

			! compute current eccentricity and inclination
			IF (ij .EQ. 1000) THEN
				CALL Eccentricity(pos_a(:,ii), vel_a(:,ii), mass_a(ii), &
						pos_c(:,ii), vel_c(:,ii), mass_c(ii), ecc)
				CALL Inclination(pos_c(:,ii), vel_c(:,ii), inc)
				! Record ecc and/or inc if they are extrema
				IF (ecc .GT. ecc_max(1,ii)) THEN
					ecc_max(1,ii) = ecc
					ecc_max(2,ii) = currtime(ii)
				ENDIF
				IF (ecc .LT. ecc_min(1,ii)) THEN
					ecc_min(1,ii) = ecc
					ecc_min(2,ii) = currtime(ii)
				ENDIF

				IF (ii .EQ. 1 .AND. verbose) THEN
					WRITE(*,'(6E15.6)') currtime(ii), ecc, 180D0-inc*180D0/3.14159D0, vel_c(:,ii)
				ENDIF
				ij = 0
			ENDIF

			currtime(ii) = currtime(ii) + dt(ii)
			ij = ij + 1
		ENDDO
	ENDDO
	!$OMP END PARALLEL DO

	WRITE(*,'(A)') "Integration Complete"

	! output stats
	DO ii=1,nsystems
		! initial conditions
		WRITE(*,'(I5,I3,6E12.4$)') ii, stat(ii), mass_b(ii), inc_b(ii), ecc_b(ii), &
				mass_c(ii), inc_c(ii), ecc_c(ii)
		! inclination oscillation
		WRITE(*,'(4E12.4$)'), ecc_max(:,ii), ecc_min(:,ii)
		! eccentricity oscillation
		WRITE(*,'()')
	ENDDO

	! remember to deallocate, someday

END PROGRAM

