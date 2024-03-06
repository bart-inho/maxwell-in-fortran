module all

    use iso_fortran_env
    implicit none
    private
    public :: everything
    
    type everything

        ! Everything is a type that contains all the variables and subroutines needed for the simulation
        integer :: i, j, nx, sou_position, ice_position, rock_position, shield_position
        real(8) :: c, mu0 = 1.256e-6, eps0 = 8.854e-12, imp, sigma_max
        real(8) :: lambda_min = 12., dz, dt, t0, tau = 50.0, w0, pi = 4*TAN(1.0d0)
        real(8) :: ice_permittivity = 3.2, rock_permittivity = 5.0, air_permittivity = 1.0, shield_permittivity = 1.0e7
        real(8), allocatable :: epr(:), mur(:), Ex(:), Hz(:)

    contains

        ! Subroutines
        procedure, public :: EverythingInit, EverythingFini, WriteToText, GaussianSource, UpdateMagneticField, UpdateElectricField
        procedure, private :: BoundaryConditions

    end type everything
        
contains

    subroutine EverythingInit(this) ! Initialize everything
        implicit none
        class(everything), intent(inout) :: this

        this%nx = 5000 ! Number of points
        this%sou_position = 2500 ! Source position
        ! this%shield_position = 2485 ! Shield position CONSTRUCTIVE INTERFERENCE
        this%shield_position = 2475 ! Shield position DESTRUCTIVE INTERFERENCE
        this%c = sqrt(1/(this%mu0*this%eps0)) ! Speed of light in free space
        this%imp = sqrt(this%mu0/this%eps0) ! Impedance of free space
        this%dz = this%lambda_min/20 ! 20 points per wavelength
        this%dt = this%dz/this%c ! Time step
        
        allocate(this%epr(this%nx), this%mur(this%nx+1), this%Ex(this%nx), this%Hz(this%nx+1)) ! Allocate arrays
        
        this%epr = 1.0 ! Set all relative permittivities to air
        this%mur = 1.0 ! Set all relative permeabilities to free space
        this%Ex = 0.0 ! Set all electric fields to zero
        this%Hz = 0.0 ! Set all magnetic fields to zero
        
        this%w0 = 30000000.0
        ! this%w0 = 2.0*this%pi*this%c/this%lambda_min ! Angular frequency
        this%t0 = this%tau * 6. ! Time delay
        
        this%ice_position = 3000
        this%rock_position = 3500
        this%epr(this%shield_position) = this%shield_permittivity ! Set permittivity of shield
        ! this%epr(this%ice_position:this%nx) = this%ice_permittivity ! Set permittivity of ice
        ! this%epr(this%rock_position:this%nx) = this%rock_permittivity ! Set permittivity of rock

        open(45, file='Ex.txt', status='new') ! Open text files
        open(29, file='Hz.txt', status='new') ! Open text files

    end subroutine EverythingInit

    subroutine EverythingFini(this) ! Deallocate everything
        implicit none
        class(everything), intent(inout) :: this

        deallocate(this%epr, this%mur, this%Ex, this%Hz) ! Deallocate arrays

        close(45) ! Close text files
        close(29) ! Close text files

    end subroutine EverythingFini

    subroutine WriteToText(this) ! Write everything to text files

        implicit none
        class(everything) :: this

        write(45,'(100000(1G13.3e3))') this%Ex ! Write Ex to Ex.txt
        write(29,'(100000(1G13.3e3))') this%Hz ! Write Hz to Hz.txt

    end subroutine WriteToText

    subroutine GaussianSource(this, t) ! Gaussian source
        implicit none
        class(everything), intent(inout) :: this
        integer, intent(in) :: t

        ! Gaussian source

        this%Ex(this%sou_position) = this%Ex(this%sou_position) + exp(-(t-this%t0)**2/this%tau**2) * sin(this%w0*t*this%dt)/this%imp
        this%Hz(this%sou_position) = this%Hz(this%sou_position) - exp(-(t-this%t0)**2/this%tau**2) * sin(this%w0*t*this%dt)/this%imp

    end subroutine GaussianSource
    
    subroutine UpdateMagneticField(this) ! Update magnetic field
        implicit none
        class(everything), intent(inout) :: this
        integer :: i = 0

        do i = 2, this%nx ! Update Hz
            this%Hz(i) = this%Hz(i) + this%dt / (this%mu0 * this%mur(i))*(this%Ex(i) - this%Ex(i-1)) / this%dz
        end do
        
        call BoundaryConditions(this) ! Apply boundary conditions

    end subroutine UpdateMagneticField

    subroutine UpdateElectricField(this) ! Update electric field
        implicit none
        class(everything), intent(inout) :: this
        integer :: i = 0
        
        do i = 1, this%nx ! Update Ex
            this%Ex(i) = this%Ex(i) + this%dt/(this%eps0 * this%epr(i))*(this%Hz(i+1) - this%Hz(i)) / this%dz
        end do
        
        ! No boundary conditions for Ex due to staggered grid stencil

    end subroutine UpdateElectricField    

    subroutine BoundaryConditions(this) ! Boundary conditions
        implicit none
        class(everything), intent(inout) :: this

        ! Boundary conditions
        this%Hz(1) = this%Hz(2) ! Set Hz at the left boundary to the value of the second boundary
        this%Hz(this%nx+1)   = this%Hz(this%nx) ! Set Hz at the second to right boundary to the value of the third to right boundary

    end subroutine BoundaryConditions
    
end module all