!==============================================================================
! 1D model for simulating a turbulent-driven wind from a Red Supergiant.
! Setup following the steady-state descriptions given in Kee+ (2021), A&A, 646.
!
! Coded up by Florian Driessen (2025).
!==============================================================================
module mod_usr

  use mod_hd
  use mod_constants, only: const_c, const_G, const_LSun, const_MSun, &
       const_RSun, mp_cgs, kB_cgs, const_years, const_sigma

  implicit none

  ! User input parameters
  real(8) :: mstar_sol, rstar_sol, twind_cgs, vturb_cgs, kappa_ross_cgs

  ! Dimensionless variables for computations
  real(8) :: mstar, rstar, rhosurf, twind, csound, csoundeff, vturb
  real(8) :: gamma_rad, kappa_ross, gmstar, mdot, rcrit, tfloor, vinf

  ! Hydrogen, helium mass fractions and floor temperature for simulation
  real(8), parameter :: x_h = 1.0d0, y_he = 0.0d0, tfloor_cgs = 100.0d0

  ! Indices of extra output variables
  integer :: itemplucy_, itausph_, igrad_

contains

!==============================================================================
! This routine should set user methods and activate the physics module.
!==============================================================================
  subroutine usr_init

    call set_coordinate_system("spherical")
    call usr_params_read(par_files)

    usr_set_parameters   => initglobaldata_usr
    usr_init_one_grid    => initial_conditions
    usr_special_bc       => special_bound
    usr_gravity          => stellar_gravity
    usr_set_pthermal     => set_ptotal
    usr_process_adv_grid => compute_extra_vars

    call hd_activate()

    itemplucy_ = var_set_extravar("Lucy_temp", "Lucy_temp")
    itausph_   = var_set_extravar("tau_sph", "tau_sph")
    igrad_     = var_set_extravar("grad", "grad")

    if (hd_energy) call mpistop("ERROR: No support for energy equation!")

  end subroutine usr_init

!==============================================================================
! Read in the usr.par file with the problem specific list.
!==============================================================================
  subroutine usr_params_read(files)

    ! Subroutine argument
    character(len=*), intent(in) :: files(:)

    ! Local variable
    integer :: n
    !--------------------------------------------------------------------------

    namelist /star_list/ mstar_sol, rstar_sol, twind_cgs, vturb_cgs, &
         kappa_ross_cgs

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, star_list, end=111)
       111 close(unitpar)
    end do

  end subroutine usr_params_read

!==============================================================================
! Compute some quantities of interest (in CGS) before making unitless.
!==============================================================================
  subroutine initglobaldata_usr

    character(len=*), parameter :: label_fmt = '(A30, " = ", ES12.6)'
    real(8) :: mstar_cgs, rstar_cgs, lstar_cgs, log_g_cgs, log_geff_cgs
    real(8) :: rhosurf_cgs, mumol, csound_cgs, hscale_cgs, vesceff_cgs
    real(8) :: rcrit_cgs, csoundeff_cgs, mdot_cgs, vinf_cgs
    real(8) :: unit_ggrav, unit_mass, alpha, pthsurf
    !--------------------------------------------------------------------------

    mstar_cgs = mstar_sol * const_MSun
    rstar_cgs = rstar_sol * const_RSun

    ! Stellar wind quantities (assume neutral gas composition)
    mumol      = x_h + 0.25d0*y_he
    lstar_cgs  = 4.0d0*dpi * rstar_cgs**2.0d0 * const_sigma * twind_cgs**4.0d0
    gamma_rad  = kappa_ross_cgs * lstar_cgs &
         / (4.0d0*dpi * const_G * mstar_cgs * const_c)
    log_g_cgs   = log10(const_G * mstar_cgs / rstar_cgs**2.0d0)
    csound_cgs = sqrt(twind_cgs * kB_cgs/(mumol * mp_cgs))

    log_geff_cgs  = log_g_cgs + log10(1.0d0 - gamma_rad)
    csoundeff_cgs = sqrt(csound_cgs**2.0d0 + vturb_cgs**2.0d0)
    hscale_cgs    = csoundeff_cgs**2.0d0 / 10.0d0**log_geff_cgs
    vesceff_cgs   = sqrt(2.0d0*const_G * mstar_cgs * (1.0d0 - gamma_rad) &
         / rstar_cgs)
    alpha         = rstar_cgs / hscale_cgs

    ! Kee+ (2021) wind quantities; eqs. 5,11,14 and modified Parker v_infinity
    rcrit_cgs   = 0.5d0*const_G * mstar_cgs * (1.0d0 - gamma_rad) &
          / csoundeff_cgs**2.0d0
    rhosurf_cgs = 2.0d0 / (3.0d0 * kappa_ross_cgs * hscale_cgs) &
         * 1.0d0 / ( 1.0d0 - exp(-rstar_cgs / hscale_cgs) )
    mdot_cgs    = 4.0d0*dpi * rhosurf_cgs * csoundeff_cgs * rcrit_cgs**2.0d0 &
         * exp(-rstar_cgs / hscale_cgs * (1.0d0 - rstar_cgs / rcrit_cgs))    &
         * exp(-0.5d0)
    vinf_cgs    = 2.0d0*csoundeff_cgs * sqrt(log(xprobmax1))

    unit_length        = rstar_cgs
    unit_density       = rhosurf_cgs
    unit_velocity      = csound_cgs
    unit_numberdensity = unit_density / (mumol * mp_cgs)
    unit_pressure      = unit_density * unit_velocity**2.0d0
    unit_temperature   = unit_pressure / (unit_numberdensity * kb_cgs)
    unit_time          = unit_length / unit_velocity
    
    if (mype == 0 .and. .not.convert) then
      write(*, '(A)') repeat('=', 50)
      write(*, '(A)') '   Unity quantities   '
      write(*, '(A)') repeat('=', 50)
      write(*, label_fmt) 'unit length', unit_length
      write(*, label_fmt) 'unit density', unit_density
      write(*, label_fmt) 'unit velocity', unit_velocity
      write(*, label_fmt) 'unit numberdensity', unit_numberdensity
      write(*, label_fmt) 'unit pressure', unit_pressure
      write(*, label_fmt) 'unit temperature', unit_temperature
      write(*, label_fmt) 'unit time', unit_time
      write(*, '(A)') repeat('=', 50)
      write(*, '(A)') '   Stellar and wind parameters in CGS units    '
      write(*, '(A)') repeat('=', 50)
      write(*, label_fmt) 'L/Lsun', lstar_cgs / const_LSun
      write(*, label_fmt) 'M/Msun', mstar_cgs / const_MSun
      write(*, label_fmt) 'R/Rsun', rstar_cgs / const_RSun
      write(*, label_fmt) 'Twind', twind_cgs
      write(*, label_fmt) 'Mean molecular weight', mumol
      write(*, label_fmt) 'Hydrogen mass fraction', x_h
      write(*, label_fmt) 'Helium mass fraction', y_he
      write(*, label_fmt) 'log(g)', log_g_cgs
      write(*, label_fmt) 'log(geff)', log_geff_cgs
      write(*, label_fmt) 'heff/Rstar', hscale_cgs / rstar_cgs
      write(*, label_fmt) 'alpha = Rstar/heff', alpha
      write(*, label_fmt) 'Rosseland opacity', kappa_ross_cgs
      write(*,*)
      write(*, label_fmt) 'isothermal csound', csound_cgs
      write(*, label_fmt) 'turbulent speed', vturb_cgs
      write(*, label_fmt) 'vesc eff', vesceff_cgs
      write(*, label_fmt) 'rcrit/Rstar', rcrit_cgs / rstar_cgs
      write(*, label_fmt) 'Kee rho surface', rhosurf_cgs
      write(*, label_fmt) 'Kee vinf', vinf_cgs
      write(*, label_fmt) 'Kee Mdot', mdot_cgs
      write(*, label_fmt) 'Kee Mdot (sol/yr)', &
           mdot_cgs * const_years / const_MSun
    endif

    ! Code units
    unit_ggrav   = unit_density * unit_time**2.0d0
    unit_mass    = unit_density * unit_length**3.0d0
    unit_opacity = unit_length**2.0d0 / unit_mass

    mstar      = mstar_cgs / unit_mass
    rstar      = rstar_cgs / unit_length
    csound     = csound_cgs / unit_velocity
    vturb      = vturb_cgs / unit_velocity
    csoundeff  = csoundeff_cgs / unit_velocity
    twind      = twind_cgs / unit_temperature
    tfloor     = tfloor_cgs / unit_temperature
    rhosurf    = rhosurf_cgs / unit_density
    pthsurf    = rhosurf * twind
    gmstar     = const_G * unit_ggrav * mstar
    rcrit      = rcrit_cgs / unit_length
    vinf       = vinf_cgs / unit_velocity
    mdot       = mdot_cgs * unit_time / unit_mass
    kappa_ross = kappa_ross_cgs / unit_opacity

    hd_adiab = pthsurf / rhosurf**hd_gamma

    if (mype == 0 .and. .not.convert) then
      write(*, '(A)') repeat('=', 50)
      write(*, '(A)') '  Dimensionless computation quantities  '
      write(*, '(A)') repeat('=', 50)
      write(*, label_fmt) 'Mstar', mstar
      write(*, label_fmt) 'Rstar', rstar
      write(*, label_fmt) 'Wind temperature', twind
      write(*, label_fmt) 'Isothermal sound speed', csound
      write(*, label_fmt) 'Turbulent speed', vturb
      write(*, label_fmt) 'Effective sound speed', csoundeff
       write(*, label_fmt) 'Rosseland opacity', kappa_ross
      write(*, label_fmt) 'Gamma radiation', gamma_rad
      write(*, label_fmt) 'Wind floor temperature', tfloor
      write(*,*)
      write(*, label_fmt) 'Kee critical point', rcrit
      write(*, label_fmt) 'Surface density', rhosurf
      write(*, label_fmt) 'Surface pressure', pthsurf
      write(*, label_fmt) 'Kee terminal wind speed', vinf
      write(*, label_fmt) 'Kee mass-loss rate', mdot
      write(*,*)  
      write(*, label_fmt) 'adiabatic constant', hd_adiab
      write(*, label_fmt) 'adiabatic gamma', hd_gamma
    endif

    ! For cgs output in convert stage
    if (convert .and. saveprim) then
      w_convert_factor(rho_)       = unit_density
      w_convert_factor(mom(1))     = unit_velocity
      w_convert_factor(itemplucy_) = unit_temperature
      w_convert_factor(igrad_)     = unit_length / unit_time**2.0d0
      length_convert_factor        = unit_length
      time_convert_factor          = unit_time
    endif

  end subroutine initglobaldata_usr

!==============================================================================
! Initial conditions start from spherically symmetric 1-D isothermal atmosphere
! with a beta-law velocity.
!==============================================================================
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    real(8) :: temp_lucy(ixO^S), tau_sph(ixO^S)
    real(8) :: sfac
    !--------------------------------------------------------------------------
    
    ! Start at half the sound speed at stellar surface
    sfac = 1.0d0 - (0.5d0*csound / vinf)

    where (x(ixI^S,1) >= rstar)
      w(ixI^S,mom(1)) = vinf * ( 1.0d0 - sfac*rstar / x(ixI^S,1) )
      w(ixI^S,rho_) = mdot / (4.0d0*dpi * w(ixI^S,mom(1)) * x(ixI^S,1)**2.0d0)
    endwhere

    call compute_lucy_temperature(ixI^L,ixO^L,w,x,kappa_ross,temp_lucy,tau_sph)
    w(ixO^S,itemplucy_) = temp_lucy(ixO^S)
    w(ixO^S,itausph_)   = tau_sph(ixO^S)

    call hd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initial_conditions

!==============================================================================
! Special user boundary conditions at inner radial boundary 
!   vr (extrapolated); rho (fixed)
!==============================================================================
  subroutine special_bound(qt, ixI^L, ixB^L, iB, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixB^L, iB
    real(8), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variable
    integer :: ir
    !--------------------------------------------------------------------------

    select case (iB)
    case(1)
      call hd_to_primitive(ixI^L,ixI^L,w,x)

      w(ixB^S,rho_) = rhosurf

      ! Radial velocity field (constant slope extrapolation: d^vr/dr^2 = 0)
      do ir = ixBmax1,ixBmin1,-1
        w(ir^%1ixB^S,mom(1)) = 2.0d0*w(ir+1^%1ixB^S,mom(1)) &
             - w(ir+2^%1ixB^S,mom(1))
      enddo

      ! Prohibit ghosts to be supersonic, also avoid overloading too much
      w(ixB^S,mom(1)) = min(w(ixB^S,mom(1)),  csoundeff)
      w(ixB^S,mom(1)) = max(w(ixB^S,mom(1)), -csoundeff)

      call hd_to_conserved(ixI^L,ixI^L,w,x)

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

!==============================================================================
! Compute the total pressure from thermal pressure, updated based on current
! hydro state during (multi-stage) time step, and the turbulent pressure.
!==============================================================================
  subroutine set_ptotal(w, x, ixI^L, ixO^L, pth)

    ! Subroutine arguments
    integer, intent(in)  :: ixI^L, ixO^L
    real(8), intent(in)  :: x(ixI^S,1:ndim)
    real(8), intent(in)  :: w(ixI^S,1:nw)
    real(8), intent(out) :: pth(ixI^S)

    ! Local variables
    real(8) :: csound(ixI^S)
    real(8) :: temp_lucy(ixO^S), tau_sph(ixO^S)
    !--------------------------------------------------------------------------

    call compute_lucy_temperature(ixI^L,ixO^L,w,x,kappa_ross,temp_lucy,tau_sph)

    csound(ixO^S) = sqrt(temp_lucy(ixO^S))
    csound(ixOmax1+1:ixImax1) = csound(ixOmax1)
    csound(ixImin1:ixOmin1-1) = csound(ixOmin1)

    pth(ixI^S) = (csound(ixI^S)**2.0d0 + vturb**2.0d0) * w(ixI^S,rho_)

  end subroutine set_ptotal

!==============================================================================
! Store extra variables as output.
!==============================================================================
  subroutine compute_extra_vars(igrid, level, ixI^L, ixO^L, qt, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: igrid, level, ixI^L, ixO^L
    real(8), intent(in)    :: qt, x(ixI^S,1:ndim)
    real(8), intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    real(8) :: temp_lucy(ixO^S), tau_sph(ixO^S)
    !--------------------------------------------------------------------------

    call compute_lucy_temperature(ixI^L,ixO^L,w,x,kappa_ross,temp_lucy,tau_sph)

    w(ixO^S,itemplucy_) = temp_lucy(ixO^S)
    w(ixO^S,itausph_)   = tau_sph(ixO^S)
    w(ixO^S,igrad_)     = gmstar * gamma_rad / x(ixO^S,1)**2.0d0

  end subroutine compute_extra_vars

!==============================================================================
! Compute the radiative equilibirum temperature structure for a spherically-
! modified gray atmosphere (Lucy 1971).
!==============================================================================
  subroutine compute_lucy_temperature(ixI^L, ixO^L, w, x, kappa, temp, tau)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: x(ixI^S,1:ndim)
    real(8), intent(in)    :: w(ixI^S,1:nw)
    real(8), intent(in)    :: kappa
    real(8), intent(inout) :: temp(ixO^S), tau(ixO^S)

    ! Local variables
    real(8) :: wdil(ixO^S), tauloc(ixO^S)
    !--------------------------------------------------------------------------

    call get_mspherical_tau(ixI^L,ixO^L,w,x,kappa,tauloc)
    tau(ixO^S) = tauloc(ixO^S)

    ! Dilution factor
    wdil(ixO^S) = 0.5d0*( 1.0d0 - sqrt(1.0d0 - (rstar/x(ixO^S,1))**2.0d0) )

    ! Lucy temperature structure
    temp(ixO^S) = twind * (wdil(ixO^S) + 0.75d0*tau(ixO^S))**0.25d0

    ! Set the minimal cut-off at floor temperature
    temp(ixO^S) = max(temp(ixO^S), tfloor)

  end subroutine compute_lucy_temperature

!==============================================================================
! Compute the spherically-modified optical depth for the temperature structure.
!==============================================================================
  subroutine get_mspherical_tau(ixI^L, ixO^L, w, x, kappa, tau_out)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    real(8), intent(in)    :: x(ixI^S,1:ndim)
    real(8), intent(in)    :: w(ixI^S,1:nw)
    real(8), intent(in)    :: kappa
    real(8), intent(inout) :: tau_out(ixO^S)

    ! Local variables
    integer :: ir
    real(8) :: tau0, dtau(ixO^S), taudum(ixO^S)
    !--------------------------------------------------------------------------

    tau0 = kappa * w(ixOmax1,rho_) * rstar**2.0d0 / x(ixOmax1,1)
    dtau(ixO^S) = -kappa * w(ixO^S,rho_) * (rstar / x(ixO^S,1))**2.0d0

    taudum(ixOmax1) = tau0
    do ir = ixOmax1-1,ixOmin1,-1
      taudum(ir) = taudum(ir+1) &
           + 0.5d0 * ( dtau(ir) + dtau(ir+1) ) * ( x(ir,1) - x(ir+1,1) )
    enddo

    tau_out(ixO^S) = taudum(ixO^S)

  end subroutine get_mspherical_tau

!==============================================================================
! Compute stellar gravity.
!==============================================================================
  subroutine stellar_gravity(ixI^L, ixO^L, wCT, x, gravity_field)

    ! Subroutine arguments
    integer, intent(in)  :: ixI^L, ixO^L
    real(8), intent(in)  :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    real(8), intent(out) :: gravity_field(ixI^S,ndim)
    !--------------------------------------------------------------------------

    gravity_field(ixO^S,:) = 0.0d0
    gravity_field(ixO^S,1) = -gmstar * (1.0d0 - gamma_rad) / x(ixO^S,1)**2.0d0

  end subroutine stellar_gravity

end module mod_usr
