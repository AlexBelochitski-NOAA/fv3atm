! Implementation of the Simplified High Order Closure (SHOC) scheme 
! of Bogenschutz and Krueger (2013), J. Adv. Model. Earth Syst, 5, 195-211,
! https://doi.org/10.1002/jame.20018 (further referred to as BK13)
! in a single column form suitable for use in a GCM physics package. 

! Alex Belochitski - SHOC v1,  heavily based on the code of Peter Bogenschutz.
! S Moorthi - optimization, cleanup, improve and customize for gsm
!           - improved solution for sgs-tke equation
! S Moorthi - 05-11-17 - modified shear production term to eliminate
!                        spurious tke ove Antarctica.
! S Moorthi - 01-12-17 - added extra pressure dependent tke dissipation at
!                        pressures below a critical value pcrit
! S Moorthi - 04-12-17 - fixed a bug in the definition of hl on input
!                        replacing fac_fus by fac_sub
! S.Moorthi - 00-00-17 - added an alternate option for near boundary cek following
!                        Scipion et. al., from U. Oklahoma.
! Alex Belochitski - 2029 - SHOC v2

      

 subroutine shoc(ix, nx, ny, nzm, nz, dtn, me, lat,              &
                 prsl, delp, phii, phil, u, v, omega, tabs,      &
                 qwv, qi, qc, qpi, qpl,                          &
                 dtabsdt, dqwvdt, dqcidt, dqcldt, dqpidt,        &
                 dqpldt, rhc, supice,                            &
                 pcrit, cefac, cesfac, tkef1, dis_opt,           &
                 cld_sgs, tke, hflx, evap, prnum, tkh,           &
                 wthv_sec, lprnt, ipr,imp_phys,     ncpl, ncpi,  &
                 qw_sec, thl_sec, wthl2, wqw2,                   &
!                 detrained_varqt, detrained_varthl, shoc_version)
                 detrained_varqt, detrained_varthl,              &
                 Cek, xkzo, xkzmo, shoc_version, shoc_diag, Diag)
!                 Cek, shoc_version)

  use funcphys , only : fpvsl, fpvsi, fpvs    ! saturation vapor pressure for water & ice

! Map constants of the NCEP GFS to those of SHOC
 
  use physcons, cp    => con_cp,      & ! Specific heat of air, J/kg/K
                ggr   => con_g,       & ! Gravity acceleration, m/s2
                lcond => con_hvap,    & ! Latent heat of condensation, J/kg
                lfus  => con_hfus,    & ! Latent heat of fusion, J/kg
                rv    => con_rv,      & ! Gas constant for water vapor, J/kg/K
                rgas  => con_rd,      & ! Gas constant for dry air, J/kg/K
                pi    => con_pi,      & ! Pi    
                epsv  => con_fvirt
  
  use GFS_typedefs, only :  clear_val, GFS_diag_type

  implicit none

! Input and output variables

  logical lprnt
  integer ipr
  integer, intent(in) :: ix    ! max number of points in the physics window in the x
  integer, intent(in) :: nx    ! Number of points in the physics window in the x
  integer, intent(in) :: ny    ! and y directions
  integer, intent(in) :: me    ! MPI rank
  integer, intent(in) :: lat   ! latitude

  integer, intent(in) :: nzm   ! Number of vertical layers
  integer, intent(in) :: nz    ! Number of layer interfaces  (= nzm + 1)   
  real,    intent(in) :: dtn   ! Physics time step, s 
  integer, intent(in) :: imp_phys ! microphysics identifier
  real,    intent(in) :: pcrit ! pressure in Pa below which additional tke 
                               ! dissipation is applied
  real,    intent(in) :: cefac   ! tunable multiplier to dissipation term 
  real,    intent(in) :: cesfac  ! tunable multiplier to dissipation term for bottom level 
  real,    intent(in) :: tkef1   ! uncentering terms in implicit tke integration         
  real,    intent(in) :: dis_opt ! when > 0 use different formula for near surface dissip.

  real,    intent(in) :: hflx(nx) ! Flux of sensible heat at the surface, K*m/s
  real,    intent(in) :: evap(nx) ! Flux of latent heat at the surface, kg/kg**m/s 

  real, intent(in)    :: prsl   (ix,ny,nzm)   ! mean layer presure   
  real, intent(in)    :: phii   (ix,ny,nz )   ! interface geopotential height
  real, intent(in)    :: delp   (ix,ny,nzm)      ! layer presure depth
  real, intent(in)    :: phil   (ix,ny,nzm)   ! layer geopotential height  
  real, intent(in)    :: u      (ix,ny,nzm)   ! u-wind, m/s
  real, intent(in)    :: v      (ix,ny,nzm)   ! v-wind, m/s
  real, intent(in)    :: omega  (ix,ny,nzm)   ! omega, Pa/s
  real, intent(inout) :: tabs   (ix,ny,nzm)   ! temperature, K
  real, intent(inout) :: qwv    (ix,ny,nzm)   ! water vapor mixing ratio, kg/kg
  real, intent(inout) :: qc     (ix,ny,nzm)   ! cloud water mixing ratio, kg/kg
  real, intent(inout) :: qi     (ix,ny,nzm)   ! cloud ice   mixing ratio, kg/kg
  real, intent(inout) :: qpl    (nx,ny,nzm)   ! rain mixing ratio, kg/kg
  real, intent(inout) :: qpi    (nx,ny,nzm)   ! snow + graupel mixing ratio, kg/kg
  real, intent(inout) :: dtabsdt(ix,ny,nzm)   ! Tendency of temperature due to diffusion, K/s
  real, intent(inout) :: dqwvdt (ix,ny,nzm)   ! Tendency of water vapor mixing ratio due to diffusion, kg/kg/s
  real, intent(inout) :: dqcldt (ix,ny,nzm)   ! Tendency of cloud water mixing ratio due to diffusion, kg/kg/s
  real, intent(inout) :: dqcidt (ix,ny,nzm)   ! Tendency of cloud ice mixing ratio due to diffusion, kg/kg/s
  real, intent(inout) :: dqpldt (ix,ny,nzm)   ! Tendency of rain  mixing ratio due to diffusion, kg/kg/s
  real, intent(inout) :: dqpidt (ix,ny,nzm)   ! Tendency of snow + graupel mixing ratio due to diffusion, kg/kg/s
! Anning Cheng 03/11/2016 SHOC feedback to number concentration
  real, intent(inout) :: ncpl   (nx,ny,nzm)   ! cloud water number concentration,m**-3
  real, intent(inout) :: ncpi   (nx,ny,nzm)   ! cloud ice   number concentration,m**-3
  real, intent(inout) :: rhc    (nx,ny,nzm)   ! critical relative humidity
  real, intent(in)    :: supice               ! ice supersaturation parameter
  real, intent(inout) :: cld_sgs(ix,ny,nzm)   ! sgs cloud fraction
  real, intent(inout) :: tke    (ix,ny,nzm)   ! turbulent kinetic energy. m**2/s**2
  real, intent(inout) :: tkh    (ix,ny,nzm)   ! eddy diffusivity for heat
  real, intent(inout) :: prnum  (nx,ny,nzm)   ! turbulent Prandtl number
  real, intent(in)    :: xkzo   (nx,ny,nzm)   ! Background diffusivity for heat
  real, intent(in)    :: xkzmo  (nx,ny,nzm)   ! Background diffusivity for momentum
  real, intent(inout) :: wthv_sec (ix,ny,nzm) ! Buoyancy flux, K*m/s
  real, intent(inout) :: Cek    (nx,ny,nzm)   ! Coefficient in the TKE dissipation term 

!  real :: xkzo   (nx,ny,nzm)   ! Background diffusivity for heat
!  real :: xkzmo  (nx,ny,nzm)   ! Background diffusivity for momentum
  real :: tkh_out  (ix,ny,nzm)   ! eddy diffusivity for heat
  real :: tkh_out1  (ix,ny,nzm)   ! eddy diffusivity for heat

! SHOC v2 support  

! Second moment total water mixing ratio, kg^2/kg^2
  real, intent(inout) :: qw_sec  (nx,ny,nzm) 
! Second moment liquid/ice static energy, K^2
  real, intent(inout) :: thl_sec (nx,ny,nzm) 
! Turbulent flux of variance of liquid/ice static energy, m/s*K^2
  real, intent(inout) :: wthl2   (nx,ny,nzm) 
! Turbulent flux of variance liquid/ice static energy, K*m/s
  real, intent(inout) :: wqw2    (nx,ny,nzm) 
! Tendency of total water variance due to convective detrainment
  real, intent(inout) :: detrained_varqt  (nx,ny,nzm) 
! Tendency of MSE variance due to convective detrainment
  real, intent(inout) :: detrained_varthl (nx,ny,nzm) 
! Version of SHOC to use, 1 or 2 currently 
  integer, intent(in) :: shoc_version       
! Flag for saving SHOC diagnostic output
  logical, intent(in) :: shoc_diag       
! GFS diagnostic output
  type(GFS_diag_type), intent(inout) :: Diag   


! Arithmetic constants
  real, parameter :: zero=0.0,  one=1.0,  half=0.5, two=2.0,    eps=0.622,           &
                     three=3.0, oneb3=one/three, twoby3=two/three, sqrt2=sqrt(two),  &
                     sqrtpii = one/sqrt(pi+pi),  twoby15 = two / 15.0, pt01=0.01

! Physical constants
  real, parameter :: lsub = lcond+lfus, fac_cond = lcond/cp, fac_fus = lfus/cp,      &
                     cpolv = cp/lcond,  fac_sub = lsub/cp, ggri = 1.0/ggr,           &
                     kapa = rgas/cp, gocp = ggr/cp, rog = rgas*ggri,                 &
                     epsterm = rgas/rv, onebeps = one/epsterm, onebrvcp= one/(rv*cp),&
                     vonk=0.4, &    ! Von Karman constant Moorthi - as in GFS
                     Pr    = 1.0, & ! Prandtl number
                     scrit=2.0e-6

!---------------------------------------------------------------------------------------
!                                       SHOC's TUNABLE PARAMETERS
!---------------------------------------------------------------------------------------

!--------------------------------
! Constraints on the sub-grid PDF
!--------------------------------

! Minimum and maximum values of individual gaussians' weights
  real, parameter :: atmin=0.01, atmax=one-atmin 
  real, parameter :: atmin_damp=0.05, atmax_damp=one-atmin_damp 
! Damping coefficient in the TKE equation (and the "return to isotropy" damping time scale
! in the prognostic variance equations) will be multiplied by (at_damp_strength+1) for the
! PDFs with values of the individual gaussian weights in the ranges (atmin:atmin_damp) and
!  (atmax:atmax_damp) 
  real, parameter :: at_damp_strength =  3.
! Maximum normalized variance of individual MSE gaussians 
  real, parameter :: max_thl_norm_var = 3.
! Maximum normalized variance of individual total water gaussians   
  real, parameter :: max_qw_norm_var =  3.    
! Maximum absolute value of skewness of W SGS PDF
  real, parameter :: max_w_skw = 10.  

! Larson, V.E. and J. Golaz, 2005, Mon. Wea. Rev., 133, 1023–1042
! https://doi.org/10.1175/MWR2902.1  (further referred to as LG05)  
! introduced  a parameterization of variable normalized width of individual W gaussians
! (as opposed to a fixed value of 0.4 in the original formulations of SHOC and CLUBB).
! If the switch below is set to  .true., the code uses the parameterization.
  logical, parameter :: variable_normalized_width_w = .false. !.true.
! Paramter, controlling magnitude of normalized variances
  real,    parameter :: gamma_vnww = 0.32 ! 0 < gamma_vnww < 1

!--------------------------------------------------------------------------
! Options for parameterizations of skewness of MSE and total water SGS PDFs
!--------------------------------------------------------------------------

! In the original formulation of SHOC skewness of MSE and total water 
! SGS PDFs is set to be linearly dependent on the skewness of W SGS PDF. 
! Will not be used if any other options for skewness parameterizaion are chosen.
! Note that  Zhu and Zuidema (2009) show that Sk of scalars is not proportional to Sk W in
! boundary layer clouds. 
! Coefficient relating skewness of total water and W SGS PDFs
  real, parameter :: skew_facw=1.2
! Coefficient relating skewness of MSE and W SGS PDFs 
  real, parameter :: skew_fact=0.0   

! Alternative parameterization for skewness of MSE PDF is based on diagnostic 
! expressions for the third moment of MSE using the Canuto et al (2001) approach.  
! If the switch below is set to  .true., the code uses the expression
! based on Canuto et al (2001)
  logical, parameter :: canuto_skew  = .false. 
! Proportionality coefficient between skewness of MSE and skewness of total water, used only
! when  canuto_skew  = .true., in enssence parametrizing skewness of Qt in terms of MSE Sk. 
! The basis for this parameterization is the observation of Zhu and Zuidema (2009) that in 
! LES simulations of stratocumulus and shallow convection skewnesses of PDFs of scalars
! have no relation to the Sk of W, but rather are proportional to each other with the opposite
! sign. ( See Fig. 2 of their paper)
  real,    parameter :: canuto_factor    = -5.
! Maximun absolute value of MSE Sk, based on Zhu and Zuidema (2009)
  real,    parameter :: canuto_sk_hl_max = 20. 


! LG05 introduced parameterizations for skewness of total 
! water and MSE PDFs which lead to a "realizable" set of moments 
! unlike the original formulation. 
! Guo, Z., et al 2015, , J. Adv. Model. Earth Syst., 7, 1005–1025, 
! https://doi.org/10.1002/2014MS000405 
! find that  gamma_vnww and beta_factor are useful in controlling
! stratocumulus in global NCAR CAM simulations with CLUBB. 
! When larson_golaz_05_skew is set to .true.,  parameterization of variable normalized 
! width of individual W gaussians is also used.
  logical, parameter :: larson_golaz_05_skew        = .false. !.true. 
! Parameter, controlling relationship betwen skewness of W SGS PDFs with other Sk of scalars.
  real,    parameter :: beta_factor  = 1.75 ! 0 <= beta_factor  <= 3

!------------------------------------
! Contrains on the moments of SGS PDF
!------------------------------------
 
  real, parameter :: w_tol     = 1e-3 !2.0e-02    ! Min value of sqrt(second moment of w), m/s   
  real, parameter :: w_tol_sqd = w_tol*w_tol ! Min value of second moment of w squared
  real, parameter :: w3_tol = 1e-20         ! Min value of third moment of w
! Distance between normalized means of W gaussians below which PDF turns into double delta function
  real, parameter :: w_thresh  = 0.0     
  real, parameter :: thl_tol  = 1e-4 !1.e-2  ! Min value of sqrt(second moment of MSE), K
  real, parameter :: rt_tol = 1e-8 !1.e-4    ! Min value of sqrt(second moment of tot. water),g/g
!                    thl_tol  = 1.e-2, rt_tol = 1.e-8 !CLUBB values
! Tuning coefficients in diagnostic expressions for higher order moments
  real, parameter :: thl2tune = 1.0,  qw2tune = 1.0,  qwthl2tune = 1.0

!------------------------------------------ 
! Parameters of the eddy length formulation
!------------------------------------------

! Maximum turbulent eddy length scale, m
!  real, parameter :: max_eddy_length_scale  = 2000.0 
  real, parameter :: max_eddy_length_scale  = 3000.0 
! Empirical time scale controlling eddy length in the boundary layer, s
! The larger the tscale, the larger the eddy length 
  real, parameter :: tscale = 400.
!  real, parameter :: tscale = 300.
! Parameter controlling contribution of strength of local thermal stability
! to eddy length scale under stable conditions 
  real, parameter :: slts = pt01
! Maximum measure of PBL height, m * 1e-2
  real, parameter :: max_l_inf = 100.0
!  real, parameter :: max_l_inf = 300.0

! Maximum "return-to-isotropy" time scale, s
!  real, parameter :: max_eddy_dissipation_time_scale =  400.0
  real, parameter :: max_eddy_dissipation_time_scale =  900.0
! Paramteter controlling magnitude of "return-to-isotropy" time scale 
! under stable conditions
  real, parameter :: lambda  = 0.04
! Threshold on the minimum amount of condensate that SHOCv2 will consider to be cloud to avoid 
! situations  where SHOC uses its in-cloud length scale on levels too close to the surface, kg/kg
  real, parameter :: thresh =  1.e-8
  real, parameter :: qt_thresh =  5.e-6
  real, parameter :: qi_thresh =  1.e-8
!  Set cloud top as the level with the most negative boyuancy flux, otherwise 
!  cloud top is at the level with cloud condensate below threshold thresh
  logical, parameter :: cloud_top_at_min_negative_boyu_flux = .true.

!------------------------------------------
! Parameters of the prognostic TKE equation
!------------------------------------------

! Number of iterations in the TKE equation solver.
! First order semi-implicit backwards Euler method is used to take one step forward in time
! and then the solution is iterated for convergence. 
  integer, parameter :: nitr = 6

  real, parameter :: min_tke = 1e-8  ! Minumum TKE value, m**2/s**2 
  real, parameter :: max_tke = 100.0 ! Maximum TKE value, m**2/s**2 
  real, parameter :: epsln=1.0e-6    ! Minumum value of eddy diffusivity
  real, parameter :: tkhmax=300.0    ! Maximum value of eddy diffusivity

! Constants for the TKE dissipation term based on Deardorff (1980)
  real, parameter :: pt19=0.19, pt51=0.51

! Constants for the TKE dissipation term 
  real, parameter :: Cs  = 0.15
  real, parameter :: Ck  = 0.1     
  real, parameter :: Ce  = 4*Ck**3/Cs**4  ! Coefficient in the TKE diss. term
!  real, parameter :: Ce  = Ck**3/Cs**4  ! Coefficient in the TKE diss. term
  real, parameter :: Ces = Ce             ! Coefficient in the TKE diss. term at sfc
! Original values below, Ce ~ 2, Ces ~ 8.5
! real, parameter :: Ce  = Ck**3/Cs**4, Ces = Ce*3.0/0.7

!--------------------------------------------------------------------
! Parameters of prognostic variances of MSE and total water equations
!--------------------------------------------------------------------     

! Type of prognostic variances solver:
! 1 - first order implicit backwards Euler method  
! 2 - second order implicit backwards Euler method
  integer, parameter :: prog_var_solver_type = 2

! Number of substeps in the prognostic variances solver
  integer, parameter :: substep=10

  real, parameter :: Cv = 4*1.04   ! Damping coefficient in prognostic variance eqns
  real, parameter :: mu = 20    ! Numerical viscosity, m**2/s, curently not used

! Tuning parameters for variance source term from detrainment
  real, parameter  :: tune_varthl=0., tune_varqt=0.

!-----------------------------------
! Tie-in with a microphysical scheme
!-----------------------------------

! Minimum T below which condensate is all ice
!  real, parameter :: tbgmin = 258.16    ! Minimum temperature for cloud water., K (ZC)
  real, parameter :: tbgmin = 253.16    ! Minimum temperature for cloud water., K (ZC)
!  real, parameter :: tbgmin = 233.16    ! Minimum temperature for cloud water., K (ZC)  
! Maximum T above which condensate is all water
  real, parameter :: tbgmax = 273.16    ! Maximum temperature for cloud ice, K
! Linear interpolation in between
  real, parameter :: a_bg   = one/(tbgmax-tbgmin)
! Condesation options
  logical, parameter :: Firl_condensation = .false. !.true.  

!--------------------------------------------------------------
! Mode of interpolation from layer centers to layer interfaces:
!--------------------------------------------------------------
! 0 - Linear (half-sums in the vertical)
! 1 - Monotone piecewise cubic Hermite interpolant
  integer, parameter :: interp_mode = 1



! Local variables. Note that pressure is in millibars in the SHOC code.

  real zl      (nx,ny,nzm)  ! height of the layer centers above surface, m 
  real zi      (nx,ny,nz)   ! height of the layer interfaces, m          
  real adzl    (nx,ny,nzm)  ! layer thickness i.e. zi(k+1)-zi(k) - defined at levels, m
  real adzi    (nx,ny,nz)   ! level thickness i.e. zl(k)-zl(k-1) - defined at interface, m
 
  real hl      (nx,ny,nzm)  ! liquid/ice moist static energy , K
  real qv      (nx,ny,nzm)  ! water vapor mixing ratio, kg/kg
  real qcl     (nx,ny,nzm)  ! cloud water mixing ratio, kg/kg
  real qci     (nx,ny,nzm)  ! cloud ice   mixing ratio, kg/kg
  real w       (nx,ny,nzm)  ! z-wind, m/s
  real bet     (nx,ny,nzm)  ! ggr/tv0
  real gamaz   (nx,ny,nzm)  ! ggr/cp*z

! Moments of the trivariate double Gaussian PDF for the SGS total water mixing ratio
! SGS liquid/ice static energy, and vertical velocity

! Second moments of MSE and total water are prognostic variables in SHOC v2
!  real qw_sec   (nx,ny,nzm) ! Second moment total water mixing ratio, kg^2/kg^2
!  real thl_sec  (nx,ny,nzm) ! Second moment liquid/ice static energy, K^2
  real qwthl_sec(nx,ny,nzm) ! Covariance tot. wat. mix. ratio and static energy, K*kg/kg
  real wqw_sec  (nx,ny,nzm) ! Turbulent flux of tot. wat. mix., kg/kg*m/s
  real wthl_sec (nx,ny,nzm) ! Turbulent flux of liquid/ice static energy, K*m/s
  real wqw_sec_d  (nx,ny,nzm) ! Turbulent flux of tot. wat. mix., kg/kg*m/s
  real wthl_sec_d (nx,ny,nzm) ! Turbulent flux of liquid/ice static energy, K*m/s
  real w_sec    (nx,ny,nzm) ! Second moment of vertical velocity, m**2/s**2
  real w3       (nx,ny,nzm) ! Third moment of vertical velocity, m**3/s**3
  real wqp_sec  (nx,ny,nzm) ! Turbulent flux of precipitation, kg/kg*m/s
  real hl3      (nx,ny,nzm) ! Third moment of MSE, K**3

! Eddy length formulation 
  real smixt    (nx,ny,nzm) ! Turbulent length scale, m
  real isotropy (nx,ny,nzm) ! "Return-to-isotropy" eddy dissipation time scale, s
  real brunt    (nx,ny,nzm) ! Moist Brunt-Vaisalla frequency, s^-1
  real conv_vel2(nx,ny,nzm) ! Convective velocity scale cubed, m^3/s^3
  real tkm      (ix,ny,nzm) ! eddy diffusivity for momentum

! Coefficient in the TKE dissipation term
  real Cek_tke(nx,ny,nzm)

! Variables interpolated to layer interfaces
  real tke_int      (nx,ny,nzm) ! turbulent kinetic energy interpolated to layer interfaces, m**2/s**2  
  real w_sec_int    (nx,ny,nzm) ! Second moment of vertical velocity interpolated to layer interfaces, m**2/s**2     
  real isotropy_int (nx,ny,nzm) ! "Return-to-isotropy" eddy dissipation time scale interpolated to layer interfaces, s  
  real brunt_int    (nx,ny,nzm) ! Moist Brunt-Vaisalla frequency interpolated to layer interfaces, s^-1 
  real bet_int      (nx,ny,nzm) ! ggr/tv0 interpolated to layer interfaces

! Output of SHOC

 real diag_frac ! SGS cloud fraction
 real diag_qn   ! SGS cloud+ice condensate, kg/kg 
 real diag_qi   ! SGS ice condensate, kg/kg
 real diag_ql   ! SGS liquid condensate, kg/kg                                                                                              
!  real diag_frac, diag_qn, diag_qi, diag_ql

! real diag_frac(nx,ny,nzm) ! SGS cloud fraction
! real diag_qn  (nx,ny,nzm) ! SGS cloud+ice condensate, kg/kg
! real diag_qi  (nx,ny,nzm) ! SGS ice condensate, kg/kg
! real diag_ql  (nx,ny,nzm) ! SGS liquid condensate, kg/kg


! Horizontally averaged variables
! real conv_vel(nzm)        ! Convective velocity scale cubed, m^3/s^3
  real wqlsb   (nzm)        ! liquid water flux, kg/kg/ m/s
  real wqisb   (nzm)        ! ice flux, kg/kg m/s


! Local variables

! real, dimension(nx,ny,nzm) :: tkesbbuoy, tkesbshear, tkesbdiss, tkesbbuoy_debug   &
!                               tkebuoy_sgs, total_water, tscale1_debug, brunt2

  real, dimension(nx,ny,nzm) :: total_water, brunt2, thv, tkesbdiss, brunt_test
  real, dimension(nx,ny,nzm) :: def2
!  real, dimension(nx,ny)     :: denom, numer, l_inf, cldarr, thedz, thedz2, sflux
  real, dimension(nx,ny)     :: denom, numer, l_inf, thedz, thedz2, sflux
!  logical, dimension(nx,ny)  :: cbl_top 
  integer, dimension(nx,ny)  :: cbl_top, cldarr

  real, dimension(nx,ny,nzm) :: dthldz, dqwdz, dwthl2dz, dwqw2dz, deltathl, deltaqw

  real, dimension(nx,ny,nzm) :: avew_save, cond_w_save, z_save, w3_save, w3var_save, w_sec_save, skew_w_save


  real lstarn,    depth,    omn,         betdz,      bbb,      term,   qsatt, dqsat,        &
                  conv_var,  tkes,       skew_w,     skew_qw,  aterm,  w1_1,  w1_2,  w2_1,  &
       w2_2,      w3var,     thl1_1,     thl1_2,     thl2_1,   thl2_2, qw1_1, qw1_2, qw2_1, &
       qw2_2,     ql1,       ql2,        w_ql1,      w_ql2,                                 &
       r_qwthl_1, r_wqw_1,   r_wthl_1,   testvar,    s1,    s2,    std_s1,  std_s2, C1, C2, &
       thl_first, qw_first,  w_first,    Tl1_1,      Tl1_2, betatest, pval, pkap,           &
       w2thl,     w2qw,w2ql, w2ql_1,     w2ql_2,                                            &
       thec,      thlsec,    qwsec,      qwthlsec,   wqwsec, wthlsec, thestd,dum,           &
       cqt1,      cthl1,     cqt2,       cthl2,      qn1,    qn2, qi1, qi2, omn1, omn2,     &
       basetemp2, beta1,     beta2,      qs1,        qs2,                                   &
       esval1_1,  esval2_1,  esval1_2,   esval2_2,   om1,    om2,                           &
       lstarn1,   lstarn2,   sqrtw2,     sqrtthl,    sqrtqt,                                &
       sqrtstd1,  sqrtstd2,  tsign,      tvar,       sqrtw2t, wqls, wqis,                   &
       sqrtqw2_1, sqrtqw2_2, sqrtthl2_1, sqrtthl2_2, sm,   prespot,                         &
       corrtest1, corrtest2, wrk,  wrk1, wrk2, wrk3, onema, pfac, sfac, sfaci, dts, Cvv,    &
       hl3var, skew_hl, save_thl_sec, save_qw_sec, tkef2,  rdzw, wrku, wrkv, wrkw,          &
       grd, buoy_sgs, ratio, a_prod_sh, a_prod_bu, smix, adiss, tscale1, wtke, wtk2, rdtn,  &
       Cee, a_diss, T1_1,      T1_2

  real, allocatable :: tke_history(:)

  integer i,j,k,km1,ku,kd,kdd,ka,kb, itr

  logical, save :: first_call = .true.

  logical ::  flag

  flag = .true.

! A cludge to deal with GFS restarts that do not have SHOC's tracers
  if (first_call) then
     
     tke     = 0.
     qw_sec  = 0.
     thl_sec = 0.
     first_call = .false.

! Initialize damping coefficients for the TKE eqn on the first time step
!        if (Cek(i,j,k) == clear_val) then
!           if (k == 1) then
!              Cek(i,j,k) = ces * cesfac
!           else
!              Cek(i,j,k) = ce  * cefac
!           endif
!        endif

        Cek = 1.  !Cv

  endif

! Map GFS variables to those of SHOC - SHOC operates on 3D fields
! Here a Y-dimension is added to the input variables, along with some unit conversions

  do k=1,nz
    do j=1,ny
      do i=1,nx
! Geometric height of a layer interface
        zi(i,j,k) = phii(i,j,k) * ggri
      enddo
    enddo
  enddo
!
! Move water from vapor to condensate if the condensate is negative
!

! Old code. Treatment of vater vapor is wrong. 
!  do k=1,nzm
!    do j=1,ny
!      do i=1,nx
!        if (qc(i,j,k) < zero) then
!          wrk = qwv(i,j,k) + qc(i,j,k)
!          if (wrk >= zero) then
!            qwv(i,j,k)  = wrk
!            tabs(i,j,k) = tabs(i,j,k) - fac_cond * qc(i,j,k)
!            qc(i,j,k)   = zero
!          else
!            qc(i,j,k)   = zero
!            tabs(i,j,k) = tabs(i,j,k) + fac_cond * qwv(i,j,k)
!            qwv(i,j,k)  = zero
!          endif
!        endif
!        if (qi(i,j,k) < zero) then
!          wrk = qwv(i,j,k) + qi(i,j,k)
!          if (wrk >= zero) then
!            qwv(i,j,k)  = wrk
!            tabs(i,j,k) = tabs(i,j,k) - fac_sub  * qi(i,j,k)
!            qi(i,j,k)   = zero
!          else
!            qi(i,j,k)   = zero
!            tabs(i,j,k) = tabs(i,j,k) + fac_sub  * qwv(i,j,k)
!            qwv(i,j,k)  = zero
!          endif
!        endif
!      enddo
!    enddo
!  enddo
             

  do k=1,nzm
    do i=1,nx
      if (qc(i,1,k) < zero) then
        qwv(i,1,k)  = qwv(i,1,k) + qc(i,1,k)
        tabs(i,1,k) = tabs(i,1,k) - fac_cond * qc(i,1,k)
        qc(i,1,k)   = zero
      endif
      if (qi(i,1,k) < zero) then
        qwv(i,1,k)  = qwv(i,1,k) + qi(i,1,k)
        tabs(i,1,k) = tabs(i,1,k) - fac_sub  * qi(i,1,k)
        qi(i,1,k)   = zero
      endif
    enddo
  enddo
! If water wapor is negative, bring it up to zero with 
  do k=nzm,2,-1
    km1 = k - 1
    do i=1,nx
      if (qwv(i,1,k) < zero) then
! This does not look right, water wapor is not conserved? 
        qwv(i,1,k) = qwv(i,1,km1) + qwv(i,1,k) * delp(i,1,k) / delp(i,1,km1)
      endif
    enddo
  enddo



  do k=1,nzm
    do j=1,ny
      do i=1,nx

        qv(i,j,k)    = max(qwv(i,j,k), zero)
        qcl(i,j,k)   = max(qc(i,j,k), zero)
        qci(i,j,k)   = max(qi(i,j,k), zero)

!        if  (k .eq. 2) then 
!            if (qcl(i,j,k) > 0. ) print *,   qcl(i,j,k)
!            if (qci(i,j,k) > 0. ) print *,  qci(i,j,k)
!         end if
               

        qpl(i,j,k)     = zero  ! comment or remove when using with prognostic rain/snow
        qpi(i,j,k)     = zero  ! comment or remove when using with prognostic rain/snow
        wqp_sec(i,j,k) = zero  ! SGS flux of precipiation

! Height of a layer center
        zl(i,j,k)    = phil(i,j,k) * ggri
        wrk          = one / prsl(i,j,k)

! Virtual temperature
!        thv(i,j,k)   = tabs(i,j,k) * (one+epsv*qv(i,j,k)) 
        thv(i,j,k)=tabs(i,j,k)*(one+epsv*qv(i,j,k) - qcl(i,j,k) - qci(i,j,k) &
                                                   - qpl(i,j,k) - qpi(i,j,k) ) 

! This convesrion of omega to w assumes hydrostatic balance. Should be fixed for when 
! the non-hydrostatic model is run at resolutions of ~4km and higher. 
! Perhaps it would be safer to just use dz/dt to begin with since it's avalaible in FV3 
        w(i,j,k)     = - rog * omega(i,j,k) * thv(i,j,k) * wrk


        total_water(i,j,k) = qcl(i,j,k) + qci(i,j,k) + qv(i,j,k)


        prespot        = (100000.0*wrk) ** kapa ! Exner function
        bet(i,j,k)     = ggr/(tabs(i,j,k)*prespot)     ! Moorthi
! Virtual potential temperature
        thv(i,j,k)     = thv(i,j,k)*prespot            ! Moorthi

! Dry adiabatic lapse rate * height = reference temperature
        gamaz(i,j,k) = gocp * zl(i,j,k)

! Liquid/ice water static energy. Note that the units are degrees K
        hl(i,j,k) = tabs(i,j,k) + gamaz(i,j,k) - fac_cond*(qcl(i,j,k)+qpl(i,j,k)) &
                                               - fac_sub *(qci(i,j,k)+qpi(i,j,k))
        w3(i,j,k) = zero

        brunt(i,j,k) = zero


!        tkh(i,j,k)= tkh(i,j,k) + xkzo(i,j,k)
!        tkm(i,j,k)= tkh(i,j,k)*prnum(i,j,k) + xkzmo(i,j,k)
        
      enddo
    enddo
  enddo



! Apply background diffusivity to SHOC's diffusivity coefficients
  do k=1,nzm-1
!  do k=1,nzm
     ku = k + 1
!     ku = k
     do j=1,ny
        do i=1,nx

!           tkh(i,j,ku)= tkh(i,j,ku)  + xkzo(i,j,k)
!           tkm(i,j,ku)= tkh(i,j,ku)*prnum(i,j,ku) + xkzmo(i,j,k)
            tkh(i,j,ku)= max(tkh(i,j,ku),xkzo(i,j,k))
            tkm(i,j,ku)= max(tkh(i,j,ku)*prnum(i,j,ku),xkzmo(i,j,k))

        enddo
     enddo
  enddo

  tkh(:,1,1)= tkh(:,1,2)
  tkm(:,1,1)= tkm(:,1,2)


  do j=1,ny
     do i=1,nx
        bet_int(i,j,:) = interp_center_to_interface(interp_mode, & 
                                                    bet(i,j,:), zl(i,j,:), zi(i,j,:))
     enddo ! i
  enddo    ! j

! Define vertical grid increments for later use in the vertical differentiation

  do k=2,nzm
    km1 = k - 1
    do j=1,ny
      do i=1,nx
        adzi(i,j,k)   = zl(i,j,k) - zl(i,j,km1)
        adzl(i,j,km1) = zi(i,j,k) - zi(i,j,km1)
      enddo
    enddo
  enddo
  do j=1,ny
    do i=1,nx
      adzi(i,j,1)     = (zl(i,j,1)-zi(i,j,1)) 
      adzi(i,j,nz)    = adzi(i,j,nzm)         
      adzl(i,j,nzm)   = zi(i,j,nz) - zi(i,j,nzm)

! Boundary conditions for fluxes of MSE and total water
      wthl_sec(i,j,1) = hflx(i) 
      wqw_sec(i,j,1)  = evap(i)
      wthl_sec_d(i,j,1) = hflx(i) 
      wqw_sec_d(i,j,1)  = evap(i)
!     wqp_sec(i,j,1)  = zero
    enddo
  enddo

  do k=2,nzm
       
    km1 = k - 1
    if (k==2) then 
       ku=3
       kd=2
    elseif (k==nzm) then

       ku=nzm
       kd=nzm-1

    else

       ku=k
       kd=k-1
    endif

    do j=1,ny
      do i=1,nx

! Calculate diffusive flux of MSE from diffusive fluxes of temperature, moisture, 
! and condensate computed by the vertical diffusion solver
        wthl_sec(i,j,k) =  wthl_sec(i,j,k-1) - adzl(i,j,k-1)*(dtabsdt(i,j,k-1)     -  &
                          fac_cond*(dqcldt(i,j,k-1) + dqpldt(i,j,k-1))            -  &
                          fac_sub* (dqcidt(i,j,k-1) + dqpidt(i,j,k-1)))            -  &
                          gocp*(tkh(i,j,ku) - tkh(i,j,kd))

      end do
   end do
end do


! Integrate prognostic TKE equation forward in time 

  rdtn = one / dtn
  tkef2 = 1 - tkef1

! Calculate shear production terms at layer centers for prognostic TKE and variance
! equations

  do k=1,nzm
    ku = k+1
    kd = k-1
    if (k == 1) then
      kd = k
    elseif (k == nzm) then
      ku = k
    endif

     do j=1,ny
        do i=1,nx
           rdzw        = one /  (zl(i,j,ku) - zl(i,j,kd))
           wrku        = (u(i,j,ku)-u(i,j,kd)) * rdzw
           wrkv        = (v(i,j,ku)-v(i,j,kd)) * rdzw
!          wrkw        = (w(i,j,ku)-w(i,j,kd)) * rdzw
           def2(i,j,k) = wrku*wrku + wrkv*wrkv !+ 2*wrkw(1) * wrkw(1)

           if (shoc_version == 2) then
             
! Shear generation terms for prognostic treatment of variances of total water
! and moist static energy.

! Vertical gradient of MSE
              dthldz(i,j,k)   =  (hl(i,j,ku)-hl(i,j,kd))*rdzw
! Vertical gradient of total water
              dqwdz(i,j,k)    =  (total_water(i,j,ku)-total_water(i,j,kd))*rdzw
! Divergence of subgrid vertical flux of MSE variance
              dwthl2dz(i,j,k) =  (wthl2(i,j,ku)-wthl2(i,j,kd))*rdzw
! Divergence of subgrid vertical flux of total watera variance
              dwqw2dz(i,j,k) =   (wqw2(i,j,ku)-wqw2(i,j,kd))*rdzw

           end if

        enddo   ! i
     enddo      ! j 
  enddo         ! k  loop


! Ensure values of TKE are reasonable

  do k=1,nzm
     do j=1,ny
        do i=1,nx
           tke(i,j,k)        = max(min_tke,tke(i,j,k))
           tkesbdiss(i,j,k)  = zero
!         tkesbshear(i,j,k) = zero
!         tkesbbuoy(i,j,k)  = zero
        enddo
     enddo
  enddo


!           a_prod_bu = ggr / thv(i,j,k) * wthv_sec(i,j,k)

! Obtain Brunt-Vaisalla frequency from diagnosed SGS buoyancy flux
! (Note that there's another BV freq. calculation in  eddy_length())


  do k=1,nzm

     ku = k+1
     kd = k
      
     if(k == 1) then
        ku = 2
        kd = 2
     elseif(k == nzm) then
        ku = k
        kd = k
     endif

     do j=1,ny
        do i=1,nx
           if (k == 1) then 
! On the first level use surface byoancy flux
              sflux(i,j)=hflx(i)+evap(i)*epsv*tabs(i,j,1)*(100000.0/prsl(i,k,1))**kapa
!              brunt(i,j,1) = - 2*ggr*sflux(i,j)/(thv(i,j,k)*(tkh(i,j,ku)+tkh(i,j,kd) + 0.0001))
             brunt(i,j,k) = - 2*ggr*wthv_sec(i,j,k)/(thv(i,j,k)*(tkh(i,j,ku)+tkh(i,j,kd) + 0.0001))
           else
              brunt(i,j,k) = - 2*ggr*wthv_sec(i,j,k)/(thv(i,j,k)*(tkh(i,j,ku)+tkh(i,j,kd) + 0.0001))
           endif
        enddo
     enddo
  enddo



!  call eddy_length()   ! Find turbulent mixing length

  do j=1,ny
     do i=1,nx
        brunt_int(i,j,:) = interp_center_to_interface(interp_mode, &
!                                                      brunt(i,j,:), zl(i,j,:), zi(i,j,:),sflux(i,j))
                                                       brunt(i,j,:), zl(i,j,:), zi(i,j,:))
     enddo ! i
  enddo    ! j 

  call eddy_length()   ! Find turbulent mixing length

! test
  Cek_tke = Ce
!  Cek_tke(:,:,1) = Ces*cesfac

  do k=1,nzm      
     ku = k+1
     kd = k
      
!      Cek = Ce * cefac

     if(k == 1) then
        ku = 2
        kd = 2
!test

     elseif(k == nzm) then
        ku = k
        kd = k

     endif

!     if (dis_opt > 0) then
!        do j=1,ny
!           do i=1,nx
!              wrk = (zl(i,j,k)-zi(i,j,1)) / adzl(i,j,1) + 1.5
!              Cek(i,j,k) = 1.0 + 2.0 / max((wrk*wrk - 3.3), 0.5)
!           enddo
!        enddo
!     else
!        if (k == 1) then
!           Cek = ces * cesfac
!        else
!           Cek = ce  * cefac
!        endif
!     endif


     do j=1,ny
        do i=1,nx

           grd = adzl(i,j,k)             !  adzl(k) = zi(k+1)-zi(k)

!         wrk = zl(i,j,k) / grd + 1.5
!         cek = one + 2.0 / (wrk*wrk -3.3)

! TKE boyancy production term. wthv_sec (buoyancy flux) is calculated in
! assumed_pdf(). The value used here is from the previous time step

           a_prod_bu = ggr / thv(i,j,k) * wthv_sec(i,j,k)

! Obtain Brunt-Vaisalla frequency from diagnosed SGS buoyancy flux
! (Note that there's another BV freq. calculation in  eddy_length())
! use epsln instead of 0.0001?

!           buoy_sgs = - (a_prod_bu+a_prod_bu) / (tkh(i,j,ku)+tkh(i,j,kd) + 0.0001) 
            buoy_sgs = brunt(i,j,k)


!           if (buoy_sgs .ne. 0) then
!              
!              wrk = abs((buoy_sgs - brunt(i,j,k)) / buoy_sgs)
!              if (wrk > 0.1) print *, "(buoy_sgs -brunt(i,j,k)) /buoy_sgs", wrk, " buoy_sgs=", buoy_sgs, " brunt=", brunt(i,j,k) 
!
!           else if (brunt(i,j,k) .ne. 0) then 
!
!              print *, "buoy_sgs = 0, brunt(i,j,k))= ", brunt(i,j,k)
!
!           endif

!           if (brunt(i,j,k) .ne. 0) then
!              
!              wrk = abs((brunt(i,j,k)  - brunt_test(i,j,k)) / buoy_sgs)
!              if (wrk > 0.1) print *, "wrk=", wrk, " brunt =", brunt(i,j,k), " brunt_test=", brunt_test(i,j,k) 
!
!           else if (brunt(i,j,k) .ne. 0) then 
!
!              print *, "brunt = 0, brunt_test= ", brunt_test(i,j,k)
!
!           endif



! Compute $c_k$ (variable Cee) for the TKE dissipation term following
! Deardorff, J.W., 1980, Boundary-Layer Meteorol, 18: 495.
! https://doi.org/10.1007/BF00119502

!           if (buoy_sgs <= zero) then
!              smix = grd
!           else
!              smix = min(grd,max(0.1*grd, 0.76*sqrt(tke(i,j,k)/(buoy_sgs+1.e-10))))
!           endif

!          ratio     = smix/grd
!           ratio     = smixt(i,j,k)/grd
  
!           if (smixt(i,j,k) > smix ) &
!           if (Cek(i,j,k) > 20 .or. smixt(i,j,k)== 5000) &
!           print *, "ratio= ", ratio, " smixt(i,j,k)= ", smixt(i,j,k), ' smix=', smix, "grd= ", grd, &
!                "(pt19 + pt51*smixt/grd)= ", (pt19 + pt51*smixt(i,j,k)/grd),  "(pt19 + pt51*smix/grd)=" , (pt19 + pt51*smix/grd), &
!                " diss_smixt= ",  Cek(i,j,k) * (pt19 + pt51*ratio)/ smixt(i,j,k), &
!                " diss_smix= ",  Cek(i,j,k) * (pt19 + pt51*smix/grd)/ smixt(i,j,k), &
!                " diss_new= ",  Cek(i,j,k) / smixt(i,j,k), k,  Cek(i,j,k)

! Make TKE dissipation term pressure dependent to increase damping at TOA
!           Cee       = Cek(i,j) * (pt19 + pt51*ratio) * max(one, sqrt(pcrit/prsl(i,j,k)))
!           Cee       = Cek(i,j) * max(one, sqrt(pcrit/prsl(i,j,k)))
! test
!           Cee       = Cek(i,j,k) * max(one, sqrt(pcrit/prsl(i,j,k)))

           Cee       = Cek_tke(i,j,k) * max(one, sqrt(pcrit/prsl(i,j,k)))
!           Cee       = Cek(i,j,k)*Cek_tke(i,j,k) * max(one, sqrt(pcrit/prsl(i,j,k)))

!           Cee       = Cek(i,j) *  max(one, sqrt(pcrit/prsl(i,j,k)))

!intrp ok?
! TKE shear production term
!           a_prod_sh = half*def2(i,j,k)*(tkh(i,j,ku)*prnum(i,j,ku)   &
!                                       + tkh(i,j,kd)*prnum(i,j,kd))

           a_prod_sh = half*def2(i,j,k)*(tkm(i,j,ku) + tkm(i,j,kd))



! Explicitly integrate TKE equation forward in time
!         a_diss     = Cee/smixt(i,j,k)*tke(i,j,k)**1.5 ! TKE dissipation term
!         tke(i,j,k) = max(zero,tke(i,j,k)+dtn*(max(zero,a_prod_sh+a_prod_bu)-a_diss))

! Semi-implicitly integrate TKE equation forward in time for one time step
! and then iterate for a fully implicit solution. 

           wrk3 =  tke(i,j,k)
           allocate(tke_history(nitr))
           
           wtke = tke(i,j,k)
           wtk2 = wtke
           wrk  = (dtn*Cee) / smixt(i,j,k)
           wrk1 = wtke + dtn*(a_prod_sh+a_prod_bu)

           do itr=1,nitr                        ! iterate for implicit solution
!              wtke   = min(max(min_tke, wtke), max_tke)
              wtke   = max(min_tke, wtke)
              a_diss = wrk*sqrt(wtke)            ! Coefficient in the TKE dissipation term
              wtke   = wrk1 / (one+a_diss)
              wtke   = tkef1*wtke + tkef2*wtk2   ! tkef1+tkef2 = 1.0

 !     if (lprnt .and. i == ipr .and. k<15) write(0,*)' wtke=',wtke,' wtk2=',wtk2,&
!        ' a_diss=',a_diss,' a_prod_sh=',a_prod_sh,' a_prod_bu=',a_prod_bu,&
!        ' wrk1=',wrk1,' kdt=',kdt,' itr=',itr,' k=',k

              wtk2   = wtke


              tke_history(itr)=wtke
     
           enddo

           tke(i,j,k) = min(max(min_tke, wtke), max_tke)

           a_diss = wrk*sqrt(tke(i,j,k)) ! Use updated value of TKE!

                  

!Eq 8 in BK13. Note that there's a typo in BK13, and the correct expression is
! tau = 2*tke/eps

           tscale1 = (dtn+dtn) / a_diss 

           tkesbdiss(i,j,k) = rdtn*a_diss*tke(i,j,k) ! TKE dissipation term, epsilon


           
!           if (wrk3 > 0) then
!           if ( tke(i,j,k)/wrk3 >  5 .and. wrk3 > 10) &
!           if ( isotropy(i,j,k) > 300. .and.  wthl_sec(i,j,k)  > 5.  ) then 
!          if ( isotropy(i,j,k) > 300. .and.  k==2  .and. thl_sec(i,j,k) > 2e2 ) then 
!           if ( tabs(i,j,k) > 325. .or. tke(i,j,k) > 90. .or. tkh(i,j,ku) > 260. ) then         
!           if ( tabs(i,j,k) > 325. .or. tke(i,j,k) > 70. ) then         
! if ( tabs(i,j,k) < 150. ) then         
 if ( tke(i,j,k) > 20. ) then         
!              if (abs(zl(i,j,k) - 286.346836672232) < 1. .and. k .eq. 12) then
                  kdd = k -1
                  if (k==1) kdd=k
                  print *, "smixt=", smixt(i,j,k), " k=", k," zl=", zl(i,j,k)," zl(k-1)=", zl(i,j,kdd)," zl(k+1)=", zl(i,j,k+1), &
                  " l_inf=", l_inf(i,j),  &
                  " tke", tke(i,j,k), " tke(kdd)", tke(i,j,kdd), " tke(ku)", tke(i,j,ku),  &
                  " brunt=", brunt(i,j,k), " brunt(k-1)=", brunt(i,j,kdd),  &
                  " brunt(k+1)=", brunt(i,j,k+1), &
                  " brunt_test=", brunt_test(i,j,k), " brunt_test(k-1)=", brunt_test(i,j,kdd),&
                       " brunt_test(k+1)=", brunt_test(i,j,k+1), &
                  " qcl(i,j,k)=", qcl(i,j,k), " qci(i,j,k)=", qci(i,j,k), &
                  " isotropy=", isotropy(i,j,k) ," tkesbdiss=", tkesbdiss(i,j,k) ,  &
                  " isotropy(kdd)=", isotropy(i,j,kdd), " isotropy(ku)=", isotropy(i,j,ku), " dtn=", dtn,&
                  " a_prod_sh=", a_prod_sh, " a_prod_bu=", a_prod_bu, " Cee=", Cee, &
                  " Cek(i,j,k)=", Cek(i,j,k), &
                  " tke_orig=", wrk3, " wrk=", wrk, " wrk1=", wrk1," a_diss=", a_diss, &
                  " def2(i,j,k)=", def2(i,j,k)," tkh(i,j,ku)=", tkh(i,j,ku), "prnum(i,j,ku)=" ,prnum(i,j,ku),  &
                  " tkh(i,j,kd)=", tkh(i,j,kd), "prnum(i,j,kd)=" ,prnum(i,j,kd), &
                  " tke_history=", tke_history, &
                  " tscale1=", tscale1, " thv(i,j,k)=", thv(i,j,k), " sflux(i,j)=", sflux(i,j), &
                  " wthv_sec(i,j,k)=", wthv_sec(i,j,k), " wthv_sec(i,j,k-1)=", wthv_sec(i,j,kdd),  &
                  " wthv_sec(i,j,k+1)=", wthv_sec(i,j,k+1),  &
                  " thl_sec(i,j,k)=", thl_sec(i,j,k) , &
                  " Cek_tke(i,j,k)=", Cek_tke(i,j,k), &
                  " Cek_tke(i,j,k) real=", Cek_tke(i,j,k)* max(one, sqrt(pcrit/prsl(i,j,k))), &
                  " u(i,j,k-1)=", u(i,j,kdd), " u(i,j,k)=", u(i,j,k), " u(i,j,k+1)=", u(i,j,k+1), &
                  " v(i,j,k-1)=", u(i,j,kdd)," v(i,j,k)=", v(i,j,k)," v(i,j,k+1)=", v(i,j,k+1) , &
                  " tabs(i,j,k-1)=", tabs(i,j,kdd)," tabs(i,j,k)=", tabs(i,j,k)," tabs(i,j,k+1)=", tabs(i,j,k+1), &
                  "  rdzw=", one /  (zl(i,j,ku) - zl(i,j,kdd)), "  wrku=" , u(i,j,ku)-u(i,j,kdd), &
                  " wrkv=", v(i,j,ku)-v(i,j,kdd), " prod=" ,  ((u(i,j,ku)-u(i,j,kdd))**2+ (v(i,j,ku)-v(i,j,kdd))**2)
 !                 endif
           endif 

           
               
            deallocate(tke_history)

! Calculate "return-to-isotropy" eddy dissipation time scale, see Eq. 8 in BK13

           if (buoy_sgs <= zero) then
              isotropy(i,j,k) = min(max_eddy_dissipation_time_scale,tscale1)
           else
              isotropy(i,j,k) = min(max_eddy_dissipation_time_scale,          &
                                   tscale1/(one+lambda*buoy_sgs*tscale1*tscale1))
           endif


! TKE budget terms

!         tkesbdiss(i,j,k)       = a_diss
!         tkesbshear(i,j,k)      = a_prod_sh
!         tkesbbuoy(i,j,k)       = a_prod_bu
!         tkesbbuoy_debug(i,j,k) = a_prod_bu_debug
!         tkebuoy_sgs(i,j,k)     = buoy_sgs

! Save diagnostic output.

            if (shoc_diag) then

               Diag%SHOC_eddy_length(i,k) =  smixt(i,j,k)

               Diag%SHOC_TKE_boyancy_prod(i,k) =  a_prod_bu 
               Diag%SHOC_TKE_shear_prod(i,k)   =  a_prod_sh
               Diag%SHOC_TKE_dissipation(i,k)  =  tkesbdiss(i,j,k)
               Diag%SHOC_TKE_damping_coef(i,k) =  Cee

               Diag%SHOC_ret2iso_tscale(i,k)   = isotropy(i,j,k)
               Diag%SHOC_tkh(i,k)              = (tkh(i,j,ku)+tkh(i,j,kd))*half
               Diag%SHOC_Pr(i,k)               = (prnum(i,j,ku)+prnum(i,j,kd))*half

            endif


        enddo ! i loop
     enddo   ! j loop
  enddo     ! k


! Calculate eddy diffusivity (without background diffusion)

  tke_int = 0.
  do j=1,ny
     do i=1,nx
       isotropy_int(i,j,:) = interp_center_to_interface(interp_mode, &
                                                    isotropy(i,j,:), zl(i,j,:), zi(i,j,:)) 
!       tke_int(i,j,:)      = interp_center_to_interface(interp_mode, tke(i,j,:),      zl(i,j,:), zi(i,j,:))  
       tkh_out1(i,j,:)          = max(0.,min(tkhmax,interp_center_to_interface(interp_mode,  &
                                 ck*tke(i,j,:)*isotropy(i,j,:),zl(i,j,:), zi(i,j,:),0.)))
       tkh_out(i,j,:)          = max(0.,min(tkhmax,interp_center_to_interface(interp_mode,  &
                                 ck*tke(i,j,:)*isotropy(i,j,:),zl(i,j,:), zi(i,j,:))))
     enddo ! i      
  enddo   ! j  


! Apply background diffusivity to diffusivity values computed by SHOC
  do k=1,nzm-1
     ku = k + 1
!  do k=1,nzm
!     ku = k 

     do j=1,ny
        do i=1,nx

!           tkh(i,j,ku)= tkh_out(i,j,ku)  + xkzo(i,j,k)
!           tkm(i,j,ku)= tkh_out(i,j,ku)*prnum(i,j,ku) + xkzmo(i,j,k)
           tkh(i,j,ku)= max(tkh_out(i,j,ku), xkzo(i,j,k))
           tkm(i,j,ku)= max(tkh_out(i,j,ku)*prnum(i,j,ku), xkzmo(i,j,k))

 if (  tke(i,j,k) > 20. ) then         
! if ( tabs(i,j,k) > 325. .or. tke(i,j,k) > 70. ) then         
! if ( tabs(i,j,k) < 150. .and. k>1) then         
                  kdd = k -1
                  if (k==1) kdd=k
                  print *, "after tkh smixt=", smixt(i,j,k), " k=", k," zl=", zl(i,j,k)," zl(k-1)=", zl(i,j,kdd)," zl(k+1)=", zl(i,j,k+1), &
                  " l_inf=", l_inf(i,j),  &
                  " tke", tke(i,j,k), " tke(kdd)", tke(i,j,kdd), " tke(ku)", tke(i,j,k+1),  &
                  " brunt=", brunt(i,j,k), " brunt(k-1)=", brunt(i,j,kdd),  &
                  " brunt(k+1)=", brunt(i,j,k+1), &
                  " brunt_test=", brunt_test(i,j,k), " brunt_test(k-1)=", brunt_test(i,j,kdd),&
                       " brunt_test(k+1)=", brunt_test(i,j,k+1), &
                  " qcl(i,j,k)=", qcl(i,j,k), " qci(i,j,k)=", qci(i,j,k), &
                  " isotropy=", isotropy(i,j,k) , &
                  " isotropy(kdd)=", isotropy(i,j,kdd), " isotropy(ku)=", isotropy(i,j,k+1), " dtn=", dtn,&
                  " tkesbdiss(k)=", tkesbdiss(i,j,k) ,  " tkesbdiss(k-1)=", tkesbdiss(i,j,kdd) , &
                  " tkesbdiss(k+1)=", tkesbdiss(i,j,k+1) , &
!                  " a_prod_sh=", a_prod_sh, " a_prod_bu=", a_prod_bu, " Cee=", Cee, 
                  " Cek(i,j,k)=", Cek(i,j,k), &
!                  " tke_orig=", wrk3, " wrk=", wrk, " wrk1=", wrk1," a_diss=", a_diss, &
                  " def2(i,j,k)=", def2(i,j,k), &
                  " tkh(i,j,k)=", tkh(i,j,k), "prnum(i,j,k)=" ,prnum(i,j,k), &
                  " tkh(i,j,ku)=", tkh(i,j,k+1), "prnum(i,j,ku)=" ,prnum(i,j,k+1),  &
                  " tkh(i,j,kd)=", tkh(i,j,kdd), "prnum(i,j,kd)=" ,prnum(i,j,kdd), &
!                  " tke_history=", tke_history, &
!                  " tscale1=", tscale1, " thv(i,j,k)=",
                  thv(i,j,k), " sflux(i,j)=", sflux(i,j), &
                  " wthv_sec(i,j,k)=", wthv_sec(i,j,k), " wthv_sec(i,j,k-1)=", wthv_sec(i,j,kdd),  &
                  " wthv_sec(i,j,k+1)=", wthv_sec(i,j,k+1),  &
                  " thl_sec(i,j,k)=", thl_sec(i,j,k) , &
                  " Cek_tke(i,j,k)=", Cek_tke(i,j,k), &
                  " Cek_tke(i,j,k) real=", Cek_tke(i,j,k)* max(one, sqrt(pcrit/prsl(i,j,k))), &
                  " u(i,j,k-1)=", u(i,j,kdd), " u(i,j,k)=", u(i,j,k), " u(i,j,k+1)=", u(i,j,k+1), &
                  " v(i,j,k-1)=", u(i,j,kdd)," v(i,j,k)=", v(i,j,k)," v(i,j,k+1)=", v(i,j,k+1) , &
                  " tabs(i,j,k-1)=", tabs(i,j,kdd)," tabs(i,j,k)=", tabs(i,j,k)," tabs(i,j,k+1)=", tabs(i,j,k+1), &
                  "  rdzw=", one /  (zl(i,j,k+1) - zl(i,j,kdd)), "  wrku=" , u(i,j,k+1)-u(i,j,kdd), &
                  " wrkv=", v(i,j,k+1)-v(i,j,kdd), " prod=" ,  ((u(i,j,k+1)-u(i,j,kdd))**2+ (v(i,j,k+1)-v(i,j,kdd))**2)
 !                 endif
           endif 

        
        enddo
     enddo
  enddo


  tkh(:,1,1)= tkh(:,1,2)
  tkm(:,1,1)= tkm(:,1,2)


!intrp
! debugging code
 wrk = half * ck
 do k=2,nzm
    km1 = k - 1
     do j=1,ny
        do i=1,nx
           tke_int(i,j,k) = min(tkhmax, wrk * (isotropy(i,j,k)  * tke(i,j,k)    &
                + isotropy(i,j,km1) * tke(i,j,km1))) 
!           tkh(i,j,k) = min(tkhmax, ck * isotropy_int(i,j,k)*tke_int(i,j,k))

        enddo ! i
     enddo   ! j
  enddo     ! k

  tke_int(:,1,1)=tke_int(:,1,2)
 




  
! Diagnose second order moments of the subgrid PDF following
! Redelsperger J.L., and G. Sommeria, 1986, JAS, 43, 2619-2635 
! https://doi.org/10.1175/1520-0469(1986)043<2619:TDSOAC>2.0.CO;2
! sans the use of  stability weighting functions 
    
! Second moment of vertical velocity defined on layer centers
! Eq 6 in BK13 

  do k=1,nzm
    ku = k+1
    kd = k-1
    ka = ku
    kb = k
    if (k == 1) then
      kd = k
      kb = ka
    elseif (k == nzm) then
      ku = k
      ka = kb
    endif
    do j=1,ny
      do i=1,nx
        if (tke(i,j,k) > zero) then
!intrp ok
!          wrk  = half*(tkh(i,j,ka)+tkh(i,j,kb))*(w(i,j,ku) - w(i,j,kd)) & 
          wrk  = half*(tkm(i,j,ka)+tkm(i,j,kb))*(w(i,j,ku) - w(i,j,kd)) & 
               * sqrt(tke(i,j,k)) / (zl(i,j,ku) - zl(i,j,kd))
          w_sec(i,j,k) = max(twoby3 * tke(i,j,k) - twoby15 * wrk, zero)
        else
          w_sec(i,j,k) = zero
        endif
      enddo ! i
    enddo   ! j 
  enddo     ! k

  do j=1,ny
     do i=1,nx
        w_sec_int(i,j,:) = interp_center_to_interface(interp_mode, & 
                                                      w_sec(i,j,:), zl(i,j,:), zi(i,j,:),0.)
     enddo ! i  
  enddo    ! j   
    
  do k=2,nzm
       
    km1 = k - 1
    if (k==2) then 
       ku=3
       kd=2
    elseif (k==nzm) then

       ku=nzm
       kd=nzm-1

    else

       ku=k
       kd=k-1
    endif

    do j=1,ny
      do i=1,nx

! Define diagnostic second moments on layer interfaces (except w_sec). Thermodynamic variables, 
! TKE, and  "return-to-isotropy" time scale are given on the layer centers, so 
! interpolate them to the interfaces. 
! Diffusion coefficient in GFS is already on the interface. 

        wrk1 = one / adzi(i,j,k)        ! adzi(k) = (zl(k)-zl(km1))
!       wrk3 = max(tkh(i,j,k),pt01) * wrk1
!        wrk3 = max(tkh(i,j,k),epsln) * wrk1
        wrk3 = tkh(i,j,k) * wrk1
!intrp 
!       sm   = half*(isotropy(i,j,k)+isotropy(i,j,km1))*wrk1*wrk3 ! Tau*Kh/dz^2 
        sm   = isotropy_int(i,j,k)*wrk1*wrk3 ! Tau*Kh/dz^2 
        

! SGS vertical flux liquid/ice water static energy. Eq 1 in BK13
             
        wrk1            = hl(i,j,k) - hl(i,j,km1)
!               No rain, snow or graupel in pdf (Annig, 08/29/2018)
!                       + (qpl(i,k) - qpl(i,km1)) * fac_cond &
!                       + (qpi(i,k) - qpi(i,km1)) * fac_sub
        wthl_sec_d(i,j,k) = - wrk3 * wrk1
!        wthl_sec(i,j,k) = - wrk3 * wrk1

!        if (wthl_sec(i,j,k) > 50 ) then
        if ( tabs(i,j,k) < 150. ) then
           kdd = k -1
           if (k==1) kdd=k

           print *, " wrk1=",  one / adzi(i,j,k) , " wrk3=", wrk3, " tkh(i,j,k)=", tkh(i,j,k), &
                " wrk1 =", wrk1, " hl(i,j,k)=", hl(i,j,k)," hl(i,j,km1)=", hl(i,j,kdd), " k=", k, &
                " km1=", km1, " i=", i, &
                 " wthl_sec(i,j,k)=", wthl_sec(i,j,k), " tabs(i,j,k)=", tabs(i,j,k), " gamaz(i,j,k)=", &
                 gamaz(i,j,k), "  fac_cond*(qcl(i,j,k)+qpl(i,j,k))=",  fac_cond*(qcl(i,j,k)+qpl(i,j,k)), &
                 " fac_sub *(qci(i,j,k)+qpi(i,j,k))=", fac_sub *(qci(i,j,k)+qpi(i,j,k)), &
                 " xkzo(i,j,k) =", xkzo(i,j,k),"  tkh_out(i,j,:)=",  tkh_out(i,j,k), " tke(i,j,:)=", tke(i,j,k), &
                 " isotropy(i,j,k)=", isotropy(i,j,k)," isotropy(i,j,km1)=", isotropy(i,j,kdd), &
                  "smixt=", smixt(i,j,k), " k=", k," zl=", zl(i,j,k), &
                  " l_inf=", l_inf(i,j),  &
                  " tke", tke(i,j,k), " brunt=", brunt(i,j,k), " brunt(k-1)=", brunt(i,j,kdd),  &
                  " brunt_test=", brunt_test(i,j,k), &
                  " qcl(i,j,k)=", qcl(i,j,k), " qci(i,j,k)=", qci(i,j,k), &
                  " isotropy=", isotropy(i,j,k) ," tkesbdiss=", tkesbdiss(i,j,k) ,  &
                  " a_prod_sh=", a_prod_sh, " a_prod_bu=", a_prod_bu, " Cee=", Cee, " Cek(i,j,k)=", Cek(i,j,k), &
                  " max(one, sqrt(pcrit/prsl(i,j,k)))=", max(one, sqrt(pcrit/prsl(i,j,k))), &
                  " tke_orig=", wrk3, " wrk=", wrk, " wrk1=", wrk1," a_diss=", a_diss, &
                  " def2(i,j,k)=", def2(i,j,k)," tkh(i,j,ku)=", tkh(i,j,ku), "prnum(i,j,ku)=" ,prnum(i,j,ku),  &
                  " tkh(i,j,ku)=", tkh(i,j,kd), "prnum(i,j,kd)=" ,prnum(i,j,kd), &
!                  " tke_history=", tke_history
                  " tscale1=", tscale1, " thv(i,j,k)=", thv(i,j,k), " sflux(i,j)=", sflux(i,j), &
                  " wthv_sec(i,j,k)=", wthv_sec(i,j,k), "  tkh_out(i,j,k)=",  tkh_out(i,j,k),  &
                  "  tkh_out1(i,j,k)=",  tkh_out1(i,j,k),  &
                  " tke_int(i,j,k)=", tke_int(i,j,k),  " xkzo(i,j,k)=", xkzo(i,j,k)

        endif

! SGS vertical flux of total water. Eq 2 in BK13

        wrk2           = total_water(i,j,k) - total_water(i,j,km1)
        wqw_sec_d(i,j,k) = - wrk3 * wrk2
!        wqw_sec(i,j,k) = - wrk3 * wrk2

!        if (.false.) then
! Calculate diffusive flux of MSE from diffusive fluxes of temperature, moisture, 
! and condensate computed by the vertical diffusion solver
!        wthl_sec(i,j,k) =  wthl_sec(i,j,k-1) - adzl(i,j,k-1)*(dtabsdt(i,j,k-1)     -  &
!                          fac_cond*(dqcldt(i,j,k-1) + dqpldt(i,j,k-1))            -  &
!                          fac_sub* (dqcidt(i,j,k-1) + dqpidt(i,j,k-1)))            -  &
!                          gocp*(tkh(i,j,ku) - tkh(i,j,kd))


             if (.false.) then
!        if ( wthl_sec(i,j,k) > 1. ) then 
!           if (abs((wthl_sec_d(i,j,k)- wthl_sec(i,j,k))/ wthl_sec(i,j,k)) < 10) &
                
                print *, " wthl_sec=", wthl_sec(i,j,k), " wthl_sec_d(k)=", wthl_sec_d(i,j,k), &
                " wthl_sec(k-1)=", wthl_sec(i,j,k-1), &
                " k=", k, "adzl=", adzl(i,j,k-1), " dtabsdt=", dtabsdt(i,j,k-1), &
                " fac_cond*dqcldt=", fac_cond* dqcldt(i,j,k-1), &
                " fac_cond*dqpldt=", fac_cond* dqpldt(i,j,k-1), &
                " fac_sub*dqcidt=", fac_sub* dqcidt(i,j,k-1) , &
                " fac_sub*dqpidt=", fac_sub* dqpidt(i,j,k-1) , &
                " gocp*tkh(ku)=",  gocp*tkh(i,j,ku), &
                " gocp*tkh(kd)=",gocp*tkh(i,j,kd), &
                " tkh(ku)=",  tkh(i,j,ku), &
                " tkh(kd)=", tkh(i,j,kd), &
                " xkzo(ku)=",  xkzo(i,j,ku), &
                " xkzo(kd)=", xkzo(i,j,kd), &
                " whole thing=", adzl(i,j,k-1)*(dtabsdt(i,j,k-1)     -  &
                fac_cond*(dqcldt(i,j,k-1) + dqpldt(i,j,k-1))            -  &
                fac_sub* (dqcidt(i,j,k-1) + dqpidt(i,j,k-1)))           -  &
                gocp*(tkh(i,j,ku) - tkh(i,j,kd)), &
                " dqwvdt=", dqwvdt(i,j,k-1), " rho=", prsl(i,j,k)/(rgas*tabs(i,j,k)*(1+0.622*total_water(i,j,k))), &
                " prsl(i,j,k)=", prsl(i,j,k)," tabs(i,j,k)=", tabs(i,j,k)
           
        endif

        
                   

! Calculate diffusive flux of total water from diffusive fluxes of  moisture
! and condensate computed by the vertical diffusion solver                            
        wqw_sec(i,j,k) =  wqw_sec(i,j,k-1)  - adzl(i,j,k-1)*(dqwvdt(i,j,k-1)     +  &
                                           dqcldt(i,j,k-1) + dqcidt(i,j,k-1))


! Incorrect
!       wqp_sec (i,j,k) = wqp_sec(i,j,k-1)  - adzl(i,j,k-1)*(dqpldt(i,j,k-1) + dqpidt(i,j,k-1))

!     endif
! Diagnostic (SHOC v1) treatment of variances of total water mixing ratio and
! liquid/ice water static energy. Prognostic code is below. 

        if (shoc_version == 1) then
! Second moment of liquid/ice water static energy. Eq 4 in BK13  

           thl_sec(i,j,k) = thl2tune * sm * wrk1 * wrk1

! Second moment of total water mixing ratio.  Eq 3 in BK13
             
           qw_sec(i,j,k)  = qw2tune * sm * wrk2 * wrk2

        endif
             
! Covariance of total water mixing ratio and liquid/ice water static energy.
! Eq 5 in BK13
             
        qwthl_sec(i,j,k) = qwthl2tune * sm * wrk1 * wrk2

      enddo ! i  loop
    enddo   ! j  loop

         if (.false.) then
!         if (k==nzm .and. any(smixt(i,j,1:5)>7000)) then
!        if ( wthl_sec(i,j,k) > 1. ) then 
!           if (abs((wthl_sec_d(i,j,k)- wthl_sec(i,j,k))/ wthl_sec(i,j,k)) < 10) &
                
                print *, " wthl_sec=", wthl_sec(i,j,1:nzm), new_line ('a'), &!
                     " wthl_sec_d(k)=", wthl_sec_d(i,j,1:nzm),  new_line ('a'),&
!                " wthl_sec(k-1)=", wthl_sec(i,j,k-1), &
                " k=", k, "adzl=", adzl(i,j,1:nzm), new_line ('a')," dtabsdt=", dtabsdt(i,j,1:nzm),  new_line ('a'),&
                " fac_cond*dqcldt=", fac_cond* dqcldt(i,j,1:nzm),  new_line ('a'),&
                " fac_cond*dqpldt=", fac_cond* dqpldt(i,j,1:nzm),  new_line ('a'),&
                " fac_sub*dqcidt=", fac_sub* dqcidt(i,j,1:nzm) ,  new_line ('a'),&
                " fac_sub*dqpidt=", fac_sub* dqpidt(i,j,1:nzm) ,  new_line ('a'),&
                " gocp*tkh(ku)=",  gocp*tkh(i,j,1:nzm),  new_line ('a'),&
!                " gocp*tkh(kd)=",gocp*tkh(i,j,kd), &
                " tkh(ku)=",  tkh(i,j,ku),  new_line ('a'),&
 !               " tkh(kd)=", tkh(i,j,kd), &
                " xkzo(ku)=",  xkzo(i,j,1:nzm),  new_line ('a'),&
 !               " xkzo(kd)=", xkzo(i,j,kd), &
!                " whole thing=", adzl(i,j,k-1)*(dtabsdt(i,j,k-1)     -  &
!                fac_cond*(dqcldt(i,j,k-1) + dqpldt(i,j,k-1))            -  &
!                fac_sub* (dqcidt(i,j,k-1) + dqpidt(i,j,k-1)))           -  &
!                gocp*(tkh(i,j,ku) - tkh(i,j,kd)), &
                " dqwvdt=", dqwvdt(i,j,1:nzm),  new_line ('a'),&
                " (smixt(i,j,:)=", smixt(i,j,1:nzm)
!                " rho=", prsl(i,j,k)/(rgas*tabs(i,j,k)*(1+0.622*total_water(i,j,k))), &
!                " prsl(i,j,k)=", prsl(i,j,k)," tabs(i,j,k)=", tabs(i,j,k)
           
        endif


  enddo     ! k  loop


! Values of second moments at the surface/lowest model interface need to be specified.
! Fluxes of total water and MSE are given by the boundary conditions.
! Variances of MSE and total water in SHOC v1 are defined on interfaces, and
! their value on the lowest interface are set to be identical to their values
! on the interface above. 
! Variances of MSE and total water in SHOC v2 are defined on the layer centers
! as they are prognostic variables, and are calculated explicitly.
  do j=1,ny
    do i=1,nx
! Fluxes of MSE and total water at the surface are taken from the host model
!     wthl_sec(i,j,1)  = wthl_sec(i,j,2)
!     wqw_sec(i,j,1)   = wqw_sec(i,j,2)
      qwthl_sec(i,j,1) = qwthl_sec(i,j,2)
      if (shoc_version == 1) then
         thl_sec(i,j,1)   = thl_sec(i,j,2)
         qw_sec(i,j,1)    = qw_sec(i,j,2)
      endif
   enddo
  enddo

  if (shoc_version == 2) then

! Semi-implicit prognostic equations for variances of total water mixing ratio
! and  liquid/ice water static energy for SHOC v2. Note that prognostic variances
! are defined on layer centers, as opposed to the diagnostic ones.  
! See Eqs. 13,14, and 24b in Golaz, J., V. E. Larson, and W. R. Cotton, 2002: 
! A PDF-Based Model for Boundary Layer Clouds. Part I: Method and Model Description. 
! J. .Atmos. Sci. , 59 , 3540–3551 
! https://doi.org/10.1175/1520-0469(2002)059<3540:APBMFB>2.0.CO;2
! (referenced as GLC02 further in the comments)  

    dts = dtn /substep 

    do k=1,nzm

       km1 = k+1
       if (k == nzm) km1=k


       do j=1,ny
          do i=1,nx

! Make dissipation terms pressure dependent to increase damping at TOA
             Cvv = Cv * max(one, sqrt(pcrit/prsl(i,j,k)))

             do itr = 1, substep 

! Second moment of liquid/ice water static energy.

! Divergence of subgrid vertical flux of MSE variance 
! For now allowing this quantity to be calulated by the diffusion solver
! for consistency with the TKE diffusion. The alternative is to calculate it
! from the PDF. 
               wrk1 =  0 ! dwthl2dz(i,j,k)
! Subgrid flux of MSE interpolated from layer interfaces to layer center
               wrk2 =  (wthl_sec(i,j,k)+wthl_sec(i,j,km1))*half
!               wrk2 = (wthl_sec(i,j,k)+wthl_sec(i,j,km1))
! Vertical gradient of MSE
               wrk3 =  dthldz(i,j,k)

! test
!               wrk  =  Cvv*dts/isotropy(i,j,k)
               wrk  =  Cvv*Cek(i,j,k)*dts/isotropy(i,j,k)


               if ( prog_var_solver_type == 1 ) then

! First order implicit backwards Euler method 

!                  thl_sec(i,j,k) = (thl_sec(i,j,k) - dts* &
!                                 (wrk1 +wrk2*wrk3 - tune_varthl*detrained_varthl(i,j,k)) &
!                                 )/(1.+Cv*dts/isotropy(i,j,k))

                  thl_sec(i,j,k) = (thl_sec(i,j,k) - dts* &
                                    (wrk1 + 2*wrk2*wrk3 - tune_varthl*detrained_varthl(i,j,k)) &
                                   )/(one + wrk)

               else
! Second order implicit backwards Euler method 
                  
                  thl_sec(i,j,k) = (thl_sec(i,j,k)*(one - half*wrk) - dts* &
                                    (wrk1 + 2*wrk2*wrk3 - tune_varthl*detrained_varthl(i,j,k)) &
                                   )/(one + half*wrk)

               endif

               if ( thl_sec(i,j,k) < 0.) thl_sec(i,j,k) = 0.

!               if (thl_sec(i,j,k) > 1e2) print *, "debug thl, wrk1, wrk2, &
!                    wrk3,2*wrk2*wrk3,1 + Cv*dts/isotropy(i,j,k), dts, tune_varthl*detrained_varthl(i,j,k),  &
!                    wrk1 + 2*wrk2*wrk3 - tune_varthl*detrained_varthl(i,j,k),  &
!                    (wrk1 + 2*wrk2*wrk3 - tune_varthl*detrained_varthl(i,j,k))*dts, thl_sec(i,j,k), k
!               if (thl_sec(i,j,k) > 2e2) &
!                    if ( tabs(i,j,k) > 325. .or. thl_sec(i,j,k) > 2e2 ) then 
                    if ( tabs(i,j,k) < 150.) then 
                       kdd = k -1
                       if (k==1) kdd=k

                    print *, "debug thl", wrk1, wrk2, &
                    wrk3,2*wrk2*wrk3,1 + half*wrk, 1- half*wrk, dts, tune_varthl*detrained_varthl(i,j,k),  &
                    wrk1 + 2*wrk2*wrk3 - tune_varthl*detrained_varthl(i,j,k),  &
                    (wrk1 + 2*wrk2*wrk3 - tune_varthl*detrained_varthl(i,j,k))*dts, "thl_sec", thl_sec(i,j,k), k, "wth_sec", wthl_sec(i,j,k), wthl_sec(i,j,k+1),  "hl", hl(i,j,kdd),  hl(i,j,k),  hl(i,j,k+1), one / adzi(i,j,k), "tkh" , max(tkh(i,j,kdd),epsln), tkh(i,j,k),  tkh(i,j,k+1), max(tkh(i,j,k),epsln)/ adzi(i,j,k), "tkh_lin", tke_int(i,j,kdd),  tke_int(i,j,k),  tke_int(i,j,k+1), "tke", tke(i,j,kdd), tke(i,j,k), tke(i,j,k+1), "isotropy", isotropy(i,j,kdd),  isotropy(i,j,k),  isotropy(i,j,k+1), "isotropy_int",  isotropy_int(i,j,kdd), isotropy_int(i,j,k),  isotropy_int(i,j,k+1),  "zi", zi(i,j,kdd),  zi(i,j,k),  zi(i,j,k+1), zi(i,j,k+2), "zl", zl(i,j,kdd),  zl(i,j,k),  zl(i,j,k+1), zl(i,j,k+2), "wthl2", wthl2(i,j,kdd), wthl2(i,j,k), wthl2(i,j,k+1),  wthl2(i,j,k+2), " dwthl2dz",  dwthl2dz(i,j,kdd),  dwthl2dz(i,j,k),  dwthl2dz(i,j,k+1), "w_sec", w_sec(i,j,kdd), w_sec(i,j,k), w_sec(i,j,k+1),  w_sec(i,j,k+2),  w_sec(i,j,k+3), "w_sec_int", w_sec_int(i,j,kdd), w_sec_int(i,j,k), w_sec_int(i,j,k+1),  w_sec_int(i,j,k+2), "  tkh_out(i,j,k)=", &
 tkh_out(i,j,k), "  tkh_out(i,j,k+1)=",  tkh_out(i,j,k+1), "  tkh_out(i,j,k+2)=",  tkh_out(i,j,k+2), &
" tke_int(i,j,k)=", tke_int(i,j,k), " tke_int(i,j,k+1)=", tke_int(i,j,k+1), " tke_int(i,j,k+2)=", tke_int(i,j,k+2),& 
" xkzo(i,j,k)=", xkzo(i,j,k), " xkzo(i,j,k+1)=", xkzo(i,j,k+1), " xkzo(i,j,k+2)=", xkzo(i,j,k+2), &

"  tkh_out1(i,j,k)=",  &
 tkh_out1(i,j,k), "  tkh_out1(i,j,k+1)=",  tkh_out1(i,j,k+1), "  tkh_out1(i,j,k+2)=",  tkh_out1(i,j,k+2)


                    endif


! Second moment of total water mixing ratio.

! Divergence of subgrid vertical flux of total water variance 
! For now allowing this quantity to be calulated by the diffusion solver                                                                             ! for consistency with the TKE diffusion. The alternative is to calculate it                                                                        ! from the PDF.    
               wrk1 =  0 ! dwqw2dz(i,j,k)
! Subgrid flux of total water interpolated from layer interfaces to layer center
               wrk2 =  (wqw_sec(i,j,k) + wqw_sec(i,j,km1))*half
!               wrk2 =  (wqw_sec(i,j,k) + wqw_sec(i,j,km1))
! Vertical gradient of total water
               wrk3 = dqwdz(i,j,k)

               if ( prog_var_solver_type == 1 ) then

! First order implicit backwards Euler method 

!                  qw_sec(i,j,k) = (qw_sec(i,j,k) - dts*   &
!                                (wrk1 + wrk2*wrk3 - tune_varqt*detrained_varqt(i,j,k)) &
!                               )/(1.+Cv*dts/isotropy(i,j,k))

                  qw_sec(i,j,k) = (qw_sec(i,j,k) - dts*   &
                                   (wrk1 + 2*wrk2*wrk3 - tune_varqt*detrained_varqt(i,j,k)) &
                                  )/(one + wrk)

               else

! Second order implicit backwards Euler method 

                  qw_sec(i,j,k) = (qw_sec(i,j,k)*(one - half*wrk)  - dts*   &
                                   (wrk1 + 2*wrk2*wrk3 - tune_varqt*detrained_varqt(i,j,k)) &
                                  )/(one + half*wrk)

               endif


               if ( qw_sec(i,j,k) < 0.) qw_sec(i,j,k) = 0.

             enddo ! itr


!                    if ( tabs(i,j,k) > 325.  ) then 
!                       kdd = k -1
!                       if (k==1) kdd=k

!                    print *, "debug thl", wrk1, wrk2, &
!                    wrk3,2*wrk2*wrk3,1 + half*wrk, 1- half*wrk, dts, tune_varthl*detrained_varthl(i,j,k),  &
!                    wrk1 + 2*wrk2*wrk3 - tune_varthl*detrained_varthl(i,j,k),  &
!                    (wrk1 + 2*wrk2*wrk3 - tune_varthl*detrained_varthl(i,j,k))*dts, "thl_sec", thl_sec(i,j,k), k, "wth_sec", wthl_sec(i,j,k), wthl_sec(i,j,k+1),  "hl", hl(i,j,kdd),  hl(i,j,k),  hl(i,j,k+1), one / adzi(i,j,k), "tkh" , max(tkh(i,j,kdd),epsln), tkh(i,j,k),  tkh(i,j,k+1), max(tkh(i,j,k),epsln)/ adzi(i,j,k), "tkh_lin", tke_int(i,j,kdd),  tke_int(i,j,k),  tke_int(i,j,k+1), "tke", tke(i,j,kdd), tke(i,j,k), tke(i,j,k+1), "isotropy", isotropy(i,j,kdd),  isotropy(i,j,k),  isotropy(i,j,k+1), "isotropy_int",  isotropy_int(i,j,kdd), isotropy_int(i,j,k),  isotropy_int(i,j,k+1),  "zi", zi(i,j,kdd),  zi(i,j,k),  zi(i,j,k+1), zi(i,j,k+2), "zl", zl(i,j,kdd),  zl(i,j,k),  zl(i,j,k+1), zl(i,j,k+2), "wthl2", wthl2(i,j,kdd), wthl2(i,j,k), wthl2(i,j,k+1),  wthl2(i,j,k+2), " dwthl2dz",  dwthl2dz(i,j,kdd),  dwthl2dz(i,j,k),  dwthl2dz(i,j,k+1), "w_sec", w_sec(i,j,kdd), w_sec(i,j,k), w_sec(i,j,k+1),  w_sec(i,j,k+2),  w_sec(i,j,k+3), "w_sec_int", w_sec_int(i,j,kdd), w_sec_int(i,j,k), w_sec_int(i,j,k+1),  w_sec_int(i,j,k+2), "  tkh_out(i,j,k)=", &
! tkh_out(i,j,k), "  tkh_out(i,j,k+1)=",  tkh_out(i,j,k+1), "  tkh_out(i,j,k+2)=",  tkh_out(i,j,k+2), &
!" tke_int(i,j,k)=", tke_int(i,j,k), " tke_int(i,j,k+1)=", tke_int(i,j,k+1), " tke_int(i,j,k+2)=", tke_int(i,j,k+2),& 
!" xkzo(i,j,k)=", xkzo(i,j,k), " xkzo(i,j,k+1)=", xkzo(i,j,k+1), " xkzo(i,j,k+2)=", xkzo(i,j,k+2), &

!"  tkh_out1(i,j,k)=",  &
! tkh_out1(i,j,k), "  tkh_out1(i,j,k+1)=",  tkh_out1(i,j,k+1), "  tkh_out1(i,j,k+2)=",  tkh_out1(i,j,k+2)


!                    endif


            if (shoc_diag) then

!               Diag%SHOC_qt_qt_shear_prod(i,k)       = 
!               Diag%SHOC_qt_qt_detrainment_prod(i,k) = 
!               Diag%SHOC_qt_qt_dissipation(i,k)      = 
               
!               Diag%SHOC_hl_hl_shear_prod(i,k)       =
!               Diag%SHOC_hl_hl_detrainment_prod(i,k) =
!               Diag%SHOC_hl_hl_dissipation(i,k)      = 
               

            endif



!             if (qw_sec(i,j,k) > 1e-2) print *, "debug qw", wrk1, wrk2, wrk3, 2*wrk2*wrk3, &
!                  1 + Cv*dts/isotropy(i,j,k), dts, tune_varqt*detrained_varqt(i,j,k), &
!                  wrk1 + 2*wrk2*wrk3 - tune_varqt*detrained_varqt(i,j,k), &
!                  (wrk1 + 2*wrk2*wrk3 - tune_varqt*detrained_varqt(i,j,k))*dts, &
!                  qw_sec(i,j,k), k


          enddo    ! i  
       enddo       ! j  
    enddo          ! k  

 endif

! Diagnose the third moment of SGS vertical velocity
! Result is returned in the global variable w3

  call canuto()

! Recover parameters of the subgrid PDF using diagnosed moments
! and calculate SGS cloudiness, condensation and its effects on temeperature
! and moisture variables

  call assumed_pdf()

  tkh=tkh_out

contains

  logical function in_cloud(ql, qi)

    real, intent(in) :: ql, qi

    in_cloud = (ql + qi > qt_thresh) .or. (qi > qi_thresh)
    
  end function in_cloud

  subroutine eddy_length()

! This subroutine computes the turbulent length scale based on a
! formulation described in BK13

! Local variables
    real    wrk, wrk1, wrk2, wrk3, smixt_type(nx,ny,nzm)
    integer i, j, k, kk, kl, ku, kb, kc, kli, kui

    smixt_type = zero
    
    do j=1,ny
      do i=1,nx

        cbl_top(i,j) = one !zero
        cldarr(i,j)  = zero
        numer(i,j)   = zero
        denom(i,j)   = zero
      enddo
    enddo

 ! Reinitialize the mixing length related arrays to zero
!         smixt(i,j,k)    = one   
          smixt    = epsln 
          brunt_test    = zero

    
! Find the length scale outside of clouds
    
!    do k=1,nzm-1
      do j=1,ny
        do i=1,nx
kloop:    do k=1,nzm-1

!Eq. 11 in BK13 (Eq. 4.13 in Pete's dissertation)
!Outside of cloud, integrate from the surface to the cloud base

! old
!!!          if (qcl(i,j,k)+qci(i,j,k) <= zero) then 
!!           if (qcl(i,j,k)+qci(i,j,k) <= thresh) then
!!            tkes       = sqrt(tke(i,j,k)) * adzl(i,j,k)
!!            numer(i,j) = numer(i,j) + tkes*zl(i,j,k) ! Numerator   in Eq. 11 in BK13
!!            denom(i,j) = denom(i,j) + tkes           ! Denominator in Eq. 11 in BK13
!!          else
!!            cldarr(i,j) = one   ! Take note of the columns containing cloud.
!!          endif

!           if (cbl_top(i,j) .or. cldarr(i,j) == one ) exit kloop


!           if (qcl(i,j,k)+qci(i,j,k) >  thresh) then 
            if ( in_cloud(qcl(i,j,k),qci(i,j,k)) ) then 

!              cldarr(i,j) = one   ! Take note of the columns containing cloud(s).
               cldarr(i,j) = k   ! Take note of the columns containing cloud(s).
              exit kloop

           endif


! Look for clear convective boundary layer top
!           if (sflux(i,j) > 0 .and. &  ! Unstable BL
!                ((wthv_sec(i,j,k)<=0).and.(wthv_sec(i,j,k)<=wthv_sec(i,j,k+1)))) cbl_top(i,j)=.true.
          

!           if (qcl(i,j,k)+qci(i,j,k) <= thresh) then
!           if (cldarr(i,j) == zero ) then
!           if ( .not.  cbl_top(i,j)) then
!           if ( cbl_top(i,j) ==  zero) then
              tkes       = sqrt(tke(i,j,k)) * adzl(i,j,k)
              numer(i,j) = numer(i,j) + tkes*zl(i,j,k) ! Numerator   in Eq. 11 in BK13
              denom(i,j) = denom(i,j) + tkes           ! Denominator in Eq. 11 in BK13

!              if (sflux(i,j) > 0 .and.  &  ! Unstable BL
!!                   ((wthv_sec(i,j,k)<=0).and.(wthv_sec(i,j,k)<=wthv_sec(i,j,k+1)))) cbl_top(i,j)=.true.
!                   ((wthv_sec(i,j,k) < 0).and.(wthv_sec(i,j,k) < wthv_sec(i,j,k+1)))) cbl_top(i,j)=k
!           else
!            if ( cldarr(i,j)>zero )  exit kloop
!           endif

 
        enddo kloop
      enddo
    enddo 

! Calculate the measure of PBL depth,  Eq. 11 in BK13 
    do j=1,ny
      do i=1,nx
        if (denom(i,j) >  zero .and. numer(i,j) > zero) then
          l_inf(i,j) = min(0.1 * (numer(i,j)/denom(i,j)), max_l_inf)
!           l_inf(i,j) = min(0.1 * (numer(i,j)/denom(i,j)), 100.0)    
        else
!          l_inf(i,j) = 100.0
           l_inf(i,j) = max_l_inf     
        endif
      enddo
    enddo
    
!Calculate length scale outside of cloud, Eq. 10 in BK13 (Eq. 4.12 in Pete's dissertation)
    do k=1,nzm

      kb = k-1
      kc = k+1
      if (k == 1) then
        kb = 1
        kc = 2
        thedz(:,:) = adzi(:,:,kc)
      elseif (k == nzm) then
        kb = nzm-1
        kc = nzm
        thedz(:,:) = adzi(:,:,k)
      else
        thedz(:,:) = adzi(:,:,kc) + adzi(:,:,k) !  = (z(k+1)-z(k-1))
      endif

      do j=1,ny
        do i=1,nx

! bet = ggr/tv0 
       
          betdz = bet(i,j,k) / thedz(i,j)
           
          tkes = sqrt(tke(i,j,k))
             
! Compute local Brunt-Vaisalla frequency on layer centers
          
!          if (.false.) then

!             print *, "Shouldn't be here"

          wrk = qcl(i,j,k) + qci(i,j,k)

!          if (wrk > zero) then            ! If in the cloud
!          if (wrk > thresh) then            ! If in the cloud
           if (in_cloud(qcl(i,j,k),qci(i,j,k))) then 
    
! Find the in-cloud Brunt-Vaisalla frequency
                
             omn = qcl(i,j,k) / (wrk+1.e-20) ! Ratio of liquid water to total water

! Latent heat of phase transformation based on relative water phase content
! fac_cond = lcond/cp, fac_fus = lfus/cp

             lstarn = fac_cond + (one-omn)*fac_fus

! Derivative of saturation mixing ratio over water/ice wrt temp. based on 
! relative water phase content
             dqsat =      omn  * dtqsatw(tabs(i,j,k),prsl(i,j,k))            &
                   + (one-omn) * dtqsati(tabs(i,j,k),prsl(i,j,k))

! Saturation mixing ratio over water/ice wrt temp  based on relative water phase content

             qsatt =      omn  * qsatw(tabs(i,j,k),prsl(i,j,k))               &
                   + (one-omn) * qsati(tabs(i,j,k),prsl(i,j,k))

             bbb = (one + epsv*qsatt-wrk-qpl(i,j,k)-qpi(i,j,k)                &
                 + 1.61*tabs(i,j,k)*dqsat) / (one+lstarn*dqsat)

! Calculate Brunt-Vaisalla frequency using centered differences in the vertical
! Eq 5.22 in Pete's dissertation
             brunt_test(i,j,k) = betdz*(bbb*(hl(i,j,kc)-hl(i,j,kb))                &
                          + (bbb*lstarn - (one+lstarn*dqsat)*tabs(i,j,k))     &
                          * (total_water(i,j,kc)-total_water(i,j,kb))         & 
                          + (bbb*fac_cond - (one+fac_cond*dqsat)*tabs(i,j,k))*(qpl(i,j,kc)-qpl(i,j,kb))  &
                          + (bbb*fac_sub  - (one+fac_sub*dqsat)*tabs(i,j,k))*(qpi(i,j,kc)-qpi(i,j,kb)) )
                
          else                       ! outside of cloud
                
! Find outside-of-cloud Brunt-Vaisalla frequency
! Only unsaturated air, rain and snow contribute to virt. pot. temp. 
! Eq 5.23 in Pete's dissertation   

             bbb = one + epsv*qv(i,j,k) - qpl(i,j,k) - qpi(i,j,k)
             brunt_test(i,j,k) = betdz*( bbb*(hl(i,j,kc)-hl(i,j,kb))                        &
                          + epsv*tabs(i,j,k)*(total_water(i,j,kc)-total_water(i,j,kb)) &
                          + (bbb*fac_cond-tabs(i,j,k))*(qpl(i,j,kc)-qpl(i,j,kb))       &
                          + (bbb*fac_sub -tabs(i,j,k))*(qpi(i,j,kc)-qpi(i,j,kb)) )
          endif
           
!       endif
             
! Reduction of mixing length in the stable regions (where B.-V. freq. > 0) is required.
! Here we find regions of Brunt-Vaisalla freq. > 0 for later use. 

          if (brunt(i,j,k) >= zero) then
            brunt2(i,j,k) = brunt(i,j,k)
          else
            brunt2(i,j,k) = zero
          endif
             
! Calculate turbulent length scale in the boundary layer.
! See Eq. 10 in BK13 (Eq. 4.12 in Pete's dissertation)

! Keep the length scale small near the surface following Blackadar (1984)
! Note that this is not documented in BK13 and was added later for SP-CAM runs

!         if (k == 1) then
!           term = 600.*tkes
!           smixt(i,j,k) = term + (0.4*zl(i,j,k)-term)*exp(-zl(i,j,k)*0.01)
!         else

! tscale is the eddy turnover time scale in the boundary layer and is 
! an empirically derived constant 

            if (tkes > zero .and. l_inf(i,j) > zero) then
              wrk1 = one / (tscale*tkes*vonk*zl(i,j,k))
!              wrk2 = one / (tscale*tkes*l_inf(i,j))

              if (k < cldarr(i,j)) then 
                 wrk2 = one / (tscale*tkes*l_inf(i,j))
              else
                 wrk2 = one / (tscale*tkes*max_l_inf)
              endif
              wrk1 = wrk1 + wrk2 + slts * brunt2(i,j,k) / tke(i,j,k)
!              wrk1 = wrk1 + wrk2 + pt01 * brunt2(i,j,k) / tke(i,j,k)
              wrk1 = sqrt(one / max(wrk1,1.0e-8)) * (one/0.3)
!             smixt(i,j,k) = min(max_eddy_length_scale, 2.8284*sqrt(wrk1)/0.3)
              smixt(i,j,k) = min(max_eddy_length_scale, wrk1)

!              if ( smixt(i,j,k) > 2950.) &
!              if ( tabs(i,j,k) > 325) &
!              if ( tke(i,j,k) >  40) &
!              if ( brunt2(i,j,k) >  1 .and. (wrk <= thresh)) &
!              if ( brunt2(i,j,k) >  0 .and. (k <= 3)) &
!              if (  l_inf(i,j) < max_l_inf  .and. (k <= 3)) &
!              if (  smixt(i,j,k) > 2* max_l_inf ) &
!              if (  (cbl_top(i,j)  .or.  cldarr(i,j) == one) .and. k == 1) &

!              if (  (cbl_top(i,j) > zero  .or.  cldarr(i,j) > zero) .and. k == 1) &
!                   print *, "Outside cloud smixt=", smixt(i,j,k),  " k=", k," zl=", zl(i,j,k), &
!                   "wrk1=", one / (tscale*tkes*vonk*zl(i,j,k)),  &
!                   " l_inf=", l_inf(i,j), " wrk2=", wrk2, &
!                   " sqrt(tke)=", sqrt(tke(i,j,k)), " brunt2=", brunt2(i,j,k),  &
!                   " qcl(i,j,k)=", qcl(i,j,k), " qci(i,j,k)=", qci(i,j,k), &
!                   " denom(i,j)=", denom(i,j)," numer(i,j)=", numer(i,j), &
!                   " sflux(i,k)=", sflux(i,j), &
!                   " brunt_test(i,j,k)=", brunt_test(i,j,k),  " wthv_sec(i,j,:)=", wthv_sec(i,j,k), &
!                   " cbl_top(i,j)=", cbl_top(i,j), " cldarr(i,j)=", cldarr(i,j), "wrk > thresh", wrk > thresh!, &


!              if (  (cbl_top(i,j) > zero  .and.  cldarr(i,j) > zero) .and. k == 1) &
!                   print *, "First layer smixt=", smixt(i,j,k),  " k=", k," zl=", zl(i,j,k), &
!                   " denom(i,j)=", denom(i,j)," numer(i,j)=", numer(i,j), " l_inf=", l_inf(i,j),&
!                   " sqrt(tke)=", sqrt(tke(i,j,k)), " brunt2=", brunt2(i,j,k), &
!                   " sflux(i,k)=", sflux(i,j), &
!                   " brunt_test(i,j,k)=", brunt_test(i,j,k), " wthv_sec(i,j,k)=", wthv_sec(i,j,k), &
!                   " qcl(i,j,k)=", qcl(i,j,k), " qci(i,j,k)=", qci(i,j,k), &!


!                   "CBL top smixt=", smixt(i,j,cbl_top(i,j)),  " cbl_top(i,j) =", cbl_top(i,j),  &
!                   " zl(cbl_top(i,j) )=", zl(i,j,cbl_top(i,j)), &
!                   " sqrt(tke)=", sqrt(tke(i,j,cbl_top(i,j))), " brunt2=", brunt2(i,j,cbl_top(i,j)),  &
!                   " brunt_test(i,j,k)=", brunt_test(i,j,cbl_top(i,j)), " wthv_sec(i,j,k)=", wthv_sec(i,j,cbl_top(i,j)), &
!                   " qcl(i,j,k)=", qcl(i,j,cbl_top(i,j)), " qci(i,j,k)=", qci(i,j,cbl_top(i,j)), &


!                   "Cloud bottom smixt=", smixt(i,j,cldarr(i,j)),  " cldarr(i,j) =", cldarr(i,j),  &
!                   " zl(cldarr(i,j) )=", zl(i,j,cldarr(i,j)), &
!                   " sqrt(tke)=", sqrt(tke(i,j,cldarr(i,j))), " brunt2=", brunt2(i,j,cldarr(i,j)),  &
!                   " brunt_test(i,j,k)=", brunt_test(i,j,cldarr(i,j)), " wthv_sec(i,j,k)=", wthv_sec(i,j,cldarr(i,j)), &
!                   " qcl(i,j,k)=", qcl(i,j,cldarr(i,j)), " qci(i,j,k)=", qci(i,j,cldarr(i,j))
                   


              
!              if (wrk1 > max_eddy_length_scale) print *, "Bndry lr", smixt(i,j,k), wrk1, wrk2, brunt2(i,j,k), tke(i,j,k), one / (tscale*tkes*vonk*zl(i,j,k))

!           smixt(i,j,k) = min(max_eddy_length_scale,(2.8284*sqrt(1./((1./(tscale*tkes*vonk*zl(i,j,k))) & 
!                  + (1./(tscale*tkes*l_inf(i,j)))+0.01*(brunt2(i,j,k)/tke(i,j,k)))))/0.3)
!           else
!             smixt(i,j,k) = zero
!            endif
             endif
             
           enddo
          
        enddo
     enddo


    do k=1,nzm

      kb = k-1
      kc = k+1
      if (k == 1) then
        kb = 1
        kc = 2
        thedz(:,:) = adzi(:,:,kc)
      elseif (k == nzm) then
        kb = nzm-1
        kc = nzm
        thedz(:,:) = adzi(:,:,k)
      else
        thedz(:,:) = adzi(:,:,kc) + adzi(:,:,k) !  = (z(k+1)-z(k-1))
      endif

      do j=1,ny
        do i=1,nx

           if (.false.) &
!           if (  (cbl_top(i,j) > zero  .and.  cldarr(i,j) > zero) .and. k == 1 &
!                .and. (( cldarr(i,j) - cbl_top(i,j)) < 11)) &
                print *, "First layer smixt=", smixt(i,j,k),  " k=", k," zl=", zl(i,j,k), &
                " denom(i,j)=", denom(i,j)," numer(i,j)=", numer(i,j), " l_inf=", l_inf(i,j),&
                " sqrt(tke)=", sqrt(tke(i,j,k)), " brunt2=", brunt2(i,j,k), &
                " sflux(i,k)=", sflux(i,j), &
                " brunt_test(i,j,k)=", brunt_test(i,j,k), " wthv_sec(i,j,k)=", wthv_sec(i,j,k), &
                " qcl(i,j,k)=", qcl(i,j,k), " qci(i,j,k)=", qci(i,j,k), &
                
                new_line('a') , & 

                "CBL top smixt=", smixt(i,j,cbl_top(i,j)),  " cbl_top(i,j) =", cbl_top(i,j),  &
                " zl(cbl_top(i,j) )=", zl(i,j,cbl_top(i,j)), &
                " sqrt(tke)=", sqrt(tke(i,j,cbl_top(i,j))), " brunt2=", brunt2(i,j,cbl_top(i,j)),  &
                " brunt_test(i,j,k)=", brunt_test(i,j,cbl_top(i,j))," wthv_sec(i,j,k)=", wthv_sec(i,j,cbl_top(i,j)), &
                " qcl(i,j,k)=", qcl(i,j,cbl_top(i,j)), " qci(i,j,k)=", qci(i,j,cbl_top(i,j)), &

                new_line('a') ,& 

                "Cloud bottom smixt=", smixt(i,j,cldarr(i,j)),  " cldarr(i,j) =", cldarr(i,j),  &
                " zl(cldarr(i,j) )=", zl(i,j,cldarr(i,j)), &
                " sqrt(tke)=", sqrt(tke(i,j,cldarr(i,j))), " brunt2=", brunt2(i,j,cldarr(i,j)),  &
                " brunt_test(i,j,k)=", brunt_test(i,j,cldarr(i,j)), " wthv_sec(i,j,k)=", wthv_sec(i,j,cldarr(i,j)), &
                " qcl(i,j,k)=", qcl(i,j,cldarr(i,j)), " qci(i,j,k)=", qci(i,j,cldarr(i,j)), &

                new_line('a') , & 

                 " wthv_sec(i,j,:)=", wthv_sec(i,j,1:cldarr(i,j))
                   

             
           enddo
        enddo
     enddo

    
! Now find the in-cloud turbulence length scale 
! See Eq. 13 in BK13 (Eq. 4.18 in Pete's disseration)  
    
    do j=1,ny
      do i=1,nx
          
!        if (cldarr(i,j) == 1) then ! If there's a cloud in this column 
        if (cldarr(i,j) > zero) then ! If there's a cloud in this column 
             
          kl = 0
          ku = 0
          do k=2,nzm-3
                
! Look for the cloud base in this column  
            wrk = qcl(i,j,k) + qci(i,j,k)
 !           if (wrk > thresh .and. kl == 0) then
!            if (in_cloud(qcl(i,j,k),qci(i,j,k)) .and. kl == 0) then
            if (in_cloud(qcl(i,j,k),qci(i,j,k))  .and. kl == 0) then
              kl = k
            endif
                
! Look for the cloud top in this column

            if (cloud_top_at_min_negative_boyu_flux) then   
! Set cloud top as  the minimum negative boyuancy flux level or cloud water below threshold
               if (kl > zero .and. (wthv_sec(i,j,k)<0 .and. wthv_sec(i,j,k)<wthv_sec(i,j,k+1)) .or. &
                    (in_cloud(qcl(i,j,k),qci(i,j,k)) .and. .not. in_cloud(qcl(i,j,k+1),qci(i,j,k+1)))) ku = k
!               if (kl > zero .and. (wthv_sec(i,j,k)<0 .and. wthv_sec(i,j,k)<wthv_sec(i,j,k+1)) .or. &
!                    (wrk > thresh .and. qcl(i,j,k+1)+qci(i,j,k+1) <= thresh) ) ku = k
            else 
! Set cloud top as the level with cloud condensate below threshold
!             if (wrk > thresh .and. qcl(i,j,k+1)+qci(i,j,k+1) <= thresh) ku = k
                 if   ( in_cloud(qcl(i,j,k),qci(i,j,k)) .and. .not. in_cloud(qcl(i,j,k+1),qci(i,j,k+1)) ) ku = k

            endif

                
! Compute the mixing length scale for the cloud layer that we just found
! Why are single layer clouds (ku-kl=1) not recognized? 
!            if (kl > 0 .and. ku > 0 .and. ku-kl > 1) then
!               if (kl > 0 .and. ku > 0 .and. ku-kl > 0) then
               if (kl > 0 .and. ku >= kl ) then
               

! Determine cubed convective velocity scale (conv_vel2) inside the cloud
! See Eq. 16 in BK13 (Eq. 4.21 in Pete's dissertation)  


               conv_var=zero
               do kk=kl,ku
                  conv_var = conv_var + 2.5*adzi(i,j,kk)*bet(i,j,kk)*wthv_sec(i,j,kk)
               enddo
               conv_var = conv_var**oneb3

               smixt_type(i,j,kl:ku) = conv_var
           
                
              if (conv_var > 0) then 

!Probably should use adzi(i,j,k1) instead.
!                depth = (zl(i,j,ku)-zl(i,j,kl)) + adzl(i,j,kl)
                depth = (zl(i,j,ku)-zl(i,j,kl)) + adzi(i,j,kl)

                      
                     
                do kk=kl,ku
! In-cloud turbulence length scale, Eq. 13 in BK13 (Eq. 4.18)

!                  wrk = conv_var/(depth*sqrt(tke(i,j,kk)))
!                  wrk = wrk * wrk + slts*brunt2(i,j,kk)/tke(i,j,kk)
!                  wrk = wrk * wrk + pt01*brunt2(i,j,kk)/tke(i,j,kk)
! This is what's actually in the paper. 
! Both seem like plausible parameterizations. 
! This one has L grow slower with TKE
                  wrk = conv_var/(depth*depth*sqrt(tke(i,j,kk)))
                  wrk = wrk  + slts*brunt2(i,j,kk)/tke(i,j,kk)

!               if (wrk > 0) then
!              if ( (one/0.3)*sqrt(one/wrk) > max_eddy_length_scale) print *, "Cloud", wrk,  (one/0.3)*sqrt(one/wrk), conv_var/(depth*sqrt(tke(i,j,kk))), depth, brunt2(i,j,kk), tke(i,j,kk), slts*brunt2(i,j,kk)/tke(i,j,kk), conv_var, kk, kl, ku
!           endif
   
                  if (wrk> zero) then

!                     smixt_type(i,j,k) = one
!                     if ((one/0.3)*sqrt(one/wrk) > depth .and. tke(i,j,kk) > min_tke) & 
!                          if ( depth > max_eddy_length_scale.and. (one/0.3)*sqrt(one/wrk) > max_eddy_length_scale .and. tke(i,j,kk) >  40) &


!                     smixt(i,j,kk) = min(max_eddy_length_scale, (one/0.3)*sqrt(one/wrk))
                     smixt(i,j,kk) = min(depth, (one/0.3)*sqrt(one/wrk))

!!                         if ( smixt(i,j,kk) >  2950. .and. kl < 4) &
!                         if ( tabs(i,j,kk) >  325) &
!!                         if ( tke(i,j,kk) >  40) &
!!                            if ( brunt2(i,j,kk) >  1) &
!                          print *, "In-cloud smixt=",  smixt(i,j,kk) , " k=", k, " kk=", kk, &
!                          " kl=",kl," ku=",ku,  " depth=", depth," conv_var=", conv_var, &
!                          " sqrt(tke)=", sqrt(tke(i,j,kk)), " brunt2=", brunt2(i,j,kk),  &
!                          " qcl(i,j,k)=", qcl(i,j,kk), " qci(i,j,k)=", qci(i,j,kk), &
!                          " adzi(i,j,kl)=", adzi(i,j,kl), " adzl(i,j,kl)=", adzl(i,j,kl), &
!                          " thv(i,j,kk)=", thv(i,j,kk), " tabs(i,j,kk)=", tabs(i,j,kk), &
!                          " new smixt+", depth*sqrt(sqrt(tke(i,j,kk))/conv_var)





                  endif

               enddo
                      
              endif ! If convective vertical velocity scale > 0


!           if (  flag .and. (conv_var > 0.) .and. (cbl_top(i,j) == one  .and.  cldarr(i,j) > zero) .and. k == ku &

           if ( .false. & !  flag &
!           if (  k==nzm-3 &
!                .and. (( cldarr(i,j) - cbl_top(i,j)) < 11)) 
!                .and. (ku - kl > 0 ) .and. kl < 5 ) then 
!               .and. any(tke(i,j,:)>5.)) then 
!              .and. any(tabs(i,j,:)<150.)) then 

              .and. ((lprnt .and. i == ipr .and. any(tabs(i,j,:)<220.)) &
               .or. any(tke(i,j,:)>20.))) then 

                ku = nzm -1
                print *, "First layer smixt=", smixt(i,j,1),  " k=", k," zl=", zl(i,j,1), &
                " denom(i,j)=", denom(i,j)," numer(i,j)=", numer(i,j), " l_inf=", l_inf(i,j),&
!                " sqrt(tke)=", sqrt(tke(i,j,1)), " brunt2=", brunt2(i,j,1), &
                " sflux(i,k)=", sflux(i,j), &
                " brunt_test(i,j,k)=", brunt_test(i,j,1), & !" wthv_sec(i,j,k)=", wthv_sec(i,j,k), &
!                " qcl(i,j,k)=", qcl(i,j,k), " qci(i,j,k)=", qci(i,j,k), &
                
                new_line('a') , & 

!                "CBL top smixt=", smixt(i,j,cbl_top(i,j)),  " cbl_top(i,j) =", cbl_top(i,j),  &
!                " zl(cbl_top(i,j) )=", zl(i,j,cbl_top(i,j)), &
!                " sqrt(tke)=", sqrt(tke(i,j,cbl_top(i,j))), " brunt2=", brunt2(i,j,cbl_top(i,j)),  &
!                " brunt_test(i,j,k)=", brunt_test(i,j,cbl_top(i,j))," wthv_sec(i,j,k)=", wthv_sec(i,j,cbl_top(i,j)), &
!                " qcl(i,j,k)=", qcl(i,j,cbl_top(i,j)), " qci(i,j,k)=", qci(i,j,cbl_top(i,j)), &

!                new_line('a') ,& 

                "Cloud bottom smixt=", smixt(i,j,cldarr(i,j)),  " cldarr(i,j) =", cldarr(i,j),  &
                " zl(bottom )=", zl(i,j,kl), &
                 "kl=",kl," ku=",ku,  " depth=", depth," conv_var=", conv_var, &
                 
                 new_line('a') , &

                 " thv(i,j,kk)=", thv(i,j,1:ku+1), &

                 new_line('a') , &

                 " tabs(i,j,kk)=", tabs(i,j,1:ku+1), &

!                " sqrt(tke)=", sqrt(tke(i,j,cldarr(i,j))), " brunt2=", brunt2(i,j,cldarr(i,j)),  &
!                " brunt_test(i,j,k)=", brunt_test(i,j,cldarr(i,j)), " wthv_sec(i,j,k)=", wthv_sec(i,j,cldarr(i,j)), &
!                " qcl(i,j,k)=", qcl(i,j,cldarr(i,j)), " qci(i,j,k)=", qci(i,j,cldarr(i,j)), &

                new_line('a') , & 

 !                " wthv_sec(i,j,:)=", wthv_sec(i,j,1:cldarr(i,j))
                " wthv_sec(i,j,:)=", wthv_sec(i,j,1:ku+1), &

                new_line('a') , & 

                 " smixt(i,j,:)=", smixt(i,j,1:ku+1), &

!                new_line('a') , & 

!                " new smixt(i,j,:)=", depth*sqrt(sqrt(tke(i,j,kl:ku))/conv_var), &

                new_line('a') , & 

                " tke(i,j,:)=", tke(i,j,1:ku+1), &

                new_line('a') , & 

                " qt(i,j,:)=", qcl(i,j,1:ku+1)+qci(i,j,1:ku+1) 

                flag = .false.

             endif



              kl = zero
              ku = zero
           endif ! if inside the cloud layer

           if ( .false.  .and. k==nzm-3 &
!                .and. (( cldarr(i,j) - cbl_top(i,j)) < 11)) 
!                .and. (ku - kl > 0 ) .and. kl < 5 ) then 
!               .and. any(tke(i,j,:)>5.)) then 
!              .and. any(tabs(i,j,:)<150.)) then 

!              .and. ((lprnt .and. i == ipr .and. any(tabs(i,j,:)<220.)) &
!               .or. any(tke(i,j,:)>20.))) then 
!           .and. any(tke(i,j,:)>20.)) then 
           .and. any(smixt(i,j,1:5)>7000.)) then 
              
                kl =1
                ku = nzm -1
                print *, "First layer smixt=", smixt(i,j,1),  " k=", k," zl=", zl(i,j,1), &
                " denom(i,j)=", denom(i,j)," numer(i,j)=", numer(i,j), " l_inf=", l_inf(i,j),&
!                " sqrt(tke)=", sqrt(tke(i,j,1)), " brunt2=", brunt2(i,j,1), &
                " sflux(i,k)=", sflux(i,j), &
                " brunt_test(i,j,k)=", brunt_test(i,j,1), & !" wthv_sec(i,j,k)=", wthv_sec(i,j,k), &
!                " qcl(i,j,k)=", qcl(i,j,k), " qci(i,j,k)=", qci(i,j,k), &
                
                new_line('a') , & 

!                "CBL top smixt=", smixt(i,j,cbl_top(i,j)),  " cbl_top(i,j) =", cbl_top(i,j),  &
!                " zl(cbl_top(i,j) )=", zl(i,j,cbl_top(i,j)), &
!                " sqrt(tke)=", sqrt(tke(i,j,cbl_top(i,j))), " brunt2=", brunt2(i,j,cbl_top(i,j)),  &
!                " brunt_test(i,j,k)=", brunt_test(i,j,cbl_top(i,j))," wthv_sec(i,j,k)=", wthv_sec(i,j,cbl_top(i,j)), &
!                " qcl(i,j,k)=", qcl(i,j,cbl_top(i,j)), " qci(i,j,k)=", qci(i,j,cbl_top(i,j)), &

!                new_line('a') ,& 

                "Cloud bottom smixt=", smixt(i,j,cldarr(i,j)),  " cldarr(i,j) =", cldarr(i,j),  &
                " zl(bottom )=", zl(i,j,kl), &
                 "kl=",kl," ku=",ku,  " depth=", depth," conv_var=", conv_var, &
                 
!                 new_line('a') , &

!                 " thv(i,j,kk)=", thv(i,j,1:ku+1), &

                 new_line('a') , &

                 " tabs(i,j,kk)=", tabs(i,j,1:ku+1), &

!                " sqrt(tke)=", sqrt(tke(i,j,cldarr(i,j))), " brunt2=", brunt2(i,j,cldarr(i,j)),  &
!                " brunt_test(i,j,k)=", brunt_test(i,j,cldarr(i,j)), " wthv_sec(i,j,k)=", wthv_sec(i,j,cldarr(i,j)), &
!                " qcl(i,j,k)=", qcl(i,j,cldarr(i,j)), " qci(i,j,k)=", qci(i,j,cldarr(i,j)), &

                new_line('a') , & 

                " wthv_sec(i,j,:)=", wthv_sec(i,j,1:ku+1), &


                 new_line('a') , & 

                " brunt(i,j,:)=", brunt(i,j,1:ku+1), &

                
                 new_line('a') , & 

                " brunt_test(i,j,:)=", brunt_test(i,j,1:ku+1), &


                new_line('a') , & 

                " smixt(i,j,:)=", smixt(i,j,1:ku+1), &

                new_line('a') , & 

                " smixt_type(i,j,:)=", smixt_type(i,j,1:ku+1), &


!                new_line('a') , & 

!                " new smixt(i,j,:)=", depth*sqrt(sqrt(tke(i,j,kl:ku))/conv_var), &

                new_line('a') , & 

                " tke(i,j,:)=", tke(i,j,1:ku+1), &

                new_line('a') , & 

                " qt(i,j,:)=", qcl(i,j,1:ku+1)+qci(i,j,1:ku+1) 

                flag = .false.

              kl = zero
              ku = zero
           endif 

                
        enddo   ! k=2,nzm-3
     endif     ! if in the cloudy column
  enddo       ! i=1,nx
enddo         ! j=1,ny
    
! Constraints on eddy length scale values

    do k=1,nzm

       if (k == nzm) then
          kb = k
       else
          kb = k+1
       endif
       
       do j=1,ny
          do i=1,nx
             
             wrk = 0.1*adzl(i,j,k)
          
! Minimum 0.1 of local dz                      
!             smixt(i,j,k) = max(wrk, min(max_eddy_length_scale,smixt(i,j,k)))
             smixt(i,j,k) = max(wrk, smixt(i,j,k))
!             smixt(i,j,k) = max(epsln, smixt(i,j,k))


! If chracteristic grid dimension in the horizontal< 1000m, set lengthscale to 
! be not larger that that.
     
!         if (sqrt(dx*dy) .le. 1000.) smixt(i,j,k)=min(sqrt(dx*dy),smixt(i,j,k)) 

! If at the cloud top and atmosphere is stable, set to  0.1 of local dz   
! to ensure that inversion is represented.  
                

!             if (qcl(i,j,kb) == 0 .and. qcl(i,j,k) > 0 .and. brunt2(i,j,k) > 1.e-4) then
             if ( in_cloud(qcl(i,j,k),qci(i,j,k)) .and. .not. in_cloud(qcl(i,j,kb),qci(i,j,kb)) .and. &
                   brunt2(i,j,k) > 1e-4 ) then  
!             if (qcl(i,j,kb) + qci(i,j,kb) < thresh .and. qcl(i,j,k) + qci(i,j,k) >= thresh .and. brunt2(i,j,k) > 1.e-4) then
!             if (qcl(i,j,kb)  < thresh .and. qcl(i,j,k) >= thresh .and. brunt2(i,j,k) > 1.e-4) then

                smixt(i,j,k) = wrk
!                smixt = epsln
             endif

!             if (smixt(i,j,k) > max_eddy_length_scale .and. zl(i,j,k) <=max_eddy_length_scale) &
!                  smixt(i,j,k) = max_eddy_length_scale

             if ( tke(i,j,k) >  50) &
                  print *, "smixt=", smixt (i,j,k), " k=", k," zl=", zl(i,j,k), &
                  " l_inf=", l_inf(i,j),  &
                  " sqrt(tke)=", sqrt(tke(i,j,k)), " brunt2=", brunt2(i,j,k),  &
                  " qcl(i,j,k)=", qcl(i,j,k), " qci(i,j,k)=", qci(i,j,k)

          enddo ! i
       enddo   ! j
    enddo     ! k     
    
  end subroutine eddy_length


  subroutine canuto()

! Subroutine impements an analytic expression for the third moment of SGS vertical 
! velocity based on Canuto et al, 2001, JAS, 58, 1169-1172 
! https://doi.org/10.1175/1520-0469(2001)058<1169:NTOMFT>2.0.CO;2
! (further referred to as C01)
! This allows to avoid having a prognostic equation for the third moment.
! Result is returned in a global variable w3 defined at the layer interfaces
    
! Local variables
    integer i, j, k, kb, kc

    real bet2,   f0,     f1,     f2,  f3,    f4,   f5,  iso, isosqr,         &
         omega0, omega1, omega2, X0,  Y0,    X1,   Y1,  AA0, AA1, buoy_sgs2, &
         cond_w, cond_hl,  wrk, wrk1,  wrk2, wrk3, avew
    real  z, Z0, Z1 

! See Eq. 7 in C01 (B.7 in Pete's dissertation)
    real, parameter :: c=7.0, a0=0.52/(c*c*(c-2.)), a1=0.87/(c*c),      &
                       a2=0.5/c, a3=0.6/(c*(c-2.)), a4=2.4/(3.*c+5.),   &
                       a5=0.6/(c*(3.*c+5))
!Moorthi               a5=0.6/(c*(3.+5.*c))
    
!   do k=1,nzm
    do k=2,nzm

      kb = k-1
      kc = k+1

      if(k == 1) then
!        kb = 1
!        kc = 2
!        do j=1,ny
!          do i=1,nx
!            thedz(i,j)  = one / adzl(i,j,kc)
!            thedz2(i,j) = thedz(i,j)
!          enddo
!        enddo
      elseif(k == nzm) then
        kb = nzm-1
        kc = nzm
        do j=1,ny
          do i=1,nx

            thedz(i,j)  = one / adzi(i,j,k)
            thedz2(i,j) = one / adzl(i,j,k-1)

          enddo
        enddo
      else
        do j=1,ny
          do i=1,nx
            thedz(i,j)  = one / adzi(i,j,k)
            thedz2(i,j) = one / (adzl(i,j,k)+adzl(i,j,kb))
          enddo
        enddo
      endif

       
      do j=1,ny
        do i=1,nx
      
!intrp       
! Interpolate "return-to-isotropy" time scale to layer interfaces
!         iso       = half*(isotropy(i,j,k)+isotropy(i,j,kb))
          iso       = isotropy_int(i,j,k)
          isosqr    = iso*iso 
!intrp    
! Interpolate BV frequency  to layer interfaces 
!         buoy_sgs2 = isosqr*half*(brunt(i,j,k)+brunt(i,j,kb))
          buoy_sgs2 = isosqr*brunt_int(i,j,k)
!intrp    
!         bet2      = half*(bet(i,j,k)+bet(i,j,kb))  
          bet2      = bet_int(i,j,k)

        
! Compute functions f0-f5, see Eq, 8 in C01 (B.8 in Pete's dissertation)
        
!intrp
!         avew = half*(w_sec(i,j,k)+w_sec(i,j,kb)) ! Interpolate second moment of w to the interface
          avew = w_sec_int(i,j,k)
          avew_save(i,j,k)=avew
!          cond_w = 1.2*sqrt(max(1.0e-20,2.*avew*avew*avew)) 
!          cond_w_save(i,j,k)=cond_w
!          if (shoc_version == 1) then
!             cond_hl= 1.2*sqrt(max(1.0e-20,2.*(thl_sec(i,j,k)*thl_sec(i,j,k)*thl_sec(i,j,k))))
!          else
!             cond_hl= 1.2*sqrt(max(1.0e-20,2.*(thl_sec(i,j,k)+thl_sec(i,j,kb))**3))
!          endif

          wrk1 = bet2*iso
!         wrk2 = thedz2(i,j)*wrk1*wrk1*iso

          wrk  = wthl_sec(i,j,kc) - wthl_sec(i,j,kb)

! In SHOC v1 MSE variances are defined on layer interface
          if (shoc_version == 1) then
             wrk3 = thl_sec(i,j,kc) - thl_sec(i,j,kb)
             wrk2 = thedz2(i,j)*wrk1*wrk1*iso 
             f0   = wrk2 * wrk1 * wthl_sec(i,j,k) * wrk3
             f1   = wrk2 * (wrk*wthl_sec(i,j,k) + half*avew*wrk3)
          else
! but in SHOC v2 they are on layer centers
             wrk3 = thl_sec(i,j,k) - thl_sec(i,j,kb)
             wrk2 = wrk1*wrk1*iso
!             f0   = wrk2 * wrk1 * wthl_sec(i,j,k) * wrk3 / adzi(i,j,k)
!             f1   = wrk2 * (wrk*wthl_sec(i,j,k)/(adzl(i,j,k)+adzl(i,j,kb)) + half*avew*wrk3/adzi(i,j,k))
             f0   = wrk2 * wrk1 * wthl_sec(i,j,k) * wrk3 * thedz(i,j)
             f1   = wrk2 * (wrk*wthl_sec(i,j,k)*thedz2(i,j) + half*avew*wrk3*thedz(i,j)) 
          endif

!          f0   = wrk2 * wrk1 * wthl_sec(i,j,k) * wrk3

!          wrk  = wthl_sec(i,j,kc) - wthl_sec(i,j,kb)
             
!          f1   = wrk2 * (wrk*wthl_sec(i,j,k) + half*avew*wrk3)
             
          wrk1 = bet2*isosqr

          f2   = thedz(i,j)*wrk1*wthl_sec(i,j,k)*(w_sec(i,j,k)-w_sec(i,j,kb))     &
               + (thedz2(i,j)+thedz2(i,j))*bet(i,j,k)*isosqr*wrk


          f3   = thedz2(i,j)*wrk1*wrk + thedz(i,j)*bet2*isosqr*(wthl_sec(i,j,k)*(tke(i,j,k)-tke(i,j,kb)))


          wrk1 = thedz(i,j)*iso*avew
          f4   = wrk1*(w_sec(i,j,k)-w_sec(i,j,kb) + tke(i,j,k)-tke(i,j,kb))
 

          f5   = wrk1*(w_sec(i,j,k)-w_sec(i,j,kb))
             
       
! Compute the "omega" terms, see Eq. 6 in C01 (B.6 in Pete's dissertation)

          omega0 = a4 / (one-a5*buoy_sgs2)
          omega1 = omega0 / (c+c)
          omega2 = omega1*f3+(5./4.)*omega0*f4
 
! Compute the X0, Y0, X1, Y1 terms,  see Eq. 5 a-b in C01  (B.5 in Pete's dissertation)

          wrk1 = one / (one-(a1+a3)*buoy_sgs2)
          wrk2 = one / (one-a3*buoy_sgs2)
          X0   = wrk1 * (a2*buoy_sgs2*(one-a3*buoy_sgs2))
          Y0   = wrk2 * (two*a2*buoy_sgs2*X0)
          X1   = wrk1 * (a0*f0+a1*f1+a2*(one-a3*buoy_sgs2)*f2)
          Y1   = wrk2 * (two*a2*(buoy_sgs2*X1+(a0/a1)*f0+f1))

          Z0   = 3.*buoy_sgs2/(2*(c-2))
          Z1   = 3.*f0/(2*(c-2))

! Compute the A0, A1 terms,  see Eq. 5d in C01 (B.5 in Pete's dissertation)

          AA0 = omega0*X0 + omega1*Y0
          AA1 = omega0*X1 + omega1*Y1 + omega2

! Finally, we have the third moment of w, see Eq. 4c in C01 (B.4 in Pete's dissertation)
! cond is an estimate of the third moment from the second moment. If the third moment is larger
! than the estimate - limit w3.

           z = (AA1-1.2*X1-1.5*f5)/(c-1.2*X0+AA0)
           z_save(i,j,k) = z
!          w3(i,j,k) = max(-cond_w, min(cond_w, (AA1-1.2*X1-1.5*f5)/(c-1.2*X0+AA0)))
!           if (z < -2 ) print *, " z=",z," cond_w=",cond_w," avew=", avew
!           w3(i,j,k) = max(-cond_w, min(cond_w,z))    
           w3(i,j,k) = z 
           w3_save(i,j,k)= w3(i,j,k)

! Find the third moment of hl, see Eq 4a and 4b in C01

           hl3(i,j,k)= (Z0*(Y0*z-Y1)-Z1)*(bet2*iso)**3
!           hl3(i,j,k)= max(-cond_hl, min(cond_hl, Z0*(Y0*z-Y1)-Z1))

! Implemetation of the C01 approach in this subroutine is nearly complete
! (the missing part are Eqs. 5c and 5e which are very simple)
! therefore it's easy to diagnose other third order moments obtained in C01 using this code. 
! In particular, we can calculate hl3 and use it to find skewness of MSE SGS distribution 
! instead of assuming that it's zero in assumed_pdf()


        enddo
      enddo
    enddo
    do j=1,ny
      do i=1,nx
        w3(i,j,1) = w3(i,j,2)
        hl3(i,j,1) = hl3(i,j,2)
      enddo
    enddo
    
  end subroutine canuto

  subroutine assumed_pdf()

! Compute SGS buoyancy flux, SGS cloud fraction, and SGS condensation 
! using assumed analytic double-gaussian PDF for SGS vertical velocity, 
! moisture, and  liquid/ice water static energy, based on the 
! general approach of  Larson et al 2002, JAS, 59, 3519-3539, 
! https://doi.org/10.1175/1520-0469(2002)059<3519:SSAMVI>2.0.CO;2
! and Golaz et al 2002, JAS, 59, 3540-3551
! https://doi.org/10.1175/1520-0469(2002)059<3540:APBMFB>2.0.CO;2 
! A detailed derivation of relationships between PDF parameters and its 
! higher order moments can be found in Larson and Golaz, 2005, MWR, 133, 1023–1042
! https://doi.org/10.1175/MWR2902.1 
! References in the comments in this code are given to 
! the Appendix A of Pete Bogenschutz's dissertation. 
    
! Local variables

    integer i,j,k,ku,kd, iter
    real wrk, wrk1, wrk2, wrk3, wrk4, bastoeps, eps_ss1, eps_ss2

    real c_w_thl, c_w_qw, c_thl_qw, cond_w, cond_hl, lqs1, lqs2, beta1no_l, beta2no_l, max_w_skw_var

!   bastoeps = basetemp / epsterm


! Initialize for statistics
    do k=1,nzm
      wqlsb(k) = zero
      wqisb(k) = zero
    enddo

!   sfac  = scrit
!   sfaci = one / sfac
    
    DO k=1,nzm
      
      kd = k
      ku = k + 1
!     if (k == nzm) ku = k
      
      DO j=1,ny
        DO i=1,nx

! Initialize cloud variables to zero  
          diag_qn   = zero
          diag_frac = zero
          diag_ql   = zero
          diag_qi   = zero

          pval  = prsl(i,j,k)
          pfac  = pval * 1.0e-5
          pkap  = pfac ** kapa

!         sfac  = scrit * sqrt(pfac)
!         sfac  = scrit
          sfac  = scrit * pfac * pfac
          sfaci = one / sfac
             
! PDF code expects input moments to on layer centers. 
! In GFS  liquid/ice static energy, total water mixing ratio, 
! and vertical velocity are already there. 
          thl_first = hl(i,j,k) ! + fac_cond*qpl(i,k) + fac_sub*qpi(i,k) ! Don't include precip in the PDF
          qw_first  = total_water(i,j,k)
          w_first   = w(i,j,k)
            
             
! Second moments of PDF (except w_sec and prognostic thl_sec and qw_sec) and w3 are 
! calaculated as vertical differences of layer centered variables, therefore
! placing higher order  moments on the layer interfaces.
! We need to place them back on the layer centers

!intrp ok
          if (k < nzm) then
            w3var    = half*(w3(i,j,kd) + w3(i,j,ku))
            hl3var   = half*(hl3(i,j,kd)+ hl3(i,j,ku))
            if (shoc_version == 1) then 
               thlsec   = max(zero, half*(thl_sec(i,j,kd)+thl_sec(i,j,ku)) )
               qwsec    = max(zero, half*(qw_sec(i,j,kd)+qw_sec(i,j,ku)) )
            endif
            if (shoc_version == 2) then 
               thlsec   = max(zero,thl_sec(i,j,k))
               qwsec    = max(zero,qw_sec(i,j,k))
            endif
            qwthlsec = half * (qwthl_sec(i,j,kd) + qwthl_sec(i,j,ku))
            wqwsec   = half * (wqw_sec(i,j,kd)   + wqw_sec(i,j,ku))
            wthlsec  = half * (wthl_sec(i,j,kd)  + wthl_sec(i,j,ku))   
          else          
! At the model top interface assume zero values for moments
            w3var    = half*w3(i,j,k)
            hl3var   = half*hl3(i,j,k)
            if (shoc_version == 1) then
               thlsec   = max(zero, half*thl_sec(i,j,k))
               qwsec    = max(zero, half*qw_sec(i,j,k))
            endif
            if (shoc_version == 2) then
               thlsec   = max(zero, thl_sec(i,j,k))
               qwsec    = max(zero, qw_sec(i,j,k))
            endif
            qwthlsec = half * qwthl_sec(i,j,k)
            wqwsec   = half * wqw_sec(i,j,k)
            wthlsec  = half * wthl_sec(i,j,k)
          endif



          
          w3var_save(i,j,k) = w3var
          w_sec_save(i,j,k) = w_sec(i,j,k)

! Compute square roots of some variables so we don't have to do it again
          if (w_sec(i,j,k) > zero) then
            sqrtw2  = sqrt(w_sec(i,j,k))
          else
            sqrtw2  = zero
          endif
          if (thlsec > zero) then
            sqrtthl = sqrt(thlsec)
          else
            sqrtthl = zero
          endif
          if (qwsec > zero) then
            sqrtqt  = sqrt(qwsec)
          else
            sqrtqt  = zero
          endif



          if (variable_normalized_width_w == .true. .or. &
               larson_golaz_05_skew == .true.               )  then
!! Correlation of w and thl
!             c_w_thl  = max(-one,min(one,wthlsec/max(sqrtw2,w_tol)/max(sqrtthl,thl_tol)))
!! Correlation of w and qw
!             c_w_qw   = max(-one,min(one,wqwsec  /max(sqrtw2,w_tol)/max(sqrtqt, rt_tol)))
!! Correlation of w and qw
!             c_thl_qw = max(-one,min(one,qwthlsec/max(sqrtthl,thl_tol)/ &
!                                                     max(sqrtqt, rt_tol)))

! Correlation of w and thl
             c_w_thl  = wthlsec/max(sqrtw2,w_tol)/max(sqrtthl,thl_tol)
! Correlation of w and qw
             c_w_qw   = wqwsec  /max(sqrtw2,w_tol)/max(sqrtqt, rt_tol)
! Correlation of w and qw
             c_thl_qw = qwthlsec/max(sqrtthl,thl_tol)/ max(sqrtqt, rt_tol)

          endif
             

! Find parameters of the double Gaussian PDF of vertical velocity


          IF (w_sec(i,j,k) <= w_tol_sqd) THEN ! If variance of w is too small then
                                              ! PDF is a sum of two delta functions
!          if (sqrtw2 <= w_tol) then 

             Skew_w = zero
             w1_1   = w_first
             w1_2   = w_first
             w2_1   = zero
             w2_2   = zero
             aterm  = half
             onema  = half
             sqrtw2t = one
          ELSE


! Calculate normalized variances of w gaussians

            if (variable_normalized_width_w == .true. .or. &
               larson_golaz_05_skew == .true.             )  then 
! Variable normalized variances of w gaussians following Larson and Golaz (2005)
               w2_1 = max(zero, gamma_vnww * ( 1 - max(c_w_thl**2,c_w_qw**2)))
            else
! In the original formulation normalized variances of w gaussians are constant
               w2_1 = 0.4
            endif

            w2_2 = w2_1



            wrk  = 1 - 2*atmin
            wrk1 = 1 - w2_1*w2_1
            wrk1 = 4*wrk1*wrk1*wrk1
            max_w_skw_var= wrk * sqrt(wrk1/(1-wrk*wrk))
!            if (max_w_skw_var<max_w_skw) print *, "max_w_skw_var=", max_w_skw_var

!            cond_w = max_w_skw*max(w3_tol,sqrtw2*sqrtw2*sqrtw2)
!             cond_w = max_w_skw*sqrtw2*sqrtw2*sqrtw2
             cond_w = max_w_skw_var*sqrtw2*sqrtw2*sqrtw2
             cond_w_save(i,j,k)=cond_w
             w3var= max(-cond_w, min(cond_w,w3var))
!             w3var_save(i,j,k) = w3var  
! Skewness of vertical velocity PDF
             Skew_w = w3var / (sqrtw2*sqrtw2*sqrtw2)     ! Moorthi
!            if (Skew_w < -2 ) print *, " Skew_w=",Skew_w," w3var=", w3var," sqrtw2=", sqrtw2
             Skew_w_save(i,j,k) = Skew_w


                
! Compute relative weight of the first PDF "plume" 
! See Eq A4 in Pete's dissertaion -  Ensure 0.01 < a < 0.99

            wrk   = one - w2_1
            aterm = max(atmin,min(half*(one-Skew_w*sqrt(one/(4.*wrk*wrk*wrk+Skew_w*Skew_w))),atmax))
! Relative weight of the second PDF "plume"   
            onema = one - aterm
                
            sqrtw2t = sqrt(wrk)

! Normalized means of w gaussians                
! Eq. A.5-A.6
! A.6 in the text is missing a -1 factor on the right hand side.
            wrk  =   sqrt(onema/aterm)
            w1_1 =   sqrtw2t * wrk
            w1_2 = - sqrtw2t / wrk

! (Non-normalized) Variances of first and second plumes of vertical velocity
            w2_1 = w2_1 * w_sec(i,j,k)
            w2_2 = w2_2 * w_sec(i,j,k)



          ENDIF


          if (variable_normalized_width_w == .true. .or. &
               larson_golaz_05_skew == .true.               )  then
! Normalize correlation of w and thl
!             c_w_thl  = sign(max(-one,min(one,abs(c_w_thl)/sqrtw2t)), c_w_thl)
             c_w_thl  = max(-one,min(one,c_w_thl/sqrtw2t))
! Normalize correlation of w and qw
!             c_w_qw   = sign(max(-one,min(one,abs(c_w_qw)/sqrtw2t)), c_w_qw)
             c_w_qw   = max(-one,min(one,c_w_qw/sqrtw2t))
! Normalize correlation of w and qw
!             c_thl_qw = max(-one,min(one,c_thl_qw))
                                         

          endif

             
!  Find parameters of the  PDF of liquid/ice static energy

          IF (thlsec <= thl_tol*thl_tol .or. abs(w1_2-w1_1) <= w_thresh) THEN

!  If variance of h is too small then the PDF is a sum of two delta functions

            thl1_1     = thl_first
            thl1_2     = thl_first
            thl2_1     = zero
            thl2_2     = zero
            sqrtthl2_1 = zero
            sqrtthl2_2 = zero
          ELSE

! Correlation of w and thl
            corrtest1 = max(-one,min(one,wthlsec/(sqrtw2*sqrtthl)))

! Normalized means of thl gaussians 
            thl1_1 = -corrtest1 / w1_2                 ! A.7
            thl1_2 = -corrtest1 / w1_1                 ! A.8
                
            if ( larson_golaz_05_skew == .false.) then
!            if ( larson_golaz_05_skew == .true.) then
       

! Skewness of MSE PDF

! In original SHOC skew_fact=0.
               Skew_hl = -skew_fact*Skew_w

! Parametrize skewness of MSE PDF using Canuto et al (2001)
               if (canuto_skew) then
!!                  if (shoc_version == 1) then
!                     cond_hl= canuto_sk_hl_max*sqrt(max(1.0e-20,2.*(thl_sec(i,j,k)*thl_sec(i,j,k)*thl_sec(i,j,k))))
!                  else
!! intrp Fix interpolation from centers to interfaces here
!                     cond_hl= canuto_sk_hl_max*sqrt(max(1.0e-20,2.*(thl_sec(i,j,k)+thl_sec(i,j,kb))**3))
!                  endif

                  cond_hl= canuto_sk_hl_max*sqrtthl*sqrtthl*sqrtthl
                  hl3var= max(-cond_hl, min(cond_hl,hl3var))
                  Skew_hl = hl3var / (sqrtthl*sqrtthl*sqrtthl)
               endif

               wrk1   = thl1_1 * thl1_1
               wrk2   = thl1_2 * thl1_2
               wrk3   = three * (one - aterm*wrk1 - onema*wrk2)
               wrk4   = Skew_hl - aterm*wrk1*thl1_1 - onema*wrk2*thl1_2
!            wrk4   = -skew_fact*Skew_w - aterm*wrk1*thl1_1 - onema*wrk2*thl1_2  ! testing - Moorthi
! Skewness is assumed to be zero for thl in the original SHOC
!           wrk4   =     - aterm*wrk1*thl1_1 - onema*wrk2*thl1_2
! Eq. A.9-A.10, plus the condition A.11
               wrk    = three * (thl1_2-thl1_1)
               if (wrk /= zero) then
! Variances for MSE for first and second plumes
! A.9  and A.11
                  thl2_1 = thlsec * min(max_thl_norm_var, max(zero,( thl1_2*wrk3-wrk4)/(aterm*wrk))) 
! A.10 and A.11 (Note that A.10 contans a typo, but the code here is correct, 
! see A.22 in Larson et al 2002, JAS, 59, 3519-3539)
                  thl2_2 = thlsec * min(max_thl_norm_var ,max(zero,(-thl1_1*wrk3+wrk4)/(onema*wrk))) 
               else
                  thl2_1 = zero
                  thl2_2 = zero
               endif

            else !  larson_golaz_05_skew == .true.
               
               ! Eqns 34 and 35 in LG05
!               thl2_1= thlsec * min(max_thl_norm_var, max(zero,(one - (c_w_thl/sqrtw2t)**2)/aterm * &
!                       (beta_factor*oneb3 + aterm*(1-twoby3*beta_factor))))
!               thl2_2= thlsec * min(max_thl_norm_var, max(zero,(one - (c_w_thl/sqrtw2t)**2)/onema * &
!                       (one - (beta_factor*oneb3 + aterm*(1-twoby3*beta_factor)))))

               thl2_1= thlsec * min(max_thl_norm_var, max(zero,(one - (c_w_thl)**2)/aterm * &
                       (beta_factor*oneb3 + aterm*(1-twoby3*beta_factor))))
               thl2_2= thlsec * min(max_thl_norm_var, max(zero,(one - (c_w_thl)**2)/onema * &
                       (one - (beta_factor*oneb3 + aterm*(1-twoby3*beta_factor)))))


            endif

! Convert from normalized means of MSE  gaussians to the means of first and second plumes on the
! MSE axis, see Eqs. A.7 and A.8
! These are (non-normalized) means of MSE gaussians
            thl1_1 = thl1_1*sqrtthl + thl_first
            thl1_2 = thl1_2*sqrtthl + thl_first

! Standard deviations of MSE for first and second plumes 

            sqrtthl2_1 = sqrt(thl2_1)
            sqrtthl2_2 = sqrt(thl2_2)

          ENDIF

!  FIND PARAMETERS FOR TOTAL WATER MIXING RATIO

          IF (qwsec <= rt_tol*rt_tol .or. abs(w1_2-w1_1) <= w_thresh) THEN
! If variance of total water is too small  then the  PDF is a sum of two delta functions
            qw1_1     = qw_first
            qw1_2     = qw_first
            qw2_1     = zero
            qw2_2     = zero
            sqrtqw2_1 = zero
            sqrtqw2_2 = zero
          ELSE

! Correlation of qw and thl
            corrtest2 = max(-one,min(one,wqwsec/(sqrtw2*sqrtqt)))

! Normalized means of total water  gaussians 
            qw1_1 = - corrtest2 / w1_2            ! A.7
            qw1_2 = - corrtest2 / w1_1            ! A.8

            tsign = abs(qw1_2-qw1_1)

            if ( larson_golaz_05_skew == .false.) then

               
               if (canuto_skew) then

                  Skew_qw = canuto_factor*Skew_hl

               else ! Original parameterization of qw skewness

                  IF (tsign > 0.4) THEN
                     Skew_qw = skew_facw*Skew_w
                  ELSEIF (tsign <= 0.2) THEN
                     Skew_qw = zero
                  ELSE
! Linear interpolation between 0 and 1.2*Skew_w
                     Skew_qw = (skew_facw/0.2) * Skew_w * (tsign-0.2)
                  ENDIF

               endif


               wrk1  = qw1_1 * qw1_1
               wrk2  = qw1_2 * qw1_2
               wrk3  = three * (one - aterm*wrk1 - onema*wrk2)
               wrk4  = Skew_qw - aterm*wrk1*qw1_1 - onema*wrk2*qw1_2
               wrk   = three * (qw1_2-qw1_1)
               
! Eq. A.9-A.10, plus the condition A.11
! Variances of total water for first and second plumes

               if (wrk /= zero) then
! A.9  and A.11
                  qw2_1 = qwsec * min(max_qw_norm_var, max(zero,( qw1_2*wrk3-wrk4)/(aterm*wrk))) 
! A.10 and A.11 (Note that A.10 contans a typo, but the code here is correct, 
! see A.22 in Larson et al 2002, JAS, 59, 3519-3539) 
                  qw2_2 = qwsec * min(max_qw_norm_var, max(zero,(-qw1_1*wrk3+wrk4)/(onema*wrk))) 
               else
                  qw2_1 = zero
                  qw2_2 = zero
               endif

            else !  larson_golaz_05_skew == .true.

               ! Eqns 34 and 35 in LG05
               qw2_1= qwsec *  min(max_qw_norm_var, max(zero,(one - (c_w_qw/sqrtw2t)**2)/aterm * &
                       (beta_factor*oneb3 + aterm*(1-twoby3*beta_factor))))
               qw2_2= qwsec *  min(max_qw_norm_var, max(zero,(one - (c_w_qw/sqrtw2t)**2)/onema * &
                       (one - (beta_factor*oneb3 + aterm*(1-twoby3*beta_factor)))))


            endif


! (Non-normalized) Means of first and second plumes on the total water axis,  
! see Eqs. A.7 and A.8
            qw1_1 = qw1_1*sqrtqt + qw_first
            qw1_2 = qw1_2*sqrtqt + qw_first

! Standard deviations for total water for first and second plumes
            sqrtqw2_1 = sqrt(qw2_1)
            sqrtqw2_2 = sqrt(qw2_2)

          ENDIF

! (Non-normalized) Means of first and second plumes on the vertical velocity  axis,  
! see Eqs. A.7 and A.8


          w1_1 = w1_1*sqrtw2 + w_first
          w1_2 = w1_2*sqrtw2 + w_first

!  FIND WITHIN-PLUME CORRELATIONS 

          if ( larson_golaz_05_skew == .false.) then

!  See Eq. A.12-A.13
             
             testvar = aterm*sqrtqw2_1*sqrtthl2_1 + onema*sqrtqw2_2*sqrtthl2_2

             IF (testvar == 0) THEN
                r_qwthl_1 = zero
             ELSE
                r_qwthl_1 = max(-one,min(one,(qwthlsec-aterm*(qw1_1-qw_first)*(thl1_1-thl_first) &
                     -onema*(qw1_2-qw_first)*(thl1_2-thl_first))/testvar)) ! A.12
             ENDIF

          else  ! larson_golaz_05_skew == .true.
             if (sqrtw2t == 0. .or. (one-(c_w_thl)**2)*(one-(c_w_qw)**2)<0.) print *, "sqrtw2t=", sqrtw2t, "(one-(c_w_thl)**2)*(one-(c_w_qw)**2) =", (one-(c_w_thl)**2)*(one-(c_w_qw)**2), " c_w_thl=", c_w_thl, " c_w_qw=", c_w_qw, &
             " wthlsec=", wthlsec, " sqrtw2=", sqrtw2, " w_tol=",w_tol, " sqrtthl=",sqrtthl," thl_tol=", thl_tol
!             testvar =  sqrt((one-(c_w_thl/sqrtw2t)**2)*(one-(c_w_qw/sqrtw2t)**2))
             testvar =  sqrt((one-(c_w_thl)**2)*(one-(c_w_qw)**2))
             
             if (testvar == 0) then
                r_qwthl_1 = zero
             else
                ! Eqn 36 in  Larson and Golaz (2005)   
!                r_qwthl_1 = max(-one,min(one,(c_thl_qw - c_w_thl*c_w_qw/sqrtw2t**2)/testvar))
                r_qwthl_1 = max(-one,min(one,(c_thl_qw - c_w_thl*c_w_qw)/testvar))
                  
             endif

          endif

! Save PDF moments and parameters to the diagnostic output file.
          if (shoc_diag) then
     
             Diag%SHOC_qt(i,k)       = qw_first
             Diag%SHOC_hl(i,k)       = thl_first
             Diag%SHOC_w(i,k)        = w_first

             Diag%SHOC_qt_qt(i,k)    = qwsec
             Diag%SHOC_hl_hl(i,k)    = thlsec
             Diag%SHOC_w_w(i,k)      = w_sec(i,j,k)
             
             Diag%SHOC_w_qt(i,k)     = wqwsec
             Diag%SHOC_w_hl(i,k)     = wthlsec
             Diag%SHOC_qt_hl(i,k)    = qwthlsec

             Diag%SHOC_w_w_w(i,k)    = w3var
             Diag%SHOC_hl_hl_hl(i,k) = hl3var 

             Diag%SHOC_w_Sk(i,k)     = Skew_w 
             Diag%SHOC_qt_Sk(i,k)    = Skew_qw
             Diag%SHOC_hl_Sk(i,k)    = Skew_hl

             Diag%SHOC_G1_weight(i,k)     = aterm 

             Diag%SHOC_G1_w_nmean(i,k)    = w1_1 
             Diag%SHOC_G1_w_nvar(i,k)     = sqrt(w2_1)
             Diag%SHOC_G1_qt_nmean(i,k)   = qw1_1
             Diag%SHOC_G1_qt_nvar(i,k)    = sqrtqw2_1
             Diag%SHOC_G1_hl_nmean(i,k)   = thl1_1
             Diag%SHOC_G1_hl_nvar(i,k)    = sqrtthl2_1
             Diag%SHOC_G1_qt_hl_ncov(i,k) = r_qwthl_1 

             Diag%SHOC_G2_w_nmean(i,k)    = w1_2
             Diag%SHOC_G2_w_nvar(i,k)     = sqrt(w2_2) 
             Diag%SHOC_G2_qt_nmean(i,k)   = qw1_2
             Diag%SHOC_G2_qt_nvar(i,k)    = sqrtqw2_2
             Diag%SHOC_G2_hl_nmean(i,k)   = thl1_2
             Diag%SHOC_G2_hl_nvar(i,k)    = sqrtthl2_2
     
          endif




! Calculate condensation and cloud fraction
          if ( Firl_condensation == .true. ) then 


          wrk1  = gamaz(i,j,k)  - fac_cond*qpl(i,j,k) - fac_sub*qpi(i,j,k)


! Amounts of frozen and liquid condesates in each individual gaussian are not known,
! so we will assume zeros initially and find the solution by iteration

! Liquid condensate in each gaussian
          ql1 = zero
          ql2 = zero

! Frozen condensate in each gaussian
          qi1 = zero
          qi2 = zero

! Obtain liquid/ice temperature from liquid/ice static energies of the first and second "plumes"
          Tl1_1 = thl1_1 - wrk1
          Tl1_2 = thl1_2 - wrk1


          do iter=1,3

! Obtain absolute temperature from liquid/ice static energies of the first and second "plumes"
          T1_1 =  thl1_1 - fac_cond*ql1 - fac_sub*qi1
          T1_2 =  thl1_2 - fac_cond*ql2 - fac_sub*qi2


! Now compute qs

          esval1_1 = zero
          esval2_1 = zero
          om1      = one
          eps_ss1  = eps
          eps_ss2  = eps
             
!     Partition saturation vapor pressure based on temperature for the first "plume"

!          IF (Tl1_1 >= tbgmax) THEN       ! Use temperature, not liq/ice temp 
          IF (T1_1 >= tbgmax) THEN       ! Use temperature, not liq/ice temp 
            esval1_1 = fpvsl(Tl1_1) ! Saturation vapor pressure over water
!           esval1_1 = esatw(Tl1_1)
            lstarn1  = lcond
!          ELSE IF (Tl1_1 < tbgmin) THEN    ! Use temperature, not liq/ice temp 
          ELSE IF (T1_1 < tbgmin) THEN    ! Use temperature, not liq/ice temp 
            esval1_1 = fpvsi(Tl1_1) ! Saturation vapor pressure over ice
!           esval1_1 = esati(Tl1_1)
            lstarn1  = lsub
            eps_ss1   = eps * supice ! Require supersaturation over ice
          ELSE
            esval1_1 = fpvsl(Tl1_1)  ! Saturation vapor pressure over water
            esval2_1 = fpvsi(Tl1_1)  ! Saturation vapor pressure over ice
!           esval1_1 = esatw(Tl1_1)
!           esval2_1 = esati(Tl1_1)
!            om1      = max(zero, min(one, a_bg*(Tl1_1-tbgmin)))  ! Use temperature, not liq/ice temp 
            om1      = max(zero, min(one, a_bg*(T1_1-tbgmin)))  ! Use temperature, not liq/ice temp 
            lstarn1  = lcond + (one-om1)*lfus
            eps_ss2   = eps * supice
          ENDIF

!   Convert  from saturation pressure to mixing ratio wrt. water/ice for the first "plume"

          qs1   =      om1  * (eps_ss1*esval1_1/max(esval1_1,pval-0.378*esval1_1))      &
                + (one-om1) * (eps_ss2*esval2_1/max(esval2_1,pval-0.378*esval2_1))

          lqs1   =      om1  * lcond *(eps_ss1*esval1_1/max(esval1_1,pval-0.378*esval1_1))      &
                + (one-om1) * lsub * (eps_ss2*esval2_1/max(esval2_1,pval-0.378*esval2_1))

!         beta1 = (rgas/rv)*(lstarn1/(rgas*Tl1_1))*(lstarn1/(cp*Tl1_1))
          beta1 = (lstarn1*lstarn1*onebrvcp) / (Tl1_1*Tl1_1)              ! A.18
          beta1no_l = (lstarn1*onebrvcp) / (Tl1_1*Tl1_1)              ! A.18


! Are the temperatures of the two plumes equal?  If so then set qs and beta
! in each plume to each other 
          IF (Tl1_1 == Tl1_2) THEN
            qs2   = qs1     
            lqs2  = lqs1     
            beta2 = beta1
            beta2no_l = beta1no_l
          ELSE 

            esval1_2 = zero
            esval2_2 = zero
            om2      = one
            eps_ss1  = eps
            eps_ss2  = eps

!     Partition saturated vapor pressure based on temperature for the second "plume"


!            IF (Tl1_2 < tbgmin) THEN    ! Use temperature, not liq/ice temp 
            IF (T1_2 < tbgmin) THEN    ! Use temperature, not liq/ice temp 
              esval1_2 = fpvsi(Tl1_2)  ! Saturation vapor pressure over water
!             esval1_2 = esati(Tl1_2)
              lstarn2  = lsub
              eps_ss1   = eps * supice
!            ELSE IF (Tl1_2 >= tbgmax) THEN  ! Use temperature, not liq/ice temp 
            ELSE IF (T1_2 >= tbgmax) THEN  ! Use temperature, not liq/ice temp 
              esval1_2 = fpvsl(Tl1_2)  ! Saturation vapor pressure over ice
!             esval1_2 = esatw(Tl1_2)
              lstarn2  = lcond
            ELSE
              esval1_2 = fpvsl(Tl1_2)  ! Saturation vapor pressure over water
              esval2_2 = fpvsi(Tl1_2)  ! Saturation vapor pressure over ice
!             esval1_2 = esatw(Tl1_2)
!             esval2_2 = esati(Tl1_2)
!              om2      = max(zero, min(one, a_bg*(Tl1_2-tbgmin)))  ! Use temperature, not liq/ice temp 
              om2      = max(zero, min(one, a_bg*(T1_2-tbgmin)))  ! Use temperature, not liq/ice temp 
              lstarn2  = lcond + (one-om2)*lfus
              eps_ss2  = eps * supice
            ENDIF

!     Saturation mixing ratio wrt. water/ice for the second "plume"
            qs2   =      om2  * (eps_ss1*esval1_2/max(esval1_2,pval-0.378*esval1_2))    &
                  + (one-om2) * (eps_ss2*esval2_2/max(esval2_2,pval-0.378*esval2_2))

            lqs2   =      om2  * lcond * (eps_ss1*esval1_2/max(esval1_2,pval-0.378*esval1_2))    &
                  + (one-om2) * lsub *(eps_ss2*esval2_2/max(esval2_2,pval-0.378*esval2_2))
                
!           beta2 = (rgas/rv)*(lstarn2/(rgas*Tl1_2))*(lstarn2/(cp*Tl1_2))   ! A.18
            beta2 = (lstarn2*lstarn2*onebrvcp) / (Tl1_2*Tl1_2)              ! A.18
            beta2no_l = (lstarn2*onebrvcp) / (Tl1_2*Tl1_2)              ! A.18
                
          ENDIF

          qs1 = qs1 !* rhc(i,j,k)
          qs2 = qs2 !* rhc(i,j,k)

          lqs1 = lqs1 !* rhc(i,j,k)
          lqs2 = lqs2 !* rhc(i,j,k)

! Now compute non-precipitating cloud condensate  

! Compute SGS condensation and cloud fraction for the first "plume"

!          cqt1   = one / (one+beta1*qs1)                                    ! A.19
          cqt1   = one / (one + beta1no_l * lqs1)
          wrk    = qs1 * (one+beta1*qw1_1) * cqt1
          s1     = qw1_1 - wrk                                              ! A.17, first line
!          cthl1  = cqt1*wrk*cpolv*beta1*pkap                                ! A.20
          cthl1  = cqt1*lqs1*cp/lstarn1*beta1no_l*pkap                                ! A.20

          wrk1   = cthl1 * cthl1
          wrk2   = cqt1  * cqt1
!         std_s1 = sqrt(max(zero,wrk1*thl2_1+wrk2*qw2_1-2.*cthl1*sqrtthl2_1*cqt1*sqrtqw2_1*r_qwthl_1))
! Eq. A.17, second line
          std_s1 = sqrt(max(zero, wrk1*thl2_1+wrk2*qw2_1                        &
                                - two*cthl1*sqrtthl2_1*cqt1*sqrtqw2_1*r_qwthl_1))
             
          qn1 = zero
          C1  = zero
             
          IF (std_s1 > zero) THEN
            wrk = s1 / (std_s1*sqrt2)
! SGS cloud fraction for the first "plume" 
            C1 = max(zero, min(one, half*(one+erf(wrk))))                   ! A.15
!     if (lprnt .and. i == ipr) write(0,*)' in shoc wrk=',wrk,' s1=','std=',std_s1,&
!         ' c1=',c1*100,' qs1=',qs1,' qw1_1=',qw1_1

! SGS condensation for the first "plume" 
            IF (C1 > zero) qn1 = s1*C1 + (std_s1*sqrtpii)*exp(-wrk*wrk)     ! A.16
!            if (qn1 < sfac) then
!              c1  = min(c1, qn1*sfaci)
!            endif
          ELSEIF (s1 > zero) THEN
             c1 = one
!            C1  = min(one, max(zero,s1*sfaci))
            qn1 = s1
          ENDIF
             
! Compute SGS condensation and cloud fraction for the second "plume"  

! If two plumes exactly equal, then use results for the first "plume"

          IF (qw1_1 == qw1_2 .and. thl2_1 == thl2_2 .and. qs1 == qs2) THEN
            s2     = s1
            cthl2  = cthl1
            cqt2   = cqt1
            std_s2 = std_s1
            C2     = C1
            qn2    = qn1
          ELSE

!            cqt2   = one / (one+beta2*qs2)           ! Eq. A.19  
            cqt2   = one / (one+beta2no_l*lqs2)           ! Eq. A.19  
            wrk    = qs2 * (one+beta2*qw1_2) * cqt2
            s2     = qw1_2 - wrk                     ! Eq. A.17, first line
!            cthl2  = wrk*cqt2*cpolv*beta2*pkap       ! Eq. A.20  
            cthl2  = wrk*lqs2*cp/lstarn2*beta2no_l*pkap       ! Eq. A.20  
            wrk1   = cthl2 * cthl2
            wrk2   = cqt2  * cqt2
!           std_s2 = sqrt(max(zero,wrk1*thl2_2+wrk2*qw2_2-2.*cthl2*sqrtthl2_2*cqt2*sqrtqw2_2*r_qwthl_1))
! Eq. A.17, second line
            std_s2 = sqrt(max(zero, wrk1*thl2_2+wrk2*qw2_2                        &
                                  - two*cthl2*sqrtthl2_2*cqt2*sqrtqw2_2*r_qwthl_1))

            qn2 = zero
            C2  = zero

            IF (std_s2 > zero) THEN
              wrk = s2 / (std_s2*sqrt2)
! SGS cloud fraction for the second "plume"
              C2  = max(zero, min(one, half*(one+erf(wrk))))
! SGS condensation for the second "plume" 
              IF (C2 > zero) qn2 = s2*C2 + (std_s2*sqrtpii)*exp(-wrk*wrk)
!              if (qn2 < sfac) then
!                c2  = min(c2, qn2*sfaci)
!              endif
            ELSEIF (s2 > zero) THEN
!              C2  = min(one, max(zero,s2*sfaci))
               c2 = one
              qn2 = s2
            ENDIF
               
          ENDIF

! Total SGS cloud fraction
          diag_frac = aterm*C1 + onema*C2 ! A.14
            
! Fraction of liquid conensate in each gaussian 
!          om1 = max(zero, min(one, (Tl1_1-tbgmin)*a_bg))  ! Use temperature, not liq/ice temp 
!          om2 = max(zero, min(one, (Tl1_2-tbgmin)*a_bg))  ! Use temperature, not liq/ice temp 
          om1 = max(zero, min(one, (T1_1-tbgmin)*a_bg))  ! Use temperature, not liq/ice temp 
          om2 = max(zero, min(one, (T1_2-tbgmin)*a_bg))  ! Use temperature, not liq/ice temp 
             
! Can't condense more than grid average total water in each gaussian
          qn1 = min(qn1,qw1_1)
          qn2 = min(qn2,qw1_2)
             
! Find liquid condensate in each gaussian
          ql1 = qn1*om1
          ql2 = qn2*om2
             
! Find frozen condensate in each gaussian
          qi1 = qn1 - ql1
          qi2 = qn2 - ql2
             
!     if (lprnt .and. i == ipr) write(0,*)' in shoc qi=',qi1,qi2,' ql=',ql1,ql2,&
!        ' c1=',c1,' c2=',c2,' s1=',s1,' s2=',s2,' k=',k

! Total amount of condensate
          diag_qn = min(max(zero, aterm*qn1 + onema*qn2), total_water(i,j,k))
! Liquid condensate
          diag_ql = min(max(zero, aterm*ql1 + onema*ql2), diag_qn)
! Frozen condesate
          diag_qi = diag_qn - diag_ql

          wrk = tabs(i,j,k)
             
! Update temperature variable based on diagnosed cloud properties
!          om1         = max(zero, min(one, (tabs(i,j,k)-tbgmin)*a_bg)) ! The only partition based on abs temp
!          lstarn1     = lcond + (one-om1)*lfus
          tabs(i,j,k) = hl(i,j,k) - gamaz(i,j,k) + fac_cond*(diag_ql+qpl(i,j,k)) &
                                                 + fac_sub *(diag_qi+qpi(i,j,k)) &
                      + tkesbdiss(i,j,k) * (dtn/cp)      ! tke dissipative heating

! Update moisture fields

! Update ncpl and ncpi Anning Cheng 03/11/2016
!         ncpl(i,j,k)    = diag_ql/max(qc(i,j,k),1.e-10)*ncpl(i,j,k)
! The following commneted by Moorthi on April 26, 2017 to test blowing up
!         ncpl(i,j,k)    = (1.0-diag_ql/max(qc(i,j,k),1.e-10)) * ncpl(i,j,k)
!         ncpi(i,j,k)    = (1.0-diag_qi/max(qi(i,j,k),1.e-10)) * ncpi(i,j,k)
          qc(i,j,k)      = diag_ql
          qi(i,j,k)      = diag_qi
          qwv(i,j,k)     = total_water(i,j,k) - diag_qn
          cld_sgs(i,j,k) = diag_frac

!! Update ncpl and ncpi Moorthi  12/12/2018                                                                                  
!        if (imp_phys > 0) then
!          if (ncpl(i,k) > nmin) then
!            ncpl(i,k) = diag_ql/max(qc(i,k),1.0d-10)*ncpl(i,k)
!          else
!            ncpl(i,k) = max(diag_ql/(fourb3*pi*RL_cub*997.0d0), nmin)
!          endif
!          if (ncpi(i,k) > nmin) then
!            ncpi(i,k) = diag_qi/max(qi(i,k),1.0d-10)*ncpi(i,k)
!          else
!            ncpi(i,k) = max(diag_qi/(fourb3*pi*RI_cub*500.0d0), nmin)
!          endif
!        endif

          if ( iter > 1 .and. wrk-tabs(i,j,k) < 0.001) exit

          end do ! iter
     
         else  ! Standard SHOC condensation code

            wrk1  = gamaz(i,j,k)  - fac_cond*qpl(i,j,k) - fac_sub*qpi(i,j,k)


! Amounts of frozen and liquid condesates in each individual gaussian are not known,
! so we will assume zeros initially and find the solution by iteration

! Liquid condensate in each gaussian
          ql1 = zero
          ql2 = zero

! Frozen condensate in each gaussian
          qi1 = zero
          qi2 = zero

! Obtain liquid/ice temperature from liquid/ice static energies of the first and second "plumes"
          Tl1_1 = thl1_1 - wrk1
          Tl1_2 = thl1_2 - wrk1
! Obtain absolute temperature from liquid/ice static energies of the first and second "plumes"
!          T1_1 =  thl1_1 - fac_cond*ql1 - fac_sub*qi1
!          T1_2 =  thl1_2 - fac_cond*ql2 - fac_sub*qi2


! Now compute qs

          esval1_1 = zero
          esval2_1 = zero
          om1      = one
          eps_ss1  = eps
          eps_ss2  = eps
             
!     Partition saturation vapor pressure based on temperature for the first "plume"

          IF (Tl1_1 >= tbgmax) THEN       ! Use temperature, not liq/ice temp 
!          IF (T1_1 >= tbgmax) THEN       ! Use temperature, not liq/ice temp 
            esval1_1 = fpvsl(Tl1_1) ! Saturation vapor pressure over water
!           esval1_1 = esatw(Tl1_1)
            lstarn1  = lcond
          ELSE IF (Tl1_1 < tbgmin) THEN    ! Use temperature, not liq/ice temp 
!          ELSE IF (T1_1 < tbgmin) THEN    ! Use temperature, not liq/ice temp 
            esval1_1 = fpvsi(Tl1_1) ! Saturation vapor pressure over ice
!           esval1_1 = esati(Tl1_1)
            lstarn1  = lsub
            eps_ss1   = eps * supice ! Require supersaturation over ice
          ELSE
            esval1_1 = fpvsl(Tl1_1)  ! Saturation vapor pressure over water
            esval2_1 = fpvsi(Tl1_1)  ! Saturation vapor pressure over ice
!           esval1_1 = esatw(Tl1_1)
!           esval2_1 = esati(Tl1_1)
            om1      = max(zero, min(one, a_bg*(Tl1_1-tbgmin)))  ! Use temperature, not liq/ice temp 
!            om1      = max(zero, min(one, a_bg*(T1_1-tbgmin)))  ! Use temperature, not liq/ice temp 
            lstarn1  = lcond + (one-om1)*lfus
            eps_ss2   = eps * supice
          ENDIF

!   Convert  from saturation pressure to mixing ratio wrt. water/ice for the first "plume"

          qs1   =      om1  * (eps_ss1*esval1_1/max(esval1_1,pval-0.378*esval1_1))      &
                + (one-om1) * (eps_ss2*esval2_1/max(esval2_1,pval-0.378*esval2_1))

!         beta1 = (rgas/rv)*(lstarn1/(rgas*Tl1_1))*(lstarn1/(cp*Tl1_1))
          beta1 = (lstarn1*lstarn1*onebrvcp) / (Tl1_1*Tl1_1)              ! A.18


! Are the temperatures of the two plumes equal?  If so then set qs and beta
! in each column to each other to save computation
          IF (Tl1_1 == Tl1_2) THEN
            qs2   = qs1     
            beta2 = beta1
          ELSE 

            esval1_2 = zero
            esval2_2 = zero
            om2      = one
            eps_ss1  = eps
            eps_ss2  = eps

!     Partition saturated vapor pressure based on temperature for the second "plume"


            IF (Tl1_2 < tbgmin) THEN    ! Use temperature, not liq/ice temp 
!            IF (T1_2 < tbgmin) THEN    ! Use temperature, not liq/ice temp 
              esval1_2 = fpvsi(Tl1_2)  ! Saturation vapor pressure over water
!             esval1_2 = esati(Tl1_2)
              lstarn2  = lsub
              eps_ss1   = eps * supice
            ELSE IF (Tl1_2 >= tbgmax) THEN  ! Use temperature, not liq/ice temp 
!            ELSE IF (T1_2 >= tbgmax) THEN  ! Use temperature, not liq/ice temp 
              esval1_2 = fpvsl(Tl1_2)  ! Saturation vapor pressure over ice
!             esval1_2 = esatw(Tl1_2)
              lstarn2  = lcond
            ELSE
              esval1_2 = fpvsl(Tl1_2)  ! Saturation vapor pressure over water
              esval2_2 = fpvsi(Tl1_2)  ! Saturation vapor pressure over ice
!             esval1_2 = esatw(Tl1_2)
!             esval2_2 = esati(Tl1_2)
              om2      = max(zero, min(one, a_bg*(Tl1_2-tbgmin)))  ! Use temperature, not liq/ice temp 
!              om2      = max(zero, min(one, a_bg*(T1_2-tbgmin)))  ! Use temperature, not liq/ice temp 
              lstarn2  = lcond + (one-om2)*lfus
              eps_ss2  = eps * supice
            ENDIF

!     Saturation mixing ratio wrt. water/ice for the second "plume"
            qs2   =      om2  * (eps_ss1*esval1_2/max(esval1_2,pval-0.378*esval1_2))    &
                  + (one-om2) * (eps_ss2*esval2_2/max(esval2_2,pval-0.378*esval2_2))
                
!           beta2 = (rgas/rv)*(lstarn2/(rgas*Tl1_2))*(lstarn2/(cp*Tl1_2))   ! A.18
            beta2 = (lstarn2*lstarn2*onebrvcp) / (Tl1_2*Tl1_2)              ! A.18
                
          ENDIF

          qs1 = qs1 !* rhc(i,j,k)
          qs2 = qs2 !* rhc(i,j,k)

! Now compute non-precipitating cloud condensate  

! Compute SGS condensation and cloud fraction for the first "plume"

          cqt1   = one / (one+beta1*qs1)                                    ! A.19
          wrk    = qs1 * (one+beta1*qw1_1) * cqt1
          s1     = qw1_1 - wrk                                              ! A.17, first line
          cthl1  = cqt1*wrk*cpolv*beta1*pkap                                ! A.20

          wrk1   = cthl1 * cthl1
          wrk2   = cqt1  * cqt1
!         std_s1 = sqrt(max(zero,wrk1*thl2_1+wrk2*qw2_1-2.*cthl1*sqrtthl2_1*cqt1*sqrtqw2_1*r_qwthl_1))
! Eq. A.17, second line
          std_s1 = sqrt(max(zero, wrk1*thl2_1+wrk2*qw2_1                        &
                                - two*cthl1*sqrtthl2_1*cqt1*sqrtqw2_1*r_qwthl_1))
             
          qn1 = zero
          C1  = zero
             
          IF (std_s1 > zero) THEN
            wrk = s1 / (std_s1*sqrt2)
! SGS cloud fraction for the first "plume" 
            C1 = max(zero, min(one, half*(one+erf(wrk))))                   ! A.15
!     if (lprnt .and. i == ipr) write(0,*)' in shoc wrk=',wrk,' s1=','std=',std_s1,&
!         ' c1=',c1*100,' qs1=',qs1,' qw1_1=',qw1_1

! SGS condensation for the first "plume" 
            IF (C1 > zero) qn1 = s1*C1 + (std_s1*sqrtpii)*exp(-wrk*wrk)     ! A.16
!            if (qn1 < sfac) then
!              c1  = min(c1, qn1*sfaci)
!            endif
          ELSEIF (s1 > zero) THEN
             c1 = one
!            C1  = min(one, max(zero,s1*sfaci))
            qn1 = s1
          ENDIF
             
! Compute SGS condensation and cloud fraction for the second "plume"  

! If two plumes exactly equal, then use results for the first "plume"

          IF (qw1_1 == qw1_2 .and. thl2_1 == thl2_2 .and. qs1 == qs2) THEN
            s2     = s1
            cthl2  = cthl1
            cqt2   = cqt1
            std_s2 = std_s1
            C2     = C1
            qn2    = qn1
          ELSE

            cqt2   = one / (one+beta2*qs2)           ! Eq. A.19  
            wrk    = qs2 * (one+beta2*qw1_2) * cqt2
            s2     = qw1_2 - wrk                     ! Eq. A.17, first line
            cthl2  = wrk*cqt2*cpolv*beta2*pkap       ! Eq. A.20  
            wrk1   = cthl2 * cthl2
            wrk2   = cqt2  * cqt2
!           std_s2 = sqrt(max(zero,wrk1*thl2_2+wrk2*qw2_2-2.*cthl2*sqrtthl2_2*cqt2*sqrtqw2_2*r_qwthl_1))
! Eq. A.17, second line
            std_s2 = sqrt(max(zero, wrk1*thl2_2+wrk2*qw2_2                        &
                                  - two*cthl2*sqrtthl2_2*cqt2*sqrtqw2_2*r_qwthl_1))

            qn2 = zero
            C2  = zero

            IF (std_s2 > zero) THEN
              wrk = s2 / (std_s2*sqrt2)
! SGS cloud fraction for the second "plume"
              C2  = max(zero, min(one, half*(one+erf(wrk))))
! SGS condensation for the second "plume" 
              IF (C2 > zero) qn2 = s2*C2 + (std_s2*sqrtpii)*exp(-wrk*wrk)
!              if (qn2 < sfac) then
!                c2  = min(c2, qn2*sfaci)
!              endif
            ELSEIF (s2 > zero) THEN
!              C2  = min(one, max(zero,s2*sfaci))
               c2 = one
              qn2 = s2
            ENDIF
               
          ENDIF

! Total SGS cloud fraction
          diag_frac = aterm*C1 + onema*C2 ! A.14
            
! Fraction of liquid conensate in each gaussian 
          om1 = max(zero, min(one, (Tl1_1-tbgmin)*a_bg))  ! Use temperature, not liq/ice temp 
          om2 = max(zero, min(one, (Tl1_2-tbgmin)*a_bg))  ! Use temperature, not liq/ice temp 
!          om1 = max(zero, min(one, (T1_1-tbgmin)*a_bg))  ! Use temperature, not liq/ice temp 
!          om2 = max(zero, min(one, (T1_2-tbgmin)*a_bg))  ! Use temperature, not liq/ice temp 
             
! Can't condense more than grid average total water in each gaussian
          qn1 = min(qn1,qw1_1)
          qn2 = min(qn2,qw1_2)
             
! Find liquid condensate in each gaussian
          ql1 = qn1*om1
          ql2 = qn2*om2
             
! Find frozen condensate in each gaussian
          qi1 = qn1 - ql1
          qi2 = qn2 - ql2
             
!     if (lprnt .and. i == ipr) write(0,*)' in shoc qi=',qi1,qi2,' ql=',ql1,ql2,&
!        ' c1=',c1,' c2=',c2,' s1=',s1,' s2=',s2,' k=',k

! Total amount of condensate
          diag_qn = min(max(zero, aterm*qn1 + onema*qn2), total_water(i,j,k))
! Liquid condensate
          diag_ql = min(max(zero, aterm*ql1 + onema*ql2), diag_qn)
! Frozen condesate
          diag_qi = diag_qn - diag_ql

          wrk = tabs(i,j,k)
             
! Update temperature variable based on diagnosed cloud properties
!          om1         = max(zero, min(one, (tabs(i,j,k)-tbgmin)*a_bg)) ! The only partition based on abs temp
!          lstarn1     = lcond + (one-om1)*lfus
          tabs(i,j,k) = hl(i,j,k) - gamaz(i,j,k) + fac_cond*(diag_ql+qpl(i,j,k)) &
                                                 + fac_sub *(diag_qi+qpi(i,j,k)) &
                      + tkesbdiss(i,j,k) * (dtn/cp)      ! tke dissipative heating

! Update moisture fields

! Update ncpl and ncpi Anning Cheng 03/11/2016
!         ncpl(i,j,k)    = diag_ql/max(qc(i,j,k),1.e-10)*ncpl(i,j,k)
! The following commneted by Moorthi on April 26, 2017 to test blowing up
!         ncpl(i,j,k)    = (1.0-diag_ql/max(qc(i,j,k),1.e-10)) * ncpl(i,j,k)
!         ncpi(i,j,k)    = (1.0-diag_qi/max(qi(i,j,k),1.e-10)) * ncpi(i,j,k)
          qc(i,j,k)      = diag_ql
          qi(i,j,k)      = diag_qi
          qwv(i,j,k)     = total_water(i,j,k) - diag_qn
          cld_sgs(i,j,k) = diag_frac

 !         if (wrk-tabs(i,j,k) < 0.001) exit
!          end do ! iter
        endif ! Firl_condesation == .true


! Compute the liquid water flux, see A.14
          wqls = aterm * ((w1_1-w_first)*ql1) + onema * ((w1_2-w_first)*ql2)
          wqis = aterm * ((w1_1-w_first)*qi1) + onema * ((w1_2-w_first)*qi2)
             
! Compute statistics for the fluxes 
          wqlsb(k) = wqlsb(k) + wqls
          wqisb(k) = wqisb(k) + wqis
             
! Diagnostic buoyancy flux on layer center.
!  Includes effects from liquid water, ice condensate, liquid & ice precipitation
!         wrk = epsv * basetemp
          wrk = epsv * thv(i,j,k)

          om1         = max(zero, min(one, (tabs(i,j,k)-tbgmin)*a_bg)) ! The only partition based on abs temp
          lstarn1     = lcond + (one-om1)*lfus


          bastoeps = onebeps * thv(i,j,k)

!intrp ok
          if (k < nzm) then
            wthv_sec(i,j,k) = wthlsec + wrk*wqwsec                                    &
                            + (fac_cond-bastoeps)*wqls                                &
                            + (fac_sub-bastoeps) *wqis                                &
                            + ((lstarn1/cp)-thv(i,j,k))*half*(wqp_sec(i,j,kd)+wqp_sec(i,j,ku))
          else
            wthv_sec(i,j,k) = wthlsec + wrk*wqwsec                                    &
                            + (fac_cond-bastoeps)*wqls                                &
                            + (fac_sub-bastoeps) *wqis                                &
                            + ((lstarn1/cp)-thv(i,j,k))*half*wqp_sec(i,j,k)
          endif

!           wthv_sec(i,j,k) = wthlsec + wrk*wqwsec                                     &
!                         + (fac_cond-bastoeps)*wqls                                 &
!                         + (fac_sub-bastoeps)*wqis                                  &
!                         + ((lstarn1/cp)-basetemp)*half*(wqp_sec(i,j,kd)+wqp_sec(i,j,ku))

! Diagnose <w'qt'qt'> and <w'theta'theta'> on the layer center.
! Vertical gradient of these quantities  is used in prognostic equations for variances.                                                      
!          if ( larson_golaz_05_skew == .false.) then

! Eqs 29 in  GLC02
             wthl2(i,j,k)  = aterm*(w1_1 - w_first)*((thl1_1 - thl_first)**2 + thl2_1) &
                  + onema*(w1_2 - w_first)*((thl1_2 - thl_first)**2 + thl2_2)

             wqw2(i,j,k)   = aterm*(w1_1 - w_first)*((qw1_1 - qw_first)**2 + qw2_1) &
                  + onema*(w1_2 - w_first)*((qw1_2 - qw_first)**2 + qw2_2)

!          else  ! larson_golaz_05_skew == .true.


!             if (w_sec(i,j,k) <= w_tol_sqd) then 

!                 wthl2(i,j,k)  = zero
!                 wqw2(i,j,k)   = zero

!              else
                 
!                 wrk =  one / ( w_sec(i,j,k)*(one - w2_1) )

!! Eq 45 in LG05
!                 wthl2(i,j,k) = wrk*w3var*(beta_factor*oneb3*thlsec + &
!                                             wrk*(one-beta_factor*oneb3)*wthlsec*wthlsec)
!                 wqw2 (i,j,k) = wrk*w3var*(beta_factor*oneb3*qwsec + &
!                                             wrk*(one-beta_factor*oneb3)*wqwsec * wqwsec)
                 
!              endif
             
!          endif

             wrk =  ( w_sec(i,j,k)*(one - w2_1) )
             if (wrk .ne. zero) then 
                wrk =  one / wrk
             else 
                wrk = zero
             end if

! Modification of Eq. 26 in GLC02

!test
!        if (k == 1) then
!           Cek(i,j,k) = Ces * Cesfac
!        else
!            Cek(i,j,k) = Ce  * Cefac
!        endif


!        Cek(i,j,k) = Cv * max(one, sqrt(pcrit/prsl(i,j,k)))     
        Cek(i,j,k) = 1.
!Apply additional damping for the PDFs with values of the individual gaussian weights
! in the ranges (atmin:atmin_damp) and  (atmax:atmax_damp) 
        if (aterm < atmin_damp) Cek(i,j,k)=Cek(i,j,k) * (one + at_damp_strength*(one - & 
                                           (aterm-atmin)/(atmin_damp-atmin)))
        if (aterm > atmax_damp) Cek(i,j,k)=Cek(i,j,k) * (one + at_damp_strength*(one - & 
                                           (atmax-aterm)/(atmax-atmax_damp)))

!        if (abs(Skew_w) > 4 ) Cek(i,j,k)=Cek(i,j,k) * (one + at_damp_strength*(one - & 
!                                           (aterm-atmin)/(atmin_damp-atmin)))
!        if (aterm ) Cek(i,j,k)=Cek(i,j,k) * (one + at_damp_strength*(one - & 
!                                           (atmax-aterm)/(atmax-atmax_damp)))


 if (.false.) then
!         if (k==nzm .and. any(smixt(i,j,1:5)>7000)) then
!        if ( any(wthv_sec(i,j,1:10) > 0.5 )) then 
!           if (abs((wthl_sec_d(i,j,k)- wthl_sec(i,j,k))/ wthl_sec(i,j,k)) < 10) &
                print *, " wthv_sec=", wthv_sec(i,j,1:nzm), new_line ('a'), &!    
                " wthl_sec=", wthl_sec(i,j,1:nzm), new_line ('a'), &!
                     " wthl_sec_d(k)=", wthl_sec_d(i,j,1:nzm),  new_line ('a'),&
                     " wqw_sec=", wqw_sec(i,j,1:nzm), new_line ('a'), &!
                     " wqw_sec_d=", wqw_sec_d(i,j,1:nzm), new_line ('a'), &!
                     " epsv * thv(i,j,k) wqw_sec=", epsv*thv(i,j,1:nzm)*wqw_sec(i,j,1:nzm), new_line ('a'), &!

!                " wthl_sec(k-1)=", wthl_sec(i,j,k-1), &
 !               " k=", k, "adzl=", adzl(i,j,1:nzm), new_line ('a')," dtabsdt=", dtabsdt(i,j,1:nzm),  new_line ('a'),&
 !               " fac_cond*dqcldt=", fac_cond* dqcldt(i,j,1:nzm),  new_line ('a'),&
 !               " fac_cond*dqpldt=", fac_cond* dqpldt(i,j,1:nzm),  new_line ('a'),&
 !               " fac_sub*dqcidt=", fac_sub* dqcidt(i,j,1:nzm) ,  new_line ('a'),&
!                " fac_sub*dqpidt=", fac_sub* dqpidt(i,j,1:nzm) ,  new_line ('a'),&
!                " gocp*tkh(ku)=",  gocp*tkh(i,j,1:nzm),  new_line ('a'),&
!                " gocp*tkh(kd)=",gocp*tkh(i,j,kd), &
                " tkh(ku)=",  tkh(i,j,1:nzm),  new_line ('a'),&
 !               " tkh(kd)=", tkh(i,j,kd), &
!                " xkzo(ku)=",  xkzo(i,j,1:nzm),  new_line ('a'),&
 !               " xkzo(kd)=", xkzo(i,j,kd), &
                " thermo term=", adzl(i,j,1:nzm)*(dtabsdt(i,j,1:nzm)     -  &
                fac_cond*(dqcldt(i,j,1:nzm) + dqpldt(i,j,1:nzm))            -  &
                fac_sub* (dqcidt(i,j,1:nzm) + dqpidt(i,j,1:nzm))),  new_line ('a')          , &
                " diff term+", gocp*(tkh(i,j,2:nzm) - tkh(i,j,1:nzm-1)),  new_line ('a'), &
!                " dqwvdt=", dqwvdt(i,j,1:nzm),  new_line ('a'),&
                " (smixt(i,j,:)=", smixt(i,j,1:nzm)
!                " rho=", prsl(i,j,k)/(rgas*tabs(i,j,k)*(1+0.622*total_water(i,j,k))), &
!                " prsl(i,j,k)=", prsl(i,j,k)," tabs(i,j,k)=", tabs(i,j,k)
           
        endif




!          if  (wthl2 (i,j,k) > 10 ) &
!               if  ((thl_sec(i,j,k) > 1e2) .and. (aterm > atmin_damp .and. aterm < atmax_damp)) &
!             if  ((thl_sec(i,j,k) > 1e2) .and. (aterm < atmin_damp .or. aterm > atmax_damp)) &  
!                  if  ((abs(wthl2(i,j,k)) > 1e1) .and. (aterm < atmin_damp .or. aterm > atmax_damp) .and. abs(Skew_w) < 3. ) & 
!        if  ((abs(wthl2(i,j,k)) > 1e1) .and. (aterm > atmin_damp .and. aterm < atmax_damp)) &  

!        if  ((thl_sec(i,j,k) > 2e2)) &

             if  ((tabs(i,j,k) < 150 )) & !.and. k==1)) &

!             if  ((wthv_sec(i,j,k) > 1.95e-2)) & !.and. k==1)) &

!        if  (Cek(i,j,k)>5) &

               print *, "aterm=", aterm, "  w1_1=",w1_1,"  w_first=", w_first, "  w2_1=",w2_1, &
               "  sqrtw2=",sqrtw2,"  sqrt(w_sec)=",sqrt(w_sec(i,j,k)) ,  "  w_sec=",w_sec(i,j,k) ,&
               "  qw1_1=",qw1_1,"  qw_first=", qw_first," qw2_1=",  qw2_1, "  qw_sec=", sqrtqt, &
               "  w1_2=",w1_2,"  w_first=", w_first, "  w2_2=",w2_2,"  w_sec=",sqrtw2, &
               "  qw1_2=",qw1_2,"  qw_first=", qw_first," qw2_2=",  qw2_2, "  qw_sec=", sqrtqt, "  wqw2=", wqw2(i,j,k), &
               "  thl1_1=",thl1_1,"  thl_first=", thl_first," thl2_1=",  thl2_1, "  thl_sec=", sqrtthl, &
               "  thl1_2=",thl1_2," thl2_2=",  thl2_2, &
               " w3var=", w3var, " Skew_w=",Skew_w," Skew_hl=",Skew_hl, " Skew_qw=",Skew_qw,i,k, &
               "wrk=", wrk,"(one - w2_1)=",(one - w2_1), " wthv_sec(i,j,k)=", wthv_sec(i,j,k),  &
               "  thlsec=", thlsec,  "  wthlsec=", wthlsec,"  qwsec=", qwsec,  "  wqwsec=", wqwsec, &
               "  bastoeps =", onebeps * thv(i,j,k), " wrk =", epsv * thv(i,j,k), &
               " wqls=", wqls, " ql1", ql1,  " ql2", ql2, &
               " wqis=", wqls, " qi1", ql1,  " qi2", ql2, &
               "  w3var=",w3var, "  w_sec=",w_sec(i,j,k) , &
               "  wrk*w3var=",wrk*w3var, " beta_factor*oneb3*thlsec=", beta_factor*oneb3*thlsec, &
               " wrk*(one-beta_factor*oneb3)*wthlsec * wqwsec=", wrk*(one-beta_factor*oneb3)*wthlsec * wthlsec, &
               " wthl2 estimate LG05(i,j,k)=", wrk*w3var*(beta_factor*oneb3*thlsec + & 
                wrk*(one-beta_factor*oneb3)*wthlsec*wthlsec)  , &
               " wthl2 estimate GLC02(i,j,k)=", &
               aterm*(w1_1 - w_first)*((thl1_1 - thl_first)**2 + thl2_1) &
               + onema*(w1_2 - w_first)*((thl1_2 - thl_first)**2 + thl2_2), &
               " c_w_thl=", c_w_thl, " sqrtw2t=", sqrtw2t, &
 !              " (one - (c_w_thl/sqrtw2t)**2)/aterm=", (one - (c_w_thl/sqrtw2t)**2)/aterm, &
!               " (one - (c_w_thl/sqrtw2t)**2)/onema=", (one - (c_w_thl/sqrtw2t)**2)/onema, &
               " (beta_factor*oneb3 + aterm*(1-twoby3*beta_factor))=", (beta_factor*oneb3 + aterm*(1-twoby3*beta_factor)), &
               " Cek(i,j,k)=", Cek(i,j,k), " Ce  * Cefac=", Ce  * Cefac, &
               gamaz(i,j,k), fac_cond*(diag_ql+qpl(i,j,k)), &
                                                  fac_sub *(diag_qi+qpi(i,j,k)), &
                      tkesbdiss(i,j,k) * (dtn/cp), tke(i,j,k)

               


          if  ( k == 2 .and. Skew_w > 8. .and. thl_sec(i,j,k) > 2e2) print *, "avew=", avew_save(i,j,k-1:k+1),"  cond_w=",cond_w_save(i,j,k-1:k+1), &
              "  z=", z_save(i,j,k-1:k+1), & 
               "  w3=",w3_save(i,j,k-1:k+1),"  w3var=",w3var_save(i,j,k-1:k+1), "  w_sec=",w_sec_save(i,j,k-1:k+1), &
               "  Skew_w=",Skew_w_save(i,j,k-1:k+1)  

!          if  ( wqw2 (i,j,k) < -10 )  &
!               print *, "wrk=", wrk,"(one - w2_1)=",(one - w2_1), &
!              "  thlsec=", thlsec,  "  wthlsec=", wthlsec,"  qwsec=", qwsec,  "  wqwlsec=", wqwsec, &
!               "  w3var=",w3var, "  w_sec=",w_sec(i,j,k) , &
!               "  wrk*w3var=",wrk*w3var, " beta_factor*oneb3*qwsec=", beta_factor*oneb3*qwsec, &
!               " wrk*(one-beta_factor*oneb3)*wqwsec * wqwsec=", wrk*(one-beta_factor*oneb3)*wqwsec * wqwsec, &
!               " wqw2 (i,j,k)=", wqw2(i,j,k)


!          if  ( wthl2 (i,j,k) > 10 )  &
!               print *, "wrk=", wrk,"(one - w2_1)=",(one - w2_1), &
!              "  thlsec=", thlsec,  "  wthlsec=", wthlsec,"  qwsec=", qwsec,  "  wqwlsec=", wqwsec, &
!               "  w3var=",w3var, "  w_sec=",w_sec(i,j,k) , &
!               "  wrk*w3var=",wrk*w3var, " beta_factor*oneb3*thlsec=", beta_factor*oneb3*thlsec, &
!               " wrk*(one-beta_factor*oneb3)*wthlsec * wqwsec=", wrk*(one-beta_factor*oneb3)*wthlsec * wthlsec, &
!               " wthl2 (i,j,k)=", wthl2(i,j,k)
               


         
        ENDDO  ! i
      ENDDO    ! j
    ENDDO      ! k=1,nzm
    
  end subroutine assumed_pdf


  function interp_center_to_interface(mode, array, zl, zi, value_at_sfc)

! Type of interpolation:
! 0 - Half-sums in the vertical. Used for sanity check. 
! 1 -  monotone piecewise cubic Hermite interpolant
    integer, intent(in)                :: mode
! Location of layer centers in the vertical 
    real, dimension(nzm), intent(in)   :: zl
! Location of layer interfaces in the vertical 
    real, dimension(nz ), intent(in)   :: zi
! Profile given on layer centers
    real, dimension(nzm), intent(in)   :: array
! Optional surface (lowest interface) value of a profile given on layer centers 
    real, optional      , intent(in)   :: value_at_sfc
    real, dimension(nzm)               :: interp_center_to_interface, out

! Derivatives of the interpolant at the data points
    real, dimension(:)  , allocatable  :: array_deriv
    real, dimension(nz)                :: zl_ext, array_ext
    
    integer ierr, k

    select case(mode)

    case(0)

         do k=2,nzm
            out(k) = (array(k) + array(k-1))*half
         end do

         if (present(value_at_sfc)) then 
            out(1)=value_at_sfc
         else
            out(1)=out(2)
         endif

         interp_center_to_interface     = out

      case default
         

         if (present(value_at_sfc)) then

            allocate(array_deriv(nz))

            zl_ext(1)        = 0.
            zl_ext(2:nz)     = zl
            array_ext(1)     = value_at_sfc
            array_ext(2:nz)  = array

! Calculate derivatives of the function at the data points
            call DPCHIM(nz,zl_ext,array_ext,array_deriv, 1,ierr)
! Calculate values of the interpolant at the  evaluation points
            call DPCHFE(nz,zl_ext,array_ext,array_deriv, 1,.false.,nzm,zi(2:nz),out,ierr)
         else

            allocate(array_deriv(nzm))

            ! Calculate derivatives of the function at the data points
            call DPCHIM(nzm,zl,array,array_deriv, 1,ierr)
! Calculate values of the interpolant at the  evaluation points
            call DPCHFE(nzm,zl,array,array_deriv, 1,.false.,nzm,zi(2:nz),out,ierr)
         endif

         deallocate(array_deriv)

         if (present(value_at_sfc)) then

            interp_center_to_interface(1)     = value_at_sfc
            interp_center_to_interface(2:nzm) = out(1:nzm-1)

         else
            interp_center_to_interface(1)     = out(1)
!           interp_center_to_interface = out
            interp_center_to_interface(2:nzm) = out(1:nzm-1)
         endif


    end select


  end function interp_center_to_interface


! Below is the code from the Piecewise Cubic Hermite Interpolant Package (PCHIP)
! PCHIP is a FORTRAN90 library which can construct a piecewise cubic Hermite interpolant 
! to data, and carry out various related operations, by Fred Fritsch, DOE LLNL.
! Reference:
! Fred Fritsch, Ralph Carlson,
! Monotone Piecewise Cubic Interpolation,
! SIAM Journal on Numerical Analysis,
! Volume 17, Number 2, April 1980, pages 238-246.

! This computer code is distributed under the GNU LGPL license.

FUNCTION DPCHST (ARG1, ARG2)

!*****************************************************************************80
!
!! DPCHST: carry out a sign test.
!
!  Discussion:
!
!    The object is to do this without multiplying ARG1*ARG2, to avoid
!    possible over/underflow problems.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!     Returns:
!        -1. if ARG1 and ARG2 are of opposite sign.
!         0. if either argument is zero.
!        +1. if ARG1 and ARG2 are of the same sign.
!
  implicit none

  real ARG1, ARG2
  real dpchst

  DPCHST = SIGN ( 1.0D+00, ARG1 ) * SIGN ( 1.0D+00, ARG2 )

  IF ((ARG1 == 0.0D+00) .OR. (ARG2 == 0.0D+00)) then
    DPCHST = 0.0D+00
  end if

  RETURN
END FUNCTION DPCHST


SUBROUTINE DPCHIM ( N, X, F, D, INCFD, IERR )

!*****************************************************************************80
!
!! DPCHIM sets derivatives for a monotone piecewise cubic Hermite interpolant.
!
!  Discussion:
!
!    Set derivatives needed to determine a monotone piecewise
!    cubic Hermite interpolant to given data.  Boundary values
!    are provided which are compatible with monotonicity.  The
!    interpolant will have an extremum at each point where 
!    monotonicity switches direction.  See DPCHIC if user control
!    is desired over boundary or switch conditions.
!
!    Sets derivatives needed to determine a monotone piecewise cubic
!    Hermite interpolant to the data given in X and F.
!
!    Default boundary conditions are provided which are compatible
!    with monotonicity.  (See DPCHIC if user control of boundary con-
!    ditions is desired.)
!
!    If the data are only piecewise monotonic, the interpolant will
!    have an extremum at each point where monotonicity switches direc-
!    tion.  (See DPCHIC if user control is desired in such cases.)
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the F- and D-arrays.
!
!    The resulting piecewise cubic Hermite function may be evaluated
!    by DPCHFE or DPCHFD.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Reference:
!
!    Fred Fritsch, Ralph Carlson,
!    Monotone Piecewise Cubic Interpolation,
!    SIAM Journal on Numerical Analysis,
!    Volume 17, Number 2, April 1980, pages 238-246.
!
!    Fred Fritsch, Judy Butland,
!    A Method for Constructing Local Monotone Piecewise
!    Cubic Interpolants,
!    SIAM Journal on Scientific and Statistical Computing,
!    Volume 5, Number 2, June 1984, pages 300-304.
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer    N, IERR
!        real    X(N), F(INCFD,N), D(INCFD,N)
!
!        CALL  DPCHIM (N, X, F, D, INCFD, IERR)
!
!  Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!           If N=2, simply does linear interpolation.
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of dependent variable values to be
!           interpolated.  F(1+(I-1)*INCFD) is value corresponding to
!           X(I).  DPCHIM is designed for monotonic data, but it will
!           work for any F-array.  It will force extrema at points where
!           monotonicity switches direction.  If some other treatment of
!           switch points is desired, DPCHIC should be used instead.
! 
!     D -- (output) real*8 array of derivative values at the data
!           points.  If the data are monotonic, these values will
!           determine a monotone cubic Hermite function.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD.LT.1 .)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that IERR switches in the direction
!                 of monotonicity were detected.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!             (The D-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
  implicit none

  integer    N, INCFD, IERR
  real    X(*), F(INCFD,*), D(INCFD,*)
  integer    I, NLESS1
  real    DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, DSAVE
  real   H1, H2, HSUM, HSUMT3, W1, W2
!  real    DPCHST
!
!  CHECK ARGUMENTS.
!
  IF ( N.LT.2 )  GO TO 5001
  IF ( INCFD.LT.1 )  GO TO 5002
  DO I = 2, N
     IF ( X(I).LE.X(I-1) )  GO TO 5003
  end do
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
  IERR = 0
  NLESS1 = N - 1
  H1 = X(2) - X(1)
  DEL1 = (F(1,2) - F(1,1))/H1
  DSAVE = DEL1
!
!  SPECIAL CASE N=2.  USE LINEAR INTERPOLATION.
!
  IF ( NLESS1 .le. 1) then
    D(1,1) = DEL1
    D(1,N) = DEL1
    return
  end if
!
!  NORMAL CASE  (N .GE. 3).
!
  H2 = X(3) - X(2)
  DEL2 = (F(1,3) - F(1,2))/H2
!
!  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!  SHAPE-PRESERVING.
!
  HSUM = H1 + H2
  W1 = (H1 + HSUM)/HSUM
  W2 = -H1/HSUM
  D(1,1) = W1*DEL1 + W2*DEL2
  IF ( DPCHST(D(1,1),DEL1) .LE. 0.0D+00 )  THEN
     D(1,1) = 0.0D+00
  ELSE IF ( DPCHST(DEL1,DEL2) .LT. 0.0D+00 )  THEN
!
!  NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
!
     DMAX = 3.0D+00 *DEL1
     IF (ABS(D(1,1)) .GT. ABS(DMAX))  D(1,1) = DMAX
  end if
!
!  LOOP THROUGH INTERIOR POINTS.
!
  DO I = 2, NLESS1

     IF ( 2 < I ) then
       H1 = H2
       H2 = X(I+1) - X(I)
       HSUM = H1 + H2
       DEL1 = DEL2
       DEL2 = (F(1,I+1) - F(1,I))/H2
     end if
!
!  SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
!
     D(1,I) = 0.0D+00
     IF ( DPCHST(DEL1,DEL2) )  42, 41, 45
!
!  COUNT NUMBER OF CHANGES IN DIRECTION OF MONOTONICITY.
!
   41    CONTINUE
     IF (DEL2 == 0.0D+00 )  GO TO 50
     IF ( DPCHST(DSAVE,DEL2) .LT. 0.0D+00 )  IERR = IERR + 1
     DSAVE = DEL2
     GO TO 50

   42    CONTINUE
     IERR = IERR + 1
     DSAVE = DEL2
     GO TO 50
!
!  USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
!
   45    CONTINUE
     HSUMT3 = HSUM+HSUM+HSUM
     W1 = (HSUM + H1)/HSUMT3
     W2 = (HSUM + H2)/HSUMT3
     DMAX = MAX( ABS(DEL1), ABS(DEL2) )
     DMIN = MIN( ABS(DEL1), ABS(DEL2) )
     DRAT1 = DEL1/DMAX
     DRAT2 = DEL2/DMAX
     D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)

   50    CONTINUE

  end do
!
!  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!  SHAPE-PRESERVING.
!
  W1 = -H2/HSUM
  W2 = (H2 + HSUM)/HSUM
  D(1,N) = W1*DEL1 + W2*DEL2
  IF ( DPCHST(D(1,N),DEL2) .LE. 0.0D+00 )  THEN
     D(1,N) = 0.0D+00
  ELSE IF ( DPCHST(DEL1,DEL2) .LT. 0.0D+00 )  THEN
!
!  NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
!
     DMAX = 3.0D+00 *DEL2
     IF (ABS(D(1,N)) .GT. ABS(DMAX))  D(1,N) = DMAX
  end if

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
  IERR = -1
!  CALL XERMSG ('SLATEC', 'DPCHIM', &
!    'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  RETURN

 5002 CONTINUE
!     INCFD.LT.1 RETURN.
  IERR = -2
!  CALL XERMSG ('SLATEC', 'DPCHIM', 'INCREMENT LESS THAN ONE', IERR, 1)
!  RETURN

 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
!  CALL XERMSG ('SLATEC', 'DPCHIM', &
!    'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
  RETURN
END SUBROUTINE DPCHIM


SUBROUTINE DPCHFE ( N, X, F, D, INCFD, SKIP, NE, XE, FE, IERR )

!*****************************************************************************80
!
!! DPCHFE evaluates a piecewise cubic Hermite function at many points.
!
!  Discussion:
!
!    DPCHFE evaluates a piecewise cubic Hermite function at an array of
!    points.  It may be used by itself for Hermite interpolation,
!    or as an evaluator for DPCHIM or DPCHIC.
!
!    Evaluates the cubic Hermite function defined by  N, X, F, D  at
!    the points  XE(J), J=1(1)NE.
!
!    To provide compatibility with DPCHIM and DPCHIC, includes an
!    increment between successive values of the F and D arrays.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        integer    N, NE, IERR
!        real    X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE)
!        LOGICAL  SKIP
!
!        CALL  DPCHFE (N, X, F, D, INCFD, SKIP, NE, XE, FE, IERR)
!
!  Parameters:
!
!     N -- (input) number of data points.  (Error return if N.LT.2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) .LT. X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
!           is the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD.LT.1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in DPCHIM or DPCHIC).
!           SKIP will be set to .TRUE. on normal return.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE.LT.1 .)
!
!     XE -- (input) real*8 array of points at which the function is to
!           be evaluated.
!
!          NOTES:
!           1. The evaluation will be most efficient if the elements
!              of XE are increasing relative to X;
!              that is,   XE(J) .GE. X(I)
!              implies    XE(K) .GE. X(I),  all K.GE.J .
!           2. If any of the XE are outside the interval [X(1),X(N)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!     FE -- (output) real*8 array of values of the cubic Hermite
!           function defined by  N, X, F, D  at the points  XE.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR.GT.0  means that extrapolation was performed at
!                 IERR points.
!           "Recoverable" errors:
!              IERR = -1  if N.LT.2 .
!              IERR = -2  if INCFD.LT.1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if NE.LT.1 .
!             (The FE-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
  implicit none

  integer   incfd

  integer    N, NE, IERR
  real    X(*), F(INCFD,*), D(INCFD,*), XE(*), FE(*)
  LOGICAL  SKIP
  integer    I, IERC, IR, J, JFIRST, NEXT(2), NJ
!
!  CHECK ARGUMENTS.
!
  IF ( .not. SKIP ) then

    IF ( N.LT.2 )  GO TO 5001
    IF ( INCFD.LT.1 )  GO TO 5002
    DO I = 2, N
      IF ( X(I).LE.X(I-1) )  GO TO 5003
    end do

  end if

!  skip = .true.
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
  IF ( NE.LT.1 )  GO TO 5004
  IERR = 0
!
!  LOOP OVER INTERVALS.
!  INTERVAL INDEX IS  IL = IR-1  .
!  INTERVAL IS X(IL).LE.X.LT.X(IR) .
!
  JFIRST = 1
  IR = 2
   10 CONTINUE
!
!  SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.
!
     IF (JFIRST .GT. NE) then
       return
     end if
!
!  LOCATE ALL POINTS IN INTERVAL.
!
     DO J = JFIRST, NE
        IF (XE(J) .GE. X(IR))  GO TO 30
     end do

     J = NE + 1
     GO TO 40
!
!  HAVE LOCATED FIRST POINT BEYOND INTERVAL.
!
   30    CONTINUE
     IF (IR == N)  J = NE + 1

   40    CONTINUE
     NJ = J - JFIRST
!
!  SKIP EVALUATION IF NO POINTS IN INTERVAL.
!
     IF (NJ == 0)  GO TO 50
!
!  EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
!
    CALL DCHFEV (X(IR-1),X(IR), F(1,IR-1),F(1,IR), D(1,IR-1),D(1,IR), &
      NJ, XE(JFIRST), FE(JFIRST), NEXT, IERC)

     IF (IERC .LT. 0)  GO TO 5005

     IF (NEXT(2) == 0)  GO TO 42
!        IF (NEXT(2) .GT. 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE
!           RIGHT OF X(IR).
!
        IF (IR .LT. N)  GO TO 41
!           IF (IR == N)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
           IERR = IERR + NEXT(2)
           GO TO 42
   41       CONTINUE
!
!  ELSE WE SHOULD NEVER HAVE GOTTEN HERE.
!
           GO TO 5005
!           end if
!        end if
   42    CONTINUE
!
     IF (NEXT(1) == 0)  GO TO 49
!        IF (NEXT(1) .GT. 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
!           LEFT OF X(IR-1).
!
        IF (IR .GT. 2)  GO TO 43
!           IF (IR == 2)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
           IERR = IERR + NEXT(1)
           GO TO 49
   43       CONTINUE
!           ELSE
!              XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
!              EVALUATION INTERVAL.
!
!  FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
!
           DO I = JFIRST, J-1
              IF (XE(I) .LT. X(IR-1))  GO TO 45
           end do
!
!  CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR IN DCHFEV.
!
           GO TO 5005

   45          CONTINUE
!
!  RESET J.  THIS WILL BE THE NEW JFIRST.
!
           J = I
!
!  NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
!
           DO I = 1, IR-1
              IF (XE(J) .LT. X(I)) GO TO 47
           end do
!
!  CAN NEVER DROP THROUGH HERE, SINCE XE(J).LT.X(IR-1).
!
   47          CONTINUE
!
!  AT THIS POINT, EITHER  XE(J) .LT. X(1)
!  OR      X(I-1) .LE. XE(J) .LT. X(I) .
!  RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE CYCLING.
!
           IR = MAX(1, I-1)
!           end if
!        end if
   49    CONTINUE

     JFIRST = J

   50 CONTINUE
  IR = IR + 1
  IF (IR .LE. N)  GO TO 10

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N.LT.2 RETURN.
  IERR = -1
!  CALL XERMSG ('SLATEC', 'DPCHFE', &
!    'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  RETURN

 5002 CONTINUE
!
!  INCFD.LT.1 RETURN.
!
  IERR = -2
!  CALL XERMSG ('SLATEC', 'DPCHFE', 'INCREMENT LESS THAN ONE', IERR, 1)
  RETURN
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
!  CALL XERMSG ('SLATEC', 'DPCHFE', 'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
  RETURN
!
 5004 CONTINUE
!     NE.LT.1 RETURN.
  IERR = -4
!  CALL XERMSG ('SLATEC', 'DPCHFE', &
!    'NUMBER OF EVALUATION POINTS LESS THAN ONE', IERR, 1)
  RETURN
!
 5005 CONTINUE
!     ERROR RETURN FROM DCHFEV.
!   *** THIS CASE SHOULD NEVER OCCUR ***
  IERR = -5
!  CALL XERMSG ('SLATEC', 'DPCHFE', &
!    'ERROR RETURN FROM DCHFEV -- FATAL', IERR, 2)
  RETURN
END


SUBROUTINE DCHFEV ( X1, X2, F1, F2, D1, D2, NE, XE, FE, NEXT, IERR )

!*****************************************************************************80
!
!! DCHFEV evaluates a cubic Hermite polynomial at many points.
!
!  Discussion:
!
!    Evaluate a cubic polynomial given in Hermite form at an
!    array of points.  
!
!    Evaluates the cubic polynomial determined by function values
!    F1,F2 and derivatives D1,D2 on interval (X1,X2) at the points
!    XE(J), J=1(1)NE.
!
!    While designed for use by DPCHFE, it may
!    be useful directly as an evaluator for a piecewise cubic
!    Hermite function in applications, such as graphing, where
!    the interval is known in advance.
!
!  Modified:
!
!    05 April 2015
!
!  Author:
!
!    Fred Fritsch
!
!  Calling sequence:
!
!        integer    NE, NEXT(2), IERR
!        real    X1, X2, F1, F2, D1, D2, XE(NE), FE(NE)
!
!        CALL  DCHFEV (X1,X2, F1,F2, D1,D2, NE, XE, FE, NEXT, IERR)
!
!  Parameters:
!
!     X1,X2 -- (input) endpoints of interval of definition of cubic.
!           (Error return if  X1 == X2 .)
!
!     F1,F2 -- (input) values of function at X1 and X2, respectively.
!
!     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE.LT.1 .)
!
!     XE -- (input) real*8 array of points at which the function is to
!           be evaluated.  If any of the XE are outside the interval
!           [X1,X2], a warning error is returned in NEXT.
!
!     FE -- (output) real*8 array of values of the cubic function
!           defined by  X1,X2, F1,F2, D1,D2  at the points  XE.
!
!     NEXT -- (output) integer   array indicating number of 
!     extrapolation points:
!            NEXT(1) = number of evaluation points to left of interval.
!            NEXT(2) = number of evaluation points to right of interval.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if NE.LT.1 .
!              IERR = -2  if X1 == X2 .
!                (The FE-array has not been changed in either case.)
!
  implicit none

  integer    NE, NEXT(2), IERR
  real    X1, X2, F1, F2, D1, D2, XE(*), FE(*)
  integer    I
  real    C2, C3, DEL1, DEL2, DELTA, H, X, XMI, XMA
!
!  CHECK ARGUMENTS.
!
  IF (NE .LT. 1)  GO TO 5001
  H = X2 - X1
  IF (H == 0.0D+00 )  GO TO 5002
!
!  INITIALIZE.
!
  IERR = 0
  NEXT(1) = 0
  NEXT(2) = 0
  XMI = MIN ( 0.0D+00, H)
  XMA = MAX ( 0.0D+00, H)
!
!  COMPUTE CUBIC COEFFICIENTS EXPANDED ABOUT X1.
!
  DELTA = (F2 - F1)/H
  DEL1 = (D1 - DELTA)/H
  DEL2 = (D2 - DELTA)/H
!
!  DELTA IS NO LONGER NEEDED.
!
  C2 = -(DEL1+DEL1 + DEL2)
  C3 = (DEL1 + DEL2)/H
!
!  H, DEL1 AND DEL2 ARE NO LONGER NEEDED.
!
!  EVALUATION LOOP.
!
  DO I = 1, NE
     X = XE(I) - X1
     FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
!
!  COUNT EXTRAPOLATION POINTS.
!
     IF ( X.LT.XMI )  NEXT(1) = NEXT(1) + 1
     IF ( X.GT.XMA )  NEXT(2) = NEXT(2) + 1
!
!  NOTE REDUNDANCY.  IF EITHER CONDITION IS TRUE, OTHER IS FALSE.
!
  end do

  RETURN
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!
!  NE.LT.1 RETURN.
!
  IERR = -1
!  CALL XERMSG ('SLATEC', 'DCHFEV', &
!   'NUMBER OF EVALUATION POINTS LESS THAN ONE', IERR, 1)
  RETURN

 5002 CONTINUE
!
!  X1 == X2 RETURN.
!
  IERR = -2
!  CALL XERMSG ('SLATEC', 'DCHFEV', 'INTERVAL ENDPOINTS EQUAL', IERR, 1)

  RETURN
END SUBROUTINE DCHFEV


! Saturation vapor pressure and mixing ratio subroutines
! Based on Flatau et al (1992), J. App. Met., 31, 1507-1513
! https://doi.org/10.1175/1520-0450(1992)031<1507:PFTSVP>2.0.CO;2
! Code by Marat Khairoutdinov
! Currenty SHOC uses saturation vapor pressure subroutines from GFS, 
! so the code below is unused. 
 

  real function esatw(t)
    real t	! temperature (K)
    real a0,a1,a2,a3,a4,a5,a6,a7,a8 
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /                       &
         6.11239921,       0.443987641,     0.142986287e-1, &
         0.264847430e-3,   0.302950461e-5,  0.206739458e-7, &
         0.640689451e-10, -0.952447341e-13,-0.976195544e-15/
    real dt
    dt    = max(-80.,t-273.16)
    esatw = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
  end function esatw

  real function qsatw(t,p)
    real t	! temperature (K)
    real p	! pressure    (Pa)
    real esat
!   esat  = fpvs(t)
    esat  = fpvsl(t)
    qsatw = 0.622 * esat/max(esat,p-0.378*esat) 
!   esat  = esatw(t)
!   qsatw = 0.622 * esat/max(esat,p-esat) 
  end function qsatw
  
  
  real function esati(t)
    real t	! temperature (K)
    real a0,a1,a2,a3,a4,a5,a6,a7,a8 
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /                     &
         6.11147274,     0.503160820,     0.188439774e-1, &
         0.420895665e-3, 0.615021634e-5,  0.602588177e-7, &
         0.385852041e-9, 0.146898966e-11, 0.252751365e-14/
    real dt
!    real esatw
    if(t > 273.15) then
       esati = esatw(t)
    else if(t.gt.185.) then
       dt    = t-273.16
       esati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
    else   ! use some additional interpolation below 184K
       dt    = max(-100.,t-273.16)
       esati = 0.00763685 + dt*(0.000151069+dt*7.48215e-07)
    endif
  end function esati
        
  real function qsati(t,p)
    real t	! temperature (K)
    real p	! pressure    (Pa)
    real esat !,esati
!   esat  = fpvs(t)
    esat  = fpvsi(t)
    qsati = 0.622 * esat/max(esat,p-0.378*esat)
!   esat  = esati(t)
!   qsati = 0.622 * esat/max(esat,p-esat)
  end function qsati
  
  real function dtesatw(t)
    real t	! temperature (K)
    real a0,a1,a2,a3,a4,a5,a6,a7,a8 
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /                        &
         0.443956472,      0.285976452e-1,   0.794747212e-3, &
         0.121167162e-4,   0.103167413e-6,   0.385208005e-9, &
        -0.604119582e-12, -0.792933209e-14, -0.599634321e-17/
    real dt
    dt      = max(-80.,t-273.16)
    dtesatw = a0 + dt* (a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
  end function dtesatw
        
  real function dtqsatw(t,p)
    real t	! temperature (K)
    real p	! pressure    (Pa)
!    real dtesatw
    dtqsatw = 100.0*0.622*dtesatw(t)/p
  end function dtqsatw
  
  real function dtesati(t)
    real t	! temperature (K)
    real a0,a1,a2,a3,a4,a5,a6,a7,a8 
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /                      &
         0.503223089,     0.377174432e-1,  0.126710138e-2, &
         0.249065913e-4,  0.312668753e-6,  0.255653718e-8, &
         0.132073448e-10, 0.390204672e-13, 0.497275778e-16/
    real dt
!    real dtesatw
    if(t > 273.15) then
       dtesati = dtesatw(t)
    else if(t > 185.) then
       dt      = t-273.16
       dtesati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
    else  ! use additional interpolation below 185K
       dt      = max(-100.,t-273.16)
       dtesati = 0.0013186 + dt*(2.60269e-05+dt*1.28676e-07)
    endif
  end function dtesati
  
  
  real function dtqsati(t,p)
    real t	! temperature (K)
    real p	! pressure    (Pa)
!    real dtesati
    dtqsati = 100.0*0.622*dtesati(t)/p
  end function dtqsati
  
end subroutine shoc
