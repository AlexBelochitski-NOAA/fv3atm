       module neural
            
        use module_iounitdef, only : nonetf
        use machine,          only : kind_phys
        use module_radlw_parameters, only :  topflw_type, sfcflw_type
        use module_radsw_parameters,  only : topfsw_type, sfcfsw_type,    &
     &                                     profsw_type,cmpfsw_type,NBDSW
       
        implicit none 

        private
 
        public :: init_lw_nn_emulator, init_sw_nn_emulator, & 
                  lw_nn_emulation, sw_nn_emulation

! Number of members in the LW NN ensemble

        integer, parameter     :: lw_num_of_members = 1

! Number of members in the SW NN ensemble

        integer, parameter     :: sw_num_of_members = 1

! Files containing LW NN weights and biases

        character(*), parameter::lw_file_name(lw_num_of_members)= &
!             (/'nnrad_lw_database.txt' /)
!
! V.K.    New file with LWR NN paramerers introduced 3/10/2011 V.K.
!             (/'/scratch1/NCEPDEV/global/Alexei.A.Belochitski/data/Lwr187-100-70-24.asc' /)
             (/'/scratch1/NCEPDEV/global/Alexei.A.Belochitski/data/Lwr196-75-69-26.asc' /)

! Files containing SW NN weights and biases

        character(*), parameter::sw_file_name(sw_num_of_members)= &
!             (/'nnrad_sw_database.txt' /)
!
! V.K.    New file with SWR NN paramerers introduced 3/10/2011 V.K.
             (/'/scratch1/NCEPDEV/global/Alexei.A.Belochitski/data/Swr190-100-73-36.asc' /)
!        (/'/scratch1/NCEPDEV/global/Alexei.A.Belochitski/data/Swr201-75-73-22.asc' /)

! Internal types and variables

        type nndata_1d
           real(kind=kind_phys), allocatable :: a(:)
        end type nndata_1d

        type nndata_2d
           real(kind=kind_phys), allocatable :: a(:,:)
        end type nndata_2d

! LW Hidden and output weights
        type(nndata_2d) :: lw_w1(lw_num_of_members),lw_w2(lw_num_of_members)
! LW Hidden and output biases
        type(nndata_1d) :: lw_b1(lw_num_of_members),lw_b2(lw_num_of_members)
! LW Number of inputs, hidden neurons and outputs
        integer      :: lw_in(lw_num_of_members),lw_hid(lw_num_of_members),  & 
                        lw_out(lw_num_of_members)

! SW Hidden and output weights
        type(nndata_2d) :: sw_w1(sw_num_of_members),sw_w2(sw_num_of_members)
! SW Hidden and output biases
        type(nndata_1d) :: sw_b1(sw_num_of_members),sw_b2(sw_num_of_members)
! SW Number of inputs, hidden neurons and outputs
        integer      :: sw_in(sw_num_of_members),sw_hid(sw_num_of_members),  & 
                        sw_out(sw_num_of_members)

! A constant passed over from rlwinit

        real(kind=kind_phys) heatfac

      contains
        
        subroutine init_lw_nn_emulator(me,heatcoeff)
   
          integer, intent(in) ::  me
          real(kind=kind_phys), intent(in) ::  heatcoeff

          integer iin,ihid,iout,member

          heatfac = heatcoeff

          if (me == 0) print*,'Module NEURAL: Number of LW NN ensemble members:', &
                               lw_num_of_members

! Load NN weights and biases

          do member=1,lw_num_of_members

             open(unit=nonetf,err=10,file=trim(lw_file_name(member)),status='old')

             read(nonetf,'(3i5)') iin,ihid,iout
             lw_in(member)=iin; lw_out(member)=iout; lw_hid(member)=ihid 
             
             allocate(lw_w1(member)%a(iin,ihid),lw_w2(member)%a(ihid,iout))
             allocate(lw_b1(member)%a(ihid),lw_b2(member)%a(iout)) 

             read(nonetf,*,err=11,end=12) lw_w1(member)%a, lw_w2(member)%a, &
                                          lw_b1(member)%a, lw_b2(member)%a  
             close(nonetf) 
 
             if (me == 0)  print*,'Module NEURAL: NN LW File Loaded: ', & 
                                   lw_file_name(member),                &
                                  ' iin,ihid,iout: ',iin,ihid,iout
  
          end do

! All NNs in the ensemble must have the same number inputs and outputs 

          if (.not.all(lw_in == lw_in(1)).or..not.all(lw_out == lw_out(1))) then
             if (me == 0) print *,  "Module NEURAL: LW NN ensemble members have different number of inputs and/or outputs. Exiting."
             stop
          endif

!          if (me == 0) print *, "NEURAL IN, HID, OUT:", lw_in(1), lw_hid(1),  lw_out(1)
          return        

! Catch file opening/reading errors

10        if (me == 0) print *, "Module NEURAL: Error opening file ",    & 
                                lw_file_name(member), ". Exiting."
          stop

11        if (me == 0) print *, "Module NEURAL: Error reading file ",    &
                                lw_file_name(member), ". Exiting."
          stop

12        if (me == 0) print *, "Module NEURAL: Reached EOF too early ", &
                                lw_file_name(member), ". Exiting."
          stop 
        
        end subroutine init_lw_nn_emulator


! Initialize SW NN

        subroutine init_sw_nn_emulator(me)
         
          integer, intent(in) ::  me  

          integer iin,ihid,iout,member

          if (me == 0) print*,'Module NEURAL: Number of SW NN ensemble members:', &
                               sw_num_of_members

! Load NN weights and biases

          do member=1,sw_num_of_members

             open(unit=nonetf,err=10,file=trim(sw_file_name(member)),status='old')

             read(nonetf,'(3i5)') iin,ihid,iout
             sw_in(member)=iin; sw_out(member)=iout; sw_hid(member)=ihid 
             
             allocate(sw_w1(member)%a(iin,ihid),sw_w2(member)%a(ihid,iout))
             allocate(sw_b1(member)%a(ihid),sw_b2(member)%a(iout)) 

             read(nonetf,*,err=11,end=12) sw_w1(member)%a, sw_w2(member)%a, &
                                          sw_b1(member)%a, sw_b2(member)%a  
             close(nonetf) 
 
             if (me == 0)  print*,'Module NEURAL: NN SW File Loaded: ', & 
                                   sw_file_name(member),                &
                                  ' iin,ihid,iout: ',iin,ihid,iout
  
          end do

! All NNs in the ensemble must have the same number inputs and outputs 

          if (.not.all(sw_in == sw_in(1)).or..not.all(sw_out == sw_out(1))) then
             if (me == 0) print *,  "Module NEURAL: SW NN ensemble members have different number of inputs and/or outputs. Exiting."
             stop
          endif

!          if (me == 0) print *, "NEURAL IN, HID, OUT:", lw_in(1), lw_hid(1),  lw_out(1)
          return        

! Catch file opening/reading errors

10        if (me == 0) print *, "Module NEURAL: Error opening file ",    & 
                                sw_file_name(member), ". Exiting."
          stop

11        if (me == 0) print *, "Module NEURAL: Error reading file ",    &
                                sw_file_name(member), ". Exiting."
          stop

12        if (me == 0) print *, "Module NEURAL: Reached EOF too early ", &
                                sw_file_name(member), ". Exiting."
          stop 
        

            
        end subroutine init_sw_nn_emulator

  

        subroutine  lw_nn_emulation                                     & 
!  Inputs:
     &     ( pmid,pint,tmid,tint,qnm,o3mr,gasvmr,                       &
     &       clouds,iovr,aerosols,sfemis,                               &
     &       NPTS, NLAY, NLP1, iflip, lprnt, xlon,xlat,hour,month,year, &
! Outputs:
     &       hlwc,topflx,sfcflx)

! NN Emulation of LW RRTM1

! Inputs
          integer,  intent(in) :: NPTS, NLAY
          integer,  intent(in) :: NLP1, iovr, iflip ! Ignored

          real (kind=kind_phys), dimension(:,:), intent(in) :: pint, tint,  &
               &       pmid, tmid, qnm, o3mr

          logical,  intent(in) :: lprnt ! Ignored

          real (kind=kind_phys), dimension(:,:,:), intent(in) :: gasvmr,    &
               &       clouds

          real (kind=kind_phys), dimension(:,:,:,:),intent(in) :: aerosols

          real (kind=kind_phys), dimension(:),     intent(in) :: sfemis

         real (kind=kind_phys),intent(in)::hour,year,month,xlon(:),xlat(:)
          
!  Outputs:
          real (kind=kind_phys), dimension(:,:), intent(out) :: hlwc
          
          type (topflw_type),    dimension(:),   intent(out) :: topflx
          type (sfcflw_type),    dimension(:),   intent(out) :: sfcflx

! Local variables

          real(kind=kind_phys) ::  nn_input_vector(lw_in(1)),  nn_output_vector(lw_out(1))
      
          real(kind=kind_phys),dimension(NPTS)::QRLI,BAL,VESI,EQ,ER
          real(kind=kind_phys),dimension(NPTS,NLAY)::VES,DQ


          integer :: icol, k


          do icol=1,NPTS
!
! V.K.  Interface changed for new GFS radiation by V.K. 2/23/2011
! Old inputs for CFS NN LWR 196:HID:69
             nn_input_vector(1) =    year
             nn_input_vector(2) =    cos(month*3.14159265/6.)
             nn_input_vector(3) =    sin(month*3.14159265/6.) 
             nn_input_vector(4) =    cos(xlon(icol))
             nn_input_vector(5) =    sin(xlon(icol))
             nn_input_vector(6) =    xlat(icol)

             nn_input_vector(7)    =    pint(icol,1)
             nn_input_vector(8:28) =    pint(icol,2:42:2)
             nn_input_vector(29)   =    pint(icol,43)
             nn_input_vector( 30:31) =    tint(icol,1:2)
             nn_input_vector( 32:61) =    tint(icol,3:61:2)
             nn_input_vector( 62:65) =    tint(icol,62:65)
             nn_input_vector(66:105) =    qnm(icol,1:40)
             nn_input_vector(106:139) =    o3mr(icol,31:64)
             nn_input_vector(140:196) =    clouds(icol,1:57,1) ! cfrac

! New inputs for GFS NN LWR 187:HID:70
!
!             nn_input_vector(1)       =    cos(month*3.14159265/6.)
!             nn_input_vector(2)       =    sin(month*3.14159265/6.) 
!             nn_input_vector(3)       =    cos(xlon(icol))
!             nn_input_vector(4)       =    sin(xlon(icol))
!             nn_input_vector(5)       =    xlat(icol)

!             nn_input_vector(6)       =    pint(icol,1)
!             nn_input_vector(7:27)    =    pint(icol,2:42:2)
!             nn_input_vector(28)      =    pint(icol,43)
!             nn_input_vector( 29:30)  =    tint(icol,1:2)
!             nn_input_vector( 31:60)  =    tint(icol,3:61:2)
!             nn_input_vector( 61:64)  =    tint(icol,62:65)
!             nn_input_vector(65:104)  =    qnm(icol,1:40)
!             nn_input_vector(105:138) =    o3mr(icol,31:64)
!             nn_input_vector(139:186) =    clouds(icol,1:48,1) ! cfrac           
!             nn_input_vector(187)     =    sfemis(icol) 

             call compute_nn(nn_input_vector,nn_output_vector,lw_num_of_members,& 
                  lw_w1,lw_w2,lw_b1,lw_b2,lw_hid)

             hlwc(icol,:)       = nn_output_vector(1:64)
 
             topflx(icol)%upfxc = nn_output_vector(65)
             topflx(icol)%upfx0 = nn_output_vector(66)
             sfcflx(icol)%upfxc = nn_output_vector(67)
!
! V.K. --- one output added by V.K. for GFS NN LWR on 2/23/2011 
!
!             sfcflx(icol)%upfx0 = nn_output_vector(68)
!
! V.K. 68 => 69 and 69 => 70 below by V.K. on 2/23/2011
!
!             sfcflx(icol)%dnfxc = nn_output_vector(69) ! (68)
!             sfcflx(icol)%dnfx0 = nn_output_vector(70) ! (69)

             sfcflx(icol)%dnfxc = nn_output_vector(68)
             sfcflx(icol)%dnfx0 = nn_output_vector(69)


          enddo

!! Make sure that heating rates are consitent with fluxes at top/bottom in the
!! sense of energy conservation. Adjust the heating rates accordingly. 
!
!          QRLI=0.; VESI=0. 
!          do k=1,NLAY
!! Convert mbar to Pa
!             VES(:,k)=-(pint(:,k+1)-pint(:,k))/heatfac
!! Convert K/day to W/m**2
!             QRLI=QRLI+hlwc(:,k)*VES(:,k)
!             VESI=VESI+VES(:,k)
!          enddo
!          ER=QRLI + topflx(:)%upfxc - sfcflx(:)%upfxc + sfcflx(:)%dnfxc
!          do k=1,NLAY
!             DQ(:,k)=-ER/VESI
!          enddo
!          hlwc=hlwc+DQ

        end  subroutine  lw_nn_emulation


!        subroutine  sw_nn_emulation  &
!! Inputs:
!     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                      &
!     &       clouds,iovr,aerosols,sfcalb,                               &
!     &       cosz,solcon,NDAY,idxday,                                   &
!     &       NPTS, NLAY, NLP1, iflip, lprnt,xlon,xlat,month,year,       &
!! Outputs:
!     &       hswc,topflx,sfcflx,                                        &
!! Optional
!     &       HSWB,FDNCMP_IN                                             &
!     &     )

        subroutine  sw_nn_emulation  &
! Inputs:                                                                                                    
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,                      &
     &       clouds,sfcalb,                               &
     &       cosz,solcon,                                  &
     &       NPTS, NLAY,xlon,xlat,month,year,       &
! Outputs:                                                                                                   
     &       hswc,topflx,sfcflx,                                        &
! Optional                                                                                                   
     &       FDNCMP_IN                                             &
     &     )
! NN specific inputs
 

!  ---  inputs:
          integer, intent(in) :: NPTS, NLAY !, NLP1, iovr, iflip, NDAY

!          integer, intent(in) :: idxday(:)
!          
!          logical, intent(in) :: lprnt

          real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, tlvl,  &
               &       plyr, tlyr, qlyr, olyr, sfcalb

          real (kind=kind_phys), dimension(:,:,:),   intent(in) :: clouds !,gasvmr

!          real (kind=kind_phys), dimension(:,:,:,:), intent(in) :: aerosols
          
          real (kind=kind_phys), intent(in) :: cosz(:), solcon, year, month,xlon(:),xlat(:)
          
!  ---  outputs:
          real (kind=kind_phys), dimension(:,:), intent(out) :: hswc

          type (topfsw_type),    dimension(:),   intent(out) :: topflx
          type (sfcfsw_type),    dimension(:),   intent(out) :: sfcflx

!! ---  optional outputs:
!          real (kind=kind_phys),dimension(:,:,:),optional,intent(out):: hswb
!          real (kind=kind_phys),dimension(:,:),  optional,intent(out):: hsw0
!          type (profsw_type), dimension(:,:), optional, intent(out) :: flxprf
          type (cmpfsw_type), dimension(:),  optional, intent(out) :: fdncmp_in

! Local variables
          
          real(kind=kind_phys) ::  nn_input_vector(sw_in(1)),  nn_output_vector(sw_out(1))
      
          real(kind=kind_phys),dimension(NPTS)::QRLI,BAL,VESI,EQ,ER
          real(kind=kind_phys),dimension(NPTS,NLAY)::VES,DQ

          logical :: lhswb, lfdncmp

          integer :: icol, k
         
!          lhswb  = present(hswb)
          lfdncmp= present(fdncmp_in)

          do icol=1,NPTS
      
             if (cosz(icol) < 0.0001) then 
 
                 nn_output_vector = 0.
 
             else 
!
! V.K.   Interface changed for new GFS radiation by V.K. 2/23/2011
! Old inputs for CFS NN SWR 201:HID:73
!                nn_input_vector(1)  =  year
!                nn_input_vector(2)  =     cos(month*3.14159265/6.)
!                nn_input_vector(3)  =     sin(month*3.14159265/6.)
!                nn_input_vector(4)  =    cos(xlon(icol))
!                nn_input_vector(5)  =    sin(xlon(icol))
!                nn_input_vector(6:6)  =    xlat(icol)
!                nn_input_vector(7) =      plvl(icol,1)
!                nn_input_vector(8:28) =      plvl(icol,2:42:2)
!                nn_input_vector(29) =      plvl(icol,43)
!                nn_input_vector(30:31) =    tlvl(icol,1:2)
!                nn_input_vector(32:61) =    tlvl(icol,3:61:2)
!                nn_input_vector(62:65) =    tlvl(icol,62:65)
!                nn_input_vector(66:105) =   qlyr(icol,1:40) 
!                nn_input_vector(106:139) =   olyr(icol,31:64) 
!                nn_input_vector(140:195) =   clouds(icol,1:56,1) 
!                nn_input_vector(196:199) =   sfcalb(icol,1:4) 
!                nn_input_vector(200) =       cosz(icol) 
!                nn_input_vector(201)     =   solcon
!
! New inputs for GFS NN SWR 190:HID:73
!                               
                nn_input_vector(1)       =    cos(month*3.14159265/6.)
                nn_input_vector(2)       =    sin(month*3.14159265/6.)
                nn_input_vector(3)       =    cos(xlon(icol))
                nn_input_vector(4)       =    sin(xlon(icol))
                nn_input_vector(5)       =    xlat(icol)

                nn_input_vector(6)       =    plvl(icol,1)
                nn_input_vector(7:27)    =    plvl(icol,2:42:2)
                nn_input_vector(28)      =    plvl(icol,43)
                nn_input_vector(29:30)   =    tlvl(icol,1:2)
                nn_input_vector(31:60)   =    tlvl(icol,3:61:2)
                nn_input_vector(61:64)   =    tlvl(icol,62:65)
                nn_input_vector(65:104)  =    qlyr(icol,1:40) 
                nn_input_vector(105:138) =    olyr(icol,31:64) 
                nn_input_vector(139:185) =    clouds(icol,1:47,1) 
                nn_input_vector(186)     =    cosz(icol) 
                nn_input_vector(187:190) =    sfcalb(icol,1:4) 
                
                call compute_nn(nn_input_vector,nn_output_vector,sw_num_of_members,& 
                     sw_w1,sw_w2,sw_b1,sw_b2,sw_hid)
                 
                where(nn_output_vector < 0.) nn_output_vector=0.
  
             endif
              
             hswc(icol,:)        = nn_output_vector(1:64)
             topflx(icol)%upfxc  = nn_output_vector(65)
! CFS
!             topflx(icol)%dnfxc  = nn_output_vector(66)
!             topflx(icol)%upfx0  = nn_output_vector(67)
!             sfcflx(icol)%upfxc  = nn_output_vector(68)
!             sfcflx(icol)%dnfxc  = nn_output_vector(69)
!             sfcflx(icol)%upfx0  = nn_output_vector(70)
!             sfcflx(icol)%dnfx0  = nn_output_vector(71)

             topflx(icol)%dnfxc  = nn_output_vector(67)
             topflx(icol)%upfx0  = nn_output_vector(66)
             sfcflx(icol)%upfxc  = nn_output_vector(68)
             sfcflx(icol)%dnfxc  = nn_output_vector(70)
             sfcflx(icol)%upfx0  = nn_output_vector(69)
             sfcflx(icol)%dnfx0  = nn_output_vector(71)
             if (lfdncmp) then 
                fdncmp_in(icol)%uvbfc  = nn_output_vector(72)
                fdncmp_in(icol)%uvbf0  = nn_output_vector(73)
             endif
          
          end do

! Make sure that heating rates are consitent with fluxes at top/bottom in the
! sense of energy conservation. Adjust the heating rates accordingly. 

!          QRLI=0.; VESI=0. 
!          do k=1,NLAY
! Convert mbar to Pa
!             VES(:,k)=-(plvl(:,k+1)-plvl(:,k))/heatfac
! Convert K/day to W/m**2
!             QRLI=QRLI+hswc(:,k)*VES(:,k)
!             VESI=VESI+VES(:,k)
!          enddo
!          ER=QRLI + topflx(:)%upfxc - sfcflx(:)%upfxc + sfcflx(:)%dnfxc - topflx(:)%dnfxc
!          do k=1,NLAY
!             DQ(:,k)=-ER/VESI
!          enddo
!          hswc=hswc+DQ

                   
          end  subroutine  sw_nn_emulation


        
        subroutine  compute_nn(X,Y,num_of_members,w1,w2,b1,b2,nhid)

      
 !  Input:
 !            X(IN) NN input vector 
          integer, intent(in) :: num_of_members,nhid(num_of_members)
          real(kind=kind_phys), intent(in)::X(:)
          type(nndata_2d), intent(in) :: w1(num_of_members),w2(num_of_members)
          type(nndata_1d), intent(in) :: b1(num_of_members),b2(num_of_members)
         

 !   Ouput:
 !            Y(OUT) NN output vector (composition coefficients for SNN)

          real(kind=kind_phys), intent(out):: Y(:)

! Local variables 
          integer i, nout
          real(kind=kind_phys), allocatable :: x2(:),x3(:)
          integer member

          nout=size(Y)

          Y = 0.

          allocate(x3(nout)) 

          do member = 1,num_of_members
  
             allocate(x2(nhid(member))) 

! Calculate neurons in the hidden layer

             forall(i = 1:nhid(member)) x2(i)= tanh(sum(X*w1(member)%a(:,i))+  & 
                                               b1(member)%a(i))

! Calculate NN output 

             forall(i=1:nout) x3(i)= sum(w2(member)%a(:,i)*x2) + b2(member)%a(i)

             Y = Y + x3 
             
             deallocate(x2)
         
          end do                    ! member

          deallocate(x3)
      
          Y = Y / num_of_members
         
    end  subroutine  compute_nn
                 

  end module neural




