!
!  Copyright 2019-2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!=======================================================================
module band_dft_sub
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root
  implicit none

contains

subroutine init_band_dft(system,band)
  use structures, only: s_dft_system,s_band_dft
  implicit none
  type(s_dft_system), intent(inout)  :: system
  type(s_band_dft), intent(inout) ::band

  system%wtk(:) = 0.0d0

  call get_band_kpt( band, system )

  band%num_band_kpt = size( band%band_kpt, 2 )
  !write(*,*) "num_band_kpt=",band%num_band_kpt

  allocate( band%check_conv_esp(band%nref_band,system%nk,system%nspin) )
  band%check_conv_esp=.false.

  if ( comm_is_root(nproc_id_global) ) then
     open(100,file='band.dat')
     write(100,*) "Number_of_Bands:",system%no
     write(100,*) "Number_of_kpt_in_each_block:",system%nk
     write(100,*) "Number_of_blocks:",band%num_band_kpt/system%nk
     close(100) 
 end if


end subroutine init_band_dft

subroutine calc_band_write(iter_band_kpt,system,band,info)
  use structures, only: s_dft_system,s_band_dft,s_parallel_info
  use math_constants, only: pi
  implicit none
  type(s_dft_system), intent(inout)  :: system
  type(s_band_dft), intent(inout) ::band
  type(s_parallel_info),intent(in) :: info
  integer :: iter_band_kpt, ik
  real(8) :: k_red(3)

   band%check_conv_esp=.false.
   ! band%band_kpt is already in Cartesian coordinates (get_band_kpt converts
   ! the reduced input once); copy it into vec_k WITHOUT a second primitive_b
   ! multiplication. Set all nk entries (not only the local ik range) so the
   ! root-rank printout below is also correct under k-distributed MPI.
   do ik=1,system%nk
      system%vec_k(:,ik) = band%band_kpt(:,iter_band_kpt+ik-1)
   end do

   if( comm_is_root(nproc_id_global) ) then
      write(*,10) iter_band_kpt, iter_band_kpt+system%nk-1
      write(*,20) "kpoints (reduced)","kpoints in Cartesian"
      open(100,file='band.dat',position="append")
      do ik=iter_band_kpt,iter_band_kpt+system%nk-1
         ! reduced coordinates: k_red(i) = a_i . k_cart / (2*pi)
         k_red(:) = matmul( transpose(system%primitive_a), band%band_kpt(:,ik) )/(2.0d0*pi)
         write(*,30) ik, k_red(:), system%vec_k(:,ik-iter_band_kpt+1)
         write(100,30) ik, k_red(:), system%vec_k(:,ik-iter_band_kpt+1)
      end do
      close(100)
   end if
10 format(1x,"iter_band_kpt=",i3," to",i3)
20 format(1x,3x,2x,a30,2x,a30)
30 format(1x,i3,2x,3f10.5,2x,3f10.5)

end subroutine calc_band_write

subroutine write_band(system,energy)
  use structures, only: s_dft_system,s_dft_energy
  implicit none
  type(s_dft_system), intent(in) :: system
  type(s_dft_energy), intent(in) :: energy
  integer :: ik, iob, ispin

  if( comm_is_root(nproc_id_global) ) then
     open(100,file='band.dat',position="append")
     do ik=1,size(energy%esp,2)
     do iob=1,size(energy%esp,1)
        write(100,*) ik,iob,(energy%esp(iob,ik,ispin),ispin=1,system%nspin)
     end do
     end do
     close(100)
  end if

end subroutine write_band

subroutine read_bandcalc_param( lattice, nref_band, ndiv_segment, kpt, kpt_label )
  use salmon_global, only: lattice_nml => lattice &
                         , nref_band_nml => nref_band &
                         , num_of_segments_nml => num_of_segments &
                         , ndiv_segment_nml => ndiv_segment &
                         , kpt_nml => kpt &
                         , kpt_label_nml => kpt_label
  implicit none
  character(3),intent(out) :: lattice
  integer,intent(out) :: nref_band
  integer,allocatable,intent(inout) :: ndiv_segment(:)
  real(8),allocatable,intent(inout) :: kpt(:,:)  ! given in reduced coordinates in reciprocal space
  character(1),allocatable,intent(inout) :: kpt_label(:)
  integer,parameter :: unit=100
  integer :: num_of_segments !,i,iformat
  if ( comm_is_root(nproc_id_global) ) then
     write(*,'(a50)') repeat("-",24)//"read_bandcalc_param(start)"
  end if
  !open(unit,file='bandcalc.dat',status='old')
  !read(unit,*) lattice; write(*,*) lattice
  !read(unit,*) nref_band
  lattice=lattice_nml
  nref_band=nref_band_nml
  if ( comm_is_root(nproc_id_global) ) write(*,*) "nref_band =",nref_band
  if ( lattice == "non" ) then
  else
     if ( comm_is_root(nproc_id_global) ) then
        write(*,'(a50)') repeat("-",23)//"read_bandcalc_param(return)"
     end if
     return
  end if
  !read(unit,*) num_of_segments
  num_of_segments=num_of_segments_nml
  allocate( ndiv_segment(num_of_segments) ); ndiv_segment=0
  allocate( kpt(3,num_of_segments+1)      ); kpt=0.0d0
  allocate( kpt_label(num_of_segments+1)  ); kpt_label=""
  !read(unit,*) ndiv_segment(:)
  ! A path of N segments has N+1 end points: slice up to num_of_segments+1
  ndiv_segment=ndiv_segment_nml(1:num_of_segments)
  kpt=kpt_nml(:,1:num_of_segments+1)
  kpt_label=kpt_label_nml(1:num_of_segments+1)
  !call check_data_format( unit, iformat )
  !select case( iformat )
  !case( 0 )
  !   do i=1,num_of_segments+1
  !      read(unit,*) kpt(1:3,i)
  !   end do
  !case( 1 )
  !   do i=1,num_of_segments+1
  !      read(unit,*) kpt_label(i), kpt(1:3,i)
  !   end do
  !end select
  !close(unit)
  if ( comm_is_root(nproc_id_global) ) then
     write(*,'(a50)') repeat("-",26)//"read_bandcalc_param(end)"
  end if
end subroutine read_bandcalc_param

subroutine check_data_format( unit, iformat )
  implicit none
  integer,intent(in) :: unit
  integer,intent(out) :: iformat
  character(100) :: ccc
  character(1) :: b(4)
  read(unit,'(a)') ccc
  backspace(unit)
  read(ccc,*,END=9) b
  iformat=1 ! 4 data in one line
  return
9 iformat=0 ! 3 data
end subroutine check_data_format

subroutine get_band_kpt( band, system )
   use structures, only: s_dft_system, s_band_dft
   use parallelization, only: nproc_id_global
   use communication, only: comm_is_root
   implicit none
   type(s_band_dft),intent(inout) :: band  ! fills band_kpt, nref_band + path metadata
   type(s_dft_system),intent(in) :: system
   real(8),allocatable :: kpt(:,:)
   integer :: nref_band ! convergence is checked up to nref_band
   real(8) :: G(3),X(3),M(3),R(3),L(3),W(3) ! XYZ coordinates of high-symmetry
   real(8) :: H(3),N(3),P(3),A(3),Q(3)      ! points in the 1st Brillouin zone
   real(8) :: al,cl ! length of the real-space lattice vectors (a- and c-axis)
   real(8) :: dk(3),k0(3),k1(3),pi,c1,c2,c3
   character(3) :: lattice
   integer,allocatable :: ndiv_segment(:)
   real(8),allocatable :: kpt_(:,:)
   character(1),allocatable :: kpt_label(:)
   integer :: nk,nnk,iseg,num_of_segments,i,ik

   if ( comm_is_root(nproc_id_global) ) then
      write(*,'(a60)') repeat("-",41)//"get_band_kpt(start)"
   end if

   pi=acos(-1.0d0)

   call read_bandcalc_param( lattice, nref_band, ndiv_segment, kpt_, kpt_label )

   if ( allocated(ndiv_segment) ) then
      if ( comm_is_root(nproc_id_global) ) then
         write(*,*) "k points are generated from 'bandcalc.dat'"
      end if
      do ik=1,size(kpt_,2)
         k0(:) = matmul( system%primitive_b, kpt_(:,ik) )
         kpt_(:,ik) = k0(:)
      end do
      num_of_segments = size( ndiv_segment )
   else if ( .not.allocated(ndiv_segment) ) then ! set default
      if ( comm_is_root(nproc_id_global) ) then
         write(*,*) "k points are generated by a default setting"
      end if
      select case( lattice )
      case( "sc" , "SC"  ); num_of_segments=5
      case( "fcc", "FCC" ); num_of_segments=5
      case( "bcc", "BCC" ); num_of_segments=5
      case( "hex", "HEX" ); num_of_segments=7
      case default
         write(*,*) "lattice=",lattice
         write(*,*)"default setting is not available for this lattice" 
         stop "stop@get_band_kpt"
      end select
      allocate( ndiv_segment(num_of_segments) ); ndiv_segment=10
      allocate( kpt_(3,num_of_segments+1)     ); kpt_=0.0d0
      allocate( kpt_label(num_of_segments+1)  ); kpt_label=""
      select case( lattice )
      case( "sc" , "SC"  ) ! G -> X -> M -> R -> G -> M  (5 segments)
         al=sqrt(sum(system%primitive_a(:,1)**2))
         c1=2.0d0*pi/al
         G=c1*(/ 0.0d0, 0.0d0, 0.0d0 /)
         X=c1*(/ 0.5d0, 0.0d0, 0.0d0 /)
         M=c1*(/ 0.5d0, 0.5d0, 0.0d0 /)
         R=c1*(/ 0.5d0, 0.5d0, 0.5d0 /)
         kpt_(:,1)=G; kpt_label(1)="G"
         kpt_(:,2)=X; kpt_label(2)="X"
         kpt_(:,3)=M; kpt_label(3)="M"
         kpt_(:,4)=R; kpt_label(4)="R"
         kpt_(:,5)=G; kpt_label(5)="G"
         kpt_(:,6)=M; kpt_label(6)="M"
      case( "fcc", "FCC" ) ! G -> X -> W -> G -> L -> X  (5 segments)
         al=sqrt(sum(system%primitive_a(:,1)**2))*sqrt(2.0d0)
         c1=2.0d0*pi/al
         G=c1*(/ 0.0d0, 0.0d0, 0.0d0 /)
         X=c1*(/ 1.0d0, 0.0d0, 0.0d0 /)
         W=c1*(/ 1.0d0, 0.5d0, 0.0d0 /)
         L=c1*(/ 0.5d0, 0.5d0, 0.5d0 /)
         kpt_(:,1)=G; kpt_label(1)="G"
         kpt_(:,2)=X; kpt_label(2)="X"
         kpt_(:,3)=W; kpt_label(3)="W"
         kpt_(:,4)=G; kpt_label(4)="G"
         kpt_(:,5)=L; kpt_label(5)="L"
         kpt_(:,6)=X; kpt_label(6)="X"
      case( "bcc", "BCC" ) ! G -> H -> N -> P -> G -> N  (5 segments)
         al=sqrt(sum(system%primitive_a(:,1)**2))*2.0d0/sqrt(3.0d0)
         c1=2.0d0*pi/al
         G=c1*(/ 0.0d0, 0.0d0, 0.0d0 /)
         H=c1*(/ 0.0d0, 1.0d0, 0.0d0 /)
         N=c1*(/ 0.5d0, 0.5d0, 0.0d0 /)
         P=c1*(/ 0.5d0, 0.5d0, 0.5d0 /)
         kpt_(:,1)=G; kpt_label(1)="G"
         kpt_(:,2)=H; kpt_label(2)="H"
         kpt_(:,3)=N; kpt_label(3)="N"
         kpt_(:,4)=P; kpt_label(4)="P"
         kpt_(:,5)=G; kpt_label(5)="G"
         kpt_(:,6)=N; kpt_label(6)="N"
      case( "hex", "HEX" ) ! G -> P -> Q -> G -> A -> L -> H -> P  (7 segments)
         al=sqrt(sum(system%primitive_a(:,1)**2))
         cl=sqrt(sum(system%primitive_a(:,3)**2))
         c1=2.0d0*pi/al
         c2=1.0d0*pi/al
         c3=1.0d0*pi/cl
         G=(/ 0.0d0, 0.0d0, 0.0d0 /)
         P=c1*(/ 2.0d0/3.0d0, 0.0d0, 0.0d0 /)
         Q=c2*(/ 1.0d0, 1.0d0/sqrt(3.0d0), 0.0d0 /)
         A=c3*(/ 0.0d0, 0.0d0, 1.0d0 /)
         L=c2*(/ 1.0d0, 1.0d0/sqrt(3.0d0), c3/c2 /)
         H=c1*(/ 2.0d0/3.0d0, 0.0d0, c3/c1 /)
         kpt_(:,1)=G; kpt_label(1)="G"
         kpt_(:,2)=P; kpt_label(2)="P"
         kpt_(:,3)=Q; kpt_label(3)="Q"
         kpt_(:,4)=G; kpt_label(4)="G"
         kpt_(:,5)=A; kpt_label(5)="A"
         kpt_(:,6)=L; kpt_label(6)="L"
         kpt_(:,7)=H; kpt_label(7)="H"
         kpt_(:,8)=P; kpt_label(8)="P"
      end select
   end if

   nk=system%nk
   nnk=sum( ndiv_segment(1:num_of_segments) )
   if ( mod(nnk,nk) /= 0 ) nnk=nnk-mod(nnk,nk)+nk

   allocate( kpt(3,nnk) ); kpt=0.0d0

   i=0
   do iseg=1,num_of_segments
      k0(:)=kpt_(:,iseg)
      k1(:)=kpt_(:,iseg+1)
      dk(:)=( k1(:) - k0(:) )/ndiv_segment(iseg)
      do ik=0,ndiv_segment(iseg)-1
         i=i+1
         kpt(:,i)=k0(:)+dk(:)*ik
      end do
   end do ! iseg

   if ( i < nnk ) then
      i=i+1
      kpt(:,i)=kpt(:,i-1)+dk(:)
   end if
   if ( i < nnk ) then
      do ik=i+1,nnk
         kpt(:,ik)=kpt(:,ik-1)+dk(:)
      end do
   end if

   if ( comm_is_root(nproc_id_global) ) then
      write(*,*) "Number of computed bands:",nref_band
      write(*,*) "Whole number of bands(system%no):",system%no
      write(*,*) "array size of wf for k points(system%nk):",nk
      write(*,*) "Number of segments:",num_of_segments
      write(*,*) "Total number of k points:",nnk
      write(*,*) "k points in Cartesian coordinates:"
      do i=1,size(kpt,2)
         write(*,'(1x,i4,3f10.5)') i,kpt(:,i)
      end do
      write(*,'(a60)') repeat("-",43)//"get_band_kpt(end)"
   end if

   ! --- path metadata for SYSNAME_bandpath.data (write_bandpath_data) --------
   ! cumulative |dk| distance of every path point and of the high-symmetry
   ! nodes (kpt_ holds the node coordinates in Cartesian in both branches).
   if ( allocated(band%kpt_dist)   ) deallocate( band%kpt_dist )
   if ( allocated(band%node_label) ) deallocate( band%node_label )
   if ( allocated(band%node_dist)  ) deallocate( band%node_dist )
   allocate( band%kpt_dist(nnk) )
   band%kpt_dist(1) = 0.0d0
   do i=2,nnk
      band%kpt_dist(i) = band%kpt_dist(i-1) + sqrt(sum( (kpt(:,i)-kpt(:,i-1))**2 ))
   end do
   band%num_of_nodes = num_of_segments + 1
   allocate( band%node_label(band%num_of_nodes), band%node_dist(band%num_of_nodes) )
   band%node_dist(1) = 0.0d0
   do i=1,band%num_of_nodes
      if ( kpt_label(i) == 'G' ) then
         band%node_label(i) = 'Gamma'   ! the plotter's Γ spelling
      else
         band%node_label(i) = kpt_label(i)
      end if
      if ( i >= 2 ) band%node_dist(i) = band%node_dist(i-1) &
                                      + sqrt(sum( (kpt_(:,i)-kpt_(:,i-1))**2 ))
   end do

   band%nref_band = nref_band
   if ( allocated(band%band_kpt) ) deallocate( band%band_kpt )
   call move_alloc( kpt, band%band_kpt )

end subroutine get_band_kpt


!! export SYSNAME_bandpath.data for theory='dft_band': the REAL DFT levels
!! along the high-symmetry path, in the exact format the plotter
!! (plot_sbe_results.py::_load_bandpath) and the SBE spectral tooling read --
!! the same contract as the Python EPM bandpath writers:
!!   # spinor = 0/1
!!   # nv = <valence states per k>
!!   # nb = <states per line>
!!   # nodes: LBL dist  LBL dist ...
!!   rows: ik dist q1 q2 q3 E_1..E_nb   (q reduced, E in Hartree)
!! Called once per Band_Iteration block (system%nk path points each); the
!! header is written with the first block, later blocks append.
subroutine write_bandpath_data(iter_band_kpt,system,band,energy)
  use structures, only: s_dft_system,s_band_dft,s_dft_energy
  use salmon_global, only: base_directory,sysname,nelec,yn_spinorbit
  use filesystem, only: get_filehandle
  use math_constants, only: pi
  implicit none
  integer,intent(in) :: iter_band_kpt
  type(s_dft_system),intent(in) :: system
  type(s_band_dft),intent(in) :: band
  type(s_dft_energy),intent(in) :: energy
  integer :: fh,ik,ikg,ib,nv,i
  real(8) :: k_red(3)
  character(256) :: file_bp

  if ( .not.comm_is_root(nproc_id_global) ) return
  if ( system%nspin /= 1 ) then
     write(*,*) "write_bandpath_data: skipped (contract is nspin=1)"
     return
  end if

  file_bp = trim(base_directory)//trim(sysname)//'_bandpath.data'
  fh = get_filehandle()
  if ( iter_band_kpt == 1 ) then
     open(unit=fh, file=trim(file_bp), status='replace', action='write')
     write(fh,'(a)') '# band path from theory=dft_band (real DFT levels)'
     if ( yn_spinorbit == 'y' ) then
        nv = nelec
        write(fh,'(a)') '# spinor = 1'
     else
        nv = nelec/2
        write(fh,'(a)') '# spinor = 0'
     end if
     write(fh,'(a,i0)') '# nv = ', nv
     write(fh,'(a,i0)') '# nb = ', system%no
     write(fh,'(a,*(2x,a,1x,f0.7))') '# nodes:', &
        & ( trim(band%node_label(i)), band%node_dist(i), i=1,band%num_of_nodes )
     write(fh,'(a)') '# ik, dist, q1, q2, q3 (reduced), E_1..E_nb [Ha]'
  else
     open(unit=fh, file=trim(file_bp), status='old', action='write', position='append')
  end if
  do ik=1,system%nk
     ikg = iter_band_kpt + ik - 1
     ! q_red(i) = a_i . k_cart / (2*pi)
     k_red(:) = matmul( transpose(system%primitive_a), band%band_kpt(:,ikg) )/(2.0d0*pi)
     write(fh,'(i6,f14.7,3f10.5,*(e18.10))') ikg, band%kpt_dist(ikg), k_red(1:3), &
        & ( energy%esp(ib,ik,1), ib=1,system%no )
  end do
  close(fh)

end subroutine write_bandpath_data


end module band_dft_sub
