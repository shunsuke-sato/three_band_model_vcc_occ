!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
!-------10--------20--------30--------40--------50--------60--------70--------80--------90
subroutine current(it,jav)
  use global_variables
  implicit none
  integer,intent(in) :: it
  real(8),intent(out) :: jav
  real(8) :: jav1,jav2,jz_intra,jz_inter,jav_l
  real(8) :: mass_c1_i,mass_c2_i,mass_h_i
  integer :: ik,ikr,ikz

  kz(:) = kz0(:) + Act(it)*fact_intra
  mass_h_i = 1d0/mass_h
  mass_c1_i = 1d0/mass_c1
  mass_c2_i = 1d0/mass_c2
  jav_l = 0d0
  do ik = NKrz_s,NKrz_e
    ikr = ikr_table(ik); ikz = ikz_table(ik)

! state 1
    jz_intra = kz(ikz)*( &
        (occ(1,ik)*abs(zCt(2,1,ik))**2 &
        +occ(2,ik)*abs(zCt(2,2,ik))**2 &
        +occ(3,ik)*abs(zCt(2,3,ik))**2 &   
        )*mass_c1_i &
      +&
      (occ(1,ik)*abs(zCt(3,1,ik))**2 &
      +occ(2,ik)*abs(zCt(3,2,ik))**2 &
      +occ(3,ik)*abs(zCt(3,3,ik))**2 &
      )*mass_c2_i &
      +&
      (occ(1,ik)*abs(zCt(1,1,ik))**2 &
      +occ(2,ik)*abs(zCt(1,2,ik))**2 &
      +occ(3,ik)*abs(zCt(1,3,ik))**2 &
      )*mass_h_i )

    jz_inter = 2d0*occ(1,ik)*real( &
      piz_dc1*conjg(zCt(1,1,ik))*zCt(2,1,ik) &
     +piz_dc2*conjg(zCt(1,1,ik))*zCt(3,1,ik) &
     +piz_dcc*conjg(zCt(2,1,ik))*zCt(3,1,ik) &
      ) &
      +2d0*occ(2,ik)*real( &
      piz_dc1*conjg(zCt(1,2,ik))*zCt(2,2,ik) &
     +piz_dc2*conjg(zCt(1,2,ik))*zCt(3,2,ik) &
     +piz_dcc*conjg(zCt(2,2,ik))*zCt(3,2,ik) &
      ) &
      +2d0*occ(3,ik)*real( &
      piz_dc1*conjg(zCt(1,3,ik))*zCt(2,3,ik) &
     +piz_dc2*conjg(zCt(1,3,ik))*zCt(3,3,ik) &
     +piz_dcc*conjg(zCt(2,3,ik))*zCt(3,3,ik) &
      ) 
      

    jav_l = jav_l + (jz_intra+jz_inter)*kr(ikr)
  end do
  jav_l=jav_l*2d0/((2d0*pi)**3)*(2d0*pi*dkr*dkz) 

  call MPI_ALLREDUCE(jav_l,jav,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

  return
end subroutine current
