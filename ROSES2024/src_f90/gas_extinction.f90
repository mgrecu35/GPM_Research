subroutine get_wv_extinction(freq, nfreq, nz, airTemp,rhowv,press, kext)
implicit none
real, intent(in) :: freq(nfreq)
real, intent(in) :: airTemp(nz), rhowv(nz), press(nz)
real, intent(out) :: kext(nfreq,nz)
integer, intent(in) :: nfreq, nz
integer :: i, j, ireturn
real :: abswv, absair

do j=1,nz
    do i=1,nfreq
        call gasabsr98(freq(i),airTemp(j),rhowv(j),press(j),absair,abswv,ireturn)
        kext(i,j)=abswv+absair
    end do
end do
    
end subroutine get_wv_extinction

subroutine get_tb(nfreqg, nfreqh, nfreqc, nx, ny, nz, kext_gases,&
    kext_hyd, asym_hyd, salb_hyd, sfcBin, skTemp, lyrTemp, emis, hFreqs,gFreqs,tb, pType,umu,dz)
    implicit none
    integer:: nx,ny,nz,nfreqc,nfreqg,nfreqh,nlyr
    real :: kext_gases(nx,ny,nz,nfreqg)
    real :: kext_hyd(nx,ny,nz,nfreqh)
    real :: asym_hyd(nx,ny,nz,nfreqh)
    real :: salb_hyd(nx,ny,nz,nfreqh)
    integer :: sfcBin(nx,ny)
    integer :: pType(nx,ny)
    real :: emis(nx,ny,nfreqc)
    real, intent(out) :: tb(nx,ny,nfreqc)
    real, intent(in) :: umu
    real :: lyrtemp(nx,ny,nz+1)
    real :: lyrhgt(nz+1), fisot
    real :: skTemp(nx,ny), dz
    integer :: i,j,ic,k
    integer :: hFreqs(nfreqc), gFreqs(nfreqc)
    real :: kext1d(nz), salb1d(nz), asym1d(nz), lyrtemp1d(nz+1)
    !integer :: nlyr
    do k=1,nz+1
        lyrhgt(k)=(k-1)*dz
    end do
    fisot=2.7
    do i=1,nx
        do j=1,ny
            do ic=1,nfreqc
                nlyr=0
                do k=0,nz-sfcBin(i,j)
                    kext1d(k+1)=kext_hyd(i,j,k+sfcBin(i,j),hFreqs(ic))+kext_gases(i,j,k+sfcBin(i,j),gFreqs(ic))
                    salb1d(k+1)=salb_hyd(i,j,k+sfcBin(i,j),hFreqs(ic))*kext_hyd(i,j,k+sfcBin(i,j),hFreqs(ic))/kext1d(k+1)
                    asym1d(k+1)=asym_hyd(i,j,k+sfcBin(i,j),hFreqs(ic))
                    lyrtemp1d(k+1)=lyrtemp(i,j,k+sfcBin(i,j))
                    nlyr=nlyr+1
                    if(pType(i,j).gt.0) then
                        !print*, kext1d(k+1), salb1d(k+1), asym1d(k+1), lyrtemp1d(k+1), lyrhgt(k+1), k
                    end if
                end do
                lyrtemp1d(nlyr+1)=lyrtemp(i,j,nz+1)
                
                call radtran(umu,nlyr, tb(i,j,ic), skTemp(i,j), lyrtemp1d, lyrhgt(1:nlyr+1), &
                    kext1d(1:nlyr), salb1d(1:nlyr), asym1d(1:nlyr), fisot, emis(i,j,ic), emis(i,j,ic))
                if (pType(i,j).eq.-1) then
                    print*, 'nlyr=',nlyr, emis(i,j,ic), fisot, umu, skTemp(i,j), tb(i,j,ic)
                    !return
                end if
            end do
        end do
    end do
end subroutine get_tb
!f2py -c -m rtlib src_f90/radtran_tau_dble.f src_f90/band_dble.f90 src_f90/rosen.f src_f90/multiscatter.f90 multiscatter2_ascii.o src_f90/gas_extinction.f90 ~/multiscatter-1.2.10/lib/libmultiscatter.a

!kextKa=libSc.gasabsr98(freqKa,airTempint[i],rhowvint[i],pressint[i],ireturn)