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

!f2py -c -m rtlib src_f90/radtran_tau_dble.f src_f90/band_dble.f90 src_f90/rosen.f src_f90/multiscatter.f90 multiscatter2_ascii.o src_f90/gas_extinction.f90 ~/multiscatter-1.2.10/lib/libmultiscatter.a

!kextKa=libSc.gasabsr98(freqKa,airTempint[i],rhowvint[i],pressint[i],ireturn)