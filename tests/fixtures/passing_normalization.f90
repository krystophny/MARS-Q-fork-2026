! Exercise the unmodified production KIA_PASS and KIA_PASS0 entry points.
! Changing reference density by s changes EPSALPHA by s and every normalized
! frequency (including OMEGACI0) by sqrt(s), at fixed physical T, B and rates.
program passing_normalization
    use rcomdm
    use dimensim
    use globalm
    use kineticm
    use anisotropicm
    implicit none
    integer :: i, j, mode
    real*8 :: eps, unit_scale, center, halfwidth, transit
    complex*16 :: value
    real*8, parameter :: eps_cases(3) = [1.d0, 1.d-2, 1.d-6]

    nepk = 17
    nepk2 = 2*(nepk-1)
    allocate(epk(nepk), zepko(nepk2), zepkn(nepk2), zepk(nepk2))
    allocate(zresu(nepk2,3,5,2), rlm(1), slam0(1,2))
    allocate(nueff(3,2,2), epsalpha(3,2,2), omegae0(3,2))
    allocate(depsalphadpsi(3,2,2), cs(3), dpsids(3), t(3))
    allocate(especies_tem(3,2,2), alphaa1(3,2,2), ispecies_ek(2))
    do i = 1, nepk
        epk(i) = 0.94d0*(i-1)/(nepk-1)
    end do
    do i = 1, nepk-1
        center = 0.5d0*(epk(i+1)+epk(i))
        halfwidth = 0.5d0*(epk(i+1)-epk(i))
        zepko(2*i-1) = center-halfwidth/sqrt(3.d0)
        zepko(2*i) = center+halfwidth/sqrt(3.d0)
    end do
    zepkn = tan(0.5d0*pi*zepko)
    zresu = 0.d0
    do j = 1, 3, 2
        zresu(:,j,1,1) = -0.4d0*exp(-zepkn)
        zresu(:,j,2,1) = -exp(-zepkn)
    end do
    rlm = 1.d0
    slam0 = -1.d0
    cs = [0.d0, 0.5d0, 1.d0]
    dpsids = 1.d0
    t = 1.d0
    depsalphadpsi = 0.d0
    especies_tem = 1.d0
    alphaa1 = 1.d0
    especies_z(1:2) = [1.d0, -1.d0]
    especies_m(1:2) = [2.d0, 0.00054463d0]
    ispecies_f0(1:2) = 0
    pspecies_fowp(1:2) = 0.d0
    ifowpsi0 = 1
    rntor = -3.d0
    b0k = 1.d0
    rqk = 0.7d0
    omegab = 0.8d0
    call set_test_omega()

    do mode = 0, 2
        inutype = mode
        ispecies_ek = 0
        do i = 1, size(eps_cases)
            eps = eps_cases(i)
            unit_scale = sqrt(eps)
            epsalpha = eps
            omegaci0 = unit_scale
            nueff = 0.2d0*unit_scale
            omegae0 = 0.37d0*unit_scale
            call kia_pass(2,1,1,0,1,0.1d0,0,value)
            write(*,'(I2,1X,I2,1X,3(ES25.16E3,1X))') mode,i,eps,dble(value),aimag(value)
        end do
    end do

    ! KIA_PASS0 is a finite-birth singular-limit coefficient. Its H4 local
    ! is unused; exercise the nonzero resonance, not its Maxwellian return.
    ispecies_ek = 1
    do i = 1, size(eps_cases)
        eps = eps_cases(i)
        epsalpha = eps
        omegaci0 = sqrt(eps)
        nueff = 0.2d0*sqrt(eps)
        transit = omegab*(rntor*rqk+rlm(1))*sqrt(2.d0*eps)
        omegae0 = -transit/rntor
        call kia_pass0(2,1,1,0,1,0.1d0,0,value)
        write(*,'(I2,1X,I2,1X,3(ES25.16E3,1X))') 3,i,eps,dble(value),aimag(value)
    end do
end program passing_normalization
