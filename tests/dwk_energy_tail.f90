program dwk_energy_tail
    use globalm
    use dimensim
    use rcomdm
    use kineticm
    use anisotropicm
    use, intrinsic :: ieee_arithmetic
    implicit none
    complex*16 al0, alnorm, alam, atau, altau, omega
    real*8 epspam, epsdet
    integer ncase, nitmax, noncon
    common/pamarg/al0, alnorm, alam, atau, altau, omega, epspam, epsdet, &
        ncase, nitmax, noncon
    character(32) arg
    integer kp, j, l, opt, pitch, grid, coll
    complex*16 z
    real*8 lamtest

    call get_command_argument(1, arg)
    read (arg, *) nepk
    nepk2 = 2*nepk - 1
    allocate (epk(nepk), zepko(nepk2), zepkn(nepk2), zepk(nepk2), fepk(nepk2))
    allocate (zresu(nepk2, 3, 0:6, 2))
    allocate (cpsi(3, 2, 2), dcdpsil(3, 2, 2), dcdpsil2(3, 2, 2))
    allocate (depsalphadpsi(3, 2, 2), depsalphadpsi2(3, 2, 2))
    allocate (omegae0(3, 2), nueff(3, 2, 2), epsalpha(3, 2, 2))
    allocate (cs(3), t(3), tm(3), dpsids(3), dpsidsm(3))
    allocate (rlm(3), slam0(3, 2), sf0(3, 2, 0:3, 4), ispecies_ek(2))
    cs = [0.1d0, 0.3d0, 0.6d0]
    t = 1.2d0
    tm = 1.2d0
    dpsids = 0.8d0
    dpsidsm = 0.8d0
    omegae0 = 0.03d0
    omega = cmplx(0.0d0, 1.0d-11, kind=8)
    rntor = -3.0d0
    b0k = 1.1d0
    omegaci0 = 50.0d0
    especies_z(1:2) = [1.0d0, -1.0d0]
    especies_m(1:2) = [2.0d0, 0.00054463d0]
    ispecies_f0 = 0
    ispecies_ek = 0
    cpsi = 1.3d0
    dcdpsil = -0.7d0
    dcdpsil2 = 0.2d0
    depsalphadpsi = -0.4d0
    depsalphadpsi2 = 0.1d0
    epsalpha = 0.9d0
    nueff = 0.05d0
    pspecies_fowp = 0.3d0
    pspecies_fowt = 0.2d0
    pspecies_ndb = 1.0d0
    ifowpsi0 = 1
    rqk = 2.1d0
    tpsi0l = 0.6d0
    tpsi0dpsil = -0.2d0
    tpsi0dlaml = 0.1d0
    slam0 = 0.0d0
    sf0 = 0.0d0
    rlm = [-1.0d0, 0.0d0, 1.0d0]
    call kenergymesh
    write (*, '(a,i5,a,es24.16)') 'GRID ', nepk, ' MAX_ENERGY ', maxval(zepkn)
    do kp = 1, 2
        call kdistribf_type0(2, 2, kp, 0.5d0)
    end do
    do j = 1, nepk2 - 1
        if (zresu(j, 1, 2, 1) == 0.0d0) then
            write (*, '(a,i5,2es24.16)') 'ZERO_DERIVATIVE ', j, zepko(j), zepkn(j)
        end if
    end do
    do coll = 0, 2
        inutype = coll
        do grid = 1, 2
            do pitch = 1, 3
                lamtest = 0.2d0*pitch
                omegab = 0.4d0 + 0.1d0*pitch
                drift = -0.3d0 + 0.2d0*pitch
                do kp = 1, 2
                    do opt = 0, 3
                        do l = 1, 3
                            call kia_trap(2, grid, kp, opt, l, lamtest, 0, z)
                            write (*, '(a,6i3,2es25.16)') &
                                'TRAP ', coll, grid, pitch, kp, opt, l, &
                                real(z), aimag(z)
                            call kia_pass(2, grid, kp, opt, l, lamtest, 0, z)
                            write (*, '(a,6i3,2es25.16)') &
                                'PASS ', coll, grid, pitch, kp, opt, l, &
                                real(z), aimag(z)
                            ispecies_ek = 1
                            call kia_trap0(2, grid, kp, opt, l, lamtest, 0, z)
                            write (*, '(a,6i3,2es25.16)') &
                                'TRAP0 ', coll, grid, pitch, kp, opt, l, &
                                real(z), aimag(z)
                            ! PASS0 is defined only at the finite birth resonance.
                            omegab = abs(real(rntor*omegae0(2, grid) - omega)/ &
                                         ((rntor*rqk + rlm(l))* &
                                          sqrt(2*epsalpha(2, grid, kp)* &
                                               especies_m(1)/especies_m(kp))))
                            call kia_pass0(2, grid, kp, opt, l, lamtest, 0, z)
                            write (*, '(a,6i3,2es25.16)') &
                                'PASS0 ', coll, grid, pitch, kp, opt, l, &
                                real(z), aimag(z)
                            ispecies_ek = 0
                            omegab = 0.4d0 + 0.1d0*pitch
                        end do
                    end do
                end do
            end do
        end do
    end do
end program
