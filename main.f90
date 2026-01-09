program main
    integer, parameter :: dp = selected_real_kind(15,307)
    real(dp), parameter :: pi = 3.141592654
    integer :: n_c, n_a, n_p, u, analiza  !st vozlisc, stevilo armaturnih palic, stevilo kablov

    real(dp) :: f_c, e_c, phi_cr
    real(dp) :: f_s, e_s
    real(dp) :: f_p, e_p
    real(dp) ::  a=0.0_dp, ix=0.0_dp
    real(dp) :: med, ned

    real(dp), allocatable :: xy_c(:,:), xy_s(:,:), xy_p(:,:), r_s(:), r_p(:), p_p(:)
    real(dp) :: def_pl(3), eps_sh(3)

    open(newunit = u, file = "in.txt",status = "old")
    read(u,*) analiza, med,ned,n_c,n_s,n_p
    close(u)

    if (np == 0) then
        allocate(xy_c(2,n_c), xy_s(2,n_s), r_s(n_s))
        open(newunit = u, file = "in.txt",status = "old")
            read(u,*) analiza,med,ned,n_c,n_s,n_p,f_c,e_c,phi_cr,f_s,e_s,xy_c,eps_sh,xy_s,r_s
        close(u)
    else
        allocate(xy_c(2,n_c), xy_s(2,n_s),xy_p(2,n_p), r_s(n_s), r_p(n_p),p_p(n_p))
        open(newunit = u, file = "in.txt",status = "old")
            read(u,*) analiza, med,ned,n_c,n_s,n_p,f_c,e_c,phi_cr,f_s,e_s,xy_c,eps_sh,xy_s,r_s,f_p,e_p,xy_p,r_p,p_p
        close(u)
    end if
    !Skladnost koordinatnega sistema
    eps_sh(2) = -eps_sh(2)


    !RAČUN KARAKTERISTIK
    block
        real(dp) :: sx=0.0_dp
        do i = 1,n_c-1
            a = a + (xy_c(1,i)*xy_C(2,i+1) - xy_c(1,i+1)*xy_c(2,i))
            sx = sx + (xy_c(1,i+1)-xy_c(1,i))*(xy_c(2,i)**2 + xy_c(2,i+1)**2+ xy_c(2,i)*xy_c(2,i+1))
        end do
        sx = sx/3
        xy_c(2,:) = xy_c(2,:)+sx/a
        xy_s(2,:) = xy_s(2,:)+sx/a
        if (n_p == 0) then
        else
            xy_p(2,:) = xy_p(2,:) - sx/a
        end if

        do i = 1,n_c-1
            ix = ix+ (xy_c(1,i)*xy_c(2,i+1) - xy_c(1,i+1)*xy_c(2,i))*(xy_c(2,i)**2 + xy_c(2,i+1)**2+ xy_c(2,i)*xy_c(2,i+1))/6
        end do

    end block



    !RAČUN NERAZPOKANEGA PREREZA
    block
        real(dp) :: c_mat(2,2), f_vec(2)
        f_vec = (/ned, med/)
        c_mat(:,:) = 0

        !print*, eps_sh,a,ix,e_c

        do i =1,n_s
            c_mat(1,2) = c_mat(1,2) - 2*e_s*xy_s(2,i)*pi/4 *r_s(i)**2
            c_mat(1,1) = c_mat(1,1) + 2*e_s*pi/4 *r_s(i)**2
            c_mat(2,2) = c_mat(2,2) + 2*e_s* xy_s(2,i)**2 *pi/4 *r_s(i)**2
        end do

        c_mat(1,1) = c_mat(1,1) +  a*e_c/(1+phi_cr)
        c_mat(2,2) = c_mat(2,2) + ix*e_c/(1+phi_cr)
        c_mat(2,1) = c_mat(1,2)

        f_vec(1) = f_vec(1) +  (eps_sh(1)*a+eps_sh(3)*ix)*e_c/(1+phi_cr)
        f_vec(2) = f_vec(2) + eps_sh(3)*ix*e_c/(1+phi_cr)


        !print*, c_mat
        !print*, f_vec

        def_pl(1:2) = (/c_mat(2,2)*f_vec(1)-c_mat(1,2)*f_vec(2) , -c_mat(1,2)*f_vec(1)+c_mat(1,1)*f_vec(2)/) /(c_mat(1,1)*c_mat(2,2)-c_mat(1,2)**2)
        def_pl(3) = 0
    end block

    !PRINTANJE REZULTATOV
    block
        real(dp) ::  def_plc(3), z_extr(3)
        def_plc(:) = 0
        z_extr(1) = minval(xy_c(2,:))
        z_extr(2) = maxval(xy_c(2,:))

        def_plc = (def_pl-eps_sh)/(1+phi_cr)

        if (eps_sh(3) /= 0) then
            z_extr(3) = (-eps_sh(2)/2/eps_sh(3))
        else
            z_extr(3) = 0
        end if
        print *, " "
        print '(a)', "NERAZPOKAN PREREZ"

        print '(a)', "      1. Deformacijska ravnina"
        print * , def_pl(1:2)
        print *, " "

        print '(a)', "      2. Napetosti v armaturnih palicah"
        print * , (def_pl(1) - def_pl(2) * xy_s(2,:))*e_s
        print *, " "

        print '(a)', "      3. Napetostne deformacije"
        print *, def_plc
        print *, " "


        if ((eps_sh(3) /= 0) .and. ((z_extr(3) > z_extr(1)) .and. (z_extr(3) < z_extr(2)) )) then
            print '(a)', "      4. Napetosti v betonu     (spodaj, zgoraj, teme parabole)"
            print *, z_extr
            print * , (def_plc(1) - def_plc(2) * z_extr)*e_c
        else
            print '(a)', "      4. Napetosti v betonu     (spodaj, zgoraj)"
            print * , (def_plc(1) - def_plc(2) * z_extr(1:2))*e_c
        end if

        print *, " "
    end block


    !RAZPOKAN PREREZ
    block
        real(dp) :: a_crac = 0, i_crac = 0, z_crac(2), def_pl_crac(3)

        def_pl_crac = def_pl

        do i=1,3
            ! y koordinata razpok
            if (eps_sh(3) == 0) then
                z_crac(:) = -(def_pl_crac(1)-eps_sh(1))/(def_pl_crac(2)-eps_sh(2))
            else
                z_crac(1) = -(def_pl_crac(2)-eps_sh(2)) + sqrt((def_pl_crac(2)-eps_sh(2))**2 + 4*eps_sh(3)*(def_pl_crac(1)-eps_sh(1)))
                z_crac(2) = -(def_pl_crac(2)-eps_sh(2)) - sqrt((def_pl_crac(2)-eps_sh(2))**2 + 4*eps_sh(3)*(def_pl_crac(1)-eps_sh(1)))
                z = -0.5*z/eps_sh(3)
            end if

        end do
    end block

    print'(a)', "Pritisni Enter"
    read *,
end program
