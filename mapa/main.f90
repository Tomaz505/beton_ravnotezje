program main
    use rutine


    integer, parameter :: dp = selected_real_kind(15,307)
    real(dp), parameter :: pi = 3.141592654
    integer :: n_c, n_a, n_p, u, analiza  !st vozlisc, stevilo armaturnih palic, stevilo kablov

    real(dp) :: f_c, e_c, phi_cr,x_area
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



    !RAČUN KARAKTERISTIK IN PREMIK TEZISCA V IZODISCE
    block
        real(dp) :: sx=0.0_dp

        call area_n(xy_c,n_c,a,0)
        call area_n(xy_c,n_c,sx,1)

        xy_c(2,:) = xy_c(2,:)-sx/a
        xy_s(2,:) = xy_s(2,:)-sx/a

        if (n_p /= 0) then
            xy_p(2,:) = xy_p(2,:) - sx/a
        end if

        call area_n(xy_c,n_c,ix,2)

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
        real(dp) ::  def_plc(3), z_extr(2)
        def_plc(:) = 0
        z_extr(1) = minval(xy_c(2,:))
        z_extr(2) = maxval(xy_c(2,:))

        def_plc = (def_pl-eps_sh)/(1+phi_cr)

        print*," "
        print*,"    Račun končan."
        print*," "

        ! dopolni če je eps_sh(3) /= 0
        call write_js(xy_c,n_c,xy_s,n_s,r_s,def_plc, def_pl, z_extr)

    end block


    !NERAZPOKAN DEL PREREZA MED ymin IN ymax
    block
        real(dp) :: xy_eff(2,2*n_c), z_crack(2), def_pl_crac(3), z_extr(2), z_range(2), jac_eq(2,2), f_eq(2)

        real(dp) :: i0,i1,i2,i3


        z_range(:) = 0
        z_extr(1) = minval(xy_c(2,:))
        z_extr(2) = maxval(xy_c(2,:))

        de_pl_crac = def_pl
        if (eps_sh(3) == 0) then
            !racun ene nicle deformacij

            !
            !do while
            !   i0 = 0
            !   i1 = 0
            !   i2 = 0
            !   i3 = 0

            !   !dolocitev z koordinate razpoke preveri predznak kappa
            !   z_crac(:) == -(def_pl_crac(1)-eps_sh(1))/(def_pl_crac(2)-eps_sh(2))
            !   if ((z_crac(1) > z_extr(1)) .and. (z_crac(1) < z_extr(2)) then
            !       if (def_pl_crac(2) >0) then
            !           z_range = (/z_crac(1),z_extr(2)+1/)
            !       else if (def_pl_crac(2) < 0)
            !            z_range = (/z_extr(1)-1,z_crac(1)/)
            !    else
            !        zrange = (/z_extr(1)-1,z_extr(2)+1/)
            !   end if
            !   call eff_section(xy_c,n_c,z_range(1),z_range(2),xy_eff)
            !
            !   !Racun momentov (integral y po mnogokotniku) rabim samo i0, i1 in i2
            !   i0 = area_n(xy_eff,2*n_c,i0,0)
            !   i1 = area_n(xy_eff,2*n_c,i1,0)
            !   i2 = area_n(xy_eff,2*n_c,i2,0)
            !  !i3 = area_n(xy_eff,2*n_c,i3,0)
            !

            !
            !end do
            !



        else
            !test za stevilo nicel b**2-4*a*c >,<,= 0
            !racun ene oz. dveh nicel
            z_crac(:) = -(def_pl_crac(1)-eps_sh(1))/(def_pl_crac(2)-eps_sh(2))



        end if




        call eff_section(xy_c,n_c,-20.0_dp,12.0_dp,xy_eff)

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
end program
