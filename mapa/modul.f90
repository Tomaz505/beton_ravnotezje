module rutine
    private dp,pi
    public area_n

    integer, parameter :: dp = selected_real_kind(15,307)
    real(dp), parameter :: pi = 3.141592654

    contains

    !Izracuna deg-ti statični moment (Povrsina deg=0, Vztrajnostni mom deg=2)
    !Visje momente (3 oz 5) potrebujes pri racunu napetosti (deformacija v betonu stopnje 1 oz 2)
    !Integriras po mnogokotniku y^deg
    subroutine area_n(xy,n,a,deg)
        !n-stevilo vozlisc
        !deg-stopnja potence y^deg v integralu
        !a - vhodna količina
        !xy - koordinate vozlisc
        !ze uposteva simetrijo in ni treba mnozit z 2
        integer :: n,deg
        real(dp) :: a, xy(:,:)

        do i=1,n-1
            x = 0.0
            do j= 0,deg
                x = x + xy(2,i)**(deg-j)*xy(2,i+1)**(j)
            end do
            a = a + x*(xy(1,i)*xy(2,i+1)-xy(1,i+1)*xy(2,i))*2/(1+deg)/(2+deg)
        end do

        return
    end subroutine

    !Zapise potrebne kolicine v datoteko out.js
    subroutine write_js(xy,n,xys,ns,ds,def_c,def_s,h,an)
        !n-stevilo vozlisc
        !deg-stopnja potence y^deg v integralu
        !a - vhodna količina
        !xy - koordinate vozlisc
        integer :: n,ns,an
        real(dp) :: xy(:,:),def_c(3),def_s(2),h(2),xys(:,:),ds(:),z,e

        open (unit = 1,file = "mapa/out.js",status = "old")
        write(1,*) "const sec_coor = `"
        do i=1,n
            write(1,*) xy(:,i)
        end do
        do i=n-1,1,-1
            write(1,*) (/ -xy(1,i),xy(2,i)/)
        end do
        write(1,*) "`;"

        write(1,*) "const a_coor = `"
        do i=1,ns
            write(1,*) xys(:,i) , ds(i)
        end do
        do i=ns,1,-1
            write(1,*) (/ -xys(1,i),xys(2,i),  ds(i)/)
        end do
        write(1,*) "`;"

        write(1,*) "const eps_c = `"
        do i=1,101
            z = (h(1)*(100-(i-1))+ h(2)*(i-1))/100.0_dp
            e = def_c(1)-def_c(2)*z+def_c(3)*z**2
            if (an == 2) then
                e = min(0.0,e)
            end if
            write(1,*) e , z
        end do
        write(1,*) "`;"

        write(1,*) "const eps_s = `"
        do i=1,ns !101
            !z = (h(1)*(100-(i-1))+ h(2)*(i-1))/100.0_dp
            write(1,*) def_s(1)-def_s(2)*xys(2,i),xys(2,i)!*z, z
        end do
        write(1,*) "`;"


        close(1)
    end subroutine

    subroutine crac_linmat_linsh(xy_c,n_c,xy_s,n_s,def_pl_crac,z_extr,eps_sh,e_s,e_c,f_eq,r_s,phi_cr)
        integer :: n_c,n_s
        real(dp) :: xy_c(:,:), xy_s(:,:),def_pl_crac(3),z_extr(2),eps_sh(3),e_s,e_c,f_eq(2),r_s(:),phi_cr

        real(dp) :: i0,i1,i2, z_crac, xy_eff(2,2*n_c),z_range(2)
        z_range(:) = 0

        i0 = 0.0_dp
        i1 = 0.0_dp
        i2 = 0.0_dp

        !dolocitev z koordinate razpoke preveri predznak kappa
        z_crac = -(def_pl_crac(1)-eps_sh(1))/(def_pl_crac(2)-eps_sh(2))

        if ((z_crac > z_extr(1)) .and. (z_crac < z_extr(2))) then
            if (def_pl_crac(2)-eps_sh(2) < 0) then
                z_range = (/z_crac,z_extr(2)+1.0_dp/)
            else if (def_pl_crac(2)-eps_sh(2) > 0) then
                z_range = (/z_extr(1)-1.0_dp,z_crac/)
            end if
        else if(((z_crac < z_extr(1)) .and.  (def_pl_crac(2)-eps_sh(2) < 0)) .or. (((z_crac > z_extr(2)) .and.  (def_pl_crac(2)-eps_sh(2) > 0)))) then
                z_range = (/z_extr(1)-1.0_dp,z_extr(2)+1.0_dp/)
        else
            z_range = (/0.0_dp,0.0_dp/)

        end if


        call eff_section(xy_c,n_c,z_range(1),z_range(2),xy_eff)

        !Racun momentov (integral y po mnogokotniku) rabim samo i0, i1 in i2
        call area_n(xy_eff,2*n_c,i0,0)
        call area_n(xy_eff,2*n_c,i1,1)
        call area_n(xy_eff,2*n_c,i2,2)


        do i=1,n_s
            f_eq =f_eq + (/1.0_dp,-xy_s(2,i)/)*2*e_s*pi/4*r_s(i)**2*(def_pl_crac(1)-def_pl_crac(2)*xy_s(2,i))
        end do

        f_eq = f_eq + e_c/(1+phi_cr)*(/(i1*(def_pl_crac(2)-eps_sh(2))+i0*(def_pl_crac(1)-eps_sh(1))) ,(i2*(def_pl_crac(2)-eps_sh(2))+i1*(def_pl_crac(1)-eps_sh(1)))/)

    end subroutine





    !Doloci del prereza, ki je med y koordinatama ymin in ymax
    !xy_sub je mnogokotnik po katerem lahko integriramo napetosti
    subroutine eff_section(xy,n,ymin,ymax,xy_sub)
        logical :: state
        integer :: n,i1
        real(dp) :: xy(:,:),ymin,ymax,xy_sub(:,:),x


        ! ali je je prvo vozlisce v prerezu
        state = ((xy(2,1) >= ymin) .and. (xy(2,1) <= ymax))
        i1 = 2
        xy_sub(:,:) = 0.0_dp
        !print*,"state = ",state

        if (state) then
            xy_sub(:,1) = xy(:,1)
        else if (xy(2,1)>ymax) then
            xy_sub(:,1) = (/0.0_dp,ymax/)
        else if (xy(2,1)<ymin) then
            xy_sub(:,1) = (/0.0_dp,ymin/)
        end if


        do i=1,n-1

            if (state) then
                !sem v prerezu. preverim stanje naslednjega vozlišča
                if ((xy(2,i+1)>=ymin) .and. (xy(2,i+1)<=ymax)) then
                   xy_sub(:,i1) = xy(:,i+1)
                    i1 = i1+1
                else if (xy(2,i+1)>ymax) then
                    x = (ymax*(xy(1,i+1)-xy(1,i))+xy(2,i+1)*xy(1,i)-xy(2,i)*xy(1,i+1))/(xy(2,i+1)-xy(2,i))
                    xy_sub(:,i1) = (/x,ymax/)
                    i1 = i1+1
                    state = (.not. state)
                else if (xy(2,1)<ymin) then
                    x = (ymin*(xy(1,i+1)-xy(1,i))+xy(2,i+1)*xy(1,i)-xy(2,i)*xy(1,i+1))/(xy(2,i+1)-xy(2,i))
                    xy_sub(:,i1) = (/x,ymin/)
                    i1 = i1+1
                    state = (.not. state)
                end if
            else
                if ((xy(2,i)<ymin) .and. (xy(2,i+1)>ymin)) then
                    x = (ymin*(xy(1,i+1)-xy(1,i))+xy(2,i+1)*xy(1,i)-xy(2,i)*xy(1,i+1))/(xy(2,i+1)-xy(2,i))

                    xy_sub(:,i1) = (/x, ymin/)

                    if (xy(2,i+1)> ymax) then
                        x = (ymax*(xy(1,i+1)-xy(1,i))+xy(2,i+1)*xy(1,i)-xy(2,i)*xy(1,i+1))/(xy(2,i+1)-xy(2,i))
                        xy_sub(:,i1+1) = (/x, ymax/)
                    else
                        xy_sub(:,i1+1) = xy(:,i+1)
                        state = (.not. state)
                    end if
                    i1 = i1+2

                else if ((xy(2,i)>ymax) .and. (xy(2,i+1)<ymax)) then
                    x = (ymax*(xy(1,i+1)-xy(1,i))+xy(2,i+1)*xy(1,i)-xy(2,i)*xy(1,i+1))/(xy(2,i+1)-xy(2,i))

                    xy_sub(:,i1) = (/x, ymax/)

                    if (xy(2,i+1)< ymin) then
                        x = (ymin*(xy(1,i+1)-xy(1,i))+xy(2,i+1)*xy(1,i)-xy(2,i)*xy(1,i+1))/(xy(2,i+1)-xy(2,i))
                        xy_sub(:,i1+1) = (/x, ymin/)
                    else
                        xy_sub(:,i1+1) = xy(:,i+1)
                        state = (.not. state)
                    end if
                    i1 = i1+2
                end if
            end if
            if (i==n-1) then
                xy_sub(:,i1) = (/0.0_dp,xy_sub(2,i1-1)/)
            end if
        end do
    end subroutine



end module
