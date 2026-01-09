module rutine
    private dp
    public area_n
    integer, parameter :: dp = selected_real_kind(15,307)

    contains

    subroutine area_n(xy,n,a,deg)
        !n-stevilo vozlisc
        !deg-stopnja potence y^deg v integralu
        !a - vhodna količina
        !xy - koordinate vozlisc
        integer :: n,deg
        real(dp) :: a, xy(:,:)
        if      (deg == 0) then
            do i =1,n-1
            a = a+(xy(1,i)*xy(2,i+1)-xy(1,i+1)*xy(2,i))
            end do
        else if (deg == 1) then
            do i =1,n-1
            a = a+(xy(1,i)*xy(2,i+1)-xy(1,i+1)*xy(2,i))*(xy(2,i)+xy(2,i+1))/3.0_dp
            end do
        else if (deg == 2) then
            do i =1,n-1
            a = a+(xy(1,i)*xy(2,i+1)-xy(1,i+1)*xy(2,i))*((xy(2,i)+xy(2,i+1))**2-xy(2,i)*xy(2,i+1))/6.0_dp
            end do
        else if (deg == 3) then
            do i =1,n-1
            a = a+(xy(1,i)*xy(2,i+1)-xy(1,i+1)*xy(2,i))*(xy(2,i)+xy(2,i+1))*(xy(2,i)**2+xy(2,i+1**2))/10.0_dp
            end do
        else if (deg == 4) then
            do i =1,n-1
                a = a+(xy(1,i)*xy(2,i+1) - xy(1,i+1)*xy(2,i))*((xy(2,i)**2+xy(2,i)*xy(2,i+1)+xy(2,i+1)**2)**2-xy(2,i)*xy(2,i+1)*(xy(2,i)+xy(2,i+1))**2)/15.0_dp
            end do
        else if (deg == 5) then
            do i =1,n-1
                a = a+(xy(1,i)*xy(2,i+1) - xy(1,i+1)*xy(2,i))*(xy(2,i)+xy(2,i+1))*(xy(2,i)**2+xy(2,i+1)*xy(2,i)+xy(2,i+1)**2)*(xy(2,i)**2-xy(2,i+1)*xy(2,i)+xy(2,i+1)**2)/21.0_dp
            end do
        end if
        return
    end subroutine



end module
