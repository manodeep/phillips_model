! Hammer's mid square method
module msq_rand_mod
    implicit none
    integer, parameter :: i8 = selected_int_kind(12)
    integer(i8), parameter :: ten10 = 10_i8**10, ten5 = 10_i8**5
contains

integer(i8) function msq_rand(x)
    integer(i8), intent(in) :: x
    integer(i8) :: a, b, t1, t2, t3
    a = x / ten5
    b = mod(x, ten5)
    t1 = mod(a*a*ten5, ten10)
    t2 = mod(2*a*b, ten10)
    t3 = (b*b) / ten5
    msq_rand = mod(t1 + t2 + t3, ten10)
end function msq_rand

end module msq_rand_mod