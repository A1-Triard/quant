module quant
  use fourier
  implicit none
  integer, parameter :: f = kind(1.0e0)
  integer, parameter :: max_time = 10
  integer, parameter, dimension(3) :: univ = [128, 128, 128]
  real(f), parameter :: k = 1.0_f
  real(f), parameter :: s_r = 8.0_f
  real(f), parameter :: s_u = 10.0_f
  real(f), parameter, dimension(3) :: s_pos = [64.5_f, 64.5_f, 64.5_f]
  real(f), parameter, dimension(3) :: r0 = [10.0_f, 64.5_f, 64.5_f]
  real(f), parameter :: d0 = 8.0_f
  real(f), parameter :: k0 = 1.0_f
contains
  subroutine stop_error(msg)
    use iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: msg
    write(error_unit, '(a)') msg
    stop 1
  end subroutine

  subroutine assert(cond, msg)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: msg
    if (.not. cond) then
      call stop_error(msg)
    end if
  end subroutine assert

  function dist(a, b)
    real(f), dimension(:), intent(in) :: a, b
    real(f) :: dist
    integer :: i
    call assert(ubound(a, 1) == ubound(b, 1), 'dist: size(a) /= size(b)')
    dist = 0.0_f
    do i = 1, ubound(a, 1)
      dist = dist + (a(i) - b(i)) ** 2
    end do
    dist = sqrt(dist)
  end function dist

  function pot(r)
    integer, dimension(3), intent(in) :: r
    real(f) :: pot
    if (dist(real(r, f), s_pos) <= s_r) then
      pot = s_u
    else
      pot = 0.0_f
    end if
  end function pot

  subroutine calc_kin_u(kin_u)
    complex(f), dimension(univ(1), univ(2), univ(3)), intent(out) :: kin_u
    integer :: x, y, z
    real(f) :: h
    do x = 1, univ(1)
      do y = 1, univ(2)
        do z = 1, univ(3)
          h = (k / 4.0_f) * (real(x, f) ** 2 + real(y, f) ** 2 + real(z, f) ** 2)
          kin_u(x, y, z) = complex(cos(h), -sin(h))
        end do
      end do
    end do
  end subroutine calc_kin_u

  subroutine calc_pot_u(pot_u)
    complex(f), dimension(univ(1), univ(2), univ(3)), intent(out) :: pot_u
    integer :: x, y, z
    real(f) :: h
    do x = 1, univ(1)
      do y = 1, univ(2)
        do z = 1, univ(3)
          h = pot([x, y, z])
          pot_u(x, y, z) = complex(cos(h), -sin(h))
        end do
      end do
    end do
  end subroutine calc_pot_u

  subroutine evolve(phi, kin_u, pot_u)
    complex(f), dimension(univ(1), univ(2), univ(3)), intent(inout) :: phi
    complex(f), dimension(univ(1), univ(2), univ(3)), intent(in) :: kin_u, pot_u
    call fft3_inplace(phi)
    phi = kin_u * phi
    call ifft3_inplace(phi)
    phi = pot_u * phi
    call fft3_inplace(phi)
    phi = kin_u * phi
    call ifft3_inplace(phi)
  end subroutine evolve

  subroutine calc_probab(phi, probab)
    complex(f), dimension(univ(1), univ(2), univ(3)), intent(in) :: phi
    real(f), dimension(univ(1), univ(2), univ(3)), intent(out) :: probab
    probab = abs(phi) ** 2
  end subroutine calc_probab

  subroutine draw_frame(time, probab, pixel, ncol)
    integer, intent(in) :: time
    real(f), dimension(univ(1), univ(2), univ(3)), intent(in) :: probab
    integer, dimension(:, :, :), intent(out) :: pixel
    integer, intent(in) :: ncol
    integer :: x, y, z
    real(f) :: p
    call assert(time >= 1 .and. time <= ubound(pixel, 1), 'time < 1 .or. time > ubound(pixel, 1)')
    call assert(univ(1) == ubound(pixel, 2), 'draw_frame: ubound(pixel, 2) /= univ(1)')
    call assert(univ(2) == ubound(pixel, 3), 'draw_frame: ubound(pixel, 2) /= univ(1)')
    do x = 1, univ(1)
      do y = 1, univ(2)
        p = 0.0_f
        do z = 1, univ(3)
          p = p + probab(x, y, z)
        end do
        pixel(time, x, y) = nint(max(min(p, 1.0_f), 0.0_f) * ncol)
      end do
    end do
  end subroutine

  subroutine calc_colormap(colormap)
    integer, dimension(:, 0:), intent(out) :: colormap
    integer :: c
    call assert(ubound(colormap, 1) == 3, 'ubound(colormap, 1) /= 3')
    colormap(1, :) = 0
    colormap(2, :) = 0
    do c = 0, ubound(colormap, dim=2)
      colormap(3, c) = nint((255.0e0 * c) / ubound(colormap, dim=2))
    end do
  end subroutine calc_colormap

  subroutine init_phi(phi)
    complex(f), dimension(univ(1), univ(2), univ(3)), intent(out) :: phi
    real(f) :: amp, arg, amp_max
    integer :: x, y, z
    amp_max = 0.0_f
    do x = 1, univ(1)
      do y = 1, univ(2)
        do z = 1, univ(3)
          amp = exp(-(dist([real(x, f), real(y, f), real(z, f)], r0) / d0) ** 2)
          arg = -k0 * (real(x, f) - r0(1))
          phi(x, y, z) = amp * complex(cos(arg), sin(arg))
          amp_max = max(amp_max, amp)
        end do
      end do
    end do
    phi = phi / amp_max
  end subroutine init_phi

  subroutine app
    use gif_module
    integer, dimension(3, 0:64) :: colormap
    complex(f), allocatable, dimension(:, :, :) :: kin_u, pot_u, phi
    real(f), allocatable, dimension(:, :, :) :: probab
    integer, allocatable, dimension(:, :, :) :: pixel
    integer :: time
    allocate(pixel(max_time, univ(1), univ(2)))
    allocate(kin_u(univ(1), univ(2), univ(3)))
    allocate(pot_u(univ(1), univ(2), univ(3)))
    allocate(probab(univ(1), univ(2), univ(3)))
    allocate(phi(univ(1), univ(2), univ(3)))
    call calc_kin_u(kin_u)
    call calc_pot_u(pot_u)
    call init_phi(phi)
    call calc_probab(phi, probab)
    call draw_frame(1, probab, pixel, ubound(colormap, dim=2))
    do time = 2, max_time
      print *, time
      call evolve(phi, kin_u, pot_u)
      call calc_probab(phi, probab)
      call draw_frame(time, probab, pixel, ubound(colormap, dim=2))
    end do
    call calc_colormap(colormap)
    call write_animated_gif('quant.gif', pixel, colormap, delay=100)
    deallocate(pixel)
    deallocate(kin_u)
    deallocate(pot_u)
    deallocate(probab)
    deallocate(phi)
  end subroutine app
end module quant

program main
  use quant, only: app
  implicit none
  call app()
end program main
