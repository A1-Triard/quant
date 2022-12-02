module quant
  implicit none
contains
  subroutine app
    print *, '!!!'
  end subroutine app
end module quant

program main
  use quant, only: app
  implicit none
  call app()
end program main
