!!=====================================================================
!!
!! Ensight outputing
!!
!!=====================================================================
module ensight_var

implicit none

logical :: ENSIGHT_BIN  !- Binary format
integer :: ENSIGHT_DIM  !- Dimension of data
integer :: SLCNUM       !- Slice position

character(len=10)               :: CASEFLU !- Case name
real(kind=8), dimension(:), allocatable :: TAB_TIME !- Time of outputing

integer :: NBNODE

character(len=10), dimension(:), allocatable :: NODELIST
integer,           dimension(:), allocatable :: NODETYPE


end module ensight_var
