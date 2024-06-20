character(len=sl), dimension(:), allocatable :: label
character(len=sl) :: database = './'

namelist /atoms/ database, label
