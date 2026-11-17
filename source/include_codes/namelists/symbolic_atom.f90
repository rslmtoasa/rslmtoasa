character(len=sl), dimension(:), allocatable :: label
character(len=pl) :: database = './'

namelist /atoms/ database, label
