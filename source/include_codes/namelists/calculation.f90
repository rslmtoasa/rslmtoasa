character(len=30) :: pre_processing
character(len=30) :: processing
character(len=30) :: post_processing
character(len=30) :: gf_route
logical :: verbose

namelist /calculation/ verbose, pre_processing, processing, post_processing, gf_route
