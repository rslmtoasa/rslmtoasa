character(len=15) :: pre_processing
character(len=15) :: processing
character(len=15) :: post_processing
logical :: verbose

namelist /calculation/ verbose, pre_processing, processing, post_processing
