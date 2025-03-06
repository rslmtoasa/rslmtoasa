character(len=30) :: pre_processing
character(len=30) :: processing
character(len=30) :: post_processing
logical :: verbose

namelist /calculation/ verbose, pre_processing, processing, post_processing
