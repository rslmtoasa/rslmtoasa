! K-path (band structure) settings
logical :: auto_kpath
integer :: nk_per_segment
integer :: override_space_group
character(len=200) :: custom_kpath_spec

namelist /kpath/ auto_kpath, nk_per_segment, override_space_group, custom_kpath_spec
