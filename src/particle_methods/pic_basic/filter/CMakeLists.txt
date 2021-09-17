ADD_LIBRARY(sll_filter STATIC
  sll_m_filter_base_1d.F90
  sll_m_filter_base_3d.F90
  sll_m_binomial_filter.F90
  sll_m_fft_filter_1d.F90
  sll_m_fft_filter_3d.F90
  )

TARGET_LINK_LIBRARIES( sll_filter sll_fft
	)


