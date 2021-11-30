# cython: language_level = 3

from libcpp.vector cimport vector

cdef extern from "../EventShapes/EventShapes.h":
	cdef cppclass EventShapes:
		EventShapes(vector[vector[float]]) except +
		void calc_all() except +
		double m_thrust
		double m_thrust_major
		double m_thrust_minor
		double m_oblateness
		double m_broadening
		double m_lin_spher_S
		double m_lin_spher_A
		double m_lin_spher_C
		double m_lin_spher_D

		void lvs_t() except +
		double m_lvs_t
		double m_lvs_tmajor
		double m_lvs_tminor
		double m_lvs_obl
		double m_lvs_brd