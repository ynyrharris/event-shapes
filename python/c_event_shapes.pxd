# cython: language_level = 3

from libcpp.vector cimport vector

cdef extern from "../EventShapes/EventShapes.h":
	cdef cppclass EventShapes:
		EventShapes(vector[vector[float]], unsigned int ndims) except +
		void calc_all() except +
		void calc_thrusts() except +
		void calc_sphericities(float) except +
		double m_thrust
		double m_thrust_major
		double m_thrust_minor
		double m_oblateness
		double m_broadening
		double m_spher_S
		double m_spher_A
		double m_spher_C
		double m_spher_D

		void lvs_t() except +
		double m_lvs_t
		double m_lvs_tmajor
		double m_lvs_tminor
		double m_lvs_obl
		double m_lvs_brd