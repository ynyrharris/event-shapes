# cython: language_level = 3

cimport c_event_shapes as c
from libcpp.vector cimport vector


cdef class EventShapes:
	# Hold a C++ instance of the class we're wrapping
	cdef c.EventShapes* _EventShapes

	def __cinit__(self, momenta, ndims=2):
		self._EventShapes = new c.EventShapes(momenta, ndims)

	def __dealloc__(self):
		del self._EventShapes

	def calc_all(self):
		self._EventShapes.calc_all()

	def calc_thrusts(self):
		self._EventShapes.calc_thrusts()

	def calc_sphericities(self, r):
		self._EventShapes.calc_sphericities(r)

	@property
	def thrust(self):
		return self._EventShapes.m_thrust

	@property
	def thrust_major(self):
		return self._EventShapes.m_thrust_major

	@property
	def thrust_minor(self):
		return self._EventShapes.m_thrust_minor
	
	@property
	def oblateness(self):
		return self._EventShapes.m_oblateness
	
	@property
	def broadening(self):
		return self._EventShapes.m_broadening
	
	@property
	def sph_S(self):
		return self._EventShapes.m_spher_S
	
	@property
	def sph_A(self):
		return self._EventShapes.m_spher_A
	
	@property
	def sph_C(self):
		return self._EventShapes.m_spher_C
	
	@property
	def sph_D(self):
		return self._EventShapes.m_spher_D


	def calc_lvs_t(self):
		self._EventShapes.lvs_t()

	@property
	def lvs_t(self):
		return self._EventShapes.m_lvs_t

	@property
	def lvs_tmajor(self):
		return self._EventShapes.m_lvs_tmajor

	@property
	def lvs_tminor(self):
		return self._EventShapes.m_lvs_tminor

	@property
	def lvs_obl(self):
		return self._EventShapes.m_lvs_obl

	@property
	def lvs_brd(self):
		return self._EventShapes.m_lvs_brd
