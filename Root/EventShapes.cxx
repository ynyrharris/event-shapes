
#include "EventShapes/EventShapes.h"

#include <Eigen/Dense>

#include <iostream>
#include <math.h>

using namespace Eigen;


EventShapes::EventShapes()
	: m_ndims(2),
	  m_thrust_axis(nullptr),
	  m_thrust_major_axis(nullptr),
	  m_thrust_minor_axis(nullptr),
	  m_thrust(-1.),
	  m_thrust_major(-1.),
	  m_thrust_minor(-1.),
	  m_oblateness(-1.),
	  m_broadening(-1.),
	  m_spher_S(-1),
	  m_spher_A(-1),
	  m_spher_C(-1.),
	  m_spher_D(-1.),
	  m_lvs_t(-1.),
	  m_lvs_tmajor(-1.),
	  m_lvs_tminor(-1.),
	  m_lvs_obl(-1.),
	  m_lvs_brd(-1.) {}


EventShapes::EventShapes(const std::vector<std::vector<float>>& momenta, unsigned int m_ndims)
	: m_ndims(m_ndims),
	  m_thrust_axis(nullptr),
	  m_thrust_major_axis(nullptr),
	  m_thrust_minor_axis(nullptr),
	  m_thrust(-1.),
	  m_thrust_major(-1.),
	  m_thrust_minor(-1.),
	  m_oblateness(-1.),
	  m_broadening(-1.),
	  m_spher_S(-1),
	  m_spher_A(-1),
	  m_spher_C(-1.),
	  m_spher_D(-1.),
	  m_lvs_t(-1.),
	  m_lvs_tmajor(-1.),
	  m_lvs_tminor(-1.),
	  m_lvs_obl(-1.),
	  m_lvs_brd(-1.) {

	m_momenta.reserve(momenta.size());

	// If problem is presented in 2d
	if (m_ndims == 2) {
		for (auto& x : momenta) {
			m_momenta.push_back(Vector3f(x.at(0), x.at(1), 0.));
		}

	// Else if problem is presented in 3d
	} else {
		for (auto& x : momenta) {
			m_momenta.push_back(Vector3f(x.at(0), x.at(1), x.at(2)));
		}
	}

	m_ntracks = m_momenta.size();
	m_min_ntracks = m_ndims;
}


double EventShapes::calc_t(const std::vector<Vector3f>& X, const Vector3f& axis) {
	/*
	 Function that calculates the thrust value of a set of input vectors 
	 'X' about the thrust axis 'axis'.
	*/

	double p = 0.;
	double q = 0.;

	for (const Eigen::Vector3f& x : X) {
		p += fabs(axis.dot(x));
		q += x.norm();
	}

	double t = q > 0. ? p / q : -1.;
	return t;
}


const std::pair<Vector3f, double> EventShapes::tradThrust(const std::vector<Vector3f>& X) {
	/*
	 Calculate the thrust axis using the traditional multiple-trials
	 approach.  This implementation starts from N^{(d - 1)/2} random
	 initial axes and builds a candidate thrust axis from each.
	 The normalised resultant that has the greatest value of thrust is
	 taken as the thrust axis.
	*/

	// Get thrust axis
	Vector3f axis;
	double axis_t = 0.;

	// Start from multiple random initial axes
	float N_trials = pow(X.size(), m_ndims - 0.5);
	for (unsigned int i = 0; i < N_trials; i++) {
		Vector3f init = Vector3f::Random();

		Vector3f v = Vector3f::Zero();
		for (const auto& x : X) {
			init.dot(x) > 0 ? v += x : v -= x;
		}
		v.normalize();

		// Keep the axis if it increases the thrust value
		double v_t = calc_t(X, v);
		if (v_t > axis_t) {
			axis = v;
			axis_t = v_t;
		}
	}

	return std::pair<Vector3f, double> (axis, axis_t);
}


void EventShapes::calcThrust() {
	/* Function that calculates the thrust proper axis and values.
	 * The class member data are set with the calculated values.
	 */

	std::pair<Vector3f, double> pair = tradThrust(m_momenta);

	// Set class data members
	if (pair.second < 0.) {
		m_thrust_axis = nullptr;
		m_thrust = -1.;
	} else {
		m_thrust_axis = new Vector3f(pair.first);
		m_thrust = pair.second;
	}

}


void EventShapes::calcThrustMajor() {
	/* Function that calculates the thrust major axis and value.
	 * This is the axis in the plane perpendicular to the thrust axis 
	 * along which the energy flow is greatest in that plane.
	 */

	// Function depends on the thrust axis
	if (!m_thrust_axis) calcThrust();
	if (!m_thrust_axis) return;

	const Vector3f z = *m_thrust_axis;

	// Calculate thrust in plane perpendicular to thrust axis proper
	std::vector<Vector3f> G;
	G.reserve(m_momenta.size());

	for (const auto& x : m_momenta) {
		G.push_back(x - x.dot(z) * z);
	}

	std::pair<Vector3f, double> pair = tradThrust(G);

	// Set class data members
	if (pair.second < 0.) {
		m_thrust_major_axis = nullptr;
		m_thrust_major = -1.;
	} else {
		m_thrust_major_axis = new Vector3f(pair.first);
		m_thrust_major = pair.second;
	}

	return;
}


void EventShapes::calcThrustMinor() {
	/* Function that calculates the thrust minor axis and value.
	 * This axis is defined as being orthogonal to the thrust and 
	 * thrust major axes.
	 */

	// This function depends on the thrust and thrust major axes
	if (!m_thrust_major_axis) calcThrustMajor();
	if (!m_thrust_major_axis) return;

	// Alias the existing axes
	const Vector3f thrust_axis = *m_thrust_axis;
	const Vector3f thrust_major_axis = *m_thrust_major_axis;

	Vector3f thrust_minor_axis = thrust_axis.cross(thrust_major_axis);
	double thrust_minor = calc_t(m_momenta, thrust_minor_axis);

	// Sanity check
	float a_dot_b = fabs(thrust_axis.dot(thrust_major_axis));
	float b_dot_c = fabs(thrust_axis.dot(thrust_minor_axis));
	float epsilon = 1e2 * std::numeric_limits<float>::epsilon();

	if (a_dot_b > epsilon || b_dot_c > epsilon) {
		std::cout << "Thrust axes not orthogonal!" << std::endl;
		std::cout << "a_dot_b: " << a_dot_b << std::endl;
		std::cout << "b_dot_c: " << b_dot_c << std::endl;
		std::cout << "epsilon: " << epsilon << std::endl;
	}

	// Set class data members
	m_thrust_minor_axis = new Vector3f(thrust_minor_axis);
	m_thrust_minor = thrust_minor;
}


void EventShapes::calcOblateness() {
	/* Function that calculates the oblateness about the thrust axis.
	 */

	// This function depends on the major and minor thrust calculations
	if (!m_thrust_minor_axis) calcThrustMinor();
	if (!m_thrust_minor_axis) return;

	m_oblateness = m_thrust_major - m_thrust_minor;

	return;
}


void EventShapes::calcBrd() {
	/* Function that calculates the event broadening with respect to 
	 * the thrust axis
	 */

	// Depends on the thrust axis
	if (!m_thrust_axis) tradThrust(m_momenta);
	if (!m_thrust_axis) return;

	// Alias thrust axis
	const Vector3f axis = *m_thrust_axis;

	// Calculate broadening in Up and Down hemispheres separately
	double Bu = 0., Bu_norm = 0.;
	double Bd = 0., Bd_norm = 0.;

	for (const auto& x : m_momenta) {
		if (axis.dot(x) > std::numeric_limits<float>::epsilon()) {
			Bu += x.cross(axis).norm();
			Bu_norm += x.norm();
		} else {
			Bd += x.cross(axis).norm();
			Bd_norm += x.norm();
		}
	}

	Bu = Bu_norm > 0. ? Bu / Bu_norm : -1.;
	Bd = Bd_norm > 0. ? Bd / Bd_norm : -1.;

	m_broadening = Bu + Bd;

	return;
}


void EventShapes::calcSph(float r) {
	/* Function that calculates the sphericity tensor and its eigen
	 * vectors and values.
	 */

	// Do not consider the case of zero tracks
	if (m_ntracks < m_min_ntracks) return;

	// Construct the sphericity tensor
	float a11 = 0.; float a12 = 0.; float a13 = 0.;
	float a21 = 0.; float a22 = 0.; float a23 = 0.;
	float a31 = 0.; float a32 = 0.; float a33 = 0.;
	float denom = 0.;

	for (const auto& x : m_momenta){
		float norm = x.norm();
		float r_norm = pow(norm, 2. - r);

		a11 += x.x() * x.x() / r_norm;
		a22 += x.y() * x.y() / r_norm;
		a33 += x.z() * x.z() / r_norm;

		a12 += x.x() * x.y() / r_norm;
		a13 += x.x() * x.z() / r_norm;
		a23 += x.y() * x.z() / r_norm;

		denom += pow(norm, r);
	}

	// Fill symmetric elements of sphericity tensor
	a21 = a12; a31 = a13; a32 = a23;

	double s11 = a11 / denom; double s12 = a12 / denom; double s13 = a13 / denom;
	double s21 = a21 / denom; double s22 = a22 / denom; double s23 = a23 / denom;
	double s31 = a31 / denom; double s32 = a32 / denom; double s33 = a33 / denom;

	// Calculate the eigenvalues
	Matrix3f eigen_problem;
	eigen_problem <<
		s11, s12, s13,
		s21, s22, s23,
		s31, s32, s33;

	SelfAdjointEigenSolver<Matrix3f> eigen_solver(eigen_problem);
	auto eigen_values = eigen_solver.eigenvalues();

	// Compute sphericity variables and set class data members
	if (m_ndims == 3) {
		double S = (eigen_values[1] + eigen_values[0]) * 3./2.;
		double A = eigen_values[0] * 3./2.;
		double C = (eigen_values[2] * eigen_values[1]
					+ eigen_values[2] * eigen_values[0]
					+ eigen_values[1] * eigen_values[0]) * 3.;
		double D = 27. * eigen_values[2] * eigen_values[1] * eigen_values[0];

		m_spher_S = S;
		m_spher_A = A;
		m_spher_C = C;
		m_spher_D = D;

	} else if (m_ndims == 2) {
		double S = eigen_values[1] * 2.;
		double C = eigen_values[1] * eigen_values[2] * 4.;

		m_spher_S = S;
		m_spher_C = C;

	}

}


void EventShapes::calc_all() {
	calcThrust();
	calcThrustMajor();
	calcThrustMinor();
	calcOblateness();
	calcBrd();
	calcSph(2);
}


void EventShapes::calc_thrusts() {
	calcThrust();
	calcThrustMajor();
	calcThrustMinor();
	calcOblateness();
	calcBrd();
}


void EventShapes::calc_sphericities(float r) {
	calcSph(r);
}


EventShapes::~EventShapes() {
	delete m_thrust_axis;
	delete m_thrust_major_axis;
	delete m_thrust_minor_axis;
}