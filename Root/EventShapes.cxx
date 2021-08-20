
#include "EventShapes/EventShapes.h"

#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"

#include <iostream>
#include <math.h>

// ClassImp(EventShapes)

EventShapes::EventShapes()
	: m_thrust_axis(nullptr),
	  m_thrust_major_axis(nullptr),
	  m_thrust_minor_axis(nullptr) {}

EventShapes::EventShapes(const std::vector<std::vector<float>>& momenta) 
	: /*m_four_momenta(momenta),*/
	  m_thrust_axis(nullptr),
	  m_thrust_major_axis(nullptr),
	  m_thrust_minor_axis(nullptr) {

	m_three_momenta.reserve(momenta.size());
	for (auto& p : momenta) {
		TVector3 tvector;
		if (p.size() ==  2) {
			tvector = TVector3(p[0], p[1], 0.);
		} else {
			tvector = TVector3(p[0], p[1], p[2]);
		}
		m_three_momenta.push_back(tvector);
	}

	m_ntracks = m_three_momenta.size();
	m_min_ntracks = 2;

	m_randg = TRandom{};
}

double EventShapes::calcThrustValue(const std::vector<TVector3>& pvec, const TVector3& axis) {
	/* Function that calculates the thrust value of a set of input vectors 
	 * 'pvec' about the thrust axis 'axis'
	 */

	double t_num = 0.;
	double t_denom = 0.;
	for (const auto& p : pvec) {
		t_num += fabs(axis.Dot(p));
		t_denom += p.Mag();
	}
	double thrust = t_num > 0 ? t_num / t_denom : 0.;

	return thrust;
}

void EventShapes::compare_calcT() {

	// Calculate thrust axis by each method
	// std::pair<TVector3, double> axis_new_ret = calcT_new(m_three_momenta);
	std::pair<TVector3, double> axis_orig_ret = calcT_orig(m_three_momenta);

	// TVector3 axis_new = axis_new_ret.first;
	// TVector3 axis_orig = axis_orig_ret.first;

	//<- Commented while many differences are found
	// for (unsigned int i = 0; i < 3; i++) {
	// 	if (fabs(fabs(axis_new[i]) - fabs(axis_orig[i])) > 1e-10) {
	// 		std::cout << "Difference found between two methods of thrust calculation" << std::endl;
	// 		std::cout << "The new and original methods yield, respectively:" << std::endl;
	// 		axis_new.Print();
	// 		axis_orig.Print();
	// 		break;
	// 	}
	// }

	// Set member data variable
	if (axis_orig_ret.first == TVector3()) m_thrust_axis = nullptr;
	else m_thrust_axis = new TVector3(axis_orig_ret.first);
	m_thrust = axis_orig_ret.second;
}

const std::pair<TVector3, double> EventShapes::calcT_new(const std::vector<TVector3>& pvec) {
	/* Function that attempts to calculate the thrust axis more rigorously
	 * based on the intuition that all 2^(n - 1) combinations of input 
	 * vectors and their inverses do not need to be tested, and that the 
	 * combinations that do need to be tested comprise a subset of the 
	 * 2^(n - 1) combinations that can be found.
	 *
	 * This method partitions the set of input vectors into n hemispheres 
	 * by dot product with each of the n input vectors in turn, builds 
	 * the Longest Vector Sum in each hemisphere, and takes the axis 
	 * that corresponds to the greatest value of thrust as the thrust axis.
	 *
	 * This method has produced the same thrust axis as the original 
	 * algorithm, implemented here as calcT_orig().  But it is yet to be 
	 * proven that this method is perfectly accurate, therefore it is 
	 * disfavoured for now.
	 */

	// Do not consider the case of zero tracks
	if (m_ntracks < m_min_ntracks) {
		return std::pair<TVector3, double> (TVector3(), -1.);
	}

	std::vector<TVector3> tvecs;
	std::vector<double> tvals;

	for (unsigned int j = 0; j < pvec.size(); j++) {

		// Transform initial vector into thrust axis
		TVector3 axis (0, 0, 0);

		// Define a hemisphere by dot product with each input vector
		TVector3 init (pvec[j]);

		for (unsigned int i = 0; i < pvec.size(); i++) {
			init.Dot(pvec[i]) >= 0. ? axis += pvec[i] : axis -= pvec[i];
		}
		axis = axis.Unit();

		// Point in the direction of greatest energy flow
		double eflow = 0.;
		for (const auto& p : pvec) eflow += axis.Dot(p);
		if (eflow < 0.) axis = -axis;

		// Get value of thrust about the calculated thrust axis
		double t = calcThrustValue(pvec, axis);

		tvecs.push_back(axis);
		tvals.push_back(t);
	}

	// Find the best thrust axis
	double thrust = 0.;
	TVector3 thrust_axis;
	for (unsigned int i = 0; i < tvecs.size(); i++) {
		if (tvals[i] > thrust) {
			thrust = tvals[i];
			thrust_axis = tvecs[i];
		}
	}

	return std::pair<TVector3, double> (thrust_axis, thrust);
}

const std::pair<TVector3, double> EventShapes::calcT_orig(const std::vector<TVector3>& three_momenta) {
	/* Method based on the original, validated algorithm supplied by 
	 * Deepak Kar and Sukanya Sinha.  Based on the iterative method 
	 * described in the Pythia 6.4 Manual.
	 *
	 * This implementation starts from n random initial axes to build 
	 * a resultant vector as a candidate thrust axis.  The normalised 
	 * resultant that has the greatest value of thrust is taken as the 
	 * thrust axis.
	 */

	// Alias vector of input three-vectors
	std::vector<TVector3> pvec = three_momenta;

	// Do not consider the case of zero tracks
	if (m_ntracks < m_min_ntracks) {
		return std::pair<TVector3, double> (TVector3(), -1.);
	}

	// Get thrust axis
	TVector3 tvec;
	double best_thrust = 0.;

	// Start from multiple random initial axes
	for (unsigned int i = 0; i < pow(pvec.size(), 2); i++) {

		double x, y, z;
		double r = 1.;
		m_randg.Sphere(x, y, z, r);
		TVector3 init (x, y, z);

		// TVector3 init (m_randg.Rndm(), m_randg.Rndm(), m_randg.Rndm());
		TVector3 axis (init);

		// Iterate the axis to local maximum
		double diff = 999.;
		while (diff > 1e-5) {
			TVector3 foo(0, 0, 0);
			for (const auto& p : pvec) {
				axis.Dot(p) > 0 ? foo += p : foo -= p;
			}
			foo = foo.Unit();
			diff = (axis - foo).Mag();
			axis = foo;
		}

		// Keep the axis if it increases the thrust value
		double thrust = calcThrustValue(pvec, axis);
		if (thrust > best_thrust) {
			tvec = axis;
			best_thrust = thrust;
		}
	}

	// Point in the direction of greatest energy flow
	double eflow = 0.;
	for (const auto& p : pvec) eflow += tvec.Dot(p);
	if (eflow <= 0.) tvec = -tvec;

	return std::pair<TVector3, double> (tvec, best_thrust);
}

void EventShapes::calcThrust() {
	/* Function that calculates the thrust proper axis and values.
	 * The class member data are set with the calculated values.
	 */

	std::pair<TVector3, double> pair = calcT_orig(m_three_momenta);

	// Check validity of results
	if (pair.first == TVector3()) { // Calculation unsuccessful
		m_thrust_axis = nullptr;
		m_thrust = -1.;
	} else {
		m_thrust_axis = new TVector3(pair.first);
		m_thrust = pair.second;
	}

}

void EventShapes::calcThrustMajor() {
	/* Function that calculates the thrust major axis and value.
	 * This is the axis in the plane perpendicular to the thrust axis 
	 * along which the energy flow is greatest in that plane.
	 * Phys. Rev. Lett. 43, 830 for a nice discussion.
	 */

	// Do not consider the case of zero tracks
	if (m_ntracks < m_min_ntracks) {
		m_thrust_major_axis = nullptr;
		m_thrust_major = -1.;
		return;
	}

	const std::vector<TVector3> pvec = m_three_momenta;

	// Function depends on the thrust axis
	if (!m_thrust_axis) calcThrust();

	const TVector3 thrust_axis (*m_thrust_axis);

	// Check that the thrust major axis exists
	// i.e. that there are vectors in the plane perp to the thrust axis
	// Unnecessary except when n input vectors < 3
	if (m_ntracks < 3) {
		bool no_thrust_major (true);
		for (const auto& p : pvec) {
			if ((p - p.Dot(thrust_axis) * thrust_axis).Mag2() > 1e-10) {
				no_thrust_major = false;
				break;
			}
		}

		if (no_thrust_major) {
			m_thrust_major_axis = nullptr;
			m_thrust_major = -1.;
			return;
		}
	}

	// Project input vectors into plane perpendicular to thrust axis proper
	std::vector<TVector3> pvec_perp;
	pvec_perp.reserve(pvec.size());
	for (const auto& p : pvec) {
		pvec_perp.push_back(p - p.Dot(thrust_axis) * thrust_axis);
	}

	std::pair<TVector3, double> thrust_major = calcT_orig(pvec_perp);
	// std::cout << "Thrust major value, n, vector: " << thrust_major.second << ", " << pvec.size() << ", ";
	// thrust_major.first.Print();

	TVector3 thrust_major_axis;
	if (thrust_major.first == TVector3()) {
		m_thrust_major_axis = nullptr;
		m_thrust_major = -1.;
		return;
	} else {
		thrust_major_axis = thrust_major.first;
	}

	// Sanity check
	if (fabs(thrust_axis.Dot(thrust_major_axis)) > 1e-10) {
		std::cout << "Major axis not orthogonal to proper axis!" << std::endl;
	}

	// Set class data members
	m_thrust_major_axis = new TVector3(thrust_major.first);
	m_thrust_major = thrust_major.second;
}

void EventShapes::calcThrustMinor() {
	/* Function that calculates the thrust minor axis and value.
	 * This axis is defined as being orthogonal to the thrust and 
	 * thrust major axes.
	 * Phys. Rev. Lett. 43, 830 for a nice discussion.
	 */

	// Do not consider the case where the thrust or thrust major values were invalid
	if (m_thrust == -1. || m_thrust_major == -1.) {
		m_thrust_minor_axis = nullptr;
		m_thrust_minor = -1.;
		return;
	}

	const std::vector<TVector3> pvec = m_three_momenta;

	// This function depends on the thrust and thrust major axes
	if (!m_thrust_major_axis) calcThrustMajor();

	// Alias the existing axes
	const TVector3 thrust_axis = *m_thrust_axis;
	const TVector3 thrust_major_axis = *m_thrust_major_axis;

	TVector3 thrust_minor_axis = thrust_axis.Cross(thrust_major_axis);
	double thrust_minor = calcThrustValue(pvec, thrust_minor_axis);

	// Sanity check
	if (fabs(thrust_axis.Dot(thrust_minor_axis)) > 1e-10 || fabs(thrust_major_axis.Dot(thrust_minor_axis)) > 1e-10) {
		std::cout << "Minor axis not orthogonal to major and proper!" << std::endl;
	}

	// Set class data members
	m_thrust_minor_axis = new TVector3(thrust_minor_axis);
	m_thrust_minor = thrust_minor;
}

void EventShapes::calcOblateness() {
	/* Function that calculates the oblateness about the thrust axis.
	 */

	// This function depends on the major and minor thrust calculations
	if (!m_thrust_minor_axis) calcThrustMinor();

	m_oblateness = m_thrust_major - m_thrust_minor;
}

void EventShapes::calcBrd() {
	/* Function that calculates the event broadening with respect to 
	 * the thrust axis
	 */

	// Do not consider the case of zero tracks
	if (m_ntracks < m_min_ntracks) {
		m_broadening = -1.;
		return;
	}

	// Alias vector of input vectors
	const std::vector<TVector3> pvec = m_three_momenta;

	// Depends on the thrust axis
	if (!m_thrust_axis) calcT_orig(m_three_momenta);

	// Alias thrust axis
	const TVector3 thrust_axis = *m_thrust_axis;

	// Calculate broadening in Up and Down hemispheres separately
	double B_U (0.), B_D (0.);
	double B_norm (0.);
	// std::cout << "New event:" << std::endl;
	for (const auto& p : pvec) {
		// std::cout << "Dot: " << thrust_axis.Dot(p) << std::endl;
		B_norm += p.Mag();
		if (thrust_axis.Dot(p) > 0) {
			B_U += p.Cross(thrust_axis).Mag();
		} else {
			B_D += p.Cross(thrust_axis).Mag();
		}
	}

	// std::cout << "B_D, B_D_norm = " << B_D << ", " << B_norm << std::endl;
	// std::cout << "B_U, B_U_norm = " << B_U << ", " << B_norm << std::endl;
	// std::cout << std::endl;

	B_U = B_U > 0. ? B_U / B_norm : 0.;
	B_D = B_D > 0. ? B_D / B_norm : 0.;

	double B = B_D + B_U;

	// std::cout << "Broadening: " << B << std::endl;

	m_broadening = B;

}

void EventShapes::calcLinSph() {
	/* Function that calculates the sphericity tensor and its eigen
	 * vectors and values.
	 */

	// Alias vector of input vectors
	const std::vector<TVector3> pvec = m_three_momenta;

	// Do not consider the case of zero tracks
	if (m_ntracks < m_min_ntracks) {

		// Set sphericities to dummy values
		m_lin_spher_S = -1.;
		m_lin_spher_A = -1.;
		m_lin_spher_C = -1.;
		m_lin_spher_D = -1.;

		return;
	}

	// Construct the sphericity tensor
	double a11 = 0.; double a12 = 0.; double a13 = 0.;
	double a21 = 0.; double a22 = 0.; double a23 = 0.;
	double a31 = 0.; double a32 = 0.; double a33 = 0.;
	double norm = 0.;

	for (const auto& p : pvec){
		double mod (p.Mag());
		norm += mod;

		a11 += p.Px() * p.Px() / mod;
		a22 += p.Py() * p.Py() / mod;
		a33 += p.Pz() * p.Pz() / mod;

		a12 += p.Px() * p.Py() / mod;
		a13 += p.Px() * p.Pz() / mod;
		a23 += p.Py() * p.Pz() / mod;
	}

	// Fill symmetric elements of sphericity tensor
	a21 = a12; a31 = a13; a32 = a23;

	double s11 = a11 / norm; double s12 = a12 / norm; double s13 = a13 / norm;
	double s21 = a21 / norm; double s22 = a22 / norm; double s23 = a23 / norm;
	double s31 = a31 / norm; double s32 = a32 / norm; double s33 = a33 / norm;

	// Calculate the eigenvalues
	TMatrixD spTensor(0,2,0,2);
	spTensor[0][0] = s11; spTensor[0][1] = s12; spTensor[0][2] = s13;
	spTensor[1][0] = s21; spTensor[1][1] = s22; spTensor[1][2] = s23;
	spTensor[2][0] = s31; spTensor[2][1] = s32; spTensor[2][2] = s33;

	TMatrixDEigen eigenProblem = TMatrixDEigen(spTensor);
	TMatrixD spEigen = eigenProblem.GetEigenValues();

	std::vector<double> eigenvalues;
	for (size_t k = 0; k < 3; ++k) eigenvalues.push_back(spEigen[k][k]);

	// Sort eigenvalues in decreasing order of magnitude, l1 >= l2 >= l3
	std::sort(eigenvalues.begin(), eigenvalues.end(), [](double a, double b) {
		return a > b;
	});

	// Compute sphericity variables and set class data members
	double S = (eigenvalues[1] + eigenvalues[2]) * 3./2.;
	double A = eigenvalues[2] * 3./2.;
	double C = (eigenvalues[0] * eigenvalues[1]
				+ eigenvalues[0] * eigenvalues[2]
				+ eigenvalues[1] * eigenvalues[2]) * 3.;
	double D = 27. * eigenvalues[0] * eigenvalues[1] * eigenvalues[2];

	m_lin_spher_S = S;
	m_lin_spher_A = A;
	m_lin_spher_C = C;
	m_lin_spher_D = D;

	// std::cout << "New event:" << std::endl;
	// std::cout << "S = " << m_lin_spher_S << std::endl;
	// std::cout << "A = " << m_lin_spher_A << std::endl;
	// std::cout << "C = " << m_lin_spher_C << std::endl;
	// std::cout << "D = " << m_lin_spher_D << std::endl;
	// std::cout << std::endl;

}

void EventShapes::calc_all() {
	calcThrust();
	calcThrustMajor();
	calcThrustMinor();
	calcOblateness();
	calcBrd();
	calcLinSph();
}

EventShapes::~EventShapes() {
	delete m_thrust_axis;
	delete m_thrust_major_axis;
	delete m_thrust_minor_axis;
}