
#ifndef INSTANTON_ANALYSIS__EVENT_SHAPES_H
#define INSTANTON_ANALYSIS__EVENT_SHAPES_H

#include "TObject.h"
#include "TRandom.h"

#include <Eigen/Dense>

#include <vector>
#include <utility>

class EventShapes : public TObject {
	/* Class that stores information intermediate to the calculation
	 * of event shapes per event.
	 */

	public:
		EventShapes();
		EventShapes(const std::vector<std::vector<float>>&, unsigned int ndims = 2);
		~EventShapes();

		// Class methods
		void calcSphericity();
		void calcThrustAxis();
		double calcThrustValue(const std::vector<Eigen::Vector3f>&, const Eigen::Vector3f&);
		void compare_calcT();
		const std::pair<Eigen::Vector3f, double> calcT_new(const std::vector<Eigen::Vector3f>&);
		const std::pair<Eigen::Vector3f, double> calcT_orig(const std::vector<Eigen::Vector3f>&);
		void calcThrust();
		void calcThrustMajor();
		void calcThrustMinor();
		void calcOblateness();
		void calcBrd();
		void calcLinSph();

		void calc_all();

		double get_thrust() { return m_thrust; }
		double get_thrust_major() { return m_thrust_major; }
		double get_thrust_minor() { return m_thrust_minor; }
		double get_oblateness() { return m_oblateness; }
		double get_broadening() { return m_broadening; }
		double get_lin_spher_S() { return m_lin_spher_S; }
		double get_lin_spher_A() { return m_lin_spher_A; }
		double get_lin_spher_C() { return m_lin_spher_C; }
		double get_lin_spher_D() { return m_lin_spher_D; }

	private:
		TRandom m_randg;
		// std::vector<TLorentzVector> m_four_momenta;
		// std::vector<TVector3> m_three_momenta;
		std::vector<Eigen::Vector3f> m_three_momenta;
		unsigned int m_ntracks;
		unsigned int m_min_ntracks;
		unsigned int ndims;

		Eigen::Vector3f* m_thrust_axis;
		Eigen::Vector3f* m_thrust_major_axis;
		Eigen::Vector3f* m_thrust_minor_axis;

	public:
		double m_thrust;
		double m_thrust_major;
		double m_thrust_minor;
		double m_oblateness;

		double m_broadening;

		double m_lambda1;
		double m_lambda2;
		double m_lambda3;
		double m_lin_spher_S;
		double m_lin_spher_A;
		double m_lin_spher_C; // 3 jet structure
		double m_lin_spher_D; // 4 jet structure

};

#endif