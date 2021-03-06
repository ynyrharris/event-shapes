
#include "EventShapes/EventShapes.h"

#include <Eigen/Dense>

#include <iostream>
#include <math.h>


void m1_base() {

	// Any random vector will do when m = 1!
	// q1 = random
	// q2 = -q1
}


void k2_base(std::vector<Eigen::Vector3f>& V,
	const std::vector<Eigen::Vector3f>& X, unsigned int m
	) {

	// Solve the problem in the plane orthogonal to X[m]
	Eigen::Vector3f zm = X[m].normalized();

	// Construct the linear span E of the plane orthogonal to X[m], H(m)
	Eigen::Vector3f e1 = Eigen::Vector3f::Random();
	e1 = e1 - (e1.dot(zm)) * zm;
	e1 = e1.normalized();

	Eigen::Vector3f e2 = zm.cross(e1);

	// Construct unit vectors orthogonal to the vectors of X in H(m)
	std::vector<Eigen::Vector3f> L;
	L.reserve(2 * (m + 1));
	Eigen::Vector3f hi, li, ri;

	for (unsigned int i = 0; i < m; i++) {
		hi = X[i] - (zm.dot(X[i])) * zm;
		li = zm.cross(hi);
		ri = -li;

		L.push_back(li);
		L.push_back(ri);
	}

	// Sort the rays by angle from e1
	std::sort(L.begin(), L.end(), [&](Eigen::Vector3f a, Eigen::Vector3f b) {
		double phi_a = acos(a.dot(e1));
		phi_a = zm.dot(e1.cross(a)) > 0. ? phi_a : phi_a + M_PI;

		double phi_b = acos(b.dot(e1));
		phi_b = zm.dot(e1.cross(b)) > 0. ? phi_b : phi_b + M_PI;

		return phi_a < phi_b;
	});

	// Take the bisectors of rays as representatives of the regions
	unsigned int L_size = L.size();
	Eigen::Vector3f vi;
	std::vector<Eigen::Vector3f> V_extra;
	V_extra.reserve(L_size * 2);

	for (unsigned int i = 0; i < L_size - 1; i++) {
		vi = L[(i + 1) % L_size] + L[i];

		auto check_redundance = [&] (const Eigen::Vector3f& v) {
			for (unsigned int j = 0; j < V.size(); j++) {

				bool redundant = true;
				for (unsigned int k = 0; k < m; k++) {
					if (v.dot(X[k]) * V[j].dot(X[k]) < std::numeric_limits<float>::epsilon()) {
						redundant = false;
						break;
					}
				}
		
				if (redundant) {
					V.erase(V.begin() + j);
					j--;
				}
			}
		};

		check_redundance(vi);

		V_extra.push_back(vi + 1e-5 * zm);
		V_extra.push_back(vi - 1e-5 * zm);
	}

	for (auto v : V_extra) {
		V.push_back(v);
	}
}


void remove_redundant(std::vector<Eigen::Vector3f>& V,
	const std::vector<Eigen::Vector3f>& X,
	unsigned int m
	) {
	// Remove all but one vectors belonging to the same region.
	// Has to be done to prevent making a list 2^n elements long.

	// Vectors (belong to same region)/(are equivalent) when
	// (u xi) (u' xi) > 0 for all i

	for (unsigned int i = 0; i < V.size() - 1; i++) {
		for (unsigned int j = i + 1; j < V.size(); j++) {

			bool redundant = true;
			for (unsigned int k = 0; k < m; k++) {
				if (V[i].dot(X[k]) * V[j].dot(X[k]) < 1e-7) {
					redundant = false;
					break;
				}
			}

			// Remove V[j] when it is equivalent to V[i]
			if (redundant) {
				V.erase(V.begin() + j);
				j--;
				continue;
			}
		}
	}
}


void split_region() {
	// Split regions C in A into C+ and C-

	// epsilon = ... (smallest separation in X)
	// vi+ = vi + epsilon gm
	// vi- = vi - epsilon gm
}


double calc_lvs_t(const std::vector<Eigen::Vector3f>& V,
	const std::vector<Eigen::Vector3f>& X) {

	double t = 0.;
	for (const Eigen::Vector3f& v : V) {
		double p = 0.;
		double q = 0.;

		Eigen::Vector3f axis = Eigen::Vector3f::Zero();
		for (const Eigen::Vector3f& x : X) {
			v.dot(x) > 0. ? axis += x : axis -= x;
		}

		axis.normalize();

		for (const Eigen::Vector3f& x : X) {
			p += fabs(axis.dot(x));
			q += x.norm();
		}

		double res = p / q;
		if (res > t) {
			t = res;
		}
	}

	return t;
}


void EventShapes::lvs_t() {

	std::vector<Eigen::Vector3f> X = m_momenta;

	// Variables to modify:
	std::vector<Eigen::Vector3f> V;
	V.reserve(pow(X.size(), 2));

	for (unsigned int m = 0; m < X.size(); m++) {

		// Apply initialiser case
		if (m == 0) {
			Eigen::Vector3f q = Eigen::Vector3f::Random();
			q.normalize();
			V.push_back(q);
			V.push_back(-q);
		} else {
			k2_base(V, X, m);
			// remove_redundant(V, X, m);
		}
	}

	double t = calc_lvs_t(V, X);
	m_lvs_t = t;
	m_lvs_lenV = V.size();

	// std::cout << "LVS thrust: " << t << std::endl;

	return;
}