
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
	const std::vector<Eigen::Vector3f>& X,
	unsigned int m
	) {

	// Construct the linear span E of the plane
	// e1 = unit vector in the plane othogonal to gm
	// e2 = e1 x gm

	// Construct unit vectors orthogonal to the vectors of X in G
	// gi = xi - (zm xi) zm
	// l(gi) = gm x gi
	// r(gi) = -l(gi)
	// L = {li, ri}

	// Sort the Li by angle from e1
	// cos(theta(Li, e1)) = Li e1

	// Take the bisectors of rays as representatives of A
	// vi = (v_{i + 1} + v_i).normalize()

	// Return the bisectors of rays as the solution to the problem in k = 2

	Eigen::Vector3f zm = X[m].normalized();

	// Construct the linear span E of the plane orthogonal to X[m], H(m)
	Eigen::Vector3f e1 = Eigen::Vector3f::Random();
	e1 = e1 - (e1.dot(zm)) * zm;
	e1 = e1.normalized();

	Eigen::Vector3f e2 = zm.cross(e1);

	// Construct unit vectors orthogonal to the vectors of X in H(m)
	std::vector<Eigen::Vector3f> L;
	for (unsigned int i = 0; i < m; i++) {
		Eigen::Vector3f hi = X[i] - (zm.dot(X[i])) * zm;
		Eigen::Vector3f li = zm.cross(hi);
		Eigen::Vector3f ri = -li;

		L.push_back(li.normalized());
		L.push_back(ri.normalized());
	}

	// Sort the rays by angle from e1
	std::sort(L.begin(), L.end(), [&](Eigen::Vector3f a, Eigen::Vector3f b) {
		double phi_a = acos(a.dot(e1));
		phi_a = zm.dot(e1.cross(a)) > 0. ? phi_a : phi_a + M_PI;

		double phi_b = acos(b.dot(e1));
		phi_b = zm.dot(e1.cross(b)) > 0. ? phi_b : phi_b + M_PI;

		return phi_a < phi_b;
	});

	std::cout << "Check sorting of " << L.size() << " elements:" << std::endl;
	for (auto l : L) {
		double phi = acos(l.dot(e1));
		phi = zm.dot(e1.cross(l)) > 0. ? phi : phi + M_PI;
		std::cout << "l dot e1: " << phi << std::endl;
	}

	// Take the bisectors of rays as representatives of the areas
	unsigned int L_size = L.size();
	for (unsigned int i = 0; i < L_size; i++) {
		std::cout << "Indices: " << (i + 1) % L_size << ", " << i << std::endl;
		Eigen::Vector3f vi = (L[(i + 1) % L_size] + L[i]) / 2.;
		vi.normalize();
		V.push_back(vi + 0.000001 * zm);
		V.push_back(vi - 0.000001 * zm);
	}
}


void remove_redundant(std::vector<Eigen::Vector3f>& V,
	const std::vector<Eigen::Vector3f>& X,
	unsigned int m
	) {
	// Remove all but one vectors belonging to the same region.
	// Has to be done to prevent making a list 2^n elements long!

	// Vectors are equivalent when
	// (u xi) (u' xi) > 0 for all i

    // vector<uint> remove_indices
	// for u in v, i: 0 --> n
	//     for u' in v, j: i + 1 --> n
	//         if u' = u continue
	//
	//         for xi in x
	//             if (u xi) (u' xi) < 0
	//                 remove_indices(j)
	//                 break

	std::vector<unsigned int> remove_indices;
	unsigned int V_size = V.size();
	for (unsigned int i = 0; i < V.size() - 1; i++) {
		std::cout << "V size at " << i << ": " << V.size() << std::endl;
		for (unsigned int j = i + 1; j < V.size(); j++) {

			bool redundant = true;
			for (unsigned int k = 0; k < m; k++) {
				if (V[i].dot(X[k]) * V[j].dot(X[k]) < 0) {
					redundant = false;
					break;
				}
			}

			if (redundant) {
				V.erase(V.begin() + j);
				j--;
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


double calc_t(const std::vector<Eigen::Vector3f>& V,
	const std::vector<Eigen::Vector3f>& X) {

	double t = 0.;
	unsigned int i = 0;
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
		std::cout << "Thrust axis w.r.t. region "  << i << ": " << res << std::endl;
		if (res > t) {
			std::cout << "New best thrust axis: " << res << std::endl;
			t = res;
		}
		i++;
	}

	return t;
}


void EventShapes::lvs_t() {

	std::cout << "Beginning LVS calculation of thrust" << std::endl;

	std::vector<Eigen::Vector3f> X = m_three_momenta;

	// Variables to modify:
	std::vector<Eigen::Vector3f> V;

	// for m in range(X)
	// solve problem P(span(xm), m):
	// m = 1: m1_base()
	// k2_base()
	// split_region()
	// remove_redundant()

	for (unsigned int m = 0; m < X.size(); m++) {

		G.push_back(X[m]);

		// Apply m = 1 case when G is still empty
		if (m == 0) {
			Eigen::Vector3f q = Eigen::Vector3f::Random();
			q.normalize();
			V.push_back(q);
			V.push_back(-q);
		} else {
			k2_base(V, X, m);
			remove_redundant(V, X, m);
		}
	}

	double t = calc_t(V, X);

	std::cout << "LSV thrust: " << t << std::endl;

	return;
}