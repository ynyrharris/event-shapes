
#include "EventShapes/EventShapes.h"

#include <Eigen/Dense>

#include <iostream>
#include <math.h>


void m1_base() {

	// Any random vector will do when m = 1!
	// q1 = random
	// q2 = -q1
}


void k2_base() {

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
}


void remove_redundant() {
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
}


void split_region() {
	// Split regions C in A into C+ and C-

	// epsilon = ... (smallest separation in X)
	// vi+ = vi + epsilon gm
	// vi- = vi - epsilon gm
}


void EventShapes::lvs_t() {

	std::cout << "Beginning LVS calculation of thrust" << std::endl;

	// Variables to modify:
	// vector<Eigen::Vector> V

	// for m in range(X)
	// solve problem P(span(xm), m):
	// m = 1: m1_base()
	// k2_base()
	// split_region()
	// remove_redundant()

	return;
}