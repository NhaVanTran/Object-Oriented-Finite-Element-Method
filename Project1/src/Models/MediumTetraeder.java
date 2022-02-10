package Models;

import FEM.*;
import static java.lang.Math.PI;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import inf.v3d.view.Viewer;

/**
 * The medium Bottrop-Tetraeder.
 * 
 * M. Baitsch
 */
public class MediumTetraeder {

	public static Structure createStructure() {
		Structure struct = new Structure();
		double lb = 15.0;
		double e = 2.1e11;
		double r = 457.2 / 2000;
		double t = 10.0 / 1000;
		double a = PI * (pow(r, 2) - pow(r - t, 2));
		Constraint c1 = new Constraint(false, false, false);

		Node n1 = struct.addNode(0.0, 0.0, lb * sqrt(2.0 / 3.0));
		Node n2 = struct.addNode(0.0, lb / sqrt(3), 0);
		Node n3 = struct.addNode(-lb / 2.0, -lb / sqrt(12.0), 0);
		Node n4 = struct.addNode(lb / 2.0, -lb / sqrt(12.0), 0);
		Node n5 = struct.addNode(0.0, lb / sqrt(3) + (lb / sqrt(3)), 0 - lb
				* sqrt(2.0 / 3.0));
		Node n6 = struct.addNode(-lb / 2.0 - lb / 2.0, -lb / sqrt(12.0) - lb
				/ sqrt(12.0), 0 - lb * sqrt(2.0 / 3.0));
		Node n7 = struct.addNode(lb / 2.0 + lb - lb / 2.0, -lb / sqrt(12.0)
				- lb / sqrt(12.0), 0 - lb * sqrt(2.0 / 3.0));
		Node n8 = struct.addNode(0.0 - lb / 2.0,
				lb / sqrt(3) - lb / sqrt(12.0), 0 - lb * sqrt(2.0 / 3.0));
		Node n9 = struct.addNode(-lb / 2.0 + lb - lb / 2.0, -lb / sqrt(12.0)
				- lb / sqrt(12.0), 0 - lb * sqrt(2.0 / 3.0));
		Node n10 = struct.addNode(0.0 + lb - lb / 2.0, lb / sqrt(3) - lb
				/ sqrt(12.0), 0 - lb * sqrt(2.0 / 3.0));

		Force f2 = new Force(-45.4, -71.93, -90.09);
		Force f3 = new Force(139.42, 48.92, -205.98);
		Force f4 = new Force(-19.91, 19.36, -21.37);
		Force f8 = new Force(134.01, 44.59, -109.7);
		Force f9 = new Force(-74.59, -25.81, -82.2);
		Force f10 = new Force(-146.26, 80.42, -144.34);

		n2.setForce(f2);
		n3.setForce(f3);
		n4.setForce(f4);
		n8.setForce(f8);
		n9.setForce(f9);
		n10.setForce(f10);

		n5.setConstraint(c1);
		n6.setConstraint(c1);
		n7.setConstraint(c1);

		struct.addElement(e, a, 0, 1);
		struct.addElement(e, a, 0, 2);
		struct.addElement(e, a, 0, 3);
		struct.addElement(e, a, 1, 2);
		struct.addElement(e, a, 2, 3);
		struct.addElement(e, a, 3, 1);
		struct.addElement(e, a, 2, 7);
		struct.addElement(e, a, 2, 5);
		struct.addElement(e, a, 2, 8);
		struct.addElement(e, a, 7, 5);
		struct.addElement(e, a, 5, 8);
		struct.addElement(e, a, 8, 7);
		struct.addElement(e, a, 3, 9);
		struct.addElement(e, a, 3, 8);
		struct.addElement(e, a, 3, 6);
		struct.addElement(e, a, 9, 8);
		struct.addElement(e, a, 8, 6);
		struct.addElement(e, a, 6, 9);
		struct.addElement(e, a, 1, 4);
		struct.addElement(e, a, 1, 7);
		struct.addElement(e, a, 1, 9);
		struct.addElement(e, a, 4, 7);
		struct.addElement(e, a, 7, 9);
		struct.addElement(e, a, 9, 4);

		return struct;
	}

	public static void main(String[] args) {
		Viewer viewer = new Viewer();
		Structure struct = createStructure();
		Visualizer viz = new Visualizer(struct, viewer);

		viz.radiusScale(3);
		viz.constrainScale(2);
		viz.forceScale(4 * 3e-2);
		viz.forceRadiusScale(.1);
		viz.displacementScale(4e+5);
		viz.elementForceScale(1e-2);

		//
		viz.drawElement();
		viz.drawConstraints();
		viz.drawForce();

		viewer.setVisible(true);
		struct.solveLinear();
		// struct.SolveNonLinearArcLengthControl(5);
		viz.drawDisplacement();// Draw Displacement after solving problem
		viz.drawElementForce();
		struct.printResult();

	}
}
