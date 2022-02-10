package Models;

import java.util.ArrayList;

import inf.v3d.view.Viewer;
import FEM.*;

public class SmallTetraeder {

	public static Structure createStructure() {
		Structure struct = new Structure();
		double lb = 15.0;
		double r = 457.2 / 2000;
		double t = 10.0 / 1000;
		double a = Math.PI * (Math.pow(r, 2) - Math.pow(r - t, 2));
		double e = 2.1e11;
		Constraint c1 = new Constraint(false, false, false);
		Constraint c2 = new Constraint(true, true, false);
		Force f = new Force(0, -20e3, -100e3);

		// create nodes
		Node n1 = struct.addNode(0.0, 0.0, lb * Math.sqrt(2.0 / 3.0));
		Node n2 = struct.addNode(0.0, lb / Math.sqrt(3), 0);
		Node n3 = struct.addNode(-lb / 2, -lb / Math.sqrt(12.0), 0);
		Node n4 = struct.addNode(lb / 2, -lb / Math.sqrt(12.0), 0);
		// apply BCs
		n1.setForce(f);
		// n2.setForce(f); To Test Load vector of Nha
		// n4.setForce(f);
		n2.setConstraint(c1);
		n3.setConstraint(c1);
		n4.setConstraint(c2);
		// create elements
		struct.addElement(e, a, 0, 1);
		struct.addElement(e, a, 0, 2);
		struct.addElement(e, a, 0, 3);
		struct.addElement(e, a, 1, 2);
		struct.addElement(e, a, 2, 3);
		struct.addElement(e, a, 3, 1);
		// return the new structure
		return struct;
	}

	public static void main(String[] args) {
		Viewer viewer = new Viewer();
		Structure struct = createStructure();
		Visualizer viz = new Visualizer(struct, viewer);

		viz.radiusScale(1);
		viz.constrainScale(1);
		viz.forceScale(1 * 3e-5);
		viz.forceRadiusScale(.1);
		viz.displacementScale(3e+3);
		viz.elementForceScale(3e-5);

		viz.drawElement();
		viz.drawConstraints();
		viz.drawForce();

		viewer.setVisible(true);
		//struct.solveLinear();
		//struct.solveArclengthControl(1);
		struct.solveNewtonRaphson();
		//viz.drawDisplacement();// Draw Displacement after solving problem
		//viz.drawElementForce();
		struct.printResult();

	}

}
